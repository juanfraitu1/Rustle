#!/usr/bin/env python3
"""Truth builders: annotation node tables, the adjudicated two-annotation truth, protein-space families, and the
protein referee (wave 7, 2026-09-24; old scripts at git tag `notebook-2026-09-24`).

Old -> new (every old command keeps its arguments):
  annotation_nodes.py {refseq|cat|ensembl} GFF CONTIGS OUT   -> truth.py nodes {refseq|cat|ensembl} GFF CONTIGS OUT
  adjudicated_truth.py build --out DIR ...                   -> truth.py adjudicated --out DIR ...
                                                                (B13 fixed: the old build raised TypeError in mcl_port
                                                                since 21d6c5c9 — int node ids; output now equals the
                                                                numpy-era build byte for byte)
  protein_families.py build --nodes N --genome FA ...        -> truth.py protein --nodes N --genome FA ...
  (inline in rna_truth_from_protein.py / soto_vs_us_referee.py)
                                                             -> truth.py protein-referee --gff G --genome FA --chrom C
                                                                --out PREFIX [--threads N]  (writes PREFIX.families.tsv,
                                                                `Gene Name` / `Family ID`, ids PF<i>)
  adjudicated_truth.py score / protein_families.py score     -> score.py adjudicated / score.py protein
Library names:
  protein_families.excluded / pair_hsps / edges_from / load_genes -> truth.excluded / pair_hsps / edges_from / load_genes
  protein_edge_gap.protein_edges                                  -> truth.protein_edges
  the two blastp all-vs-all blocks (protein_edge_gap, protein_families) -> truth.blastp_all_vs_all

⚠ Two protein sets are both called "§6ko" and are NOT the same (audit D13/D14): `truth.py protein` translates the
annotation_nodes `.cds.tsv` with `lib.translate_phased` and keeps per ORDERED pair the greedy HSP union on the longer
protein (`edges_from`); the referee (`protein-referee`, `score.py referee|rna-ceiling|edge-gap`) takes the longest
CDS per `gene=` symbol, translates with `lib.translate_refseq` and uses `protein_edges` (unordered pairs, HSPs whose
QUERY is the longer protein, closed intervals). Both are kept as they were.

Only the standard library is imported at module top; pysam is imported inside the functions that need it.
"""
import argparse
import bisect
import collections
import csv
import itertools
import os
import re
import subprocess
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import lib  # noqa: E402

BLAST = os.environ.get('BLAST_BIN', '/home/juanfra/miniforge3/envs/blast/bin')
MIN_COV = 0.30


# ================================================================ protein space (was protein_families / protein_edge_gap)
def excluded(biotype, rule):
    """r1: pseudogenes; r2: r1 + V(D)J recombining antigen-receptor segments (AN revision r2)."""
    if rule >= 1 and "pseudogene" in biotype:
        return True
    if rule >= 2 and (any(x in biotype for x in ("V_segment", "D_segment", "J_segment", "C_region"))
                      or biotype.startswith(("IG_", "TR_"))):
        return True
    return False


def load_genes(nodes, contigs):
    """annotation_nodes tables -> {idx: gene dict with 0-based (start, end, phase) CDS segments}."""
    names = {r["idx"]: (r["name"], r["biotype"]) for r in csv.DictReader(open(nodes + ".names.tsv"), delimiter="\t")}
    chrom = {r["idx"]: r["chrom"] for r in csv.DictReader(open(nodes), delimiter="\t")}
    genes = {}
    for r in csv.DictReader(open(nodes + ".cds.tsv"), delimiter="\t"):
        if chrom[r["idx"]] not in contigs:
            continue
        segs = [(int(x.split("-")[0]), int(x.split("-")[1].split(":")[0]), int(x.split(":")[1])) for x in r["cds"].split(",")]
        genes[r["idx"]] = {"idx": r["idx"], "chrom": chrom[r["idx"]], "strand": r["strand"], "cds": segs,
                           "name": names[r["idx"]][0], "biotype": names[r["idx"]][1]}
    return genes


def blastp_all_vs_all(faa, prefix, db_suffix, threads):
    """All-vs-all blastp -evalue 1e-5 of `faa`, cached at PREFIX.blastp.tsv (on the file's EXISTENCE, as before).
    db_suffix was '_db' in protein_edge_gap and '_protdb' in protein_families; everything else is the same command."""
    bl = prefix + '.blastp.tsv'
    if not os.path.exists(bl):
        subprocess.run([BLAST + '/makeblastdb', '-dbtype', 'prot', '-in', faa, '-out', prefix + db_suffix],
                       stdout=subprocess.DEVNULL, check=True)
        with open(bl + '.tmp', 'w') as fh:
            subprocess.run([BLAST + '/blastp', '-query', faa, '-db', prefix + db_suffix, '-evalue', '1e-5',
                            '-max_target_seqs', '100000', '-num_threads', str(threads), '-outfmt',
                            '6 qseqid sseqid nident length qstart qend sstart send bitscore'],
                           stdout=fh, check=True)
        os.replace(bl + '.tmp', bl)
    return bl


def pair_hsps(path, plen):
    """(q, s) -> list of (bitscore, nident, length, q0, q1, s0, s1) from outfmt 6 lines."""
    hs = collections.defaultdict(list)
    for line in open(path):
        q, s, nid, ln, q0, q1, s0, s1, bits = line.rstrip("\n").split("\t")
        if q != s:
            hs[(q, s)].append((float(bits), int(nid), int(ln), int(q0) - 1, int(q1), int(s0) - 1, int(s1)))
    return hs


def edges_from(hs, plen, min_ident):
    """§6ko edges of `truth.py protein` (ordered pairs; greedy non-overlapping HSPs on the longer protein, half-open
    union >= 0.30 of the longer, identity >= min_ident; weight identity x coverage)."""
    best = {}
    for (q, s), rows in hs.items():
        longer_is_q = plen[q] >= plen[s]
        taken, L, N = [], 0, 0
        for bits, nid, ln, q0, q1, s0, s1 in sorted(rows, reverse=True):
            iv = (q0, q1) if longer_is_q else (s0, s1)
            if any(iv[0] < y and x < iv[1] for x, y in taken):
                continue
            taken.append(iv)
            L += ln
            N += nid
        cov = sum(y - x for x, y in lib.merge(taken)) / max(plen[q], plen[s])
        ident = N / L if L else 0.0
        if cov >= MIN_COV and ident >= min_ident:
            k = (min(q, s, key=int), max(q, s, key=int))
            w = ident * min(cov, 1.0)
            if w > best.get(k, (0.0,))[0]:
                best[k] = (w, ident, cov)
    return best


def protein_edges(faa, out, plen, threads):
    """§6ko edges of the protein REFEREE (protein_edge_gap.protein_edges): unordered pairs; greedy non-overlapping
    HSPs by bitscore whose QUERY is the longer protein, closed intervals, edge iff they cover >= 0.30 of it."""
    bl = blastp_all_vs_all(faa, out, '_db', threads)
    hs = collections.defaultdict(list)
    for line in open(bl):
        q, s, nid, ln, q0, q1, s0, s1, bits = line.rstrip('\n').split('\t')
        if q != s:
            hs[tuple(sorted((q, s)))].append((float(bits), int(q0), int(q1), q, s))
    ed = set()
    for (a, b), v in hs.items():
        longer = max(plen.get(a, 0), plen.get(b, 0))
        if not longer:
            continue
        # greedy non-overlapping HSPs by bitscore, projected onto the LONGER protein
        taken = []
        for bits, q0, q1, q, s in sorted(v, reverse=True):
            if plen.get(q, 0) != longer:
                continue
            lo, hi = min(q0, q1), max(q0, q1)
            if all(hi < t0 or lo > t1 for t0, t1 in taken):
                taken.append((lo, hi))
        cov = sum(hi - lo + 1 for lo, hi in taken) / longer
        if cov >= 0.30:
            ed.add((a, b))
    return ed


def write_proteins(fa, chrom, cds, faa):
    """Translate each gene's longest CDS (lib.translate_refseq), keep >= 10 aa; returns {gene: protein length}."""
    plen = {}
    with open(faa, 'w') as fh:
        for g, (st, segs) in cds.items():
            p = lib.translate_refseq(fa, chrom, st, segs)
            if len(p) >= 10:
                plen[g] = len(p); fh.write(f'>{g}\n{p}\n')
    return plen


def protein_referee(gff, fa, chrom, out, threads, reuse_faa=False):
    """The protein REFEREE on one chromosome (§6ko rule, not re-tuned): longest CDS per gene, r2 exclusions
    (pseudogenes + V(D)J segments), translated, all-vs-all blastp, `protein_edges`, MCL I = 2.8.
    Returns {'PF<i>': sorted members} for MCL clusters with >= 2 members (i = MCL cluster index).

    `fa` is an open pysam.FastaFile. reuse_faa=True reuses an existing OUT.proteins.faa (what soto_vs_us_referee did,
    reading protein lengths back from it); otherwise the .faa is always rewritten (rna_truth_from_protein). The blastp
    table is cached on existence either way."""
    cds = lib.longest_cds(gff, chrom)
    # §6ko rule r2, applied verbatim: exclude pseudogenes AND V(D)J recombining segments. Without it an
    # IGKV cluster merges into one 75-member "family" and dominates every pair-weighted statistic.
    bt = lib.gene_biotypes(gff, chrom)
    cds = {g: v for g, v in cds.items() if not excluded(bt.get(g, ''), 2)}
    faa = out + '.proteins.faa'
    plen = {}
    if not (reuse_faa and os.path.exists(faa)):
        plen = write_proteins(fa, chrom, cds, faa)
    if reuse_faa and not plen:
        n = None
        for line in open(faa):
            if line.startswith('>'):
                n = line[1:].strip()
            elif n:
                plen[n] = len(line.strip()); n = None
    pedges = protein_edges(faa, out, plen, threads)
    fams = lib.mcl({(x, y): 1.0 for x, y in pedges}, inflation=2.8)   # {(a,b): weight}
    truth = {}
    for i, mem in enumerate(fams):
        m = sorted(set(mem))
        if len(m) >= 2:
            truth[f'PF{i}'] = m
    return truth


def cmd_protein_referee(a):
    import pysam
    fa = pysam.FastaFile(a.genome)
    truth = protein_referee(a.gff, fa, a.chrom, a.out, a.threads)
    with open(a.out + '.families.tsv', 'w') as fh:
        fh.write('Gene Name\tFamily ID\n')
        for fid, mem in truth.items():
            for g in mem:
                fh.write(f'{g}\t{fid}\n')
    print(f'{a.chrom}: protein referee {len(truth)} families over {sum(map(len, truth.values()))} genes '
          f'-> {a.out}.families.tsv')


def cmd_protein(a):
    """Prereg Addendum AN build (was `protein_families.py build`)."""
    import pysam
    contigs = set(a.contigs.split(","))
    genome = pysam.FastaFile(a.genome)
    genes = load_genes(a.nodes, contigs)
    rule = max(a.rule, 1 if a.no_pseudogenes else 0)
    genes = {k: g for k, g in genes.items() if not excluded(g["biotype"], rule)}
    faa = a.out + ".proteins.faa"
    plen = {}
    with open(faa, "w") as fh:
        for k, g in genes.items():
            p = lib.translate_phased(genome, g["chrom"], g["strand"], g["cds"])
            if len(p) >= 10:
                plen[k] = len(p)
                fh.write(f">{k}\n{p}\n")
    bl = blastp_all_vs_all(faa, a.out, "_protdb", a.threads)
    E = edges_from(pair_hsps(bl, plen), plen, a.min_ident)
    fams = [c for c in lib.mcl({k: v[0] for k, v in E.items()}) if len(c) >= 2]
    fams.sort(key=len, reverse=True)
    tag = a.out if a.min_ident == 0 else f"{a.out}.i{a.min_ident:.2f}"
    with open(tag + ".families.tsv", "w") as fh:
        fh.write("family_id\tidx\tname\tbiotype\tchrom\tstrand\tcds\n")
        for i, c in enumerate(fams):
            for k in sorted(c, key=int):
                g = genes[k]
                fh.write(f"P{i}\t{k}\t{g['name']}\t{g['biotype']}\t{g['chrom']}\t{g['strand']}\t"
                         f"{','.join(f'{x}-{y}' for x, y, _ in g['cds'])}\n")
    with open(tag + ".edges.tsv", "w") as fh:
        fh.write("u\tv\tweight\tidentity\tcoverage\n")
        for (u, v), (w, i, c) in sorted(E.items(), key=lambda kv: (int(kv[0][0]), int(kv[0][1]))):
            fh.write(f"{u}\t{v}\t{w:.4f}\t{i:.4f}\t{c:.4f}\n")
    print(f"{tag}: proteins {len(plen)}; edges {len(E)}; families {len(fams)} ({sum(map(len, fams))} genes, "
          f"largest {len(fams[0]) if fams else 0})")


# ================================================================ node tables (was annotation_nodes.py)
def nodes_load_genes(gff, contigs):
    """RefSeq gene/pseudogene `Name=` -> (chrom, start0, end, strand), last record wins; exon union from exon `gene=`."""
    genes, exons = {}, collections.defaultdict(list)
    for line in open(gff):
        if line.startswith("#"):
            continue
        f = line.rstrip("\n").split("\t")
        if len(f) < 9 or f[0] not in contigs:
            continue
        if f[2] in ("gene", "pseudogene"):
            n = re.search(r"(?:^|;)Name=([^;]+)", f[8])
            if n:
                genes[n.group(1)] = (f[0], int(f[3]) - 1, int(f[4]), f[6])
        elif f[2] == "exon":
            g = re.search(r"(?:^|;)gene=([^;]+)", f[8])
            if g:
                exons[g.group(1)].append((int(f[3]) - 1, int(f[4])))
    return genes, {n: lib.merge(exons.get(n) or [(g[1], g[2])]) for n, g in genes.items()}


def longest(cds_by_tx):
    """transcript -> [(start0, end, phase)] -> the longest CDS's segments (ties: first transcript id)."""
    best = None
    for t in sorted(cds_by_tx):
        segs = cds_by_tx[t]
        L = sum(e - s for s, e, _ in segs)
        if best is None or L > best[0]:
            best = (L, sorted(segs))
    return best[1] if best else None


def nodes_refseq(gff, contigs):
    genes, exons = nodes_load_genes(gff, contigs)
    biotype, n_records = {}, collections.Counter()
    cds = collections.defaultdict(lambda: collections.defaultdict(list))
    for line in open(gff):
        f = line.rstrip("\n").split("\t")
        if len(f) < 9 or f[0] not in contigs:
            continue
        if f[2] in ("gene", "pseudogene"):
            a = lib.gff_attrs(f[8])
            if "Name" in a:
                n_records[a["Name"]] += 1
                biotype[a["Name"]] = a.get("gene_biotype", f[2])
        elif f[2] == "CDS":
            a = lib.gff_attrs(f[8])
            if "gene" in a:
                cds[a["gene"]][a.get("Parent", "?")].append((int(f[3]) - 1, int(f[4]), int(f[7]) if f[7] in "012" else 0))
    dup = sum(1 for v in n_records.values() if v > 1)
    print(f"refseq: {len(genes)} named genes ({dup} names on > 1 record; last record kept, as load_genes)")
    return [(c, s, e, st, exons[n], n, biotype.get(n, "?"), longest(cds[n]) if n in cds else None)
            for n, (c, s, e, st) in genes.items()]


def nodes_gencode_like(gff, contigs, gene_types, rename=None):
    genes, parent, exons = {}, {}, collections.defaultdict(list)
    cds = collections.defaultdict(list)
    for line in open(gff):
        if line.startswith("#"):
            continue
        f = line.rstrip("\n").split("\t")
        if len(f) < 9:
            continue
        chrom = rename(f[0]) if rename else f[0]
        if chrom not in contigs:
            continue
        a = lib.gff_attrs(f[8])
        if f[2] in gene_types:
            genes[a["ID"]] = (chrom, int(f[3]) - 1, int(f[4]), f[6], a.get("Name", a.get("gene_id", a["ID"])),
                              a.get("gene_biotype", a.get("gene_type", a.get("biotype", "?"))))
        elif f[2] == "exon":
            exons[a.get("Parent")].append((int(f[3]) - 1, int(f[4])))
        elif f[2] == "CDS":
            cds[a.get("Parent")].append((int(f[3]) - 1, int(f[4]), int(f[7]) if f[7] in "012" else 0))
        elif "ID" in a and "Parent" in a:
            parent[a["ID"]] = a["Parent"]
    gex = collections.defaultdict(list)
    gcds = collections.defaultdict(dict)
    for t, ex in exons.items():
        if parent.get(t) in genes:
            gex[parent[t]] += ex
    for t, segs in cds.items():
        if parent.get(t) in genes:
            gcds[parent[t]][t] = segs
    print(f"{len(genes)} genes, {sum(1 for g in genes if gex.get(g))} with exons, {len(gcds)} with CDS")
    return [(c, s, e, st, lib.merge(gex.get(g) or [(s, e)]), name, bt, longest(gcds[g]) if g in gcds else None)
            for g, (c, s, e, st, name, bt) in genes.items()]


def nodes_cat(gff, contigs):
    return nodes_gencode_like(gff, contigs, ("gene",))


def nodes_ensembl(gff, contigs):
    return nodes_gencode_like(gff, contigs, ("gene", "ncRNA_gene", "pseudogene"),
                              rename=lambda c: c if c.startswith("chr") else "chr" + c)


def cmd_nodes(a):
    """Prereg Addenda AI/AJ/AK gene-level node tables (was annotation_nodes.py).
    refseq:  gene+pseudogene records, exon union from exon lines' `gene=`.
    cat:     CAT/Liftoff GENCODE `gene` records, exon union from the exons of the gene's transcripts (exon -> Parent
             transcript -> Parent gene); a gene with no exon keeps its span.
    ensembl: Ensembl gene / ncRNA_gene / pseudogene records (contigs named N, renamed chrN), exon union from the exons
             of every feature whose Parent is the gene.
    Writes the node table (idx order = sorted by chrom, start), `<OUT>.names.tsv` (idx, name, biotype) and
    `<OUT>.cds.tsv` (idx, strand, CDS segments `start-end:phase` of the gene's transcript with the longest CDS)."""
    kind, gff, contigs, out = a.kind, a.gff, a.contigs, a.out
    rows = sorted({"refseq": nodes_refseq, "cat": nodes_cat, "ensembl": nodes_ensembl}[kind](gff, set(contigs.split(","))),
                  key=lambda r: r[:7])
    with open(out, "w") as fh, open(out + ".names.tsv", "w") as fn, open(out + ".cds.tsv", "w") as fc:
        fh.write("idx\tchrom\tstart\tend\tstrand\tn_exon\tn_reads\texons\n")
        fn.write("idx\tname\tbiotype\n")
        fc.write("idx\tstrand\tcds\n")
        for i, (c, s, e, st, ex, name, bt, cd) in enumerate(rows):
            fh.write(f"{i}\t{c}\t{s}\t{e}\t{st}\t{len(ex)}\t0\t{','.join(f'{x}-{y}' for x, y in ex)}\n")
            fn.write(f"{i}\t{name}\t{bt}\n")
            if cd:
                fc.write(f"{i}\t{st}\t{','.join(f'{x}-{y}:{p}' for x, y, p in cd)}\n")
    print(f"wrote {out}: {len(rows)} nodes")


# ================================================================ adjudicated truth (was adjudicated_truth.py build)
MIN_ID, AK_MIN_COV, SD_SLACK = 0.70, 0.30, 1000


def load_records(name, spec, contigs):
    nodes_tsv, tag = spec.split(":", 1)
    names = {r["idx"]: (r["name"], r["biotype"]) for r in csv.DictReader(open(nodes_tsv + ".names.tsv"), delimiter="\t")}
    cds = {r["idx"]: (r["strand"], [(int(a), int(b.split(":")[0]), int(b.split(":")[1]))
                                    for a, b in (x.split("-") for x in r["cds"].split(","))])
           for r in csv.DictReader(open(nodes_tsv + ".cds.tsv"), delimiter="\t")}
    fam = {}
    for r in csv.DictReader(open(tag + ".clusters.tsv"), delimiter="\t"):
        fam[f"{r['chrom']}:{r['start']}-{r['end']}"] = r["cluster_id"]
    rep = {}
    if os.path.exists(tag + ".loci.tsv"):
        for r in csv.DictReader(open(tag + ".loci.tsv"), delimiter="\t"):
            rep[r["annotation"]] = r["representative"]
    recs = []
    for r in csv.DictReader(open(nodes_tsv), delimiter="\t"):
        if r["chrom"] not in contigs:
            continue
        key = f"{r['chrom']}:{int(r['start']) + 1}-{r['end']}"
        f = fam.get(key) or fam.get(rep.get(key, ""))
        ex = [tuple(map(int, b.split("-"))) for b in r["exons"].split(",")]
        recs.append({"ann": name, "idx": r["idx"], "chrom": r["chrom"], "start": int(r["start"]), "end": int(r["end"]),
                     "exons": ex, "name": names[r["idx"]][0], "biotype": names[r["idx"]][1], "family": f,
                     "cds": cds.get(r["idx"])})
    return recs


def load_sedef(path, contigs):
    rows = collections.defaultdict(list)
    for line in open(path):
        f = line.rstrip("\n").split("\t")
        if f[0] in contigs and f[9] in contigs:
            rows[f[0]].append((int(f[1]), int(f[2]), f[9], int(f[10]), int(f[11]), f[13]))
    idx = {}
    for c, v in rows.items():
        v.sort()
        idx[c] = (v, [x[0] for x in v], max(x[1] - x[0] for x in v))
    return idx


def exon_hull_in(exons, s, e):
    parts = [(max(a, s), min(b, e)) for a, b in exons if a < e and s < b]
    return (min(p[0] for p in parts), max(p[1] for p in parts)) if parts else None


def sd_evidence(sd, u, v):
    if u["chrom"] not in sd:
        return False
    rows, starts, maxlen = sd[u["chrom"]]
    lo = bisect.bisect_left(starts, u["start"] - maxlen)
    hi = bisect.bisect_left(starts, u["end"])
    for a0, a1, oc, b0, b1, strand in rows[lo:hi]:
        if a1 <= u["start"] or oc != v["chrom"] or b1 <= v["start"] or b0 >= v["end"]:
            continue
        hu = exon_hull_in(u["exons"], a0, a1)
        hv = exon_hull_in(v["exons"], b0, b1)
        if not hu or not hv:
            continue
        ratio = (b1 - b0) / max(1, a1 - a0)
        if strand == "-":
            y = sorted((b1 - (hu[0] - a0) * ratio, b1 - (hu[1] - a0) * ratio))
        else:
            y = (b0 + (hu[0] - a0) * ratio, b0 + (hu[1] - a0) * ratio)
        tol = abs((b1 - b0) - (a1 - a0)) + SD_SLACK
        if y[0] - tol < hv[1] and hv[0] < y[1] + tol:
            return True
    return False


def blast_evidence(a, loci, involved, genome):
    """X(u -> v): dc-megablast of u's exon-union sequence (both annotations; soft-masked repeats do not seed) against v's
    genomic span; non-overlapping HSPs (on the query, best nident first) sum to >= 300 aligned bp at identity >= 0.70 and
    cover >= 0.30 of u's exonic length. Returns the set of (u, v) with X."""
    q_fa, t_fa, db, out = (f"{a.out}/exons.fa", f"{a.out}/spans.fa", f"{a.out}/spans_db", f"{a.out}/blast.tsv")
    qlen = {}
    with open(q_fa, "w") as fq, open(t_fa, "w") as ft:
        for k in involved:
            l = loci[k]
            seq = "".join(genome.fetch(l["chrom"], s, e) for s, e in l["exons"])
            qlen[k] = len(seq)
            fq.write(f">L{k}\n{seq}\n")
            ft.write(f">L{k}\n{genome.fetch(l['chrom'], l['start'], l['end'])}\n")
    if not os.path.exists(out):
        subprocess.run([a.blast_bin + "/makeblastdb", "-dbtype", "nucl", "-in", t_fa, "-out", db],
                       stdout=subprocess.DEVNULL, check=True)
        with open(out + ".tmp", "w") as fh:
            subprocess.run([a.blast_bin + "/blastn", "-task", "dc-megablast", "-query", q_fa, "-db", db, "-lcase_masking",
                            "-evalue", "1e-5", "-max_target_seqs", "100000", "-num_threads", str(a.threads),
                            "-outfmt", "6 qseqid sseqid qstart qend nident length"], stdout=fh, check=True)
        os.replace(out + ".tmp", out)
    hsps = collections.defaultdict(list)
    for line in open(out):
        qs_, ss_, q0, q1, nid, ln = line.split("\t")
        if qs_ != ss_:
            q0, q1 = sorted((int(q0), int(q1)))
            hsps[(int(qs_[1:]), int(ss_[1:]))].append((int(nid), int(ln), q0 - 1, q1))
    hit = set()
    for (u, v), hs in hsps.items():
        taken, L, N = [], 0, 0
        for nid, ln, q0, q1 in sorted(hs, reverse=True):
            if any(q0 < y and x < q1 for x, y in taken):
                continue
            taken.append((q0, q1))
            L += ln
            N += nid
        cov = sum(y - x for x, y in lib.merge(taken)) / max(1, qlen[u])
        if L >= 300 and N / L >= MIN_ID and cov >= AK_MIN_COV:
            hit.add((u, v))
    return hit


def overlap_groups(idx_list, recs):
    """Union-find over records whose exon unions share >= 1 bp (the `--merge-overlapping-loci` rule)."""
    uf = lib.UF()
    by_chrom = collections.defaultdict(list)
    for i in idx_list:
        uf.find(i)
        for s, e in recs[i]["exons"]:
            by_chrom[recs[i]["chrom"]].append((s, e, i))
    for blocks in by_chrom.values():
        blocks.sort()
        cur_end, cur_i = -1, None
        for s, e, i in blocks:
            if cur_i is not None and s < cur_end:
                uf.union(cur_i, i)
            if e > cur_end:
                cur_end, cur_i = e, i
    groups = collections.defaultdict(list)
    for i in idx_list:
        groups[uf.find(i)].append(i)
    return list(groups.values())


def joint_groups(recs, A, Bn, mode):
    """union (AL): one union-find over both annotations — chains distinct genes through overlapping models of the other
    annotation (§6kn). matched (AM): each annotation's own loci, then a RefSeq locus and a GENCODE locus are one joint
    locus iff each is the other's best exonic-overlap partner; everything else stays alone."""
    if mode == "union":
        return overlap_groups(list(range(len(recs))), recs)
    per = {n: overlap_groups([i for i, r in enumerate(recs) if r["ann"] == n], recs) for n in (A, Bn)}
    blocks = collections.defaultdict(list)
    for n in (A, Bn):
        for k, members in enumerate(per[n]):
            for s, e in lib.merge([b for i in members for b in recs[i]["exons"]]):
                blocks[recs[members[0]]["chrom"]].append((s, e, n, k))
    ovc = collections.Counter()
    for bl in blocks.values():
        bl.sort()
        active = []
        for s, e, n, k in bl:
            active = [x for x in active if x[1] > s]
            for s2, e2, n2, k2 in active:
                if n2 != n:
                    key = (k, k2) if n == A else (k2, k)
                    ovc[key] += min(e, e2) - s
            active.append((s, e, n, k))
    best_a, best_b = {}, {}
    for (ka, kb), bp in ovc.items():
        if bp > best_a.get(ka, (0, None))[0]:
            best_a[ka] = (bp, kb)
        if bp > best_b.get(kb, (0, None))[0]:
            best_b[kb] = (bp, ka)
    matched_b, groups = set(), []
    for ka, members in enumerate(per[A]):
        kb = best_a.get(ka, (0, None))[1]
        if kb is not None and best_b.get(kb, (0, None))[1] == ka:
            groups.append(members + per[Bn][kb])
            matched_b.add(kb)
        else:
            groups.append(members)
    groups += [m for kb, m in enumerate(per[Bn]) if kb not in matched_b]
    return groups


def cmd_adjudicated(a):
    """Prereg Addendum AK: adjudicated two-annotation ground truth (was `adjudicated_truth.py build`).
    NODES_TSV from `truth.py nodes` (with .names.tsv / .cds.tsv); TAG_PREFIX = the E1 construction's `<tag>` (reads
    `<tag>.clusters.tsv`, `<tag>.loci.tsv`, `<tag>.graph.tsv`). Exactly two --ann. Joint loci: records of both
    annotations whose exon unions share >= 1 bp. Opinion per annotation: SAME / DIFF / NONE. AGREED TRUE = SAME in
    both; DISPUTED = SAME in exactly one -> TRUE if evidence (dc-megablast X or SEDEF S), FALSE if both loci coding and
    no evidence, else UNSCORED. AK-0 gates are printed. Writes DIR/loci.tsv, DIR/edges.tsv, DIR/pairs.tsv,
    DIR/clusters.tsv. The dc-megablast table is cached at DIR/blast.tsv (on existence)."""
    import pysam
    os.makedirs(a.out, exist_ok=True)
    contigs = set(a.contigs.split(","))
    genome = pysam.FastaFile(a.genome)
    anns = [x.split("=", 1) for x in a.ann]
    assert len(anns) == 2, "exactly two --ann"
    (A, specA), (Bn, specB) = anns
    recs = [r for name, spec in anns for r in load_records(name, spec, contigs)]
    groups = joint_groups(recs, A, Bn, a.joint_loci)
    loci, rec_locus = [], {}
    for members in groups:
        rs = [recs[i] for i in members]
        loci.append({"chrom": rs[0]["chrom"], "start": min(r["start"] for r in rs), "end": max(r["end"] for r in rs),
                     "exons": lib.merge([b for r in rs for b in r["exons"]]), "members": members,
                     "coding": any(r["cds"] for r in rs), "has": {n: any(r["ann"] == n for r in rs) for n in (A, Bn)},
                     "names": sorted({r["name"] for r in rs})})
    loci.sort(key=lambda l: (l["chrom"], l["start"], l["end"]))
    for k, l in enumerate(loci):
        l["id"] = k
        for i in l["members"]:
            rec_locus[(recs[i]["ann"], f"{recs[i]['chrom']}:{recs[i]['start'] + 1}-{recs[i]['end']}")] = k
    # E1 graph edges of each annotation, on joint loci
    E = {}
    unmapped = collections.Counter()
    for name, spec in anns:
        tag = spec.split(":", 1)[1]
        E[name] = {}
        for line in open(tag + ".graph.tsv"):
            x, y, w = line.rstrip("\n").split("\t")
            if x.rsplit(":", 1)[0] not in contigs or y.rsplit(":", 1)[0] not in contigs:
                continue
            u, v = rec_locus.get((name, x)), rec_locus.get((name, y))
            if u is None or v is None:
                unmapped[name] += 1
                continue
            if u != v:
                k = (min(u, v), max(u, v))
                E[name][k] = max(E[name].get(k, 0.0), float(w))
    agreed = set(E[A]) & set(E[Bn])
    disputed = set(E[A]) ^ set(E[Bn])
    # HGNC hard negatives: coding loci sharing a gene group, both annotated in both, no edge in either graph
    groups_of = collections.defaultdict(set)
    for r in csv.DictReader(open(a.hgnc), delimiter="\t"):
        for g in r["gene_group_id"].split("|"):
            if g:
                groups_of[r["symbol"]].add(g)
    grp_loci = collections.defaultdict(set)
    for l in loci:
        if l["coding"] and l["has"][A] and l["has"][Bn]:
            for nm in l["names"]:
                for g in groups_of.get(nm, ()):
                    grp_loci[g].add(l["id"])
    hard_neg = set()
    for ids in grp_loci.values():
        for p in itertools.combinations(sorted(ids), 2):
            if p not in E[A] and p not in E[Bn]:
                hard_neg.add(p)
    involved = sorted({x for p in (agreed | disputed | hard_neg) for x in p})
    X = blast_evidence(a, loci, involved, genome)
    sd = load_sedef(a.sedef, contigs)

    def ev(p):
        u, v = p
        return ((u, v) in X or (v, u) in X), (sd_evidence(sd, loci[u], loci[v]) or sd_evidence(sd, loci[v], loci[u]))
    evd = {p: ev(p) for p in sorted(agreed | disputed | hard_neg)}
    w = lambda p: sum(E[n][p] for n in (A, Bn) if p in E[n]) / sum(1 for n in (A, Bn) if p in E[n])
    strict = {p: w(p) for p in agreed | {p for p in disputed if any(evd[p])}}
    permissive = {p: w(p) for p in agreed | disputed}
    label = {}
    # B13: mcl_port (the Rust-bin shim since 21d6c5c9, r1047) takes STRING node ids and the old build passed int locus
    # ids, so `adjudicated_truth.py build` raised TypeError there. Zero-padded ids sort (byte order, as the bin does)
    # exactly like the ints the numpy mcl_port sorted, so this reproduces the pre-21d6c5c9 clusters.
    width = len(str(max(len(loci) - 1, 0)))
    for tag, graph in (("strict", strict), ("permissive", permissive)):
        sgraph = {(f"{u:0{width}d}", f"{v:0{width}d}"): wt for (u, v), wt in graph.items()}
        for k, c in enumerate(lib.mcl(sgraph, a.inflation, a.prune)):
            for x in c:
                label[(tag, int(x))] = (k, len(c))

    def co(tag, u, v):
        lu, lv = label.get((tag, u)), label.get((tag, v))
        return lu is not None and lv is not None and lu[0] == lv[0] and lu[1] >= 2
    cand = set()
    for tag in ("strict", "permissive"):
        cl = collections.defaultdict(list)
        for (t, x), (k, n) in label.items():
            if t == tag and n >= 2:
                cl[k].append(x)
        for ms in cl.values():
            cand.update(itertools.combinations(sorted(ms), 2))
    status = {p: ("TRUE" if co("strict", *p) and co("permissive", *p) else "UNSCORED") for p in cand}
    with open(f"{a.out}/loci.tsv", "w") as fh:
        fh.write(f"locus\tchrom\tstart\tend\tcoding\thas_{A}\thas_{Bn}\tnames\texons\n")
        for l in loci:
            fh.write(f"L{l['id']}\t{l['chrom']}\t{l['start']}\t{l['end']}\t{int(l['coding'])}\t{int(l['has'][A])}\t"
                     f"{int(l['has'][Bn])}\t{','.join(l['names'])[:500]}\t{','.join(f'{x}-{y}' for x, y in l['exons'])}\n")
    with open(f"{a.out}/edges.tsv", "w") as fh:
        fh.write(f"u\tv\tin_{A}\tin_{Bn}\tblast\tsd\tclass\n")
        for p in sorted(evd):
            cls = "agreed" if p in agreed else ("disputed" if p in disputed else "hard_negative")
            fh.write(f"L{p[0]}\tL{p[1]}\t{int(p in E[A])}\t{int(p in E[Bn])}\t{int(evd[p][0])}\t{int(evd[p][1])}\t{cls}\n")
    with open(f"{a.out}/pairs.tsv", "w") as fh:
        fh.write("u\tv\tstatus\n")
        for p in sorted(status):
            fh.write(f"L{p[0]}\tL{p[1]}\t{status[p]}\n")
    uf2 = lib.UF()
    for p, s in status.items():
        if s == "TRUE":
            uf2.union(p[0], p[1])
    comps = collections.defaultdict(list)
    for x in uf2.p:
        comps[uf2.find(x)].append(x)
    with open(f"{a.out}/clusters.tsv", "w") as fh:
        fh.write("cluster_id\tlocus\tchrom\tstart\tend\n")
        for k, (root, ms) in enumerate(sorted(comps.items())):
            for x in sorted(ms):
                l = loci[x]
                fh.write(f"T{k}\tL{x}\t{l['chrom']}\t{l['start'] + 1}\t{l['end']}\n")
    n_true = sum(1 for s in status.values() if s == "TRUE")
    n_uns = len(status) - n_true
    sens = sum(1 for p in agreed if any(evd[p])) / max(1, len(agreed))
    fpr = sum(1 for p in hard_neg if any(evd[p])) / max(1, len(hard_neg))
    dis_kept = sum(1 for p in disputed if any(evd[p]))
    print(f"records {len(recs)}; joint loci {len(loci)}; graph edges not mapped to a locus {dict(unmapped)}")
    print(f"edges {A} {len(E[A])}, {Bn} {len(E[Bn])}; agreed {len(agreed)}; disputed {len(disputed)} "
          f"(with evidence {dis_kept}); hard negatives {len(hard_neg)}")
    print(f"TRUE pairs {n_true}; UNSCORED {n_uns}; truth clusters {len(comps)} (largest {max(map(len, comps.values())) if comps else 0})")
    ok = sens >= 0.80 and fpr <= 0.20 and n_uns <= 0.50 * (n_true + n_uns)
    print(f"GATE: evidence sensitivity on agreed edges {sens:.3f} (>= 0.80); hard-negative evidence rate {fpr:.3f} "
          f"(<= 0.20); unscored share {n_uns / max(1, n_true + n_uns):.3f} (<= 0.50) -> {'VALID' if ok else 'NOT VALID'}")


# ================================================================ CLI
def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split('\n\n')[0],
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd", required=True)

    p = sub.add_parser("nodes", help="gene-level node tables from RefSeq / CAT / Ensembl (was annotation_nodes.py)",
                       description=cmd_nodes.__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("kind", choices=("refseq", "cat", "ensembl"))
    p.add_argument("gff")
    p.add_argument("contigs", help="comma-separated contig names")
    p.add_argument("out", help="node table path; also writes OUT.names.tsv and OUT.cds.tsv")
    p.set_defaults(func=cmd_nodes)

    p = sub.add_parser("adjudicated", help="AK two-annotation truth build (was adjudicated_truth.py build)",
                       description=cmd_adjudicated.__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    for k in ("--out", "--contigs", "--genome", "--sedef", "--hgnc"):
        p.add_argument(k, required=True)
    p.add_argument("--ann", action="append", required=True, help="NAME=NODES_TSV:TAG_PREFIX (exactly two)")
    p.add_argument("--blast-bin", default="/home/juanfra/miniforge3/envs/blast/bin")
    p.add_argument("--threads", type=int, default=4)
    p.add_argument("--inflation", type=float, default=2.8)
    p.add_argument("--joint-loci", choices=("union", "matched"), default="matched")
    p.add_argument("--prune", type=float, default=1e-9)
    p.set_defaults(func=cmd_adjudicated)

    p = sub.add_parser("protein", help="protein-space families, Addendum AN (was protein_families.py build)",
                       description="One protein per gene (longest CDS from `<NODES_TSV>.cds.tsv`, >= 10 aa); all-vs-all "
                       "blastp -evalue 1e-5 (cached at PREFIX.blastp.tsv); per ordered pair, non-overlapping HSPs greedy "
                       "by bitscore on the longer protein; EDGE iff their union covers >= 0.30 of the longer protein (and "
                       "identity >= --min-ident); weight identity x coverage; MCL I=2.8 prune 1e-9. Writes "
                       "PREFIX[.i<min-ident>].families.tsv and .edges.tsv.")
    for k in ("--nodes", "--genome", "--contigs", "--out"):
        p.add_argument(k, required=True)
    p.add_argument("--min-ident", type=float, default=0.0)
    p.add_argument("--no-pseudogenes", action="store_true", help="AN r1: drop records whose biotype contains 'pseudogene'")
    p.add_argument("--rule", type=int, default=0, help="0 AN, 1 r1 (no pseudogenes), 2 r2 (r1 + V(D)J segments)")
    p.add_argument("--threads", type=int, default=4)
    p.set_defaults(func=cmd_protein)

    p = sub.add_parser("protein-referee", help="the protein referee as a Gene Name / Family ID table (new)",
                       description=protein_referee.__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    for k in ("--gff", "--genome", "--chrom", "--out"):
        p.add_argument(k, required=True)
    p.add_argument("--threads", default="4")
    p.set_defaults(func=cmd_protein_referee)

    a = ap.parse_args(argv)
    a.func(a)


if __name__ == "__main__":
    main()
