#!/usr/bin/env python3
"""PREREG r1820 (docs/PREREG_undecidable_copies_2026-11-19.md): do undecidable reads depict the copies they cannot be
placed on? Four subcommands, one workspace directory OUT:

  undecidable.py inputs  --truth READS.tsv --tie-split T.tsv --bam SEL.bam --out OUT
      OUT/P.bam      primary records only
      OUT/PU.bam     P + the U* records (a SEVERAL COPIES read's records on its truth origin loci, other than its primary's)
      OUT/PUms.bam   P + the Ums records (every non-primary record at the best ms of a read whose two best places tie exactly
                     on ms; no truth used)
      OUT/copies.tsv landing copies (merged truth loci holding a U* record), with U reads / U* records / Ums overlap
      OUT/inputs.tsv counts per read kind
  undecidable.py classes --truth READS.tsv --tie-split T.tsv --bam SEL.bam --annot ANNOT.gff3[.gz] --out OUT
                         [--own 2] [--shadow 5] [--read-family RF.tsv]
      OUT/classes.tsv  each landing copy: gene, biotype group, Dd origin reads, Dd power reads (any score / within 10 ms),
                       class own / own_1 / shadow / unjudgeable;  OUT/q1.tsv, OUT/q2.tsv the descriptive tables
  undecidable.py proj    --fasta GENOME.fa --p-gtf P.gtf --truth READS.tsv --tie-split T.tsv --bam SEL.bam --out OUT
                         [--threads 8] [--min-cov 0.9] [--min-id 0.85] [--pad 5000]
      OUT/PROJ.gtf  P's multi-exon chains at each U read's primary place, projected (minimap2 splice:hq) onto that read's
                    landing copies; kept at cov >= 0.9 and identity >= 0.85, as row 1519
  undecidable.py score   --out OUT --annot ANNOT.gff3[.gz] --arm P=GTF --arm PU=GTF --arm PUms=GTF --arm D=GTF
                         --arm PROJ=OUT/PROJ.gtf [--arm NAME=GTF ...]
      OUT/depicted.tsv per copy x arm; OUT/bars.tsv C1, C1b, C2, C3 and the cost lines, ms first (the U set is the truth's)

Conventions: intervals 0-based half-open; an intron chain is the tuple of (exon end, next exon start); a transcript is
multi-exon when it has >= 1 intron. Wilson 95% intervals (records.wilson). No path is fixed here."""
import argparse, collections, csv, gzip, os, re, subprocess, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import records as R

GROUPS = {
    "protein_coding": "protein_coding",
    "unprocessed_pseudogene": "unprocessed_pseudogene", "transcribed_unprocessed_pseudogene": "unprocessed_pseudogene",
    "unitary_pseudogene": "unprocessed_pseudogene", "transcribed_unitary_pseudogene": "unprocessed_pseudogene",
    "polymorphic_pseudogene": "unprocessed_pseudogene", "translated_unprocessed_pseudogene": "unprocessed_pseudogene",
    "processed_pseudogene": "processed_pseudogene", "transcribed_processed_pseudogene": "processed_pseudogene",
    "translated_processed_pseudogene": "processed_pseudogene",
    "pseudogene": "pseudogene_unsplit",          # NCBI and GENCODE's bare 'pseudogene': type not given
    "lncRNA": "lncRNA",
}
CLASSES = ("own", "own_1", "shadow", "unjudgeable")
CLAIM = ("shadow", "unjudgeable")                 # the copies the claim is about

def group(biotype):
    return GROUPS.get(biotype, "none" if biotype is None else "other")

# ---------------------------------------------------------------- read side

def u_reads(truth):
    """reads of kind several (the truth's SEVERAL COPIES ties) = U"""
    return {q for q, t in truth.items() if t.kind == "several"}

def primary_of(rs):
    return next((r for r in rs if r.primary), rs[0])

def u_star(rs, t):
    """ordinals of a U read's records on its origin loci, other than the primary's locus (and not the primary itself)"""
    p = primary_of(rs)
    pi = R.locus_index(p, t)
    out = []
    for i, r in enumerate(rs):
        if r.primary:
            continue
        li = R.locus_index(r, t)
        if li is not None and li in t.origin and li != pi:
            out.append(i)
    return out

def u_ms(rs):
    """ordinals of the non-primary records at the best ms of a read whose two best places tie exactly on ms"""
    if R.lead(rs, "ms") != 0:
        return []
    best = max(r.ms for r in rs)
    return [i for i, r in enumerate(rs) if not r.primary and r.ms == best]

def merge_spans(spans):
    """[(chrom, strand, s, e, payload)] -> merged copies [(chrom, strand, s, e, [payloads])], sorted"""
    out = []
    for c, st, s, e, x in sorted(spans, key=lambda v: (v[0], v[1], v[2], v[3])):
        if out and out[-1][0] == c and out[-1][1] == st and s < out[-1][3]:
            out[-1][3] = max(out[-1][3], e); out[-1][4].append(x)
        else:
            out.append([c, st, s, e, [x]])
    return [tuple(v) for v in out]

def landing_copies(truth, recs):
    """-> list of copies (id, chrom, strand, start, end, U reads, U* records); the span of a copy is the union of the truth
    loci (not the records) that hold a U* record"""
    spans = []
    for q in u_reads(truth):
        rs, t = recs.get(q), truth[q]
        if not rs:
            continue
        for i in u_star(rs, t):
            c, s, e, st, _ = t.loci[R.locus_index(rs[i], t)]
            spans.append((c, st, s, e, q))
    out = []
    for k, (c, st, s, e, qs) in enumerate(merge_spans(spans)):
        out.append(dict(id=f"c{k + 1}", chrom=c, strand=st, start=s, end=e, u_reads=len(set(qs)), u_records=len(qs)))
    return out

class CopyIndex:
    """overlap lookup: (chrom, strand, start, end) -> copy ids"""
    def __init__(self, copies, bin_size=100000):
        self.b, self.n = collections.defaultdict(list), bin_size
        for c in copies:
            for k in range(c["start"] // bin_size, (c["end"] - 1) // bin_size + 1):
                self.b[(c["chrom"], c["strand"], k)].append(c)
    def hits(self, chrom, strand, s, e):
        seen = {}
        for k in range(s // self.n, max(s, e - 1) // self.n + 1):
            for c in self.b.get((chrom, strand, k), ()):
                if c["start"] < e and s < c["end"]:
                    seen[c["id"]] = c
        return list(seen.values())

def write_bams(bam, plans, out):
    import pysam
    tmp = os.environ.get("TMPDIR", out)
    with pysam.AlignmentFile(bam, check_sq=False) as fh:
        w = {x: pysam.AlignmentFile(f"{tmp}/undecidable.{os.getpid()}.{x}.bam", "wb", template=fh) for x in ("P", "PU", "PUms")}
        seen = collections.Counter()
        for r in fh.fetch(until_eof=True):
            if r.is_unmapped or r.is_supplementary:
                continue
            i = seen[r.query_name]; seen[r.query_name] += 1
            keep = plans.get(r.query_name)
            if keep is None:
                continue
            for x in ("P", "PU", "PUms"):
                if i in keep[x]:
                    w[x].write(r)
        for x, f in w.items():
            f.close()
            pysam.sort("-o", f"{out}/{x}.bam", f.filename.decode()); pysam.index(f"{out}/{x}.bam"); os.remove(f.filename)

def cmd_inputs(a):
    os.makedirs(a.out, exist_ok=True)
    truth, recs = R.load_truth(a.truth, a.tie_split), R.load_records(a.bam)
    kinds, plans, n = collections.Counter(), {}, collections.Counter()
    for q, rs in recs.items():
        t = truth.get(q)
        kinds[t.kind if t else "no_truth"] += 1
        prim = [i for i, r in enumerate(rs) if r.primary][:1] or [0]
        us = u_star(rs, t) if t is not None and t.kind == "several" else []
        um = u_ms(rs)
        plans[q] = dict(P=set(prim), PU=set(prim) | set(us), PUms=set(prim) | set(um))
        n["U*_records"] += len(us); n["Ums_records"] += len(um); n["both_records"] += len(set(us) & set(um))
        n["U_reads_with_U*"] += bool(us); n["reads_with_Ums"] += bool(um)
    copies = landing_copies(truth, recs)
    write_bams(a.bam, plans, a.out)
    with open(f"{a.out}/copies.tsv", "w") as fh:
        fh.write("copy\tchrom\tstrand\tstart\tend\tu_reads\tu_records\n")
        for c in copies:
            fh.write("\t".join(str(c[k]) for k in ("id", "chrom", "strand", "start", "end", "u_reads", "u_records")) + "\n")
    with open(f"{a.out}/inputs.tsv", "w") as fh:
        fh.writelines(f"kind_{k}\t{v}\n" for k, v in sorted(kinds.items()))
        fh.writelines(f"{k}\t{v}\n" for k, v in n.items())
        fh.write(f"landing_copies\t{len(copies)}\n")
    print(open(f"{a.out}/inputs.tsv").read(), end="")

def read_copies(out):
    with open(f"{out}/copies.tsv") as fh:
        return [dict(r, start=int(r["start"]), end=int(r["end"]), id=r["copy"]) for r in csv.DictReader(fh, delimiter="\t")]

# ---------------------------------------------------------------- annotation

def opn(p):
    return gzip.open(p, "rt") if p.endswith(".gz") else open(p)

def load_annot(path):
    """GFF3 (CAT / Liftoff slim, NCBI) or GTF -> (transcripts, genes)
    transcripts: {tid: dict(chrom, strand, exons, gene)}; genes: {gid: dict(chrom, strand, start, end, biotype, name)}.
    An exon whose parent is a gene-level feature (NCBI pseudogenes) makes that gene its own transcript."""
    genes, parent, ex = {}, {}, collections.defaultdict(list)
    for l in opn(path):
        if l[0] == "#":
            continue
        f = l.rstrip("\n").split("\t")
        if len(f) < 9:
            continue
        if "=" in f[8] and '"' not in f[8].split(";")[0]:
            at = dict(x.split("=", 1) for x in f[8].split(";") if "=" in x)
            fid, par = at.get("ID"), at.get("Parent", "").split(",")[0]
        else:
            at = dict(re.findall(r'(\S+) "([^"]*)"', f[8]))
            fid = at.get("transcript_id") if f[2] in ("transcript", "exon") else at.get("gene_id")
            par = at.get("gene_id") if f[2] == "transcript" else at.get("transcript_id", "")
        if "gene_biotype" in at or f[2] in ("gene", "pseudogene"):
            genes[fid] = dict(chrom=f[0], strand=f[6], start=int(f[3]) - 1, end=int(f[4]),
                              biotype=at.get("gene_biotype") or at.get("gene_type"),
                              name=at.get("gene_name") or at.get("Name") or at.get("gene") or fid)
        elif f[2] == "exon":
            ex[par].append((f[0], f[6], int(f[3]) - 1, int(f[4])))
        elif fid:
            parent[fid] = par
    tx = {}
    for t, es in ex.items():
        g = parent.get(t, t) if t not in genes else t
        es.sort(key=lambda v: v[2])
        tx[t] = dict(chrom=es[0][0], strand=es[0][1], exons=[(s, e) for _, _, s, e in es], gene=g)
    return tx, genes

def chain(exons):
    return tuple((a[1], b[0]) for a, b in zip(exons, exons[1:]))

def contiguous_part(small, big):
    """small is a contiguous run of >= 2 introns of big (or equal to it)"""
    k = len(small)
    return k >= 2 and any(big[i:i + k] == small for i in range(len(big) - k + 1))

def copy_annotation(copies, tx, genes):
    """-> {copy id: dict(gene, biotype, group, chains={annotated multi-exon chains}, exonic=bp of the gene's exon union)}"""
    idx = CopyIndex(copies)
    best = {}
    for gid, g in genes.items():
        for c in idx.hits(g["chrom"], g["strand"], g["start"], g["end"]):
            ov = min(c["end"], g["end"]) - max(c["start"], g["start"])
            if ov > best.get(c["id"], (0, None))[0]:
                best[c["id"]] = (ov, gid)
    chains, union = collections.defaultdict(set), collections.defaultdict(list)
    for t in tx.values():
        s, e = t["exons"][0][0], t["exons"][-1][1]
        for c in idx.hits(t["chrom"], t["strand"], s, e):
            ch = chain(t["exons"])
            if ch:
                chains[c["id"]].add(ch)
            if best.get(c["id"], (0, None))[1] == t["gene"]:
                union[c["id"]] += t["exons"]
    out = {}
    for c in copies:
        gid = best.get(c["id"], (0, None))[1]
        g = genes.get(gid, {})
        out[c["id"]] = dict(gene=g.get("name", ""), biotype=g.get("biotype"), group=group(g.get("biotype")) if gid else "none",
                            chains=chains[c["id"]], exonic=sum(e - s for s, e in merge_iv(union[c["id"]])), exons=merge_iv(union[c["id"]]))
    return out

def merge_iv(iv):
    out = []
    for s, e in sorted(iv):
        if out and s <= out[-1][1]:
            out[-1][1] = max(out[-1][1], e)
        else:
            out.append([s, e])
    return [tuple(v) for v in out]

# ---------------------------------------------------------------- classes (Q1, Q2)

def classify(copies, truth, recs, own=2, shadow=5):
    """-> {copy id: dict(dd_origin, dd_power_any, dd_power_10, cls)}"""
    idx = CopyIndex(copies)
    st = {c["id"]: dict(dd_origin=0, dd_power_any=0, dd_power_10=0) for c in copies}
    for q, t in truth.items():
        if t.kind != "one":
            continue
        o = t.loci[next(iter(t.origin))]
        home = {c["id"] for c in idx.hits(o[0], o[3], o[1], o[2])}
        for cid in home:
            st[cid]["dd_origin"] += 1
        rs = recs.get(q, [])
        if not rs:
            continue
        best = max(r.ms for r in rs)
        anyc, near = set(), set()
        for r in rs:
            for c in idx.hits(r.chrom, r.strand, r.start, r.end):
                if c["id"] not in home:
                    anyc.add(c["id"])
                    if best - r.ms <= 10:
                        near.add(c["id"])
        for cid in anyc:
            st[cid]["dd_power_any"] += 1
        for cid in near:
            st[cid]["dd_power_10"] += 1
    for v in st.values():
        v["cls"] = ("own" if v["dd_origin"] >= own else "own_1" if v["dd_origin"] == 1 else
                    "shadow" if v["dd_power_any"] >= shadow else "unjudgeable")
        v["cls_10"] = ("own" if v["dd_origin"] >= own else "own_1" if v["dd_origin"] == 1 else
                       "shadow" if v["dd_power_10"] >= shadow else "unjudgeable")
    return st

def read_family(path):
    if not path:
        return {}
    with open(path) as fh:
        return dict(l.rstrip("\n").split("\t")[:2] for l in fh if "\t" in l)

def cmd_classes(a):
    truth, recs = R.load_truth(a.truth, a.tie_split), R.load_records(a.bam)
    copies = read_copies(a.out)
    tx, genes = load_annot(a.annot)
    ann = copy_annotation(copies, tx, genes)
    st = classify(copies, truth, recs, a.own, a.shadow)
    fam = read_family(a.read_family)
    cfam = collections.defaultdict(collections.Counter)
    if fam:
        idx = CopyIndex(copies)
        for q in u_reads(truth):
            rs = recs.get(q)
            if not rs:
                continue
            for i in u_star(rs, truth[q]):
                for c in idx.hits(rs[i].chrom, rs[i].strand, rs[i].start, rs[i].end):
                    cfam[c["id"]][fam.get(q, "-")] += 1
    with open(f"{a.out}/classes.tsv", "w") as fh:
        fh.write("copy\tchrom\tstrand\tstart\tend\tu_reads\tu_records\tgene\tbiotype\tgroup\tannotated_chains\tfamily\t"
                 "dd_origin\tdd_power_any\tdd_power_10\tclass\tclass_10\n")
        for c in copies:
            s, n = st[c["id"]], ann[c["id"]]
            f = cfam[c["id"]].most_common(1)[0][0] if cfam[c["id"]] else "-"
            fh.write("\t".join(map(str, (c["id"], c["chrom"], c["strand"], c["start"], c["end"], c["u_reads"], c["u_records"],
                                          n["gene"], n["biotype"] or "-", n["group"], len(n["chains"]), f, s["dd_origin"],
                                          s["dd_power_any"], s["dd_power_10"], s["cls"], s["cls_10"]))) + "\n")
    q1, q2 = collections.Counter(), collections.Counter()
    for c in copies:
        g = ann[c["id"]]["group"]
        q1[(g, "copies")] += 1; q1[(g, "u_reads")] += int(c["u_reads"]); q1[(g, "u_records")] += int(c["u_records"])
        q2[(st[c["id"]]["cls"], g)] += 1
    gs = sorted({k[0] for k in q1})
    with open(f"{a.out}/q1.tsv", "w") as fh:
        fh.write("group\tcopies\tu_reads\tu_records\n")
        fh.writelines(f"{g}\t{q1[(g, 'copies')]}\t{q1[(g, 'u_reads')]}\t{q1[(g, 'u_records')]}\n" for g in gs)
    with open(f"{a.out}/q2.tsv", "w") as fh:
        fh.write("class\t" + "\t".join(gs) + "\ttotal\n")
        for k in CLASSES:
            fh.write(k + "\t" + "\t".join(str(q2[(k, g)]) for g in gs) + f"\t{sum(q2[(k, g)] for g in gs)}\n")
        fh.write(f"# shadow under the within-10-ms power rule (printed beside): {sum(v['cls_10'] == 'shadow' for v in st.values())}\n")
    print(open(f"{a.out}/q1.tsv").read() + open(f"{a.out}/q2.tsv").read(), end="")

# ---------------------------------------------------------------- PROJ

def read_gtf_tx(path):
    """an assembler's GTF -> {tid: dict(chrom, strand, exons)}"""
    ex = collections.defaultdict(list)
    for l in opn(path):
        f = l.rstrip("\n").split("\t")
        if len(f) < 9 or f[2] != "exon":
            continue
        t = re.search(r'transcript_id "([^"]+)"', f[8]).group(1)
        ex[t].append((f[0], f[6], int(f[3]) - 1, int(f[4])))
    return {t: dict(chrom=v[0][0], strand=v[0][1], exons=sorted((s, e) for _, _, s, e in v)) for t, v in ex.items()}

def revcomp(s):
    return s.translate(str.maketrans("ACGTNacgtn", "TGCANtgcan"))[::-1]

def sam_hits(sam_lines):
    """minimap2 SAM (with --eqx) -> [(query, target, reverse, target_start, exons, cov, ident)]"""
    out = []
    for l in sam_lines:
        if l[0] == "@":
            continue
        f = l.split("\t")
        flag = int(f[1])
        if flag & 4 or flag & 2048:
            continue
        pos, cig = int(f[3]) - 1, re.findall(r"(\d+)([MIDNSHP=X])", f[5])
        qlen = sum(int(n) for n, o in cig if o in "MI=XSH")
        aq, match, cols, exons, cur = 0, 0, 0, [], pos
        s0 = pos
        for n, o in cig:
            n = int(n)
            if o in "=XM":
                aq += n; cols += n; match += n if o in "=M" else 0; cur += n
            elif o == "I":
                aq += n; cols += n
            elif o == "D":
                cols += n; cur += n
            elif o == "N":
                exons.append((s0, cur)); cur += n; s0 = cur
        exons.append((s0, cur))
        out.append((f[0], f[2], bool(flag & 16), exons, aq / qlen if qlen else 0, match / cols if cols else 0))
    return out

def cmd_proj(a):
    import pysam
    truth, recs = R.load_truth(a.truth, a.tie_split), R.load_records(a.bam)
    copies = read_copies(a.out)
    cidx = CopyIndex(copies)
    ptx = {t: v for t, v in read_gtf_tx(a.p_gtf).items() if len(v["exons"]) > 1}
    tb = collections.defaultdict(list)
    for t, v in ptx.items():
        tb[(v["chrom"], v["strand"], v["exons"][0][0] // 100000)].append(t)
    def p_at(chrom, strand, s, e):
        hits = set()
        for k in range(s // 100000 - 1, (e - 1) // 100000 + 1):
            for t in tb.get((chrom, strand, k), ()):
                v = ptx[t]
                if v["exons"][0][0] < e and s < v["exons"][-1][1]:
                    hits.add(t)
        return hits
    pairs = collections.defaultdict(set)              # template -> copies
    for q in u_reads(truth):
        rs, t = recs.get(q), truth[q]
        if not rs:
            continue
        p = primary_of(rs)
        lands = {c["id"] for i in u_star(rs, t) for c in cidx.hits(rs[i].chrom, rs[i].strand, rs[i].start, rs[i].end)}
        for tid in p_at(p.chrom, p.strand, p.start, p.end):
            pairs[tid] |= lands
    n, nc = project(pairs, ptx, copies, a.fasta, f"{a.out}/PROJ.gtf", a.threads, a.min_cov, a.min_id, a.pad)
    print(f"proj: {len(pairs)} templates, {nc} copies, {n} projected chains kept")

def project(pairs, ptx, copies, fasta, out_gtf, threads=8, min_cov=0.9, min_id=0.85, pad=5000):
    """templates (ptx[tid]) projected onto their paired copies (pairs: tid -> copy ids); -> (chains kept, copies used)"""
    import pysam
    fa = pysam.FastaFile(fasta)
    tmp = os.environ.get("TMPDIR", os.path.dirname(os.path.abspath(out_gtf)))
    qf, tf = f"{tmp}/proj.{os.getpid()}.q.fa", f"{tmp}/proj.{os.getpid()}.t.fa"
    byid = {c["id"]: c for c in copies}
    used = sorted({c for v in pairs.values() for c in v})
    with open(qf, "w") as fh:
        for tid in sorted(pairs):
            v = ptx[tid]
            s = "".join(fa.fetch(v["chrom"], s, e) for s, e in v["exons"])
            fh.write(f">{tid}\n{revcomp(s) if v['strand'] == '-' else s}\n")
    with open(tf, "w") as fh:
        for cid in used:
            c = byid[cid]
            s0 = max(0, c["start"] - pad)
            fh.write(f">{cid}|{c['chrom']}|{s0}\n{fa.fetch(c['chrom'], s0, c['end'] + pad)}\n")
    cmd = ["minimap2", "-ax", "splice:hq", "-uf", "--eqx", "-N", "50", "-p", "0.5", "--secondary=yes", "-t", str(threads), tf, qf]
    sam = subprocess.run(cmd, check=True, capture_output=True, text=True).stdout.splitlines()
    n, kept = 0, set()
    with open(out_gtf, "w") as fh:
        for qn, tn, rev, exons, cov, ident in sam_hits(sam):
            cid, chrom, s0 = tn.split("|"); s0 = int(s0)
            if cid not in pairs.get(qn, ()) or cov < min_cov or ident < min_id or len(exons) < 2:
                continue
            strand = "-" if rev else "+"             # the template is in transcript orientation, the target forward
            if strand != byid[cid]["strand"]:
                continue
            g = [(s + s0, e + s0) for s, e in exons]
            key = (cid, chain(g))
            if key in kept:
                continue
            kept.add(key); n += 1
            tid = f"proj.{n}"
            attr = f'gene_id "{cid}"; transcript_id "{tid}"; template "{qn}"; cov "{cov:.3f}"; ident "{ident:.3f}";'
            fh.write(f"{chrom}\tproj\ttranscript\t{g[0][0] + 1}\t{g[-1][1]}\t.\t{strand}\t.\t{attr}\n")
            for s, e in g:
                fh.write(f"{chrom}\tproj\texon\t{s + 1}\t{e}\t.\t{strand}\t.\t{attr}\n")
    os.remove(qf); os.remove(tf)
    return n, len(used)

# ---------------------------------------------------------------- score (C1, C1b, C2, C3)

def arm_chains(gtf, copies):
    """-> {copy id: (set of multi-exon chains, list of single-exon (s, e))}"""
    idx = CopyIndex(copies)
    out = collections.defaultdict(lambda: (set(), []))
    for v in read_gtf_tx(gtf).values():
        s, e = v["exons"][0][0], v["exons"][-1][1]
        for c in idx.hits(v["chrom"], v["strand"], s, e):
            if len(v["exons"]) > 1:
                out[c["id"]][0].add(chain(v["exons"]))
            else:
                out[c["id"]][1].append((s, e))
    return out

def depicted(chains, ann):
    """(exact, part) for a set of emitted chains against a copy's annotated chains"""
    exact = bool(chains & ann)
    part = exact or any(contiguous_part(c, b) for c in chains for b in ann)
    return exact, part

def matches(c, ann):
    """'exact', 'part' or None"""
    if c in ann:
        return "exact"
    return "part" if any(contiguous_part(c, b) for b in ann) else None

def single_cover(singles, exons, total):
    if not total:
        return False
    for s, e in singles:
        cov = sum(max(0, min(e, b) - max(s, a)) for a, b in exons)
        if cov >= 0.8 * total:
            return True
    return False

def fmt(k, n):
    lo, hi = R.wilson(k, n)
    return f"{k}/{n} = {k / n if n else 0:.3f} [{lo:.3f}, {hi:.3f}]"

def bars(cls, ann, A):
    """the readings; cls {copy: class}, ann {copy: annotation}, A {arm: {copy: (chains, singles)}} -> list of lines"""
    claim = [c for c, k in cls.items() if k in CLAIM]
    dep = {x: {c for c in claim if depicted(A[x][c][0], ann[c]["chains"])[0]} for x in A}
    added_dep = dep["PU"] - dep["P"]
    added = [(c, ch) for c in claim for ch in A["PU"][c][0] - A["P"][c][0]]
    k = sum(matches(ch, ann[c]["chains"]) == "exact" for c, ch in added)
    lo, hi = R.wilson(k, len(added))
    c1 = "HOLDS" if len(added_dep) >= 20 and lo >= 0.5 else "FAILS" if len(added_dep) < 20 or hi < 0.5 else "BETWEEN"
    L = [f"claim copies (shadow + unjudgeable)\t{len(claim)}",
         f"C1 added depicted copies (PU not P)\t{len(added_dep)}",
         f"C1 depicted share of chains PU adds\t{fmt(k, len(added))}",
         f"C1\t{c1}"]
    notproj = added_dep - dep["PROJ"] if "PROJ" in dep else set()
    kp = sum(matches(ch, ann[c]["chains"]) == "exact" for c in added_dep for ch in A["PROJ"][c][0]) if "PROJ" in A else 0
    np_ = sum(len(A["PROJ"][c][0]) for c in added_dep) if "PROJ" in A else 0
    L += [f"C1b added depicted copies not depicted by PROJ\t{len(notproj)}",
          f"C1b PROJ depicted share on the same copies\t{fmt(kp, np_)}",
          f"C1b\t{'HOLDS' if c1 == 'HOLDS' and len(notproj) >= 10 else 'NOT MET'}"]
    if "D" in dep:
        nd = len(added_dep - dep["D"])
        lo2, _ = R.wilson(nd, len(added_dep))
        L += [f"C2 share of C1's copies D does not depict\t{fmt(nd, len(added_dep))}",
              f"C2\t{'NEW OVER DEFAULTS' if lo2 >= 0.3 else 'NOT NEW'}"]
    if "PUms" in dep:
        reach = len(added_dep & dep["PUms"])
        addm = [(c, ch) for c in added_dep for ch in A["PUms"][c][0] - A["P"][c][0]]
        addu = [(c, ch) for c in added_dep for ch in A["PU"][c][0] - A["P"][c][0]]
        km = sum(matches(ch, ann[c]["chains"]) == "exact" for c, ch in addm)
        ku = sum(matches(ch, ann[c]["chains"]) == "exact" for c, ch in addu)
        sm, su = (km / len(addm) if addm else 0), (ku / len(addu) if addu else 0)
        ok = c1 == "HOLDS" and reach >= 0.5 * len(added_dep) and abs(sm - su) <= 0.1
        L += [f"C3 PUms reach of C1's copies\t{reach}/{len(added_dep)}",
              f"C3 depicted share PUms / PU on those copies\t{fmt(km, len(addm))} / {fmt(ku, len(addu))}",
              f"C3\t{'HOLDS' if ok else 'NOT MET' if c1 == 'HOLDS' else 'not read (C1 does not hold)'}"]
    for x in [x for x in ("PU", "PUms") if x in A]:
        for k2 in ("shadow", "own"):
            cs = [c for c, v in cls.items() if v == k2]
            un = sum(matches(ch, ann[c]["chains"]) is None for c in cs for ch in A[x][c][0] - A["P"][c][0])
            L.append(f"cost {x} unmatched chains added on {k2} copies\t{un}")
    return L

def cmd_score(a):
    copies = read_copies(a.out)
    with open(f"{a.out}/classes.tsv") as fh:
        cls = {r["copy"]: r["class"] for r in csv.DictReader(fh, delimiter="\t")}
    tx, genes = load_annot(a.annot)
    ann = copy_annotation(copies, tx, genes)
    A = {}
    for spec in a.arm:
        name, gtf = spec.split("=", 1)
        A[name] = arm_chains(gtf, copies)
    for need in ("P", "PU"):
        if need not in A:
            sys.exit(f"undecidable score: --arm {need}=GTF is required")
    with open(f"{a.out}/depicted.tsv", "w") as fh:
        fh.write("copy\tclass\tgroup\tannotated_chains\t" + "\t".join(f"{x}_chains\t{x}_exact\t{x}_part\t{x}_single" for x in A) + "\n")
        for c in copies:
            cid, n = c["id"], ann[c["id"]]
            row = [cid, cls[cid], n["group"], len(n["chains"])]
            for x in A:
                ch, sg = A[x][cid]
                e, p = depicted(ch, n["chains"])
                row += [len(ch), int(e), int(p), int(n["group"] == "processed_pseudogene" and single_cover(sg, n["exons"], n["exonic"]))]
            fh.write("\t".join(map(str, row)) + "\n")
    L = bars(cls, ann, A)
    with open(f"{a.out}/bars.tsv", "w") as fh:
        fh.write("\n".join(L) + "\n")
    print("\n".join(L))

# ---------------------------------------------------------------- arrays (descriptive side table, no truth, no bar)

STEM = re.compile(r"^([A-Za-z][A-Za-z0-9-]*[A-Za-z])(\d+[A-Z]?\d*P?)$")
ARRAY_GROUPS = ("protein_coding", "unprocessed_pseudogene", "processed_pseudogene", "pseudogene_unsplit", "lncRNA")

def stem(name):
    """gene name -> family stem (TSPY10, TSPY6P -> TSPY; USP17L12 -> USP17L; RBMY1A1 -> RBMY); None for accession-like names"""
    if not name or "." in name or name.startswith(("LOFF_", "CHM13_", "ENSG", "LOC")) or "-" in name[:3]:
        return None
    m = STEM.match(name)
    return m.group(1).rstrip("-") if m and len(m.group(1).rstrip("-")) >= 3 else None

def find_arrays(genes, min_units=5, gap=200000):
    """clusters of >= min_units annotated genes sharing a stem on one contig, neighbours <= gap apart (either strand)
    -> [(array name, chrom, start, end, [gene ids])]"""
    by = collections.defaultdict(list)
    for gid, g in genes.items():
        st = stem(g["name"])
        if st and group(g["biotype"]) in ARRAY_GROUPS:
            by[(st, g["chrom"])].append((g["start"], g["end"], gid))
    out = []
    for (st, chrom), v in sorted(by.items()):
        v.sort(); cur = [v[0]]
        for x in v[1:] + [None]:
            if x is not None and x[0] - max(e for _, e, _ in cur) <= gap:
                cur.append(x); continue
            if len(cur) >= min_units:
                out.append((f"{st}@{chrom}:{cur[0][0]}", chrom, cur[0][0], max(e for _, e, _ in cur), [g for _, _, g in cur]))
            if x is not None:
                cur = [x]
    return out

def cmd_arrays(a):
    """slice the BAM at the arrays and write the three arm inputs; units.tsv lists every unit"""
    import pysam
    os.makedirs(a.out, exist_ok=True)
    tx, genes = load_annot(a.annot)
    arrays = find_arrays(genes, a.min_units)
    with open(f"{a.out}/units.tsv", "w") as fh:
        fh.write("array\tunit\tgene\tchrom\tstrand\tstart\tend\tbiotype\tgroup\n")
        for name, chrom, s, e, gids in arrays:
            for gid in gids:
                g = genes[gid]
                fh.write(f"{name}\t{gid}\t{g['name']}\t{chrom}\t{g['strand']}\t{g['start']}\t{g['end']}\t{g['biotype'] or '-'}\t{group(g['biotype'])}\n")
    recs = collections.defaultdict(list)
    with pysam.AlignmentFile(a.bam) as fh:
        tmp = os.environ.get("TMPDIR", a.out)
        allw = pysam.AlignmentFile(f"{tmp}/arrays.{os.getpid()}.all.bam", "wb", template=fh)
        seen = set()
        for name, chrom, s, e, _ in arrays:
            for r in fh.fetch(chrom, max(0, s - a.pad), e + a.pad):
                if r.is_unmapped or r.is_supplementary:
                    continue
                key = (r.query_name, r.reference_name, r.reference_start, r.flag, r.cigarstring)
                if key in seen:
                    continue
                seen.add(key); allw.write(r)
                recs[r.query_name].append(R.Rec(r.reference_name, r.reference_start, r.reference_end, "-" if r.is_reverse else "+",
                                                (), r.get_tag("ms"), r.get_tag("AS"), not r.is_secondary))
        allw.close()
    pysam.sort("-o", f"{a.out}/ALL.bam", allw.filename.decode()); pysam.index(f"{a.out}/ALL.bam"); os.remove(allw.filename)
    plans = {}
    for q, rs in recs.items():
        prim = [i for i, r in enumerate(rs) if r.primary][:1]
        plans[q] = dict(P=set(prim), PU=set(), PUms=set(prim) | set(u_ms(rs)))
    # the slice holds a read's records at the arrays only: its ties are judged among those (documented in the PREREG)
    with pysam.AlignmentFile(f"{a.out}/ALL.bam") as fh:
        w = {x: pysam.AlignmentFile(f"{a.out}/{x}.unsorted.bam", "wb", template=fh) for x in ("P", "PUms")}
        idx = collections.Counter()
        order = {}
        for q, rs in recs.items():
            order[q] = rs
        for r in fh.fetch(until_eof=True):
            q = r.query_name
            rs = order[q]
            key = (r.reference_name, r.reference_start, "-" if r.is_reverse else "+", not r.is_secondary)
            i = next(j for j, x in enumerate(rs) if (x.chrom, x.start, x.strand, x.primary) == key)
            for x in ("P", "PUms"):
                if i in plans[q][x]:
                    w[x].write(r)
        for x, f in w.items():
            f.close()
            pysam.sort("-o", f"{a.out}/{x}.bam", f"{a.out}/{x}.unsorted.bam"); pysam.index(f"{a.out}/{x}.bam"); os.remove(f"{a.out}/{x}.unsorted.bam")
    n = collections.Counter(len(rs) > 1 for rs in recs.values())
    print(f"arrays: {len(arrays)} arrays, {sum(len(x[4]) for x in arrays)} units; reads {len(recs)} "
          f"({n[True]} with several records in the slice); PUms adds {sum(len(p['PUms']) - len(p['P']) for p in plans.values())} records")

def read_units(out):
    with open(f"{out}/units.tsv") as fh:
        return [dict(r, start=int(r["start"]), end=int(r["end"]), id=r["unit"]) for r in csv.DictReader(fh, delimiter="\t")]

def cmd_arrays_proj(a):
    units = read_units(a.out)
    ptx = {t: v for t, v in read_gtf_tx(a.p_gtf).items() if len(v["exons"]) > 1}
    uidx = CopyIndex(units)
    arr = collections.defaultdict(list)
    for u in units:
        arr[u["array"]].append(u["id"])
    pairs = collections.defaultdict(set)
    for tid, v in ptx.items():
        for u in uidx.hits(v["chrom"], v["strand"], v["exons"][0][0], v["exons"][-1][1]):
            pairs[tid] |= set(arr[u["array"]]) - {u["id"]}
    n, nc = project(pairs, ptx, units, a.fasta, f"{a.out}/PROJ.gtf", a.threads)
    print(f"arrays proj: {len(pairs)} templates, {nc} units, {n} projected chains kept")

def cmd_arrays_score(a):
    units = read_units(a.out)
    tx, genes = load_annot(a.annot)
    # a unit's annotated chains: its own gene's transcripts only (units of an array may overlap their neighbours' spans)
    own = collections.defaultdict(set)
    for t in tx.values():
        ch = chain(t["exons"])
        if ch:
            own[t["gene"]].add(ch)
    A = {spec.split("=", 1)[0]: arm_chains(spec.split("=", 1)[1], units) for spec in a.arm}
    tot = collections.defaultdict(collections.Counter)
    with open(f"{a.out}/arrays_depicted.tsv", "w") as fh:
        fh.write("array\tunit\tgene\tgroup\tannotated_chains\t" + "\t".join(f"{x}_chains\t{x}_exact\t{x}_part\t{x}_unmatched" for x in A) + "\n")
        for u in units:
            ann = own[u["id"]]
            row = [u["array"], u["id"], u["gene"], u["group"], len(ann)]
            for x in A:
                ch = A[x][u["id"]][0]
                e, p = depicted(ch, ann)
                un = sum(matches(c, ann) is None for c in ch)
                row += [len(ch), int(e), int(p), un]
                if ann:
                    k = (u["array"].split("@")[0], u["group"])
                    tot[k][x + "_exact"] += e; tot[k][x + "_unmatched"] += un
            if ann:
                tot[(u["array"].split("@")[0], u["group"])]["units"] += 1
            fh.write("\t".join(map(str, row)) + "\n")
    with open(f"{a.out}/arrays_summary.tsv", "w") as fh:
        fh.write("family\tgroup\tunits_with_annotated_chains\t" + "\t".join(f"{x}_depicted\t{x}_unmatched" for x in A) + "\n")
        for k in sorted(tot):
            fh.write(f"{k[0]}\t{k[1]}\t{tot[k]['units']}\t" + "\t".join(f"{tot[k][x + '_exact']}\t{tot[k][x + '_unmatched']}" for x in A) + "\n")
    print(open(f"{a.out}/arrays_summary.tsv").read(), end="")

def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sp = ap.add_subparsers(dest="cmd", required=True)
    for name in ("inputs", "classes", "proj"):
        p = sp.add_parser(name)
        p.add_argument("--truth", required=True); p.add_argument("--tie-split", required=True)
        p.add_argument("--bam", required=True); p.add_argument("--out", required=True)
        if name == "classes":
            p.add_argument("--annot", required=True); p.add_argument("--own", type=int, default=2)
            p.add_argument("--shadow", type=int, default=1); p.add_argument("--read-family")
        if name == "proj":
            p.add_argument("--fasta", required=True); p.add_argument("--p-gtf", required=True)
            p.add_argument("--threads", type=int, default=8); p.add_argument("--min-cov", type=float, default=0.9)
            p.add_argument("--min-id", type=float, default=0.85); p.add_argument("--pad", type=int, default=5000)
    p = sp.add_parser("score")
    p.add_argument("--out", required=True); p.add_argument("--annot", required=True); p.add_argument("--arm", action="append", required=True)
    p = sp.add_parser("arrays")
    p.add_argument("--annot", required=True); p.add_argument("--bam", required=True); p.add_argument("--out", required=True)
    p.add_argument("--min-units", type=int, default=5); p.add_argument("--pad", type=int, default=5000)
    p = sp.add_parser("arrays-proj")
    p.add_argument("--out", required=True); p.add_argument("--fasta", required=True); p.add_argument("--p-gtf", required=True)
    p.add_argument("--threads", type=int, default=8)
    p = sp.add_parser("arrays-score")
    p.add_argument("--out", required=True); p.add_argument("--annot", required=True); p.add_argument("--arm", action="append", required=True)
    a = ap.parse_args(argv)
    dict(inputs=cmd_inputs, classes=cmd_classes, proj=cmd_proj, score=cmd_score, arrays=cmd_arrays,
         **{"arrays-proj": cmd_arrays_proj, "arrays-score": cmd_arrays_score})[a.cmd](a)

if __name__ == "__main__":
    main()
