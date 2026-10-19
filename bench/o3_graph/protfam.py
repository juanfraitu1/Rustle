#!/usr/bin/env python3
"""Gene families and the multigene families inside them (2026-10-18; user: "the definition of multigene families should be the
ones the Makova group uses in both IsoCon and Yoo's papers" -- then: "those are not multigene families, they are gene families").

Level 1, the GENE FAMILY -- Yoo et al. 2025 (Nature, Supplementary VIII, the Canzar lab's section): blastp of the protein
sequences of all protein-coding genes (longest isoform per gene), homology at >= 50% identity with the alignment covering >= 35%
of BOTH proteins (the user's reading of "35% protein"), families the connected components (transitive closure). Here CAT's
protein-coding genes, each its longest CDS translated; every coding gene is a node (no span cap, no one-gene-per-locus drop).
Yoo's second step -- paralogous loci missing from the annotation or pseudogenized, found by read bundles -- is taken by the
near-identical level (annotated pseudogenes join their coding copies there) and, at run time, by O1's genome search.
(Yoo's supplement calls these "multicopy gene families"; at 50% protein identity they are homology gene families.)
Level 2, the MULTIGENE FAMILY in IsoCon's sense (Sahlin et al. 2018: the highly similar duplicate copies -- the nine Y-ampliconic
families; reads among its copies are ambiguous, and O2 and O3 act there): seqfamilies.py's rule (sequence >= 0.90 over >= 0.8 and the carried intron chain, every
chain kept, non-coding members). Level 1's components include level 2's links, so every multigene family nests inside one
gene family.

    protfam.py build OUT.json LEVEL2_POOL.json [--threads N] [--protein-links EARLIER.json.links.tsv]
        -> OUT.json (gene families: genes, the multigene families nested in each), OUT.json.proteins.fa, OUT.json.links.tsv
    protfam.py check OUT.json
        -> the mainstay and IsoCon's nine Y-ampliconic families at both levels (member_report), OUT.json.check.json
"""
import collections
import csv
import importlib.util
import json
import os
import subprocess
import sys
import tempfile

HERE = os.path.dirname(os.path.abspath(__file__))
MIN_IDENT, MIN_COVER = 50.0, 0.35
BLAST = "/home/juanfra/miniforge3/envs/blast/bin"
CAT_FULL = "/mnt/linuxdisk/home/juanfraitu/winloci_data/gencode_chm13/chm13v2.0_gencode.gff3"
_COMP = str.maketrans("ACGTNacgtn", "TGCANtgcan")
_CODE = {a + b + c: aa for (a, b, c), aa in zip(
    ((x, y, z) for x in "TCAG" for y in "TCAG" for z in "TCAG"),
    "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG")}


def translate_cds(blocks, strand, fetch, phase=0):
    """The protein of CDS blocks (0-based half-open, any order) on `strand`; fetch(start, end) -> genome bases. phase: bases
    to skip at the CDS's 5' end. A trailing stop is dropped and an internal one becomes X (blastp's unknown residue)."""
    seq = "".join(fetch(s, e) for s, e in sorted(blocks)).upper()
    if strand == "-":
        seq = seq[::-1].translate(_COMP)
    seq = seq[phase:]
    aa = "".join(_CODE.get(seq[i:i + 3], "X") for i in range(0, len(seq) - 2, 3))
    return aa[:-1].replace("*", "X") if aa.endswith("*") else aa.replace("*", "X")


def longest_cds(records):
    """records (gene, transcript, (start, end)) of CDS blocks -> {gene: the transcript with the longest CDS} (ties: the first
    transcript id)."""
    tot = collections.Counter()
    gene_of = {}
    for g, t, (s, e) in records:
        tot[t] += e - s
        gene_of[t] = g
    best = {}
    for t in sorted(tot):
        g = gene_of[t]
        if g not in best or tot[t] > tot[best[g]]:
            best[g] = t
    return best


def protein_links(hits, lens, min_ident=MIN_IDENT, min_cover=MIN_COVER):
    """{(a, b), a < b} of blastp hits (query, subject, % identity, qstart, qend, sstart, send; 1-based) at >= min_ident whose
    alignment covers >= min_cover of both proteins."""
    out = set()
    for q, s, ident, qs, qe, ss, se in hits:
        if q == s or ident < min_ident:
            continue
        if (abs(qe - qs) + 1) / lens[q] >= min_cover and (abs(se - ss) + 1) / lens[s] >= min_cover:
            out.add(tuple(sorted((q, s))))
    return out


def families(prot_links, near_links, level2=None):
    """Level-1 families: connected components (>= 2 genes) of the protein links and the near-identical links, each with its
    genes and the level-2 subfamilies inside it ({name: genes})."""
    par = {}

    def find(x):
        par.setdefault(x, x)
        while par[x] != x:
            par[x] = par[par[x]]
            x = par[x]
        return x
    for a, b in list(prot_links) + list(near_links):
        par[find(a)] = find(b)
    comp = collections.defaultdict(list)
    for x in list(par):
        comp[find(x)].append(x)
    out = []
    for genes in sorted((sorted(c) for c in comp.values() if len(c) >= 2), key=lambda c: c[0]):
        gs = set(genes)
        subs = sorted(n for n, m in (level2 or {}).items() if set(m) <= gs)
        out.append(dict(genes=genes, subfamilies=subs))
    return out


KNOWN = {"NPIP": r"^NPIP[AB]", "TBC1D3": r"^TBC1D3([A-L]|P[0-9]+)?$", "GOLGA8": r"^GOLGA8", "GOLGA6L": r"^GOLGA6L",
         "CT47A": r"^CT47A", "XAGE": r"^XAGE", "SPATA31": r"^SPATA31",                           # the mainstay set
         "BPY2": r"^BPY2[BC]?$", "CDY": r"^CDY[12][AB]?$", "DAZ": r"^DAZ[1-4]$", "HSFY": r"^HSFY[12]$", "PRY": r"^PRY2?$",
         "RBMY": r"^RBMY1[A-J][0-9]?$", "TSPY": r"^TSPY[0-9]+$", "VCY": r"^VCY(1B)?$", "XKRY": r"^XKRY2?$"}   # IsoCon's nine


def member_report(known, genes, l1, l2, spans=None):
    """For each known family ({name: regex on gene names}; names only identify the truth here): its members ({gene id: name}
    matching), how the two levels place them -- level 1: the families holding members, the most members in one, the members
    in none, the other genes beside the most members; level 2 ({subfamily: gene ids}): the same for the subfamilies. spans
    ({gene id: (contig, start, end)}): a member in none whose span overlaps a placed member's is a second annotation of that
    locus (CAT names RBMY1A1's locus RBMY1B too), reported as a duplicate, not a miss."""
    import re

    def split(missing, placed):
        if not spans:
            return sorted(genes[g] for g in missing), []
        dup = {g for g in missing if any(spans[g][0] == spans[p][0] and spans[g][1] < spans[p][2] and spans[p][1] < spans[g][2]
                                         for p in placed)}
        return sorted(genes[g] for g in missing - dup), sorted(genes[g] for g in dup)
    out = {}
    for fam, pat in known.items():
        mem = {g for g, n in genes.items() if re.search(pat, n)}
        hold = [(len(mem & set(f["genes"])), f) for f in l1 if mem & set(f["genes"])]
        best = max(hold, key=lambda x: x[0]) if hold else (0, dict(genes=[]))
        in1 = set().union(*(set(f["genes"]) for _, f in hold)) if hold else set()
        subs = [(len(mem & set(v)), k) for k, v in l2.items() if mem & set(v)]
        in2 = set().union(*(set(l2[k]) for _, k in subs)) if subs else set()
        m1, d1 = split(mem - in1, mem & in1)
        m2, d2 = split(mem - in2, mem & in2)
        out[fam] = dict(members=len(mem), level1_families=len(hold), level1_largest=best[0], level1_missing=m1,
                        level1_duplicates=d1, level1_other_genes=len(set(best[1]["genes"]) - mem),
                        level2_subfamilies=len(subs), level2_largest=max((n for n, _ in subs), default=0),
                        level2_missing=m2, level2_duplicates=d2)
    return out


def cds_records(gff=CAT_FULL):
    """Protein-coding CDS blocks of CAT's full GFF: [(gene id, transcript id, chrom, strand, (start0, end), phase)]."""
    out = []
    with open(gff) as fh:
        for ln in fh:
            f = ln.split("\t", 8)
            if len(f) < 9 or f[2] != "CDS" or "transcript_biotype=protein_coding" not in f[8]:
                continue
            at = dict(x.split("=", 1) for x in f[8].rstrip("\n").split(";") if "=" in x)
            out.append((at.get("gene_id"), at.get("transcript_id"), f[0], f[6], (int(f[3]) - 1, int(f[4])),
                        int(f[7]) if f[7] in "012" else 0))
    return out


def cmd_build(out_p, level2_pool, threads=4, links_path=None):
    import pysam
    spec = importlib.util.spec_from_file_location("o3_graph_families", os.path.join(HERE, "families.py"))
    F = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(F)
    if links_path:                                        # the protein level of an earlier build (blastp: ~35 min)
        links = {(r["a"], r["b"]) for r in csv.DictReader(open(links_path), delimiter="\t")}
        return _families_out(out_p, level2_pool, links)
    recs = cds_records()
    best = longest_cds([(g, t, b) for g, t, _, _, b, _ in recs])
    keep = {t for t in best.values()}
    blocks, info = collections.defaultdict(list), {}
    for g, t, ch, st, b, ph in recs:
        if t in keep:
            blocks[t].append((b, ph))
            info[t] = (g, ch, st)
    fa = pysam.FastaFile(F.GENOME)
    prot = {}
    for t, bl in blocks.items():
        g, ch, st = info[t]
        bl.sort()
        first = bl[-1] if st == "-" else bl[0]                       # the CDS's 5' block carries the phase
        aa = translate_cds([b for b, _ in bl], st, lambda s, e, ch=ch: fa.fetch(ch, s, e), phase=first[1])
        if len(aa) >= 10:
            prot[g] = aa
    with open(f"{out_p}.proteins.fa", "w") as fh:
        fh.writelines(f">{g}\n{aa}\n" for g, aa in sorted(prot.items()))
    print(f"[protfam] {len(prot)} protein-coding genes translated", flush=True)
    with tempfile.TemporaryDirectory(dir=os.environ.get("TMPDIR")) as d:
        subprocess.run([f"{BLAST}/makeblastdb", "-in", f"{out_p}.proteins.fa", "-dbtype", "prot", "-out", f"{d}/db"],
                       check=True, capture_output=True)
        res = subprocess.run([f"{BLAST}/blastp", "-query", f"{out_p}.proteins.fa", "-db", f"{d}/db", "-evalue", "1e-5",
                              "-max_target_seqs", "5000", "-num_threads", str(threads),
                              "-outfmt", "6 qseqid sseqid pident qstart qend sstart send"],
                             check=True, capture_output=True, text=True).stdout
    hits = [(x[0], x[1], float(x[2]), int(x[3]), int(x[4]), int(x[5]), int(x[6]))
            for x in (ln.split("\t") for ln in res.splitlines())]
    lens = {g: len(aa) for g, aa in prot.items()}
    links = protein_links(hits, lens)
    with open(f"{out_p}.links.tsv", "w") as fh:
        fh.write("a\tb\n")
        fh.writelines(f"{a}\t{b}\n" for a, b in sorted(links))
    return _families_out(out_p, level2_pool, links)


def _families_out(out_p, level2_pool, links):
    l2 = json.load(open(level2_pool))
    (pool,) = l2["sets"].values()
    level2 = {f["family"]: [c["gene_id"] for c in f["copies"]] for f in pool}
    near = [(r["a"], r["b"]) for r in csv.DictReader(open(f"{level2_pool}.links.tsv"), delimiter="\t") if r["struct_link"] == "1"]
    fams = families(links, near, level2)
    json.dump(dict(rule="protfam: blastp >= 50% identity over >= 35% of both proteins (Yoo et al. 2025), components with the "
                        "near-identical links of " + os.path.basename(level2_pool), min_ident=MIN_IDENT, min_cover=MIN_COVER,
                   level2_pool=level2_pool, families=fams), open(out_p, "w"), indent=1)
    sizes = collections.Counter(len(f["genes"]) for f in fams)
    print(f"[protfam] {len(links)} protein links; {len(fams)} families ({sum(sizes.values())} with >= 2 genes; largest "
          f"{max(sizes) if sizes else 0}); {sum(len(f['subfamilies']) for f in fams)} near-identical subfamilies nested -> {out_p}",
          flush=True)
    return 0


def cmd_check(out_p):
    """The known families (KNOWN) at both levels, against CAT's genes of every biotype."""
    spec = importlib.util.spec_from_file_location("o3_graph_families", os.path.join(HERE, "families.py"))
    F = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(F)
    d = json.load(open(out_p))
    rows = list(csv.DictReader(open(f"{F.CAT}/chm13v2.0_CAT_Liftoff.genes.tsv"), delimiter="\t"))
    genes = {r["gene_id"]: r["gene_name"] for r in rows}
    spans = {r["gene_id"]: (r["chrom"], int(r["start0"]), int(r["end"])) for r in rows}
    (pool,) = json.load(open(d["level2_pool"]))["sets"].values()
    l2 = {f["family"]: [c["gene_id"] for c in f["copies"]] for f in pool}
    rep = member_report(KNOWN, genes, d["families"], l2, spans=spans)
    for fam, r in rep.items():
        print(f"{fam:8s} members {r['members']:3d} | multigene family: {r['level2_largest']:3d} in its largest of "
              f"{r['level2_subfamilies']}, missing {len(r['level2_missing'])} {r['level2_missing'][:6]}, duplicate annotations "
              f"{r['level2_duplicates']} | gene family: "
              f"{r['level1_largest']:3d} in its largest of {r['level1_families']} (+{r['level1_other_genes']} other genes), "
              f"missing {len(r['level1_missing'])}", flush=True)
    json.dump(rep, open(f"{out_p}.check.json", "w"), indent=1)
    return 0


def build_args(argv):
    """`build OUT.json LEVEL2_POOL.json [--threads N] [--protein-links LINKS.tsv]` (after `build`) -> (OUT, LEVEL2_POOL,
    threads, protein links to reuse or None)."""
    out, pool = argv[0], argv[1]
    threads = int(argv[argv.index("--threads") + 1]) if "--threads" in argv else 4
    links = argv[argv.index("--protein-links") + 1] if "--protein-links" in argv else None
    return out, pool, threads, links


def main(argv):
    if argv and argv[0] == "check" and len(argv) == 2:
        return cmd_check(argv[1])
    if argv and argv[0] == "build" and len(argv) >= 3:
        out, pool, threads, links = build_args(argv[1:])
        return cmd_build(out, pool, threads=threads, links_path=links)
    raise SystemExit(__doc__)


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
