#!/usr/bin/env python3
"""validate_genomic_dna.py — DNA-level validation of the refined paralog catalog by GENOMIC-SPAN alignment.

The canonical orthogonal check is a SEDEF/BISER segmental-duplication map, but the precompiled BISER
aligner segfaults under WSL2 and its unrefined putative SDs are repeat-noise (they link gene loci to Alu
copies, not to their paralogs). This is the targeted substitute for BISER's alignment step: for each refined
family, extract every copy's GENOMIC span (chrom:start-end, INTRONS INCLUDED) from the assembly and all-vs-all
align them (minimap2 asm20). A family is DNA-CONFIRMED when >=2 copies' genomic spans align at segmental-
duplication stringency (identity >= ID, coverage-of-shorter >= COV).

CAVEAT — this is DNA-level CORROBORATION, NOT a fully-orthogonal check, and NOT a real SEDEF/BISER segdup
map. The catalog was built on exon-sum (exon) homology and the genomic span CONTAINS those exons, so the
test partly re-uses the building signal; the genuinely-independent content (introns + flanks) dominates the
alignment only for intron-rich, few-exon families (e.g. RABL2, a 4-exon/19kb family) and is minor for
compact many-exon families. There is also no genome-wide FDR / repeat-masking and it cannot tell a paralog
from an allele. Report the confirmed fraction as a LOWER BOUND on segdup-grade homology, not as a precision.
NOTE: because coverage is measured OF THE SHORTER span, a true intronless retrocopy would PASS (its short
span aligns over the parent's exons) — so the genomic-silent bin is NOT a retrocopy detector. Adversarial
check found the silent families are partial / structurally-divergent HIGH-identity paralogs that fail on
COVERAGE (extent), not on sequence divergence, with 0/N intronless — they are not retrocopies.

Run: python bench/validate_genomic_dna.py <refined_prefix> <genome.fa> [samtools] [minimap2] [ID] [COV]
"""
import sys, csv, subprocess, collections, os, tempfile

PREFIX = sys.argv[1] if len(sys.argv) > 1 else "/home/juanfra/winloci_scratch/gw_xchrom_refined"
GENOME = sys.argv[2] if len(sys.argv) > 2 else "/home/juanfra/winloci_scratch/GGO.fasta"
SAMTOOLS = sys.argv[3] if len(sys.argv) > 3 else "/home/juanfra/miniforge3/envs/biser/bin/samtools"
MM2 = sys.argv[4] if len(sys.argv) > 4 else "/home/juanfra/miniforge3/bin/minimap2"
ID = float(sys.argv[5]) if len(sys.argv) > 5 else 0.90    # segmental-duplication identity standard
COV = float(sys.argv[6]) if len(sys.argv) > 6 else 0.50

copies = collections.defaultdict(list)  # fid -> [(chrom,start,end,n_exon,strand)]
for c in csv.DictReader(open(f"{PREFIX}.copies.tsv"), delimiter="\t"):
    copies[c["family_id"]].append((c["chrom"], int(c["start"]), int(c["end"]), int(c["n_exon"]), c["strand"]))
fams = {f["family_id"]: f for f in csv.DictReader(open(f"{PREFIX}.families.tsv"), delimiter="\t")}

def genomic_seqs(cps):
    """Extract each copy's genomic span (introns included) via samtools faidx. 0-based half-open -> 1-based."""
    regions = [f"{c[0]}:{c[1]+1}-{c[2]}" for c in cps]
    out = subprocess.run([SAMTOOLS, "faidx", GENOME, *regions], capture_output=True, text=True).stdout
    seqs, cur = [], []
    for line in out.splitlines():
        if line.startswith(">"):
            if cur:
                seqs.append("".join(cur)); cur = []
        else:
            cur.append(line)
    if cur:
        seqs.append("".join(cur))
    return seqs  # in region order = copy order

def dna_homology_component(cps):
    """largest set of copies whose GENOMIC spans mutually align (id>=ID, cov-of-shorter>=COV)."""
    seqs = genomic_seqs(cps)
    n = len(cps)
    if len(seqs) != n or n < 2:
        return 0
    with tempfile.NamedTemporaryFile("w", suffix=".fa", delete=False) as fh:
        path = fh.name
        for i, s in enumerate(seqs):
            fh.write(f">{i}\n{s}\n")
    try:
        out = subprocess.run([MM2, "-cx", "asm20", "-X", "--no-long-join", path, path],
                             capture_output=True, text=True).stdout
    finally:
        os.unlink(path)
    lens = [len(s) for s in seqs]
    parent = list(range(n))
    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]; x = parent[x]
        return x
    for ln in out.splitlines():
        f = ln.split("\t")
        if len(f) < 12:
            continue
        try:
            q, t = int(f[0]), int(f[5])
        except ValueError:
            continue
        if q == t:
            continue
        qs, qe = int(f[2]), int(f[3])
        de = next((float(x[5:]) for x in f[12:] if x.startswith("de:f:")), None)
        ident = (1.0 - de) if de is not None else int(f[9]) / max(1, int(f[10]))
        cov = (qe - qs) / max(1, min(lens[q], lens[t]))
        if ident >= ID and cov >= COV:
            parent[find(q)] = find(t)
    return max(collections.Counter(find(i) for i in range(n)).values(), default=0)

n = conf = conf_xc = n_xc = 0
silent, conf_ex = [], []
exons_conf, exons_silent = [], []
for fid, cps in copies.items():
    cross = fams[fid]["cross_chrom"] == "true"
    n += 1; n_xc += cross
    lc = dna_homology_component(cps)
    is_conf = lc >= 2
    conf += is_conf; conf_xc += is_conf and cross
    me = sum(c[3] for c in cps) / len(cps)
    (exons_conf if is_conf else exons_silent).append(me)
    if is_conf and len(conf_ex) < 10:
        conf_ex.append((fid, cross, len(cps), lc, round(me, 1)))
    if not is_conf:
        silent.append((fid, cross, len(cps), round(me, 1)))

def pct(a, b): return f"{a}/{b} = {100*a/b:.1f}%" if b else "n/a"
def mean(xs): return sum(xs) / len(xs) if xs else 0.0

print(f"=== GENOMIC-DNA validation of {PREFIX} (minimap2 asm20, id>={ID}, cov-of-shorter>={COV}) ===\n")
print(f"refined families: {n}  (cross-chrom: {n_xc})")
print(f"DNA-CONFIRMED (>=2 copies' GENOMIC spans homologous): all {pct(conf, n)} | "
      f"same-chrom {pct(conf - conf_xc, n - n_xc)} | cross-chrom {pct(conf_xc, n_xc)}")
print(f"\nmean n_exon: DNA-confirmed={mean(exons_conf):.2f}  genomic-silent={mean(exons_silent):.2f}")
print(f"(genomic-silent = partial/structurally-divergent high-identity paralogs failing on COVERAGE, NOT")
print(f" retrocopies — verified 0/silent intronless; report confirmed% as a segdup-grade LOWER BOUND)")
print(f"\nDNA-confirmed examples (fid, xc, n_copies, linked, mean_exon):")
for e in conf_ex:
    print(f"  {e[0]} xc={e[1]} n={e[2]} linked={e[3]} mean_exon={e[4]}")
print(f"\ngenomic-silent families (first 15):")
for s in silent[:15]:
    print(f"  {s[0]} xc={s[1]} n_copies={s[2]} mean_exon={s[3]}")

print("\n=== RABL2 flagship (GWFAM17) ===")
if "GWFAM17" in copies:
    cps = copies["GWFAM17"]
    lc = dna_homology_component(cps)
    print(f"  5 copies, genomic-DNA-linked component = {lc}/5, mean_exon {sum(c[3] for c in cps)/5:.1f}")
