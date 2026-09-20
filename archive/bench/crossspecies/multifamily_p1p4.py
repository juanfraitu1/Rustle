#!/usr/bin/env python3
"""P1/P2/P4 on families OTHER than NPIP -- the still-unverified parts of the definition.

NPIP established: P1 seed-invariance (19/19 seeds -> identical gene set), P2 complete graph
(density 1.000), P4 substrate convergence (Jaccard 0.963, zero RNA-only edges). All three were
measured on ONE family. This runs the same three on PCDHB, NBPF, TUBA, TBC1D3, GOLGA.

PRE-REGISTERED EXPECTATION: NPIP is unusually clean (a ~16 kb cassette with no structural variation
between copies). PCDHB is a tandem array; GOLGA has 38 members with known sub-structure. P1 is
EXPECTED TO PARTLY FAIL. That is not a defeat -- it would mean F must be a CLOSURE (iterate seeding to
a fixed point) rather than a single pass, which is still constructive. Recording this before the run so
the outcome cannot be reinterpreted after the fact.

⚠⚠ TIER CORRECTION, 2026-08-10. The paragraph that used to sit here claimed "Edge rule = the SHIPPED
one: asm20 id>=0.80 UNION sensitive(-k11 -w5) id>=0.60, coverage >=0.50 AGGREGATED over records".
THAT WAS FALSE ON FOUR OF FIVE AXES, and every P1-P4 number this script has ever printed was computed
under the wrong rule. It ran `-N 200 -p 0.02` with NO `-X`; it kept an asm20 leg the shipped default
SKIPS; it used nmatch/blocklen where the binary uses 1-de; and "AGGREGATED over records" is
RUSTLE_ER_SUM_COVERAGE, which is DEFAULT OFF in the binary and carries a measured cross-family
precision cost. It also carried defect M1 (a query-axis numerator over a possibly-target denominator).

The rule now comes from ONE place -- bench/soto/rustlib.py, the Python mirror of `nucleotide_edges_scored`
-- so this script cannot drift from the binary again. Read that file for the full axis-by-axis table.
EVERY NUMBER THIS SCRIPT PRINTS MUST NAME ITS TIER; it now prints the argv it ran.

⚠ All seeds go through ONE minimap2 genome pass -- one index build, not one per seed.
⚠ The GENOME pass below is deliberately NOT changed: seed->genome mapping legitimately keeps
  `-N 200 -p 0.02` and no `-X`. Only the ALL-VS-ALL (E_r) pass is the shipped tier.
⚠ RNA nodes come from the FULL A119b BAM. Both Soto caches are -M -L subsets.
⚠ -F 2308 for the per-read CIGAR walk; N in an RNA CIGAR is an intron spliced OUT.
"""
import os
import subprocess
import sys
from collections import defaultdict

sys.path.insert(0, os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "soto"))
import rustlib as er_tier  # noqa: E402  -- bench/soto/rustlib.py is THE definition of the E_r rule

W = "/home/juanfra/winloci_scratch/seedfam/mf"
HS = "/mnt/c/Users/jfris/Desktop/Reference/chm13v2.0.fa"
BAM = "/mnt/linuxdisk/home/juanfraitu/winloci_data/A119b.t2t.bam"
GENES = "/home/juanfra/winloci_scratch/seedfam/hs_genes.bed"
HERE = os.path.dirname(os.path.abspath(__file__))
FAMS = sys.argv[1].split(",") if len(sys.argv) > 1 else \
    ["PCDHB", "NBPF", "TUBA", "TBC1D3", "GOLGA"]
T = "4"
MIN_BLOCK, MERGE_D, MIN_ALN, MIN_SPAN, MIN_READS = 1000, 10000, 5000, 5000, 3
os.makedirs(W, exist_ok=True)


def sh(cmd, out=None):
    if out:
        with open(out, "w") as fh:
            return subprocess.run(cmd, stdout=fh, stderr=subprocess.DEVNULL).returncode
    return subprocess.run(cmd, stdout=subprocess.DEVNULL,
                          stderr=subprocess.DEVNULL).returncode


genes = []
for line in open(GENES):
    f = line.rstrip("\n").split("\t")
    if len(f) >= 4:
        genes.append((f[0], int(f[1]), int(f[2]), f[3]))

members = defaultdict(list)
for c, s, e, g in genes:
    hit = [f for f in FAMS if g.startswith(f)]
    if hit:
        members[max(hit, key=len)].append((c, s, e, g))
# ⚠⚠ THE MEMBER SET IS BUILT BY PREFIX MATCH (`g.startswith(f)`) AND THAT IS NOT A FAMILY.
# Measured 2026-08-11 on CHM13: the `TBC1D3` panel picks up TBC1D30, TBC1D31 and TBC1D32 — three
# unrelated TBC-domain genes that merely share the string — and each of them seeds a DIFFERENT part
# of the genome, which is most of why P1 reads "NO" and why |V| reaches 205 for a 9-copy family.
# `TUBA` likewise picks up TUBA1B-AS1 (an antisense RNA) and TUBAL3; `PCDHB` picks up PCDHB1-AS1 and
# the cluster record PCDHB@. The seed set is an ANNOTATION-STRING artifact, it is INDEPENDENT of the
# alignment tier, and a P1 failure here must not be read as a failure of the definition until the
# member set is curated. Left as-is deliberately: changing it silently would move every number this
# script has ever printed. The offending symbols are printed below so no row can be quoted blind.
# ⚠ No classifier is attempted here: "is TBC1D30 a TBC1D3 paralog?" is not decidable from the string
#   (TBC1D30 is not; TBC1D3B is), so the symbols are simply PRINTED and the reader adjudicates.
for f in FAMS:
    members[f].sort()
    syms = sorted({g for _, _, _, g in members[f]})
    print(f"{f:<10} {len(members[f]):>3} members   {', '.join(syms)}")

# ---- ONE genome pass for every seed of every family -------------------------
seedfa = f"{W}/allseeds.fa"
if not os.path.exists(seedfa) or os.path.getsize(seedfa) == 0:
    with open(seedfa, "w") as fh:
        for f in FAMS:
            for c, s, e, g in members[f]:
                out = subprocess.run(["samtools", "faidx", HS, f"{c}:{s+1}-{e}"],
                                     capture_output=True, text=True).stdout.splitlines()
                if len(out) > 1:
                    fh.write(f">{f}|{g}\n" + "\n".join(out[1:]) + "\n")
n_seeds = sum(1 for l in open(seedfa) if l[0] == ">")
print(f"\nseeds: {n_seeds}   -> ONE minimap2 genome pass")
paf = f"{W}/allseeds.paf"
if not os.path.exists(paf) or os.path.getsize(paf) == 0:
    sh(["minimap2", "-x", "asm20", "-c", "--eqx", "-N", "200", "-p", "0.02",
        "-I", "2G", "-t", T, HS, seedfa], out=paf)
print(f"records: {sum(1 for _ in open(paf))}\n")

hits = defaultdict(list)
for line in open(paf):
    f = line.rstrip("\n").split("\t")
    if len(f) < 12:
        continue
    q, t, ts, te, nm, bl = f[0], f[5], int(f[7]), int(f[8]), int(f[9]), int(f[10])
    if bl >= MIN_BLOCK and nm / bl >= 0.80:
        hits[q].append((t, ts, te, nm))


def loci_for(seed):
    by = defaultdict(list)
    for t, ts, te, nm in hits.get(seed, []):
        by[t].append((ts, te, nm))
    out = []
    for c, v in by.items():
        v.sort()
        cs = ce = None
        aln = 0
        for s, e, nm in v:
            if cs is None:
                cs, ce, aln = s, e, nm
            elif s <= ce + MERGE_D:
                ce = max(ce, e)
                aln += nm
            else:
                if aln >= MIN_ALN and ce - cs >= MIN_SPAN:
                    out.append((c, cs, ce))
                cs, ce, aln = s, e, nm
        if cs is not None and aln >= MIN_ALN and ce - cs >= MIN_SPAN:
            out.append((c, cs, ce))
    return sorted(out)


def genes_in(L, fam):
    hit = set()
    for c, a, b in L:
        for gc, gs, ge, g in members[fam]:
            if gc == c and gs < b and a < ge:
                hit.add(g)
    return hit


# ---- P1 --------------------------------------------------------------------
print("=" * 78)
print("P1 SEED-INVARIANCE  (does every member seed return the same family?)")
print("=" * 78)
print(f"{'family':<10}{'seeds':>6}{'|F(s)| min-max':>16}{'genes recovered':>18}{'identical set?':>16}")
p1 = {}
for f in FAMS:
    res = {g: loci_for(f"{f}|{g}") for _, _, _, g in members[f]}
    gs = {g: genes_in(L, f) for g, L in res.items()}
    sizes = sorted(len(v) for v in res.values())
    allsame = len({frozenset(v) for v in gs.values()}) == 1
    full = sum(1 for v in gs.values() if len(v) == len(members[f]))
    p1[f] = (res, gs, allsame)
    print(f"{f:<10}{len(res):>6}{f'{sizes[0]}-{sizes[-1]}':>16}"
          f"{f'{full}/{len(res)} seeds -> all {len(members[f])}':>18}"
          f"{('YES' if allsame else 'NO'):>16}")
for f in FAMS:
    res, gs, allsame = p1[f]
    if not allsame:
        u = set().union(*gs.values()) if gs else set()
        worst = sorted(gs.items(), key=lambda kv: len(kv[1]))[:3]
        print(f"  {f}: union {len(u)}/{len(members[f])}; weakest seeds "
              + ", ".join(f"{k}({len(v)})" for k, v in worst))

# ---- node set per family = union over seeds (the CLOSURE, one iteration) ----
print("\n" + "=" * 78)
print("P2 GRAPH STRUCTURE  (SHIPPED rule -- tier and floors printed below)")
print("=" * 78)


def er_edges(prefix, ns):
    """E_r at the SHIPPED tier: ONE sensitive pass, single-record rule, M1-correct coverage."""
    return er_tier.er_edges([(f"{prefix}.paf", er_tier.SENSITIVE_IDENTITY)], nodes=set(ns),
                           min_coverage=er_tier.MIN_COVERAGE, denom="min")


comps = er_tier.components


def ava(src, prefix):
    """All-vs-all at the SHIPPED tier -- ONE pass, from bench/soto/rustlib.py."""
    er_tier.ava(src, f"{prefix}.paf", threads=int(T))


nodesets = {}
# PROVENANCE. Every number below is a statement about ONE tier; print which, so no row can be
# quoted without it. (Defect B1: the rows this script printed before 2026-08-10 were a DIFFERENT
# tier from the binary, and nothing in the output said so.)
print(f"\nE_r TIER: {er_tier.er_provenance(int(T))}")
print(f"E_r RULE: identity 1-de >= {er_tier.SENSITIVE_IDENTITY:.2f}, "
      f"coverage-of-shorter >= {er_tier.MIN_COVERAGE:.2f}, ANY single record")
# γ is reported at BOTH values because the spec and the code disagree (docs say 0.40,
# gw_family_catalog hardcodes 0.20). Reporting one alone has repeatedly produced unquotable rows.
print(f"{'family':<10}{'|V|':>5}{'edges':>8}{'density':>9}{'comps':>7}{'largest':>9}"
      f"{'g=0.40':>9}{'g=0.20':>9}")
for f in FAMS:
    res, gs, _ = p1[f]
    U = sorted({x for L in res.values() for x in L})
    regfile = f"{W}/{f}.regions"
    with open(regfile, "w") as fh:
        for c, a, b in U:
            fh.write(f"{c}:{a+1}-{b}\n")
    sh(["samtools", "faidx", HS, "-r", regfile], out=f"{W}/{f}.fa")
    sh(["samtools", "faidx", f"{W}/{f}.fa"])
    ns = [l.split("\t")[0] for l in open(f"{W}/{f}.fa.fai")]
    nodesets[f] = ns
    ava(f"{W}/{f}.fa", f"{W}/{f}_dna")
    E = er_edges(f"{W}/{f}_dna", set(ns))
    n = len(ns)
    cs = comps(ns, E)
    big = cs[0]
    d = er_tier.density(big, E)
    print(f"{f:<10}{n:>5}{len(E):>8}{d:>9.3f}{len(cs):>7}{len(big):>9}"
          f"{d/0.40:>8.2f}x{d/0.20:>8.2f}x")

# ---- P4 --------------------------------------------------------------------
print("\n" + "=" * 78)
print("P4 SUBSTRATE CONVERGENCE  (same nodes; genomic sequence vs read-supported exons)")
print("=" * 78)
print(f"{'family':<10}{'|V|':>5}{'DNA E':>8}{'RNA E':>8}{'shared':>8}{'DNAonly':>9}{'RNAonly':>9}{'Jaccard':>9}")
for f in FAMS:
    ns = nodesets[f]
    rna = f"{W}/{f}_rna.fa"
    with open(rna, "w") as fh:
        for reg in ns:
            c, rng = reg.rsplit(":", 1)
            a, b = rng.split("-")
            a, b = int(a) - 1, int(b)
            p = subprocess.run(["samtools", "view", "-F", "2308", BAM, f"{c}:{a+1}-{b}"],
                               capture_output=True, text=True)
            ex = subprocess.run(["python3", f"{HERE}/read_exons.py"], input=p.stdout,
                                capture_output=True, text=True).stdout
            iv = defaultdict(int)
            for line in ex.splitlines():
                x = line.split("\t")
                if len(x) >= 3:
                    iv[(int(x[1]), int(x[2]))] += 1
            blocks = []
            cur = None
            for (s, e) in sorted(iv):
                s2, e2 = max(s, a), min(e, b)
                if e2 - s2 < 20:
                    continue
                if cur and s2 <= cur[1]:
                    cur = (cur[0], max(cur[1], e2), cur[2] + iv[(s, e)])
                else:
                    if cur and cur[2] >= MIN_READS:
                        blocks.append(cur[:2])
                    cur = (s2, e2, iv[(s, e)])
            if cur and cur[2] >= MIN_READS:
                blocks.append(cur[:2])
            if not blocks:
                continue
            seq = []
            for s, e in blocks:
                o = subprocess.run(["samtools", "faidx", HS, f"{c}:{s+1}-{e}"],
                                   capture_output=True, text=True).stdout.splitlines()
                seq.append("".join(o[1:]))
            fh.write(f">{reg}\n" + "".join(seq) + "\n")
    have = [l[1:].strip() for l in open(rna) if l[0] == ">"]
    if len(have) < 2:
        print(f"{f:<10}{len(ns):>5}   (only {len(have)} nodes have read support -- skipped)")
        continue
    ava(rna, f"{W}/{f}_rna")
    Ed = er_edges(f"{W}/{f}_dna", set(have))
    Er = er_edges(f"{W}/{f}_rna", set(have))
    inter = Ed & Er
    print(f"{f:<10}{len(have):>5}{len(Ed):>8}{len(Er):>8}{len(inter):>8}"
          f"{len(Ed-Er):>9}{len(Er-Ed):>9}{len(inter)/max(len(Ed|Er),1):>9.3f}")
print("\nDONE")
