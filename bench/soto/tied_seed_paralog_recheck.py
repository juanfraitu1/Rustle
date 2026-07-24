#!/usr/bin/env python3
"""Fix the tied-seed eval's paralog-count metric.

The tied_seed_eval.md verification scored each --tied-seed over-seed by aligning its consensus to the
genome + ANNOTATION and eyeballing paralogs; it labelled os1 (chr11:89,746,906) "0 paralogs, weak
homology -> likely phantom". That was wrong: os1's true paralog (os4, chr11:89,934,235) is itself an
UNANNOTATED gap copy, so an annotation-anchored search cannot see it.

This re-scores paralogy with the PIPELINE'S OWN homology oracle: minimap2 asm20 (`-c --eqx -x asm20
--secondary=no`, the exact call in copy_assign_pipeline.rs) among the emitted copies, family edge =
identity >= 0.80 AND query-coverage >= 0.50 (the --refine family-edge oracle). A copy's paralog count is
its degree in that homology graph.

Read-only; no copy_assign. Usage: python3 bench/soto/tied_seed_paralog_recheck.py
"""
import subprocess, tempfile, os, sys
from collections import defaultdict

FASTA = "/mnt/linuxdisk/home/juanfraitu/winloci_data/chm13v2.0.fa"
QUANT = "/home/juanfra/winloci_scratch/ts_on.quant.tsv"      # --tied-seed copy set (17 copies)
ID_MIN, COV_MIN = 0.80, 0.50
CHROM = "chr11"

# the 6 copies that appear ONLY with --tied-seed (baseline->tied diff), with their eval labels
NEW = {
    89746906: "os1 (eval:0-paralog/PHANTOM)",
    89785154: "TRIM64B (Soto member)",
    89832458: "TRIM49D1 (gene)",
    89849042: "TRIM49D2 (gene)",
    89893563: "TRIM64 (Soto member)",
    89934235: "os4 (eval:real-ish)",
}

def sh(cmd):
    return subprocess.run(cmd, capture_output=True, text=True).stdout

# 1. load emitted copies
copies = []  # (start,end,label)
with open(QUANT) as f:
    next(f)
    for line in f:
        c = line.rstrip("\n").split("\t")
        s, e = int(c[4]), int(c[5])
        copies.append((s, e))
copies.sort()
def name(s): return f"{s}"

# 2. write each copy's reference sequence to one multi-fasta
tmp = tempfile.mkdtemp()
fa = os.path.join(tmp, "copies.fa")
with open(fa, "w") as out:
    for s, e in copies:
        seq = sh(["samtools", "faidx", FASTA, f"{CHROM}:{s+1}-{e}"])
        seq = "".join(l for l in seq.splitlines() if not l.startswith(">"))
        out.write(f">{name(s)}\n{seq}\n")

# 3. all-vs-all asm20 (the pipeline oracle), parse PAF -> homology edges
paf = sh(["minimap2", "-c", "--eqx", "-x", "asm20", "-t", "4", "-p", "0.1", "-N", "50", fa, fa])
edges = defaultdict(dict)   # a -> {b: (id,cov)}
best = defaultdict(dict)
for line in paf.splitlines():
    p = line.split("\t")
    if len(p) < 12: continue
    q, qlen, qs, qe = p[0], int(p[1]), int(p[2]), int(p[3])
    t, tlen = p[5], int(p[6])
    nmatch, alnlen = int(p[9]), int(p[10])
    if q == t or alnlen == 0: continue
    idt = nmatch / alnlen
    cov = (qe - qs) / qlen
    # keep the best (by id*cov) per unordered pair, symmetric
    for a, b, c_ in ((q, t, cov), (t, q, (qe-qs)/tlen)):
        cur = best[a].get(b, (0, 0))
        if idt * c_ > cur[0] * cur[1]:
            best[a][b] = (idt, c_)

def paralogs(a):
    out = []
    for b, (idt, cov) in best[a].items():
        if idt >= ID_MIN and cov >= COV_MIN:
            out.append((b, idt, cov))
    return sorted(out, key=lambda x: -x[1])

# 4. report — focus on the 6 new copies
print(f"# tied-seed paralog re-check via pipeline homology oracle (asm20, id>={ID_MIN}, cov>={COV_MIN})")
print(f"# copies: {len(copies)}  (all-vs-all)\n")
print(f"{'copy_start':>11}  {'label':<28}{'n_paralogs':>10}  paralog_partners(id)")
for s, e in copies:
    if s not in NEW: continue
    pl = paralogs(name(s))
    partners = ", ".join(f"{b}({idt:.3f})" for b, idt, cov in pl) or "-"
    flag = ""
    if s == 89746906:
        flag = "  <== eval said 0 paralogs; oracle finds " + str(len(pl))
    print(f"{s:>11,}  {NEW[s]:<28}{len(pl):>10}  {partners}{flag}")

print("\n# ALL copies (paralog degree in the emitted homology graph):")
for s, e in copies:
    pl = paralogs(name(s))
    tag = f"  [{NEW[s]}]" if s in NEW else ""
    print(f"  {s:>11,}: {len(pl)} paralog(s){tag}" + (("  -> "+", ".join(f'{b}' for b,_,_ in pl)) if pl else ""))

# 5. persist the 6 over-seeds' corrected paralog counts
tsv = os.path.join(os.path.dirname(os.path.abspath(__file__)), "tied_seed_paralog_recheck.tsv")
with open(tsv, "w") as out:
    out.write("copy_start\tlabel\toracle_n_paralogs\tparalog_partner\tpartner_id\teval_verdict\tcorrected\n")
    for s, e in copies:
        if s not in NEW: continue
        pl = paralogs(name(s))
        partner = pl[0][0] if pl else "-"; pid = f"{pl[0][1]:.4f}" if pl else "-"
        ev = "0-paralog/PHANTOM" if s == 89746906 else ("real-ish" if s == 89934235 else "real")
        corrected = "REAL (1 paralog)" if pl else "still 0"
        out.write(f"{s}\t{NEW[s]}\t{len(pl)}\t{partner}\t{pid}\t{ev}\t{corrected}\n")
print(f"\nwrote {tsv}")
