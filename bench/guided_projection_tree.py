#!/usr/bin/env python3
"""Prereg Addendum Q: subfamily trees from a reference-projected alignment (minimap2 asm20 of every member onto one
reference member, CIGAR-stacked on reference columns) instead of MAFFT, on Addendum P1's members; compares each
literature-group call with P1's MAFFT tree for the same run.

usage: guided_projection_tree.py <workdir> <refseq gff> <genes.gff.gz> <genome.fa> <iqtree3>
"""
import os
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))
src = open(os.path.join(HERE, "guided_tree_width.py")).read().split("mismatch = 0")[0]
argv = list(sys.argv)
exec(compile(src, "guided_tree_width.py(setup)", "exec"))
sys.argv = argv

PD = f"{W}/tree_proj"
os.makedirs(PD, exist_ok=True)


def rc(s):
    return s.translate(COMP)[::-1]


def project_alignment(tag, members):
    names = list(members)
    fa = f"{PD}/{tag}.fa"
    with open(fa, "w") as fh:
        for n in names:
            fh.write(f">{n}\n{members[n]}\n")
    paf = subprocess.run(["minimap2", "-c", "-x", "asm20", "-X", "-N", "50", "-p", "0.1", "-t", "4", fa, fa],
                         capture_output=True, text=True, check=True).stdout
    aligned = collections.Counter()
    for line in paf.splitlines():
        f = line.split("\t")
        if f[0] == f[5]:
            continue
        aligned[f[0]] += int(f[3]) - int(f[2])
        aligned[f[5]] += int(f[8]) - int(f[7])
    ref = max(sorted(names), key=lambda n: aligned[n])
    ref_fa, oth_fa = f"{PD}/{tag}.ref.fa", f"{PD}/{tag}.others.fa"
    with open(ref_fa, "w") as fh:
        fh.write(f">{ref}\n{members[ref]}\n")
    with open(oth_fa, "w") as fh:
        for n in names:
            if n != ref:
                fh.write(f">{n}\n{members[n]}\n")
    paf = subprocess.run(["minimap2", "-c", "-x", "asm20", "-N", "50", "-p", "0.1", "-t", "4", ref_fa, oth_fa],
                         capture_output=True, text=True, check=True).stdout
    by_q = collections.defaultdict(list)
    for line in paf.splitlines():
        f = line.split("\t")
        AS = next((int(x[5:]) for x in f[12:] if x.startswith("AS:i:")), 0)
        cg = next((x[5:] for x in f[12:] if x.startswith("cg:Z:")), "")
        by_q[f[0]].append((AS, int(f[1]), int(f[2]), int(f[3]), f[4], int(f[7]), cg))
    Lr = len(members[ref])
    rows = {ref: list(members[ref])}
    for n in names:
        if n == ref:
            continue
        row, filled = ["-"] * Lr, bytearray(Lr)
        seq = members[n]
        for AS, qlen, qs, qe, strand, ts, cg in sorted(by_q.get(n, []), key=lambda r: -r[0]):
            s = seq if strand == "+" else rc(seq)
            q = qs if strand == "+" else qlen - qe
            t = ts
            for k, op in re.findall(r"(\d+)([MID])", cg):
                k = int(k)
                if op == "M":
                    for j in range(k):
                        if not filled[t + j]:
                            row[t + j] = s[q + j]
                            filled[t + j] = 1
                    t += k
                    q += k
                elif op == "D":
                    t += k
                else:
                    q += k
        rows[n] = row
    keep = [j for j in range(Lr) if sum(1 for n in names if rows[n][j] != "-") > 0.5 * len(names)]
    aln = f"{PD}/{tag}.proj.fa"
    with open(aln, "w") as fh:
        for n in names:
            fh.write(f">{n}\n{''.join(rows[n][j] for j in keep)}\n")
    return aln, ref, Lr, len(keep)


def iqtree(tag, aln):
    t0 = time.time()
    subprocess.run([IQTREE, "-s", aln, "-m", "MFP", "-B", "1000", "-alrt", "1000", "-T", "4", "--seed", "1", "--prefix",
                    f"{PD}/{tag}", "-redo", "-quiet"], check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL,
                   timeout=600)
    model = next((l.split(":", 1)[1].strip() for l in open(f"{PD}/{tag}.iqtree") if l.startswith("Best-fit model")), "?")
    return parse_newick(open(f"{PD}/{tag}.treefile").read()), time.time() - t0, model


out = []
for tag, fam, label_of, members in runs:
    if len(members) < 4:
        continue
    t0 = time.time()
    aln, ref, Lr, kept = project_alignment(tag, members)
    t_aln = time.time() - t0
    tree, t_tree, model = iqtree(tag, aln)
    q = tree_scores(tree, fam, label_of)
    p1 = None
    if os.path.exists(f"{W}/tree/{tag}.treefile"):
        p1 = tree_scores(parse_newick(open(f"{W}/tree/{tag}.treefile").read()), fam, label_of)
    out.append((tag, fam, ref, Lr, kept, t_aln, t_tree, model, q, p1, len(members), sum(1 for m in members if label_of.get(m))))
    print(f"[proj] {tag}: reference {label_of.get(ref, ref)} ({Lr} bp) -> {kept} columns; align {t_aln:.1f} s, tree {t_tree:.0f} s, {model}",
          file=sys.stderr, flush=True)

print("## Q — reference-projected alignment trees (SH-aLRT > 75)")
agree = collections.Counter()
summary = collections.defaultdict(collections.Counter)
for tag, fam, ref, Lr, kept, t_aln, t_tree, model, q, p1, nleaves, ntruth in out:
    res, positional, part, mid = q
    ps, pp, br, bp, ex, nt = part
    groups = "; ".join(f"{g}: {v[0]}" + (f" (SH {v[1]:.1f}/UFB {v[2]:.0f})" if v[1] is not None else "") for g, v in res.items())
    extra = "" if positional is None else f" | positional split supported: {positional} -> {'WRONG' if positional else 'CORRECT'}"
    print(f"{tag}: {nleaves} leaves ({ntruth} truth), reference {ref} {Lr} bp -> {kept} columns, {model}, "
          f"{t_aln:.1f}+{t_tree:.0f} s | {groups}{extra}")
    print(f"      midpoint-root L1 partition: pairwise sens {ps:.3f} prec {pp:.3f} | bipartite R {br:.3f} P {bp:.3f} exact {ex}/{nt}")
    print(f"      root sides: {' '.join(mid[0])}  ||  {' '.join(mid[1])}")
    kind = "reference" if tag.startswith("ref") else "leave-out"
    for g, v in res.items():
        summary[(fam, kind, g)][v[0]] += 1
        if p1 is not None and v[0] != "absent":
            agree["agree" if (v[0] == "RECOVERED") == (p1[0][g][0] == "RECOVERED") else "disagree"] += 1
    if positional is not None:
        summary[(fam, kind, "positional split")]["WRONG" if positional else "CORRECT"] += 1
print("\n## Q summary (counts over runs)")
for k in sorted(summary):
    print(f"{k[0]:6s} {k[1]:9s} {k[2]:28s} {dict(summary[k])}")
print(f"\n## agreement with P1 (MAFFT) on RECOVERED calls, over (run, group) cells: {dict(agree)}")


def _ms(v):
    v = [x for x in v if not (isinstance(x, float) and math.isnan(x))]
    return "nan" if not v else (f"{statistics.mean(v):.3f}±{statistics.stdev(v):.3f}" if len(v) > 1 else f"{v[0]:.3f}")


print("\n## P2 — hybrid width on the G1 candidates (median over replicates' medians; counts mean±SD)")
for level in ("half", "keep1"):
    for fam in families + ["ALL"]:
        for scope in ("all", "like-for-like"):
            for arm in ("M0", "G1-clip", "H", "H-union"):
                R = [r for r in width_rows if r["level"] == level and r["family"] == fam and r["arm"] == arm and r["scope"] == scope]
                g = lambda k: _ms([r[k] for r in R])
                print(f"{level:5s} {fam:6s} {scope:13s} {arm:8s} n {g('n')} tx-boundary {g('tx_frac')} | J {g('jaccard')} "
                      f"5' {g('off5')} 3' {g('off3')} trunc {g('truncated')} overext {g('overextended')}")
