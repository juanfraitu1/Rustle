#!/usr/bin/env python3
"""r1696 step 0 (docs/PREREG_bubble_vote_2026-11-13.md): per-bubble vote on row 1649's family graphs.
A locus (copies with overlapping spans) carries a node when any of its paths holds it; an informative node of a read's
walk is carried by some but not all loci; assign to the locus with the fewest contradictions iff the runner-up has
>= MARGIN more.
shortcut: a node visit is a vote, no CIGAR or base-quality check; add one before any registered row.
Usage: bubble_vote.py pick GRAPH.gfa ALN.gaf COPIES.tsv MARGIN > picks.tsv
       bubble_vote.py score PICKS KMER_PICKS COPIES.tsv DUMP [DUMP...]"""
import collections, csv, re, sys

def spans(copies):
    return {(r["family_id"], r["copy_idx"]): (r["chrom"], int(r["start"]), int(r["end"])) for r in csv.DictReader(open(copies), delimiter="\t")}

def ov(span, a, b):
    return span[a][0] == span[b][0] and span[a][1] < span[b][2] and span[b][1] < span[a][2]

def pick(gfa, gaf, copies, margin, out=sys.stdout):
    span = spans(copies); loci = []; carry = collections.defaultdict(set)  # node -> locus indices
    for l in open(gfa):
        if l[0] != "P": continue
        f = l.rstrip("\n").split("\t"); k = tuple(f[1].split("|")[:2])
        i = next((i for i, m in enumerate(loci) if any(ov(span, k, x) for x in m)), None)
        if i is None: loci.append({k}); i = len(loci) - 1
        else: loci[i].add(k)
        for x in f[2].split(","): carry[x[:-1]].add(i)
    best = {}
    for l in open(gaf):
        f = l.rstrip("\n").split("\t"); q = f[0].split()[0]; q = q[2:] if q.startswith("R|") else q
        cov = (int(f[3]) - int(f[2])) / int(f[1]); idn = next((float(t[5:]) for t in f[12:] if t.startswith("id:f:")), 0.0)
        if q not in best or (cov, idn) > best[q][:2]: best[q] = (cov, idn, set(re.findall(r"[<>]([^<>]+)", f[5])))
    n = collections.Counter(aligned=len(best))
    for q, (cov, _, nodes) in best.items():
        if cov < 0.90: n["abstain_cover"] += 1; continue
        info = [carry[v] for v in nodes if len(carry[v]) < len(loci)]
        contra = sorted((sum(i not in c for c in info), i) for i in range(len(loci)))
        if len(contra) > 1 and contra[1][0] - contra[0][0] < margin: n["abstain_margin"] += 1; continue
        k = sorted(loci[contra[0][1]])[0]; n["assigned"] += 1; print(q, k[0], k[1], sep="\t", file=out)
    print(dict(n), file=sys.stderr)

def score(picks, kmer, copies, *dumps):
    span = spans(copies); rd = lambda p: {l.split("\t")[0]: tuple(l.rstrip("\n").split("\t")[1:3]) for l in open(p)}
    pk, km = rd(picks), rd(kmer); mol = {}
    for d in dumps:
        for l in open(d):
            q, fam, st, v, al, _, t = l.rstrip("\n").split("\t")
            m = mol.setdefault(q, dict(al=al, o2=False, T=None if t == "OUTSIDE" else {tuple(x.split(":")) for x in t.split(";")}))
            m["o2"] |= st == "assigned"
    fams = {f for f, _ in pk.values()}; c = collections.Counter(); per = collections.defaultdict(collections.Counter)
    right = lambda p, q: any(p == t or ov(span, p, t) for t in mol[q]["T"])
    # the panel = molecules whose truth or pick is in the graphs' families (row 1649 built per-family read sets)
    for q, m in mol.items():
        if q not in pk and not (m["T"] and {t[0] for t in m["T"]} & fams): continue
        if not m["T"]:
            c["negative assigned"] += q in pk; continue
        c["decidable"] += 1; k = q in km; g = q in pk
        c["kmer assigned"] += k; c["kmer right"] += k and right(km[q], q)
        if not g: continue
        r = right(pk[q], q); c["assigned"] += 1; c["right" if r else "wrong"] += 1; per[pk[q][0]]["right" if r else "wrong"] += 1
        if not k and not m["o2"]:
            c["only (kmer and O2 abstain)"] += 1; c["only right" if r else "only wrong"] += 1
            c["only: primary " + ("right" if m["al"] == "al_right" else "wrong" if m["al"] == "al_wrong" else "NA")] += 1
            c["only: vote right, primary wrong"] += r and m["al"] != "al_right"; c["only: vote wrong, primary right"] += (not r) and m["al"] == "al_right"
    for k in sorted(c): print(k, c[k], sep="\t")
    print({f: dict(v) for f, v in per.items()})

def _selfcheck():
    import tempfile, io
    d = tempfile.mkdtemp(); open(d + "/c.tsv", "w").write("family_id\tcopy_idx\tchrom\tstart\tend\nF\t0\tc\t0\t10\nF\t1\tc\t100\t110\n")
    open(d + "/g.gfa", "w").write("P\tF|0|cat|0\t1+,2+,4+,5+,7+\t*\nP\tF|1|cat|0\t1+,3+,4+,6+,7+\t*\n")
    g = lambda w: f"\t100\t0\t100\t+\t{w}\t9\t0\t9\t9\t9\t60\tid:f:1\n"
    # r1: two copy-0 alleles; r2: one each (a containment rule finds no path, the vote ties); r3: one vote only
    open(d + "/a.gaf", "w").write("r1" + g(">1>2>4>5>7") + "r2" + g(">1>2>4>6>7") + "r3" + g(">1>3>4"))
    o = io.StringIO(); e, sys.stderr = sys.stderr, io.StringIO()
    try: pick(d + "/g.gfa", d + "/a.gaf", d + "/c.tsv", 2, o)
    finally: sys.stderr = e
    assert o.getvalue() == "r1\tF\t0\n", o.getvalue()

if __name__ == "__main__":
    _selfcheck()
    if sys.argv[1] == "pick": pick(*sys.argv[2:5], int(sys.argv[5]))
    else: score(*sys.argv[2:])
