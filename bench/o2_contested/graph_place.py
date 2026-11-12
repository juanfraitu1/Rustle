#!/usr/bin/env python3
"""r1649 step 0 (docs/STEP0_graph_placement_2026-11-12.md): picks from GraphAligner GAF on a seqwish family graph whose
paths are the input sequences (FAMILY|COPY|KIND|ID). Best alignment per read (longest aligned query span, then
identity), >= 0.90 of the read aligned; the read is compatible with every path that contains all its nodes; assign iff
the compatible sequences are copies of one locus.
shortcut: node-set containment ignores order and orientation along the path; upgrade to walk matching before any row.
Usage: graph_place.py GRAPH.gfa ALN.gaf COPIES.tsv > OUT.picks.tsv"""
import collections, csv, re, sys

def main(gfa, gaf, copies):
    span = {}
    for r in csv.DictReader(open(copies), delimiter="\t"): span[(r["family_id"], r["copy_idx"])] = (r["chrom"], int(r["start"]), int(r["end"]))
    paths = {}
    for l in open(gfa):
        if l[0] == "P":
            f = l.rstrip("\n").split("\t"); paths[f[1]] = {x[:-1] for x in f[2].split(",")}
    best = {}
    for l in open(gaf):
        f = l.rstrip("\n").split("\t"); q = f[0].split()[0]; q = q[2:] if q.startswith("R|") else q
        cov = (int(f[3]) - int(f[2])) / int(f[1]); idn = next((float(t[5:]) for t in f[12:] if t.startswith("id:f:")), 0.0)
        if q not in best or (cov, idn) > best[q][:2]: best[q] = (cov, idn, set(re.findall(r"[<>]([^<>]+)", f[5])))
    same = lambda a, b: span[a][0] == span[b][0] and span[a][1] < span[b][2] and span[b][1] < span[a][2]
    n = collections.Counter(aligned=len(best))
    for q, (cov, _, nodes) in best.items():
        if cov < 0.90: n["abstain_cover"] += 1; continue
        ks = sorted({tuple(p.split("|")[:2]) for p, ns in paths.items() if nodes <= ns})
        if not ks: n["abstain_no_path"] += 1; continue
        if not all(same(ks[0], k) for k in ks[1:]): n["abstain_multi_locus"] += 1; continue
        n["assigned"] += 1; print(q, ks[0][0], ks[0][1], sep="\t")
    print(dict(n), file=sys.stderr)

def _selfcheck():
    import tempfile, os
    d = tempfile.mkdtemp(); open(d + "/c.tsv", "w").write("family_id\tcopy_idx\tchrom\tstart\tend\nF\t0\tc\t0\t10\nF\t1\tc\t100\t110\n")
    open(d + "/g.gfa", "w").write("P\tF|0|cat|0\t1+,2+,4+\t*\nP\tF|1|cat|0\t1+,3+,4+\t*\n")
    open(d + "/a.gaf", "w").write("r1\t100\t0\t100\t+\t>1>2>4\t9\t0\t9\t9\t9\t60\tid:f:1\nr2\t100\t0\t100\t+\t>1>4\t9\t0\t9\t9\t9\t60\tid:f:1\nr3\t100\t0\t50\t+\t>3\t9\t0\t9\t9\t9\t60\tid:f:1\n")
    import io, contextlib; o = io.StringIO(); e = io.StringIO()
    with contextlib.redirect_stdout(o), contextlib.redirect_stderr(e): main(d + "/g.gfa", d + "/a.gaf", d + "/c.tsv")
    assert o.getvalue() == "r1\tF\t0\n", o.getvalue()  # r2 fits both copies, r3 is under 0.90

if __name__ == "__main__":
    _selfcheck(); main(*sys.argv[1:4])
