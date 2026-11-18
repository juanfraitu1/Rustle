#!/usr/bin/env python3
"""Tool table of the cross-copy splicing design (docs/superpowers/specs/2026-11-18-cross-copy-splicing-design.md, section 5):
for COPY-SPECIFIC calls, which arms emit the junction at the copy that uses it (U), and at the one that does not (S)?

    tool_table.py [--dev GENOME=DEV.sites.tsv] GENOME,FAM.psi.tsv,SITES.tsv,ARM=GTF[,ARM=GTF ...] [...]      one argument per sample

SITES = site_cause.py's table of the sample (DEV.sites.tsv: of the development sample): the units, the MISKEYED calls
(left out) and the verdict at tau 1 for the split. A unit (genome, U's junction, S's anchor site) counts once: where two
samples of a genome call it, the one with more certain reads at U speaks, among the samples that have the arm (for a
pair of arms: that have both), and the unit counts when it is COPY-SPECIFIC in that sample. Units of the development
sample count nowhere. An arm is "behind" another when its share at U is more than 0.05 lower on >= 30 shared units,
else "level". Descriptive: no bar."""
import collections, itertools, os, sys
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path[:0] = [HERE, os.path.join(HERE, "..", "..", "tools")]
import calls as K
from psi_bunches import gtf_introns

BEHIND, MIN_UNITS = 0.05, 30

class Table:
    def __init__(self, dev, specs):
        read = lambda path: [dict(zip(ls[0].split("\t"), l.split("\t"))) for ls in [open(path).read().splitlines()] for l in ls[1:]]
        self.seen = set()
        if dev:
            g, path = dev.split("=", 1)
            self.seen = {(g, r["unit"]) for r in read(path)}
        self.samples = []
        for spec in specs:
            g, psi, sites, *arms = spec.split(",")
            verdict = {(r["family"], r["anchor"], r["col_anchor"], r["col_other"]): r["verdict_1"] for r in read(sites)}
            calls = [c for c in K.load_calls(psi) if verdict[c[:4]] != "MISKEYED" and (g, K.unit_name(c)) not in self.seen]
            self.samples.append((g, calls, {a.split("=", 1)[0]: set(gtf_introns(a.split("=", 1)[1])) for a in arms}, verdict))
        self.arms = sorted({a for s in self.samples for a in s[2]})

    def units(self, *arms):
        """{unit: (call, the sample's arms, its A verdict)}: the copy-specific units over the samples that have every arm of `arms`"""
        out = {}
        for g, calls, have, verdict in self.samples:
            if all(a in have for a in arms):
                for c in calls:
                    u = (g, K.unit_name(c))
                    if u not in out or c.reads_u > out[u][0].reads_u:
                        out[u] = (c, have, verdict[c[:4]])
        return {u: v for u, v in out.items() if v[0].copy_specific}

    def shares(self, arm, among=None):
        """(units, share with the junction emitted at U, share with its image emitted at S)"""
        U = [v for v in self.units(*(among or [arm])).values()]
        n = len(U)
        return (n, sum(c.u_intron in have[arm] for c, have, _ in U) / n, sum(c.s_intron in have[arm] for c, have, _ in U) / n) if n else (0, 0.0, 0.0)

    __getitem__ = shares

    def reading(self, arm, other):
        n, a, _ = self.shares(arm, [arm, other]); _, b, _ = self.shares(other, [arm, other])
        return "behind" if n >= MIN_UNITS and round(b - a, 6) > BEHIND else "level"

table = Table

def main(argv):
    dev = argv.pop(1) if argv and argv[0] == "--dev" else None
    t = Table(dev, argv[1:] if dev else argv)
    print(f"samples {len(t.samples)}; copy-specific units {len(t.units())} ({len(t.seen)} calls of the development sample left out wherever they recur)")
    print("arm\tsamples\tunits\temitted at U\temitted at S")
    for a in t.arms:
        n, u, s = t[a]
        print(f"{a}\t{sum(a in x[2] for x in t.samples)}\t{n}\t{round(u * n)} = {u:.3f}\t{round(s * n)} = {s:.3f}")
    print("pairs (arm, other, shared units, share at U of each, reading of arm against other)")
    for a, b in itertools.permutations(t.arms, 2):
        n, x, _ = t.shares(a, [a, b]); _, y, _ = t.shares(b, [a, b])
        print(f"{a}\t{b}\t{n}\t{x:.3f}\t{y:.3f}\t{t.reading(a, b)}")
    print("by the sequence test's verdict (arm, verdict, units, emitted at U, emitted at S)")
    for a in t.arms:
        by = collections.defaultdict(list)
        for c, have, v in t.units(a).values():
            by[v].append((c.u_intron in have[a], c.s_intron in have[a]))
        for v, xs in sorted(by.items()):
            print(f"{a}\t{v}\t{len(xs)}\t{sum(x for x, _ in xs)}\t{sum(y for _, y in xs)}")

if __name__ == "__main__":
    main(sys.argv[1:])
