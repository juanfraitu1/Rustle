#!/usr/bin/env python3
"""copy_phasing (correlation clustering, Strainy / SDA lineage) against HG002's own-haplotype truth on CHM13.

  cp_truth.py inputs --bam SEL.bam --out PREFIX [--species human] [--pad 2000]
      PREFIX.copies.tsv  the truth families (famtruth recent paralogue components, >= 2 genes) as copy_phasing's copies
      PREFIX.net.tsv     every read with a non-supplementary record on a copy of the family (real_net.py's convention)
  cp_truth.py score --truth READS.tsv --tie-split T.tsv --bam SEL.bam --copies PREFIX.copies.tsv --run RUN_PREFIX [--out FILE]
      O2: for molecules the truth decides (kind one) with a copy_phasing assignment, assigned copy vs origin, beside minimap2's
          primary on the same molecules; split by near-tied on ms (two best places within 10) or clear.
      O3: families with >= 6 `outside` molecules (origin at an HG002 place no CHM13 locus of the read represents: the copy may be
          one CHM13 lacks) against families copy_phasing flags (certified candidates or putative copies)."""
import argparse, collections, csv, os, sys
import pysam
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE); sys.path.insert(0, os.path.join(HERE, "..", "denovo_eval"))
import famtruth as F
import records as R
from member_completeness import families

def cmd_inputs(a):
    G = F.load_genes(a.species)
    fams, _ = families(G, a.species)
    with open(a.out + ".copies.tsv", "w") as fh:
        fh.write("family_id\tcopy_idx\ttid\tchrom\tstart\tend\tstrand\texons\n")
        for i, f in enumerate(fams):
            for j, g in enumerate(sorted(f)):
                ex = sorted(G[g]["exons"])
                fh.write(f"F{i}\t{j}\t{g}\t{G[g]['chrom']}\t{ex[0][0]}\t{ex[-1][1]}\t{G[g]['strand']}\t" + ",".join(f"{s}-{e}" for s, e in ex) + "\n")
    n = 0
    with pysam.AlignmentFile(a.bam) as b, open(a.out + ".net.tsv", "w") as w:
        w.write("family\tread\tused\tvia\n")
        for i, f in enumerate(fams):
            reads = set()
            for g in f:
                ex = G[g]["exons"]
                if G[g]["chrom"] not in b.references:
                    continue
                for r in b.fetch(G[g]["chrom"], max(0, min(s for s, _ in ex) - a.pad), max(e for _, e in ex) + a.pad):
                    if not r.is_unmapped and not r.is_supplementary:
                        reads.add(r.query_name)
            for q in sorted(reads):
                w.write(f"F{i}\t{q}\t1\tB\n")
            n += len(reads)
    print(f"{len(fams)} families, {sum(len(f) for f in fams)} copies, {n} net rows")

def cmd_score(a):
    T = R.load_truth(a.truth, a.tie_split)
    recs = R.load_records(a.bam)
    cp = {}
    with open(a.copies) as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            cp[r["tid"]] = (r["chrom"], int(r["start"]), int(r["end"]), r["strand"], r["family_id"])
    def at_copy(locus, g):
        c, s, e, st = locus[:4]
        x = cp.get(g)
        return x is not None and x[0] == c and x[1] < e and s < x[2]
    def primary_locus(q, t):
        p = next((r for r in recs.get(q, ()) if r.primary), None)
        return R.locus_index(p, t) if p else None
    o2 = collections.Counter()
    with open(a.run + ".o2.tsv") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            q, g = r["read"], r["assigned"]
            t = T.get(q)
            if t is None or t.kind != "one" or q not in recs:
                continue
            o = next(iter(t.origin))
            near = "near_tied" if R.near_tied(recs[q], "ms") else "clear"
            pl = primary_locus(q, t)
            pr_right = pl == o
            kind = ("new" if g.startswith("new:") else "abstain" if g not in cp       # ABSTAIN, HOPELESS (exact ties it refuses)
                    else "right" if at_copy(t.loci[o], g) else "wrong")
            o2[(near, "n")] += 1
            o2[(near, kind)] += 1
            o2[(near, "primary_right")] += pr_right
            o2[(near, "primary_right_on_placed")] += pr_right and kind in ("right", "wrong")
            o2[(near, "fixes")] += kind == "right" and not pr_right
            o2[(near, "breaks")] += kind == "wrong" and pr_right
    L = ["O2 (truth decides the molecule; copy_phasing assigned it):"]
    for near in ("near_tied", "clear"):
        n = o2[(near, "n")]
        placed = o2[(near, "right")] + o2[(near, "wrong")]
        L.append(f"  {near}: n {n}; copy_phasing right {o2[(near, 'right')]} / wrong {o2[(near, 'wrong')]} / abstain {o2[(near, 'abstain')]} / new copy {o2[(near, 'new')]}; "
                 f"primary on the {placed} placed: right {o2[(near, 'primary_right_on_placed')]}; primary on all {n}: right {o2[(near, 'primary_right')]}; "
                 f"fixes {o2[(near, 'fixes')]}, breaks {o2[(near, 'breaks')]}")
    # O3
    fam_of_read = collections.defaultdict(set)
    with open(a.run.rsplit(".", 0)[0] + ".reads.tsv") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            fam_of_read[r["read"]].add(r["family"])
    outside = collections.Counter()
    for q, t in T.items():
        if t.kind != "outside":
            continue
        for f in fam_of_read.get(q, ()):
            outside[f] += 1
    flagged = set()
    with open(a.run + ".candidates.tsv") as fh:
        flagged |= {r["family"] for r in csv.DictReader(fh, delimiter="\t") if r["flagged"] == "1"}
    putative = set()
    with open(a.run + ".putative.tsv") as fh:
        putative |= {r["family"] for r in csv.DictReader(fh, delimiter="\t")}
    with open(a.run + ".families.tsv") as fh:
        ran = {r["family"] for r in csv.DictReader(fh, delimiter="\t")}
    pos = {f for f in ran if outside[f] >= 6}
    neg = {f for f in ran if outside[f] == 0}
    L.append(f"O3 (families copy_phasing ran: {len(ran)}; with >= 6 outside molecules: {len(pos)}; with none: {len(neg)}):")
    for name, S in (("certified", flagged), ("putative", putative), ("either", flagged | putative)):
        L.append(f"  {name}: flags {len(S & ran)}; in >= 6-outside families {len(S & pos)}/{len(pos)}; in no-outside families {len(S & neg)}/{len(neg)}")
    out = "\n".join(L) + "\n"
    if a.out:
        open(a.out, "w").write(out)
    print(out, end="")

def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sp = ap.add_subparsers(dest="cmd", required=True)
    p = sp.add_parser("inputs"); p.add_argument("--bam", required=True); p.add_argument("--out", required=True)
    p.add_argument("--species", default="human"); p.add_argument("--pad", type=int, default=2000)
    p = sp.add_parser("score"); p.add_argument("--truth", required=True); p.add_argument("--tie-split", required=True)
    p.add_argument("--bam", required=True); p.add_argument("--copies", required=True); p.add_argument("--run", required=True); p.add_argument("--out")
    a = ap.parse_args(argv)
    dict(inputs=cmd_inputs, score=cmd_score)[a.cmd](a)

if __name__ == "__main__":
    main()
