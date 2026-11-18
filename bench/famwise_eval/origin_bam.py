"""PREREG r1794 (docs/PREREG_famwise_origin_2026-11-18.md), row b: the assembly inputs of the oracle comparison, from one
alignment.

  origin_bam.py --truth READS.tsv [--tie-split T.tsv] --bam SEL.bam --out DIR [--promoted-mapq Q]

DIR/D.bam  a read with a truth entry (kinds one, several, tie_unsplit) and a record at its origin keeps every record;
           any other read keeps its primary record only.
DIR/O.bam  the same, but a read of the first kind keeps only the records on its copy or copies of origin; if its primary
           is not among them, the kept record with the highest ms becomes the primary (flag 0x100 cleared).
DIR/regions.bed  merged spans of the origin places, and of the primary's place, of every read with a record at its origin
                 that is tied (kinds several, tie_unsplit) or whose primary is off its origin.
DIR/inputs.tsv   reads per case, and what the default seeding band (AS >= 0.98 x best) leaves out of D's exact ms ties.
The three inputs: P = D.bam, primary alignments only (assemble --seed-pool primary; StringTie: samtools view -b -F 256);
D = D.bam with the defaults; O = O.bam with every record (assemble --seed-pool all; StringTie: O.bam)."""
import argparse, collections, os, sys
import pysam
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import records as R

ORACLE = ("one", "several", "tie_unsplit")

def plan(recs, t):
    """-> (case, kept ordinals for D, kept ordinals for O, ordinal to promote in O or None)"""
    prim = [i for i, r in enumerate(recs) if r.primary] or [0]
    if t is None or t.kind not in ORACLE:
        return "primary_only", prim, prim, None
    home = [i for i, r in enumerate(recs) if R.locus_index(r, t) in t.origin]
    if not home:
        return "no_home", prim, prim, None
    every = list(range(len(recs)))
    if prim[0] not in home:
        return "moved", every, home, max(home, key=lambda i: recs[i].ms)
    return ("kept_as_is" if home == every else "restricted"), every, home, None

def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--truth", required=True); ap.add_argument("--tie-split"); ap.add_argument("--bam", required=True)
    ap.add_argument("--out", required=True); ap.add_argument("--promoted-mapq", type=int, help="MAPQ written on a promoted record (default: unchanged)")
    a = ap.parse_args()
    os.makedirs(a.out, exist_ok=True)
    truth, recs = R.load_truth(a.truth, a.tie_split), R.load_records(a.bam)
    plans, cases, spans, below = {}, collections.Counter(), collections.defaultdict(list), 0
    for q, rs in recs.items():
        t = truth.get(q)
        plans[q] = case, d, o, _ = plan(rs, t)
        cases[case] += 1
        if case == "moved" or (case in ("kept_as_is", "restricted") and t.kind != "one"):    # tied, or primary off its origin; has a record at home
            for i in t.origin:
                spans[t.loci[i][0]].append(t.loci[i][1:3])
            p = next((r for r in rs if r.primary), rs[0])
            i = R.locus_index(p, t)
            spans[p.chrom].append(t.loci[i][1:3] if i is not None else (p.start, p.end))
        if len(d) > 1:
            bm, ba = max(r.ms for r in rs), max(r.AS for r in rs)
            below += sum(not r.primary and r.ms == bm and r.AS < 0.98 * ba for r in rs)
    tmp = os.environ.get("TMPDIR", a.out)
    with pysam.AlignmentFile(a.bam, check_sq=False) as fh:
        outs = {x: pysam.AlignmentFile(f"{tmp}/origin_bam.{os.getpid()}.{x}.bam", "wb", template=fh) for x in "DO"}
        seen = collections.Counter()
        for r in fh.fetch(until_eof=True):
            if r.is_unmapped or r.is_supplementary:
                continue
            i = seen[r.query_name]; seen[r.query_name] += 1
            _, d, o, promote = plans[r.query_name]
            if i in d:
                outs["D"].write(r)
            if i in o:
                if i == promote:
                    r.flag &= ~0x100
                    if a.promoted_mapq is not None:
                        r.mapping_quality = a.promoted_mapq
                outs["O"].write(r)
        for x, w in outs.items():
            w.close()
            pysam.sort("-o", f"{a.out}/{x}.bam", w.filename.decode()); pysam.index(f"{a.out}/{x}.bam"); os.remove(w.filename)
    with open(f"{a.out}/regions.bed", "w") as fh:
        for c in sorted(spans):
            cur = None
            for s, e in sorted(spans[c]):
                if cur and s <= cur[1]:
                    cur[1] = max(cur[1], e)
                else:
                    if cur:
                        fh.write(f"{c}\t{cur[0]}\t{cur[1]}\n")
                    cur = [s, e]
            fh.write(f"{c}\t{cur[0]}\t{cur[1]}\n")
    with open(f"{a.out}/inputs.tsv", "w") as fh:
        fh.write("# reads per case; promoted records keep their MAPQ" + ("" if a.promoted_mapq is None else f" -> {a.promoted_mapq}") + "\n")
        fh.writelines(f"{c}\t{cases[c]}\n" for c in ("kept_as_is", "restricted", "moved", "no_home", "primary_only"))
        fh.write(f"D_records_tied_on_ms_below_the_AS_band\t{below}\n")
    print(open(f"{a.out}/inputs.tsv").read(), end="")

if __name__ == "__main__":
    main()
