#!/usr/bin/env python3
"""Does O1 lay down ALL the possibilities O2 needs?

THE ASYMMETRY. O2 assigns a read to a copy. A copy O1 MISSES can never receive a read — the
possibility is not on the table. A copy O1 wrongly ADDS, O2 can decline. So for the thesis's core
objective, O1's RECALL matters more than its precision — and essentially all effort has gone to
precision (seven closed repair routes on the coverage clause).

THE MEASUREMENT. Segmental duplications predict where copies should be. For every catalog copy, take
its SD partner interval and ask whether the catalog has anything there. If not, that is a possibility
O2 never sees. Then split those by whether they are RECOVERABLE:

  partner interval has READS      -> EXPRESSED, so a node COULD have been built  -> a genuine miss
  partner interval has ~no reads  -> unexpressed, structurally invisible to RNA  -> not recoverable

⚠ An SD partner need not contain a gene at all, so this is an UPPER bound on missed copies; the read
split is what makes it interpretable. ⚠ Reads are counted `-F 2308` (primary only), the standing
invariant before any per-read statistic.
"""
import bisect, collections, csv, subprocess, sys

G   = "/mnt/linuxdisk/home/juanfraitu/o1_gw"
D   = "/mnt/linuxdisk/home/juanfraitu/winloci_data"
SD  = "/mnt/c/Users/jfris/Desktop/final.bed"
BAM = f"{D}/GGO_ds.bam"
MIN_ID, MIN_UNIT = 0.90, 1000
SAMPLE = int(sys.argv[1]) if len(sys.argv) > 1 else 400


def main():
    loc, fam = {}, {}
    for r in csv.DictReader(open(f"{G}/ggo_gw.copies.tsv"), delimiter="\t"):
        k = f"{r['family_id']}~{r['copy_idx']}"
        loc[k] = (r["chrom"], int(r["start"]), int(r["end"])); fam[k] = r["family_id"]
    by = collections.defaultdict(list)
    for k, (c, s, e) in loc.items():
        by[c].append((s, e, k))
    for c in by:
        by[c].sort()
    st = {c: [x[0] for x in v] for c, v in by.items()}

    def covered(c, s, e):
        v = by.get(c)
        if not v:
            return False
        i = bisect.bisect_left(st[c], e)
        for j in range(max(0, i - 500), i):
            s0, e0, _ = v[j]
            if e0 > s:
                return True
        return False

    def hits(c, s, e):
        v = by.get(c)
        if not v:
            return []
        i = bisect.bisect_left(st[c], e); out = []
        for j in range(max(0, i - 500), i):
            s0, e0, k = v[j]
            if e0 > s:
                out.append(k)
        return out

    predicted = {}        # partner interval -> the catalog copy that predicts it
    for line in open(SD):
        f = line.rstrip("\n").split("\t")
        if len(f) < 34:
            continue
        try:
            idn = float(f[-1])
        except ValueError:
            continue
        a0, a1, b0, b1 = int(f[1]), int(f[2]), int(f[4]), int(f[5])
        if idn < MIN_ID or a1 - a0 < MIN_UNIT or b1 - b0 < MIN_UNIT:
            continue
        for (ca, sa, ea, cb, sb, eb) in ((f[0], a0, a1, f[3], b0, b1), (f[3], b0, b1, f[0], a0, a1)):
            src = hits(ca, sa, ea)
            if src and not covered(cb, sb, eb):
                predicted.setdefault((cb, sb, eb), src[0])

    print(f"catalog: {len(loc)} copies in {len(set(fam.values()))} families")
    print(f"SD-PREDICTED copies the catalog has NOTHING at: {len(predicted)}", flush=True)

    keys = sorted(predicted)
    step = max(1, len(keys) // SAMPLE)
    samp = keys[::step][:SAMPLE]
    print(f"sampling {len(samp)} of them for read support (-F 2308) ...", flush=True)

    expressed = collections.Counter(); depths = []
    for (c, s, e) in samp:
        p = subprocess.run(["samtools", "view", "-c", "-F", "2308", BAM, f"{c}:{s+1}-{e}"],
                           capture_output=True, text=True)
        try:
            n = int(p.stdout.strip())
        except ValueError:
            n = 0
        depths.append(n)
        expressed["reads >= 3 (a node COULD be built)" if n >= 3 else
                   "reads 1-2" if n >= 1 else "NO reads (structurally invisible)"] += 1

    print("\n=== are the missed possibilities RECOVERABLE? ===")
    n = len(samp)
    for k, v in expressed.most_common():
        print(f"  {k:<40} {v:>5}/{n}  = {v/n:.4f}")
    import statistics as stt
    print(f"\n  median primary reads in a predicted-but-absent interval: {stt.median(depths):.0f}")
    rec = expressed["reads >= 3 (a node COULD be built)"]
    print(f"\n  ⟹ scaled to all {len(predicted)} predicted-absent intervals: "
          f"~{int(len(predicted)*rec/n)} are EXPRESSED and were still not laid down.")


if __name__ == "__main__":
    main()
