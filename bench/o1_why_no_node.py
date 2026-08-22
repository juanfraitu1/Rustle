#!/usr/bin/env python3
"""WHY did the expressed, duplication-supported locations not become nodes?

CONTEXT. SDs predict 3,928 copy locations the catalog has nothing at; ~35% of them carry >= 3 primary
reads, i.e. ~1,374 EXPRESSED candidates that O2 never gets to consider. Is that a design choice the
node builder makes correctly, or a gap?

THE NODE BUILDER'S OWN CRITERIA are the yardstick, not my opinion of them:
  * the shipped catalog has 0/1415 SINGLE-EXON nodes -> an unspliced read pile is correctly rejected;
  * MIN_READS = 3 applies to the supporting chain, not to the interval;
  * the mis-chain filter drops a giant intron (> 50 kb) supported by < 3 reads.

SO EACH EXPRESSED-BUT-ABSENT INTERVAL FALLS INTO ONE OF:
  UNSPLICED           no read carries an N -> correctly not a node
  NO CHAIN >= 3       spliced, but no single intron chain reaches MIN_READS -> correctly not a node
  ⭐ SHOULD BE A NODE  a consistent intron chain with >= 3 reads -> a genuine RECALL GAP

⚠ `-F 2308` before any per-read CIGAR statistic — the standing invariant.
"""
import bisect, collections, csv, re, subprocess, sys

G   = "/mnt/linuxdisk/home/juanfraitu/o1_gw"
D   = "/mnt/linuxdisk/home/juanfraitu/winloci_data"
SD  = "/mnt/c/Users/jfris/Desktop/final.bed"
BAM = f"{D}/GGO_ds.bam"
MIN_ID, MIN_UNIT, MIN_READS = 0.90, 1000, 3
CIG = re.compile(r"(\d+)([MIDNSHP=X])")
SAMPLE = int(sys.argv[1]) if len(sys.argv) > 1 else 250


def main():
    loc = {}
    for r in csv.DictReader(open(f"{G}/ggo_gw.copies.tsv"), delimiter="\t"):
        loc[f"{r['family_id']}~{r['copy_idx']}"] = (r["chrom"], int(r["start"]), int(r["end"]))
    by = collections.defaultdict(list)
    for k, (c, s, e) in loc.items():
        by[c].append((s, e, k))
    for c in by:
        by[c].sort()
    st = {c: [x[0] for x in v] for c, v in by.items()}

    def any_hit(c, s, e):
        v = by.get(c)
        if not v:
            return False
        i = bisect.bisect_left(st[c], e)
        return any(v[j][1] > s for j in range(max(0, i - 500), i))

    predicted = set()
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
            if any_hit(ca, sa, ea) and not any_hit(cb, sb, eb):
                predicted.add((cb, sb, eb))

    keys = sorted(predicted)
    step = max(1, len(keys) // SAMPLE)
    samp = keys[::step][:SAMPLE]
    print(f"SD-predicted, catalog-absent intervals: {len(predicted)}; sampling {len(samp)}", flush=True)

    cls = collections.Counter(); gaps = []
    for (c, s, e) in samp:
        p = subprocess.run(["samtools", "view", "-F", "2308", BAM, f"{c}:{s+1}-{e}"],
                           capture_output=True, text=True)
        chains = collections.Counter(); n_primary = n_spliced = 0
        for line in p.stdout.splitlines():
            f = line.split("\t")
            if len(f) < 6:
                continue
            n_primary += 1
            pos = int(f[3]); introns = []
            for num, op in CIG.findall(f[5]):
                num = int(num)
                if op in "MD=XN":
                    if op == "N":
                        introns.append((pos, pos + num))
                    pos += num
            if introns:
                n_spliced += 1
                chains[tuple(introns)] += 1
        if n_primary < MIN_READS:
            cls["(below MIN_READS overall)"] += 1
        elif n_spliced == 0:
            cls["UNSPLICED — correctly not a node"] += 1
        elif not chains or max(chains.values()) < MIN_READS:
            cls["NO CHAIN >= 3 — correctly not a node"] += 1
        else:
            cls["⭐ SHOULD BE A NODE — recall gap"] += 1
            gaps.append((c, s, e, n_primary, max(chains.values()), len(next(iter(
                sorted(chains.items(), key=lambda x: -x[1])))[0])))

    print("\n=== why is there no node here? (unit = SD-predicted absent interval) ===")
    n = len(samp)
    for k, v in cls.most_common():
        print(f"  {k:<38} {v:>4}/{n}  = {v/n:.4f}")
    g = cls["⭐ SHOULD BE A NODE — recall gap"]
    print(f"\n  ⟹ scaled to all {len(predicted)}: ~{int(len(predicted)*g/n)} intervals carry a"
          f" >= {MIN_READS}-read spliced chain and still have no node.")
    print("\n  examples (chrom, start, end, primaries, best-chain reads, introns in that chain):")
    for x in gaps[:10]:
        print(f"    {x[0]}:{x[1]}-{x[2]}  primaries={x[3]}  chain_reads={x[4]}  introns={x[5]}")


if __name__ == "__main__":
    main()
