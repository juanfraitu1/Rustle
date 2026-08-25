#!/usr/bin/env python3
"""What margin should the read-strand vote require before it is allowed to overrule the placeholder?

WHY. RUSTLE_READ_STRAND calls strand by a STRICT MAJORITY of read orientation, measured at 0.9650 against
junction-determined truth. Those 3.5% of wrong calls are not cosmetic: `family_detect.rs:670` collapses
co-located reps on strand, so one wrong call merges two distinct loci or splits one. Measured cost
(ledger 4p): of 35 dissolved families, 16 = 0.4571 were antisense artefacts (correctly dissolved) but
18 = 0.5143 were still same-strand -- a genuine loss.

THE PROPOSAL. Do not flip on a bare majority. ABSTAIN -- keep the '+' placeholder -- when the vote is not
decisive, so a wrong call degrades to the status quo ante instead of corrupting a locus. This is the
shape O2 is already defended on (`n_decisive >= 1`, assign-or-abstain).

THE MEASUREMENT. SPLICED reps have junction-determined strand = ground truth. For each, tally read
orientation, compute the margin, and sweep a threshold: what accuracy does a given margin buy, and at
what abstention cost? Then read the same margin distribution off the SINGLE-EXON reps, which is where
the rule would actually apply.

⚠ T8 / APPROXIMATE. The shipped vote is over the SKELETON's reads; this queries every primary read at the
rep's locus, which is a superset. It gives the shape of the accuracy-vs-margin trade, not the exact
operating point the binary would see. Confirm the chosen threshold through the real binary before
shipping it.
⚠ COMPARATOR: a constant '+' predictor scores 0.4867 on spliced reps. Accuracy must be read against that,
never against 0.5.
"""
import collections, csv, random, subprocess, sys

BAM = "/mnt/linuxdisk/home/juanfraitu/winloci_data/GGO_ds.bam"
N = "/mnt/linuxdisk/home/juanfraitu/o1_strand/off/dump/ggo.nodes.tsv"
NSPL, NSE = 1200, 800


def vote(chrom, start, end):
    p = subprocess.run(["samtools", "view", "-F", "2308", BAM, f"{chrom}:{start+1}-{end}"],
                       capture_output=True, text=True)
    f = r = 0
    for line in p.stdout.splitlines():
        if int(line.split("\t", 2)[1]) & 0x10:
            r += 1
        else:
            f += 1
    return f, r


rows = list(csv.DictReader(open(N), delimiter="\t"))
rng = random.Random(101)
spl = rng.sample([x for x in rows if int(x["n_exon"]) >= 2], NSPL)
se = rng.sample([x for x in rows if int(x["n_exon"]) == 1], NSE)

lab = []
for x in spl:
    f, r = vote(x["chrom"], int(x["start"]), int(x["end"]))
    if f + r:
        lab.append((f, r, x["strand"]))
sev = []
for x in se:
    f, r = vote(x["chrom"], int(x["start"]), int(x["end"]))
    if f + r:
        sev.append((f, r))

print(f"labelled (spliced, junction-determined) n={len(lab)}   single-exon n={len(sev)}\n")
print(f"{'margin>=':>9} {'called':>8} {'abstain':>8} {'ACCURACY':>9}   {'SE flip':>8} {'SE abstain':>11}")
for thr in (0.0, 0.60, 0.70, 0.80, 0.90, 0.95, 1.00):
    called = [(f, r, t) for f, r, t in lab if max(f, r) / (f + r) >= thr]
    ok = sum(1 for f, r, t in called if ("-" if r > f else "+") == t)
    ab = 1 - len(called) / len(lab)
    sec = [(f, r) for f, r in sev if max(f, r) / (f + r) >= thr]
    flip = sum(1 for f, r in sec if r > f)
    print(f"{thr:>9.2f} {len(called):>8} {ab:>8.4f} {ok/max(1,len(called)):>9.4f}   "
          f"{flip:>8} {1-len(sec)/len(sev):>11.4f}")
print(f"\n  COMPARATOR: a constant '+' predictor scores {sum(1 for _,_,t in lab if t=='+')/len(lab):.4f} on this labelled set.")
print("  'SE flip' = single-exon reps the rule would set to '-'; the rest keep the '+' placeholder.")
