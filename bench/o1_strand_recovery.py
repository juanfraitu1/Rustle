#!/usr/bin/env python3
"""Can a single-exon rep's strand be MEASURED instead of guessed?

`denovo_assemble.rs:1010` stamps `strand.unwrap_or('+')` when no canonical junction determines
orientation. Result: ALL 5,928 single-exon reps carry '+' (zero '-'), while spliced reps split
0.4867/0.5133. Because a rep's `seq` is stored in its `strand` orientation, a wrongly-'+' rep is stored
reverse-complemented and aligns MINUS to its true paralogue -- which the orientation guard discards.
Measured: 3,951/4,009 = 0.9855 of guard-blocked pairs involve such a rep (vs 0.3943 among kept pairs,
2.50x enrichment).

THIS TEST. Use SPLICED reps as a LABELLED SET: their strand IS junction-determined, so it is ground
truth. Ask whether the reads' own orientation recovers it. If it does, the same signal determines strand
for single-exon reps and the guard keeps its precision while losing its artefact cost.

⚠ The BAM is minimap2 `-ax splice:hq -uf`, so spliced alignments carry `ts:A:+/-` (transcript strand) and
FLAG 0x10 gives read orientation. `-uf` means ts is relative to the FORWARD transcript, so the effective
transcript strand is ts XOR reverse-flag. Both signals are reported separately below; do not merge them
without looking.

⚠ COMPARATOR: a constant predictor ("always +") scores 0.4867 on spliced reps. Any recovery rule must
beat that, not merely exceed 0.5.
"""
import collections, csv, random, subprocess

BAM = "/mnt/linuxdisk/home/juanfraitu/winloci_data/GGO_ds.bam"
N = "/mnt/linuxdisk/home/juanfraitu/o1_reps2/dump/ggo.nodes.tsv"

rows = [r for r in csv.DictReader(open(N), delimiter="\t")]
spliced = [r for r in rows if int(r["n_exon"]) >= 2]
rng = random.Random(101)
samp = rng.sample(spliced, 400)          # random, not the biggest or the most interesting
print(f"labelled set: {len(samp)} spliced reps sampled at random from {len(spliced)}")

agree_ts = tot_ts = agree_fl = tot_fl = 0
for r in samp:
    reg = f'{r["chrom"]}:{int(r["start"])+1}-{r["end"]}'
    p = subprocess.run(["samtools", "view", "-F", "2308", BAM, reg], capture_output=True, text=True)
    ts_votes, fl_votes = collections.Counter(), collections.Counter()
    for line in p.stdout.splitlines():
        f = line.split("\t")
        rev = (int(f[1]) & 0x10) != 0
        fl_votes["-" if rev else "+"] += 1
        ts = next((x[5:] for x in f[11:] if x.startswith("ts:A:")), None)
        if ts:
            # -uf: ts is relative to the forward transcript; flip it when the read is reverse-mapped
            eff = ts if not rev else ("+" if ts == "-" else "-")
            ts_votes[eff] += 1
    truth = r["strand"]
    if ts_votes:
        tot_ts += 1; agree_ts += (ts_votes.most_common(1)[0][0] == truth)
    if fl_votes:
        tot_fl += 1; agree_fl += (fl_votes.most_common(1)[0][0] == truth)

print(f"\n  ts:A: majority vs junction-determined strand : {agree_ts}/{tot_ts} = {agree_ts/max(1,tot_ts):.4f}")
print(f"  FLAG 0x10 majority vs junction-determined     : {agree_fl}/{tot_fl} = {agree_fl/max(1,tot_fl):.4f}")
print(f"  COMPARATOR — constant '+' predictor           : 0.4867  (the spliced base rate)")
