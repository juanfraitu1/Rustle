#!/usr/bin/env python3
"""How often could/does an EM change a paralog-copy assignment vs plain AS-argmax?

For each family multimapper (read mapping to >=2 distinct paralog loci) we look at its best-two
placements and classify:
  - AS-DECISIVE (gap > DEC): the alignment score picks a clear winner. ANY likelihood EM is
    dominated by the same mismatch evidence -> it MUST agree with AS-argmax. EM cannot change these.
  - NEAR-TIE (0 < gap <= DEC): the score margin is within noise. This is the ONLY place an EM's
    abundance prior can legitimately flip the call.
  - FULLY-TIED (gap == 0 AND NM/de/ms equal): unassignable by anything (covers zero PSVs).
We also count AS-argmax vs NM-argmin DISAGREEMENT = reads where the EM's mismatch-likelihood
(edit distance, what the fingerprint-EM scores) points to a DIFFERENT copy than the alignment score.
That disagreement rate is the maximum the EM's *evidence* (not its prior) can move vs argmax.

Usage: python3 em_vs_argmax.py <contig> [<contig> ...]
"""
import sys, subprocess, collections

BAM = "/mnt/c/Users/jfris/Desktop/GGO.bam"
DEC = 5          # AS gap above which the score is "decisive"
DIST = 50000     # min distance between placements to count as distinct copies

def aln_tags(fields):
    AS = NM = ms = None; de = None
    for t in fields[11:]:
        if t.startswith("AS:i:"): AS = int(t[5:])
        elif t.startswith("NM:i:"): NM = int(t[5:])
        elif t.startswith("ms:i:"): ms = int(t[5:])
        elif t.startswith("de:f:"): de = t[5:]
    return AS, NM, de, ms

def run(contigs):
    grand = collections.Counter()
    for contig in contigs:
        reads = collections.defaultdict(list)
        p = subprocess.Popen(["samtools", "view", BAM, contig], stdout=subprocess.PIPE,
                             text=True, bufsize=1 << 20)
        for line in p.stdout:
            f = line.rstrip("\n").split("\t")
            flag = int(f[1])
            if flag & 0x800 or flag & 0x4:   # skip supplementary / unmapped
                continue
            AS, NM, de, ms = aln_tags(f)
            if AS is None:
                continue
            reads[f[0]].append((int(f[3]), AS, NM, de, ms))
        p.wait()

        c = collections.Counter()
        for rn, places in reads.items():
            if len(places) < 2:
                continue
            ps = sorted(places, key=lambda x: -x[1])  # by AS desc
            kept = []
            for pl in ps:
                if all(abs(pl[0] - k[0]) > DIST for k in kept):
                    kept.append(pl)
            if len(kept) < 2:
                continue
            c["family_multimappers"] += 1
            b0, b1 = kept[0], kept[1]
            as_gap = b0[1] - b1[1]
            # classification
            if as_gap == 0 and b0[2] == b1[2] and b0[3] == b1[3] and b0[4] == b1[4]:
                c["fully_tied"] += 1
            elif as_gap > DEC:
                c["as_decisive"] += 1
            else:
                c["near_tie"] += 1
            # EM-evidence (NM) vs AS argmax: does min-NM disagree with max-AS?
            # b0 = max-AS placement. NM-argmin among the best-two:
            if b0[2] is not None and b1[2] is not None:
                if b1[2] < b0[2]:           # the lower-AS copy has FEWER edits -> evidence conflict
                    c["nm_disagrees_with_as"] += 1
                # also flag if NM ties but AS decides (AS using gaps/clip beyond edit distance)
                if b1[2] == b0[2] and as_gap > 0:
                    c["nm_tie_as_breaks"] += 1
        # report
        N = c["family_multimappers"] or 1
        print("=" * 72)
        print(f"CONTIG {contig}: {c['family_multimappers']:,} family multimappers")
        print(f"  AS-DECISIVE (gap>{DEC}, EM MUST agree w/ argmax) : {c['as_decisive']:6d}  ({100*c['as_decisive']/N:5.1f}%)")
        print(f"  NEAR-TIE   (0<gap<={DEC}, EM prior COULD flip)   : {c['near_tie']:6d}  ({100*c['near_tie']/N:5.1f}%)")
        print(f"  FULLY-TIED (unassignable by anything)           : {c['fully_tied']:6d}  ({100*c['fully_tied']/N:5.1f}%)")
        print(f"  --")
        print(f"  NM-argmin DISAGREES with AS-argmax (evidence)   : {c['nm_disagrees_with_as']:6d}  ({100*c['nm_disagrees_with_as']/N:5.2f}%)")
        print(f"  NM tied but AS breaks it (AS uses gaps/clip)     : {c['nm_tie_as_breaks']:6d}  ({100*c['nm_tie_as_breaks']/N:5.2f}%)")
        for k, v in c.items():
            grand[k] += v
    # grand total
    N = grand["family_multimappers"] or 1
    print("=" * 72)
    print(f"ALL CONTIGS: {grand['family_multimappers']:,} family multimappers")
    print(f"  AS-DECISIVE (EM forced to agree with argmax)        : {100*grand['as_decisive']/N:5.1f}%")
    print(f"  NEAR-TIE  (only place EM's prior can legitimately act): {100*grand['near_tie']/N:5.1f}%")
    print(f"  FULLY-TIED (no method can assign)                   : {100*grand['fully_tied']/N:5.1f}%")
    print(f"  EM-evidence (NM) ever contradicts AS-argmax         : {100*grand['nm_disagrees_with_as']/N:5.2f}%")

if __name__ == "__main__":
    run(sys.argv[1:] if len(sys.argv) > 1 else ["NC_073248.2", "NC_086018.1", "NC_073224.2"])
