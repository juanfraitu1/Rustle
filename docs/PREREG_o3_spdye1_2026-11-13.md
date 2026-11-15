# PREREG row 1700 (2026-11-13): SPDYE1's hidden copy in HG002, descriptive, one locus

User go 2026-11-13 ("O3: SPDYE1 on HG002"). The 2026-10-23 read-level section (archive `LEDGER_O3.md` "4. Read level:
does any arm find SPDYE1's hidden copy?") already answers the question for minimap2, `candidate_copies`, Strainy and the
frozen Python port: none finds it (all 99 labelled reads have a unique primary on CHM13 copy SPDYE1.18; the hidden copy
is 9% of that pile, under every allele-frequency floor; the port's cluster "looks like a chimera"). The one detector
built since and never run on it is the shipped `copy_phasing` binary.

## Plan (descriptive; one locus is neither PASS nor FAIL)
1. Rebuild the label: HG002 reads whose alignment to the maternal haplotype lies in the donor-only insertion
   (recorded as chr7_MATERNAL 104,183,399-104,189,258, "99 / 342").
2. Run `copy_phasing` at defaults on the family holding that CHM13 copy, HG002 Kinnex on CHM13.
3. Report: the labelled reads' largest cluster (share of them, purity) and its tier (certified / putative / none).
Validity gate for step 1: the recorded interval, taken from the maternal assembly on Fritz, must be absent from CHM13
(a flank-to-flank alignment with a multi-kb gap). If not, the label cannot be rebuilt here and the row stops.

## Outcome (2026-11-13): STOPPED at the validity gate, not run
The maternal assembly on Fritz (`r1644_hg002_unspliced/mat.fa`, GenBank CM039017.1 = chromosome 7, 162,676,710 bp) holds
no insertion at the recorded interval: the interval with 10 kb flanks (25,860 bp) aligns to CHM13 chr7:102,231,314-
102,257,174 end to end, 2 mismatches, no gap. The recorded coordinates belong to another HG002 assembly version (the
workstation's), or to another coordinate frame; the 99 read names and that assembly are not on Fritz.
What stands is the 2026-10-23 section: no arm finds the copy. To run step 2 the label is needed: the 99 read names, or
the workstation's HG002 maternal FASTA. Workspace `r1700_spdye1` holds the two extracted sequences.

## Amendment 1 (2026-11-14, written before any download or alignment)
The stop above was a coordinate-frame error, not a missing insertion. The label was built on hg002v1.1 (archive
`LEDGER_O3.md` section 3: "from hg002v1.1"), where chr7_MATERNAL is 160,955,030 bp; the maternal assembly on Fritz is an
older GenBank version (CM039017.1, 162,676,710 bp), so the recorded coordinates point elsewhere in it. hg002v1.1 is
public (human-pangenomics S3, `T2T/HG002/assemblies/hg002v1.1.fasta.gz`). The user cannot supply the 99 read names and
asked on 2026-11-14 to recalculate. The row stays descriptive (one locus, neither PASS nor FAIL).

Label, archived method unchanged (`bench/famcmp/spdye1_label.py`):
1. hg002v1.1 chr7_MATERNAL 104,183,399-104,189,258 with 20 kb flanks, aligned to CHM13 (asm20, `--eqx`). Gate: one
   alignment spans both flanks and holds an insertion of at least 3 kb. The insertion's donor coordinates are taken from
   that alignment and reported beside the recorded ones (the record says "6.4 kb" over a 5,860 bp span).
2. Family F = the family of the catalog copy nearest the insertion's CHM13 anchor (catalog: A119b on CHM13,
   `~/runs/human_A119b/a1`, as in row 1651). Reads = every HG002 Kinnex read with any record on a copy of F in
   `HG002_kinnex.t2t.bam`.
3. Those reads aligned to the diploid hg002v1.1 (splice:hq, both haplotypes in one index, so the primary is the best
   donor alignment). A read is labelled when its primary has at least 50 aligned exon bases inside the insertion.
Reproduction checks, fixed now (archive: 99 reads, all primaries unique on one CHM13 copy, 9% of a pile of ~1,050, median
NM 7 against CHM13 and 2 against the donor): SAME LABEL = 80 to 120 labelled reads and at least 90% of their CHM13
primaries on one catalog copy. Outside that it is a DIFFERENT LABEL: everything below is still reported, under that
heading, and is not compared with the 2026-10-23 section. Pile share and the two median NM are reported, not barred.

Detector (step 2 of the plan): `candidate_copies --families F --net-reads` (defaults), then `copy_phasing run` at
defaults on that net. Reported: labelled reads in the net; their clusters in `reads.tsv` (largest cluster: share of the
labelled reads, purity = labelled / cluster size); that cluster's fate and whether a certified or putative candidate
holds it. Fewer than 20 labelled reads in the net: the tier result says nothing and is reported as such.

## Outcome of Amendment 1 (2026-11-14, job 4471029, workspace `r1700_spdye1`): run; DIFFERENT LABEL; no tier fires
- **Insertion: reproduced exactly.** hg002v1.1 chr7_MATERNAL 104,183,399-104,189,258 (5,860 bp, the recorded span; the
  record's "6.4 kb" is a rounding of something else) is one insertion in a flank-to-flank alignment to CHM13 chr7 at
  103,953,526. Gate passed. (The 2026-11-13 stop compared the coordinates with an older assembly; the CHM13 position it
  gave, 102.2 Mb, was wrong.)
- **Label: DIFFERENT by the fixed check.** 17,309 reads have a record on family MCL92 (23 copies, A119b catalog);
  435 are labelled (archive 99; band 80-120). 414 of them (0.952) have their CHM13 primary on one copy, MCL92:21
  (chr7:103,868,698-103,880,329), 429 at MAPQ >= 1; 15 sit on MCL687 copies, 6 off the family. Share of that copy's
  pile 0.153 (2,709 primaries; archive 9% of ~1,050); median NM 12 against CHM13, 7 against the donor (archive 7 / 2).
  The shape is the archive's (one copy, a minority of its pile, closer to the donor), the count is 4.4x. Two differences
  are known and neither was tested: the archive labelled only `no_counterpart` reads of its own family table (342), and
  this BAM holds 2.6x the reads on the copy. Nothing below is compared number by number with the 2026-10-23 section.
- **Detector.** `candidate_copies` nets cap at 1,000 reads per family (default `--max-reads`): 1,000 of 20,822, of which
  21 are labelled (floor 20). `copy_phasing` at defaults on that net: 20 piles, 4 haplotypes, **certified 0, putative 0**
  (5,420 s on one family). `candidate_copies` itself puts 20 of the 21 in clusters: 16 in an 81-read cluster and 3 in a
  32-read cluster, both linked to the reference copy as alleles (d 0.003 / 0.008), and 1 in a 10-read flagged candidate
  (`cand_MCL92_4`, purity 0.10).
- **Reading (one locus, descriptive).** The 2026-10-23 answer stands on the rebuilt label: no shipped arm finds the
  hidden copy; its reads are absorbed as alleles of MCL92:21. The cap makes this a 5% sample of the family's reads;
  a run with the net restricted to MCL92:21's pile was not registered and not run.

## Amendment 2 (2026-11-14, written before the run): the detector on the one copy's pile
Amendment 1's detector run saw 21 labelled reads because the family net is capped at 1,000 of 20,822 reads. One more
run, descriptive like the row: net = 1,000 primaries of MCL92:21's pile (2,709), the 1,000 with the smallest sha1 of the
read name (so the labelled share stays the pile's, about 15%, expected about 150 reads), all `used = 1`;
`copy_phasing run --families MCL92` at defaults. Reported: labelled reads in the net; certified / putative counts; for
each candidate, labelled reads / reads (largest labelled cluster: share of the labelled reads in the net, purity).
"Finds the copy" for this report = a certified or putative candidate holding at least half of the net's labelled reads
at purity at least 0.8. Anything else = not found. No further run on this locus.

### Amendment 2 outcome (2026-11-14, job 4471120, 5,470 s)
- Net: 1,000 of the pile's 2,709 primaries; **154 labelled** (0.154, the pile's share).
- `copy_phasing`: piles 3, clusters 28 (4 conflicting), haplotypes 3 (all `unexplained`), **certified 0, putative 1**.
- The putative candidate (`putative_MCL92_0`, on the copy's own locus, consensus 2,754 bp) holds **6 reads**
  (multi 1.000, other holders 0.409). The run does not write a putative candidate's read names, so its labelled count
  is not measured; at most 6 of 154 (0.039) against the bar of 77.
- **Not found** by the written definition, whatever the purity of the 6. With 154 of the copy's reads in front of it
  the detector still does not separate them: the cap was not the reason for Amendment 1's result. No further run.
