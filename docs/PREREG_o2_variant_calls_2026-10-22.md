# Pre-registration: O2 from the sample's own variant calls (2026-10-22)

**User (2026-10-22):** "it just needs to be something that passes from beginning to end, like if I posit that I have a gene
family I need to show it (lets say DAZ or so) then show how if there are PSVs, bam tags or other differences among the possible
placements and show that we got the right placement instead of minimap2 or that we found a possibly missing copy based on the
variants" — one family picked by hand (GTF2IRD2), with backups.

**Design:** `docs/superpowers/specs/2026-10-22-o2-variant-calls-design.md`. The spec's "How reads are scored" and "Not
evaluable" amendments are in commit c63bbd6a. Plan: `docs/superpowers/plans/2026-10-22-o2-variant-calls.md`.

**Frozen code:** `07441e7c` on branch `o3-fast`, holding `bench/o3_graph/o2vc.py`, `o2vc_data.py`, `o2vc_rule.py`. The
weak-het policy is `P1 `keep` (the record's GT alleles whatever FILTER; `--weak-het` default, commit fdb41f25)`, chosen on the development families by the rule in "What was seen before".

**Expected failures, pre-registered separately:** `docs/PREREG_o2_variant_calls_expected_failures_2026-10-22.md` (0b31688c). It
holds 16 failure modes with their signatures, and blind panel predictions.

## The rule (verbatim from the spec's method, summarised)

- **Domain.** A tied read is one whose alignments within 10 AS of its best lie on ≥ 2 copies of one family and on no locus off it.
  minimap2's primary may lie outside the tied set.
- **Teachers.**
  - A copy's teachers are reads primary on it whose every alignment on another copy of the family is ≥ 10 AS lower. A tied read
    never teaches.
  - longcallR 2.0.0 (`-p hifi-masseq --min-mapq 0 --no-bam-output`, defaults otherwise) calls SNVs on each copy's teachers.
- **Genotypes.**
  - **Known:** the copy is known at a position when depth ≥ 6, counted as longcallR's DP counts it. The genotype is the GT
    alleles whatever FILTER, with the weak-het policy `P1 `keep` (the record's GT alleles whatever FILTER; `--weak-het` default, commit fdb41f25)` applied.
  - **Unknown:** `RnaEdit`, `Dense` and depth > 50,000 are unknown.
  - **No record:** {reference base}.
- **Assign or abstain.**
  - **Separates:** a column separates two copies when both are known and their allele sets are disjoint.
  - **Votes:** the read's base, at Q ≥ 10, votes for the copy whose set holds it.
  - **Beats:** c beats d when ≥ 1 separating column votes for c and none votes for d.
  - **Decision:** assign c iff c beats every other tied copy; else abstain with the set of unbeaten copies.
  - **Frame:** all comparisons are made in the read's frame, so a copy aligned opposite to the read has its alleles complemented.
- **Possible missing copy.**
  - **Foreign:** an allele is foreign at a column where every tied copy is known and none holds it.
  - **Flag:** a flag is a connected component of recurrent foreign alleles (≥ 2 reads each) that are linked (≥ 2 reads carry
    both, at different columns).
  - **Editing:** A>G changes in the transcript frame never count toward a flag.
  - **Effect:** a flagged read is not placed.

## Arms (all on the same reads)

| arm | definition |
|---|---|
| N1 | minimap2's primary |
| N2 | N1 at MAPQ ≥ 1 |
| N3 | the same rule, every known genotype = {reference base} |
| N4 | the frozen corrector (`o2correct.py donor` v2), where its output holds the read |
| VC | this design |

## Panels and exact commands (run once each, in this order, after this file is committed)

```bash
cd /mnt/linuxdisk/home/juanfraitu/rustle_soto_subset
T=/mnt/linuxdisk/tmp/o3_graph; G=/mnt/linuxdisk/home/juanfraitu/npip_ladder/idx/target.fa; P=$T/famident/par_families.json
# A + C-human: HG002 real Kinnex reads, every pool family (chrM families out); truth = HG002's own assembly (answer key only)
python3 bench/o3_graph/o2vc.py run $T/hg002proof/chm13.bam $T/o2psv_dev/pool_copies.tsv $G $T/o2vc/hg002 --exclude-chrM
for F in S19GTF2IRD2 S19NBPF1 S19PDPR; do
  python3 bench/o3_graph/o2vc.py score $T/o2vc/hg002 $T/hg002proof/truth.tsv $T/hg002proof/truth_untied.tsv --mode donor \
    --n4 $T/o2psv_dev/hg002.v2.tsv --panel A --focus $F --exclude $P
done
python3 bench/o3_graph/o2vc.py score $T/o2vc/hg002 $T/hg002proof/truth.tsv $T/hg002proof/truth_untied.tsv --mode donor \
  --n4 $T/o2psv_dev/hg002.v2.tsv --panel C --exclude $P
# B: DAZ, A119b testis, no truth (descriptive)
python3 bench/o3_graph/o2vc.py run $T/yag_testis/yag.bam $T/o2vc/inputs/daz_copies.tsv $G $T/o2vc/daz --families S19DAZ1
python3 bench/o3_graph/o2vc.py score $T/o2vc/daz --mode none
# E: HG00733 simulated from its own haplotypes, held out (never run or scored during development)
S=$T/donorsim/HG00733/sim
python3 bench/o3_graph/o2vc.py run $S/chm13.bam $T/o2psv_dev/pool_copies.tsv $G $T/o2vc/hg00733
python3 bench/o3_graph/o2vc.py score $T/o2vc/hg00733 $S/truth.tsv --mode donor --n4 $S/second/o2correct.tsv --panel E --exclude $P
```

**Disclosed precompute.**
- **What was run early:** the longcallR calls (`OUTDIR/vcf/*.vcf`) for `hg002` and `hg00733` were computed BEFORE this file was
  committed. They come from a frozen checkout at a2098a96, in batches, while the policy was being chosen. This was to save time:
  the user has a meeting.
- **Why it is safe:** calls do not depend on the weak-het policy, and `run` reuses saved VCFs. The development-time changes after
  a2098a96 did not touch calling (`git diff a2098a96..07441e7c` touches no calling code. It changes only `parse_vcf(with_af)`, the weak-het interpretation of saved VCFs, the attempt counter and exit code, the scorer and the CLI.).
- **What was not looked at:** the truth-free placements those batches wrote were neither opened nor scored. The commands above
  regenerate every output except the VCFs with the frozen code.

**Deferred (fast path, user 2026-10-22):**
- panel C-gorilla (KB3781 NPIP unit C);
- panel D (one-genome no-harm worlds).

Each gets its own amendment before it runs.

## Scoring (spec, as amended in c63bbd6a)

- **Scored reads.** These are tied reads with a definite truth: classes unique, unique_elsewhere, elsewhere and no_counterpart.
  - A placement is right iff the placed copy is the truth copy. Every other placement on a scored read is wrong. That includes a
    truth copy outside the tied set, `elsewhere`, and `no_counterpart`.
  - Fixes, breaks and wrong/placed are counted over all scored reads, for every arm, on the same reads. A break is a read N1 has
    right that the arm places wrong; this includes outside-primary reads VC moves.
  - The decidable subset (truth among the tied copies) is reported for comparison with earlier registers.
- **Calibration counts.**
  - K1 = VC placements among `own_tie` tied reads (HG002 only).
  - K2 = flagged tied reads whose truth copy is present (unique, unique_elsewhere, own_tie; in E every read).
  - C1 = flags with ≥ 2 `no_counterpart` reads.
- **Exclusions.**
  - PAR families (`famident/par_families.json`) are out of every pooled number and reported apart.
  - chrM families are out of the HG002 run.
- **Not evaluable.** A bar is not evaluable when there is nothing to measure, or when no outcome could meet it for lack of reads:
  - a zero denominator;
  - A2 with < 10 minimap2-wrong scored reads on the focus family;
  - C1 with < 2 `no_counterpart` tied reads;
  - E2 with no VC placement.
- **Verdict precedence.** FAIL if any evaluable bar fails; else UNDERPOWERED if any bar is not evaluable; else PASS.
- **The "fix possible" bound** for A2, A3 and E2 counts every scored read that N1 places wrong, as implemented. That includes
  `elsewhere` / `no_counterpart` reads and reads whose truth copy lies outside the tied set, none of which any arm can fix. It is
  therefore an upper bound on the possible fixes (Task 6 review, minor 1).
  - At scoring, the decidable N1-wrong count (`decidable.arms.N1.wrong`) is reported next to it.
  - If a bar is evaluable only by the upper bound, the result says so.

## Bars (verbatim from the spec)

- **A (focus GTF2IRD2; the backups NBPF1 and PDPR get the same bars and are reported):**
  - A1 VC wrong ≤ 2% of VC placements;
  - A2 fixes ≥ 10;
  - A3 VC fixes > N3 fixes;
  - K1 ≤ 2% and K2 ≤ 2%, pooled over HG002's pool families, PAR out.
- **B:** descriptive.
- **C-human:**
  - C1: at least one flag holds ≥ 2 `no_counterpart` tied reads;
  - C2 (reported): every flag with its reads' truth classes. The audit checks each counted flag's donor locus for a CHM13 copy.
- **E:**
  - E1 wrong ≤ 2% of placements;
  - E2 fixes > breaks;
  - K2 ≤ 2%.

**Decision.**
- The showcase passes iff panel A for GTF2IRD2 passes.
- A failed GTF2IRD2 is a FAIL. The backups are reported, never substituted.
- C-human and E each get their own verdict, never pooled. Human and gorilla are never pooled.

## What was seen before this file

- **The showcase choice.** The corrector's HG002 per-family numbers on tied reads with donor truth chose it: GTF2IRD2 19 fixes /
  0 breaks, NBPF1 11 / 0, PDPR 14 / 0. So did samplecopy's held-out GTF2IRD2 31 / 0 (register 1333). The showcase is hand-picked,
  and hand-picked sets overstate (metric trap 09-18).
- **Development, HG002 SMN2 and MSTRG3150 only** (`docs/O2VC_DEV_2026-10-22.md`): pooled over both development families (1,440 scored tied reads, every arm on the same reads): **VC 1,330 placed, 1 wrong, 647 fixes, 0 breaks**; N1 (minimap2) 656 wrong of 1,440; N2 623 wrong of 1,358; N3 (reference only) 692 placed, 15 wrong, 8 fixes; N4 (corrector) 19 wrong, 638 fixes, 1 break; K1 0/128; K2 3/1,568; depth agreement 41/41. Policy choice: P3 (`ref`) excluded (K1 98/128); P1 (`keep`) and P2 (`unknown`) tie on every term (identical VC on all 1,849 reads), so the spec's P1 stays. ⚠On SMN2, VC's whole gain over N3 rests on ONE LowQual 0/1 record (chr5:70823918, AF 0.14, DP 79, exon 2a) that erases one reference PSV; the audit re-calls SMN2 and GTF2IRD2 five times each.
- **The SMN2 smoke** (Task 4): one LowQual 0/1 at AF 0.14 erased a reference PSV for 1,589 of 1,655 tied reads.
- **P1 vs P2.** The frozen rule does not separate them. VC is identical on all 1,849 dev tied reads; only N3 differs, on 922
  SMN2 reads. P1 is kept as the spec's own rule, the status quo; it is not derived from the rule.
- **Policy choice rule (frozen before scoring the dev runs):**
  1. Exclude any policy whose dev K1 is above 2%.
  2. Take the lowest dev wrong/placed.
  3. Break ties by more fixes.

## Caveats carried into every claim

- **The caller is not held out on HG002.** longcallR's PacBio MAS-seq model was trained on PacBio's HG002 Kinnex cDNA set.
  Provenance of our HG002 reads: the `chm13.bam` header and run IDs (m84036_230422…) look like PacBio's public HG002 Kinnex release, so the reads may be longcallR's own MAS-seq training data. The headers cannot prove it. Any HG002 result is therefore caller-optimistic..
- **SNVs only.** Indel PSVs are not used.
- **Absorbed extra copies are not flagged.** An extra donor copy absorbed by its nearest reference copy is never flagged, and is
  visible only through `no_counterpart` placements.
- **E is simulated.** It has even allele expression, no editing and no donor-extra copies, so it is optimistic.
- **Non-determinism.** longcallR is not deterministic at weak sites (4 of 41 identical runs). The saved VCFs are the record. The
  audit re-calls GTF2IRD2 five times and reports how many placements change.
- **Novelty.** This is DuploMap's reliable-PSV principle (Prodanov & Bansal 2020, DNA) carried to long-read RNA, with an
  assign-or-abstain certificate and an RNA-side missing-copy flag.

## After the runs

1. Append "Results (code `07441e7c`; run once)":
   - one table per panel, with each bar's value and verdict;
   - GTF2IRD2's worked reads from its `summary.md`.
2. Audit with the `metric-auditor` agent on a frozen copy of the outputs.
3. Add register row 1374 and update memory.
