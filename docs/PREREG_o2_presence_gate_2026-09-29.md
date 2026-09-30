# Pre-registration: a DNA PRESENCE GATE for O2 on a matched individual (HG002)

**Written 2026-09-29 (KEY=o2_presence) BEFORE any per-read truth was computed and before any gated or null arm was
run.** HUMAN only (CHM13 v2.0 reference, HG002 individual); nothing is pooled with gorilla. Nothing in `src/` is
edited; the gate is a prototype outside `src/` (BAM surgery + the shipped binary). Nothing is committed. Scratch:
`/mnt/linuxdisk/tmp/rustle_figures_dev/o2_presence/` (code in `src/`, frozen binaries in `rbin/`).

## 0. The question

O2 assigns an AS-tied read to one of the reference's copies or abstains. The reference's copy set is not the
sample's: HG002 lacks many CHM13-specific k-mers in 22% of Soto's Nearly-Fixed genes (PREREG_soto_parcn_assembly,
C2). **If a read's candidate copies are first restricted to those whose copy-specific k-mers are present in the
individual's own genome (a yes/no restriction, never a 1/k weight), do assignments become more correct, and does
the restriction ever remove the read's true copy?** The individual's own diploid assembly supplies per-read truth.

## 1. Data and what was done before this file (no truth, no arm outcome seen)

- **RNA:** PacBio Kinnex full-length RNA, HG002 (Coriell cells, Revio), `DATA-Revio-HG002-1/2-FLNC/flnc.bam`
  (58.9 GB, 37,248,737 FLNC reads) from downloads.pacbcloud.com at a measured ~12 MB/s (parallel ranges do not help;
  GIAB's Baylor/NCBI copies ran at ~7 MB/s). Only the first 4.0 GB (8 GB fetched, 4 GB used) was read:
  **2,329,643 FLNC reads** (BAM order = ZMW order, i.e. not ordered by gene).
- **Screen (input definition, fixed before any tie was counted in the main set):** a read is kept iff ≥ 30% of its
  CHM13-present canonical 30-mers (1/16 FracMinHash sample) occur ≥ 2 times in CHM13 v2.0 (`src/dupscreen.c`).
  On an unscreened random sample (first 20,000 reads, mapped with the shipped command) the screen keeps 22.6% of
  reads and **59/59 reads that are AS-tied at ≥ 2 distinct loci** (exact tie of the top two AS values, copy_assign's
  definition). 521,027 reads pass; mapping speed (~50 reads/s on 5 cores for SD-rich reads) fixed the **read set =
  the first 120,000 screened reads** (files `b00`+`b01`), decided before any tie was counted in them.
- **Mapping (BASE input):** shipped command `minimap2 -ax splice:hq -uf --eqx -Y -N 50 -p 0.1 --secondary=yes` against
  the prebuilt CHM13 v2.0 splice index (`-K 20M` added only so that 10-min foreground windows keep completed reads;
  batching does not change a read's alignments). 119,997 reads mapped.
- **O1 on this BAM (shipped driver `tools/rustle_pipeline.sh`, binaries frozen in `rbin/`, sha1 copy_assign
  57168f3f…, gw_family_catalog 8f86b2d1…, mcl_families 748ea0bd…, as_table f2abaede…):** `assemble` 2,035 transcripts;
  `families` 418 copies / 120 families; `catalog` (what the driver's `assign` stage consumes) 985 copies / 163 families.
- **BASE O2 already run** (driver `assign` command + `--union-certificate`, 3 min 45 s): AS-tied gate 30,106
  region-local tied molecule entries, 29,564 with a tied placement outside every catalog copy (never assignable);
  695 molecules in the table; contested 488 → **assigned 2 / tied 461 / ambiguous 25**. Its correctness is not known.
- **The tied population (genome-wide, `src/placements.py`):** 1,373 reads whose best AS is reached at ≥ 2 distinct
  loci (overlapping best-AS records on one chromosome = one locus). 652 of O2's 695 table molecules are among them;
  the other 43 are region-local ties only (their genome-wide best lies elsewhere) and are excluded, counted.
  Primary on chr16: 211 reads (DEV); elsewhere: 1,162 (HELD-OUT).
- **Presence instrument already computed (`src/build_q2.py`, `src/presence.py`, the soto_parcn_asm `kc30` counter
  unchanged):** for each of the 18,113 best-AS placements of the 1,373 reads, SPEC = canonical 30-mers fully inside
  the placement's CHM13 span (exons + introns) with CHM13 v2.0 count == 1; `share` = fraction of SPEC present
  (count ≥ 1) anywhere in HG002 v1.1 (both haplotypes, exact count). Seen (instrument only): 17,067 placements have
  < 100 SPEC (undetermined; mostly reads with 25–51 identical short placements on acrocentric arms), 851 present,
  **195 absent**; the share histogram of the resolvable ones is NOT bimodal at 0 (1 placement < 0.05, 19 in
  0.10–0.25, 175 in 0.25–0.50, 128 in 0.50–0.75, 723 ≥ 0.75). 194 reads have ≥ 1 absent candidate (DEV 127,
  HELD-OUT 67; 19 merged loci, chr16 EIF3CL 28.99 Mb and the NPIPB6 region 28.66 Mb carry most dev placements);
  191 of them would keep exactly one candidate. 98.2% of placements carry no SPEC 30-mer inside the read's own
  aligned blocks.

## 2. Instrument (fixed now)

**2.1 Presence rule (the gate, a priori; = the parCN machinery's zero).** Candidate = one best-AS placement of a
tied read (its CHM13 span). **ABSENT iff |SPEC(span)| ≥ 100 and share < 0.5** (median HG002 count over SPEC = 0,
i.e. parCN_HG002 = 0; the ≥ 100 floor is PREREG_soto_parcn_assembly §2.3). |SPEC| < 100 → UNDETERMINED → kept.
No edit-depth filter: QuicK-mer2's < 100-neighbour rule protects short-read depth from mismapping; an exact count
in an assembly has nothing to protect (disclosed deviation from the soto_parcn_asm SPEC; that machinery's counter
and encoding are reused unchanged).

**2.2 Resolvability (the sequence limit, per placement):** `n_spec_blocks` = SPEC 30-mers fully inside the
placement's aligned blocks (what the read covers). A read's true copy is **resolvable** iff n_spec_blocks ≥ 1.

**2.3 Arms** (same read set, same catalog `hg002.cat.copies.tsv/.fa`, same frozen binary, same command:
`copy_assign --bam B --fasta chm13v2.0.fa --regions <all contigs> --families hg002.cat.copies.tsv --copies-fa
hg002.cat.copies.fa --union-certificate`):
- **BASE** — the BAM as mapped.
- **GATE** — for every tied read, every best-AS record at an ABSENT placement is deleted (other records untouched).
  If the deleted set contains the primary, the remaining best-AS record with the lowest record index becomes primary
  (SEQ/QUAL copied, reverse-complemented if the strand differs; MAPQ 60 if it is now the only best-AS locus, else 0).
  If a catalog copy loses every read, it is removed from the catalog handed to that arm (it is ABSENT by the rule;
  copy_assign refuses readless copies) and the removal is counted.
- **NULL** — per read, the same NUMBER of best-AS candidate loci as GATE deleted, chosen uniformly at random among
  that read's candidate loci; everything else as GATE. **5 seeds** (1–5); mean and range reported.

**2.4 Decision per read and arm** (the pipeline's own semantics): candidate loci left = 0 → **abstain**
("all candidates absent", counted); exactly 1 → **assigned by elimination** to it (the read is no longer AS-tied;
O2 skips it and its unique best placement is the answer, as for every non-tied read); ≥ 2 → O2's verdict: assigned
iff a row is `assigned` with `origin_rejected = 0` (union verdict when the read is in the union certificate's
scope), else **abstain** (tied / ambiguous / no row). Two assigned rows naming different loci = **wrong** unless
both are the true locus. BASE never eliminates (every read in the population has ≥ 2 candidate loci).

**2.5 Per-read TRUTH from the individual's own genome.**
1. HG002 target set (a small index, never whole-genome): (a) every tied read's CHM13 best-AS placement spans lifted
   to HG002 through the Q100 project's one-to-one chains `CHM13v2.0_to_hg002v1.1.{mat,pat}.chain.gz` (nf-LO +
   rustybam trim, homologous chromosomes only); ∪ (b) every 1-kb HG002 bin holding ≥ 3 exact 30-mers shared with
   the tied reads whose HG002 count is ≤ 20 (a homology screen that also finds HG002 copies absent from CHM13);
   bins merged across gaps < 150 kb, padded ± 10 kb.
2. The 1,373 reads are mapped to the target set with the shipped command (splice:hq, -N 50, -p 0.1).
3. **Own-genome copy** = the read's best-AS HG002 alignments. Their aligned blocks are lifted to CHM13 through
   `hg002v1.1_to_CHM13v2.0.chain.gz` (t = HG002, q = CHM13). **Correspondence:** a best HG002 alignment corresponds
   to the read's CHM13 candidate whose aligned blocks overlap ≥ 50% of the lifted aligned bases.
4. Truth classes: **RESOLVED** (every best HG002 alignment lifts ≥ 50% of its bases and all correspond to ONE CHM13
   candidate locus — MAT and PAT copies of one locus agree by construction); **OUTSIDE** (lifts, but onto no
   candidate of the read: the syntenic copy is not among the reference's best loci — kept; any assignment is wrong,
   abstention is right); **UNLIFTED** (a best alignment lifts < 50%: HG002 sequence with no CHM13 correspondent,
   e.g. an HG002-specific copy, or trimmed by the one-to-one chains — excluded, counted); **OWN-TIE** (best HG002
   alignments correspond to ≥ 2 different CHM13 candidates: unresolvable even in the own genome — excluded,
   counted); **WORSE** (best AS in HG002 < best AS in CHM13: the target set probably misses the origin — excluded,
   counted). Correctness: an assigned catalog copy is correct iff ≥ 50% of the true placement's aligned bases lie
   inside the copy's span; an elimination is correct iff the surviving locus is the true locus.
5. **Presence-call check (independent of reads):** an ABSENT-called placement is **syntenically present** iff ≥ 50%
   of its span's bases lift to HG002 (either haplotype) through the CHM13→HG002 chains.

**2.6 NM-identical twins:** for a read, candidates whose NM equals the minimum NM among its remaining candidates.
A read with ≥ 2 such candidates after gating is **NM-identical** and must stay an abstention.

## 3. Metrics and clauses (evaluated in this order; verdict strings fixed now)

Per arm × stratum (DEV chr16 / HELD-OUT): assigned-correct, assigned-wrong, abstained, split by resolvability of the
true copy (2.2); the abstentions split by NM-identical (2.6); the truth-class counts; per absent locus (19) the
eliminations and their truth. Denominator = every RESOLVED + OUTSIDE read of the 1,373 (never conditioned on an arm).

- **C0 (gate on the truth; nothing below is claimed if it fails).** (a) RESOLVED + OUTSIDE ≥ 50% of the 1,373;
  (b) among them the own-genome best AS ≥ the CHM13 best AS for ≥ 95%. Fail → "TRUTH NOT ESTABLISHED".
- **C1 (primary, HELD-OUT).** Precision of GATE's eliminations = correct / (correct + wrong) among reads GATE
  assigns and BASE does not. ≥ 0.95 → "GATE ADDS CORRECT ASSIGNMENTS"; [0.80, 0.95) → "MOSTLY CORRECT"; < 0.80 →
  "GATE MISASSIGNS".
- **C2 (primary, against NULL, HELD-OUT).** GATE wrong < min over the 5 NULL seeds' wrong AND GATE correct > max
  over NULL correct → "GATE BEATS NULL"; else "NOT DISTINGUISHABLE FROM RANDOM REMOVAL".
- **C3 (no harm).** Reads BASE assigns correctly that GATE turns wrong or abstain: 0 → "NO HARM"; else "HARM (n)".
- **C4 (NM-identical invariant, every arm).** Assignments among reads that are NM-identical after gating: 0 →
  "NM-IDENTICAL ABSTENTIONS HOLD"; else "VIOLATED (n)".
- **C5 (presence ≠ absence?).** Share of ABSENT-called placements (by merged locus) that are syntenically present
  (2.5.5): ≥ 50% → "THE RULE FLAGS ALLELES, NOT COPIES"; ≤ 20% → "THE RULE FLAGS MISSING COPIES"; between →
  "MIXED". Reported per locus.
- **C6 (DEV chr16).** C1–C4 statistics reported on the development chromosome, no separate bar (it is where O2 was
  developed; HELD-OUT decides).
- **D (descriptive).** Resolvability split (how many true copies have ≥ 1 SPEC 30-mer in the read's blocks),
  OUTSIDE/UNLIFTED/OWN-TIE counts, NM-identical abstentions per arm, BASE O2's own accuracy on its 2–4 assignments.

## 4. Predictions (this author, before the truth)

1. C0 passes (Q100 chains are one-to-one and cover > 99% of non-SD sequence); OWN-TIE is large (≥ 20%), because
   exact CHM13 ties are mostly exact copies that are also identical in HG002 (X/Y PAR, acrocentric arrays).
2. **C5 "MIXED" or "FLAGS ALLELES"**: the share histogram piles up at 0.25–0.50, the signature of PSV alleles that
   differ between CHM13 and HG002, not of a deleted copy (which would sit near 0).
3. Following from 2, **C1 < 0.95** (I expect 0.5–0.8: a present copy with a non-CHM13 PSV haplotype is gated, and
   its reads are then eliminated onto the sibling) — "MOSTLY CORRECT" at best, possibly "MISASSIGNS".
4. C2: GATE beats NULL on correct (the rule is not random) but its wrong count is not below every seed's.
5. C3 "NO HARM" (BASE assigns ~2 reads; none is at an absent locus). C4 holds in BASE and GATE by construction of
   the certificate; eliminations never touch an NM-identical set with ≥ 2 kept twins.
6. Resolvability: ≥ 95% of RESOLVED reads' true copies have no SPEC 30-mer inside the read's blocks — the sequence
   limit is the binding one, which is why BASE assigns almost nothing.

## 5. Falsifiers

- C5 "FLAGS ALLELES" falsifies the premise that k-mer absence = copy absence at the parCN = 0 cut.
- C1 < 0.80 on HELD-OUT means the gate as stated is unsafe for O2 (assign-or-abstain must not buy yield with error).
- C2 failing means the gate's gains are what any candidate removal of that size buys.
- Any C4 violation means the gate (or the certificate) assigns reads the sequence cannot resolve.

## 6. Hostile self-review (before running)

- **Tiny O2 scope.** BASE assigns 2 reads: "wrong assignments fall" cannot be tested on O2's own certificate here;
  what is testable is whether the restriction's eliminations are right. Said up front, not discovered.
- **Few loci.** 194 affected reads sit on 19 merged loci; EIF3CL alone is most of DEV. Every number is also
  reported per locus; a HELD-OUT verdict resting on one locus is flagged as such.
- **The screen and the 120k cut** shape the population; validated for ties (59/59), disclosed.
- **Truth is aligner-based.** The own-genome alignment uses the same minimap2 settings; its ties are real
  indistinguishability in HG002 and are excluded (OWN-TIE), never resolved by a rule.
- **One-to-one chains** drop overlapping alignments in SDs (rustybam trim): some true correspondences are missing
  → UNLIFTED/OUTSIDE inflated; reported, not repaired.
- **Cell line.** Kinnex HG002 RNA is from Coriell cells, the assembly from HG002 cell-line DNA; somatic CNV/EBV
  changes are possible and would appear as OUTSIDE/UNLIFTED.
- **Seen before writing:** the share histogram, the reach counts (194/191), BASE's 2/461/25. Not seen: any truth,
  any gated or null arm, any presence-vs-chain comparison.
- **Circularity check:** the presence rule uses HG002 k-mers; the truth uses HG002 alignments + chains. Both read
  the same assembly, so an assembly error at a locus would bias both the same way; the chains are the only
  non-k-mer ingredient and come from the Q100 project, not from this analysis.
- No parameter (0.5, 100, 50% overlap, 3 k-mers/bin, 150 kb merge) is changed after the truth exists.

---

# OUTCOME (2026-09-29, appended after scoring; sections 0–6 unchanged since sha1 40820f5e)

Scripts (all in `/mnt/linuxdisk/tmp/rustle_figures_dev/o2_presence/src/`): `truth.py` (truth classes + C5),
`gate_bam.py gate|null:<s>` (arms), `run_o2.sh` (the one O2 command), `score.py` (decisions §2.4, tables §3).
Per-read table: `score/per_read.tsv`. Arms: BASE, GATE, NULL seeds 1–5, same frozen `copy_assign` (sha1 57168f3f).
O2's contested counts: BASE 2 assigned / 461 tied / 25 ambiguous; GATE and every NULL 2 / 328 / 21 (the reads the
restriction touches leave the contested set in both).

## C0 — ⛔ "TRUTH NOT ESTABLISHED" (the pre-registered gate fails; everything below is DESCRIPTIVE, not a verdict)

Truth classes over the 1,373 tied reads: **RESOLVED 337, OUTSIDE 47** (together 384 = **0.280 < 0.50**, C0a fails),
UNLIFTED 497, OWN-TIE 481, WORSE 11 (C0b: own-genome best AS ≥ CHM13 best AS for 1,362/1,373 = 0.992, passes).
Why C0a fails: **462 of the 497 UNLIFTED reads sit on the acrocentric short arms** (chr13/14/15/21/22 p-arms: the
Q100 chains keep homologous-chromosome alignments only, and those arms are not lifted by construction; spot checks
of five such spans: no chain block); OWN-TIE is chr16 192 (every EIF3C/EIF3CL and NPIPB6-region read also ties in
HG002) and X/Y PAR 192. By stratum: DEV 12 RESOLVED + 6 OUTSIDE of 211; HELD-OUT 325 + 41 of 1,162. The bar was
ill-posed for this population (I did not foresee that 36% of exact ties are acrocentric-arm reads); it is reported,
not re-cut. **Where the gate acts, truth coverage is high:** of the 67 HELD-OUT reads the gate touches, 64 are
RESOLVED/OUTSIDE (61 R, 3 O; 2 WORSE, 1 OWN-TIE). DEV: all 127 touched reads are OWN-TIE (unscorable).

## C5 — "THE RULE FLAGS ALLELES, NOT COPIES" (independent of the per-read truth)

**17 of the 19 ABSENT-called loci have a syntenic HG002 counterpart** (156 of 158 placement spans lift ≥ 50% to a
haplotype through `CHM13v2.0_to_hg002v1.1.{mat,pat}`): EIF3CL chr16:28.99 Mb lifts 1.000 to PATERNAL (0.02 MAT) —
HG002 carries one EIF3CL copy; the NPIPB6 region 28.66 Mb lifts to both (0.95–0.98 / 0.999); chr12:93.62 Mb both
(1.000/1.000); chr17 47.1–47.5 Mb (LRRC37A/KANSL1, 17q21) PAT and/or MAT; X-PAR1 to chrX_MATERNAL (0.96–0.999).
The two without a counterpart: **chr7:76.59 Mb (share 0.000 — the one genuinely absent copy)** and chr22:5.55 Mb
(an acrocentric p-arm the chains cannot lift, so uninformative). 17/19 = 0.89 ≥ 0.50.

## C1 — descriptive string "GATE MISASSIGNS" (HELD-OUT)

GATE assigns 61 reads BASE does not, all **by elimination**: **13 correct, 48 wrong → precision 0.213**. In all 48
wrong cases the gate deleted the read's TRUE copy (`true_removed` = 48). DEV: nothing scorable (above).

## C2 — descriptive string "NOT DISTINGUISHABLE FROM RANDOM REMOVAL" — in fact WORSE than every seed

| HELD-OUT, 64 touched reads | correct | wrong | abstain |
|---|---|---|---|
| BASE | 0 | 0 | 64 |
| **GATE** | **13** | **48** | 3 |
| NULL seeds 1–5 | 21 / 28 / 30 / 27 / 23 | 40 / 33 / 31 / 34 / 38 | 3 |

GATE has fewer correct and more wrong than **every** NULL seed. Per locus (correct/wrong, GATE vs NULL range):
**chr12:93.62 Mb 0/24 vs 5–12/12–19** (the gate removes the true copy every time); chr17 LRRC37A/KANSL1 7/12 vs
7–10/9–12; X-PAR1 6/9 vs 5–9/6–10; three single-read loci (chr15:57.57, chr9:40.11, chr9:42.06 Mb) 0/3. **Robustness:** without chr12, GATE 13/24 vs NULL
14–19 / 18–23; restricted to reads whose own-genome AS margin (best vs next distinct AS in HG002) is ≥ 10, GATE 11/21
vs NULL 11–16 / 16–21. The gate is at best random removal, and systematically wrong at one locus.

## C3 — "NO HARM" (n = 1)

BASE assigns 3 reads in the denominator: 1 correct (HELD-OUT, chr5) and **2 wrong (DEV chr16:21.34/21.68 Mb, both
OUTSIDE: the reads' syntenic copy is not among their CHM13 candidates; NM 157 and 283)**. The one correct read stays
correct under GATE.

## C4 — "VIOLATED (4)" in EVERY arm, the same 4 reads — not caused by the gate

O2 as shipped (`--union-certificate`, `origin_drop_indels` on) assigns 4 reads whose candidates tie on BOTH AS and NM
(n_decisive 3–12): substitution-only PSV evidence differs when the NM tie is split between substitutions and indels.
Outcomes: 1 correct, 2 wrong (OUTSIDE), 1 OWN-TIE. The gate changes none of them. Reading the intended concern (NM-
identical twins must stay abstentions) on the BASE state: of 242 denominator reads that are NM-identical in BASE,
GATE resolves 60 by elimination — **13 correct, 47 wrong**; NULL seeds 21–30 correct, 30–39 wrong.

## D — descriptive

- Resolvability: of 337 RESOLVED reads, 135 (40%) have ≥ 1 CHM13-unique 30-mer inside the true placement's aligned
  blocks, 202 (60%) none. BASE's one correct assignment is resolvable (its 2 wrong ones are OUTSIDE, no true locus among the
  candidates); eliminations fall in both (12/43 and 2/5 correct/wrong among non-resolvable / resolvable).
- Own-genome margins of the eliminated reads: median 10 AS (≈ 2 mismatches), min 5 — HG002 separates the copies the
  reference ties.
- The only truly absent copy (chr7:76.59 Mb, share 0.000, no chain counterpart) is the case the gate was designed
  for; it touches 1 read (an abstention, truth not scorable).

## The mechanism (why the gate is anti-correlated with the truth)

A read from HG002's copy A ties EXACTLY between CHM13's A and B when HG002's A carries non-CHM13 alleles at A-vs-B
paralogous sites (making the read equidistant from the two reference copies). Those same alleles erase A's
CHM13-specific 30-mers from HG002, so A's share falls below 0.5 and the gate calls A absent. **The condition that
sends a read to O2 (the tie) selects exactly the copies whose reference version the individual does not carry, and
k-mer absence cannot tell "this copy is missing" from "this copy carries different alleles".** At chr12:93.62 Mb (a
1.4-kb span with 119 CHM13-unique 30-mers) ~3 HG002 variants take the share to 0.25, and all 24 reads come from that
very copy (spot check: AS 1444/1444 in CHM13, 1454 at the chr12:93.6 counterpart on BOTH haplotypes vs 1449/1444 at
the other locus in HG002). Short spans make it worse (one variant kills up to 30 k-mers).

## Predictions vs outcome

1. C0 passes ✗ (0.280; acrocentric arms). OWN-TIE ≥ 20% ✓ (35%).
2. C5 "MIXED" or "FLAGS ALLELES" ✓ ("FLAGS ALLELES", 17/19).
3. C1 < 0.95, expected 0.5–0.8 ✓ direction, ✗ magnitude (0.213).
4. GATE beats NULL on correct but not on wrong ✗ — worse than every seed on both.
5. C3 "NO HARM" ✓ (n = 1); C4 holds in BASE and GATE ✗ (4 certificate assignments on NM-tied reads in every arm).
6. ≥ 95% non-resolvable ✗ (60%).

Falsifiers triggered: C5 (k-mer absence ≠ copy absence at the parCN = 0 cut), C1 < 0.80 (unsafe), C2 (no better than
random removal). **Decision: the DNA presence gate as stated must NOT be adopted for O2.** What the data support
instead: presence must be judged by SYNTENIC correspondence (the chain/assembly test used here as truth), not by the
survival of reference-specific k-mers; and the copy set O2 needs is the individual's own (HG002 separates 61/64 of
these reads by ≥ 1 mismatch), i.e. assign against the sample's copies, not gate the reference's.

## Hostile self-review of the outcome

- C0 failed; every number above is descriptive. The core result nonetheless rests on reads with a resolved truth
  (64/67 touched held-out reads) and on C5, which does not use the per-read truth.
- One locus (chr12, 24 reads) carries the "worse than random" contrast; without it the gate is merely ≈ random
  (13/24 vs 14–19/18–23). Loci, not reads, are the independent units: 6 held-out loci with eliminations, 0 where
  the gate beats every NULL seed.
- Truth uses HG002 alignments + the Q100 chains; the presence rule uses HG002 k-mers. An assembly error would bias
  both alike; it cannot produce the observed anti-correlation (truth says the read comes from A; k-mers say A's
  reference form is missing — consistent, not circular).
- Small O2 scope: BASE assigns 3 denominator reads (1 correct, 2 wrong); nothing here measures the PSV certificate's
  accuracy at scale. The 2 wrong BASE calls are OUTSIDE reads (reference bias), a finding for O3 rather than O2.
- Population: first 120,000 screened reads of 2.33M (0.26% of raw reads are exact distinct-locus ties); a 4× larger
  sample would add reads, mostly at the same loci.

## Draft register rows (suffix E, NOT appended to `docs/NEGATIVE_RESULTS_REGISTER.md`)

See the session report `figs/o2_presence.md` (rows 1165E–1168E).
