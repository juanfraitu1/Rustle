# Pre-registration — gorilla TBC1D3 as a HELD-OUT family for the NPIP loss hypothesis

**Written 2026-09-29 01:38 PDT, before any read at a TBC1D3 copy is opened** (no BAM record, no assembly, no GTF line
at the TBC1D3 contigs has been looked at in this study). GORILLA ONLY: OR6737 testis (`GGO_mm.bam`) and KB3781
fibroblast (`GCA_029281585.2_flnc_mm.bam`) are analysed and reported separately and never pooled; no human file is
used. Dev diagnostic: nothing in `src/` is edited, no default changes, nothing is committed. Scratch
`/mnt/linuxdisk/tmp/rustle_figures_dev/ggo_tbc1d3_holdout/`.

## 1. The hypothesis under test and why TBC1D3

**NPIP hypothesis** (`ggo_npip_loss` 2026-09-29; memory `project_ggo_npip_loss`): gorilla NPIP copies that have reads
get no locus *because the individual's reads are divergent and carry non-canonical junctions*, not because of an
assembler rule. Its evidence was developed on NPIP only:

- all 31 absences (OR 14 / KB 17 of 25) are decided at three steps: pass-1 floor 8, gate 16 (strict junctions and/or
  the single-exon '+' placeholder), polish mono floor 7; 0 at seeding, dedupe, other gate clauses or merging;
- 70% (OR) / 73% (KB) of the copies' own spliced primaries carry ≥ 1 non-canonical junction (control 2.9% / 3.4%);
- median `de` of canonical-only own reads 0.0059 / 0.0032 vs control 0.0015 / 0.0013 (3.9× / 2.5×);
- completeness is lost before the assembler (49/50 copy-samples have no whole-chain read; 0 complete shipped);
- every single switch that rescues a copy maps to a closed register row.

The 2×2 variant/truncation simulation (`PREREG_ggo_npip_variant_sim_2026-09-29.md`, sibling study, outcome NOT known
when this was written) asks whether the individuals' copy sequences and/or 5′-truncated read lengths reproduce that loss.

**Why a held-out family** (memory `feedback_hold_a_substrate_back`): a mechanism established on the family it was
developed on is not established. TBC1D3 is the other core-duplicon family of the thesis (NPIP/TBC1D3 layer-order,
certificates, nested lattice), it is gorilla-annotated, and this study has not opened its reads.

## 2. Truth — from annotation only (built before any read; `code/build_truth.py` 34d27c47)

- Source: `winloci_data/GGO_genomic.gff` (gorilla RefSeq; the two contigs subset to `ann/ggo2.gff`).
- Records: every `gene`/`pseudogene` whose `description` matches `^TBC1 domain family member 3[A-Z]` (3G, 3G-like,
  3K-like, 3B-like). TBC1D30/31/32 (`member 30/31/32`) are excluded by the pattern.
- ⚠ **14 records, not 15.** The brief said 13 on NC_073228.2 + 2 on NC_073224.2. A grep of the WHOLE GFF for the
  description or a TBC1D3 `product=` finds 12 gene records on NC_073228.2 (31.89-62.53 Mb) and 2 on NC_073224.2
  (224.19-224.20 Mb), and nothing elsewhere besides TBC1D30/31/32. All 14 are used; none is dropped.

| cid | record | contig:start | str | type | tx | distinct canonical chains | primary junctions / length | annotated non-canonical intron |
|---|---|---|---|---|---|---|---|---|
| t00 | LOC101151653 (3G-like) | 228:31,893,545 | − | protein coding | 1 | 1 | 13 / 2,119 | — |
| t01 | LOC129533458 (3G-like) | 228:32,009,446 | + | protein coding | 1 | 1 | 13 / 2,142 | — |
| t02 | LOC101144080 (3G) | 228:44,678,474 | − | protein coding | 1 | 1 | 13 / 2,970 | 43 bp CA..CT (in the only model) |
| t03 | LOC115933306 (3G-like) | 228:54,239,509 | + | protein coding | 4 | 3 | 14 / 2,141 | 43 bp TG..AG (1 of 4 models) |
| t04 | LOC129533792 (3G-like) | 228:54,292,843 | − | protein coding | 3 | 2 | 13 / 2,144 | 43 bp CT..CA |
| t05 | LOC101125558 (3G) | 228:54,400,033 | − | protein coding | 2 | 1 | 13 / 2,126 | 43 bp CT..CA |
| t06 | LOC129533808 (3G-like) | 228:54,880,881 | + | protein coding | 2 | 1 | 13 / 2,125 | 43 bp TG..AG |
| t07 | LOC129533806 (3G-like) | 228:54,930,303 | − | protein coding | 2 | 1 | 13 / 2,119 | 43 bp CT..CA |
| t08 | LOC115934662 (3G-like) | 228:56,462,389 | + | protein coding | 2 | 1 | 13 / 2,131 | 43 bp TG..AG |
| t09 | LOC129533797 (3G-like) | 228:56,531,258 | − | protein coding | 3 | 2 | 15 / 2,497 | 43 bp CT..CA |
| t10 | LOC129533813 (3G-like) | 228:56,834,995 | + | protein coding | 2 | 1 | 13 / 2,119 | 43 bp TG..AG |
| t11 | LOC109026840 (3G) | 228:62,522,473 | − | protein coding | 1 | 1 | 10 / 1,277 | — |
| t12 | LOC115931404 (3B-like) | 224:224,191,304 | − | **pseudogene** (exons parented to the gene) | 1 | 1 | 4 / 1,036 | 31 bp CG..AA |
| t13 | LOC134759231 (3K-like) | 224:224,196,681 | − | protein coding | 1 | 1 | 3 / 423 | — |

(228 = NC_073228.2, 224 = NC_073224.2; motifs are +-strand genomic.)

- **Construction mirrors `ggo_npip_sim/build_tx.py` (FL arm):** territory = merged union of the record's raw annotated
  exons (0-based half-open; the record is its own "native", nothing clipped); every non-canonical annotated intron is
  merged (10 introns: nine 43-bp last introns of alternative models, one 31-bp pseudogene intron — RefSeq
  indel-correction introns, `ann/nc_introns.tsv`); primary chain = most junctions, then longest, then id; junction =
  (last exon base, first exon base), 1-based. 26 transcripts → **18 distinct FL chains** (`ann/tx.tsv`).
- Territories do not overlap. t12/t13 are adjacent (313 bp apart, same strand; plausibly one copy split into a 3′
  pseudogene part and a 5′ 3K-like part): a read on both counts as an own read of both (the diagnostic's definition),
  and is reported.
- **Prior read-level knowledge, disclosed:** (a) memory `project_crossspecies_expansions` (07-16): gorilla testis
  `copy_assign` gave "TBC1D3 ~100 reads/locus × 3" — so OR6737 very likely has reads and loci at ≥ 3 copies;
  (b) r810: gorilla TBC1D3 is a 0.997-identity clique with a 0.883 halo; (c) r822/r932: RefSeq names 0 gorilla
  TBC1D3 by symbol (LOC ids). Nothing about non-canonical junctions, divergence or loss steps at TBC1D3 is known.

## 3. Part A — the loss diagnostic on TBC1D3, and the transfer clauses

**Instruments (the NPIP diagnostic's code, path-only copies; §6).** Binaries `cc_bin_frozen/` (= the NPIP studies':
`copy_assign` 57168f3f, `as_table` f2abaede, `mcl_families` 748ea0bd). Regions = the two contigs WHOLE
(`tbc_contigs.txt`: NC_073228.2:0-195332687, NC_073224.2:0-243847345), as the NPIP run used whole contigs. AS table =
each sample's genome-wide `molecules.tsv` (byte-identical to the 09-25 runs' tables, checked). Arms = the diagnostic's
set via `arm_diag.sh`: ship, nopol, maj, maj_nopol, keepdup(_nopol), noseed(_nopol), readstrand(_nopol),
allsec(_nopol), majrs(_nopol), nofloor, nofloor_rs (light lock); floor1(_nopol) (heavy lock).

**Instrument checks (must pass before Part A is interpreted):**
- I1: the region-restricted `ship` GTF equals the 09-25 genome-wide (unrestricted) GTF on both contigs, line by line,
  cov/TPM/FPKM removed (`bytecmp.py`), both samples. If it differs, the difference is stated and the regional run is
  what is scored.
- I2: `validate.py` — the emulated gate survivors equal the binary's no-polish transcripts on all 14 territories in
  strict / majority / read-strand modes, both samples.
- I3: `polish_port.py` — ALL IDENTICAL on every contig for nopol→ship, maj_nopol→maj, readstrand_nopol→readstrand.

**Transfer clauses** (per sample; `transfer.py` 7b60b4d3; numbers from the diagnostic's own definitions):

- **T1 — same three losing steps.** Absent copies WITH reads (first losing step ≠ "0 no reads") are all lost at the
  pass-1 floor, the gate (junction and/or placeholder label) or the polish mono floor, with ≤ 1 exception per
  sample. Power: ≥ 3 absent-with-reads copies in the sample, else NO POWER.
- **T2 — non-canonical junctions far above background.** Own spliced primaries at the copies with ≥ 1
  non-canonical junction: fraction ≥ 0.30 AND ≥ 10 × the control fraction (the NPIP control: 40,000 spliced
  primaries, NC_073242.2 / NC_073244.2 40-60 Mb, no TBC1D3 copy there). Power: ≥ 20 own spliced primaries.
- **T2b — not only the annotated indel introns.** The T2 rule with junctions that exactly equal an annotated
  non-canonical intron (§2) excluded (`recurrent_tbc.py` `nc_unannotated`).
- **T3 — divergence elevated.** Median `de` of own canonical-only spliced primaries ≥ 2 × the control's
  canonical-only median (NPIP 3.9× / 2.5×). Power: ≥ 20 reads.
- **T4 — completeness lost before the assembler.** Of copies with own reads, ≥ 80% have no whole-chain read, and
  ≤ 1 copy is complete in the shipped output.
- Descriptive, not decisive: T6 — the non-canonical fraction of absent vs present copies; the single switches that
  rescue each absence and their register rows; whether present holders are 3′ fragments; merges of two copies into
  one locus (tandem pairs t03/t04 and t06/t07 are 38 kb apart).

**Diagnostic verdict per sample:** TRANSFERS = T2, T2b, T3, T4 pass and T1 passes or has no power; DOES NOT
TRANSFER = T2 fails (with power); PARTIAL = anything else (e.g. T2 passes only through the annotated 43-bp introns,
or T1 fails with power); NO POWER = T2 has no power. Overall TRANSFERS only if both samples transfer.

**Predictions (committed now): TRANSFERS in both samples** — T1 pass or no power; T2 ≥ 0.30 and ≥ 10× control; T2b
pass; T3 ≥ 2×; T4 pass.

**What falsifies transfer:**
- T2 fails in either sample: TBC1D3's reads are not atypically non-canonical, so the NPIP cause is NPIP-specific.
- T1 fails with power: absences with reads arise at another step (e.g. seeding when every read ties to another copy
  — a multi-mapping/O2 cause, not divergence).
- T2 passes only through annotated non-canonical introns (T2b fails): the "signature" is a known RefSeq feature.
- T4 fails: whole-chain reads exist and the assembler still does not emit complete copies (a method loss).

## 4. Part B — the 2×2 variant/truncation simulation, mirrored from the NPIP prereg

Everything below is the NPIP prereg's design with the TBC1D3 truth, contigs and real depths substituted. Code =
copies of the sibling's frozen scripts (verified against its sha1s: `callvar.py` 592846d1, `build_reads.py` d000f3a8,
`common.py` 52a2b6fa, `recurrent.py` 382864df, `map_batch.sh` be09bd0f, `spike.sh` 8d459c39, `arm.sh` c6ed4cbb,
`trace.py` a556b01f, `polish_port.py` 854c4d78, `arms_score.py` 970356be, `attribute.py` d0b6b5f7, `noncanon.py`
4b7d3eff, `dediv.py` 48dae80d) with path-only edits (§6).

- **Step 1 — the individuals' copy sequences:** `callvar.py` unchanged: own primaries (primary, orientation = copy
  strand, ≥ 1 aligned base on the territory, no MAPQ filter); pileup over =/X/M/D; call at coverage ≥ 4 when an
  allele has ≥ 3 reads and > 50%; SNV / deletion / insertion; `N` gaps never called. Reported per copy: calls by
  class, near-splice calls, the callable fraction, and R_TBC explainability.
- **Step 2 — arms**, n = the copy's real own-primary count, read *i* paired across arms: **C** reference, full
  length (0-30 bp end jitter); **T** reference, 3′-anchored, length = the aligned query length of the copy's *i*-th
  real own primary (seeded permutation); **V** variant copy, full length; **VT** variant copy, real lengths. HiFi
  errors sub 0.001 / indel 0.0003. Mapping `minimap2 2.30 -ax splice:hq -uf --eqx -Y -N 50 -p 0.1 --secondary=yes
  -t 4` on `winloci_data/GGO.splice.mmi` (whole genome), batches of ≤ 1,200 reads under the heavy lock.
- **Step 3 — spike-in assembly:** background = the sample's real records on the two contigs minus every record of
  every own-primary read; plus the arm's simulated records; as_table = real genome-wide table + the arm's simulated
  table. Arm **E** = background only. Runs ship / nopol / maj (`arm.sh`).
- **Supplementary, reported SEPARATELY (the brief's zero-depth rule):** copies with 0 real own primaries in a sample
  get 0 reads in C/T/V/VT (as NPIP). Separately, each such copy gets a **fixed depth of 10 reads** (`build_zreads.py`):
  ZC = reference full length, ZT = reference, 3′-anchored, length drawn from the sample's pooled own-primary lengths.
  No variant is callable at a 0-read copy, so ZV = ZC and ZVT = ZT. Spike-ins CZ = C + ZC, TZ = T + ZT, VZ = V + ZC,
  VTZ = VT + ZT (`spikez.sh`); only the zero-read copies' presence / first losing step is reported from them. They
  never enter the verdict.
- **Step 4 — scoring:** the same code per arm: presence (`arms_score.py`), first losing step (`trace.py` →
  `polish_port.py` → `attribute.py`), completeness, `noncanon.py`, `dediv.py`, `recurrent.py` (R_TBC exact + the NPIP
  classes, expected 0 here), `recurrent_tbc.py`, placement. Families are NOT run (secondary in NPIP; the question
  here is loss).

**R_TBC and J2 for TBC1D3 (the one design element that cannot be copied).** NPIP's J2 names two NPIP-specific motif
classes. Here R_TBC = the NPIP rule verbatim (`common.R_junctions()`: a non-canonical junction carried by ≥ 1 read of a
≥ 2-read skeleton on a copy territory in BOTH libraries), classes = R_TBC by +-strand motif (upper-cased) with band
[min length − 100, max length + 100] (NPIP's hand bands lie inside this rule's bands), J2 classes = classes with ≥ 2
distinct R_TBC junctions (else all classes), **J2 = every J2 class is carried by ≥ 2 reads at the copies.** R_TBC
empty ⇒ J2 not evaluable, J judged on J1 alone and flagged.

**Decision rules (NPIP's, both samples; `verdict.py` a3467c80):**
- J1: non-canonical read fraction within real ± 0.15. J = J1 ∧ J2.
- L1: |absent − real absent| ≤ 3 over the 14 copies. L2: half the L1 distance between the step-category counts
  {none/seeding/dedupe, pass-1 floor, gate, polish, other} ≤ 3. L = L1 ∧ L2. ⚠ 3 of 14 is looser than 3 of 25; the
  size-scaled bar (≤ 2) is reported next to it (not decisive).
- **Power (added, fixed now):** L is NO POWER in a sample whose real data has < 3 absent copies with ≥ 1 own primary;
  L is judged on the powered samples only.
- Verdict mapping (NPIP's table; my deterministic reading, written into `verdict.py` before any arm exists):
  sequence sufficient = V has J and L; interaction = only VT has J and L; truncation sufficient = T has L and V does
  not; structural / not explained = neither V nor VT satisfies J2. L and J verdicts are also reported separately
  (L: C has L ⇒ "depth / isoform spread", then T-not-V ⇒ truncation, V ⇒ sequence, VT only ⇒ interaction; J: V ⇒
  sequence, VT ⇒ interaction, neither satisfies J2 ⇒ structural). If the sibling's Outcome reads the NPIP table
  differently, both readings are given.

**Predictions (committed now):**
- C: non-canonical fraction ≤ 0.05; J2 absent; ≤ 2 absences among copies with real reads per sample; complete ≥
  half of the copies with reads.
- T: non-canonical ≤ 0.05; J2 absent; complete ≤ C's; L passes in ≥ 1 powered sample.
- V: non-canonical < 0.15 (J1 fails); J2 absent; |absent(V) − absent(C)| ≤ 2.
- VT: close to T (|absent(VT) − absent(T)| ≤ 2); J1 and J2 fail.
- Arm E: presence at < 5 copies (NPIP's falsifier bar; ≥ 3 reported as the size-scaled flag).
- **Predicted verdict = NPIP's predicted verdict: loss "truncation sufficient"; junction signature "structural / not
  explained".** ⚠ TBC1D3 transcripts are ~2.1 kb (NPIP's up to 8.9 kb), so real reads may be near full length and T ≈ C;
  that would make this prediction fail, which is the point of stating it.
- Supplementary Z: a zero-read copy given 10 reference full-length reads is present in CZ for ≥ 80% of such copies.

**2×2 transfer:** TBC1D3's L- and J-verdict labels equal the NPIP study's OUTCOME labels (compared when both exist;
if the NPIP outcome is not available at write-up, the comparison is marked pending and only the prediction is
judged).

**2×2 falsifiers (NPIP's, transposed):** V's non-canonical fraction ≥ real − 0.15 in both samples (junctions are an
aligner response to point/small-indel divergence); V or VT produces J2 (the recurrent junctions are sequence-driven);
C reproduces L (the loss is depth and isoform spread); arm E gives presence at ≥ 5 copies (read the per-copy
comparison net of E).

## 5. Hostile self-review

1. **Near-identical copies change the mechanism space.** TBC1D3 is a 0.997 clique (r810). Reads may tie across copies
   (MAPQ 0), so a copy's "own primaries" are partly an aligner tie-break, and absence may be a PLACEMENT effect (O2)
   rather than divergence. T1 would catch it as losses at seeding or as absences with ties; the 2×2 reproduces
   placement ties by construction (simulated reads are placed by the same aligner).
2. **Variant calls are blends.** With near-identical copies, the reads at one copy come from several; the called
   "individual copy" is a majority blend, hets are dropped. More acute than for NPIP.
3. **The annotated 43-bp non-canonical introns** (nine models) are RefSeq indel corrections. If the individuals'
   reads carry them, T2 could pass on a known reference feature; T2b guards this. The simulation merges them (as NPIP),
   and `N` gaps are never called, so V cannot reproduce them: "structural" is then partly by construction; R_TBC is
   split into annotated / unannotated in the report.
4. **Known prior:** the 07-16 testis census (~100 reads/locus × 3). OR6737 likely has few absences with reads; T1 and
   L may have NO POWER there. The power rules were fixed now for that reason.
5. **Different truth construction from NPIP.** NPIP's territory is the T_member landing (liftoff + identity), with
   lifted exons; TBC1D3's is the annotated record alone. "Own read" therefore differs slightly in definition.
6. **14 copies, 2 of them atypical** (a 423-bp 3-junction 3K-like model next to a 4-junction pseudogene). Thresholds
   copied from a 25-copy family are looser here; scaled bars are reported.
7. **Held-out family, not held-out library.** Same two individuals, same libraries, same aligner and reference as the
   NPIP study. A pass shows the mechanism is not NPIP-specific in these libraries, not that it generalises across
   individuals or technologies.
8. **J2 class rule is new** (NPIP's classes were hand-set after seeing R). The rule is fixed now, before R_TBC exists.
9. **My reading of the NPIP verdict table is mine** (§4); the sibling has not frozen a verdict script.
10. **n = 2 samples, dev only, thresholds are mine.** Sufficiency, not necessity: an arm reproducing the loss shows a
    sufficient cause in this model, not how the real data arose. minimap2 2.30 (sim) vs 2.31 (real BAMs).
11. Metric traps checked: the denominators (own spliced primaries, all 14 copies) are fixed before scoring and are not
    conditioned on presence; the scored set is all 14 records, chosen from annotation; records are addressed by
    ID + coordinates, never by name alone.

## 6. Frozen code (sha1, first 8; `code/`, recorded 01:38 before any read is opened)

**Copies with path-only edits** (truth path, study root, contigs, regions file; diffs recorded in `code/orig/`):
- from `ggo_npip_varsim/code/`: `trace.py` e0659be3, `polish_port.py` 854c4d78 (unchanged), `arms_score.py`
  4991bed8, `attribute.py` 64d24e3c, `dediv.py` c2eec75b, `callvar.py` 1eaa05f6, `build_reads.py` ae01b0e0,
  `common.py` 3bb20fdc, `recurrent.py` faac7d78, `map_batch.sh` c257ad55, `spike.sh` 36782967, `arm.sh` 4a0b6ebe;
- from `ggo_npip_loss/`: `arm_diag.sh` 96d07c04, `validate.py` 29e3b8a2, `recur.py` 299161dd, `waterfall.py`
  4060b004, `table.py` c936f0f8;
- `noncanon.py` b497aae5: the truth path, plus ONE non-path edit — a `max(·, 1)` guard on the control denominators,
  because the spike-in BAMs hold no control contig (the real-data result is unchanged).

**New for this study:** `build_truth.py` 34d27c47, `build_reads_lib.py` bdaba758 (build_reads.py's three functions
verbatim), `build_zreads.py` 1014cc32, `spikez.sh` a2c19c0f, `recurrent_tbc.py` 7081ae59, `bytecmp.py` 7032fbd9,
`transfer.py` 7b60b4d3, `verdict.py` a3467c80.

**Annotation products:** `ann/truth.json` 7df0f9de, `ann/tx.tsv` d130163e, `ann/nc_introns.tsv` 0c9f4e75,
`ann/records.tsv` fdbb60b8, `tbc_contigs.txt` 5672cae8.

## Outcome

*(Appended 2026-09-29 ~03:45 PDT. The text above this heading is byte-identical to the frozen version, sha1
0b3e207e, recorded in `FROZEN.sha1` at 01:40 before any read was opened. Everything was run as registered except
the deviations listed in §O5.)*

### O1. Answer first

- **The NPIP hypothesis does NOT transfer to gorilla TBC1D3** in OR6737, the only sample with power. **KB3781 has no
  power**: the fibroblast library has 4 own primaries over the 14 copies.
- **The TBC1D3 loss is real, but its cause is not the NPIP cause.**
  - Where copies with their own reads get no locus (OR: 4 copies), the reads are canonical: 0/15 of those copies'
    spliced primaries carry a non-canonical junction.
  - Instead, 1-6 reads per copy each carry a DIFFERENT intron chain, so the exact-chain pass-1 floor drops them.
- **Completeness is not lost before the assembler.** 6 of OR's 7 present copies are complete (the full 13/15-junction
  primary chain), against 0/50 for NPIP.
- **2×2 (my pre-registered reading of NPIP's table): loss = "truncation sufficient"; junctions = "structural / not
  explained".**
  - The junction verdict is vacuous: there is no recurrent non-canonical junction to reproduce (R_TBC is empty).
  - Pre-registered 2×2 transfer test: TBC1D3's labels vs **NPIP's OUTCOME** (loss reproduced by no arm; junctions
    structural). **The loss label differs, so the 2×2 does not transfer either.** TBC1D3 matches NPIP's *prediction*,
    not its outcome.
- **Zero-depth supplement: every zero-read copy is present and complete with 10 reference reads** (OR 5/5, KB 10/10,
  in all four Z arms). This includes the pseudogene LOC115931404 and the 423-bp 3K-like model.

### O2. Part A — the diagnostic on TBC1D3 (real BAMs)

**Instrument checks: all pass.**
- I1: the region-restricted `ship` run equals the 09-25 genome-wide (unrestricted) GTF line for line, cov/TPM/FPKM
  removed, on BOTH contigs in BOTH samples: OR 65,989 + 91,500 lines; KB 64,107 + 89,796 lines.
- I2: `validate.py` gives 14/14 copies identical in strict, majority and read-strand modes, in both samples.
- I3: `polish_port.py` is ALL IDENTICAL for all three polish pairs in both samples.
- Contig mono floors: OR 14 / 14, KB 17 / 17 (NC_073224.2 / NC_073228.2).

**Per copy, shipped command** (`diag/score/table.md`; prim / kept sec / all same-orientation sec → first losing step):

| sample | present | complete | absent: pass-1 floor | absent: seeding | other steps |
|---|---|---|---|---|---|
| OR | 7 | 6 (LOC101125558, LOC129533808, LOC129533806, LOC115934662, LOC129533797, LOC129533813) + LOC115933306 at 12/14 junctions | 4 (LOC101151653 1/2/75, LOC129533458 2/1/74, LOC101144080 6/0/85, LOC129533792 6/0/89) | 3 (LOC109026840 0/0/13, LOC115931404 0/0/5, LOC134759231 0/0/4) | 0 |
| KB | 0 | 0 | 8 (1-2 seeded reads each; 4 primaries in total) | 6 | 0 |

- **Present copies are held by their own reads, with no merges.** No holder holds a second copy, including the tandem
  pairs 38 kb apart. Two of OR's present copies have 0 primaries (LOC129533806, LOC129533813) and are held by 20
  seeded secondaries each: near-identical copies share reads.
- **Why OR's four own-read absences hit the pass-1 floor** (`trace.OR.json` single chains). Every read is its own
  chain, and every chain passes the strict junction clause:
  - LOC101151653 / LOC129533458: fragments with 12 and 3 junctions;
  - LOC129533792: 6 reads carrying 5-13 junctions (0-2 unannotated);
  - LOC101144080: 3 of its 6 reads carry 33-34 junctions, 26-32 of them unannotated (chains that run past the
    copy).
- **Single switches** (polish on): all-secondaries +4 and floor-1 +4 (the same four copies, OR); all-secondaries +2 (KB:
  LOC115933306 through 1,054 low-AS secondaries, LOC129533797 through 65). Majority, keep-duplicates, read-strand,
  no-mono-floor and no-polish: +0. Both effective switches are closed rows (r1060/r1100; r1067-1073/r1082/r1083).
- **Waterfall** (OR): 90 primaries + 123 of 1,317 secondaries seeded → 205 after dedupe → 109 in 1-read chains →
  96 in ≥ 2-read skeletons (all spliced) → 96 same-strand gate survivors (**0 removed by strict junctions, 0 by the '+'
  placeholder**) → 64 in shipped transcripts. KB: 9 seeded, all 9 in 1-read chains.

**Transfer clauses** (`transfer.py`, `diag/score/transfer.json`):

| clause | OR6737 | KB3781 |
|---|---|---|
| T1 same three steps (≤ 1 exception, power ≥ 3) | **FAIL**: 7 absences with reads; 4 pass-1 floor, **3 seeding** | FAIL: 8 pass-1 floor, 6 seeding |
| T2 non-canonical ≥ 0.30 and ≥ 10× control | **FAIL: 0.022 (2/90) vs control 0.029 (0.8×)** | NO POWER (3/4; n < 20) |
| T2b the same without the annotated indel introns | FAIL: 0.022 (0 reads carry an annotated 43/31-bp intron) | NO POWER |
| T3 canonical-only `de` ≥ 2× control | **FAIL: 0.0027 vs 0.0015 (1.8×)** | NO POWER (n = 1) |
| T4 ≥ 80% no whole-chain read and ≤ 1 complete | **FAIL: 7/14 without a whole-chain read, 6 complete** | PASS (trivially: ≤ 2 reads per copy) |
| per-sample verdict | **DOES NOT TRANSFER** | **NO POWER** |

- T6 (descriptive, the direction NPIP implies): the absent copies' spliced primaries are 0/15 non-canonical; the
  present copies' are 2/75.
- Junction instances: non-canonical 0.42% (control 0.32%); annotated 65% (NPIP 9% / 25%).
- **R_TBC is empty.** No non-canonical junction recurs in both libraries (`recur.py`, `recurrent_tbc.py R`), and 0
  reads carry the NPIP classes.
- KB's 4 own primaries: 3 carry ≥ 1-kb non-canonical gaps at `de` 0.070, 5-50× the family's own divergence. Most
  likely they are foreign or unreferenced molecules, but at n = 4 this is anecdote.

**Predictions:** "TRANSFERS in both samples" — **refuted in OR** (T2, T2b, T3, T4 all fail) and untestable in KB.
The only NPIP feature TBC1D3 shares is the step class of the low-depth losses (pass-1 floor) and the fact that every
rescuing switch is a closed row.

### O3. Part B — the 2×2 (`score/verdict.json`, `score/verdict.out`)

**Step 1, the individuals' copies.**
- OR: 33 calls (30 SNV, 3 insertions, 0 deletions) at 5 copies:
  - LOC129533797 17, LOC115933306 8, LOC129533792 5, LOC101144080 2, LOC115934662 1;
  - 17 exonic, 4 within 20 bp of an annotated boundary, 0 in a splice dinucleotide;
  - callable fraction 0.55-1.00 at those copies, 0 elsewhere.
- KB: 0 calls (≤ 1 read per copy). **KB's V and VT reads are therefore identical in sequence to C and T.**
- R_TBC is empty, so no explainability table.

**Reads and mapping.**
- 376 arm reads (OR 90 × 4, KB 4 × 4) + 300 Z reads = 676, mapped in one batch (59 s, 15.0 GB).
- Primaries on the source copy: C 90/90, T 87/90, V 64/90, VT 62/90. So 26-28 V reads land on ANOTHER copy: the
  called "individual copy" is partly a blend of copies (hostile review 2).
- Canonical-only `de`: C 0.0014, T 0.0014, **V 0.0024, VT 0.0027**, real 0.0027. The called point variants account
  for the whole of TBC1D3's mild excess divergence.

| arm | OR absent (real 7) | L1 / L2 dist | OR steps floor/early | OR complete (real 6) | KB absent (real 14) | non-canonical OR / KB (real .022 / .750) | L | J |
|---|---|---|---|---|---|---|---|---|
| C | 3 | 4 ✗ / 2.0 | 0 / 3 | 11 | 14 (L1 0, L2 1.0) | 0 / 0 | ✗ | ✗ (KB J1) |
| T | 5 | 2 ✓ / 1.0 ✓ | 2 / 3 | 8 | 14 | 0 / 0 | **✓** | ✗ (KB J1) |
| V | 3 | 4 ✗ / 2.0 | 0 / 3 | 11 | 14 | 0 / 0 | ✗ | ✗ |
| VT | 5 | 2 ✓ / 1.0 ✓ | 2 / 3 | 8 | 14 | 0 / 0 | **✓** | ✗ |
| E | 14 | 7 ✗ | — | 0 | 14 | — | ✗ | — |

- **L power:** both samples meet the registered power rule (4 absent copies with ≥ 1 own primary each).
- **Scaled bar** (≤ 2): T and VT still pass in OR (distance 2).
- **J2 is not evaluable** (R_TBC empty), so J = J1:
  - J1 passes in OR in every arm (reference reads match an absent signal);
  - J1 fails in KB in every arm (0 vs 3/4).

**Per copy (OR):**
- T and VT lose LOC101151653 (1 primary) and LOC129533458 (2) at the pass-1 floor, as the real data does.
- No arm loses LOC101144080 or LOC129533792 (6 primaries each). Their real reads carry unannotated junctions and ends
  that annotated-isoform reads cannot have (NPIP hostile review 4).
- The 3 seeding copies are absent in every arm, E included.

**Verdict** (`verdict.py`): L = **truncation sufficient**; J = **structural / not explained (on J1; R_TBC empty)**;
combined = **truncation sufficient**. The C falsifier did not fire (C lacks L). V did not produce J. E is 0 present in
both samples (bar < 5).

**Predictions vs outcome:**

| prediction | outcome |
|---|---|
| C: NC ≤ 0.05; J2 absent; ≤ 2 absences among copies with reads; complete ≥ half | ✓ 0 / 0; n/a (R empty); ✓ OR 0 of 9, **✗ KB 4 of 4** (1-read copies); ✓ OR 11, **✗ KB 0** |
| T: NC ≤ 0.05; J2 absent; complete ≤ C; L in ≥ 1 powered sample | ✓; n/a; ✓ 8 ≤ 11; ✓ (both) |
| V: NC < 0.15; J2 absent; \|V − C\| ≤ 2 | ✓ 0; n/a; ✓ 0 / 0 |
| VT: \|VT − T\| ≤ 2; J1 and J2 fail | ✓ 0 / 0; J fails, **but only through KB's 4 reads** (J1 passes in OR) |
| E < 5 present | ✓ 0 / 0 |
| verdict = NPIP's predicted (truncation sufficient; structural) | ✓ as labels; the J label is vacuous |
| Z: ≥ 80% of zero-read copies present in CZ | ✓ 15/15 copy-samples, all complete, in CZ/TZ/VZ/VTZ |
| **2×2 transfer: labels equal NPIP's outcome** | **✗ on the loss** (NPIP: not reproduced by any arm; TBC1D3: truncation sufficient); J label equal but vacuous |

- **Z side effect:** in KB, the Z reads also give LOC129533458 and LOC101125558 a locus in every Z arm. Their 1-read
  chains join the Z reads' chains through seeded secondaries. This is the near-identical-copy read sharing again.

### O4. What this says about the NPIP hypothesis

- **On the held-out family the mechanism is family-specific.** At TBC1D3 the absences with reads are
  depth-and-chain-diversity losses at the exact-chain pass-1 floor:
  - 1-6 canonical reads per copy;
  - different 5′ truncation points and unannotated extensions;
  - in a 0.997-identity family whose copies share reads through seeded secondaries.
- They are not losses to divergent reads or non-canonical junctions. So the NPIP explanation — divergent reads with
  non-canonical junctions, an alignment/reference (O3) question — is **NPIP-specific**. It should be quoted as
  "NPIP (and not TBC1D3)", not as a general account of multi-copy loss.
- **What does generalise:** in both families, every switch that restores a copy is a closed register row (all
  secondaries, floor 1).
- **The present TBC1D3 copies show the method is not the limit when reads are full-length:** 6/7 complete in OR, and
  15/15 zero-read copies complete at 10 reads. This matches the NPIP FL sim, now on real reads.
- **KB3781 carries no information** about either mechanism for TBC1D3 (4 reads).

### O5. Deviations (all recorded before the verdict was computed, except where stated)

1. **14 records, not 15** (§2). The brief's count was off by one; no record was dropped.
2. **KB `floor1` was not run.** `copy_assign` reached 23 GB RSS on NC_073228.2 alone and the machine had 0 GB
   available. The process was killed by PID after a `readlink /proc/<pid>/cwd` check (rc 143).
   - KB's single-switch column therefore lacks floor-1.
   - The trace's floor-1 gate emulation says only 2 of KB's 8 floor-lost copies have a single read that would pass
     strict (LOC101151653, LOC109026840).
   - OR `floor1` / `floor1_nopol` ran (245 / 113 s, 14.5 GB). Their `arms_score` entries were merged in through
     `arms_score_or.py` (the sample loop limited to OR).
3. **`spike.sh` / `spikez.sh` contig order** (sha1 37c7dddb / 1f0c6630).
   - `samtools view` region extraction keeps argument order, so `NC_073228.2 NC_073224.2` produced unsorted spike-ins.
   - This showed up only where simulated reads map to NC_073224.2 (T/VT failed at `samtools index`).
   - Fixed to header order; T/VT rebuilt. C/V had no NC_073224.2 records and were valid.
4. **`verdict.py`**: `nc.get('reads', 0)` (sha1 ea349ec3), because arm E has no own reads and the key is absent. The
   crash happened AFTER the four arms' L/J lines had printed; the logic is unchanged.
5. Spike-in BAMs and per-arm table copies were deleted after scoring (disk at 97%). They can be rebuilt from
   `aln/sim.bam` + `spike/*.bg.bam` with `spike.sh` / `spikez.sh`.
6. The copied scorers keep NPIP labels (`NPIP_pooled`, `npip_can`); here they mean the 14 TBC1D3 copies.

**Post-hoc readings (not registered, labelled as such):**
- T1 fails only on the 3 seeding copies, whose 4-13 secondaries all belong to reads placed ≥ 2% better elsewhere (0
  primaries). Read as unexpressed, OR's four own-read absences are 4/4 at the pass-1 floor. That would pass T1, but
  T2-T4 fail regardless.
- KB meets the L power rule by its letter, but every arm reproduces KB exactly, because ≤ 2 reads per copy cannot pass
  the floor in any arm. The rule should have required reads per copy, not copies with reads.

### O6. Register rows (appended 2026-09-29 to `docs/NEGATIVE_RESULTS_REGISTER.md` as 1155-1157)

| N | date | area | claim | verdict |
|---|---|---|---|---|
| 1155 | 2026-09-29 | multi-copy loss (held-out family) | The gorilla NPIP loss mechanism — copies with reads get no locus because the individual's reads are divergent and carry non-canonical junctions — holds for gorilla TBC1D3 (PREREG_ggo_tbc1d3_holdout_2026-09-29, sha1 0b3e207e) | ⛔ **Does not transfer (OR6737; KB3781 no power, 4 own primaries).** OR: non-canonical 2/90 = 0.022 vs control 0.029 (NPIP 0.70); canonical `de` 1.8× control (NPIP 3.9×); 6/7 present copies complete (NPIP 0/19); R_TBC empty. The 4 own-read absences are canonical 1-6-read copies whose reads each form a distinct chain (pass-1 floor), in a 0.997 family that shares reads through seeded secondaries. The NPIP account is NPIP-specific. Frozen instruments reproduce the genome-wide GTF byte for byte. |
| 1156 | 2026-09-29 | multi-copy loss (2×2 sim, held-out family) | The 2×2 variant/truncation verdict for gorilla NPIP (loss reproduced by no arm; junctions structural) is the same for TBC1D3 | ⛔ **Loss label differs.** TBC1D3: T/VT reproduce L in OR (absent 5 vs real 7, L2 1.0; C/V 3), so "truncation sufficient". They reproduce 2 of 4 pass-1 losses (the 1-2-read copies); the 6-read copies carry unannotated structure no arm simulates. The J label ("structural") is vacuous: R_TBC is empty and J1 fails only on KB's 4 reads. V reproduces the real divergence (de 0.0024-0.0027 vs 0.0027). |
| 1157 | 2026-09-29 | multi-copy loss (depth) | Unexpressed TBC1D3 copies (0 own primaries) would stay absent even with modest expression | ⛔ **No: 10 reference reads give presence and completeness at 15/15 zero-read copy-samples** (OR 5, KB 10, all four Z arms), including a pseudogene and a 423-bp model. Their real absence is expression, not method. |
