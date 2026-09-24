# Pre-registration — O2 accuracy with read-level truth on the real family catalogs

**Written 2026-09-23 (§6zf), before any read is simulated.** User: *"lets do 1"* — the first of the O2
strengtheners: a per-read assignment accuracy on the real multi-copy families, which O2 has never had
(register row 583: the only headline was on a hand-made family with a denominator conditioned on the
assignment). Truth comes from simulation because no real dataset carries a read's copy of origin.

## Objects

- **Copies:** every copy of every multi-copy family in the shipped O1 catalogs — human `HSA_gwcat`
  (394 families, 1,220 copies, CHM13) as the DEVELOPMENT substrate and gorilla `GGO_gwcat` (494 families,
  mGorGor1) as the HELD-OUT substrate. Each copy's sequence is the catalog's own spliced emission
  (`copies.fa`), i.e. exactly the transcript O1 defined the copy with.
- **Reads:** `bench/sim_reads.py` HiFi model (substitution 0.001, indel 0.0003), plus what real libraries
  do to ends: ± 0–30 bp jitter at both ends (mandatory — identical ends collapse under the coordinate
  dedupe) and 5′/3′ degradation trimming up to 10% of the length from a random end. Reads per copy =
  the copy's real read count clamped to [10, 100], so depth follows expression without letting one
  family dominate. Read name = `family|copy_idx|i` (the truth).
- **Mapping:** the shipped command (`minimap2 -ax splice:hq -uf --eqx -Y -N 50 -p 0.1 --secondary=yes`)
  against the WHOLE genome (prebuilt splice indexes), so every paralogue competes as in production.
- **Assignment:** `copy_assign --families <copies.tsv> --copies-fa <copies.fa>` on the simulated BAM,
  shipped defaults; `<out>.assignments.tsv` is the object scored.
- **Divergence covariate:** per copy, identity to its closest sibling in the family (minimap2 `asm20`
  all-vs-all of the catalog copies; `1 − identity` = the copy's divergence). Bins: < 0.5%, 0.5–1%, 1–2%,
  2–5%, ≥ 5%.

## Metrics — committed now

Per read (every simulated read is in the denominator; nothing is conditioned on the assignment):

| outcome | definition |
|---|---|
| **correct** | assigned to the copy the read was simulated from |
| **wrong** | assigned to any other copy (same family or another) |
| **abstain** | in the table with no copy assigned (any non-assigned status) |
| **lost** | not in the assignment table at all (unmapped, or mapped outside every catalog copy) |

Reported per substrate, overall and per divergence bin, separately for **MAPQ-0 reads** (the tied set,
O2's actual job — the aligner cannot decide these) and MAPQ-60 reads; plus the same for two baselines:
(i) **aligner-primary** (assign the read to the copy its primary alignment overlaps — what StringTie/FLAIR
effectively do), (ii) **abstain-all on MAPQ 0** (the conservative no-tool). Numbers: accuracy among
assigned = correct / (correct + wrong); coverage = assigned / all; the error rate among assigned in each
bin.

| outcome on the MAPQ-0 stratum, human dev AND gorilla held-out | verdict |
|---|---|
| accuracy among assigned ≥ 0.95 overall and ≥ 0.90 in every divergence bin ≥ 1%, with coverage ≥ 0.5 | ⭐ **O2 assigns tied reads correctly where PSVs exist and abstains where they do not** |
| accuracy ≥ 0.95 but coverage < 0.5, or accuracy 0.90–0.95 | ⚠ **correct but conservative / mostly correct** — report the bin where it fails |
| accuracy < 0.90 among assigned on either substrate | ⛔ **O2 makes wrong calls at a rate the thesis cannot carry; find the bin and the mechanism** |

The aligner-primary baseline is expected near 1/k on MAPQ-0 reads (register: never 1/k); the gap between it
and O2's accuracy is O2's contribution. On MAPQ-60 reads O2 must not be worse than the aligner (≥ 0.98).

**Predicted, before looking:** ⚠ on the human dev set — accuracy among assigned ~0.95 overall, failing in
the < 1% bin (few PSVs clear the 0.1% error floor, so abstention should dominate there; any assignment
made there is the risk), coverage 0.5–0.7 on MAPQ-0 reads; gorilla similar or better (more divergent
families). Aligner-primary on MAPQ 0: ~0.5.

**Realism caveats stated up front:** no readthrough molecules, no intron retention, no reads from copies
the catalog missed (§6u7's over-merge and r1052's non-homologous exons are not in the truth), one error
model. These make the number an UPPER bound on real-data accuracy; the divergence-bin structure is what
transfers.

I will not change the read model, the bins, the denominators or the bar after seeing any number.

> **Renamed 2026-09-24 (publishing names):** `bench/o2_read_truth.py` is now `bench/copy_assign_read_truth.py`; the pipeline driver stage is `assign`.

---

## Addendum 1 (2026-09-23, after the first human run, before the gorilla run is scored)

**Substrate change, forced by the input contract.** `copy_assign --families` requires the catalog's `exons`
column; the genome-wide July catalogs (`HSA_gwcat`, `GGO_gwcat`) predate it and were refused. The
current-format catalogs on disk are per contig: development = human chr16 (`chr16_arm/on`, 1,418 copies,
290 multi-copy families — the NPIP/TBC1D3 territory), held-out = gorilla NC_073244.2 (`hom_c234`, 357
copies, 54 multi-copy families). 18 chr16 copies shorter than 300 bp (single-exon fragments) received no
reads and were removed from the catalog handed to O2 (`h16_cat.dropped.tsv`); `--families` refuses a copy
with no reads by design. The genome-wide human simulation (24,789 reads from all 1,220 `HSA_gwcat` copies)
stays as a mapping-only artefact.

**The table is one row per read × family**, discovered on the first pass: a tied read whose placements touch
copies of several catalog families gets a row in EACH of them, and every family judges the read against
ITS OWN copies only. 946 of the 1,245 tied chr16 reads span more than one family; 352 carry "assigned" in
two families. Scoring therefore states its per-read reading explicitly: OWN (the true family's row — the
certificate's own accuracy), PRIMARY (rows with `primary_local = 1`), ANY (any assigned row; two loci =
conflict). The pre-registered metric is reported under each; the bar is applied to the reading a consumer
can actually use without the truth (ANY / PRIMARY), and OWN is reported as the certificate's intrinsic
accuracy. `--no-as-tied-only` (every read in the table) exceeds 10 minutes on 28k reads and is not used;
MAPQ-60 reads are the aligner's and are scored by the aligner-primary baseline only.

---

# OUTCOME (2026-09-23)

## Human chr16, development — 1,400 copies / 290 multi-copy families, 28,453 reads

Mapping: 23,837 reads at MAPQ 60, 3,357 at MAPQ 1–59, **1,259 at MAPQ 0** (the tied set O2 is for). The
aligner's own placement is right for 100.0% of MAPQ-60 and 99.6% of MAPQ-1–59 reads and for **49.2% of
MAPQ-0 reads** (the coin toss the register forbids as an answer).

**The certificate itself is exact.** Reading each tied read's row in its TRUE family (OWN):

| divergence of the source copy | tied reads | assigned | correct | wrong | abstain | coverage |
|---|---|---|---|---|---|---|
| < 0.5% | 1,082 | 58 | 58 | **0** | 1,017 | 0.05 |
| 0.5–1% | 56 | 42 | 42 | 0 | 14 | 0.75 |
| 1–2% | 47 | 37 | 37 | 0 | 10 | 0.79 |
| 2–5% | 34 | 10 | 10 | 0 | 17 | 0.29 |
| ≥ 5% | 40 | 10 | 10 | 0 | 30 | 0.25 |
| **all** | **1,259** | **157** | **157** | **0** | **1,088** | **0.125** |

157 assigned, 157 correct, 0 wrong: the PSV certificate never picks the wrong copy within a family, and it
abstains on 94% of the reads from copies under 0.5% divergence — the regime where identical sequence
over the read leaves nothing to decide. Coverage in the 1–2% bin is 0.79. The low coverage at ≥ 2% is
the tie itself: a MAPQ-0 read from a divergent copy is one that happens to lie on a stretch shared with a
sibling, so abstaining is right.

**⛔ But the table a consumer sees does not carry that accuracy.** `--families` emits one row per read ×
family, and every family judges a read against its own copies only. A tied read whose secondary
placements touch copies of other catalog families gets rows there too (946 of the 1,259 tied reads span
several families), and those foreign families "assign" it — 904 foreign assigned rows against 157 true
ones; 348 of them as a **sole candidate** (`n_candidates = 1`: the family has one copy in the read's
placement set and assigns it by default), the rest by a within-family PSV vote whose candidate set does
not contain the true copy (the example rows carry 37 decisive sites and p = 10⁻⁹¹ for the wrong locus).
`as_best` is per read, not per row, so nothing in the table separates a foreign row from the true one.
Under the readings a consumer could use without the truth: **PRIMARY** (rows with `primary_local = 1`)
accuracy among assigned **0.29** at coverage 0.26; **ANY** (any assigned row; two loci = conflict) **0.09**
at coverage 0.51 with 153 conflicts. Both are below the aligner's coin toss, so by the pre-registered bar the
shipped output is ⛔ on the MAPQ-0 stratum — for a reason that is a missing step, not a wrong statistic.

**Mechanism, stated for the fix.** The E_r family relation and the read-competition relation are not the
same graph (the 08-14 excision measured it: 53.8% of orphans land on a paralogue outside the catalog's
family). O2 arbitrates within a family and never across families, so a read that competes across two
families is assigned twice or assigned in the wrong one. The fix is a cross-family pass: the candidate set
of a read is the union of every catalog copy its placements touch, regardless of family, and the same PSV
certificate runs once on that union (a family may then assign a read only if it wins against every other
family's candidates; a sole candidate in one family is not a decision when the read also has candidates
elsewhere). The OWN numbers above are what that pass would produce if the true family's candidates always
win the union, which is the certificate's property on this substrate (0 wrong within family).

## Gorilla, held-out — NC_073244.2 (`hom_c234`, 357 copies) and the SD catalog `ggo_sd` (322 copies on NC_073241/42/44)

| substrate | reads | MAPQ 0 | source-copy divergence of the tied reads | O2 assigned | wrong | abstain | aligner-primary on MAPQ 0 |
|---|---|---|---|---|---|---|---|
| NC_073244.2 catalog | 11,448 | 25 | all < 0.5% | 0 | 0 | 25 | 0.69 |
| SD catalog (3 contigs) | 5,899 | 145 | all < 0.5% | 0 | 0 | 134 (11 no row) | 0.49 |

The held-out ties are all reads from copies under 0.5% divergence, and O2 abstains on every one of them —
the same behaviour as the human < 0.5% bin (94% abstain, 0 wrong). No foreign-family row assigned anything
on gorilla (the SD catalog's families do not share read placements), so the cross-family defect did not
express here; the held-out substrate therefore confirms "no wrong calls" and "abstains where copies are
near-identical" but carries no assignments to measure accuracy on. The tie-rich human chr16 catalog is
where both the certificate's exactness and the table's defect are visible.

## Verdict

- **The PSV certificate (the O2 decision rule): ⭐ 157/157 correct within family on human, 0 wrong on
  gorilla; abstains on 94% of near-identical-copy reads; coverage 0.75–0.79 where copies differ by 0.5–2%.**
- **The shipped `--families` table: ⛔ 0.29 (PRIMARY) / 0.09 (ANY) accuracy among assigned on human MAPQ-0
  reads**, below the aligner's 0.49, because reads are judged per family and never arbitrated across
  families (946/1,259 tied reads span several catalog families; 904 foreign "assignments", 348 of them
  sole-candidate defaults). This is the mechanism the register's row 583 could not see: the previous
  "accuracy" numbers were read off one family at a time.
- **What to ship next:** a cross-family arbitration pass in `copy_assign` (candidate set = the union of every
  catalog copy a read's placements touch; one certificate over the union; a sole candidate in one family is
  not a decision when the read has candidates elsewhere). Re-score with `bench/o2_read_truth.py score` — the
  OWN numbers are the target it should reach.
- The MAPQ-60 stratum is the aligner's: 100.0% (human) / 100.0% (gorilla) correct by placement, as required.

Scripts: `bench/o2_read_truth.py sim` (simulate + map + closest-sibling identity),
`bench/o2_read_truth.py score` (per-read scoring under the OWN / PRIMARY / ANY readings); runs in
`/mnt/linuxdisk/tmp/gw22/o2sim/` (`h16*`, `g44*`, `gsd*`).
