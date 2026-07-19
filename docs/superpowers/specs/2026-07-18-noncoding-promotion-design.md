# Non-Coding-Aware Promotion of Reference-Absent Copies — Design Spec

**Status:** approved (brainstorming), ready for implementation plan
**Date:** 2026-07-18
**Builds on:** the genome-wide reference-absent funnel (`bench/promote_genomewide.py` → `gw_reference_absent_copies.json`; the `copy_vs_allele_structural.py` / `o4_diploid_validate.py` discrimination that writes `gw_discriminated.json`). Related: the O4 reference-absent catalog and the copy-graph "overhang" work.
**Objective:** recover the credible **non-coding / novel** reference-absent collapses that the current promotion funnel silently drops because it gates on protein homology (BLASTx) rather than collapse quality.

## Goal

Add a **new, isolated, purely-additive** promotion path that promotes a reference-absent collapse onto a separate **non-coding track** when it clears a **collapse-quality bar** (divergence band + genome coverage + single/own-locus + balanced co-segregation), with **protein homology and ORF demoted from gates to labels**. The output is a new file `gw_noncoding_copies.json`; no existing output is regenerated or changed.

## Motivation (grounded, from the 2026-07-18 funnel-gap workflow `wf_2ae48d9c-212`)

The promotion funnel is a **protein-translatability gate, not a collapse-quality gate**:

- **0 of 31** protein-`NO-HIT` reference-absent records were ever promoted (all 73 promoted carry a BLASTx hit).
- Yet the `NO-HIT` records are **not weaker collapses**: NO-HIT median divergence 11.2% ≈ protein-HIT 10.2%; genome_id / genome_cov are comparable. **The only variable that separates promoted from dropped is ORF length** (median 220 aa promoted vs 101 aa dropped) — i.e. translatability, not evidence of a real collapse.
- Applying a collapse-quality criterion (divergence 5–27%, genome_id ≥ 0.74, genome_cov ≥ 0.90) flags **~19 of 31** NO-HIT as credible; **2** are high-divergence repeat/chimera artifacts (genome_id ~0.26–0.29); the rest are weak. Of the **9 strong** NO-HIT (alt_reads ≥ 20 AND alt_cols ≥ 20), **7 are credible, 2 artifact**. **Zero of the 7 are annotated protein-coding paralogs** — 2 lncRNA, 2 unannotated-intergenic (one a 1028-aa novel-ORF candidate), 3 divergent non-coding copies inside protein-coding gene bodies. None could ever produce the reference protein hit the funnel requires.
- The flagship is `NC_073236.2_139051025`: an unannotated intergenic non-coding (lncRNA-like, ORF 102 aa / 9.2% of a 3342-bp transcript) transcript with a ~7% divergent second copy that survives the full adversarial screen (not editing — 31% transversions + fixed ~90% differences; not chimera — one clean mapq60 block; not repeat/multimap — maximal entropy, flat depth CV 0.067; not an ordinary allele — 7% is 30–40× gorilla heterozygosity).

This vindicates the thesis stance — RNA-level multi-copy families defined from read-conflict/PSV topology, **not** from reference annotation or protein homology. The protein gate is structurally blind to exactly the class the RNA method is built to detect.

**Honesty caveat (binding on the design):** "credible" here is an **RNA-only** criterion. At RNA level a reference-absent *duplicated copy* (parCN > 2) and a *hyperdivergent single-copy allele* (parCN = 2) produce an **identical** signature (single reference homolog, unique mapping, two co-segregating haplotypes). Copy-vs-allele is **not resolvable from RNA** and requires DNA parCN. Therefore every promoted record is a **flagged reference-*divergent* candidate**, never a confirmed copy.

## Design

### §1 — Architecture

A new standalone script `bench/promote_noncoding.py`, structured as a **pure classifier + an I/O shell**:

- `classify_noncoding(rec) -> (promote: bool, call: str, reason: str)` — a **pure function** over a plain dict of per-candidate signals (genome_id, genome_cov, own_locus, alt_cols, alt_read_fraction, alt_reads, n_primary, orf_aa, cons_len). No I/O. This is the unit-tested seam.
- an I/O shell (`main`) that loads inputs, runs the single genome re-score, builds each `rec`, calls `classify_noncoding`, attaches labels, and writes the outputs.

The script imports the existing `orf_aa` helper (or a copy) from `promote_genomewide.py` for the ORF label, and reuses `copy_vs_allele_structural.py`'s SEDEF loader pattern for the segdup-partner label. It **does not** import or trigger any consensus-building or coding-catalog code.

### §2 — Inputs (all already on disk; no rebuild)

- `cons.fa` = `/home/juanfra/winloci_scratch/refabsent/gw_promoted/cons.fa` — **all 734 consensuses built before the ORF gate** (only 145 survived into `gw_reference_absent_copies.json`; the other 589 are where ORF<80 non-coding copies live). These are already **post the A→I-editing filter** (they were only written if the flag passed the `agtc < 0.50` & ≥5-subtype spectrum test in `promote_genomewide.py`), so editing-cleanliness is inherited, not recomputed.
- the flags JSON = `genomewide_flags_new.json` (default; overridable arg) — per-cid `chrom, start, end, n_primary_reads, n_alt_positions (=alt_cols), n_alt_reads, alt_read_fraction`.
- `GGO.fasta` (`/home/juanfra/winloci_scratch/GGO.fasta`, indexed) — genome for the re-score.
- `GGO_genomic.gff` (`/mnt/linuxdisk/home/juanfraitu/winloci_data/GGO_genomic.gff`) — biotype label; **optional** (if unreadable, biotype = `unknown`, never a gate).
- `final.bed` (SEDEF self-align) — segdup-partner label; optional.
- carried protein labels: if a cid is present in `gw_reference_absent_copies.json`, copy its `protein`/`prot_id` as a label; else `protein = "not-tested"` (BLAST binaries are absent; protein is never recomputed and never a gate).
- de-duplication set `gw_discriminated.json` (the ~73 already-promoted **coding** copies): cids present here are **excluded** from the non-coding track. The collapse-quality bar is protein-agnostic, so without this de-dup the track would re-list coding copies that already passed the funnel; excluding them makes the non-coding track exactly "the credible collapses the protein gate MISSED" (matching the ~19 expected yield). Optional; absent → no exclusion.

### §3 — Genome re-score (single process)

One `minimap2 -cx splice:hq --eqx -N1 -t <threads> GGO.fasta cons.fa` (mirrors `promote_genomewide.py:105`; splice preset because the consensuses are spliced RNA). Parse the PAF; per cid keep the best hit by `identity × coverage`:
- `genome_id = matches/blocklen`, `genome_cov = aligned/qlen`, `genome_hit = tgt:pos`, `mapq` (stored as a **label**, not gated).
- `own_locus` = best hit target == the flag's chrom AND `|pos − start| < 200_000` (the existing `own` definition).

Note on multiplicity: with `-N1` the PAF has one hit per query, so we do **not** re-test copy multiplicity here (that is the discriminator's job). "Single/own-locus" means the best hit lands on the candidate's own locus (`own_locus`); reads-independent second-locus evidence is captured by the SEDEF `sedef_partner` label (§5), and `mapq` is stored as a label. This keeps the gate honest about what a `-N1` re-score can and cannot show.

This is a **single foreground minimap2 process** (734 small queries vs the genome) — not `copy_assign`, no whole-genome re-index beyond minimap2's own, WSL2-safe.

### §4 — The collapse-quality bar (`classify_noncoding`)

Promote iff **all** hold (constants named at module top so the recall/precision knob is one edit):

| gate | constant (default) | condition |
|---|---|---|
| divergent, not artifact | `GID_LO=0.60`, `GID_HI=0.97` | `GID_LO ≤ genome_id ≤ GID_HI` |
| full-length, not chimera | `COV_MIN=0.90` | `genome_cov ≥ COV_MIN` |
| single / own-locus | — | `own_locus` is True |
| co-segregation breadth | `MIN_COLS=8` | `alt_cols ≥ MIN_COLS` |
| balanced collapse | `AF_LO=0.15`, `AF_HI=0.60` | `AF_LO ≤ alt_read_fraction ≤ AF_HI` |
| depth | `MIN_ALT=5`, `MIN_DEPTH` (from `hidden_copy_scan`) | `alt_reads ≥ MIN_ALT` and `n_primary ≥ MIN_DEPTH` |

(`AF_LO/AF_HI` gate the overall **alt-read fraction** and are deliberately named apart from the imported per-column balance constants `BAL_LO/BAL_HI` in `hidden_copy_scan`, which measure a different quantity and are not reused here.)

Rationale (workflow): the 7 strong-credible span genome_id 0.744–0.93 / genome_cov 0.91–1.07 / alt_cols 21–82 / alt_read_fraction ~0.32–0.49; both artifacts are genome_id ≤ 0.30. `GID_LO=0.60` (div 40%) sits safely below every credible member and above both artifacts, leaning toward recall as chosen. Records failing only `genome_id ≥ 0.97` are labeled `~REF` (essentially reference) and excluded; records failing `genome_id < GID_LO` are labeled `artifact` and excluded — both logged, not silently dropped.

`call` ∈ {`noncoding-candidate`, `~REF`, `artifact`, `thin`, `not-own-locus`}. Only `noncoding-candidate` is promoted.

### §5 — Labels, not gates

Computed and stored, never gating:
- `orf_aa` (from the consensus) and `coding_potential`: `noncoding` if `orf_aa*3/cons_len < 0.30`, else `coding-capable/novel-protein` (this honestly flags the 1028-aa intergenic candidate `NC_073242.2_18569043` as a possible novel *protein* family rather than an lncRNA — non-coding-track-eligible but noted).
- `biotype` from GFF overlap of the **source** locus: `lncRNA` / `protein_coding-gene-body` / `pseudogene` / `intergenic` / `unknown`.
- `sedef_partner`: a distinct segdup partner locus (reads-independent second-locus evidence) or `null`.
- `protein`: carried label or `"not-tested"`.

### §6 — Output & the honesty rail

`{OUT}/gw_noncoding_copies.json` (+ a `.tsv` mirror with the same columns) — one record per promoted locus:

```
cid, chrom, start, end, track:"noncoding",
genome_hit, genome_id, genome_cov, divergence,
alt_cols, alt_read_fraction, alt_reads, n_primary,
orf_aa, coding_potential, biotype, sedef_partner, protein,
copy_vs_allele:"candidate-DNA-needed",
status:"flagged-reference-divergent-candidate",
reason
```

- **`copy_vs_allele` is the literal constant `"candidate-DNA-needed"` on every record**, and `status` is always `"flagged-reference-divergent-candidate"`. The script discovers and flags; it never asserts a confirmed copy. This is enforced in code (the fields are set unconditionally, not derived), so no future threshold change can accidentally promote a record to "confirmed".
- A stderr/stdout summary tallies: promoted (by biotype), `~REF`, `artifact`, `thin`, `not-own-locus`, with the excluded artifacts named.

### §7 — Determinism & crash-safety

- One `minimap2` process; deterministic given fixed `cons.fa` + genome. No POA rebuild, no `copy_assign`, no background/nohup/waiter/pkill. Outputs under `winloci_scratch/refabsent/gw_promoted/`, not `/tmp`.
- Runs foreground (or harness-tracked background — it is a Python/minimap2 script, not `copy_assign`).
- Missing optional inputs (GFF, SEDEF, protein labels) degrade to `unknown`/`null`/`not-tested`, never an error and never a gate.

### §8 — Testing

- **Hermetic unit tests** on `classify_noncoding` (no minimap2/pysam): flagship-like (genome_id 0.93, cov 1.02, alt_cols 82, alt_read_fraction 0.318, own True) → `(True, "noncoding-candidate", …)`; artifact-like (genome_id 0.28) → `(False, "artifact", …)`; near-reference (genome_id 0.99) → `(False, "~REF", …)`; thin (alt_cols 2) → `(False, "thin", …)`; chimera (cov 0.5) → `(False, "not-... ", …)`; not-own-locus (own False) → `(False, "not-own-locus", …)`. Plus a boundary test at each constant.
- **Honesty-rail test:** every promoted record has `copy_vs_allele == "candidate-DNA-needed"` and `status == "flagged-reference-divergent-candidate"`.
- **Data-gated integration** on the real 734 (`cons.fa`): expect the flagship `NC_073236.2_139051025` promoted; the 7 strong-credible NO-HIT present; ~19 total flagged; the 2 high-divergence artifacts (`NC_073243.2_30120917` div 70.6, `NC_073231.2_39477379` div 73.6) **absent**. Foreground, `winloci_scratch` (WSL2 crash rule).

## Out of scope

- DNA parCN validation (the separate step that finalizes copy-vs-allele; on-disk mGorGor1 mat/pat `.mmi` + SUN density).
- The flagship window-merge fix (reading the ~174 kb span as one ~126 kb spliced gene) and the copy_graph overhang render — separate follow-ups.
- Re-running the consensus build, or any change to `promote_genomewide.py`, the discriminator, or the coding catalog.
- Recomputing protein hits (BLAST binaries are absent; protein stays a carried label).

## Global constraints

- **Purely additive:** no existing file is modified; no existing output is regenerated. The coding catalog and all prior outputs stay **byte-identical**.
- **Protein/ORF are labels, never gates.**
- **Every record carries `copy_vs_allele = "candidate-DNA-needed"` and `status = "flagged-reference-divergent-candidate"`, set unconditionally.**
- The gate is a **pure function** with an I/O shell (minimap2-free unit tests).
- Thresholds are named module-level constants (single-edit recall/precision knob).
- Deterministic; single minimap2 process; no `copy_assign`; outputs under `winloci_scratch`.
- `~REF`, `artifact`, `thin`, `not-own-locus` exclusions are **logged**, not silently dropped.
