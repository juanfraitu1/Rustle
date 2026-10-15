# PREREG — are `candidate_copies`' remaining false flags other apes' transcripts? (2026-10-14; to run on the main machine)

**Written before any of the alignments below.** It follows the F < 5 branch of `docs/PREREG_o3_fewread_flags_2026-10-14.md`:
a clean simulation (96 NULL-type units, paralogues down to 1%, present copies down to 4 reads, 5' G read ends) produced no
few-read false flag, so the A14 / H false flags come from something the simulation lacks. This file tests the candidate with
independent evidence on the same library (M3).

## Background, and what is and is not known

- The A14 (no-deletion control) and H (held-out, nothing deleted) runs use reads of the KB3781 fibroblast library. After the
  Amendment 15 fix they still carry false flags: A14 16 (2 b_allele + 14 c_unmatched, 10 families; C1' bar ≤ 8 failed), H 10
  (2 b + 8 c, 6 families; H1 bar ≤ 4 failed) (`docs/O3_CANDIDATES_AMENDMENT15_2026-10-05.md`).
- On the same library, `docs/PREREG_unmapped_rescue_2026-10-08.md` Amendment 43 found that most RNA-only wrong flags of the
  unmapped-rescue chain best match another ape (orangutan 279, siamang 244, chimpanzee 144, gorilla 29, human 2 of 698), with
  their variant k-mers absent from the animal's WGS. How such reads entered the library was NOT established (cross-sample
  contamination of multiplexed ape libraries is the stated candidate; no SMRT cell is shared with the project's own chimp /
  orangutan libraries). Those flags are not the `candidate_copies` flags; whether the A14 / H false flags are of the same kind
  has never been checked.

## Inputs (main machine; nothing here is re-run)

- Flagged unions: `iso.contigs.fa` (+ `iso_names.tsv`) of the corrected runs, `/mnt/linuxdisk/tmp/rna_allele/a14a15/` and
  `/mnt/linuxdisk/tmp/rna_allele/a15h/` (`accept_o3_candidates.sh contigs`), with each union's Amendment 9 class from the same
  runs' classify output (a / b / c / pri). The unit is a flagged union; **false** = class b, c or pri (16 + 10 = 26 expected;
  any other count is reported and used as found).
- Genomes, each aligned as the FASTA or the splice index already on disk:
  gorilla primary `GGO.fasta` (GCF_029281585.2) and KB3781's maternal / paternal assemblies (`mGorGor1.{mat,pat}.splice.mmi`)
  = **own**; human CHM13 v2.0 (`npip_ladder/idx/target.splice.mmi`), chimpanzee mPanTro3 (GCF_028858775.2), orangutan mPonPyg2
  (GCF_028885625.2), siamang mSymSyn1 v2.1 (GCF_028878055.3, the copy Amendment 43 used) = **foreign**.

## Command and score

`minimap2 -c -x splice:hq -uf -N 20 -t 4 GENOME iso.contigs.fa > iso.<genome>.paf` per genome (minimap2 2.30).
Score of a union on a genome = matching bases of its primary records / union length (capped at 1; Amendment 43c's score).
`own` = max over the three gorilla assemblies; `foreign` = max over the four others, with the species that gives it.

## Classes, fixed now (per false union)

- **FOREIGN**: foreign ≥ 0.995 and foreign − own ≥ delta/2 = 0.00479 (Addendum 1's foreign rule, 2026-09-23).
- **FOREIGN-LEANING**: foreign − own ≥ 0.00479 but foreign < 0.995.
- **OWN**: own ≥ foreign − 0.00479 (no other ape explains it better; the next suspects are real-read error structure and
  real paralogue complexity).

**Self-check, required:** the class-a flags (GWFAM175_B0 in A14 and the H class-a union) must come out OWN with own ≥ 0.999.
If either comes out FOREIGN the score is broken and nothing below is read.

## Reading, fixed now

| FOREIGN among the false unions (A14 + H pooled, n ≈ 26) | reading |
|---|---|
| ≥ half | **M3 is the main cause.** The stage's remaining false flags are other apes' transcripts in the library. A species-panel screen becomes the candidate fix, under its own pre-registration that counts its cost in real flags (Amendment 43c's screen kept only 80.5% of real flags); C1' / H1 recomputed without the FOREIGN unions are reported beside, never as a pass. |
| ≤ 5 | **M3 is not the cause.** The false flags are gorilla sequence; each is diagnosed individually (read origins against mat / pat, per-column cover of its consensus, the haplotype it best matches) before any rule is proposed. |
| in between | mixed; both follow-ups, per union. |

Reported beside: every false union with its length, reads, d, own and foreign scores and best species; FOREIGN-LEANING
counts; the same table for the class-a and true flags.

# OUTCOME

(to be appended after the run)
