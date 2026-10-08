# What happens to the reads of a copy the reference haplotype lacks, and what a reference-free step recovers (2026-10-08)

Study of `docs/PREREG_o3_maternal_reference_2026-10-08.md` (body, Amendment 1, Amendment 2); plan `docs/superpowers/plans/2026-10-08-o3-maternal-reference.md`; code `bench/o3_maternal/` (83 unit tests);
work dir `/mnt/linuxdisk/tmp/o3_mat/`; artifact https://claude.ai/artifact/R1QoJa9GLfxnt51haKG6k8 (private). Substrate: KB3781 (male) fibroblast Iso-Seq against its own maternal (`mat`) and paternal (`pat`) assemblies.
Every number below is copied from `W/<ref>/fate/fate.tsv`, `side.json`, `W/<ref>/score/*.json` and `W/compare_*.json`.

## 1. Results

**Q1, where the reads of a missing copy go.** In neither direction do they become unmapped. With the mother's genome as reference, the 7 judged large copies (1,132 reads) are 92.2% absorbed on a relative, 7.7% tied and 0.1% unmapped or partial; with the father's genome as reference, the 13 judged copies (1,217 reads) are 97.4% absorbed, 2.0% tied and 0.7% unmapped or partial. The control of all other reads is 6.4% tied (mother) and 5.0% tied (father), so the tying seen at 2% to 8% of the missing copies' reads is of the same size as the ordinary multi-copy background (the exceptions are the loci tied almost entirely). The registered claim ("missing copies hide as unmapped or multi-mapped reads") is REFUTED at 5 of 7 judged copies with the mother as reference and 11 of 13 with the father as reference; the exceptions are GWFAM205_B0 (81 of 82 reads tied with a paralog 0.65% away) and GWFAM390_B1 (20 of 20 tied, 0.58% away), plus GWFAM64_B11 (between bars) and GWFAM382_B0 (partly).
The only judged copy farther than the allele cutoff from any relative on the reference is GWFAM175_B0 (304 reads): its reads are absorbed on a copy 8.7% away, at median divergence 0.0662 on the mother's genome against 0.0007 on the father's, with MAPQ 60 on both. Three other judged copies have an identical or near-identical counterpart on the mother's genome (GWFAM208_B0 identity 1.0000, GWFAM64_B17 0.9999, GWFAM214_B0 0.9995), so they are absent by the registered synteny rule and present by sequence.
Reads that could not be mapped at all do exist, but not for this reason: of the 959 reads unmapped on the combined primary assembly, 132 map on the mother's genome and 3 on the father's, so those are sequences the primary assembly lacks.

**LRPAP1.** By the registered rule the only LRPAP1 locus absent from the mother's genome is the chrY copy (sex control, 1 read). p12 (the chr12 5′ fragment) is not absent: it lifts as class `T?` to a mat locus 97.2% identical (Amendment 2); the chr12 cluster has 4 loci on `pat` and 3 on `mat`, and the registered lifts are many-to-one there, so which locus is extra is unresolved. p12 is kept as a descriptive locus: all 83 of its reads are TIED with p14 even on `pat` (the two differ by 4 bp in 2,446), and on `mat` 66 are absorbed on the p14 ortholog (99.59%), 16 tied, median divergence 0.0025 against 0.0008 on `pat`.

**Q2, recovery from the maternal alignment alone.** IsoCon recovers GWFAM175_B0, the only copy beyond the allele cutoff that has at least 20 reads: R1 PASS, with 1 of 16 expressed copies recovered overall (9 flagged candidates, 1 matching a truth locus). R2 passes (13 within-cutoff copies, none flagged) but R3 FAILS: 6 of 24 families without a missing copy (25%, bar 20%) carry a flagged candidate that recovers nothing. The in-house stage (with the Amendment 15/15b fix) flags 1 candidate and recovers nothing: R1 FAIL, R2 PASS, R3 PASS (1 of 24). With the father's genome as reference (held-out direction), R1 cannot be tested (no judged copy beyond the cutoff), R2 passes for both arms (26 copies), R3 passes (IsoCon 2 of 21, in-house 0 of 21), and neither arm recovers any of the 30 expressed copies.

**Side by side.** On the genome that has the copy, the reads sit on it at median divergence within the allele cutoff (≤ 0.00958) for 20 of the 21 large copies; the exception is GWFAM382_B0, whose reads are 5–6% divergent from BOTH genomes (see section 6, item 7). The chain's view of GWFAM175_B0: in the mother's-genome run IsoCon has 1 transcript of the copy, kept as a new copy and recovered by a candidate with 4 transcripts; in the father's-genome run the same transcript matches the reference and is not flagged.

## 2. Fate tables (LARGE = at least 20 reads)

### Reference = mother's genome (`mat`); truth = `pat`

| locus | kind | reads | unmapped + partial | tied | on nearest relative | elsewhere | nearest-relative identity | de (reference) | de (other haplotype) | verdict | registered prediction (REFUTED) |
|---|---|---|---|---|---|---|---|---|---|---|---|
| GWFAM208_B0 | catalog | 551 | 0 | 1 | 550 | 0 | 1.0000 | 0.0008 | 0.0008 | REFUTED | held |
| GWFAM175_B0 | catalog | 304 | 0 | 2 | 17 | 285 | 0.9126 | 0.0662 | 0.0007 | REFUTED | held |
| GWFAM214_B0 | catalog | 122 | 0 | 0 | 122 | 0 | 0.9995 | 0.0010 | 0.0010 | REFUTED | held |
| LRPAP1_p12 | lrpap1_desc | 83 | 1 | 16 | 66 | 0 | 0.9959 | 0.0025 | 0.0008 | descriptive |  |
| GWFAM205_B0 | catalog | 82 | 1 | 81 | 0 | 0 | 0.9935 | – | 0.0006 | SUPPORTED | FAILED |
| GWFAM64_B17 | catalog | 31 | 0 | 0 | 31 | 0 | 0.9999 | 0.0015 | 0.0015 | REFUTED | held |
| GWFAM227_B0 | catalog | 21 | 0 | 0 | 21 | 0 | 0.9911 | 0.0057 | 0.0013 | REFUTED | held |
| GWFAM64_B11 | catalog | 21 | 0 | 3 | 18 | 0 | 0.9980 | 0.0008 | 0.0024 | between bars | FAILED |

Control, all other reads (29,807): 193 partial, 1,910 tied, 27,704 absorbed, median de 0.0016.

### Reference = father's genome (`pat`); truth = `mat`

| locus | kind | reads | unmapped + partial | tied | on nearest relative | elsewhere | nearest-relative identity | de (reference) | de (other haplotype) | verdict | registered prediction (REFUTED) |
|---|---|---|---|---|---|---|---|---|---|---|---|
| GWFAM309_B1 | catalog | 482 | 2 | 0 | 480 | 0 | 0.9991 | 0.0010 | 0.0011 | REFUTED | held |
| GWFAM70_B2 | catalog | 195 | 1 | 0 | 194 | 0 | 0.9999 | 0.0007 | 0.0007 | REFUTED | held |
| GWFAM402_B0 | catalog | 100 | 2 | 1 | 97 | 0 | 0.9991 | 0.0018 | 0.0019 | REFUTED | held |
| GWFAM70_B4 | catalog | 97 | 0 | 0 | 97 | 0 | 0.9927 | 0.0013 | 0.0090 | REFUTED | held |
| GWFAM70_B3 | catalog | 65 | 0 | 0 | 65 | 0 | 0.9998 | 0.0012 | 0.0012 | REFUTED | held |
| GWFAM390_B0 | catalog | 59 | 0 | 1 | 58 | 0 | 0.9996 | 0.0010 | 0.0010 | REFUTED | held |
| GWFAM169_B1 | catalog | 46 | 0 | 0 | 46 | 0 | 1.0000 | 0.0010 | 0.0010 | REFUTED | held |
| GWFAM162_B0 | catalog | 41 | 0 | 0 | 41 | 0 | 0.9993 | 0.0008 | 0.0011 | REFUTED | held |
| GWFAM309_B4 | catalog | 35 | 0 | 0 | 35 | 0 | 1.0000 | 0.0009 | 0.0009 | REFUTED | held |
| GWFAM70_B1 | catalog | 29 | 0 | 0 | 29 | 0 | 0.9967 | 0.0043 | 0.0026 | REFUTED | held |
| GWFAM169_B0 | catalog | 26 | 0 | 0 | 1 | 25 | 0.9992 | 0.0008 | 0.0012 | REFUTED | held |
| GWFAM382_B0 | catalog | 22 | 3 | 2 | 17 | 0 | 0.9948 | 0.0565 | 0.0570 | PARTLY | FAILED |
| GWFAM390_B1 | catalog | 20 | 0 | 20 | 0 | 0 | 0.9942 | – | 0.0013 | SUPPORTED | FAILED |

Control, all other reads (29,148): 189 partial, 1,461 tied, 27,498 absorbed, median de 0.0017.

The "absorbed" label in the control row has no nearest-relative split (control reads have no registered nearest relative), so "elsewhere" there means "mapped untied on a primary".

## 3. Selection caveat

Every read of the 34-family and LRPAP1 sets was mapped on the combined primary assembly before this test, so UNMAPPED on `mat`/`pat` here means "mapped on `_pri`, lost on this genome" and is 0 by construction in the fate table; reads unmapped on `_pri` are the separate set of 959 reads (132 map on `mat`, 3 on `pat`, query coverage ≥ 0.8). No pooled percentage across loci is a rate: the judged loci are 7 (mat) and 13 (pat), from one animal and one tissue.

## 4. Amendment 1 predictions

- **S1** (median de lower on the haplotype that has the copy, within the allele range): holds for 4 of 8 large loci with the mother as reference and 2 of 13 with the father as reference. The failures are loci whose counterpart on the reference is identical or within 0.1% (equal medians) and GWFAM64_B11, whose paternal median is higher; GWFAM382_B0 is divergent on both. With the mother as reference S1 holds for exactly the four copies whose nearest relative there is not a near-identical copy (GWFAM175_B0, GWFAM205_B0, GWFAM227_B0, p12); its second clause (≤ 0.00958 on the present haplotype) holds for 20 of 21 large copies, the exception being GWFAM382_B0.
- **S2** (the transcripts exist in both runs; kept as new and recovered only in the absent run): holds for GWFAM175_B0 (1 output on the locus in each run; kept new and recovered only in the `mat` run). For the other loci IsoCon has outputs on the locus in both runs (for example GWFAM208_B0 9 and 5) and none is kept as new, because they link to a reference locus within the cutoff.
- **S3** (the `pat` run applies the registered bars unchanged): reported above; R1 not testable, R2 and R3 pass, nothing recovered.

## 5. Arms

### Reference = `mat`

| arm | R1 | R2 | R3 (bar ≤ 20%) | candidates / flagged | expressed copies recovered | flagged candidates by class |
|---|---|---|---|---|---|---|
| IsoCon | PASS (1 of 1) | PASS (0 of 13 flagged) | FAIL (6 of 24) | 21 / 9 | 1 of 16 | {'c_unmatched': 2, 'b_other': 6, 'a_recovered': 1} |
| in-house stage | FAIL (0 of 1) | PASS (0 of 13 flagged) | PASS (1 of 24) | 20 / 1 | 0 of 16 | {'c_unmatched': 1} |

p12 line (IsoCon outputs with identity × coverage ≥ 0.999): on p12 (paternal) 1, on p14 (paternal) 1, on p14 (maternal) 1, on a p12 counterpart in `mat` 0. Three sequence types are present in the reads at the p12/p14 locus; the p12-type output links to p14 by the registered allele cutoff and is not kept as a new copy.

### Reference = `pat` (held-out direction)

| arm | R1 | R2 | R3 (bar ≤ 20%) | candidates / flagged | expressed copies recovered | flagged candidates by class |
|---|---|---|---|---|---|---|
| IsoCon | NOT TESTABLE (0 of 0) | PASS (0 of 26 flagged) | PASS (2 of 21) | 22 / 7 | 0 of 30 | {'b_other': 4, 'c_unmatched': 3} |
| in-house stage | NOT TESTABLE (0 of 0) | PASS (0 of 26 flagged) | PASS (0 of 21) | 16 / 1 | 0 of 30 | {'c_unmatched': 1} |

## 6. What did not go as predicted

1. S5 prediction "REFUTED at every LARGE locus" fails at GWFAM205_B0 (SUPPORTED) and GWFAM64_B11 (between bars) on `mat`, and at GWFAM390_B1 (SUPPORTED) and GWFAM382_B0 (PARTLY) on `pat`.
2. p12 is not a mother-absent locus (Amendment 2); the LRPAP1 autosomal difference is "4 loci vs 3 on chr12, which one is extra is unresolved", and the only absent locus is the sex-chromosome copy.
3. Several catalog "absent" loci have identical counterparts on the reference (208_B0, 214_B0, 64_B17): the synteny-based truth rule overstates absence.
4. S1 fails for 4 of 8 (mat) and 11 of 13 (pat) loci.
5. IsoCon R3 FAILS on `mat` (25%); the in-house stage fails R1.
6. The in-house arm is not the registered one: the binary contains the Amendment 15/15b fix (commit a13b817f), contrary to the prereg's "not applied" (Amendment 2); its acceptance re-run was not run.
7. Observation, not a result: GWFAM382_B0's 22 reads are 5–6% divergent from both haplotype assemblies (de 0.0526 on `pat`, 0.0570 on `mat`), so they come from a copy present in neither genome, the class of the 2026-08-13 reference-absent candidates; this study did not test it.

## 7. Provenance

minimap2 2.30-r1287 (baseline BAM 2.31); reads aligned with `-ax splice:hq -uf --eqx -Y -N 50 -p 0.1 --secondary=yes -K 100M`; nearest-relative identities asm20 `-N 50 -p 0.1`; IsoCon 0.3.3 (`IsoCon pipeline --nr_cores 4`, defaults); `o3_candidates` binary built 2026-10-06 09:09 (sha256 prefix 1e01f514e387cadf), after commit a13b817f; delta 0.00958; the `mat` verdicts were written to `W/mat/VERDICTS.txt` at 2026-10-08 08:50:06 and the first `pat` file at 08:50:20, the `pat` run applying the same code and constants. Label cross-check against Amendment 10 (`express`, both haplotypes jointly): same order for GWFAM175_B0 304/281, GWFAM205_B0 82/77, GWFAM227_B0 21/13, but GWFAM208_B0 551/0, GWFAM214_B0 122/0, GWFAM64_B17 31/0, GWFAM64_B11 21/0 (those reads tie across haplotypes in Amendment 10 and are untied on `pat` here). Execution rulings are listed in the session ledger and in Amendment 2.

## 8. Not covered

Between-individual differences (one animal); the testis library; copies absent from both haplotypes (GWFAM382_B0 above); the two `_pri` copies with a lift below 50% (GWFAM175:2, GWFAM491:1); which chr12 LRPAP1 locus is the extra one; reads unmapped on `_pri` beyond the 959; any threshold tuning (none); the Amendment 15 acceptance re-run of the in-house stage; independent replication.
