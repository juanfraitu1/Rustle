# PREREG: which tied copies can sequence differences decide, and does StringTie read secondaries (row 1870, 2026-11-21)

User (2026-11-21): the claim must be airtight for the advisor; StringTie in long-read mode is understood to use primaries
only; "is there anything more clever that would allow to disambiguate which ones really belong via PSVs?" Written before
either measurement. Earlier PSV rows (1778, 1783, 1829, 1844) asked where a READ belongs and found no gain over the primary;
this asks whether a COPY has evidence, from all reads of its family at once.

## A. StringTie and secondary alignments
`stringtie -L` (the command of `apetools/job.sbatch`) on chimp testis `PTR_mm.bam` as it is (primaries + secondaries, the
input of every StringTie assembly compared so far) and on the same BAM without secondaries (`samtools view -F 0x100`).
- IGNORES SECONDARIES: the two GTFs have the same set of transcripts (chrom, strand, exon coordinates).
- USES THEM: otherwise; then StringTie's counts in rows 1859-1869 already include secondaries, and the sentence "tools that
  discard secondaries" is wrong for StringTie. Reported: transcripts only in one, and named copies modeled by each.

## B. Copy-level evidence from positions where copies differ (`r1870_copy_psv/psv.py`)
Sample A119b, the mode's assembly (`r1865_multicopy_oneread/a1/rescue.gtf` statuses as in `r1863_copy_evidence/mode.a1.tsv`),
named families, copies with reads. A read aligned to several copies of its family (records with `=` / `X` CIGARs) shows, at
each of its bases, which copies it matches. Per copy C:
- private support at a site: the read matches C there and mismatches EVERY other family copy it is aligned to at that base.
- contradiction at a site: the read mismatches C there and matches another copy.
- two copies are IDENTICAL OVER THE READS when no read aligned to both differs between them at any base (>= 2 shared reads).
Verdict per copy: SUPPORTED (>= 1 site with >= 2 private-support reads); else IDENTICAL (to at least one other copy);
else REFUTED (>= 1 contradiction site with >= 2 reads); else UNRESOLVED. A sequencing error mismatches every copy and
counts for none.

## How it is read
- Calibration: the method is trusted only if >= 90% of the copies with own reads (row 1863's OWN) come out SUPPORTED.
- Then the tie-only copies are reported by verdict. IDENTICAL = no sequence-based method can tell the copies apart from
  these reads (irreducible; report the group). REFUTED = the copy can be dropped. SUPPORTED = it can be kept as a copy.
- OR6737 beside, same script.

## Predictions (mine)
A: StringTie ignores secondaries. B: calibration holds; most tie-only TSPY copies are IDENTICAL; a minority of tie-only
copies is REFUTED; few are SUPPORTED.

## Limits
Reference copies of CHM13 against another individual's reads: a private allele of the sample's copy looks like an error.
Sites are seen only through reads (no multiple alignment of the copies). Annotated copies of fifteen families.

## Result (2026-11-21)
Job 4491428; `r1870_copy_psv/{A.txt,st.*.sig,a1.tsv,or.tsv}`. StringTie 3.0.3.

### A. StringTie USES secondaries (by the registered rule), without effect on the named families
Chimp testis: 1,509,097 primary and 2,443,626 secondary records. `stringtie -L` on the full BAM: 62,418 transcripts,
identical to the existing `st_dn.gtf` (the run is deterministic). Without secondaries: 62,494; 61,268 shared, 1,226 only
without, 1,150 only with (2%). So StringTie does not ignore secondary records, and my prediction was wrong. On the named
families the two runs are the same: 121 copies / 0.73 / 424 transcripts, TSPY 6, NPIP 35. What StringTie does with a
secondary record was not established here (source not read).
Consequence for the claim: every StringTie assembly compared in rows 1855-1869 HAD the secondaries in its input. "Tools
that discard secondary alignments" is not a correct description of StringTie; "StringTie, given the same alignments, does
not turn them into copies" is.

### B. Copy-level evidence at differing positions, A119b
Calibration holds: 81 of 85 copies with own reads are SUPPORTED (95%; the others: TSPY2 and one TSPY3 REFUTED, SPATA31E1
and VCY1B UNRESOLVED). OR6737: 29 of 36 (81%; 6 REFUTED), no bar there.
The 38 tie-only copies: SUPPORTED 7 (BPY2B, CT47A5, CT47A6, HSFY1, HSFY2, two XAGE1A), IDENTICAL 14 (all TSPY),
REFUTED 17 (12 TSPY, CT47A1-4, XAGE1B). Predictions: most tie-only TSPY IDENTICAL held (14 of 26); "a minority REFUTED"
was wrong (17 of 38, the largest class); "few SUPPORTED" 7.
A119b copies of the fifteen families, mode's assembly: 88 SUPPORTED, 14 IDENTICAL (TSPY), 19 REFUTED, 2 UNRESOLVED.
Caution on REFUTED: it means no read carries a base only this copy explains while reads contradict it; a copy whose every
difference is shared with some other copy, or that differs only by insertions / deletions (not scored), is expressed and
still lands there. TSPY2 (12 own reads by score, REFUTED) is such a case. Read it as "not needed to explain the reads".
