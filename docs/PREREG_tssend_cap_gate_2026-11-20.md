# Row 1834: does the assembler's cap signal say where the median 5' end helps? (pre-registered 2026-11-20)

User (2026-11-20): test the median 5' end layer gated on a 5' signal. Row 1831: the layer of `final.gtf`
(`--polish-tss-end median`, row 1487, set on A119b) moves the 5' end from the most upstream read start to the median read
start; it helps A119b (+25 points of correct-chain transcripts with the 5' end within 50 bp of the annotated start) and mouse
ES (+17), and is worse on five datasets (OR6737 -23, the others -1 to -2).

**The signal.** The binary already tests each contig for a cap signature (non-templated G at read 5' ends: capped share at the
proven first-exon clusters against exon bodies; `polish_tss_cap_signal` in `PREFIX.params.tsv`) and gates the TSS rescue on it.
No annotation is used. A library "has a cap signal" when more than half of its tested contigs say yes.

**Development table (row 1831 runs, looked at before this registration).**
| dataset | contigs yes / no | capped at first-exon clusters | layer effect, points |
|---|---|---|---|
| fly sim | 0 / 36 | 0.006 | -1.5 |
| mouse sim | 0 / 26 | 0.008 | -1.2 |
| zebrafish sim | 0 / 143 | 0.006 | -1.5 |
| mouse ES | 0 / 110 | 0.038 | +17.4 |
| mouse heart | 0 / 93 | 0.094 | -2.0 |
| OR6737 | 0 / 26 | 0.040 | -23.4 |
| A119b | 24 / 1 | 0.693 | +25.1 |
"Cap signal = layer helps" holds on 6 of 7; mouse ES is the exception (no signal, layer helps). A gate would keep A119b's gain,
avoid OR6737's loss and forgo mouse ES's gain.

**Test, on six libraries not looked at for this layer:** mouse cortex (ENCSR862DMK, GENCODE vM27), KB3781 gorilla fibroblast,
chimpanzee testis (`PTR_mm`), bonobo testis (`PPA_mm`), chimpanzee LCL (`PTR_AG18354_lcl_mm`), bonobo fibroblast
(`PPA_PR00251_fib_mm`); RefSeq per species. Arms `base` and `tssend_off` (`RUSTLE_POLISH_TSS_END=off`), stage `assemble`,
`RUSTLE_PTC=0`; scorer of row 1831. Workspace `r1834_tssend_cap_gate`. Registered before any arm ran.
- Layer effect = share of correct-chain transcripts of `final.gtf` with the 5' end within 50 bp of the nearest annotated start of
  that chain, base minus `tssend_off`. *Helps* = at least +1 point; *hurts* = at most -1 point.
- **Prediction per library, fixed by the rule:** cap signal -> helps; no cap signal -> does not help.
- **Bar G1:** the prediction is right on at least 5 of 6. **G2:** no library with a cap signal is hurt.
- Reported: the gated layer's effect per library (the layer's where there is a signal, 0 elsewhere) against always on and
  always off, over all thirteen libraries.
**Decision.** G1 and G2 = the gate is a candidate (the layer applied only on contigs with a cap signal, as the rescue already
is); the user rules. Otherwise reported; the layer's default stays the user's call with the thirteen-library table.

**Limits, stated in advance.** Annotated starts are weak truth (rows 1718 / 1742 used CAGE on human). The mouse ES exception
means the cap signal cannot be the whole explanation: the layer helps when the most upstream read overshoots the annotated
start, which a cap test only partly captures. Five of the six test libraries are apes. A library can pass G1 trivially when the
layer does nothing on it.

## Result (2026-11-20): FAIL; the cap signal does not say where the layer helps, and the layer hurts five of six test libraries

Jobs 4484231-4484236, every arm exit 0; `r1834_tssend_cap_gate/score.NAME.txt`. Repo `9956b32f`, binaries of 2026-11-20.

| library | cap signal (contigs yes / no; capped share at clusters) | 5' within 50 bp, with -> without the layer | effect, points | predicted | right |
|---|---|---|---|---|---|
| mouse cortex | no (0 / 93; 0.079) | 0.626 <- 0.567 | +5.8 | no help | no |
| KB3781 gorilla fibroblast | no (0 / 26; 0.651) | 0.367 <- 0.615 | -24.7 | no help | yes |
| chimpanzee testis | yes (25 / 1; 0.556) | 0.427 <- 0.578 | -15.1 | helps | no |
| bonobo testis | yes (25 / 1; 0.547) | 0.312 <- 0.463 | -15.1 | helps | no |
| chimpanzee LCL | no (0 / 26; 0.748) | 0.423 <- 0.651 | -22.8 | no help | yes |
| bonobo fibroblast | no (0 / 26; 0.768) | 0.287 <- 0.483 | -19.6 | no help | yes |

**G1 FAIL** (3 of 6 right). **G2 FAIL** (both libraries with a cap signal are hurt, by 15 points). No gate proposed.

Over the thirteen libraries (rows 1831 + 1834), effect of the layer in points: A119b +25.1, mouse ES +17.4, mouse cortex +5.8;
mouse heart -2.0, the three simulated sets -1.2 to -1.5; OR6737 -23.4, KB3781 -24.7, chimpanzee testis -15.1, bonobo testis
-15.1, chimpanzee LCL -22.8, bonobo fibroblast -19.6. Mean: always on -6.0, gated on the cap signal -0.4, always off 0.
With the layer `final.gtf` also holds 9-21 fewer correct chains on the six test libraries and 0.5-1.3 points less base
sensitivity (base precision +0.1 to +1.4).

**Reading (my interpretation).** The layer was set on A119b and is right there and on two of three real mouse libraries; on
every other ape library it moves starts away from the annotated ones by 15-25 points. Two explanations cannot be separated
here: (a) the layer fits human / mouse libraries and not the others; (b) the annotations differ: human RefSeq and GENCODE starts
are curated, while RefSeq models of the other apes are largely predicted from RNA evidence and may sit at the most upstream
read, which is what the assembler emits without the layer. Row 1487 judged the layer on A119b with RefSeq and row 1718 with
CAGE; no CAGE exists for the other apes. Either way the layer is not general as measured, and no reference-free signal tested
here predicts where it helps (the cap signal also reads "no" on libraries where 65-77% of the reads at first-exon clusters are
capped: its verdict depends on the exon-body rate too).
Candidate for the user: the median 5' end as an opt-in layer (`RUSTLE_POLISH_TSS_END=median` already exists; the default of
`final.gtf` would change), stated with explanation (b) as the open question.
