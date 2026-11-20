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
