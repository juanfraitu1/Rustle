# PREREG END ARTIFACTS, PAS EVIDENCE (register 1502, 2026-11-02, before code or numbers)

Why: r1501 rule J lost 0.83% of `=` (bar 0.5%); 184 of the 325 lost were A-run-only drops. Reads carry no poly(A) tails (FLNC
trimmed: 1/1,639 A119b reads, 0 in OR6737 / GGO / HG002). Genome evidence instead: a polyadenylation signal (PAS) upstream of a real
cleavage site.
PAS = any of AATAAA ATTAAA AGTAAA TATAAA CATAAA GATAAA AATATA AATACA AATAGA ACTAAA AATGAA (Tian 2005 hexamers, AAAAAG left out:
an A-run matches it) wholly inside the 40 bp before the 3' end, transcript orientation.
Rule K = [RTS part of rule J AND reads <= 2] OR [A share > 60% AND no PAS AND no support model 3' end within 10 bp] (no read cap on
the A part). Constants otherwise as r1501. Input, reference and bars as r1501: chain precision >= +0.5, `=` loss <= 0.5%,
multicopy/single <= 1.5x. Pass -> RUSTLE_FINAL_END_ARTIFACTS layer proposed (separate commit). Report: PAS share of A-run models.
