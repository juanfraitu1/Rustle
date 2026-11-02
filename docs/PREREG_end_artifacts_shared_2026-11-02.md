# PREREG END ARTIFACTS, SHARED-EVIDENCE EXEMPTION (register 1501, 2026-11-02, before code or numbers)

Why: r1500 rule E (RTS or A-run, reads <= 2) lost 2.3% of `=`; the oracle arm that spared reference matches lost 0.008%. The
exemption must come from our own evidence. Rule J, on final.gtf, locus = (gene_id, contig, strand), support model = a model of the
same locus with reads >= 3:
- RTS part: some RTS junction of the model is used by no support model;
- A part: A share > 60% AND no support model has its 3' end within 10 bp (the --polish-subchain tolerance) of the model's 3' end;
- flag = (RTS part OR A part) AND reads <= 2. No tuning; constants as r1500 (SQANTI3 6.1 definitions, `utilities end-artifacts`).
Input A119b gw.LS final.gtf; gffcompare vs chains/ref.gtf. Bars: chain precision >= +0.5 points, `=` loss <= 0.5%, multicopy
(SD90) / single-copy removal <= 1.5x. Pass -> RUSTLE_FINAL_END_ARTIFACTS layer proposed as the default of PREFIX.final.gtf (a
separate commit, after the check). Fail -> recorded, tool stays opt-in. RTS-only and A-only parts reported, not decided on.
