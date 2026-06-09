## Task 4 integration (NC_073224.2 rna-XM_063708549.1)

OFF has chain: True
ON  has chain: True

union-baseline log: (no [VG-UNION-BASELINE] log lines emitted)

### Interpretation: DONE_WITH_CONCERNS

Both OFF and ON outputs are **byte-identical** (0 diff lines, 334 lines each).
The flag `RUSTLE_VG_UNION_BASELINE=1` did not trigger on this locus slice — no
`[VG-UNION-BASELINE]` log lines appeared in `/tmp/un.log`.

The slice (NC_073224.2:39996000-40059000) already emits a transcript matching the
reference intron chain in both modes (`RSTL.1.5`, 16 exons, chain_len=15). The
over-collapse regression that drops `rna-XM_063708549.1` in the genome-wide run
requires the full secondary/family cross-mapping context from surrounding
paralog copies, which is absent in this narrow BAM slice.

The fix was not exercised on this slice. The real validation must come from
the genome-wide run (next task, T12/T13) where the full over-collapse context
is present.

### Run details

- Control (OFF): exit 0, RSTL.1.5 has chain
- Treatment (ON):  exit 0, RSTL.1.5 has chain, outputs byte-identical to OFF
- Flags used: RUSTLE_VG_TANDEM=1 RUSTLE_VG_TANDEM_PRIMARY_JUNCTIONS=1
              RUSTLE_VG_DECISIVE_GATE=1 RUSTLE_VG_DECISIVE_GATE_MIN_PRIM=4
