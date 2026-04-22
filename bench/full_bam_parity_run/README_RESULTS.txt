Full-BAM parity run (GGO_19.bam)
================================

Inputs
  StringTie: GGO_19_stringtie.gtf  (precomputed: stringtie -L GGO_19.bam, v3.0.3)
  Rustle:    bench/full_bam_parity_run/GGO_19_rustle.gtf  (./target/release/rustle -L -p 8)

gffcompare: StringTie = reference, Rustle = query
  cmp_rs.*  (--stats / .tmap used by scripts/full_bam_rustle_vs_stringtie_audit.py)

Headline (cmp_rs.stats)
  Transcript level:  Sensitivity 92.8%  |  Precision 9.7%
  Matching transcripts: 1707 / 1839 reference mRNAs (stats file; tmap lists 1807 ref rows)
  Query side: 17540 transcripts in 652 loci (multi-exon-heavy)

Systemic miss breakdown (100 StringTie refs without an '=' match to any Rustle tx)
  per_ref_class among misses:
    c  51   (contained / partial match class — predcluster overlap, isofrac, fragmentation)
    j  28   (intron boundary / splice shift vs StringTie — junction snap, -E window, longtrim)
    k   9   (single-exon / local class — SE hull, singlethr, runoff)
    n   5   (novel splice / local)
    o   4   (overlap / exon boundary)
    m   4   (fragmentation / multi-hit)

Reverse comparison (Rustle = reference, StringTie = query) — cmp_st.stats
  Transcript level:  Sn 11.8%  |  Pr 92.8%
  (Expected: Rustle emits far more loci/transcripts than StringTie; StringTie is a small query set.)

Artifacts
  systemic_audit_rs_vs_st.tsv          — aggregated class counts
  systemic_audit_rs_vs_st.miss_refs.tsv — 100 non-exact StringTie transcript IDs + resolved class

Priority for StringTie parity (by volume)
  1) Precision / extra Rustle isoforms  — pairwise, isofrac, dedup, junction_path caps (not sensitivity)
  2) c-class                             — containment / subset vs dominant isoform
  3) j-class                             — splice-site agreement with StringTie (-E, RUSTLE_EMIT_REF_SSDIST)

Next steps
  - Long-read (`-L`) uses StringTie-like rules: higher `-f` floors, junction-path cap, per-splice
    pruning, containment drop, etc. (`RunConfig::apply_stringtie_parity_preset`). Re-run full BAM
    after changes; `--max-sensitivity` is the escape hatch for diagnostics.
  - Short-read: pass `--stringtie-parity` to apply the same preset without `-L`.
  - Trace miss refs: ./scripts/run_full_bam_stringtie_trace.sh GGO_19.bam GGO_19_stringtie.gtf trace_out
    or: rustle --trace-reference GGO_19_stringtie.gtf --trace-output trace.txt ...
  - Cross-check hotspots with rlink.cpp line refs in src/rustle/*.rs.
