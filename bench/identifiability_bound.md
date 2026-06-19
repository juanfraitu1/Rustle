# The identifiability bound on copy detection (3 independent confirmations)

The method recovers every copy the read data can DISTINGUISH; the rest are unidentifiable by definition.
Confirmed on three independent fronts this session:

1. DIVERGENT PARALOGS (Phase 1, divergent_phase1_features.py): on a gold Compara set (781 paralog /
   1294 non-paralog), no RNA sequence feature separates true divergent paralogs from domain-sharers
   (best AUC 0.629, protein; DNA <=0.56). Distinguishable only with phylogenetic-level input.

2. COLLAPSED COPIES (--recover-copies, recover_collapsed_copies): feeding AS-tied secondary reads +
   running copy_split PAST the family gate recovers 0 NEW copies in no-family collapsed arrays
   (DAZ1/RABL2A/GOLGA6L7/PRAMEF). Distinguishable copies form families (in-family split already gets
   them); collapsed copies are not distinguishable (that is why they collapsed). Exhaustive + complementary.

3. PRIOR VG work: "decisive secondary signal self-contradictory"; multimapper advantage confined to
   tied co-located copies which are by definition non-separable.

THESIS STATEMENT: copy detection/assignment is COMPLETE up to identifiability; identifiability is
information-theoretic (does a copy carry a copy-specific PSV/junction?), not algorithmic. Engineering
(secondary reads, longer k-mers, masking, protein-ORF) cannot cross it. The in-family PSV split realizes
the bound. Transcripts themselves are validated real (bench/transcript_validation.md: 86% intron precision,
98.9% gene-associated), so the bound is about copy SEPARABILITY, not transcript quality.
