# DNA-derived PSV identifiability catalog (Phase 1)

- co-located classified pairs: **1387** -> resolvable **331** (24%), genuine-K=0 **1056** (76%). NOTE: this is the DNA reference universe (all aligned co-located pairs, including unexpressed identical tandem copies the RNA census never observes); on the **137** pairs the RNA census actually classifies, DNA and RNA agree **86%** and both find that expressed subset mostly resolvable.
- pairs excluded from K: **14262** unaligned (copy did not align to ref0 — divergent/short paralog), **14** unannotated (no overlapping GFF exon)
- **cross-check DNA-K=0 vs RNA-K0** on 137 census-classified pairs: concordance **86%** (confusion {(True, False): 14, (True, True): 8, (False, False): 110, (False, True): 5})
- discordant DNA-K=0 ∧ RNA-resolvable: **14** (candidate: indel / splice-shift pseudo-K=0 — substitution-only PSVs miss it; Phase-2 private_exon_bp will test this)
- discordant DNA-K≥1 ∧ RNA-tied: **5** (candidate: PSV in a poorly-expressed exon — reference identifiability ≥ read identifiability)

