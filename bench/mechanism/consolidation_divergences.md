# Consolidation divergences (deliverable B)

Consolidations attempted whose byte-identity gate FAILED — kept both impls, deferred to
deliverable C. Empty means every attempted consolidation was byte-identical.

| B2 item | canonical target | what diverged (families/rows) | deferred to |
|---|---|---|---|
| B2.2 EM | em_assign_family | gstm.quant.tsv abundances differ (soft_quantify_em QUANT_ERROR=0.01/100-iter vs em_assign_family error_rate=0.003/eps=1e-6/200-iter — different error model AND convergence; entangled with the 5-epsilon inconsistency). Observed on the gate corpus: CAFAM0 copy_index 1 abundance 0.0564 (soft_quantify_em) vs 0.0565 (em_assign_family), all other rows unchanged — small but gate-detectable (md5 strict). | deliverable C |
