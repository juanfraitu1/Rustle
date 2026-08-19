# Expanded O1 known-family purity audit

**Scope:** 19 post-hoc family-bearing graphs from frozen gorilla and human catalogs. This
measures within-emitted-graph typing purity and orientation sensitivity. It is not discovery
recall: choosing a catalog family after emission conditions on Rustle having found it.

Primary graphs retain **108/133** audit nodes and all **75**
independently named target-family members. They withhold **16** nodes named as
unrelated genes and report **1** broad-family/recent-subfamily outgroup separately.
The orientation guard removes **15** within-family edges. Named
members remain connected in **15/15 testable** cases; the remaining
**4** cases have only one independently named member and are reported as
annotation-limited rather than forced passes.

| species | family | audit→primary nodes | named | related outgroups | conflicts | R0→forward edges | primary components | status |
|---|---|---:|---:|---:|---:|---:|---:|---|
| GGO | C4 | 2→2 | 2 | 0 | 0 | 1→1 | 2 | PASS_NAMED_CONNECTED |
| GGO | TSPYL | 2→2 | 2 | 0 | 0 | 1→1 | 2 | PASS_NAMED_CONNECTED |
| GGO | GSTM | 4→4 | 4 | 0 | 0 | 5→5 | 4 | PASS_NAMED_CONNECTED |
| GGO | RGPD | 6→6 | 1 | 0 | 0 | 9→9 | 6 | ANNOTATION_LIMITED |
| GGO | ANKRD18 | 4→4 | 1 | 0 | 0 | 4→4 | 4 | ANNOTATION_LIMITED |
| GGO | HERC2 | 2→2 | 1 | 0 | 0 | 1→1 | 2 | ANNOTATION_LIMITED |
| GGO | MAGEA | 11→10 | 5 | 0 | 1 | 47→47 | 10 | PASS_NAMED_CONNECTED |
| GGO | DAZ | 2→2 | 2 | 0 | 0 | 1→1 | 2 | PASS_NAMED_CONNECTED |
| GGO | APOBEC3 | 2→2 | 2 | 0 | 0 | 1→1 | 2 | PASS_NAMED_CONNECTED |
| HSA | GOLGA6_8 | 22→19 | 11 | 1 | 1 | 86→86 | 19 | PASS_NAMED_CONNECTED |
| HSA | TBC1D3 | 6→6 | 5 | 0 | 0 | 15→15 | 6 | PASS_NAMED_CONNECTED |
| HSA | RANBP2_RGPD | 9→8 | 8 | 0 | 1 | 31→31 | 8 | PASS_NAMED_CONNECTED |
| HSA | RABL2 | 2→2 | 2 | 0 | 0 | 1→1 | 2 | PASS_NAMED_CONNECTED |
| HSA | AMY1 | 6→6 | 2 | 0 | 0 | 15→15 | 6 | PASS_NAMED_CONNECTED |
| HSA | PCDHB | 2→2 | 2 | 0 | 0 | 1→1 | 2 | PASS_NAMED_CONNECTED |
| HSA | MAGEA | 10→7 | 7 | 0 | 3 | 25→25 | 7 | PASS_NAMED_CONNECTED |
| HSA | RBMY1 | 3→3 | 3 | 0 | 0 | 3→3 | 3 | PASS_NAMED_CONNECTED |
| HSA | NBPF_core | 21→20 | 14 | 0 | 1 | 87→87 | 20 | PASS_NAMED_CONNECTED |
| HSA | NBPF_repeat_bridge | 17→1 | 1 | 0 | 9 | 36→21 | 1 | ANNOTATION_LIMITED |

`graphs/<species>.<family>.gene_family.gfa` is the purified view. Its sibling `.audit.gfa`
retains every emitted node and reverse-only edge. Purple nodes are documented broad-family
homologs outside the recent-copy subfamily; orange nodes are independently named unrelated
genes; grey nodes are untyped; blue nodes are unnamed candidates connected to a named member.

The panel deliberately contains difficult mixtures (GOLGA6/8, RANBP2/RGPD, MAGEA, NBPF) and
the known NBPF/TTC6-DNAH14 repeat-bridge component. These are more informative for purity than
adding only clean two-copy pairs. Family stems are an evaluation typing device, not a production
coordinate blacklist or a substitute for O1 discovery.
