# KEY=ourhomology
edges: Soto's exon map-back 12,231; ours 10,731 gene pairs (86 loci folding several annotation records)

| homology | copy numbers | ARI all / dev / held-out | exact all (dev / held-out) | nesting | bipartite sens / prec |
|---|---|---|---|---|---|
| H0 Soto map-back | none (sequence only) | 0.7307 / 0.6418 / 0.8693 | 345 | 440/444 | - |
| H0 Soto map-back | S1C (Soto's) | 0.9698 / 0.9708 / 0.9681 | 479 (216 / 263) | 440/444 | 0.980 / 1.000 |
| H0 Soto map-back | ours | 0.9277 / 0.9227 / 0.9343 | 411 (177 / 234) | 440/444 | 0.935 / 0.975 |
| H1 ours | none (sequence only) | 0.7648 / 0.7146 / 0.8242 | 186 | 257/444 | - |
| H1 ours | S1C (Soto's) | 0.8235 / 0.7675 / 0.8844 | 249 (113 / 136) | 257/444 | 0.690 / 0.966 |
| H1 ours | ours | 0.7856 / 0.7323 / 0.8446 | 221 (91 / 130) | 257/444 | 0.656 / 0.947 |

VERDICT (H1 with S1C famCN, held-out ARI 0.8844 vs H0 0.9681): PARTIAL
