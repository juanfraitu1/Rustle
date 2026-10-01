# KEY=unitcovercn
multi genes 149, clean 2071; location baseline mean Jaccard dev 0.221 / held-out 0.277

| arm | famCN | units | Jaccard dev / held-out | p dev / held-out | exact | recall | precision | clean: own family / 2+ families / empty |
|---|---|---|---|---|---|---|---|---|
| A | - | 1546 | 0.439 / 0.473 | 0.0010 / 0.0010 | 0.107 | 0.555 | 0.718 | 0.546 / 0.303 / 0.020 |
| B | - | 1241 | 0.440 / 0.464 | 0.0010 / 0.0010 | 0.121 | 0.496 | 0.800 | 0.607 / 0.201 / 0.023 |
| C | s1c | 2255 | 0.548 / 0.624 | 0.0010 / 0.0010 | 0.215 | 0.917 | 0.617 | 0.644 / 0.224 / 0.045 |
| C | ours | 2400 | 0.492 / 0.601 | 0.0010 / 0.0010 | 0.161 | 0.817 | 0.608 | 0.626 / 0.225 / 0.053 |
| D | s1c | 1715 | 0.693 / 0.660 | 0.0010 / 0.0010 | 0.349 | 0.869 | 0.734 | 0.697 / 0.157 / 0.052 |
| D | ours | 1837 | 0.540 / 0.654 | 0.0010 / 0.0010 | 0.195 | 0.768 | 0.739 | 0.668 / 0.165 / 0.061 |

VERDICT (arm D, S1C famCN, held-out): REFINES  [Jaccard 0.660 vs A 0.473; clean false-multi 0.157 vs A 0.303; p 0.0010]

NPIP side, arm D (S1C | ours famCN):
- PKD1P6-NPIPP1: Soto ID_149, ID_151, ID_152, ID_153, ID_154, ID_155 | D/s1c ID_149, ID_151, ID_152, ID_153, ID_154, ID_155 | D/ours ID_149, ID_153, ID_154, ID_155
- AC126755.6: Soto ID_149, ID_151, ID_152, ID_153, ID_154, ID_155 | D/s1c ID_149, ID_151, ID_152, ID_153, ID_154, ID_155 | D/ours ID_149, ID_152, ID_153, ID_154, ID_155
- MSTRG.2119: Soto ID_41, ID_151, ID_152, ID_153, ID_169 | D/s1c ID_41, ID_151, ID_152, ID_153, ID_169 | D/ours ID_41, ID_152, ID_153, ID_169
- PDXDC2P-NPIPB14P: Soto ID_149, ID_151, ID_152, ID_153, ID_154, ID_155 | D/s1c ID_149, ID_151, ID_152, ID_153, ID_154, ID_155 | D/ours ID_149, ID_151, ID_152, ID_153, ID_154, ID_155
- AP001120.2: Soto ID_149, ID_151, ID_152, ID_153, ID_154, ID_155 | D/s1c ID_149, ID_151, ID_152, ID_153, ID_154, ID_155 | D/ours ID_149, ID_153, ID_154
