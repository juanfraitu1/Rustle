# KEY=unitcover
genes 2334: clean 2071, multi 149; owned duplicons 1546 (of 1546 touched by clean genes); families owning >= 1 duplicon 392

| half | multi genes | mean Jaccard (structure) | mean Jaccard (location baseline) | permutation p | verdict |
|---|---|---|---|---|---|
| dev | 81 | 0.439 | 0.221 | 0.0010 | HOLDS |
| heldout | 68 | 0.473 | 0.277 | 0.0010 | HOLDS |

VERDICT (held-out decides): HOLDS  [dev agrees: HOLDS]
- structure, all 149 multi genes: exact set 0.107, recall 0.555, precision 0.718, empty 1
- location baseline, all 149 multi genes: exact set 0.054, recall 0.286, precision 0.426, empty 0
- clean genes (leave-one-out), 2071: predicted exactly their family 0.546, two or more families 0.303, empty 0.020

NPIP side (ID_149-ID_155):
- PKD1P6-NPIPP1: Soto ID_149, ID_151, ID_152, ID_153, ID_154, ID_155 | structure ID_149, ID_154 | location ID_28, ID_149, ID_154
- AC126755.6: Soto ID_149, ID_151, ID_152, ID_153, ID_154, ID_155 | structure ID_149, ID_154 | location ID_28, ID_149, ID_154
- MSTRG.2119: Soto ID_41, ID_151, ID_152, ID_153, ID_169 | structure ID_41 | location ID_16
- PDXDC2P-NPIPB14P: Soto ID_149, ID_151, ID_152, ID_153, ID_154, ID_155 | structure ID_28, ID_151, ID_154 | location ID_154
- AP001120.2: Soto ID_149, ID_151, ID_152, ID_153, ID_154, ID_155 | structure ID_154 | location ID_155
