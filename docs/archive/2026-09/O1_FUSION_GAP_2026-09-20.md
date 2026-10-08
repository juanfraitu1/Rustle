# Where O1 is incomplete, and what fusions specifically cost it

Measured 2026-09-20 on the 9-chromosome panel built for the held-out test (chr2/5/6/7/8/10/16/17/21).
Descriptive: no rule changed, nothing adopted.

## 1. The structural fact: O1 emits a PARTITION

**0 of 2,670 clustered loci appear in more than one cluster.** The DNA catalog is a strict partition —
every gene gets exactly one family. A readthrough/fusion gene is, by construction, a member of **two**
families, so the output type cannot represent it. The method must pick a side, and it does.

## 2. What that costs, measured

RefSeq curates the label (`description=readthrough`): **209 such genes genome-wide, 98 on this panel.**
Of those, 19 land in a cluster at all, in 11 clusters — and **8 of the 11 also contain a normal gene**,
i.e. the fusion has pulled in (or been pulled into) one parent family:

| chromosome | cluster | the fusion(s) | clustered with |
|---|---|---|---|
| chr10 | MCL1 | BMS1P2-AGAP9, BMS1P4-AGAP5, PARGP1-AGAP4, TIMM23B-AGAP6 | AGAP7P/10P/11/12P/13P/14P |
| chr16 | MCL2 | PKD1P3-NPIPA1, PKD1P4-NPIPA8, PKD1P6-NPIPP1, PKD1P5-… | PKD1P2 |
| chr16 | MCL11 | BOLA2-SMG1P6 | BOLA2B |
| chr17 | MCL39 | TBC1D3P1-DHX40P1 | DHX40 |
| chr17 | MCL13 | ZNF286A-TBC1D26 | TBC1D28 |
| chr2 | MCL9 | INO80B-WBP1 | WBP1P1, WBP1P2 |
| chr7 | MCL58 | GIMAP1-GIMAP5 | GIMAP3P |

**The discarded half is often the bigger one.** Genomic alignment (PAF match bases, ≥ 300 bp) from each
fusion into its own cluster vs the cluster it was not placed in:

| fusion | placed in | own-cluster bp | bp into the OTHER parent's cluster |
|---|---|---|---|
| PKD1P6-NPIPP1 | MCL2 (PKD1P) | 676,288 | **1,712,037 → MCL0 (NPIP)** — 2.5× more |
| PKD1P4-NPIPA8 | MCL2 | 616,987 | **1,200,730 → MCL0** |
| PKD1P3-NPIPA1 | MCL2 | 664,337 | **899,219 → MCL0** |
| BOLA2-SMG1P6 | MCL11 (BOLA2) | 2,988 | **105,974 → MCL1 (SMG1)** — 35× more |

Each fusion has ≥ 300 bp alignment into **5–15 different clusters**; exactly one survives.

⚠ Read this carefully: the placement is not *wrong*. Under a partition, keeping PKD1P6-NPIPP1 out of
NPIP is the defensible choice — it stops fusions polluting the NPIP family, which is the §6m1/§6q
boundary story. The point is that **the choice is forced**: you must either pollute NPIP or discard the
fusion's NPIP membership. There is no third option inside a partition.

## 3. Why this is not simply "re-propose the cover"

Two adjacent things were already pre-registered and refuted, and both close doors:

- **r846 / §6m3** — splitting a chimera into two nodes at the parent boundary, one per family. Safe but
  **harmful**: it DOUBLES rather than separates, and the short pieces become hubs.
- **r845 / §6m2** — a `--chimera-policy` lever including **`multi` (a chimera belongs to both parent
  families)**. `multi` never cleared the bar (max +0.018) and was **negative at NPIP L3 (−0.006)**.

⭐ But r845's own root cause is the opening, and it is on the **prediction** side, not the truth side:
its verdict reads *"a second truth label **the prediction cannot match** costs more than it buys"*, and
it attributes the null to set size — *"only 1 of 27 NPIP members and 0 of 19 TBC1D3 members are flagged
— too small a set to move a family metric even with a perfect label"*. It closes with an explicit
condition: **"Do not re-propose as a default without a family whose chimera fraction is large."**

That family now exists, and it was found on a near-zero-exposure chromosome:

| family | members | curated readthroughs | chimera fraction |
|---|---|---|---|
| **chr10 MCL1 (AGAP)** | 10 | **4** | **40.0%** |
| NPIP | 27 | 1 | 3.7% ← r845's root cause |
| TBC1D3 | 19 | 0 | 0.0% |

r845 changed what the **truth** may assert while the prediction stayed a partition; nothing has yet
tested letting the **prediction** carry dual membership, and the AGAP cluster is the first substrate
where the label set is large enough for that to be measurable at all.

## 4. The other, larger incompleteness — not fusions

Fusions are a small population. The dominant O1 gap is elsewhere and already on record:

- ⭐⭐⭐ **§6o8: the RNA edge graph caps pairwise recall at 0.052 (human) / 0.229 (gorilla), because
  > 50% of truth families have NO EDGE AT ALL.** Grouping is saturated; no clustering operator can fix
  it. **Edge construction, then node completeness, is the priority** — not the grouping rule.
- **De novo purity**: recall is 1.000 (21/21 NPIP genes covered) but the median gene is split across
  **2 de novo loci**, and 37 clustered loci touching 21 genes land in 9 clusters.
- **Subfamilies**: all 21 NPIP genes land in ONE cluster, so **NPIPA is not separated from NPIPB**
  (advisor Q9), and under a per-subfamily truth NPIPA scores as wholly missed.
- **Precision is the weak side on held-out data**: chr8 sensitivity 0.950 but precision 0.583.
