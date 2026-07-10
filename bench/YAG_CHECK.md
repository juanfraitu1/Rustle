# Y ampliconic genes (IsoCon's benchmark): how many copies do we recover?

**Date:** 2026-07-09. **Substrate:** `GGO_mm.bam`, gorilla chrY = **`NC_073248.2`** (confirmed by SRY, ZFY, UTY,
USP9Y). Run with `--min-copies 2 --skip-poa-diagnostic --homology-primary`, one region at a time, foreground.

IsoCon (Sahlin 2018) was built on the human Y ampliconic gene families (TSPY, RBMY, DAZ, CDY, HSFY, BPY2 …) and
reported more copies/isoforms than were previously known. Those families are annotated in gorilla only by
`description=`, not by symbol, so they were recovered from the GFF description field.

## Coverage first — three of six families have no data

| family | annotated copies | window | primary reads |
|---|---|---|---|
| DAZ | 2 | `42778133-42950552` | 220 |
| BPY2 | 2 | `42640394-43088291` | 224 |
| RBMY (distal cluster) | 6 | `19597754-19735926` | 87 |
| RBMY (proximal) | 6 | `15472440-16459625` | 77 |
| TSPY | 6 | `34731504-34847734` | 34 |
| HSFY | 5 | `17247746-17494787` | **10** |
| CDY | 1 | `19341920-19354346` | **0** |

CDY and HSFY carry no usable data in this sample, so nothing can be said about them. TSPY is being asked to
resolve **6 copies from 34 reads.**

## Result

| family | annotated copies in window | de-novo reps | **copies recovered** | assigned / tied / ambiguous |
|---|---|---|---|---|
| **RBMY (distal cluster)** | 6 | 5 | **6 — exact match** | 503 / 215 / **170** |
| **TSPY** | 6 | 4 | **5** (4 assembled + 1 rescued `RC_`) | 105 / 113 / 0 |
| RBMY (proximal) | 6 | 2 | 2 | 243 / 23 / 0 |
| **DAZ** | 2 | 2 | **0 families** | — |
| **BPY2** | 2 | 2 | **0 families** | — |

**RBMY's distal cluster is an exact hit**: six annotated copies, six recovered, from 87 reads — including copies
the assembler never built (5 reps → 6 copies, the sixth recovered by collapsed/rescue detection). This is the
regime the method was designed for.

**TSPY recovers 5 of 6.** The missing copy is `34777900-34780649`. One of the five (`RC_..._34757175`) was
recovered by rescue, not assembly. 52% of placements are tied — with 34 reads over 6 near-identical copies, most
reads simply never reach a distinguishing column. That is the coverage cause of abstention, not the K=0 wall.

**RBMY proximal recovers 2 of 6** — 77 reads spread over a 1 Mb span; four copies are simply not transcribed
deeply enough to assemble.

⚠ **RBMY distal reports 170 ambiguous reads (19%).** Genome-wide the ambiguous rate is 0.05% (20 reads in
43,239). RBMY is the one place the gate genuinely cannot decide, rather than certifying a tie. Unexplained;
worth a look.

## DAZ and BPY2: a readthrough artifact is *hiding a real copy*

Both windows return **zero families**, and for the same reason. Every assembled spliced transcript in the DAZ
window starts at ~42783130 — they are all isoforms of **DAZ1** (the `-` copy). The second DAZ copy
(`42879918-42945552`, `+`) yields **no spliced model at all**. The only `+`-strand object is a **164 kb
single-exon readthrough** (`42781428-42945549`, 12 reads) spanning *both* DAZ copies. Because a readthrough's
exon-sum is genomic rather than mRNA, the homology oracle finds **0 edges**, and the family collapses.

The BPY2 window is the same picture, worse: its only `+`-strand model is a **297.8 kb single-exon readthrough**
(`42781428-43079270`, 14 reads) running from DAZ through BPY2. Both annotated BPY2 copies produce nothing.

So the readthrough artifact is not merely noise that inflates abstention — **it stands in the place of a real
copy and suppresses the family.** Dropping it (below) does not by itself recover DAZ2 or BPY2; that needs
assembly at low coverage, which is precisely IsoCon's contribution.

## Rule work: R4 confirmed on two new artifacts; a new candidate rejected

The readthrough rule validated in [READTHROUGH_RULE_VALIDATION](READTHROUGH_RULE_VALIDATION.md) is:

> a single-exon transcript is the unspliced form iff **≥ 5 distinct junctions, each with ≥ 2 reads, lie entirely
> inside its span**.

Both YAG readthroughs are caught, and not marginally:

| artifact | span | distinct junctions | R4 (≥5) |
|---|---|---|---|
| DAZ-window readthrough | 164.1 kb | **56** | flagged |
| BPY2-window readthrough | 297.8 kb | **56** | flagged |

That raises the artifact set to **15, minimum 14 distinct junctions**, against a control maximum of 4. The margin
holds on an independent chromosome and gene set.

### A new rule was derived and it FAILS: exonic coverage breadth

Candidate **R8**: *a real intronless mRNA is uniformly covered; a readthrough spans introns where coverage
collapses.* Measure the fraction of positions in the single exon covered by ≥ 2 primary reads.

| object | breadth (≥2×) |
|---|---|
| DAZ readthrough | 0.086 |
| RFPL 128 kb | 0.160 |
| artifact (`182489940-182666627`) | 0.332 |
| artifact (`72892769-72912934`) | 0.426 |
| **GSTM readthrough** | **0.785** |
| **TSPYL1** (real intronless gene) | **0.464** |
| GSPT2 (real intronless gene) | 0.981 |
| ATXN7L3B (real intronless gene) | 0.875 |

The GSTM readthrough (0.785) is better covered than the real intronless gene TSPYL1 (0.464). **The classes
overlap; breadth cannot separate them at any cutoff.** Rejected.

Five candidate rules have now been tested (assembled-intron, any-read-junction, span/read ratio, contains-a-
transcript, coverage breadth). **Only distinct-junction count survives.** The structural `catalog_overlaps`
check (`Containment` / `SharedAcrossFamilies`) remains the belt-and-braces safety net for anything the rule lets
through.

## Still open

- **No retrocopy positive control.** CDY — the one intronless, retrogene-derived YAG, and the object R4 is most
  likely to destroy — has **0 reads** here. The gap flagged in `READTHROUGH_RULE_VALIDATION.md` is still open;
  chrY does not close it.
- DAZ2 / BPY2 need assembly at low coverage, not a better filter.
- RBMY's 19% ambiguous rate is unexplained.

## Reproduce

```
copy_assign --bam GGO_mm.bam --fasta GGO.fasta --region NC_073248.2:19597754-19735926 \
    --min-copies 2 --skip-poa-diagnostic --homology-primary --out RBMY
```
Related: `bench/READTHROUGH_RULE_VALIDATION.md`, `bench/FAMILY_SPOT_CHECK.md`, `reference_isocon_sahlin`.
