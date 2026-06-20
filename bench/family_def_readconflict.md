# Defining a multi-copy family at the RNA level: read-conflict vs sequence similarity

## The question
The de-novo pipeline currently defines a family by **sequence similarity** — an edge between two rep
transcripts iff `contiguous_core_coverage ≥ T_CORE = 0.13` (or, as a contrast, minimizer-Jaccard), then dense
communities of that graph. This conflates true recent paralogs with **domain-sharers** (genes sharing one
exon/domain but not paralogous), uses an arbitrary threshold, and measures a *proxy* for the thing that
actually matters downstream. Do we have a good definition, and can we do better?

## The proposed definition (operational / identifiability-grounded)
A **family is a connected component of the read-conflict graph**: two loci are linked iff reads actually
cross-map between them with tied alignment scores — a read with a primary alignment in one locus and an
**AS-tied secondary** in the other (a genuine alternative placement at a *distinct* locus, not one read
overlapping two nested gene annotations). This is Canzar's multimapping-conflict graph, and it is exactly the
unit the copy-assignment problem operates on. Properties:
- **No arbitrary similarity threshold** — the boundary is the alignment-score tie (a property of the data).
- **Provable decomposition** — by construction reads never cross-map outside their component, so the
  assignment problem separates exactly across families with no information lost.
- It distinguishes paralogs from domain-sharers by the *extent* of cross-mapping (whole-transcript vs a domain).

## The experiment (`family_def_readconflict.py`, `family_def_comparison_figure.py`)
On a Compara-labelled set — 5 confirmed recent-paralog pairs (APOBEC3D/F, RFPL1/2/3, RABL2A/B) and 7
domain-sharer pairs (ASDURF/ASNSD1, CASP8/FLACC1, …, all nested opposite-strand genes sharing one exon) — we
measure, per pair, minimizer-Jaccard (similarity) and the read-conflict count + mutual coverage from `GGO.bam`.

| pair | class | Jaccard | read-conflict | MAPQ-0 signal |
|---|---|---|---|---|
| RABL2A~RABL2B | paralog | 0.313 | **190** (mut_cov 1.00) | 170/200 secondary |
| APOBEC3D~APOBEC3F | paralog | 0.144 | **6** (mut_cov 0.92) | within-family secondaries |
| RFPL1/2/3 (×3) | paralog | 0.13–0.17 | **0** | 0 secondary, 0 MAPQ-0 |
| all 7 domain-sharers | domain-sharer | 0.10–**0.41** | **0** | 0 |

## Findings
1. **Sequence similarity does NOT separate the classes.** The domain-sharer CREB1~METTL21A (Jaccard 0.406)
   scores *higher* than every true paralog (max 0.313). The domain-sharer range (0.10–0.41) encloses the
   paralog range (0.13–0.31). No threshold cuts them — the current definition's core failure, on labelled data.
2. **Read-conflict has zero false positives on domain-sharers** (7/7 = 0). The thing the similarity definition
   gets wrong (merging domain-sharers), read-conflict gets right — because nested genes sharing one exon
   produce no *alternative placement*; a read maps to one locus and merely overlaps two annotations.
3. **Read-conflict flags exactly the families where copy-assignment is needed** (RABL2 190, APOBEC3 6) and is
   **silent on RFPL** — a true paralog whose reads place *uniquely* (0 secondary, 0 MAPQ-0). RFPL is paralogous
   evolutionarily but poses no copy-assignment problem, so excluding it from the resolution-needed set is the
   correct operational behaviour.

## Interpretation
Read-conflict is **not** a paralog classifier and should not claim to be — RNA cannot robustly recover
evolutionary paralogy (established separately; needs protein/genomic orthology). It measures
**identifiability-RELEVANCE**: the set of loci among which reads are genuinely confused. That is a strictly
better definition for this thesis than sequence similarity, because (a) it never false-alarms on
domain-sharers, (b) it has no tuned threshold, and (c) it picks out precisely the families the downstream
machinery (PSV/junction assignment, the identifiability theorem) exists to resolve.

## Secondary-cap exposure (`secondary_cap_exposure.py`)
The BAM was built with `minimap2 -ax splice:hq` (no `-N`) → default **N=5 secondaries (≤6 placements/read)**.
For a co-located array all of a read's placements fall in one region, so records-per-read = multimapping degree;
a spike at 6 = saturation (a read may have MORE unreported placements → potentially missed conflict edges).

| region | multimappers | saturated (≥6) | max degree |
|---|---|---|---|
| TSPY array (>6 copies) | 44 | **36 (82%)** | 6 |
| RBMY1 array | 127 | 9 (7%) | 6 |
| chrY amplicon | 8 | 0 | 2 |
| APOBEC3 / RABL2 / RFPL / DSFAM38 (≤7 copies) | ≤15 | **0 (0%)** | ≤2 |

Verdict: the cap bites only **>6-copy co-located arrays** (TSPY-class). Even there the family is a connected
*component* — each read still carries 6 placements, over-connecting the array — so missing edges don't fragment
it; only edge *weights* and the inclusion of a perpetually-marginal copy are at risk. For every ≤6-copy family
(all validated ones), exposure is zero. → safe to port; flag edge weights as lower bounds in mega-arrays.

## Next step (implementation)
Port the criterion into `detect_edges`: keep POA contiguous-core as a cheap homology *prefilter* (necessary to
cross-map), but **confirm and weight each edge by read cross-mapping** (the AS-tied secondary index already
exists: `tied_secondary_reads`). Families become the conflict-graph components. The portable kernel
(`read_conflict.rs`: AS-tied conflict edges + connected-component families) is TDD'd in the Rust tree; the
remaining integration is plumbing per-locus secondary placements into the detection stage.
