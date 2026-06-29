# VG coherence check + non-linear structure (are we using the real VG?)

`family_def_vg_coherence.py` — per family, align every member to the longest member (minimap2 -c -N 10),
parse ALL blocks, and compute (A) member-side backbone-threading coherence and (B) the NON-LINEAR
signatures the collinear pipeline is blind to. 826 recent-coding families, 2087 member alignments.

## (A) Coherence — mostly redundant with the protein bar
Member-side forward threading <50% in 14% of families — but this is largely **UTR-driven** (cDNA
includes UTRs; the protein whole-gene bar already validated the CDS). So coherence at the cDNA level
adds little over the protein bar; it is not a clean independent bridge detector. (The ~10 true
low-density bridges from the FP audit remain the residual; protein homology already gates them.)

## (B) NON-LINEAR structure — the real finding (VERIFIED)
The whole family-definition pipeline is COLLINEAR (cDNA/protein all-vs-all, co-threading, and we even
dropped reverse-strand edges). A real variation graph models inversions (reversing traversal), tandem
repeats (CYCLES), and exon shuffling (reordered shared nodes). We use none of it. How often do these
actually occur in our families?

| non-linear case | families | % | status |
|---|--:|--:|---|
| INVERSION (≥10% reverse block within family) | 0 | 0% | none at this level (inverted gene-dups are sense-sense in cDNA) |
| TANDEM REPEAT (backbone region matched ≥2×) | 20 | 2% | **REAL — verified** |
| EXON SHUFFLE (≥20% out-of-order blocks) | 16 | 2% | real rearrangement + some repeat-driven |

**Verification (RCF_13, all-LOC tandem-repeat family):**
- Self-alignment of members shows genuine INTERNAL tandem repeats (6 and 4 off-diagonal self-repeat
  blocks) — the repeat is in the sequence, not a family-alignment artifact.
- Building the **actual `vg msga` graph** for the family yields a **CYCLIC graph** (198 nodes, 218
  edges, `cyclic`) — a cycle *is* tandem-repeat structure. The real VG represents what our collinear
  homology flattens.

## Implication
~4% of families carry non-linear (repeat / rearrangement) structure that our collinear definition
flattens to a single backbone. This matters for the thesis: a tandem-repeat COPY-NUMBER difference
between copies is a structural distinguishing feature — a PSV-like signal our SNP/indel-bubble PSV
model misses. Modeling repeats/CNV in the family VG would add copy-distinguishing power (move the
identifiability frontier), and the `vg` toolkit (installed, cyclic-graph-capable) is the right engine.

NEXT: a full non-linear census (vg graph per flagged family, classify cycle vs rearrangement vs
artifact); fold structural (repeat-CNV) variants into the PSV set for copy-assignment.
