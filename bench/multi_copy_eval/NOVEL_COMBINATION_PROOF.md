# Can we PROVE rustle finds novel cross-copy exon combinations "impossible otherwise"?

Short answer: **for paralog copies, no — and the reason is fundamental (a measured data-level
barrier), not a rustle limitation.** The graph's genuinely-unique power is on the *sequence /
read-assignment* axis (gene-conversion detection, multimapper rescue), **not** the
novel-junction axis. This document records the discriminator, the measurements, and the
mechanism.

## The discriminator (what "impossible otherwise" must mean)

An exon combination E (an ordered intron chain) is graph-unique — unrecoverable by any
linear / primary-only full-length assembler — only if **all** hold:

- **D1** no single read's aligned chain spans E (else a full-length/read-chain linear tool emits it);
- **D2** E is assemblable only from evidence a primary-only tool mis-routes or discards;
- **D3** E is *not* a parsimony walk over one copy's own splice graph — at least one junction of E
  must come from a **different physical copy** (a cross-copy chain).

Against this bar, every prior result falls short:

| result | graph-unique? | why |
|---|---|---|
| 1668 class-`j` novel junctions | **no** | 83% spanned by a *primary* read (fail D1); the rest are ordinary path enumeration (fail D3). A real *discovery* story, not impossibility. |
| gene-conversion mosaic (CONFIRMED) | **yes, but** | genuinely graph-unique, but it is *detection* (a GTF attribute, `gtf.rs:242-246`) and has **zero novel junctions** — same intron chain. Not an exon combination. |
| RABL2 9-vs-4 | **no** | real graph-level multimapper *rescue*, but each isoform is structurally single-copy and read-spanned (fail D1/D3). Proves copy-disambiguation. |

## The measurement: cross-copy recombinants never produce an observable novel junction

A cross-copy chain (D3) needs a molecule that joins copy-A exons to copy-B exons. We aligned such
recombinants (`minimap2 -ax splice:hq`, secondaries+supplementary on) across every regime:

| regime | setup | result |
|---|---|---|
| **colinear** (substitution divergence) | 2 copies 200 kb apart, identity 0.95 / 0.97 / 0.98 / 0.99 | **0/40 split** at every identity. All align contiguously to ONE locus with that locus's standard 4-intron chain `400M2600N…450M`; the foreign half is just mismatches → a **SNP mosaic, no novel junction**. |
| **inversion** | copy1 an inverted duplication (− strand) | **0/40 split.** Each read gets *two full contiguous* alignments (one per locus, opposite strands); minimap2 picks one, never splits. |
| **tandem read-through** | copy0 e0,e1 → copy1 e3,e4, copies 9.5 kb apart | **0/20 split.** Aligns contiguously with a long intron `…300M11900N450M` — one colinear alignment a linear splice graph also sees. |

**Mechanism (why it's fundamental):** paralog copies are *homologous* — every copy's locus
contains a counterpart of every exon. So a cross-copy mRNA recombinant **always has a contiguous,
colinear alignment available at a single locus**, and the aligner always prefers it over a
split (a 200 kb "intron" + supplementary is far more expensive than a few % mismatch). The
structural recombination is therefore **never observable as a novel junction** — it collapses to
a sequence (PSV) mosaic on top of one copy's ordinary junction chain.

**The formal trap (D1 ∧ D3 is self-defeating for paralogs):** the only thing that could assert a
cross-copy junction adjacency is a molecule spanning it. If such a molecule exists, it aligns
contiguously and a read *does* span the chain → **fails D1** (linearly findable). If you fragment
reads so none spans it → **nothing asserts the adjacency** → it is unassertable. There is no
middle. So a cross-copy *novel junction* is either linearly recoverable or unobservable.

(The only true escape — a forced split with a genuinely novel junction — needs **non-homologous**
partners, i.e. a gene fusion / trans-splice. That is fusion detection from a split alignment,
which a fusion-aware *linear* tool also does; it is neither "paralog multimapper" nor graph-*unique*.)

## What IS graph-unique and realizable (the honest win)

The variation graph's unique power is **not** inventing structurally-novel junctions — any
observable junction lives in some read, so a linear/fusion-aware tool can reach it. Its unique
power is on the **sequence / assignment axis**:

1. **Gene-conversion / cross-copy *sequence* recombinant detection** — `RUSTLE_VG_MOSAIC_ON`,
   CONFIRMED on synthetic + realistic badread HiFi (breakpoint 9655, 0 bp dispersion, invisible
   to StringTie). Graph-unique because a linear assembler must assign a read to one copy and
   cannot represent a within-read copy switch. **But:** detect-and-flag, **no novel junction**.
2. **Multimapper copy-rescue** — RABL2 9-vs-4: reassigning discarded secondaries by edit-distance
   ownership recovers copy-specific isoforms StringTie drops. Real, verified, but single-copy
   structures.

Both resolve *which copy* the evidence belongs to. Neither manufactures a novel exon combination.

## Answering "any other way to model alternate exon donors and acceptors?"

Rustle **already models alternate donors/acceptors as first-class graph variation *per locus***
(`graph_build.rs:957-1015` group-by-acceptor for alt-donors, group-by-donor for alt-acceptors;
every distinct splice coordinate is a node split-point + distinct edge), plus an opt-in
within-locus enumerator (`alt_donor_acceptor_rescue`, `RUSTLE_ALT_SPLICE_RESCUE`).

The **cross-copy** gap is real (the family graph's `JunctionEdge` carries no donor/acceptor
coordinate; it rounds boundaries to 10 bp and keys edges by node pair, collapsing distinct
boundaries). Ranked options to model cross-copy alt-boundaries:

1. **Promote the gene-conversion detector to EMIT** the recombinant as a flagged isoform (option 3)
   — the one with a *demonstrated* graph-unique signal. **Sequence-channel, not structural.**
2. Coordinate-labelled family edges (anchor-relative alt-boundary bubbles) + seq-to-graph read
   voting — makes cross-copy alt-boundaries representable and read-confirmable.
3. Base-level POA variation graph over the homologous exon — the full vg-style endgame.

**Honest bound:** by the barrier above, none of these recovers a *structurally* novel cross-copy
junction for paralogs (none is observable); and the borrow/completion machinery is measured INERT
on real loci (divergence that starves a copy also kills the linkage). The realizable cross-copy
payoff is the **sequence channel** (option 1), and the unbounded payoff is the **per-locus**
alt-donor/acceptor modeling rustle already has.

## Verdict

> "Find unique exon combinations impossible otherwise" conflates two axes. On the **structural
> (novel-junction)** axis, it is **not provable for paralogs** — homology guarantees a contiguous
> colinear alignment, so a cross-copy novel junction is either linearly findable or unobservable
> (measured: 0 splits across colinear/inverted/tandem, identity 0.95–0.99). On the **sequence
> (copy-assignment)** axis, rustle's graph **does** do what linear tools cannot — detect cross-copy
> gene-conversion recombinants and rescue copy-specific isoforms from discarded multimappers — but
> these carry no novel junction. The graph's edge is *resolution*, not *invention*.
