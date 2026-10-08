# PREREG — read-level node splitting, two triggers tested in parallel

**Written 2026-09-21 before any score.** Both options prototyped as a Python read-splitting
preprocessor over the real chr16 BAM (`/mnt/linuxdisk/tmp/o1denovo/chr16.bam`), producing a modified
BAM fed unchanged through the shipped `copy_assign --assemble-only` (JUNCTION_MAJORITY now default) and
the shipped `mcl_families` family-definition pipeline — no Rust source touched at this stage. If either
clears its bar, a proper Rust implementation is the natural next step; if neither does, this closes the
line without having risked a code change.

## Why prototype at the READ level, not as a post-hoc node cut

Register 846 (post-hoc parent-boundary node cut) failed for two specific, mechanistic reasons:
**(a) cutting DOUBLES rather than separates** — 7 of 8 chimeras kept both halves in the same family
component anyway, because the "other half" was independently homologous to the same family through a
different route; **(b) a shorter node is an easier node** — re-scoring the cut piece's coverage against
its own (shorter) length inflated its edges. Splitting at the READ level, BEFORE assembly ever builds one
locus out of the chimeric reads, structurally avoids both: the two resulting fragments never exist as one
object to begin with (no doubling to detect), and each fragment's reads flow through the assembler's
ordinary, unmodified coverage/support accounting from scratch (no denominator to get wrong).

Register 815 (junction-based bridge cut, `RUSTLE_LOCUS_BRIDGE_CUT`) won on development and failed the
gorilla holdout on every metric. Both triggers here are tested on a HELD-OUT chromosome before any
conclusion, per that lesson.

## Trigger 1 — chimeric-bridge junction (reuses `is_chimeric_bridge`'s exact test)

The assembler already has `is_chimeric_bridge` (`denovo_assemble.rs`): a SKELETON is chimeric if it
shares a junction with two other skeletons whose genomic spans are disjoint. Measured this session: this
guard excludes a chimera from READ-SUPPORT POOLING only — it does not stop the chimeric skeleton itself
from being assembled and emitted as one fused transcript, which is exactly what produced NPIPB4+RRN3P1's
shared 65 kb locus, CDR2's 91 kb locus, and PKD1P6-NPIPP1's 36 kb locus this session.

    for a primary read r with intron chain [...,(d,a),...]:
        LET L = the set of OTHER reads sharing a junction with r whose spans lie entirely LEFT of (d,a)
        LET R = the set of OTHER reads sharing a junction with r whose spans lie entirely RIGHT of (d,a)
        IF L and R are both non-empty AND span(L) and span(R) are DISJOINT:
            cut r's alignment at (d,a) into two independent sub-alignments (mirrors `split_mischained_reads`,
            which already does exactly this cut-and-keep-both-flanks mechanic for a different trigger)

Exact, not thresholded, matching `is_chimeric_bridge`'s own design principle.

## Trigger 2 — read-identity turnover ("coverage cliff")

    for a primary read r with intron chain [...,(d,a),...]:
        LET S_before = the set of OTHER reads whose alignment blocks overlap r's exon immediately BEFORE (d,a)
        LET S_after  = the set of OTHER reads whose alignment blocks overlap r's exon immediately AFTER (d,a)
        turnover = 1 - |S_before ∩ S_after| / |S_before ∪ S_after|
        IF turnover >= TURNOVER_FLOOR (0.90, chosen before looking -- near-total read-identity replacement):
            cut r at (d,a), same mechanic as Trigger 1

⚠ Related to, but distinct from, §6u7/§6v0's bridge-FRACTION metric (which scored whether two ALREADY
SEPARATE loci should MERGE, and topped out at AUC ~0.65-0.67, never clearing a positive bar). This is the
mirror question — whether one ALREADY-FUSED skeleton should split — scored on a different statistic
(read-set turnover across a single intron, not bridge-read count between two candidate loci). The two are
not the same test and a prior negative result for the merge direction does not presume this one; stated
here so a negative result is not read as a surprise.

## Substrates

- **Development**: chr16 (heavily used this session; every over-merge case that motivated this thread was
  found here).
- **HELD OUT, run last: chr9** (16 doc mentions this session, never used for any decision; still has a
  built PAF/graph pair from §6u9's cover-and-jn work, so a fresh chr9 build is not required if the same
  loci/genes need re-scoring, though the split test itself needs a fresh chr9 BAM pass).

## Scoring

Both triggers scored against:
1. **§6u7's per-copy node-construction endpoint** (universe fixed on the unmerged output, one-to-one
   max-overlap claim, collisions counted as misses) — this is what over-merge actually costs.
2. **Family-definition F on both truths** (Soto cover, protein referee), using the EXISTING max-overlap
   resolver (never the refuted multi-label one from §6v9).
3. **Register 846's specific failure modes, checked explicitly, not assumed away:**
   - **Doubling check**: for each split candidate, do the two fragments' resulting loci still end up in
     the SAME family cluster anyway (via a different transitive path)? Report the rate.
   - **Coverage-inflation check**: report each fragment's own length vs the original locus's length, and
     confirm no denominator is computed against a fragment's shorter length in a way that inflates its
     apparent identity/coverage relative to the original (this is naturally avoided by construction here,
     since fragments are assembled from scratch rather than re-scored from a pre-existing node, but is
     checked, not assumed).

## Bars — stated before looking

1. **PRIMARY: per-copy correctness must rise by >= +2.0 percentage points on chr16** — the same bar
   §6u7 used, for direct comparability to that refutation.
2. **GUARD: false-merge rate must not rise.**
3. **GUARD: NPIP dominant-cluster coverage must stay >= 20/21** (the exact number the §6p1 distance merge
   broke).
4. **HELD-OUT (chr9): per-copy correctness must also rise, by any positive amount, at the same trigger
   parameters chosen on chr16.** A held-out failure means NOT ADOPTED at the Python-prototype stage,
   full stop — no Rust implementation is proposed from a rule that fails its own holdout.
5. If a trigger clears bars 1-4, doubling and coverage-inflation are reported as a further, non-blocking
   diagnostic (they would only block a subsequent Rust implementation decision, not this prototype's
   verdict).

## What a negative result means

If both triggers fail, it strengthens the case that node-level over-merge from readthrough genuinely
cannot be resolved by ANY local, per-junction rule (post-hoc cut, pre-assembly read split, or merge-side
bridge-fraction) — three structurally different attack angles, all refuted, would be a strong, specific
finding in its own right, not just an absence of success.

---

## OUTCOME (appended 2026-09-21, after scoring)

**REFUTED on the primary bar for both triggers.** Per-copy correctness on chr16 (universe ~1,327-1,332,
baseline 705/1327 = 0.5313):

| arm | correct | rate | delta | collisions | false-merge |
|---|---|---|---|---|---|
| baseline | 705 | 0.5313 | — | 212 | 26.6% |
| Trigger 1 (chimeric-bridge) | 705 | 0.5313 | **+0.00pp** | 212 | 26.6% |
| Trigger 2 (turnover) | 718 | 0.5390 | **+0.78pp** | 207 | 26.4% |

Bar was **>= +2.0pp**. Trigger 1 moved nothing at all; Trigger 2 moved in the right direction on every
guard (correctness up, collisions down, false-merge rate down) but well short of the bar. Pooled
family-definition F (both truths) was byte-identical to baseline for both triggers — neither split
touched a locus that changes which cluster a truth gene resolves to.

**Trigger 1 flagged 63 reads genome-wide, Trigger 2 flagged 1,926 — and NEITHER flagged a single read at
any of the three known over-merge sites this whole investigation is about** (NPIPB4/RRN3P1's shared
locus, CDR2's 91 kb engulfing locus, PKD1P6-NPIPP1). Checked directly, not assumed.

## Why Trigger 1 cannot work here — confirmed mechanistically, not just observed

Checked one of the 66 individual reads that independently span CDR2's full 91,627 bp engulfing region
(8-10 introns each, not an outlier — 66 separate reads show this same wide structure). Its own junctions
are shared by 7 to 777 other reads each. Merging all 560 distinct neighbour spans collapses to **ONE**
group, not two: every "side" of the apparent bridge is itself connected to the other side by yet more
reads of intermediate, overlapping span (tandem, physically adjacent genes in a segmental duplication).
`is_chimeric_bridge`'s design assumes an OUTLIER read bridges two otherwise-separate, cleanly-spanning
populations — that assumption fails whenever the "bridge" itself is the well-replicated, dominant
structure at a locus, which register 852/[[project_pkd1p6_npipp1_is_real]] already established is the
normal case here, not the exception (PKD1P6-NPIPP1: 110 MAPQ-60 reads, a real fusion, not rare artifact).

**This is register 846's finding (a) — transitivity defeats separation — reproduced at a completely
different pipeline stage.** 846 found post-hoc node-cutting fails because both cut halves stay in the
same family component via transitive edges elsewhere. Here, PRE-assembly read-level bridge detection
fails for the identical structural reason one stage earlier: there is no clean disjoint pair to detect in
the first place, because overlapping intermediate reads chain the whole region into one transitively
connected population. Two structurally different attack angles (post-hoc node cut, pre-assembly read
split), at two different pipeline stages, both defeated by the same underlying biology.

## Held-out substrate

Not run. Bar 4 required clearing the primary development bar first ("at the SAME trigger parameters
chosen on chr16"); neither trigger did, so there is no rule to carry to chr9.

## Conclusion

Three independent, structurally different attempts at fixing readthrough-driven over-merge are now
refuted on this substrate: register 846 (post-hoc parent-boundary node cut), register 815 (junction-based
bridge cut, won on dev/failed gorilla holdout), and this session's read-level chimeric-bridge and
turnover triggers (failed on development itself, with a confirmed mechanism). The common thread across
all three: **whatever signal is used to identify "this should be split" is looking for an outlier or a
sharp discontinuity, and the dominant over-merge cases on this substrate are neither** — they are
well-supported, densely-replicated readthrough transcription through tandemly duplicated, physically
overlapping genes, which is nearly indistinguishable, read by read, from one long real transcript.
**A future attempt needs a signal that does not depend on rarity or population-level asymmetry** — for
example, external evidence (annotation-independent expression breaks, or genuinely orthogonal molecular
evidence) rather than anything derivable from the read population's own internal structure, since that
structure is, at these loci, exactly what a genuine multi-gene readthrough looks like.
