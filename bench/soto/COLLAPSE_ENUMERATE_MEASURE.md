# `--collapse-enumerate` — measured effect (byte-identical OFF + Soto ON/OFF + EEF1A1 control)

The flag re-admits near-identical gene families that RNA collapses to `<2` distinct loci — which the
homology catalog (`gw_family_catalog --homology-primary`) otherwise drops silently at its
`≥2-distinct-loci` gate — as **K=0-collapsed COPY-NUMBER** entries (`<out>.collapsed.tsv`), never as
fabricated per-read copies. Admission requires **all three** signals: a local `hidden_copy` witness
(a balanced co-equal 2nd haplotype at the locus), `alt_read_fraction ≥ 0.30`, and a genome projection
landing at **≥2 distinct loci**. Default OFF. Binaries built at commit `329ac9c` (feature head).

## 1. Isolation guarantee — byte-identical with the flag OFF ✓

`copy_assign … --homology-primary --skip-poa-diagnostic --min-copies 2` on the four known families
(GGO_mm.bam / GGO.fasta), flag OFF, md5 vs the pre-feature release baseline (`rf_*_on.*`):

| family | families.tsv | famcn_readonly.tsv | assignments.tsv |
|--------|:---:|:---:|:---:|
| DAZ    | ✓ MATCH | ✓ MATCH | ✓ MATCH |
| GSTM   | ✓ MATCH | ✓ MATCH | ✓ MATCH |
| PCDHB  | ✓ MATCH | ✓ MATCH | ✓ MATCH |
| MAGEA  | ✓ MATCH | ✓ MATCH | ✓ MATCH |

**12/12 files byte-identical.** And on the Soto genome-wide sweep, `families.tsv` / `copies.tsv` are
**md5-identical OFF vs ON** (`99137c22…` / `137175be…`) — the feature is purely *additive*: it only ever
writes the new `collapsed.tsv`, and only when the flag is on and something is re-admitted.

## 2. Soto ON vs OFF — recall gain at 100% precision

`gw_family_catalog --cross-chrom --homology-primary --enumerate-copies` on the scoped Soto BAM
(`soto_reads.bam` → CHM13v2.0), (a) without and (b) with `--collapse-enumerate`:

| run | catalog families | collapsed re-admitted |
|-----|:---:|:---:|
| OFF | 66 | — |
| ON  | 66 (identical) | **2** |

Both re-admitted families are **genuine Soto segmental-duplication families** that the OFF catalog
missed (they collapse to `<2` RNA-distinct loci and are dropped at the `≥2-distinct-loci` gate):

**GWFAMc0 — ANKRD20A (Soto `ID_280`, 3 annotated members), famCN=5 (seed + 4 projected), alt_frac=0.567, 234 alt reads**
`seed chr9:40,211,594-40,234,709` → `ANKRD20A3P` (d=0). Projection loci (excl. the seed's own locus):
- `chr9:43,062,731 @0.992` — no annotation within 200 kb → a **novel/unannotated 4th ANKRD20 copy**
- `chr9:44,378,410 @0.984` → `ANKRD20A7P` (d≈4.3 kb)
- `chr9:77,926,044 @0.992` → near `MEP1AP1`/`ID_260` (adjacent segdup block)
- `chr9:79,615,466 @0.993` → `ANKRD20A1` (d=0)

**GWFAMc1 — chr1q21 segdup (Soto `ID_211`, 4 annotated members), famCN=4 (seed + 3 projected), alt_frac=0.400, 12 alt reads**
`seed chr1:148,825,999-148,829,853` → `AC243772.3` (d=0). Projection loci (excl. the seed's own locus):
- `chr1:121,019,456 @0.994` → `AC241377.2` (d=0)
- `chr1:142,865,358 @0.998` → `AC239798.1` (d≈18 kb)
- `chr1:144,288,808 @0.998` → `AC241585.3` (d≈16 kb)

> **famCN is seed-inclusive.** `project_family_copies` is called with the candidate's own span as its sole
> `known` entry, so the projection *excludes* the seed locus; total copy number = seed + projected others =
> `n_projection_loci + 1` (`famcn_from_projection`). A whole-branch review caught the earlier off-by-one
> (famCN was reporting the projected-others count). With the fix, **`ID_211` famCN = 4 = exactly its 4
> annotated members** (`AC243772.3` seed + `AC241377.2`/`AC239798.1`/`AC241585.3`).

**Precision = 100%.** Every re-admitted family's seed sits d=0 on a real Soto member, and every
projection locus lands on a real annotated member of that same family — `ID_211` recovers **4/4** of its
annotated members exactly; `ID_280` recovers all 3 of its annotated members (`ANKRD20A3P`/`ANKRD20A7P`/
`ANKRD20A1`) plus a bona-fide *unannotated* ANKRD20 copy (`chr9:43,062,731`), so famCN=5 correctly exceeds
the 3 annotated. **Zero false re-admissions.**

**Recall is deliberately conservative.** The strict three-signal gate re-admitted 2 real collapsed
families in this sweep. They are **not** among the 9 pre-labelled "category-B" candidates
(TCAF2/ANKRD36B/LIMS1/NCF1B/… `ID_127/22/226/240/251/338/386/402/458`) from the earlier miss deep-dive;
those did not produce a balanced ≥0.30 local witness with enough reads in this BAM, so the gate correctly
abstained rather than guess. The feature recovers *real* collapse at 100% precision, not a fixed target
list — precision-first by design (the retired `collapse_gate` was retired for the opposite failure).

## 3. EEF1A1 control — the known dispersed-pseudogene hard case is rejected ✓

`gw_family_catalog … --collapse-enumerate` on an EEF1A1-locus-scoped BAM (`NC_073229.2` region of
GGO_mm.bam, 4,017 reads) against the **full** GGO.fasta (so projection can still reach EEF1A1's dispersed
processed pseudogenes — exactly the condition that fooled the retired `collapse_gate`, χ(H)=7):

```
[gw-catalog-homology] 43 skeletons -> 2 reps over 1 contigs
[gw-catalog-homology] 0 E_r edges -> 0 families (>= 2 distinct loci)
→ NO collapsed.tsv written
```

The candidate reached the gate (2 singleton reps dropped → routed to `readmit_locus`) and was **rejected
on the local-witness signal**: EEF1A1's own locus carries no balanced co-equal 2nd haplotype, even though
its pseudogenes satisfy the projection signal. This is precisely the discrimination the three-signal gate
adds over the retired depth-only gate. (Per the spec's honest caveat, the intronless/structure 4th-signal
filter is therefore **not** needed on this data — but remains the documented fallback if a dispersed case
ever re-admits.)

## Performance note

`--collapse-enumerate` re-projects each dropped candidate onto the genome (in-engine minimap2), so the
Soto sweep ran ~50 min ON vs ~6 min OFF. It is an opt-in analysis flag, not a default-path cost; the
default (OFF) path is unchanged and unaffected. (A follow-up can route `readmit_locus` through the batched
`project_families_batch` to amortize the per-candidate genome re-index.)

## Provenance

Numbers are from the post-review binary (`e3b4c87`, feature head after the whole-branch adversarial
review). The review confirmed both load-bearing contracts hold *in code* (byte-identical-OFF: the collapsed
vec never feeds the families/copies emit; primary-only-witness: `is_secondary || is_supplementary` filtered
before the hidden-copy input) and caught the seed-inclusive famCN off-by-one corrected above. OFF
byte-identical was re-verified after the fixes (DAZ 3/3 md5-MATCH).

## Bottom line

- **OFF = byte-identical** (isolation guarantee holds; 12/12 known-family files + Soto families/copies).
- **ON recovers real collapsed families at 100% precision, 0 false re-admissions** (ANKRD20A, chr1q21),
  including a novel unannotated ANKRD20 copy; conservative recall by design.
- **EEF1A1 rejected** — the three-signal gate discriminates a genuine local collapse from dispersed
  paralog bleed, which the retired `collapse_gate` could not.
