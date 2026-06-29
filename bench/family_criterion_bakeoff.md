# Conflict-edge tie criterion: AS vs NM vs de (the advisor's "don't trust AS" test)

## Question
The conflict-graph family definition links two loci iff `>= MIN_READS` reads cross-map between them with a
**tied** placement. The tie was defined on minimap2's **AS** (alignment score). The advisor's objection: AS is
a length + gap-model composite — the aligner's best guess, not raw evidence — and should not be taken at face
value. Do we get a better, threshold-defensible family definition by tying on raw-er signals: **NM** (edit
distance) or **de** (minimap2's gap-compressed per-base divergence)?

## Method (`family_criterion_bakeoff.py`)
Per labelled pair, collect every read's alignments in both member spans, attribute each **record** to its
best-overlap member (the same fix that killed the nested-gene FP in the Rust pipeline), and test the read's two
placements for a tie under each criterion:
- **AS-tie**: `min(AS) >= 0.9·max(AS)` (the shipped Rust default)
- **NM-tie**: `|NM_a/alnlen_a − NM_b/alnlen_b| <= D` AND `max <= CEIL`
- **de-tie**: `|de_a − de_b| <= D` AND `max(de_a,de_b) <= CEIL`

Edge iff conflict-read count `>= MIN_READS`. Panel: 5 Compara paralog pairs + 7 domain-sharers (labelled
ground truth) + 7 cross-chromosome paralogs (5–90 % identity) + 3 co-located arrays (MAGEA/PRAMEF/XAGE).
Verified by 5 adversarial sub-investigations (workflow `wf_01631729-ab7`); script audited sound.

## Headline finding — AS has a systematic false-positive mode on retrocopies
On **diverged processed-pseudogene / retrocopy** pairs, AS gives *near-identical scores* to placements that NM
and de show are 10–20× apart in quality, because the read's AS at the spliced **parent** pays per-intron
gap-open penalties that the single-exon **retrocopy** alignment never incurs — so AS ranks the worse copy
at-or-above the true parent and fires a spurious conflict edge:

| pair (parent ~ retrocopy) | AS edges | de/NM edges | median de gap | parent favoured |
|---|---|---|---|---|
| EEF1A1 ~ LOC109023808 | **3347 (FP)** | 0 | 0.0168 (~16×) | 3347/3347 |
| CNN2 ~ LOC129524764 | **121 (FP)** | 0 | 0.0084 (~10×) | 121/121 |
| ADH5 ~ LOC101125574 | **3 (FP)** | 0 | 0.0296 | 3/3 |

AS scores the *worse* copy strictly higher in 71 % (EEF1A1) / 90 % (CNN2) of reads. These reads are cleanly
resolvable to the intron-bearing parent (parent reads carry 6–8 CIGAR `N`-gaps the copy lacks); the AS tie is
an artifact, not ambiguity. This is the advisor's thesis proven on real data — and it matters because
retrocopies/pseudogene copies are a large slice of the copy landscape (1,313 pseudogene + 665 unannotated
copies genome-wide).

## Corrected confusion matrix (operating point D=0.005, CEIL=0.05, MIN_READS=3)
Negatives (should NOT fire): 7 domain-sharers + 3 retrocopies + APOBEC3 (asymmetric, minimap2-resolved 15/15) = 11.
Positives (should fire): RABL2A~RABL2B, AK6~LOC (jac .73), CCDC196~LOC (jac .90) = 3. RFPL excluded (0 cross-mapping reads — silent by absence).

| criterion | TP | FN | TN | FP | precision | recall |
|---|---|---|---|---|---|---|
| **AS** | 3/3 | 0 | 7/11 | **4** (EEF1A1, CNN2, ADH5, APOBEC3) | **0.43** | 1.0 |
| **NM** | 3/3 | 0 | 11/11 | 0 | 1.0 | 1.0 |
| **de** | 3/3 | 0 | 11/11 | 0 | **1.0** | 1.0 |

**`de` wins** (NM corroborates, AS demoted to audit). de is clip-robust where NM is not (APOBEC3D NM-rate
inflates to 0.16–0.22 on soft-clip/poly-A tails while de stays ~0.02), and has the wider safety margin
(CCDC196 fires 24 under de vs 16 under NM). Crucially **de-tie ⊂ AS-tie** on every pair — de removes exactly
AS's 4 false positives — so `de-tie is a strict subset of AS-tie` is the portable regression check.

## Adjudications
- **APOBEC3 not-firing is correct, not a FN.** Its cross-mapping reads fit APOBEC3D 2–3× better than APOBEC3F
  (median |de gap| 0.024; 0/15 reads fit both ≤5 %; minimap2 resolves 15/15 at MAPQ60-vs-MAPQ0). Contrast the
  genuine RABL2 fire: median |de gap| 0.006, 187/195 reads fit both ≤5 %, 118/195 MAPQ-0-at-both. APOBEC3 sits
  with the resolvable cases. It only fires under de at the loosest D=0.02/CEIL≥0.07 corner — spurious.
- **Arrays:** annotated MAGEA/PRAMEF/XAGE members show 0 cross-member conflict — a **test-design artifact**, not
  a miss. The pipeline's real MAGEA conflicts are between **unannotated de-novo sub-loci** 24–207kb from any
  annotated gene, and those cross-map (103/69/311 reads). The port's conflict graph must be built on the
  pipeline's coordinate sub-loci, not annotated-gene vertices.

## Operating point & stability
`D_DE=0.005, CEIL_DE=0.05, MIN_READS=3`. The tight corner is the unique fully-clean point on the 9-config
D×CEIL grid: D is the load-bearing knob (D→0.01 re-admits CNN2; D→0.02 re-admits EEF1A1 catastrophically), and
0.005 sits below the smallest genuine-conflict de gap (RABL2/CCDC196 ~0.004–0.006) and below the smallest
retrocopy gap (CNN2 0.0084). Domain-sharers stay silent across the entire grid (that axis is easy).

## Port plan (into `read_conflict.rs` + `build_read_placements`)
1. `BamRead`/placement carries per record: **de** (f32, `de:f` tag), **nm** (u32, `NM:i`), **aln_len**
   (query-aligned length = read len − clips). Keep **AS** for audit only.
2. Tie predicate becomes **de-tie** (`|de_a−de_b| ≤ 0.005 && max ≤ 0.05`); NM-tie computed in parallel for the
   audit log; AS-tie kept only to assert `de-tie ⊆ AS-tie` on real data.
3. Edge iff de-tied count `≥ MIN_READS=3` (raise from 2 to clear the CNN2 NM noise-floor; de alone is 0-FP even at 2).
4. **Replace the 200 bp distinct-record guard** with a CIGAR/identity comparison before trusting this on
   genuinely collapsed tandems (<200 bp-separated copies) — the guard is panel-safe only because every same-chrom
   copy here is ≥7 kb apart. (Open item.)
5. Make `RUSTLE_CONFLICT_{D_DE,CEIL_DE,D_NM,CEIL_NM,MIN_READS}` env-overridable, defaulting to the operating point.

## Open questions
The 200 bp guard is unprincipled (untested on collapsed tandems); AK6/CCDC196 should-fire labels are soft
(near-error-floor de ties minimap2 itself resolves by MAPQ); de is minimap2's own tag (a re-align could shift
it); the 0-FP result is a 14-pair curated panel (genome-wide unverified); the de-tie behaviour on the
pipeline's de-novo sub-loci vertices (vs annotated-gene vertices) is not yet measured.

---

## Tag-discriminator dig (does any minimap2 tag beat or augment de?)
Full tag vocabulary in GGO.bam: NM ms AS nn ts tp cm s1 s2 de rl. Dumped the per-read tag vector for every
cross-mapping read on the panel (`tag_discriminator_dump.py` → `.tsv`; genuine=246, artifact_retro=3540,
artifact_asym=15; domain-sharers=0 rows, pre-killed by best-overlap attribution) and evaluated 12 single tags +
combinations + a confound audit (workflow `wf_09afe78c-599`, 5 analyses + synthesis).

**Verdict: `de` wins unchanged; nothing beats it.** Two failure modes, and `de` already covers the one that
survives attribution:

| discriminator | independent of de? | per-pair result | verdict |
|---|---|---|---|
| **de-tie** `\|Δde\|≤0.005 & max≤0.05` | — (the axis) | **7/7, 0 FP/0 FN** | **WINNER** |
| `\|Δnmr\|` (NM/alnlen) | no — `nm/alnlen ≈ de` | false-fires CNN2 | de-in-disguise |
| `Δms`, `Δcm`, `Δs1` | no — corr w/ de −0.63/−0.53/−0.28 AND length-confounded (corr w/ alnlen .58/.79/.92) | false-fire CNN2/APOBEC3 | circular + length artifact |
| `s2/s1` (chaining uniqueness) | yes | inverted on retros (~1.0), 98% empty on copy side | **dud** (blind to retrocopies) |
| `nintron` (CIGAR N-gaps) | yes (structural) | FN on genuine AK6 (copy is itself an intronless retrogene → same signature) | independent but **label-confounded** |
| **`mapq`** | **yes (placement confidence)** | both-mapq0 = genuine; retro reads mapq60-at-true | **independent corroborator** |
| `ts`, `rl`, `min_mapq`, `ms_raw` | — | AUC ≈ 0.5 | inert |

**Net value of the dig:** (a) no tag catches a case `de` misses — `de`'s tie clause already rejects APOBEC3
(all 15 reads `|Δde|` 0.013–0.038 > 0.005), the case we hoped `s2/s1` or `mapq` would need to catch;
(b) `mapq` is the one genuinely-independent, non-confounded axis — worth **logging** per edge (fraction of
supporting reads with both-side `mapq==0` = genuine-multimapper signature) and as an **optional default-off
AND-gate** (`mapq_a≤40 & mapq_b≤40`) that widens the safe-`D` window 0.005→0.014 (2.8×) for noisier substrates;
(c) the `max(de)≤0.05` ceiling is redundant on this panel (the tie clause carries the discrimination) but kept
as cheap insurance. `D=0.005` is a knife-edge against CNN2 (false-fires 29 reads at `D=0.008`) — do not loosen
without the mapq gate.

**Refined port delta:** BamRead carries **`de` + `aln_len` + `mapq`** per record (DROP `nm` — it is
de-in-disguise and added nothing); keep `AS` only to assert `de-tie ⊆ AS-tie`; exclude **supplementary**
(chimeric/split ≠ multimapping ambiguity; the shipped Rust path currently includes it). Named consts
`DE_TIE_DELTA=0.005, DE_MAX=0.05, MIN_CONFLICT_READS=3` (no buried thresholds). Per-edge report logs the
both-mapq0 fraction as a confidence field. Do NOT port `nintron`/`ms`/`cm`/`s1`/`s2`/`ts` into the gate.
