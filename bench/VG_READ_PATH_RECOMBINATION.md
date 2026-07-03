# VG-native read-path recombination — do real crossover / gene-conversion reads exist in the gorilla RNA families?

**Object.** In the O2 copy-assignment variation graph (`bench/o2_vg_visualization.py`) a multi-copy
family is one graph: **copies are paths, PSVs are bubbles**. Thread a spanning read through the ordered
bubbles → a per-bubble allele-vector. A read whose vector is consistent with **no single copy-path** but
is a clean **concatenation of two copies' paths** is a *recombinant* read — the concrete, per-read form of
the theory's **K-frontier recombination obstruction** (`bench/THEORY.md` §5 Proposition: `sep+link`
recovers copies at `K=2` but **fails at `K≥3` through cross-copy recombination**; the tight condition is
*recombination-freeness*). Two kinds, **bubble-bounded** so endpoints are exact graph objects:
`CROSSOVER` = 1 switch `A…A|B…B`; `GENE-CONVERSION` = 2 switches `A…A|B…B|A…A` (first/last switched
bubble = tract endpoints).

**Detector:** `bench/vg_read_path_recombination.py` (build + eval). Reuses the materialized family VG
(`materialize_family`), the real gate (`copy_assign.assign_read`), and the shipped RT/template-switch
discriminator (`mosaic.rs::classify_event` + `genome.rs::breakpoint_microhomology`, ported byte-faithfully
with the same constants). Outputs `vg_read_path_recombination.tsv` (per read) + `.json`
(`{summary, families, evaluation}`). Deterministic (`PYTHONHASHSEED=0`, sorted, no RNG). Substrate:
119 co-located gorilla families with ≥2 copies **and** ≥2 PSV bubbles; reads = `GGO_mm.bam` via the cache.

---

## Headline verdict (honest, after adversarial review)

**VG-native read-path recombination is REAL but VERY NARROW, and the *biological* signal is a hypothesis,
not a finding. The *operational* payoff is the robust, useful result.**

- The **raw** switch signal is dominated by **systematic near-identical-copy ambiguity**, not
  recombination: **1447/2318 rows (62%)** are whole-molecule copy-splits (recomb-fraction median **0.74**).
  **Crossovers yield ZERO credible biological recombination** (a confirmed crossover is by construction
  copy-ambiguity, not a local tract).
- The recurrent localized tracts (674 rows) **shrink again** under a **full-pattern cleanliness test**
  added for this review: only **134 reads / 134 molecules across 10 distinct loci** are clean bi-copy
  tracts. The remaining **540** are noisy A/B mosaics whose "clean 2-switch tract" is a **column-subsetting
  artifact** of unreliable asm20 copy-vectors.
- The previously-flagship **GSTM2** "textbook conversion hotspot" (fam13, 292 localized reads) **collapses
  to 1 clean read** (291 noisy, full-switch median **12**). **The GSTM2 gene-conversion claim is retracted.**
- The RT-switch microhomology veto fired **0/2318 (INERT)** and DNA is **0/2318 (unchecked)** — so
  RT-template-switch is **not** actually excluded, and the surviving tract columns are ordinary
  heterozygous sites RNA **cannot** separate from allelic SNPs without DNA parCN.
- **The robust conclusion is OPERATIONAL:** the O2 gate **force-assigns 2214/2318 (95.5%)** recombinant
  reads to a single copy — including **134/134 (100%)** of the cleanest credible tracts — when every one of
  them conflicts with **every** single copy and **should ABSTAIN**. This is a concrete per-read gate leak
  the `mosaic_discriminator` cannot see, and it is the empirical realization of the K≥3 obstruction.

---

## (1) Do real crossover / gene-conversion reads exist? — counts, named, tract lengths, post-veto

| stage | rows | distinct molecules | note |
|---|---|---|---|
| raw recombinant (61 families) | 2318 | **2227** | 91 double-counted via duplicate-locus materializations |
| — crossover / conversion | 1528 / 790 | | |
| systematic_copy_split (recFrac ~0.74) | 1447 | | **artifact class**: near-identical-copy ambiguity / unrepresented recombinant copy |
| sporadic_chimera (non-recurrent) | 197 | | per-molecule RT/PCR chimera or error |
| localized_tract (recurrent, recFrac ~0.11) | 674 | 598 | recurrent bounded tracts — *candidate* conversion |
| **↳ clean_full (full-pattern clean)** | **134** | **134** | **credible-real** clean bi-copy tracts |
| ↳ noisy A/B mosaic (subset-clean only) | 540 | | clean tract is a column-subsetting artifact |

**The full-pattern test (the key review fix).** The detector fits its "2-switch tract" on the subset of
columns where *both* copies' alleles are read-supported. Over **all** read-supported discriminating columns
(`full_switches`) the same read often carries many more switches. Flagship read
`m64076…/10225148` (GSTM2): detector = **2 switches**; full pattern = **17 switches** (`BBABBAABB…`).
The dropped columns are **asm20 copy-vector miscalls** — **96/193 (49%)** of GSTM2's copy0-vs-copy1
discriminating columns carry a copy allele with **~0 read support** (e.g. copy1='T' at a column piled
`{C:689, T:0}`). A read is credible only if it is clean over the full pattern (`clean_full=1`,
`full_switches == n_switches`). Full-switch medians: **credible = 2**, **noisy = 10**, all-localized = 7.

**Credible clean gene-conversion survivors — 134 reads / 10 distinct loci** (dedup fam10≡fam132):

| locus | copies | K | clean reads | (noisy) | tract bubbles | full-sw median |
|---|---|---|---|---|---|---|
| **THOC3** (fam73) | LOC101131460 / THOC3 | 2 | **53** | 0 | 2 | 2 |
| **RABL2** (fam39) | LOC134756389 / RABL2B | 5 | **36** | 21 | 2 / 5 / 6 | 2 |
| LOC101146782 (fam17) | LOC101146782 / LOC101141192 | 3 | 13 | 0 | 2 | 2 |
| **ANKRD18A** (fam22) | LOC101142783 / ANKRD18A | 4 | 7 | 0 | 2 | 2 |
| LOC101138917 (fam135) | LOC101138917 / LOC115932771 | 2 | 7 | 0 | 2 | 2 |
| **LARP1B** (fam123) | LOC129533124 / LARP1B | 2 | 6 | 0 | 2 | 2 |
| LOC109025132 (fam31) | LOC109025132 / LOC129523505 | 3 | 5 | 1 | 2 | 2 |
| LOC101139703 (fam47) | LOC101139703 / LOC129526951 | 2 | 4 | 0 | 2 | 2 |
| **SMG1** (fam12) | SMG1 / LOC101128221 | 5 | 2 | 3 | 3 | 2 |
| **GSTM2** (fam13) | GSTM2 / LOC134756922 | 3 | **1** | **291** | 2 | 2 (clean) / 12 (noisy) |

**Crossover survivors: 0.** Every confirmed crossover is a whole-molecule copy-split, not a bounded tract,
so the VG-native *crossover* channel yields **no** biological recombination.

**Nuance:** the credible clean tracts are **not** exclusively a K≥3 phenomenon — the largest (THOC3, 53
reads) and several others (LOC101138917, LARP1B, LOC101139703) are **K=2** families. The K-enrichment
below is a property of the *raw recombinant read-mass*, not of this small clean subset; RABL2 (K=5) is the
one clean locus with genuinely long tracts (5–6 bubbles).

**Post-veto / artifact separation is NOT delivered by the veto or by DNA.** RT-switch microhomology fired
**0/2318** (the shipped leg targets distal splice-like donor/acceptor direct repeats; it is inert for
adjacent-PSV paralog switches, only firing inside a ≥3-periodic microsatellite — none present). DNA is
**0/2318** (the DNFAM catalog uses an independent id-space; no reliable join). So **RT-template-switch is
not excluded**, and the only working filters are recurrence + recomb-fraction + full-pattern cleanliness.
A recurrent RT/PCR chimera at a homology hotspot — exactly where these reads live — would pass them.
Moreover the surviving tract columns are ordinary het sites (GSTM2 col `{C:259, T:18}`) that **RNA cannot
distinguish from a heterozygous SNP at one copy's locus** without DNA parCN. **"Gene conversion" is a
recurrent-structural-tract hypothesis, not a DNA-verified finding.**

---

## (2) K-frontier concentration — real but partly a detection-power confound

- **Enrichment (real).** K≥3 families carry recombinants **25/34** vs K==2 **36/85** (Fisher odds
  **3.78, p = 0.0024**). K≥3 are **29%** of processed families but carry **73.6%** of recombinant read-mass
  (1706/2318). This is consistent with recombination as the K≥3 obstruction: more copies → more copy-pairs →
  more spoofable combinations.
- **Confound (honest caveat).** Part of the enrichment is **mechanical detection power**, not biology:
  per-family `n_recombinant` correlates with `n_bubbles` (Pearson **r = 0.392**) and with #copy-pairs
  `C(K,2)` (**r = 0.356**). K≥3 families simply give the pair-fitting detector more bubbles and pairs to fit.
- **Read-level K-frontier (the sharp point).** The credible tracts do **not** sit on the *copy-level*
  SUN-Tier-2 (hap-unique-only, recombination-vulnerable) copies — the sole genome-wide Tier-2 family
  (fam42, LOC129529768) carries **zero** recombinant reads. Instead **127/134** credible reads bridge two
  **Tier-1 SUN-identifiable** copies (7 bridge two Tier-3 collapse copies). A single recombinant read
  carries copy A's private **SUN** in one arm and copy B's in the other, satisfying **two Strong-Separation
  witnesses at once** — impossible for any real single copy. SUN identifiability guarantees a *non*-recombinant
  read is taggable; it gives **no** protection against a recombinant read, which belongs to no single copy
  and must abstain. This is precisely the theory's cross-copy recombinant that spoofs an alternative
  minimum cover at K≥3.

---

## (3) O2 gate: does it force-assign or abstain them? — the operational payoff

The significance gate (`copy_assign.assign_read`, the real O2 engine) **force-assigns 2214/2318 (95.5%)**
recombinant reads to a single copy (104 already abstain). By structural class:

| class | total | force-assigned (leak) | abstained |
|---|---|---|---|
| systematic_copy_split | 1447 | ~1357 (93.8%) | ~90 |
| sporadic_chimera | 197 | ~183 (92.9%) | ~14 |
| **localized_tract** | 674 | **674 (100%)** | 0 |
| **↳ clean credible tracts** | **134** | **134 (100%)** | **0** |

**100% of even the cleanest credible tracts are force-assigned.** The gate matches whichever copy's private
SUN falls in the read's covered arm — exactly the recombination failure mode. These reads spoof a copy
combination and belong to no single copy; **they should ABSTAIN, not be counted toward one copy's quant.**
This payoff is **artifact-independent**: whether a switched read is real allelic conversion, an
unrepresented recombinant paralog, or a recurrent chimera, it conflicts with every single copy by ≥2
bubbles and must not be force-assigned.

---

## (4) Value over the shipped mosaic_discriminator

- **Different axis.** `mosaic.rs::classify_event` (with `--vg` `RtSwitchArtifact` suppression) gates
  isoform **emission** — which junctions to emit. It never touches copy **assignment**, so it cannot see
  that a recombinant read is being force-assigned to a copy. VG-native read-path detection lives on the
  **assignment** axis — a genuinely new axis.
- **New discriminators.** Run naked here the mosaic MH leg fires 0/2318 and would label all confirmed reads
  `GeneConversion`, **missing** the 1447 systematic copy-ambiguity reads **and** the 540 noisy mosaics.
  Bubble-bounded **recomb-fraction** and the **full-pattern switch count** are the new legs that separate
  them.
- **Exact objects.** Bubble-bounding makes the tract endpoints **exact graph objects** (first/last switched
  PSV column) and the fraction/full-switch count computable, vs the mosaic's ±`breakpoint_tol=50` bp
  donor/acceptor coords with no fraction.

---

## Is this a real, useful use of the VG for O2?

**Yes — but only as an ABSTAIN trigger, not as a gene-conversion discovery instrument.** The biological
recombination signal is too rare and too artifact-/het-ambiguous to headline: 62% of raw switches are copy
ambiguity, crossovers yield nothing, and after the full-pattern test only 134 reads on 10 loci remain — of
which even the flagship GSTM2 collapses to a single clean read, RT-switch is not excluded, and the tract
columns are unresolvable from allelic SNPs by RNA alone. What **is** robust and useful is the operational
finding: bubble-bounded VG-native threading identifies a per-read population (2214/2318 force-assigned;
100% of even the cleanest tracts) that the current O2 gate mis-assigns and that the `mosaic_discriminator`
cannot catch. These reads are the concrete data realization of the theory's **K≥3 cross-copy recombination
obstruction** (recombination-freeness fails at the read level), and the correct O2 behavior is to **flag and
abstain** them. That is the deliverable: a VG-native abstain lever for the copy-assignment gate.

---

## (5) SHIPPED — the recombinant-abstain gate leg is WIRED (default-on)

The abstain lever is now **wired into the O2 copy-assignment gate** as a **default-on** leg with an
opt-out, mirroring `--no-repeat-gate` / `--no-split-recombinants`.

- **Assign path wired:** `o2_vg_visualization.materialize_family`, which threads every read through the
  significance gate `copy_assign.assign_read` (the min_p / margin certificate) and labels the per-read
  `status`. That per-read `status` is exactly what force-assigned the 2214. The abstain leg runs on the
  materialized VG (`recombinant_abstain.apply_abstain_to_vg(vg)`) at both `materialize_family` return
  points. `assign_read` itself is **unchanged** (byte-identity for `validate_sim5x` / `crossval` /
  `assign_family`), and the o2vg **cache stays the pure-gate assignment** — the leg is applied on read,
  never baked into the JSON. All Python; the 2214 empirics run through this Python path (no Rust path).
- **Detection logic reused, not re-derived:** `bench/recombinant_abstain.py` is the single source of
  `detect_read_recombination` + `full_pattern_switches` (moved verbatim out of the detector, which now
  IMPORTS them) and adds the reusable `is_recombinant(obs, supported, copy_vecs, names) -> (bool, kind,
  arms)` predicate + `apply_abstain_to_vg`. The diagnostic detector runs with the leg opted OUT so it
  always characterizes the raw pure-gate leak.
- **Opt-out (default-on):** `--no-recombinant-abstain` (sets `RUSTLE_NO_RECOMBINANT_ABSTAIN=1`) makes
  `apply_abstain_to_vg` a no-op; the assignment output is then byte-identical to the prior pure gate.

**Reads moved assigned → abstain (default-on): 2214, by structural class**

| class | moved (assigned→abstain) | of which clean_full=1 | of which noisy (clean_full=0) |
|---|---|---|---|
| systematic_copy_split | 1357 | 425 | 932 |
| localized_tract | 674 | **134** | 540 |
| sporadic_chimera | 183 | 62 | 121 |
| **total** | **2214** | **621** | **1593** |

- The **674 localized_tract** reads (incl **GSTM2 fam13 = 292 localized**, `1 clean / 291 noisy`) all
  move to abstain. **Honest label:** these are **belongs-to-no-copy reads carrying a subset-clean bi-copy
  recombinant pattern** — NOT "credible gene-conversion". Only **134/674** localized (and **621/2214**
  overall) are clean over the full read-supported pattern; the **GSTM2 gene-conversion claim is
  retracted** (§1 above, same commit).
- **The leg deliberately does NOT gate on `clean_full`.** Abstaining is correct for the full 2214
  because **every** abstained read — clean or noisy — conflicts with every single copy by **≥2**
  read-supported bubbles (`min_single_mm ≥ 2`) = belongs to no single copy. Gating on `clean_full` would
  wrongly **re-admit** the 540 noisy belongs-to-no-copy reads to force-assignment.

**BYTE-IDENTITY on non-recombinant reads: PASS.** Independent leg-OFF vs leg-ON rerun over the cache
(**102,224** reads across the flippable families): **2214 flipped**, **0 non-recombinant reads changed**,
**0 monotonicity violations**. Every flip is exactly `assigned → recombinant` with `best_copy=None`;
every non-recombinant read (and every already-abstaining read) is bit-for-bit unchanged. It is a **pure
monotone addition of abstentions** — no copy-swap, no abstain→assign, ever. Opt-out recovers the
pure-gate baseline byte-for-byte (0 violations). Deterministic under `PYTHONHASHSEED=0` (two leg-ON
passes, 0 mismatches).

**Coverage-vs-correctness trade (honest).** Across the flippable families (102,224 reads) assigned
**40,534 → 38,320** (39.65% → 37.49%). Coverage cost = **2214 abstentions = 5.46% of previously-assigned
reads (2.17 pts)**. Correctness gain = **2214 belongs-to-no-copy force-assignments removed**, incl 100%
(674/674) of the localized tracts. RNA cannot separate allelic gene-conversion from an unrepresented
recombinant paralog (needs DNA parCN), but **both classes abstain**, so the leg is correct regardless.

**Tests:** `bench/test_recombinant_abstain.py` — 5 pass: (1) predicate flags crossover not clean
single-copy; (2) monotone (only the force-assigned recombinant read flips; all others byte-identical);
(3) opt-out no-op; (4) determinism; (5) data-backed 2214 move (674 localized incl GSTM2, 1357 systematic,
183 sporadic; 0 byte-identity violations).

---

### Reproduce
```
PYTHONHASHSEED=0 /home/juanfra/miniforge3/bin/python bench/vg_read_path_recombination.py        # build+eval
PYTHONHASHSEED=0 /home/juanfra/miniforge3/bin/python bench/vg_read_path_recombination.py eval    # eval only
```
Outputs: `bench/vg_read_path_recombination.tsv` (per-read; `full_switches`, `clean_full`, `recomb_fraction`,
`structural_call`, `gate_status`, `gate_leak`, `should_abstain`), `bench/vg_read_path_recombination.json`
(`summary` / `families` / `evaluation`).

### Review fixes applied to the detector
1. **`full_switches` + `clean_full` columns** (per read): switch count over all read-supported discriminating
   columns; credible-real is restricted to `clean_full=1` (674→134). Exposes and demotes the column-subsetting
   artifact (GSTM2 292→1).
2. **Distinct-molecule dedup** of duplicate-locus families (fam10≡fam132 share all reads): 2318→2227 total,
   674→598 localized, 17→16→(clean) 10 loci reported as loci not families.
3. **RT-switch veto reframed as INERT** (fired 0/2318): the report no longer claims the veto separates real
   from artifact; recurrence + fraction + full-pattern do, and RT-template-switch is explicitly *not* excluded.
4. **Allele-vs-conversion disclaimer**: tract columns are het sites RNA cannot separate from allelic SNPs
   without DNA parCN.
5. **Detection-power confound** reported: `n_recombinant` ~ `n_bubbles` r=0.392, ~`C(K,2)` r=0.356.
6. **`flag_filter=0xF04`** in `psv_graph_genomewide.column_alleles` to exclude supplementary (0x800) from the
   pileup (source-corrected; effect ≤1 read on the cached run, immaterial).
