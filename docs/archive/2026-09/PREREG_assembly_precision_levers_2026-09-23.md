# Pre-registration — two precision levers for `--assemble-only`, and whether they become the defaults

**Written 2026-09-23 (§6za), before any arm is scored.** User: *"yes lets try those, and if they look good lets
change the assembler defaults."* Both levers come from r1069/r1070 (the tool-strata table and the simulation).

## Lever 1 — strict junctions for the transcript product

The §6m8 majority mode (default since 2026-09-21) tolerates a minority of non-canonical junctions ≤ 10 kb when
a canonical majority fixes the strand. It was adopted for FAMILY recovery (NPIPB12: 9 junctions, 8 canonical).
r1069 measured its cost on the transcript product: **1,391 chains (9.7% of human chr20-22 output) at 1.2%
annotation-exact; 217 (5.1%) at 0.0% on gorilla NC_073244.2** — against 16.0% / 36.6% for the rest.

Arm: the shipped polished assembly with `RUSTLE_JUNCTION_MAJORITY=0` (strict: every junction canonical, one
strand). If adopted, it ships as `copy_assign --assembly-junctions strict|majority` with `strict` the default
**under `--assemble-only` only**; every family path keeps majority (their catalog numbers stay byte-identical).

## Lever 2 — retained-intron chains spanned by a dominant junction

r1070 (simulation): 3.75% of transcripts yield a wrong chain carried by ≥ 2 reads, 71% of them a "missing
junction" — the aligner reading through a short exon, consistently across reads of one molecule. In the real
output such a chain is a transcript T with an exon that fully contains the intron of a junction J carried by
another transcript U at the same locus. Rule (polish stage, reads-only, annotation-free, hence legal de novo):

> drop T if some junction J of another transcript at T's locus (same `gene_id`, same strand) lies strictly
> inside one of T's exons and `support(J) ≥ R × reads(T)`, where `support(J)` = reads of all emitted
> transcripts at the locus that carry J.

R is swept on human over {1, 2, 5, 10}; the selection rule is the largest R whose human intron-chain precision
gain is ≥ 0.5 pt (a looser R drops less), else the R with the best precision gain at ≤ 1% chain loss. One R
goes to gorilla. Emulated in Python on the polished GTF first; if it passes, ported into the Rust polish
(`--polish-retained-ratio R`, default R if adopted) with byte-parity against the emulation.

## Substrates and metrics

- **Development: human A119b chr20/21/22**, RefSeq CHM13, gffcompare v0.12.10 on the polished output
  (baseline 15,876 transcripts, intron chain 23.3 / 16.0, 2,295 matching chains).
- **Held out: gorilla `GGO_mm.bam`, all 26 contigs** (baseline 78,602, 27.0 / 33.1, 25,829 chains).
- Also reported: the annotation-exact rate of the DROPPED set vs the kept set (is the lever removing junk?),
  the 15-cell tool comparison after the change, and the tool-consensus share (r1065) after the change.

## The bar — judged on held-out gorilla, each lever alone, then both together

| outcome | verdict |
|---|---|
| intron-chain precision **+ ≥ 1.0 pt** and matching chains **− ≤ 1.0%** | ⭐ **ADOPT as default** |
| precision + ≥ 1.0 pt and chains − 1.0 to 3.0% | ⚠ **TRADE** — ship as a flag, user's call on the default |
| precision + < 1.0 pt, or chains − > 3.0% | ⛔ **NO** |

**Predicted, before looking:** Lever 1 ⭐ — the tolerated chains are 5-10% of the output at ~0-1% precision, so
dropping them should add ~1.5-2 pts of precision for a loss of a few chains (the 1-2% of them that are
annotated); gorilla's 217 at 0.0% makes it near-free there. Lever 2 ⚠/⭐ — the retained-intron class is smaller
(the taxonomy's `retained` was 3.3% of misses, and genuine retained-intron isoforms exist in RefSeq), so a
precision gain under 1 pt at R = 5-10 is likely; the risk is dropping real minor isoforms at low R.

I will not change the rules, the substrates or the bar after seeing any number.

---

# OUTCOME (2026-09-23)

## Lever 1 — strict junctions (`RUSTLE_JUNCTION_MAJORITY=0` on the shipped polished assembly)

| substrate | arm | transcripts | intron SN / PR | **intron chain SN / PR** | **matching chains** | Δ precision | Δ chains |
|---|---|---|---|---|---|---|---|
| human chr20-22 (dev) | majority (shipped) | 15,876 | 61.0 / 56.3 | 23.3 / 16.0 | 2,295 | | |
| | **strict** | 14,349 | 60.9 / **62.1** | 23.2 / **17.5** | 2,280 | **+1.5** | −15 (−0.65%) |
| gorilla, 26 contigs (held out) | majority (shipped) | 78,602 | 59.6 / 81.7 | 27.0 / 33.1 | 25,829 | | |
| | **strict** | 75,686 | 59.3 / **82.7** | 27.0 / **34.3** | 25,797 | **+1.2** | −32 (−0.12%) |

⭐ **ADOPT** on held-out: precision + 1.2 pt, chains − 0.12%. Intron-level precision moves the most
(+5.8 / +1.0), which is the class of junction the tolerance was admitting.

## Lever 2 — retained-intron chains spanned by a dominant junction (emulated, then ported)

Human sweep of R (dropped-set annotation-exact rate in brackets): R = 1 → 2,159 chains at 18.6 (4.9%);
R = 2 → 2,217 at 18.1 (3.7%); R = 5 → 2,260 at 17.5 (2.4%); **R = 10 → 2,273 at 17.1 (2.0%)**. Selection
rule (largest R with ≥ 0.5 pt gain) picks **R = 10**: human precision +1.1, chains −22 (−0.96%).

| substrate | arm | transcripts | intron chain SN / PR | matching chains | Δ precision | Δ chains | dropped set annotated / kept set annotated |
|---|---|---|---|---|---|---|---|
| gorilla (held out) | shipped | 78,602 | 27.0 / 33.1 | 25,829 | | | |
| | **R = 10** | 75,427 | 26.9 / **34.3** | 25,721 | **+1.2** | −108 (−0.42%) | **3.4% vs 34.2%** |
| | R = 5 (info) | 74,606 | 26.8 / 34.6 | 25,628 | +1.5 | −0.78% | 5.0% vs 34.5% |

⭐ **ADOPT** on held-out at R = 10: precision + 1.2 pt, chains − 0.42%, and the dropped set is ten times less
often annotated than what stays. Ported to the Rust polish as `--polish-retained-ratio` (one pass over the
survivors of the other steps, support and candidate junction sets fixed before any drop).

## Both together, parity, and the default flip

Rust port parity: the polished transcript sets under `--polish-retained-ratio 10` are identical to the
emulation's on human chr20-22 (14,792 = 14,792) and gorilla genome-wide (75,427 = 75,427).

| substrate | arm | transcripts | intron SN / PR | **intron chain SN / PR** | **matching chains** | Δ precision | Δ chains |
|---|---|---|---|---|---|---|---|
| human chr20-22 (dev) | shipped (09-22) | 15,876 | 61.0 / 56.3 | 23.3 / 16.0 | 2,295 | | |
| | **strict + retained 10** | 13,386 | 60.7 / 62.6 | 23.0 / **18.7** | 2,259 | **+2.7** | −36 (−1.6%) |
| gorilla, 26 contigs (held out) | shipped (09-22) | 78,602 | 59.6 / 81.7 | 27.0 / 33.1 | 25,829 | | |
| | **strict + retained 10** | 72,690 | 59.3 / 82.8 | 26.9 / **35.6** | 25,710 | **+2.5** | −119 (−0.46%) |

⭐ **ADOPT (held out: + 2.5 pt, − 0.46%)** — the two levers are additive (1.2 + 1.2 ≈ 2.5). 15-cell tool
comparison: gorilla **10 → 12/15** (StringTie now 5/5, its precision 34.6 vs our 35.6); human chr20-22
12 → 10/15 (two FLAIR cells lost by a hair: chain sensitivity 23.0 vs 23.1, matching chains 2,259 vs
2,268 — the dev-side cost of −36 chains). Tool consensus (r1065) is unaffected in kind: the dropped
chains are 1-3% annotated.

**Defaults changed (user, 2026-09-23):** under `--assemble-only`, `--assembly-junctions strict` and
`--polish-retained-ratio 10` are the defaults; every family path (`gw_family_catalog`, `mcl_families`,
`copy_assign` without `--assemble-only`) is untouched — the junction switch only runs under `--assemble-only`
and an explicit `RUSTLE_JUNCTION_MAJORITY` always wins. Verified on full chr21: a default run is
**byte-identical** to the combined arm, and `--assembly-junctions majority --polish-retained-ratio 0` is
**byte-identical** to the 2026-09-22 output. `bench/assembly_polish.py` mirrors the new step
(`--retained-ratio`). `REPRODUCE.md` and `tools/genome_wide_sweep.sh` carry both recipes.

---

# ADDENDUM (2026-09-23) — the partial matches in OUR output: is a "full transcript" rule missing?

Every emitted multi-exon transcript under the new defaults, by gffcompare class (`partials.py`):

| class | human chr20-22 (12,076) | gorilla genome-wide (72,271) |
|---|---|---|
| `=` exact | 18.7% | 35.6% |
| `j` shares junctions, different combination | 51.7% | 47.3% |
| `c` contained in the reference (fragment) | 5.9% | 8.3% |
| `k` contains the reference (extra junctions) | 2.9% | 3.2% |
| `m` / `n` retained intron | 2.3% / 3.6% | 1.4% / 0.9% |
| `u` / `i` / `x` / `o` (unannotated, intronic, antisense, other overlap) | 6.1 / 3.6 / 2.2 / 3.0% | 1.4 / 0.3 / 1.0 / 0.5% |

**Fragments (`c`)** are 5′-truncated in 82% / 79%; median 8 / 5 reads. Reads carrying the FULL reference chain:
0 in 61% / 62%, 1 in 10% / 12%, ≥ 2 in 29% / 26%. Attribution on human (721 fragments): **75% have no buildable
full chain (fewer than 2 distinct-coordinate reads, r1067/r1058 territory); 23% are emitted alongside the
full chain** (the ISM ratio kept them as a shorter isoform with ≥ 70% of the container's support); only **2%**
lost the full chain to the fraction rule and 1% to another polish step — where it happens, the fragment has
median 18 reads and the full chain 4. So the polish is not the reason we emit fragments.

**The obvious rule was tested and is worse**: replace an emitted fragment by its best-supported raw container
(a chain of ≥ 2 or ≥ 3 reads at the same locus that contains it, extending at the 5′ side only or at any
side), human chr20-22, new defaults as baseline:

| arm | intron chain SN / PR | matching chains |
|---|---|---|
| new defaults | 23.0 / 18.7 | 2,259 |
| extend 5′ to a ≥ 2-read container | 22.5 / 18.3 | 2,209 |
| extend 5′ to a ≥ 3-read container | 22.8 / 18.6 | 2,245 |
| extend any side, ≥ 2 reads | 21.4 / 17.4 | 2,107 |
| extend any side, ≥ 3 reads | 22.4 / 18.2 | 2,198 |

⛔ Every variant loses chains and precision: the fragment we emit is the annotated form more often than the
longer read-supported chain that contains it — r1072's finding again, now with ≥ 2-read containers.

**Over-extensions (`k`)**: extra junctions at the 5′ side in 59% / 72%, median 4 reads, and the shorter
reference chain is itself emitted in only 15.5% / 9.0% of cases (it rarely has 2 exact reads of its own).
Without annotation there is no signal to trim by; the reads carry the extension. **`j`** is half of the output
on both species, with the closest reference chain also emitted in 50% / 65% of cases — alternative junction
combinations at annotated genes, which the annotation cannot adjudicate.

**Conclusion: no further "full transcript" rule is missing at the chain level.** What is partial in our output
is partial because the full chain is not carried by two reads, or because a shorter isoform was deliberately
kept beside the full one; extending fragments makes the output worse.

# ADDENDUM 2 (2026-09-23, written before running) — trim 5′ over-extensions by junction support?

`k` transcripts (2.9% / 3.2% of output) extend the annotated chain, 59-72% at the 5′ side, on a median of 4
reads. Their own reads carry the extension, so the only annotation-free signal is RELATIVE junction support:
the extension's junction is carried by few reads at the locus while the interior junctions are carried by many.

**Rule (one parameter f):** for every emitted multi-exon transcript with ≥ 2 introns, let e = reads carrying
its 5′-most junction (any read at the locus) and m = median reads carrying its other junctions. If e < f × m,
remove the 5′-terminal exon and intron (one trim). If the trimmed chain is already emitted, the transcript
is dropped instead (merged). Variants: 5′ side only (pre-registered decision arm) and both ends (info).
f ∈ {0.05, 0.10, 0.20} swept on human chr20-22; selection = the f with the largest intron-chain precision gain
at ≤ 1% chain loss; that f goes to gorilla genome-wide. Bar as the levers': ⭐ precision + ≥ 1.0 pt and chains
− ≤ 1%; ⚠ + ≥ 1.0 pt at − 1-3%; ⛔ otherwise.

**Prediction: ⛔.** The extension's own reads are real evidence; low relative support mostly means a minor
alternative-TSS isoform, and trimming manufactures chains carried by no read (isoseq's 0-read class sits at
0% precision, r1069). Expect chains to fall faster than precision rises.

## ADDENDUM 2 — OUTCOME (2026-09-23): ⛔ refuted on development; gorilla not spent

Human chr20-22, new defaults as baseline (`trim_emul.py`, junction support from all primary reads):

| arm | transcripts | intron SN / PR | intron chain SN / PR | matching chains |
|---|---|---|---|---|
| new defaults | 13,386 | 60.7 / 62.6 | **23.0 / 18.7** | **2,259** |
| trim 5′, f = 0.05 (931 trimmed, 544 merged) | 12,842 | 59.2 / 65.1 | 20.8 / 17.8 | 2,049 |
| trim 5′, f = 0.10 | 12,576 | 58.4 / 65.6 | 19.5 / 17.0 | 1,914 |
| trim 5′, f = 0.20 | 12,267 | 57.6 / 66.0 | 17.9 / 16.1 | 1,764 |
| trim both ends, f = 0.05 | 12,636 | 58.8 / 66.4 | 20.4 / 17.7 | 2,003 |
| trim both ends, f = 0.20 | 11,869 | 56.8 / 68.1 | 16.7 / 15.5 | 1,640 |

Every setting loses matching chains (−9% at the gentlest) AND chain-level precision, while intron-level
precision rises — the signature of manufacturing chains no read carries: the trimmed transcript keeps its
good interior junctions (intron precision up) but is no longer any annotated chain (chain precision down).
The low-support terminal junction is, more often than not, a real minor 5′ isoform. No point reaches the
bar on development, so the held-out substrate was not spent. Prediction confirmed.

# ADDENDUM 3 (2026-09-23, written before running) — restricted trims, stronger retained-intron rules, unlikely exon combinations

User: *"what if we tried that rule only on certain cases? Also should there be rules to avoid retained introns?
or can we model unlikely exon combinations?"* Three families, each one or two parameters, emulated on human
chr20-22 (new defaults as baseline), the best of each family taken to gorilla genome-wide only if it reaches the
levers' bar on development (⭐ precision + ≥ 1.0 pt at chains − ≤ 1%; ⚠ + ≥ 1.0 at − 1-3%; ⛔ otherwise).

**T — restricted 5′ trims** (the r1078 trim only where it cannot manufacture a new chain):
- T1 *merge-only*: trim the 5′ terminal exon+intron only when the trimmed chain is ALREADY an emitted
  transcript (the extended form is then dropped as a 5′-extended duplicate), with the r1078 support condition
  at f ∈ {0.05, 0.10, 1.0} (1.0 = no support condition, pure merge).
- T2 *unique-extension*: as T1, but the condition is that the terminal junction is carried by no read other
  than the transcript's own (support == own reads).

**RI — retained introns beyond the adopted ratio-10 filter:** R ∈ {5, 7} (finer than the adopted sweep), and
RI-n: also drop a transcript whose exon overlaps another transcript's junction PARTIALLY (the `n` class:
exon end inside the intron) when that junction's support ≥ R × its reads.

**S — unlikely exon combinations (exon skipping the aligner can manufacture, r1070's "other" class):** drop T
if an intron of T strictly contains an exon of another emitted transcript at the locus of length ≤ L, and both
of that exon's flanking junctions carry ≥ R × reads(T). L ∈ {50, 100, 200}, R ∈ {5, 10}.

**Prediction:** T1 at f = 1.0 (pure merge of 5′-extended duplicates into the shorter emitted chain) ⚠ at best —
it removes only transcripts whose shorter form already exists, so chains cannot fall, but r1072 says the
extension is sometimes the annotated one, so precision moves little; T2 ⛔ (own-read-only junctions are
the minor-isoform signature). RI at R = 5/7: a trade below the bar (R = 5 was + 0.3 pt for − 0.4% more chains
on gorilla); RI-n ⛔ (partial overlaps are alternative splice sites more than artefacts). S ⛔/⚠: microexon
skipping is real biology as often as an alignment miss; expect < 0.5 pt at any L.

## ADDENDUM 3 — OUTCOME (2026-09-23): ⛔ all three families on development; gorilla not spent

Human chr20-22, new defaults 13,386 transcripts / intron chain 23.0 / **18.7** / **2,259** chains (`rules3_emul.py`):

| family / setting | dropped | intron chain SN / PR | matching chains | Δ precision | Δ chains |
|---|---|---|---|---|---|
| T1 merge-only trim, no support condition | 1,083 | 21.3 / 19.0 | 2,091 | +0.3 | **−7.4%** |
| T1 merge-only, terminal junction < 10% of the others | 337 | 22.4 / 18.8 | 2,204 | +0.1 | −2.4% |
| T1 merge-only, < 5% | 249 | 22.6 / 18.8 | 2,224 | +0.1 | −1.5% |
| T2 unique-extension merge (terminal junction carried only by its own reads) | 104 | 22.7 / 18.7 | 2,235 | 0.0 | −1.1% |
| RI-n partial retained intron, R = 10 | 3,863 | 18.5 / 22.1 | 1,819 | +3.4 | **−19.5%** |
| RI-n, R = 5 | 5,053 | 16.6 / 23.2 | 1,628 | +4.5 | **−27.9%** |
| S exon-skip, exon ≤ 50 bp, flanks ≥ 5× | 235 | 22.6 / 18.8 | 2,227 | +0.1 | −1.4% |
| S, ≤ 100 bp, ≥ 5× | 1,226 | 21.0 / 19.1 | 2,069 | +0.4 | −8.4% |
| S, ≤ 100 bp, ≥ 10× | 884 | 21.7 / 19.1 | 2,137 | +0.4 | −5.4% |
| S, ≤ 200 bp, ≥ 5× | 2,203 | 19.6 / 19.5 | 1,927 | +0.8 | −14.7% |

- **Restricted trims (T):** even the pure merge of a 5′-extended chain into its already-emitted shorter form
  costs 7.4% of matching chains for +0.3 pt — the extended form is the annotated one about as often as the
  shorter one when both are read-supported (r1072's 2.5:1 held only for single-read extensions). The
  support-conditioned variants do almost nothing.
- **Partial retained introns (RI-n):** the class is mostly alternative splice sites, not artefacts — dropping
  it buys +3-4 pts of precision at −20-28% chains, far past the trade band. The adopted strict-containment
  rule (r1075) is the right boundary.
- **Exon skipping (S):** at ≤ 50 bp (true microexons) the rule finds 235 transcripts and nothing changes; wider
  windows remove real skipping isoforms 10× faster than they raise precision. The aligner's "other" class in
  the simulation (1.9% of reads) does not concentrate in emitted ≥ 2-read chains.

⛔ No setting reaches + 1.0 pt at ≤ 3% chain loss. Predictions held (T1 was predicted ⚠, measured ⛔).

# ADDENDUM 4 (2026-09-23) — StringTie's rule: one full-length read per chain, junctions covered ≥ J

Emulated as: add every single-read all-canonical chain not already emitted whose junctions each carry ≥ J
primary reads (StringTie `-L`: per-bp coverage ≥ 1 for the chain, `-j` junction coverage; J = 2 is the
"1.5 reads per junction" reading), on top of the new-default output (`stringtie_rule.py`):

| substrate | arm | transcripts | intron SN / PR | **intron chain SN / PR** | **matching chains** |
|---|---|---|---|---|---|
| human chr20-22 | new defaults | 13,386 | 60.7 / 62.6 | 23.0 / **18.7** | 2,259 |
| | + 1-read chains, junctions ≥ 2 (35,199 added) | 48,585 | 69.9 / 45.6 | 27.9 / **5.8** | 2,742 |
| | junctions ≥ 3 (30,184) | 43,570 | 67.5 / 51.3 | 26.8 / 6.2 | 2,640 |
| | junctions ≥ 5 (24,923) | 38,309 | 64.8 / 56.8 | 25.9 / 6.9 | 2,544 |
| gorilla NC_073244.2 (held out) | new defaults | 3,862 | 54.8 / 83.3 | 28.3 / **40.8** | 1,571 |
| | + 1-read chains, junctions ≥ 2 (6,314) | 10,176 | 61.7 / 75.7 | 34.1 / **18.7** | 1,897 |
| | junctions ≥ 3 (5,037) | 8,899 | 58.4 / 80.1 | 32.9 / 20.5 | 1,825 |
| | junctions ≥ 5 (3,985) | 7,847 | 56.0 / 82.3 | 31.0 / 22.0 | 1,723 |

⛔ **It is StringTie's operating point, not a rule that holds ours:** +21% matching chains for a 3.2× (human)
/ 2.2× (gorilla) drop in intron-chain precision, and the added chains are annotation-exact at ~1.4% / ~5%
(483 true of 35,199; 326 of 6,314). Raising the junction floor to 5 removes a third of the additions and
recovers a fifth of the precision. Under the 15-cell comparison every precision cell against StringTie and
FLAIR would flip to a loss. Consistent with r1067 (the same predicate per chain), r1069 (StringTie's own
1-read stratum at 4.8% / 14.3%) and r1073 (no learned cut lifts this class).

# ADDENDUM 5 (2026-09-23) — StringTie's isoform fraction (0.01) in two places

**(a) as the gate for one-read chains** (Addendum 4's additions kept only where 1 read ≥ F × the locus maximum):

| substrate | arm | transcripts | intron chain SN / PR | matching chains | F |
|---|---|---|---|---|---|
| human chr20-22 | new defaults | 13,386 | 23.0 / 18.7 | 2,259 | 20.6 |
| | + one-read chains, junctions ≥ 2 | 48,585 | 27.9 / 5.8 | 2,742 | 9.6 |
| | … kept only if ≥ 1% of locus (22,274 kept) | 35,660 | 26.8 / 7.7 | 2,637 | 12.0 |
| | … ≥ 2% (15,841 kept) | 29,227 | 25.9 / 9.1 | 2,549 | 13.5 |
| gorilla NC_073244.2 | new defaults | 3,862 | 28.3 / 40.8 | 1,571 | 33.4 |
| | + one-read chains, junctions ≥ 2 | 10,176 | 34.1 / 18.7 | 1,897 | 24.2 |
| | … ≥ 1% (4,945 kept) | 8,807 | 33.7 / 21.3 | 1,871 | 26.1 |
| | … ≥ 2% (3,755 kept) | 7,617 | 32.9 / 24.0 | 1,825 | 27.8 |

The fraction removes the singletons at deep loci (a third to half of them) but the survivors are still ~2-6%
annotated; F stays 7-8 points below the defaults. ⛔.

**(b) as our polish's isoform fraction** (`--polish-isoform-fraction 0.01` in place of 0.02, everything else the
new defaults):

| substrate | fraction | transcripts | intron chain SN / PR | matching chains | F |
|---|---|---|---|---|---|
| human chr20-22 | 0.02 (default) | 13,386 | 23.0 / 18.7 | 2,259 | 20.6 |
| | **0.01** | 15,048 | 24.0 / 17.2 | 2,358 (+4.4%) | 20.0 |
| gorilla, 26 contigs (held out) | 0.02 (default) | 72,690 | 26.9 / 35.6 | 25,710 | 30.6 |
| | **0.01** | 81,743 | 27.7 / 32.5 | 26,432 (+2.8%) | 29.9 |

A genuine recall/precision dial with a far better exchange rate than any singleton rule (+2.8-4.4% chains for
−1.5 / −3.1 pts), but net negative by F on both substrates and the wrong direction for the levers' bar. The
0.02 default stands (it was validated held-out in §6p9 and again here); 0.01 is available as the flag for a
recall-leaning run.
