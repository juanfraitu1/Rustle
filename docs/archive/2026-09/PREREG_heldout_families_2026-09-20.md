# Pre-registration — HELD-OUT FAMILIES: does the shipped O1 rule generalise, or is it NPIP-shaped?

**Written 2026-09-20, BEFORE running anything on the held-out chromosomes.** User: *"test O1, O2 and
frequently asked questions from my advisor to test if everything holds and test more families to prove it
is not a fluke or overfitting."*

This answers `docs/ADVISOR_QUESTIONS.md` **§1.3 "Everything you have is one family"** and **Q6 "Does it
port to other families?"**, whose current answer concedes human n = 3 (NPIP, TBC1D3, AMY) and gorilla
n = 1. Nothing is tuned here. The rule below is the one already shipped; if it needs adjusting to pass,
it has failed.

## 1. What is held out, and why these chromosomes

Counting `\bchrN\b` mentions in `docs/o1_ledger.md` — the record of every O1 decision ever made:

| chromosome | ledger mentions | status |
|---|---|---|
| chr16 | 115 | **development** (NPIP) — every early O1 decision was scored here |
| chr5 / chr20 / chr15 / chr1 | 51 / 46 / 42 / 35 | exposed (AMY on chr1; assembly panel) |
| chr11 / chrY / chr14 / chr21 / chr7 | 25 / 21 / 16 / 15 / 13 | exposed |
| chr22 / chr19 / chr17 / chrX / chr9 | 8 / 7 / 6 / 5 / 4 | lightly exposed (TBC1D3 on chr17) |
| ⭐ **chr2** | **0** | **HELD OUT — never mentioned in the ledger** |
| ⭐ **chr6** | **0** | **HELD OUT — never mentioned in the ledger** |

**chr2 and chr6 are the primary test.** They are the only two chromosomes with zero ledger exposure, so
no parameter, threshold or exception in the shipped rule can have been chosen by looking at them. Every
other chromosome is reported as a secondary, weaker arm.

## 2. The rule — frozen, identical to the chr16 run

Exactly the guided recipe re-verified today in `docs/SLIM_REGRESSION_2026-09-20.md`:

```sh
# 1. gene + pseudogene bodies from RefSeq CHM13 (never HSA_genomic.gff — it drops 29.1% of loci)
zcat Reference/chm13v2.0_RefSeq_full.gff.gz \
  | awk -F'\t' '$1=="chrN" && ($3=="gene"||$3=="pseudogene"){print $1":"$4"-"$5}' | sort -u > chrN.regions
samtools faidx chm13v2.0.fa -r chrN.regions > chrN.bodies.fa      # keep faidx's native names
# 2. all-vs-all
minimap2 -x asm20 -c --eqx -P -t 4 chrN.bodies.fa chrN.bodies.fa > chrN.paf
# 3. THE SHIPPED RULE — no per-chromosome flags, no exceptions
mcl_families --paf chrN.paf --gff chm13.gff --min-exonic-bp 1 --min-shared-exon-frac 0.60 --out chrN
```

⚠ `--min-shared-exon-frac` is inert without `--min-exonic-bp 1`. ⚠ `mcl_families` keys the PAF on
`chrom:start-end` (`parse_gene_key`); a PAF named by gene symbol silently yields 0 nodes.

## 3. The truth — committed now, and its known weakness

**A truth family is a set of ≥ 3 RefSeq genes on that chromosome whose `Name=` shares a symbol root.**
The root is the symbol with a trailing copy-suffix stripped by ONE fixed regex, applied blind:

```
root = re.sub(r'(?:P\d+|\d+|[A-Z])$', '', symbol)      applied ONCE, not repeatedly
```

> ⚠ **AMENDMENT, 2026-09-20, before any held-out chromosome was scored** (chr2's all-vs-all was still
> running; no cluster file existed yet, and none had been opened). The regex as first committed said
> *"applied repeatedly until stable"*, which **over-strips**: `TBC1D3K → TBC1D3 → TBC1D → TBC1 → TBC`
> and `GOLGA6L2 → GOLGA6L → GOLGA6 → GOLGA → GOLG`. That would have merged genuinely distinct families
> into one truth set and made the truth, not the method, the thing under test. Changed to a **single
> pass**. Disclosed here rather than silently: the amendment is to the TRUTH CONSTRUCTION, made before
> any score was computed, and it makes the test **harder**, not easier — a single pass splits NPIP into
> NPIPA and NPIPB, so the method must now recover finer families than "one NPIP blob".

- Genes whose symbol begins `LOC` are **excluded from truth** (unnamed loci carry no nomenclature claim)
  but are **left in the input**, so they can still cost us precision. This is deliberate: it is the
  harder direction.
- Roots shorter than 3 characters are dropped (avoids collapsing unrelated symbols).
- Families with < 3 members are dropped (a pair is not evidence of a multi-copy family).

⚠ **Stated weakness, before the result.** Nomenclature is not a phylogeny. Register rows already record
that a *name-based* subgroup truth for TBC1D3 was **wrong and cost 16 retracted claims** — that failure
was about **subgroups within a family** (AE/CDKL), not about family membership, which is what is scored
here. Even so, symbol roots will (a) split genuine families that were named inconsistently and (b) merge
paralogous-but-distinct groups. Both directions are recorded as failures below, not explained away.

## 4. Metrics — sensitivity, precision, bipartite F (the user's standing reporting rule)

For each truth family, and pooled over all of them:

- **bipartite matching** between predicted clusters and truth families, one-to-one, maximising total
  overlap (`scipy.optimize.linear_sum_assignment`);
- **sensitivity** = matched members / truth members; **precision** = matched members / predicted members
  of the matched cluster; **F** = harmonic mean.
- Reported **per chromosome** and **pooled**, and additionally as **families fully recovered**
  (F = 1.000) out of all truth families.

## 5. The bar — committed BEFORE the run

The honest overall development-set number on record is **F = 0.628 across all 72 families**
(`project_overall_family_metrics`; the widely quoted 0.881 is 11 hand-picked families and is NOT the
comparison). The ground-truth ceiling against a single annotation is ~0.8 (§6kl/§6km).

| outcome on chr2 + chr6 pooled | verdict |
|---|---|
| **F ≥ 0.55** | ⭐ **HOLDS** — within ~0.08 of the development-set 0.628; the rule ports |
| 0.40 ≤ F < 0.55 | ⚠ **PARTIAL** — degrades off the development set; report as a limit, do not claim generality |
| **F < 0.40** | ⛔ **FAILS** — the rule is NPIP-shaped. Say so plainly; do not retune and re-run |

Secondary, also committed now:

- **precision ≥ 0.70 pooled.** Families we emit should not be junk even where recall is hard.
- **≥ 1 truth family per held-out chromosome recovered at F = 1.000.** If the method never nails a single
  held-out family exactly, "it finds families" is too strong a claim.
- ⛔ **A chromosome yielding 0 families is not automatically a failure** — chr20 correctly yields 0
  (median `cov_longer` 0.0026, below the 0.30 floor). It counts as a failure **only if** the truth set
  for that chromosome is non-empty, which is checked and reported either way.

## 6. What I will not do

- Not tune any flag per chromosome, and not introduce a new flag.
- Not drop a held-out chromosome after seeing its score.
- Not substitute a different truth definition after seeing the result. If this truth is judged unfit, the
  whole test is reported as void — not replaced with one that scores better.
- Not quote the pooled number without the per-chromosome table, including the worst chromosome.

---

# TEST 1 IS VOID — and TEST 2, pre-registered before any Soto scoring

## Test 1 verdict: the truth was unfit, exactly as §6 warned

Scored on the **development** chromosomes only (chr2/chr6 were never opened for this):

| | truth families | touched | mean F | exact |
|---|---|---|---|---|
| chr16 (development, NPIP) | 55 | 19 | **0.215** | 1 |
| chr17 (development, TBC1D3) | 86 | 32 | **0.220** | 3 |

The development set scores **below this pre-registration's own FAIL bar (0.40)**, which means the truth,
not the method, is what failed. The evidence is in what it counts as a "family": of chr16's 36 wholly
missed truth families the largest are **MIR (46 members), RNA5S (31), ZNF (29), C16orf (12), SNORA,
SNORD, RN7SK, CDH, PRSS** — microRNA and rRNA gene classes, ancient protein-coding families, and
`C16orf`, which is not a family at all but the naming convention "open reading frame on chromosome 16".
The method's stated scope is *recent, sequence-similar* multi-copy families (identity ≥ 0.7,
cov_longer ≥ 0.3); it is correct to leave those alone. Meanwhile the families it does recover are exactly
the real recent duplications: **NOMO F = 1.000**, HERC2, SMG1 (7/7), HBA, CLEC18, ABHD17A.

⚠ Filtering by GFF biotype (external annotation, not our output) removes only **9 of the 36** misses
(6 miRNA, 2 snoRNA, 1 tRNA) and leaves 21 protein-coding nomenclature classes, so **no symbol-root truth
can be rescued by scoping**. Per §6 of this pre-registration the test is reported **VOID**, not replaced
with a better-scoring variant of itself.

⭐ One real finding survives it: with the amended single-pass root, **NPIPA (7 members) is scored as
wholly missed** — the method puts all 21 NPIP genes in ONE cluster, so bipartite matching gives that
cluster to NPIPB and leaves NPIPA unmatched. That is a genuine limitation and it is advisor question
**Q9 ("NPIPA and NPIPB should be distinct subfamilies")** reappearing as a measurement.

## Test 2 — external, duplication-aware truth on chromosomes outside the threshold's selection set

**Truth: Soto et al. 2025 published families** (`bench/soto/soto_famCN_S1C.tsv`, `Family ID` column) —
an external, published, SD-derived family assignment on the same CHM13 assembly. A truth family is one
with **≥ 3 members on the chromosome being scored**, matched to RefSeq by `Gene Name`.

**Why these chromosomes.** `--min-shared-exon-frac 0.60` was chosen by LORO on **chr5 / chr7 / chr21**
(§6o2–§6o6, `bench/FEX_SWEEP_LORO.md`), scoring **F 0.845 / precision 0.898** — and the ledger itself
records the caveat *"three regions is a weak LORO"*. So:

| arm | chromosomes | Soto families | why |
|---|---|---|---|
| **development** | chr5, chr7, chr21 | 3 + 12 + ? | the LORO set the 0.60 threshold was selected on |
| ⭐ **held out** | **chr2, chr8** | 7 + 5 | **zero mentions in both the ledger and the register** |
| ⭐ **held out** | chr10 | 8 | 1 ledger mention, 0 register — effectively unexposed |
| reference | chr16, chr17 | 15 + 15 | development families (NPIP, TBC1D3) |

**Both arms are scored by the SAME scorer** (`bench/heldout_family_score.py`, bipartite
sensitivity/precision/F) so the comparison is not against a differently-computed published number.

**The bar, committed now:**

| held-out pooled F (chr2 + chr8 + chr10) vs development pooled F (chr5/7/21) | verdict |
|---|---|
| held-out ≥ development − 0.10 | ⭐ **HOLDS** — the threshold was not fitted to its selection set |
| development − 0.25 ≤ held-out < development − 0.10 | ⚠ **PARTIAL** — real degradation; report as a limit |
| held-out < development − 0.25 | ⛔ **FAILS** — overfitted to chr5/7/21. Say so; do not retune |

Secondary: **≥ 1 exact (F = 1.000) family recovery on a zero-exposure chromosome** (chr2 or chr8).

⚠ **Stated limitation before the run.** Soto is the truth-source the 0.60 threshold was tuned against,
so what is held out here is the **chromosome**, not the truth-source. LORO is designed against exactly
the per-family overfitting at issue, but this is not an independent-truth test and will not be described
as one.
