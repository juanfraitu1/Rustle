# Two candidate definitions tested: cover prediction, and J_N as an MCL weight

**§6u8, 2026-09-21.** Pre-registration `docs/PREREG_cover_and_jn_weight_2026-09-21.md` (md5
`eed4b0fc`), committed `6f18304b` before any score. Tool `bench/cover_and_jn_definition.py`.
**NEITHER IS ADOPTED.** ⚠HUMAN substrate (CHM13/Soto) — do not pool with gorilla.

## 0. A correction to how the truth has been read

Soto's S1C carries an explicit **`No. Assigned Families`** column, and **149 of 2,334 CHM13 gene IDs
(6.4%) are assigned to more than one family**. The derived family-set size equals that column for
**2,334 / 2,334** records — so the published truth is a **COVER, authored as one**.

Every scorer in this project applies `if len(fids) != 1: continue` (*"a partition needs one label per
gene"*). Because dropping a member can push a family under the >= 3-member floor, that removes whole
families, not just genes:

| chromosome | cover truth | partition truth (what was used) |
|---|---|---|
| chr16 | **15 families / 70 members** | 8 / 43 |
| chr7 | 12 | 9 |
| chr9 | 7 | 5 |
| chr5 | 3 | 2 |

Everything below is scored against the **cover** truth. ⚠Keying is by Gene Name because no CHM13 CAT
annotation (`CHM13_G*`) is available to reach loci; the pre-registered sensitivity arm (drop the 21
names that are multi-family only by gene-ID collision) agreed in sign on both rules.

## 1. Setup

Graph = the shipped `mcl_families --dump-graph` (edge set and exon conjunct shipped). Comparator =
**mcl_port-MCL on the same graph**, never the shipped Rust F (register 917: mcl_port is not
bit-identical). Six new chromosome graphs were built for this test (chr4, chr7, chr9, chr15, chr17,
chr22), taking the development truth from **35 to 76 families**.

Development: chr2, chr8, chr10, chr16, chr15, chr17, chr22, chr7 — **76 cover families**.
Held-out: chr4 + chr9 (3 and 16 doc mentions; chr19 was rejected — it has ONE family >= 3 members).

## 2. Result — development (76 families)

| arm | sens | prec | F | dF | better/worse/unchanged |
|---|---|---|---|---|---|
| **baseline MCL** | 0.7029 | 0.6150 | **0.6197** | — | — |
| cover k>=2 | 0.7379 | 0.5870 | 0.6152 | **-0.0045** | 6 / **17** / 53 |
| cover k>=3 | 0.7239 | 0.5940 | 0.6138 | -0.0059 | 3 / 13 / 60 |
| cover k>=4 | 0.7060 | 0.6014 | 0.6142 | -0.0055 | 0 / 7 / 69 |
| J_N gate>=3 | 0.7037 | 0.6309 | 0.6299 | **+0.0102** | 5 / 2 / 69 |
| J_N gate>=5 | 0.7037 | 0.6309 | 0.6299 | +0.0102 | 5 / 2 / 69 |
| J_N gate>=8 | 0.7037 | 0.6269 | 0.6278 | +0.0081 | 3 / 2 / 71 |

Bar was **> +0.02**. Neither clears.

## 3. Test 1 (cover) — REFUTED, and the small-sample version was a sign flip

On the first four chromosomes (35 families) cover k>=2 scored **+0.0098**. Adding three chromosomes took
it to **-0.0064** (64 families) and the full set to **-0.0045** (76). The sign flipped with power, and
throughout, **more families get WORSE than better** (6 vs 17 at k>=2). Dual membership buys sensitivity
(0.7029 -> 0.7379) and pays more for it in precision (0.6150 -> 0.5870).

⟹ The partition is the right shape for the prediction, even where the truth is a cover. That answers
register 845's re-open condition (chr10 AGAP at 40% chimeras): the condition was met, the test was run,
and dual membership still does not pay.

## 4. Test 2 (J_N weight) — the whole gain is ONE family

+0.0102 pooled, precision +0.016 with sensitivity flat — but only 7 of 76 families move at all, and:

| chrom | family | n | baseline sens/prec/F | J_N sens/prec/F | dF |
|---|---|---|---|---|---|
| chr2 | **ID_176** | 4 | 0.750 / **0.130** / 0.222 | 0.750 / **1.000** / 0.857 | **+0.635** |
| chr10 | ID_270 | 5 | 0.400 / 0.400 / 0.400 | 0.600 / 0.429 / 0.500 | +0.100 |
| chr10 | ID_409 | 5 | 1.000 / 0.833 / 0.909 | 1.000 / 1.000 / 1.000 | +0.091 |
| chr15 | ID_78 | 4 | 0.250 / 0.200 / 0.222 | 0.250 / 0.333 / 0.286 | +0.063 |
| chr15 | ID_113 | 16 | 1.000 / 0.457 / 0.627 | 1.000 / 0.485 / 0.653 | +0.026 |
| chr16 | ID_154 | 16 | 0.938 / 0.556 / 0.698 | 0.938 / 0.536 / 0.682 | -0.016 |
| chr2 | ID_65 | 7 | 0.571 / 1.000 / 0.727 | 0.429 / 1.000 / 0.600 | -0.127 |

**ID_176 alone is +0.635 of the +0.772 total.** Without it the other six net +0.137, i.e. **+0.002
pooled** — nothing. This is register 917's trap in its purest form (*"it is TWO FAMILIES"*).

⭐The mechanism is real and worth keeping: ID_176's baseline cluster had **precision 0.130** — an
over-merged hub — and J_N cleaned it to **1.000**. That is exactly the large-component regime where
§6u3 measured J_N at AUC 0.924. **J_N fixes over-merged hubs; there is only one such family in 76.**

## 5. Held-out (chr4 + chr9) — structurally uninformative

Every arm is an **exact no-op**: dF 0.0000, 0 of 9 families moved, for both rules at every parameter.
The held-out arm could not have decided anything — the same failure mode as register 845, whose
held-out arm moved by exactly 0.000 because *"0 of the 15 chimeras is inside the scored universe"*.
⚠The decision was made on development; this arm is reported for completeness and cannot resurrect
either rule.

## 6. Consequences

- The definition is unchanged: `mcl_families --min-exonic-bp 1 --min-shared-exon-frac 0.60`.
- **J_N is confirmed a SCORER, not a definitional ingredient** — closing the last live lead from §6u3.
  It remains worth reporting as a per-cluster certificate flagging over-merged hubs (precision 0.130
  clusters), which is where its AUC 0.924 lives.
- The cover truth should replace the partition-restricted truth in future scoring: it is what Soto
  published, and the restriction was suppressing ~47% of chr16's families.
- Six reusable chromosome graphs now exist at `/mnt/linuxdisk/tmp/famgraph/`; the Soto truth universe is
  ~97 families genome-wide at >= 3 members, of which 76 are now in the development set — so the ceiling
  on statistical power for ANY future family-definition test is now known and nearly exhausted.

---

# 7. Soto is not the only truth — both candidates re-scored on an independent referee

**User, 2026-09-21:** *"Soto should not be the only ground truth, they have some epic fails in their
logic too and to my knowledge some of the family members are overmerged/undermerged, we should use it
as a middle ground."* Correct, and §6u5 (r934) had already flagged that every Soto-scored number is
**agreement**, not precision. Everything in §§1-6 above is therefore re-scored here.

**Second truth: the §6ko protein-family referee** — longest CDS per gene, translated in-house, all-vs-all
blastp e <= 1e-5, edge iff non-overlapping HSPs cover >= 0.30 of the longer protein, MCL I=2.8, r2
biotype exclusions. Independent of our genomic alignment gate AND of Soto's SD/WSSD construction (it is
amino-acid evidence about the product), though not of the RefSeq annotation it reads CDS from.
Symbol-root truth was NOT used: register 902 voided it (it calls `MIR`, `RNA5S`, `ZNF`, `C16orf`
naming classes "families").

It is also far larger: **289 families / 1,718 members** over the 8 development chromosomes, against
Soto's 76 / 342 — so it relieves the power ceiling of row 944.

## 7.1 Both candidates are negative on the independent truth too

| arm | sens | prec | F | dF | better/worse/unchanged |
|---|---|---|---|---|---|
| **baseline MCL** | 0.2759 | 0.4146 | **0.3028** | — | — |
| cover k>=2 | 0.2777 | 0.4079 | 0.2996 | **-0.0032** | 2 / 11 / 276 |
| cover k>=3 | 0.2769 | 0.4128 | 0.3025 | -0.0003 | 2 / 7 / 280 |
| cover k>=4 | 0.2767 | 0.4133 | 0.3027 | -0.0001 | 1 / 3 / 285 |
| J_N gate>=3 | 0.2759 | 0.4129 | 0.3015 | **-0.0013** | 2 / 4 / 283 |
| J_N gate>=5 | 0.2759 | 0.4129 | 0.3015 | -0.0013 | 2 / 4 / 283 |
| J_N gate>=8 | 0.2759 | 0.4146 | 0.3028 | -0.0000 | 1 / 2 / 286 |

⚠Absolute F is low by construction — protein families include ancient homology that our recent-duplication
gate (identity >= 0.7, cov >= 0.3) deliberately excludes — so this table ranks arms, it does not measure us.

⭐**J_N's +0.0102 Soto gain does NOT replicate: it is -0.0013 here.** The one family carrying that gain
(chr2 `ID_176`) is a Soto artifact of scale. **Both candidates are now refuted on two independent truths**,
which is a stronger result than §§2-5 alone.

⚠A bug was found and fixed en route: `protein_referee` numbers families `PF{i}` independently per
chromosome, so per-family diffs keyed on the bare id compared different families across chromosomes
(reported sum dF +10.0 against a pooled -0.0032). Keys are now chromosome-qualified. Pooled figures were
never affected.

## 7.2 Auditing Soto: under-merge YES, over-merge NO

On the **118 genes both truths place** (Soto pairs 343, referee pairs 590, our pairs 471):

| relation | pairs | reading |
|---|---|---|
| S ∩ R | **343** | every Soto pair survives the referee |
| **S \ R** (Soto joins, referee splits) | **0** | **no over-merge detected — 0 of 24 Soto families** |
| **R \ S** (referee joins, Soto splits) | **247** | under-merge candidates |

The over-merge half of the concern is **not supported**: the referee never splits a Soto family.

The under-merge half **is**, once the referee's coarseness is controlled for. Of the 247:

- **164 (66.4%) carry no sequence homology edge at all** — ancient paralogy only the protein referee
  sees. Soto is *right* to exclude these; its families are recent duplications.
- **83 (33.6%) carry a direct homology edge**, median aggregated identity **0.858**, with **71 at >= 0.80,
  37 at >= 0.90, 30 at >= 0.95** — recent duplicates by any standard, which Soto splits.
  (Control: pairs both truths join have median identity 0.962.)

## 7.3 What that costs us

| scored against | precision | recall |
|---|---|---|
| Soto alone | **0.7113** (TP 335, FP 136) | 0.9767 |
| Soto, not charging the 83 referee-corroborated pairs | **0.8875** | 0.9767 |

**61% of our apparent false positives against Soto (83 of 136) are corroborated by an independent
protein family AND a direct sequence edge.** 53 remain uncorroborated.

⚠**This is a bounded correction, not a new scoring truth.** The sequence edge is our own gate's evidence,
so folding these pairs into a truth we then score ourselves on would condition the denominator on the
prediction (register 770). It is a LOWER BOUND on how much Soto understates our precision, and it is the
right thing to quote beside any Soto number — never instead of one.

⟹ **Soto is a middle ground, as asked: a high-precision, low-coverage reference.** Its pairs are
trustworthy (0 over-merges, and we recover 97.7% of them), its omissions are not evidence of our error.
