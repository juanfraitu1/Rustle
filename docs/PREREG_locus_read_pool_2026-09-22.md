# Pre-registration — which alignments define the INITIAL loci: primaries only, all alignments, or primaries + near-tied secondaries?

**Written 2026-09-22 (§6z7), before any arm is scored.** Advisor (via the user): building the initial loci from
PRIMARY alignments only is unsafe, because for AS-tied multi-mappers the primary placement is a coin toss. The
pipeline answers the *missing-locus* half later (O2 assigns tied reads by PSV; `--tied-seed` seeds starved
copies), but the other half — **a coin-toss primary that seeded a locus it had no business in** — is not
addressed anywhere. User's two asks: **(1)** measure how neutral or detrimental it is to define the initial
loci from ALL alignments; **(2)** since minimap2 ran with `-p 0.1`, many secondaries are bad — test an early
filter that defines loci from primaries + GOOD secondaries.

## What is already known (not predictions)

- §6n2 / `RUSTLE_GTF_SECONDARY=1`: admitting every secondary into pass-1 lifted complete NPIP chains 10→20/26
  at a **9× transcript cost** (echoes of one molecule at several paralogues). "O1 only, never quantification."
- r850 / `RUSTLE_GTF_SECONDARY_AS_RATIO`: the per-REGION AS-tie filter is inert (transcripts 13,132 at every
  width) because 83.4% of molecules contribute one record per region, so their in-scope best is themselves.
  **Verdict then: needs a genome-wide best-AS pre-pass per molecule.** That pre-pass is what ask (2) builds.
- ⚠ Found while wiring it: `copy_assign` always takes `BamIndexCache::reads_in_region` (the cache opens
  whenever a `.bai` exists), and that path applied **no** AS-tie filter at all — so the ratio env never
  reached the assembler through `copy_assign`. Fixed (both paths now filter; ratio 0 = unchanged).

## Arms — only the read pool that seeds pass-1 changes

| arm | pass-1 read pool |
|---|---|
| **P** (shipped) | primary records only (`-F 2308`) |
| **ALL** | primary + every secondary (`RUSTLE_GTF_SECONDARY=1`) — ask (1) |
| **GOOD** | primary + secondaries with `AS ≥ 0.98 × best AS of that molecule over the WHOLE BAM` (`RUSTLE_GTF_SECONDARY_AS_TABLE` from one scan of all records; ratio 0.95 and 0.90 reported as a descriptive sweep, **0.98 is the decision arm** — the same width `copy_assign --as-ratio` uses for "tied") — ask (2) |

Everything downstream is identical and is the shipped de novo O1 chain: `copy_assign --assemble-only` with the
shipped polish → one locus per `gene_id` (span = its transcripts' hull, exons = the rep's, `pick_locus_rep`)
→ `minimap2 -x asm20 -c -X -N 50 -p 0.1` all-vs-all → `mcl_families --min-exonic-bp 1 --min-shared-exon-frac
0.60` → `family_score` (sensitivity / precision / one-to-one bipartite F / collapsed).

## Substrates

- **Development: human A119b chr20** (all 1,104,846 records on the contig: 414,711 primary, 680,620
  secondary), truth = the chr20 protein referee (42 families / 152 genes, `/mnt/linuxdisk/tmp/sedef/truth/chr20.tsv`).
- **Held out: gorilla `GGO_mm.bam` NC_073244.2** (473,231 records), truth = a protein referee built the same way
  from `GGO_genomic.gff` + `GGO.fasta` **before any gorilla arm is run**. Reported with no re-tuning; never
  pooled with human.

## Metrics — committed now

**Primary (family level, the objective):** referee sensitivity, precision, bipartite F, collapsed, per arm.

**Locus level (the mechanism):** loci per arm; loci with ≥ 1 annotated gene/pseudogene overlap; **fill-in** =
ALL/GOOD loci absent from P (no P locus overlapping) and how many of those overlap an annotated
gene/pseudogene; **echo** = an ALL/GOOD locus whose supporting reads are all secondary records; **lost** = P
loci with no ALL/GOOD locus overlapping.

**Coin-toss audit on P (no arm needed):** for every P locus, the fraction of its supporting primaries whose
molecule is genome-wide AS-tied (`second-best ≥ 0.98 × best`). Report the share of loci that are
majority-tied, and of those, the share with **no** annotated gene/pseudogene overlap — the loci the coin toss
could have put in the wrong place. Also the same share restricted to loci that end up in a family.

**Transcript level (reported, not decisive):** gffcompare vs RefSeq, matching chains and precision — the
9× echo cost of §6n2 will show here and is the reason ALL is O1-only.

## The bar

Judged on **held-out gorilla**, family-level F vs P:

| outcome | verdict |
|---|---|
| ALL or GOOD: F up by ≥ 0.02 and collapsed not up | ⭐ **BETTER** — the read pool should change |
| \|ΔF\| < 0.02 | ⚠ **NEUTRAL** — then GOOD is preferred over ALL only if it cuts ALL's locus/transcript inflation by ≥ 50% at the same F |
| F down by ≥ 0.02 | ⛔ **DETRIMENTAL** |

**Predicted, before looking:** ⚠ NEUTRAL for both on the family level. Multi-copy families are exactly where
tied reads live, so ALL will add paralogue nodes (§6n2's +10 NPIP copies), but the graph already reaches most
of those copies from DNA homology of the loci that DO have a primary, so F should barely move. Locus
inflation under ALL will be large (echo loci at every tied paralogue, most of them at annotated
genes/pseudogenes — that IS the fill-in) and GOOD at 0.98 should retain most of the fill-in with far fewer
echoes, because `-p 0.1` secondaries are mostly far below the best. The coin-toss audit should show
majority-tied P loci to be a small minority overall and concentrated in family members.

I will not change the arms, the substrates, the metrics or the bar after seeing any number.

---

# OUTCOME, part 1 (2026-09-22) — ask (1), ALL alignments vs primaries only

## Family level (the objective; `family_score`, protein referee)

| substrate | arm | nodes / edges | clusters | sens | prec | **F** | collapsed | truth genes with no clustered locus |
|---|---|---|---|---|---|---|---|---|
| chr20 (dev) | **P** | 503 / 550 | 6 | 0.039 | **0.857** | 0.075 | **1** | 143 / 152 |
| chr20 (dev) | ALL | 3,378 / 79,845 | 39 | **0.164** | 0.568 | **0.255** | 18 | 92 / 152 |
| NC_073244.2 (held out) | **P** | 216 / 192 | 23 | 0.015 | **1.000** | 0.030 | **12** | 705 / 775 |
| NC_073244.2 (held out) | ALL | 881 / 5,779 | 64 | **0.071** | 0.724 | **0.129** | 99 | 559 / 775 |

**ALL is neither neutral nor simply detrimental: it is a sensitivity-for-precision trade** — F +0.18 (dev) and
**+0.10 (held out)**, precision −0.29 / −0.28, and collapsed truth genes 1 → 18 / 12 → 99. By the bar the ⭐
row needs "collapsed not up", which fails on both substrates, so ALL is **not adoptable as written**; it is
also not ⚠ (|ΔF| ≥ 0.02) nor ⛔ (F is up). ⚠ Both F values are low in absolute terms because this pure de novo
chain reaches few referee genes (P leaves 143/152 and 705/775 without a clustered locus, the alignment
ceiling of r1026); the relative comparison is what the arms test.

## Where ALL's gain comes from — mostly NOT invisible loci

| substrate | referee genes newly covered by a clustered locus under ALL | P had NO locus there ("invisible via primaries") | P had a locus there that was simply not clustered |
|---|---|---|---|
| chr20 | 52 (9 → 60, lost 1) | 12 (**23%**) | 40 (**77%**) |
| NC_073244.2 | 147 (70 → 216, lost 1) | 65 (**44%**) | 82 (**56%**) |

The advisor's "loci invisible because their reads' primaries went elsewhere" is real but the minority
(23–44% of the gain). The majority is **connectivity**: a secondary-seeded locus carries the SAME intron chain
as the primary locus of the same molecules, placed at the paralogue, so its edge to that locus passes the
exon conjunct trivially — where the paralogue's own primary-built locus (different rep, different exons) had
failed it. Edges go 550 → 79,845 and 192 → 5,779; clustered loci 69 → 850 and 63 → 389, of which 376 (44%)
and 117 (30%) overlap **no** annotated gene or pseudogene.

## Locus level (mechanism)

| substrate | arm | loci | median span | overlapping an annotated gene/pseudogene | loci absent from the other arm | … of which annotated | echo loci (zero primary records over the span) | P loci lost |
|---|---|---|---|---|---|---|---|---|
| chr20 | P | 1,856 | 6,947 bp | 82.3% | 0 | – | – | – |
| chr20 | ALL | 4,956 | 4,616 bp | 65.4% | 1,494 (30.1%) | 403 (27.0%) | 192 (3.9%) | 0 |
| NC_073244.2 | P | 1,113 | 15,416 bp | 99.0% | 0 | – | – | – |
| NC_073244.2 | ALL | 1,920 | 10,198 bp | 90.2% | 392 (20.4%) | 226 (57.7%) | 166 (8.6%) | 0 |

Only 3.9% / 8.6% of ALL's loci are true echoes (no primary record anywhere over their span). Most of the new
loci sit where primaries already were but had not reached a locus (support below the floor, or removed by the
polish); under ALL the secondaries lift them over it. **P loses no locus in either species.**

## Transcript level (reported, not decisive — the §6n2 echo cost)

| substrate | arm | transcripts | intron chain SN/PR | matching chains |
|---|---|---|---|---|
| chr20 | P | 6,121 | 24.9 / **19.2** | 1,069 |
| chr20 | ALL | 13,523 | 25.3 / 8.5 | 1,084 |
| NC_073244.2 | P | 4,279 | 28.3 / **36.8** | 1,573 |
| NC_073244.2 | ALL | 8,197 | 28.9 / 19.9 | 1,606 |

+15 / +33 matching chains for a halved precision, and the chr20 all-vs-all took 2,783 s against P's 102 s.

*Part 2 (ask 2, the GOOD arm) and the coin-toss audit follow once the genome-wide AS table is built.*

---

# OUTCOME, part 2 (2026-09-23) — ask (2), GOOD secondaries, and the coin-toss audit

Genome-wide table: one scan of every record (`samtools view -F 2052 | awk | sort`, 113 min for the 96 GB A119b
BAM, 17 min for gorilla) → per molecule best AS, second-best AS, primary contig. **Tie population: 2.5% of
A119b molecules are AS-tied at 0.98 (5.2% at 0.90; 27.3% have a second placement at all); gorilla 1.1% / 3.9%.**

## Family level

| substrate | arm | nodes / edges | clusters | sens | prec | **F** | collapsed | no clustered locus |
|---|---|---|---|---|---|---|---|---|
| chr20 (dev) | P | 503 / 550 | 6 | 0.039 | 0.857 | 0.075 | 1 | 143 |
| chr20 (dev) | **GOOD 0.98** | 538 / 801 | 7 | 0.053 | **0.889** | **0.099** | **0** | 142 |
| chr20 (dev) | GOOD 0.95 / 0.90 | 586 / 1,327 · 627 / 1,960 | 7 · 7 | 0.053 | 0.889 | 0.099 | 0 | 142 |
| chr20 (dev) | ALL | 3,378 / 79,845 | 39 | 0.164 | 0.568 | 0.255 | 18 | 92 |
| NC_073244.2 (held out) | P | 216 / 192 | 23 | 0.015 | 1.000 | 0.030 | 12 | 705 |
| NC_073244.2 (held out) | **GOOD 0.98** | 245 / 424 | 24 | 0.046 | **1.000** | **0.089** | **12** | 682 |
| NC_073244.2 (held out) | GOOD 0.95 / 0.90 | 256 / 450 · 301 / 1,993 | 25 · 25 | 0.048 · 0.046 | 1.000 | 0.091 · 0.089 | 11 · 10 | 680 |
| NC_073244.2 (held out) | ALL | 881 / 5,779 | 64 | 0.071 | 0.724 | 0.129 | 99 | 559 |

⭐ **GOOD 0.98 meets the ⭐ row on the held-out substrate: F +0.059 (≥ 0.02), collapsed 12 = 12, and precision
stays at 1.000** (chr20: F +0.024, precision UP 0.857 → 0.889, collapsed 1 → 0). It keeps roughly half of
ALL's recall gain at none of its cost, and the ratio does not matter at the family level (0.98–0.90 give the
same chr20 score and F .089/.091/.089 on gorilla).

## Locus and transcript level — GOOD is nearly free

| substrate | arm | loci | new vs P (annotated) | echo loci | transcripts | intron chain SN/PR | matching chains |
|---|---|---|---|---|---|---|---|
| chr20 | P | 1,856 | – | – | 6,121 | 24.9 / 19.2 | 1,069 |
| chr20 | GOOD 0.98 | 1,896 (+2.2%) | 35 (8) | 13 (0.7%) | 6,175 (+0.9%) | 24.9 / 19.1 | 1,069 |
| chr20 | ALL | 4,956 (+167%) | 1,494 (403) | 192 (3.9%) | 13,523 (+121%) | 25.3 / 8.5 | 1,084 |
| NC_073244.2 | P | 1,113 | – | – | 4,279 | 28.3 / 36.8 | 1,573 |
| NC_073244.2 | GOOD 0.98 | 1,145 (+2.9%) | 30 (**29**) | 17 (1.5%, all annotated) | 4,329 (+1.2%) | 28.7 / **37.0** | **1,597** |
| NC_073244.2 | ALL | 1,920 (+72%) | 392 (226) | 166 (8.6%) | 8,197 (+92%) | 28.9 / 19.9 | 1,606 |

On gorilla the 30 loci GOOD adds are annotated genes/pseudogenes in 29 cases — the genuinely invisible copies —
and the transcript set gains 24 matching chains at higher precision. All-vs-all: 174 s (GOOD) vs 2,783 s (ALL)
on chr20.

## The coin-toss audit (P loci, genome-wide ties at 0.98)

| substrate | P loci | majority-tied (≥ 50% of supporting primaries tied) | … unannotated | fully tied | loci in a family | majority-tied among them |
|---|---|---|---|---|---|---|
| chr20 | 1,856 | **5 (0.3%)** | **0** | 0 | 69 | 0 |
| NC_073244.2 | 1,113 | **6 (0.5%)** | **0** | 1 | 63 | 4 (all annotated) |

Median tied fraction per locus is 0.000 on both. ⛔ **The feared failure — a locus seeded by coin-toss
primaries where nothing belongs — does not occur at a measurable rate**: a locus is never built from tied reads
alone, and the handful that are majority-tied all sit on annotated genes. The coin toss decides *which*
paralogue receives a molecule, so its only cost is the invisible-locus half, which is what GOOD repairs.

## Verdict

- **Ask (1), ALL alignments:** a sensitivity-for-precision trade with a collapse explosion and 27× compute — not
  neutral, not adoptable as written (register 1059).
- **Ask (2), primaries + near-tied secondaries at the genome-wide best:** ⭐ better on held-out by the bar,
  free at the locus and transcript level (register 1060). Shipped as
  `RUSTLE_GTF_SECONDARY=1 RUSTLE_GTF_SECONDARY_AS_RATIO=0.98 RUSTLE_GTF_SECONDARY_AS_TABLE=<molecules.tsv>`
  (byte-identical unset). **Not the default yet** because it needs the genome-wide pre-pass; the follow-up is a
  Rust `as_table` bin (one BAM scan → the table) so the recipe has no `samtools|awk|sort` step, and then the
  user's call on the flip.
- **The advisor's mis-seeding concern:** negligible at the locus level (register 1061).
