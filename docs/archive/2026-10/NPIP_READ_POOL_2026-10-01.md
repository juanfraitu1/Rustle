# Which alignments build the loci, on NPIP (human chr16, A119b): result, 2026-10-01

Pre-registration: `docs/PREREG_npip_read_pool_2026-10-01.md` (commit e5530d37, before any arm ran). Scorers: `bench/npip_read_pool/score.py`
(pre-registered numbers), `bench/npip_read_pool/figdata.py` (figure windows), `bench/npip_read_pool/cointoss.py` (tied primaries, added after the
prereg). Work dir `/mnt/linuxdisk/tmp/readpool_npip/` (all outputs; drivers `bench/npip_read_pool/arm.sh`, `posthoc.py`, `pagedata.py`, page template `page.tmpl.html`). NPIP copies = the 25 chr16 copies of the
CAT/Liftoff v2.0 truth (Amendment 2). Human only; nothing here is pooled with gorilla.

## Arms and cost

| | P (primaries) | GOOD (+ secondaries AS >= 0.98 best; default since 09-24) | ALL (+ every secondary) |
|---|---|---|---|
| transcripts on chr16 | 8,673 | 9,473 | 14,183 |
| loci on chr16 | 2,550 | 2,802 | 5,826 |
| all-vs-all records | 72,969 | 95,530 | 1,604,462 |
| all-vs-all wall time | 39 s | 46 s | ~20 min in 13 shards (one process did not finish in 10 min) |
| multi-member clusters | 70 | 110 | 277 |

GOOD reproduces the 2026-09-29 copy-recovery run exactly (9,473 chr16 transcripts).

## Locus definition at the 25 NPIP copies

| | P | GOOD | ALL |
|---|---|---|---|
| loci on copies (rep exons overlap a copy's exons, same strand) | 31 | 46 | 169 |
| copies with >= 1 locus | 24 | 22 | **25** |
| loci per covered copy, median / max | 1 / 3 | 1 / 8 | 5 / 25 |
| loci on >= 2 copies by rep exons (pre-registered "fused") | 0 | 0 | 6 |
| loci whose span covers exons of >= 2 copies (post hoc) | 0 | 5 | 29 |
| echo loci on copies (zero primaries over the span) | 0 | 0 | 0 |

- P misses only NPIPB2 (an assembly-polish loss, registered 09-29). ALL builds a locus at every copy, NPIPB2 included.
- GOOD loses NPIPA6 and NPIPB6 as separate loci: secondary alignments join them to a neighbour. At NPIPB4/NPIPB5 a 2-read transcript
  has one 429,587-bp "intron" from an NPIPB4 exon to an NPIPB5 exon, a read split across two copies by its secondary alignment. At
  NPIPA6 the locus becomes one 76 kb locus whose representative is on NPIPA7.
- ALL piles up loci: 25 at NPIPB4, 20 at LOC128966608, 16 at NPIPB5, plus antisense and intronic pieces at most copies.
- The pre-registered echo test (zero primaries over the span) never fires at NPIP: every copy has primaries, most of them tied
  (below), so that test cannot see secondary-built loci here.

## Does the all-vs-all + clustering remove what ALL adds?

**Pre-registered rule: f = 61 / 2,370 = 0.026 -> "claim holds".** Of ALL's 2,431 loci that overlap no GOOD locus, 2,370 are echoes
or off-copy (by the rule's definition), and 61 of them end up in the largest NPIP cluster.

**The rule measured the wrong thing.** Its denominator is every extra locus on chr16, almost all of them far from NPIP, and those can
never join an NPIP cluster. It does not ask whether the NPIP clusters are polluted. Described without a rule:

- **What the step does remove:** 1,573 of the 2,431 extra loci (65%) align to nothing and end as singletons; folding (records on top of
  a member become that member) absorbs 76 loci into the largest NPIP cluster.
- **What it keeps:** every cluster that holds a locus on an NPIP copy, counted as nodes after folding:

| | P | GOOD | ALL |
|---|---|---|---|
| clusters | 5 | 5 | 12 |
| nodes | 38 | 51 | 147 |
| on an NPIP copy | 27 | 29 | 60 |
| inside a copy's span, off its exons | 2 | 3 | 23 |
| antisense to a copy | 2 | 8 | 25 |
| elsewhere | 7 | 11 | 39 |
| copies with a node | 23 | 21 | 24 |

  Largest NPIP cluster: P 22 nodes (20 on copies), GOOD 27 (22), ALL 27 (17 on copies, 9 with no gene, 1 on another gene; 76 loci
  folded in). With ALL, fewer copies keep a node of their own (17 vs 20) because loci spanning two copies fold copies together.
- **Why homology cannot do it:** a secondary alignment exists only where the read's sequence is homologous, so a locus built from
  secondaries is homologous to the locus built from the same reads' primaries by construction. The all-vs-all finds that homology and
  keeps the edge. It removes loci with no homolog (the 65% singletons), never the echoes of a real family.

## The coin toss is real at NPIP (`bench/npip_read_pool/cointoss.py`)

Primaries whose aligned blocks overlap a copy's exons, and whether the read is AS-tied genome-wide (second AS >= 0.98 x best AS, from
the `as_table` molecules table): **6,046 of 10,305 (58.7%) are tied** (genome-wide A119b: 2.5%). **15 of 25 copies are majority-tied**
(NPIPA6 97.9%, NPIPA7 98.8%, NPIPA9 98.5%, NPIPB15 98.7%, LOC124907807 99.0%, LOC124907808 100%). **LOC124907808 (0 untied) and
LOC124907807 (1 untied) have primaries-only loci built from coin tosses.** NPIPB2 (0.5%), PKD1P6-NPIPP1 (1.5%), NPIPA2 (4.2%) and
NPIPB7 (7.5%) are nearly untied. This revises r1061 for this family: genome-wide 0.3-0.5% of loci are majority-tied, NPIP 60%.

## Reading

- The advisor's premise holds at NPIP: primaries-only loci at most copies rest on coin tosses.
- His remedy does not: ALL keeps every tied copy (25/25), but it also builds 5 loci per copy, 87 non-copy nodes in the NPIP clusters, and
  29 loci spanning two copies. The all-vs-all keeps them all, since they are homologous by construction. It costs 22x the all-vs-all records and ~30x the time.
- GOOD (the default) is not clean on NPIP either: 21 copies with a node vs P's 23, 22 non-copy nodes vs 11, and five two-copy loci from
  split secondary alignments. Its genome-scale held-out gain (§6z7 r1060, gorilla) is not overturned by one family, but NPIP is where
  it costs.
- Where a read came from among tied copies is decided downstream, by PSVs (O2: assign or abstain), not by which pool built the loci.
