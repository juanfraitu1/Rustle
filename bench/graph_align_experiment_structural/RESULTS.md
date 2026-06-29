# Graph vs Linear paralog copy-assignment — STRUCTURAL regimes

**Question.** A prior experiment (`../graph_align_experiment/RESULTS.md`) found
sequence-to-graph alignment **TIES** linear minimap2 (both 100%) on a SNP-only
2-copy paralog family — because clean substitution PSVs are signal minimap2
already exploits. The open question this experiment answers: does graph alignment
**WIN** when the two copies differ by **STRUCTURAL** variation (inversion / large
indel), the regime where linear seed-chaining is supposed to struggle and produce
split/chimeric alignments?

## VERDICT (TL;DR)

**No regime here shows graph alignment materially beating linear for paralog copy
assignment.** On both structural regimes, **all three methods correctly assign
100% of the *resolvable* reads** (graph, linear-AS, and even naive-primary). The
hypothesized failure mode — structural variants forcing minimap2 into
split/chimeric/soft-clipped alignments that break coordinate-based assignment —
**does not occur in the assignment setting**, because the read's *correct* copy is
present in the reference and aligns to it collinearly. minimap2 always has a clean
copy to land on, so it never needs to split.

The split DOES occur — but only in a control where the read's true copy is
**absent** from the reference (copyB reads forced onto a copyA-only reference):
then a 600 bp inversion produces 367–530 bp soft-clips and supplementary records.
That is not the copy-assignment task (both copies are in the reference by
definition), so it doesn't help the hypothesis.

The **one** place the methods differ is on **info-theoretically unresolvable
reads** (indel regime, reads that span neither the indel nor any SNP): graph and
linear-AS both honestly **abstain** (tie / zero votes), while naive-primary
**fabricates a confident wrong call** (9/12 wrong). That is an argument for
*AS-based or graph-based abstention over naive-primary*, **not** for graph over
linear — linear-AS abstains identically.

---

## Tools (verified)
vg v1.73.0 · GraphAligner 1.0.20 · minigraph 0.21 · minimap2 2.30-r1287 ·
samtools 1.22.1 · python3 3.14 + pysam 0.23.3.

## Benchmarks (single-exon transcript copies, ~1.8–2.3 kb; `generate.py`, seed 1234)

**Regime 1 — INVERSION** (`inversion/`): copyA = `[800bp left][200bp seg][800bp right]`;
copyB identical except the 200 bp internal segment is **reverse-complemented**,
plus 2 flanking SNPs (left@200, right@1600). 76 reads: 30 full-length, 30
inversion-spanning (±120–300 bp flanks), 16 non-spanning-but-carry-a-SNP. Read
names `uniqA_*`/`uniqB_*` encode true copy; class in name (`full`/`span`/`nospan_snp`).

**Regime 2 — LARGE INDEL + sparse SNPs** (`indel/`): copyA = `[1000][300bp INSERTION][1000]`
(2300 bp); copyB lacks the insertion (2000 bp); only 2 SNPs elsewhere. 72 reads:
24 full, 24 indel-spanning, 12 non-spanning-with-SNP, **12 non-spanning-NO-SNP
(info-theoretically UNRESOLVABLE — carry zero diagnostic signal)**.

## Methods
- **GRAPH**: `vg msga -f copies.fa` → GFA with named paths copyA/copyB →
  `GraphAligner -x vg` → GAF. Assign by (a) **exclusive-node vote** (nodes unique
  to one path = the SNPs / the insertion) **+** (b) for the inversion,
  **orientation through the flipped segment** (nodes shared between paths but
  traversed `+` in copyA vs `-` in copyB — the inversion bubble). Tie/zero → ambiguous.
- **LINEAR-AS**: `minimap2 -a -x map-pb -N 5 -p 0.1` to copyA and copyB
  **separately**; assign to the copy with higher `AS`. Exact tie → ambiguous.
- **LINEAR-NAIVE-PRIMARY**: `minimap2 -a -x map-pb` to a **concatenated** reference
  (both copies); assign by which copy the **primary** alignment lands on — what a
  coordinate-ingesting tool (rustle) effectively sees. Forces a call (never abstains).

---

## Did the graphs encode the structural variant? YES (both)

**Inversion graph** (`vg msga`, despite its band-width warning) modeled the
inversion correctly as an **orientation-flipped bubble**:

| metric | value |
|---|---|
| nodes | 65 |
| A-exclusive nodes / bp | 4 / 18 bp (the 2 SNPs + boundary) |
| B-exclusive nodes / bp | 3 / 18 bp |
| **shared-but-orientation-FLIPPED nodes** | **6 nodes / 184 bp** |

copyA path: `…29,30+,31+,32+,33+,34+,35+,38…`  copyB path:
`…29,35-,34-,33-,32-,31-,30-,38…` — i.e. copyB traverses the same 6 nodes in
**reverse**. That is a textbook inversion bubble, and it is the signal the graph
assignment uses.

**Indel graph**: 310 bp of **A-exclusive** sequence (the ~300 bp insertion +
1 SNP, fragmented into ≤32 bp nodes by msga's band width) and 10 bp B-exclusive
(the 2 SNPs). The insertion is captured as a copyA-only super-bubble.

---

## 1. Per-regime accuracy (n correct / n total), by read class

### Regime 1 — INVERSION (76 reads)
| read class (n) | GRAPH | LINEAR-AS | LINEAR-NAIVE-PRIMARY |
|---|---|---|---|
| full (30) | **30/30** | **30/30** | **30/30** |
| span (30) | **30/30** | **30/30** | **30/30** |
| nospan_snp (16) | **16/16** | **16/16** | **16/16** |
| **resolvable total (76)** | **76/76** | **76/76** | **76/76** |

### Regime 2 — INDEL (72 reads)
| read class (n) | GRAPH | LINEAR-AS | LINEAR-NAIVE-PRIMARY |
|---|---|---|---|
| full (24) | **24/24** | **24/24** | **24/24** |
| span (24) | **24/24** | **24/24** | **24/24** |
| nospan_snp (12) | **12/12** | **12/12** | **12/12** |
| **resolvable total (60)** | **60/60** | **60/60** | **60/60** |
| nospan_nosnp (12) — UNRESOLVABLE | 0 wrong, **12 abstain** | 0 wrong, **12 abstain** | **9 WRONG**, 3 lucky, 0 abstain |

**Decision margins** (so a 100% tie is distinguishable from a knife-edge — it is
NOT a knife-edge for any method):

| regime / class | GRAPH node-vote gap (min/mean/max) | LINEAR-AS AS gap (min/mean/max) |
|---|---|---|
| inv full | 9 / 9.5 / 10 | 586 / 916 / 1938 |
| inv span | 7 / 7.5 / 8 | 526 / 628 / 842 |
| inv nospan_snp | 1 / 1.0 / 1 | 6 / 6 / 6 |
| indel full | 3 / 7.4 / 12 | 336 / 627 / 936 |
| indel span | 1 / 5.5 / 10 | 288 / 615 / 918 |
| indel nospan_snp | 1 / 1.0 / 1 | 6 / 6 / 6 |

Every correct call clears its decision boundary comfortably. The single-SNP
classes (gap = 1 node / 6 AS) are the tightest but still unambiguous — and graph
and linear are tight *together*, not one rescuing the other. On
inversion-spanning reads the linear AS gap is **+526 to +842** (huge), because
aligning a read to the wrong copy turns the 200 bp inverted segment into a block
of ~200 noisy substitutions/indels → ~600+ AS penalty.

---

## 2. Do inversion-spanning reads get SPLIT by minimap2? Does it break naive-primary? Does graph recover?

**No, and therefore there is nothing to recover.**

For an inversion-spanning read aligned to the **wrong** copy, minimap2 keeps a
**single primary** alignment and represents the inverted middle as a run of
mismatches/indels (real example, true-copy-A read vs copyB):
`27M1I14M1I216M12D…<noisy inverted region>…251M`, `AS:i:770` vs `AS:i:1386` on the
correct copy. There is no chimeric split because the **correct** copy is also in
the reference and aligns collinearly with the full read.

In the **naive concat** reference (both copies present), across **all** inversion
span reads (and an adversarial set of 30 short-flank reads, 40–90 bp flanks +
full 200 bp inversion, and a 600 bp-inversion set):

| naive-primary diagnostic (inversion span, n=30) | count |
|---|---|
| primary mis-assigned to wrong copy | **0** |
| primary soft-clip ≥ 50 bp | **0** (max soft-clip = 0 bp) |
| reads with a **supplementary** record | **0** |
| reads with a secondary record | 1 |
| adversarial short-flank reads (n=30): primary wrong / suppl / heavy-clip | **0 / 0 / 0** |
| 600 bp inversion reads (n=24): primary wrong / suppl / heavy-clip | **0 / 0 / 0** |

So naive-primary is **not** broken by the inversion in the assignment setting, and
graph has nothing extra to recover — both are already 100%.

**Where the split actually happens (control):** force copyB reads (inverted
middle, 600 bp) onto a **copyA-only** reference (true copy absent). Now minimap2
**must** break: 12/12 reads get a **367–530 bp soft-clip** on the primary, and 2
also emit **supplementary** records (the inverted middle aligning in reverse). The
split is real — it just requires the read's matching copy to be **missing**, which
never happens in copy *assignment* (the task is precisely to choose between two
present copies).

For the **indel** regime: an insertion-spanning read on the wrong copy aligns as a
single clean `…147M 300I 60M…` (the 300 bp insertion = one I operation),
`AS:i:574` vs `1484` correct. No split, no soft-clip; correct copy wins by +910 AS.

---

## 3. Verdict per regime

- **INVERSION — TIE (graph does NOT beat linear).** Graph 76/76, linear-AS 76/76,
  linear-naive-primary 76/76. The inversion never forces a split/chimeric primary
  in the assignment setting, so naive-primary does not fail and graph buys nothing.
  The graph *correctly* encodes the inversion as a fwd/rev bubble (real, verified)
  — but linear already resolves every read with a +526…+842 AS margin.

- **INDEL — TIE on resolvable reads (graph does NOT beat linear).** Graph 60/60,
  linear-AS 60/60, naive-primary 60/60. The 300 bp insertion is a single CIGAR `I`
  op, not a split. The **only** divergence is on the 12 truly-unresolvable reads:
  graph and linear-AS both **abstain** (0 wrong); naive-primary forces calls and is
  **9/12 wrong**. This separates *abstaining methods* (graph **and** linear-AS) from
  *forced-call naive-primary* — it is **not** a graph-vs-linear distinction. The
  takeaway is "**use AS-tie / node-vote abstention instead of trusting the primary
  position**," and AS abstention is a linear technique.

---

## 4. Overall — is there ANY regime where graph is materially better?

**No.** Across SNP-only (prior experiment) and both structural regimes here, graph
alignment is **equal**, never better, for paralog copy assignment:

1. **The premise (structural SV → linear split → naive-primary fails) is false in
   the assignment setting.** A read's correct copy is in the reference and aligns
   collinearly; minimap2 lands the primary there with a large AS margin. SV only
   forces split/soft-clip when the matching copy is *absent* — not the assignment
   task.
2. **Inversions/indels make linear discrimination EASIER, not harder:** mapping to
   the wrong copy turns the SV into a wall of penalties (≈ inversion length of
   substitutions, or a full-length insertion op), inflating the AS gap to
   hundreds–thousands of points — a far wider margin than the SNP-only regime.
3. **The graph's structural modeling is correct but redundant:** the inversion
   bubble (6 flipped nodes / 184 bp) and the insertion super-bubble (310 bp) are
   built and used, yet only reproduce assignments linear already makes.
4. **The genuine win is abstention, not graph:** on unresolvable reads, both graph
   and linear-AS abstain; only naive-primary fabricates wrong calls. So the
   actionable lesson for a coordinate-ingesting tool (rustle) is to treat
   AS-ties / zero-diagnostic reads as ambiguous rather than trusting the primary
   coordinate — achievable purely with linear AS, no graph required.

**Honest bottom line:** for 2-copy paralog assignment with the copies present in
the reference, graph alignment does not beat linear seed-chaining — not on SNPs,
not on a 200/600 bp inversion, not on a 300 bp indel. Structural divergence, if
anything, *helps* linear by widening the AS margin. A graph advantage would
require a different problem (matching copy genuinely absent/unassembled, or
cross-read PSV phasing of individually-ambiguous reads), not better
discrimination between two represented copies.

---

## Exact commands & tool quirks (reproducible)

```bash
cd bench/graph_align_experiment_structural
python3 generate.py                 # -> inversion/{copies,reads}.fa, indel/{copies,reads}.fa

for reg in inversion indel; do
  cd $reg
  # graph
  vg msga -f copies.fa > graph.vg 2> msga.log   # prototype; "length>128 band" warning is HARMLESS
  vg view graph.vg > graph.gfa
  vg validate graph.vg                            # "graph: valid"
  GraphAligner -g graph.gfa -f reads.fa -a aln.gaf -x vg
  # linear-AS (per-copy)
  samtools faidx copies.fa copyA > copyA.fa; samtools faidx copies.fa copyB > copyB.fa
  minimap2 -a -x map-pb -N 5 -p 0.1 copyA.fa reads.fa > linA.sam
  minimap2 -a -x map-pb -N 5 -p 0.1 copyB.fa reads.fa > linB.sam
  # linear-naive (concat, primary position)
  cat copyA.fa copyB.fa > concat.fa; samtools faidx concat.fa
  minimap2 -a -x map-pb concat.fa reads.fa > naive.sam
  cd ..
  python3 evaluate.py $reg            # -> $reg/eval_output.txt
done
```

**Tool quirks observed**
- `vg msga` is a deprecated prototype but builds correct 2-path graphs, **including
  the inversion as an orientation-flipped bubble** (nodes shared between paths,
  `+` in copyA / `-` in copyB). The "sequence length > 128 non-chunked band" warning
  is harmless for these short copies; alignments validated clean.
- **GraphAligner GAF path orientation IS the inversion signal.** Assignment must
  (i) normalize each read's global strand from non-flipped shared nodes, then
  (ii) compare its orientation at the flipped nodes against copyA's — otherwise a
  reverse-strand read inverts every vote. (Handled in `evaluate.py`.)
- **`minimap2 -x map-pb` does not split structural-variant reads when a collinear
  copy exists.** It absorbs a 200/600 bp inversion or a 300 bp indel into a single
  primary alignment (noisy block / one big I op) rather than emitting a
  supplementary. Splitting only appeared when the matching copy was *removed* from
  the reference (`absent.sam`: 367–530 bp soft-clips, some supplementaries).
- **`-N 5 -p 0.1`** is required for linear-AS so the second-best (wrong-copy) hit is
  reported and the AS gap can be measured; without it some reads report only the
  winning copy (still correct, but no margin).
- Linear-AS abstention = exact AS tie; verified to fire on all 12
  info-theoretically unresolvable indel reads (AS_A == AS_B to the point).

## Artifacts
`generate.py`, `evaluate.py`, this `RESULTS.md`, and per-regime `inversion/` &
`indel/` dirs each with `copies.fa`, `reads.fa`, `graph.{vg,gfa}`, `aln.gaf`,
`linA.sam`, `linB.sam`, `naive.sam`, `meta.txt`, `eval_output.txt`. Inversion dir
also has the stress-test controls: `adv_reads.fa`/`adv_naive.sam` (short-flank),
`big_copies.fa`/`big_reads.fa`/`big_naive.sam` (600 bp inversion), `absent.sam`
(true-copy-absent split control).
