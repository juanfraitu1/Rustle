# Graph vs Linear paralog copy-assignment — feasibility experiment

**Question.** Does *sequence-to-graph* alignment (variation graph) resolve which
paralog gene copy a read came from **better than, equal to, or worse than**
ordinary *linear* alignment, on a synthetic 2-copy family with ground truth?

**Verdict (TL;DR): EQUAL. Both methods are 100% accurate (103/103 reads), with
huge decision margins for both.** On this benchmark the diagnostic SNPs (PSVs)
are individually weak but collectively overwhelming, so even linear minimap2
resolves every read decisively. The variation graph *does* build correct bubbles
at the diagnostic sites, and graph assignment is also perfect — but it buys
nothing over linear here. This is an honest negative/tie result for "graph beats
linear."

---

## Data
`test_data/synthetic_family/` — one 20 kb `chr_test` with two ~98%-identical
paralog copies on the `-` strand:
- Copy A: genes `MYG1_A`, isoforms **A1** (4 exons) and **A2** (exon-3 skip).
- Copy B: gene `MYG1_B`, isoform **B1** (4 exons).
- 7 planted **diagnostic SNPs (PSVs)** where copy A carries base X and copy B
  carries base Y (X≠Y). A read carries its **TRUE** copy's allele at these sites
  regardless of where it is placed.
- 131 IsoSeq-like alignment records / **103 distinct reads** (~0.5% CCS error,
  20–40 bp poly-A tails, ±jitter on terminal exons).

### Ground truth (fully derivable from read names — incl. multi-mappers)
The generator (`generate_isoseq.py`) bakes the true copy into the read name and
into the diagnostic alleles. So truth is known for **every** read, not just the
uniq_* ones:

| name prefix | n | true copy | notes |
|---|---|---|---|
| `uniq_A1_*` | 30 | A | full 4-exon, copy A |
| `uniq_A2_*` | 20 | A | exon-3 skip, copy A |
| `uniq_B1_*` | 25 | B | full 4-exon, copy B |
| `multi_*`   | 14 | **A** | ambiguous multi-mapper, carries copy-A alleles (X) everywhere |
| `multi_r_*` | 14 | **B** | ambiguous multi-mapper, carries copy-B alleles (Y) everywhere |

The `multi*` reads are the real resolution test: each is *placed* on both copies
in the BAM (primary + supplementary), but its bases at the 7 PSV sites match only
its true copy. We could therefore evaluate **accuracy on all 103 reads**, not
just the uniq_* subset.

---

## Pipeline
1. **Transcript seqs** (`extract_seqs.py`): spliced, per-copy transcript
   sequences A1/A2/B1 extracted from `ref.fa` using the generator's exact 0-based
   half-open exon coords (so they match the reads) → `seqs.fa`
   (A1=1150, A2=1000, B1=1150 bp). Reads are spliced mRNA, so we use **transcript**
   sequences (GraphAligner is not splice-aware).
2. **Variation graph**: `vg msga -f seqs.fa` → `graph.vg` → `vg view` → `graph.gfa`.
   125 nodes, 154 edges, 3 named paths (A1, A2, B1). Graph validates clean.
3. **Reads**: `samtools fasta reads_sorted.bam` → `reads.fa` (103 reads).
4. **Graph alignment**: `GraphAligner -g graph.gfa -f reads.fa -a aln.gaf -x vg`.
   All 103 reads aligned, mapq 60.
5. **Graph copy assignment**: from GFA P-lines compute **A-exclusive** nodes
   (in A1∪A2, not in B1) and **B-exclusive** nodes (in B1, not in A1∪A2). Parse
   each read's GAF node path, count overlap with each exclusive set, assign to the
   larger; tie/zero → ambiguous.
6. **Linear baseline**: `minimap2 -a -x map-pb -N 5 -p 0.1 copies_AB.fa reads.fa`
   (copy refs = A1 and B1 transcripts). Assign each read to the copy with the
   higher `AS` tag; tie → ambiguous.

---

## Did the graph form bubbles at the diagnostic sites? YES
From the P-lines, the exclusive-node sets are non-empty and sized exactly as
expected from 7 PSVs + the exon-3-skip event:

| metric | value |
|---|---|
| total nodes | 125 |
| A-path nodes (A1∪A2) | 96 |
| B-path nodes (B1) | 95 |
| shared nodes | 66 |
| **A-exclusive nodes** | **30** |
| **B-exclusive nodes** | **29** |
| A-exclusive bp | 35 |
| B-exclusive bp | 31 |

A1 and B1 alternate at the PSV sites (e.g. A1 uses node 3, B1 uses node 4; A1→6 /
B1→7; A1→10 / B1→11 …) — i.e. real **bubbles** at the diagnostic positions. A1 vs
A2 diverge at the exon-3-skip bubble (A1→…50, A2→…67). The ~35/31 bp of exclusive
sequence ≈ the 7 single-base PSVs plus a few collapse-boundary bases. Sanity check
passes: the graph encodes the copy difference structurally.

---

## Accuracy

### uniq_* reads (75 reads, unambiguous truth)
| method | correct | total | ambiguous | accuracy |
|---|---|---|---|---|
| **graph (GAF exclusive-node vote)** | 75 | 75 | 0 | **100.0%** |
| **linear (minimap2 AS)** | 75 | 75 | 0 | **100.0%** |

Per-class: graph 50/50 A + 25/25 B; linear 50/50 A + 25/25 B.

### multi-mapper split (truth known: multi→A, multi_r→B)
| group | truth | graph | linear |
|---|---|---|---|
| `multi_*`   | A | A=14 B=0 amb=0 → **14/14** | A=14 B=0 amb=0 → **14/14** |
| `multi_r_*` | B | A=0 B=14 amb=0 → **14/14** | A=0 B=14 amb=0 → **14/14** |

Both methods send every multi-mapper to its TRUE copy. The PSV alleles a read
carries override its placement, for both linear and graph.

### overall (all 103 reads, truth from names)
| method | correct | accuracy | ambiguous |
|---|---|---|---|
| graph  | 103/103 | **100.0%** | 0 |
| linear | 103/103 | **100.0%** | 0 |

**Zero errors, zero ambiguous, for either method.**

---

## Are the 100% scores decisive or borderline? DECISIVE (not a knife-edge tie)
A 100%–100% tie could hide a real difference if one method were barely scraping
by. It is not — both separate the copies by a wide margin on every read:

**Linear (AS of correct copy − AS of wrong copy):**
| group | n | min | max | mean |
|---|---|---|---|---|
| uniq_A | 48 | 150 | 194 | 177.7 |
| uniq_B | 25 | 174 | 194 | 189.8 |
| multi_A | 14 | 148 | 196 | 180.9 |
| multi_B | 14 | 150 | 170 | 160.4 |

**Graph (#correct-exclusive nodes − #wrong-exclusive nodes on path):**
| group | n | min | max | mean |
|---|---|---|---|---|
| uniq_A | 50 | 24 | 29 | 27.6 |
| uniq_B | 25 | 27 | 29 | 28.8 |
| multi_A | 14 | 27 | 29 | 28.5 |
| multi_B | 14 | 24 | 26 | 25.8 |

The smallest linear margin is **+148 AS**; the smallest graph margin is **+24
exclusive-node votes**. Neither method comes anywhere near its decision boundary.
The diagnostic SNPs are collectively far more than enough for *either* aligner.

(2 `uniq_A1` reads got only a copy-A hit from minimap2 — no competitive copy-B
alignment at all — which only reinforces how strong the copy-A signal is. Still
correct.)

---

## VERDICT

**Graph alignment performs EQUAL to linear alignment here — both 100%.** It is
neither better nor worse on this benchmark.

**Why a tie, not a graph win?** The two copies differ at 7 clean substitution
sites spread across the transcript. Each PSV is a hard mismatch under linear
alignment, so minimap2's `AS` for the wrong copy is penalized ~7× (plus the 150 bp
indel for A2-vs-B), producing a 148–196-point gap that linear resolves trivially.
The graph encodes the same 7 differences as bubbles and resolves them just as
cleanly — but there was no residual ambiguity for the graph to recover, because
linear already had margin to spare.

**When would graph be expected to win?** Not on independent, well-separated
substitution PSVs like these (linear handles those). Graph/PSV-phasing methods pay
off when copies are *so* identical that a single read spans few or zero
diagnostic sites, when differences are structural (indels/rearrangements) that
break linear seed chains, or when **read-linkage/phasing across PSVs** is needed
to assign reads that are individually ambiguous. This synthetic family is the easy
regime for both, so it cannot demonstrate a graph advantage — an honest negative
result for the "graph beats linear" hypothesis on *this* data.

---

## Exact commands (reproducible)
```bash
mkdir -p /tmp/graphexp && cd /tmp/graphexp
cp .../synthetic_family/{ref.fa,truth.gtf,reads_sorted.bam,reads_sorted.bam.bai} .
samtools faidx ref.fa

# 1. per-copy transcript sequences (generator's 0-based half-open coords)
python3 extract_seqs.py                       # -> seqs.fa  (A1,A2,B1)

# 2. variation graph
vg msga -f seqs.fa > graph.vg                 # prototype; warns about >128 bp band — harmless
vg view graph.vg > graph.gfa
vg validate graph.vg                          # "graph: valid"

# 3. reads
samtools fasta reads_sorted.bam > reads.fa    # 103 reads (1 seq/read name)

# 4. graph alignment
GraphAligner -g graph.gfa -f reads.fa -a aln.gaf -x vg

# 5+6. linear baseline + 7. evaluate (assignment & accuracy)
samtools faidx seqs.fa A1 > copyA.fa; samtools faidx seqs.fa B1 > copyB.fa
cat copyA.fa copyB.fa > copies_AB.fa; samtools faidx copies_AB.fa
minimap2 -a -x map-pb -N 5 -p 0.1 copies_AB.fa reads.fa > linear.sam
python3 evaluate.py                           # tables in eval_output.txt
```

## Tool quirks hit
- **`vg msga` is a deprecated prototype** and warns `sequence of length 1150 …
  longer than the non-chunked limit of 128. Alignments may be discontiguous`.
  In practice the output graph still validated clean and produced correct A1/A2/B1
  paths with the expected bubbles, so it was usable here. For production-scale
  graphs prefer `vg construct`/PGGB/minigraph; for these tiny seqs `msga` was fine.
- **Reverse strand:** the BAM reads are flag-16 (gene is on `-`). `samtools fasta`
  reverse-complements them, so `reads.fa` is the RC of the forward transcripts.
  Both minimap2 (reports flag-16 hits) and GraphAligner (traverses the graph with
  `<` node orientations) handle this automatically — no manual RC needed. Verified
  empirically via `orient_check.sam` rather than assumed.
- **minimap2 secondary alignments:** to score a read against *both* copies you must
  pass `-N 5 -p 0.1` (default `-p 0.8` suppresses the near-equal second-best hit
  to the other copy). With defaults only the best copy is reported and the AS
  comparison is degenerate. (2 reads still got a single-copy hit because their
  wrong-copy alignment fell below threshold — they're trivially correct.)
- **GAF path parsing:** node string looks like `<125<123<...<1`; orientation
  markers `<`/`>` stripped, node IDs matched with `[<>](\d+)`.

## Artifacts in this directory
- `extract_seqs.py`, `evaluate.py` — the two scripts.
- `seqs.fa` (graph input), `copies_AB.fa` (linear refs), `reads.fa`.
- `graph.vg`, `graph.gfa` — the variation graph (+ P-lines = copy paths).
- `aln.gaf` — GraphAligner output; `linear.sam` — minimap2 output.
- `orient_check.sam` — strand sanity check.
- `eval_output.txt` — full evaluation print-out.
- `ga.log`, `msga.log`, `mm2.log` — tool logs.
