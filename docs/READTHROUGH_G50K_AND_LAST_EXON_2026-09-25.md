# Readthroughs: minimap2 `-G 50k` vs a "reaches past the last exon" rule

## §0 Pre-registration (2026-09-25, written before any number below was computed)

**Advisor:** the readthroughs are pipeline sloppiness; use minimap2's levers, specifically `-G 50k` (max intron
50 kb; the `splice` preset default is 200 kb). **User:** that cuts real genes with >50 kb introns; prefer a rule
that removes a readthrough when it reaches past the last exon of a gene.

**Instrument** `bench/mechanism/readthrough_rules.py`: primary reads (`-F 2308`), one contig at a time. Junction
classes from the annotation (never from the reads):
- **RT** (readthrough junction): donor base in the exon union of gene A, acceptor base in gene B, A ≠ B, same
  strand, A and B spans do not overlap, and the junction is not an annotated intron of any gene. RefSeq records
  whose description says `readthrough` are removed from the gene set first, so they cannot hide an RT.
- **ANN**: junction exactly equal to an annotated intron of a (non-readthrough) gene; split `<50 kb` / `≥50 kb`.
Only junctions with ≥2 reads are scored.

**Q1 (`-G 50k`, descriptive + a direct realignment):** (a) fraction of RT junctions and RT reads with intron
> 50 kb (what `-G 50k` can remove at most); (b) ANN introns ≥ 50 kb with reads: junctions, reads, genes (what it
costs). (c) Realign on chr20 the primary reads that carry either an RT junction or an ANN ≥ 50 kb junction with
`minimap2 -ax splice:hq -uf` at default `-G` and at `-G 50k`, and count how many keep the junction.
**Prediction:** most RT introns are short (PKD1P6-NPIPP1's fusion intron is 6.6 kb), so `-G 50k` removes a
minority of RT junctions while breaking every expressed ≥50 kb real intron.

**Q2 (annotation-free "past the last exon" rule).** The last exon of gene A is where its mature transcripts end
(the polyA site). A readthrough splice has A's polyA site **inside its intron** (it splices out of or over A's
terminal exon); a real intron rarely does. Feature per junction J: `T` = spliced primary reads (≥1 intron) whose
3′ end (strand-aware) falls strictly inside J's intron, `S` = reads using J; `R = T / (T + S)`.
Rule: J is a readthrough iff `R ≥ θ`. DEV = human chr16 (θ picked there: largest θ with RT recall ≥ 0.5,
reported with the full ROC); HELD-OUT = human chr20 and gorilla NC_073244.2 (no re-tuning, never pooled).
**Bar:** held-out AUC ≥ 0.80 AND at θ, ANN ≥ 50 kb false-removal rate below `-G 50k`'s (which is 100% of
those junctions by construction) and ANN overall FPR ≤ 5%. Otherwise the rule is refuted.

### §0b Refinement on DEV, frozen before held-out
The first feature (`R_any`: any spliced read's 3′ end inside the intron) had dev AUC .962 but counted reads of
genes NESTED inside long introns. Refined once, on dev only: `R_up` counts only reads whose 5′ end lies
UPSTREAM of the donor (transcripts of the donor's own gene that terminate inside the intron). Dev AUC .961
(vs ANN ≥ 50 kb .821 → .890). Pre-registered θ (largest θ with dev RT recall ≥ 0.5): **θ = 0.9524**
(dev: recall .500, ANN FPR .0044, precision .716). Frozen for chr20 and gorilla.

## §1 Verdict

- ⛔ **`-G 50k` is the wrong lever.** Readthrough introns are SHORT (median 17.8 / 18.0 / 18.9 kb; only
  12-17% exceed 50 kb). In a direct realignment it keeps **100%** of the ≤50 kb readthrough read-junctions and
  breaks **99.6%** of the real ≥50 kb intron read-junctions, **92% of them into supplementary split
  alignments**. That creates exactly the chimeric-looking records it was meant to remove.
- ⭐ **The "reaches past the last exon" rule works and passes the pre-registered bar on both held-out
  substrates**, with no annotation. It removes about half of the readthrough junctions and ~0.1-0.5% of real
  junctions, and it never breaks a real long intron.

## §2 Q1: what `-G 50k` removes vs costs

| | human chr16 (dev) | human chr20 | gorilla NC_073244.2 |
|---|---|---|---|
| readthrough (RT) junctions / reads | 242 / 2,368 | 123 / 1,146 | 103 / 612 |
| RT intron length median (IQR) | 17.8 kb (9.3-35.3) | 18.0 kb (10.4-37.2) | 18.9 kb (8.8-38.9) |
| **RT junctions `-G 50k` can remove** | 29 (12.0%), 268 reads | 16 (13.0%), 64 reads | 17 (16.5%), 91 reads |
| **real annotated introns ≥ 50 kb it breaks** | **135 junctions, 25,240 reads, 68 genes** | **94, 7,455 reads, 53 genes** | **43, 2,248 reads, 25 genes** |

Direct realignment (chr20, the 8,437 primary reads carrying an RT or a ≥50 kb real junction, `minimap2 -ax
splice:hq -uf`, chr20 only):

| read-junctions | kept at default `-G 200k` | kept at `-G 50k` | lost at 50k → supplementary split |
|---|---|---|---|
| RT, intron ≤ 50 kb (1,082) | 1.000 | **1.000** | 0 |
| RT, intron > 50 kb (64) | 1.000 | 0.016 | 61 |
| real annotated, ≥ 50 kb (7,455) | 1.000 | **0.004** | 6,829 |

## §3 Q2: the rule

**Rule (annotation-free):** for a splice junction J, count `S` = reads using J and `U` = spliced reads whose 5′ end
is upstream of J's donor and whose 3′ end falls inside J's intron, i.e. transcripts of the same gene that
polyadenylate where J keeps going. Remove J iff `U / (U + S) ≥ 0.95`. In words: **drop a splice if at least 95%
of the gene's transcripts that reach its donor end before its acceptor.** This is "reaches further than the last
exon", with the last exon read off the reads' own polyA sites.

| (θ = 0.9524 frozen from dev) | AUC (vs all real / vs real ≥50 kb) | RT recall | real-junction FPR | real ≥ 50 kb removed | precision |
|---|---|---|---|---|---|
| human chr16 (dev) | .961 / .890 | .500 | .0044 | 2 / 135 | .716 |
| **human chr20 (held-out)** | **.960 / .895** | **.520** | **.0039** | **1 / 94** | **.719** |
| **gorilla (held-out)** | **.964 / .937** | **.417** | **.0012** | **0 / 43** | **.782** |

Bar (AUC ≥ .80, FPR ≤ 5%, fewer long real introns lost than `-G 50k`): ⭐ **passes on both held-out substrates**.
Read-weighted at θ: removes 566 / 459 / 145 RT reads (24-40%) against 267 / 100 / 55 reads on real junctions,
vs `-G 50k`'s 268 / 64 / 91 RT reads against 25,240 / 7,455 / 2,248 real ones.

**θ sweep (descriptive).** θ is a ratio, and it moves recall against precision smoothly: 0.80 → recall .68-.72,
FPR .005-.020; 0.90 → .54-.59, FPR .003-.009; 0.99 → .07-.29, FPR < .001. There is no cliff to tune.

**What it removes and what it keeps (and why that is the right split).**
- **Removed:** minor run-on past a strong polyA site. Example: ARL6IP1→RPS15A, where 57 reads splice through but
  12,201 of ARL6IP1's transcripts end first (R = .995).
- **Kept:** dominant, constitutive fusions whose donor is not a strong terminator, e.g. SLX1B→SULT1A4 (120
  through, 10 end), BOLA2→SMG1P6 and SMG1P2→SMG1P6. SLX1-SULT1A and BOLA2-SMG1P6 are curated RefSeq readthrough
  genes, and PKD1P6-NPIPP1 (110 MAPQ-60 reads) is the same class. This is the biology register 967/968 kept
  hitting: a well-supported fusion is a real molecule, so a rule that removed it would be wrong.
- **False removals are alternative last exons** (CLEC16A, IFT140, BBS2): a gene whose short isoform terminates
  inside the intron that its long isoform splices over. Topologically that IS a readthrough; only the
  annotation says the downstream exons belong to the same gene.

## §4 Caveats / not yet measured

- Junction-level only. The family-level effect is NOT measured here. Register r1017's oracle (a perfect split of
  every fused chr16 locus) bounds the gain at +0.021 referee F with NPIP worse, so this rule cleans up transcript
  models and loci more than it moves family F.
- Positives are defined geometrically (donor in gene A, acceptor in non-overlapping gene B). Some are unannotated
  isoform extensions into lncRNA/LOC models rather than readthroughs in the strict sense.
- Needs full-length reads with trustworthy 3′ ends (Kinnex/Iso-Seq; 3′ dispersion ±5-25 bp, §6w3). Not for
  short-read or 3′-truncated libraries.
- Not implemented in the assembler. Natural place: a junction filter before pass-1 in `denovo_assemble.rs`,
  default OFF, with the family-level scoring as its own prereg.

## Reproduce
```bash
O=/mnt/linuxdisk/tmp/readthrough
python3 bench/mechanism/readthrough_rules.py hsa16 /mnt/linuxdisk/tmp/advisor_jaccard/chr16.bam /mnt/linuxdisk/tmp/sedef/chr16.genes.gff chr16 $O
python3 bench/mechanism/readthrough_rules.py hsa20 /mnt/linuxdisk/tmp/ppar/chr20.bam /mnt/linuxdisk/tmp/sedef/chr20.genes.gff chr20 $O
python3 bench/mechanism/readthrough_rules.py ggo44 /mnt/linuxdisk/tmp/gw22/sec/ggo44.bam /mnt/linuxdisk/tmp/gw22/sec/ref/NC_073244.2.genes.gff NC_073244.2 $O
# -G realignment: $O/hsa20.sel.fa (primaries carrying RT / >=50 kb junctions) vs $O/chr20.fa at -G 200k and 50k -> $O/sel_G{200k,50k}.bam
```
`<tag>.junctions.tsv`: class, intron, reads, ends inside, `R` (= `R_up`), gene pair.

## §5 Pre-registration: 3′-end DENSITY to find the real ends (2026-09-25, before computing)

User: *"could we use density to infer where the real ends are?"* Current rule counts EVERY upstream-starting
spliced read ending inside the intron (`R_up`). Density version: cluster spliced reads' 3′ ends per strand (a new
cluster when the gap to the previous end exceeds 25 bp). A cluster is a **real end** iff it holds ≥ 3 reads AND
is not internal priming (the 20 genomic bases downstream of the cluster mode, strand-aware, are < 60% A, the
SQANTI3 criterion). `R_peak` = the same ratio counting only ends that fall in real-end clusters. θ is picked on
dev chr16 exactly as before (largest θ with RT recall ≥ .5) and frozen for chr20 and gorilla.
**Bar:** density helps iff on BOTH held-out substrates, at the recall the current rule reaches, the real-junction
FPR is ≥ 25% lower (relative), and AUC does not drop by > .01. Descriptive: canonical PAS (AATAAA/ATTAAA within
50 bp upstream) rate of the dominant in-intron peak for RT vs falsely removed real junctions. If the false removals
are alternative last exons (a real end), their PAS rate should match the readthroughs', and density cannot remove
them. **Prediction:** a small gain on long introns (scattered degraded ends drop out), no gain on the
alternative-last-exon false positives.

## §6 Result: density does NOT improve the rule ⛔ (it confirms what the false removals are)

θ frozen on dev: `R` 0.9524, `R_peak` 0.918.

| | AUC `R` → `R_peak` | FPR at matched recall, `R` → `R_peak` | `R_peak` at its own θ: recall / FPR / prec |
|---|---|---|---|
| chr16 (dev) | .961 → .897 | .0044 → .0040 (−8%) | .500 / .0040 / .733 |
| chr20 (held-out) | .960 → **.894** | .0039 → .0028 (**−28%**) | .569 / .0052 / .680 |
| gorilla (held-out) | .964 → **.924** | .0012 → .0010 (**−17%**) | .456 / .0017 / .734 |

Bar (≥25% lower FPR on BOTH held-out AND AUC drop ≤ .01): ⛔ **fails twice**. Gorilla gains only 17%, and AUC
drops .04-.07 everywhere. The AUC loss is readthroughs with sparse end evidence: fewer than 3 ends in one peak
gives `R_peak` = 0, so thin but correct signal is thrown away.

**Why density cannot help: the false removals ARE real ends.** A canonical polyA signal (AATAAA/ATTAAA ≤ 50 bp
upstream) sits at the dominant in-intron peak of **73% / 76% / 82%** of the falsely removed real junctions,
vs **81% / 59% / 76%** of the flagged readthroughs. Both are genuine, dense, PAS-marked polyA sites. The false
positives are alternative last exons, where a real polyA site sits inside an intron that a longer isoform of the
SAME gene splices over. No end-based signal (count, density or PAS) separates "the next exons belong to another
gene" from "the next exons belong to this gene". Only the annotation, or homology of the downstream exons, can.
⟹ Keep `R` (all ends). Density is a sound way to call polyA sites (PAS-validated) but adds nothing to this rule.

## §7 Pre-registration: 5′-START density (an independent promoter for the downstream gene) + canonical splicing

(2026-09-25, written before computing.) The remaining false removals are alternative last exons (§6). A true
readthrough goes into a gene B that has **its own promoter**: B's own transcripts start inside J's intron (at B's
first exon, which the readthrough usually skips) and run on through J's acceptor exon. In an alternative last
exon, the downstream exons are reached only from A's promoter.

**Feature.** `V` = spliced primary reads whose 5′ end lies inside J's intron and whose 3′ end lies beyond J's
acceptor (strand-aware), i.e. transcripts from an internal promoter that reach the downstream exons without using
J. **Density version (primary):** cluster spliced reads' 5′ ends per strand (new cluster when the gap exceeds
100 bp, since 5′ ends are dispersed ±150-250 bp, §6w3); a **real start** = a cluster of ≥ 3 reads. `V_peak` counts
only reads whose 5′ end is in a real start. Ratio `Q = V_peak / (V_peak + S)`.

**Canonical splicing (user):** every junction is motif-checked (GT-AG, GC-AG, AT-AC, strand-aware) and **only
canonical junctions are scored**, which is also what `--assemble-only`'s default `--assembly-junctions strict`
admits. The canonical fraction of RT vs real junctions is reported.

**Rule:** remove J iff `R ≥ 0.9524` (frozen, §3) AND `Q ≥ θ_Q`. θ_Q is picked on dev chr16 as the largest value
keeping ≥ 80% of the RT junctions that `R` alone removes. Frozen for chr20 and gorilla.
**Bar (both held-out):** real-junction FPR ≥ 50% lower (relative) than `R` alone, with RT recall ≤ 0.10 lower
(absolute). **Prediction:** it works on multi-gene readthroughs into expressed genes, and fails where B is a
lncRNA/LOC model with little own expression, or where an alternative last exon also has an internal promoter.

## §8 Result: 5′-start density — ⚠ a real precision gain that fails the recall clause; canonical splicing changes nothing

**Canonical splicing.** Readthrough junctions are canonical at **98.3% / 100% / 97.1%** (chr16 / chr20 /
gorilla), against 99.9-100% for real introns. Enforcing GT-AG/GC-AG/AT-AC removes 4, 0 and 3 readthrough junctions.
**Readthroughs are properly spliced**, which is more evidence that they are real processing and not alignment
artifacts. All numbers below are on canonical junctions only.

θ_Q frozen on dev = **0.6667** (keeps ≥ 80% of the 120 canonical RT junctions that `R` alone removes).

| | `R` alone: recall / FPR / prec | `R & Q`: recall / FPR / prec | ΔFPR (rel) | Δrecall | bar |
|---|---|---|---|---|---|
| chr16 (dev) | .504 / .0044 / .714 | .403 / .0017 / .835 | −61% | −.101 | — |
| **chr20 (held-out)** | .520 / .0040 / .719 | .374 / **.0011** / **.868** | **−72%** ✓ | **−.146** ✗ | ⛔ |
| **gorilla (held-out)** | .420 / .0012 / .778 | .250 / **.0001** / **.962** | **−92%** ✓ | **−.170** ✗ | ⛔ |

Among junctions `R` flags, the independent-promoter share `Q` separates cleanly: median **.96 / .89 / .85 for
readthroughs vs .31 / .33 / .00 for real junctions**. The prediction held. **Verdict per the bar: not adopted as
specified.** The false-removal half passed with room to spare (FPR −72% / −92%, 7 and 1 remaining false removals)
but recall fell by more than .10 on both held-out substrates.

**The recall loss is not a threshold choice** (descriptive sweep, not used for the verdict): θ_Q from 0.2 to 0.67
gives chr20 recall .40 → .37 and gorilla .25 flat. The readthroughs it gives up are ones where **gene B has no
transcripts of its own in this library**: HAGH→IGFALS, CASC16→TOX3, PTPRA→GNRH2, MAG→CD22, LPAR2→PBX4 (liver /
neuronal / immune genes in testis) and lncRNA/LOC models (PSMF1→LOC…, SMOX→LINC01433). With no independent
promoter reads there, the RNA cannot tell "readthrough into a silent gene" from "an unannotated extension of A".
Keeping them is the conservative error.

**Recommendation.** Use `R & Q` where precision matters (default-safe filter: ~1 real junction in 1,000-10,000
removed, precision .87-.96) and `R` alone where recall matters. Either way it is one opt-in knob with two
dimensionless ratios, not `-G`. A third arm with a relaxed recall clause would be a new prereg, not a re-reading of
this one.

## §9 Pre-registration: do per-read BAM tags mark readthrough reads? (2026-09-25, before computing)

User: do readthrough reads have too-high `NM` / `de` or other tags? Readthrough reads = primaries carrying an RT
junction. Controls = primaries whose every junction is an annotated intron (random 10% sample). Features:
`de`, `NM`/aligned length, `AS`/aligned length, `s2/s1` (minimap2's second-best/best chain score), `cm`
(minimizers on the chain), MAPQ, total soft-clip, and two junction-local features for the RT junction (or a random
junction for controls): the shorter flanking exon anchor, and mismatches+indels within 20 aligned bases of the
junction on either side. Scored at read level (AUC, direction-free = max(AUC, 1−AUC)) and at junction level
(median of the junction's reads, RT vs ANN junctions). DEV chr16, HELD-OUT chr20 + gorilla.
**Bar:** a tag is useful iff its junction-level AUC is ≥ .70 on BOTH held-out substrates. **Prediction:** none
passes. Readthroughs are canonical, MAPQ-60, correctly aligned molecules, so their alignment quality should equal
normal reads'. Junction-local anchor/mismatch features have the best chance, if some readthrough splices are
aligner-forced.

## §10 Result: no BAM tag marks readthrough reads ⛔ (`bench/mechanism/readthrough_tags.py`)

Junction-level AUC (direction-free), chr16 dev / **chr20** / **gorilla**. RT reads 2,349 / 1,135 / 604 vs control
reads 34,995 / 20,068 / 11,886:

| feature | RT vs control median (chr20) | AUC chr16 | **AUC chr20** | **AUC gorilla** |
|---|---|---|---|---|
| `de` | .0020 vs .0017 | .567 | .558 | .589 |
| `NM` / aligned bp | .0021 vs .0019 | .561 | .551 | .596 |
| `AS` / aligned bp | .880 vs .893 | .584 | .534 | .614 (reversed: RT higher) |
| `s2/s1` | .016 vs 0 | .612 | .560 | .670 |
| `cm` | 671 vs 659 | .505 | .535 | .669 (RT reads are longer) |
| MAPQ | 60 vs 60 | .541 | .500 | .504 |
| soft clip | 1 vs 1 bp | .536 | .592 | .571 |
| shorter anchor at the junction | 117 vs 96 bp | .547 | .538 | .527 |
| errors within 20 bp of the junction | 0 vs 0 | .510 | .502 | .547 |

⛔ **None reaches .70 on both held-out substrates** (best: `s2/s1` .56/.67). The prediction held. Readthrough reads
align as cleanly as normal reads: divergence differs by 0.03 percentage points, MAPQ is 60, anchors are ~100 bp and
there are no errors next to the splice. The alignment is correct. Only locus-level context (where A's transcripts
end, whether B has its own promoter) separates them.

## §11 Pre-registration: combining the tags (2026-09-25, before computing)

User: the advisor wants the de novo route, so can the tags be combined? Note: `R` and `Q` (§3/§8) are computed
from reads only, with no annotation, so they are already de novo-compatible. The annotation is used only to LABEL
junctions for scoring.
**Model:** L2 logistic regression (class-balanced, features standardised on dev) on junction-level medians of the 9
tags of §9 (+ log read count), trained on chr16, frozen, applied to chr20 and gorilla. Two models: **T** = tags
only; **T+RQ** = tags + `R` + `Q`, compared with **RQ** = `R` + `Q` in the same logistic form.
**Bars:** T is useful iff held-out AUC ≥ .70 on BOTH. The tags add to the rule iff T+RQ beats RQ by ≥ .01 AUC on
BOTH held-out substrates. **Prediction:** T < .70 (no single tag carries signal and they are correlated: `de`/`NM`/`AS`
all measure divergence); T+RQ ≈ RQ.

## §12 Result: combining the tags adds nothing ⛔; the apparent win was read DEPTH (`bench/mechanism/readthrough_combined.py`)

Canonical junctions; RT 238 / 121 / 100 vs control junctions 5,245 / 3,023 / 3,414 (chr16 / chr20 / gorilla).

| model (trained chr16, frozen) | AUC chr16 | **AUC chr20** | **AUC gorilla** |
|---|---|---|---|
| T = 9 tags + log read count (as pre-registered) | .988 | .980 | .966 |
| ⚠ read count ALONE | .978 | .973 | .981 |
| **9 tags, no count** | .639 | **.527** | **.480** |
| 9 tags, controls depth-matched to RT (log2 bins, 5:1) | .835 | .773 | **.611** |
| R + Q (the locus rule) | .974 | .995 | .980 |
| R + Q, depth-matched | .958 | .968 | .971 |
| T + RQ | .995 | .996 | .989 |

- The pre-registered T "passes" (.98 / .97), but **that is the read count** (largest coefficient, −3.54; count alone
  gives the same AUC). ⚠ The control junctions are also depth-biased, because reads were sampled at 10% and so
  well-read junctions dominate. The tags themselves carry nothing: chance without the count (.53 / .48), and
  depth-matched they miss the .70 bar on gorilla (.61). **⛔ Tags: refuted.**
- Tags on top of the rule: T+RQ − RQ = +.001 (chr20) / +.009 (gorilla) AUC, under the +.01 bar on both, and that
  margin is carried by the count. **⛔ Nothing to add.**
- **Read count as a rule, measured on the full junction table** (unbiased): AUC .862 / .852 / .842. Readthrough
  junctions are rare (median 4 / 4 / 3 reads vs 68 / 63 / 28 for real introns), but a floor that catches half the
  readthroughs (≤ 4 reads) removes **12-14% of real junctions** (1,372 / 780 / 1,378). The locus rule removes 48 /
  25 / 12 at the same recall. The r1067-r1073 floor study already closed this lever: low-count junctions are mostly
  real, lowly expressed introns.
- **The locus rule is already de novo.** `R` and `Q` are computed from read 3′ and 5′ ends only, with no
  annotation, and depth-matched it keeps AUC .96-.97.

## §13 Pre-registration: enforcing B's OWN FIRST EXON (2026-09-25, before computing)

User: how can we enforce that the downstream gene has its own promoter? `Q` (§8) counts reads that start inside J's
intron in a ≥3-read start cluster and end beyond J's acceptor, but an UNSPLICED intronic read (pre-mRNA, internal
priming) also passes. **`Q1` (strict):** additionally, the read's FIRST intron (in transcript orientation) must have
its donor inside J's intron. The read then has a first exon of its own inside the region the readthrough skips and
splices out of it into the downstream exons, i.e. a transcript started at B's promoter.
Rule `R ≥ .9524 AND Q1 ≥ θ`, θ picked on dev exactly as in §8 (keep ≥ 80% of the RT that `R` removes), frozen.
**Bar vs `R & Q` on BOTH held-out:** FPR not higher, and recall ≥ `R & Q` − .02 (the same promoter evidence, stated
cleanly), OR FPR ≥ 25% lower at recall within .05. **Prediction:** nearly identical. Most start-cluster reads
that run past the acceptor are already spliced.

## §14 Result: B's own first exon (`Q1`) — ⭐ passes as the clean statement of the same evidence

θ frozen on dev: `Q` .6667, `Q1` .625. Canonical junctions.

| | `R & Q`: recall / FPR / prec | `R & Q1`: recall / FPR / prec |
|---|---|---|
| chr16 (dev) | .403 / .0017 (19) / .835 | .403 / **.0013 (14)** / **.873** |
| **chr20 (held-out)** | .374 / .0011 (7) / .868 | .358 / .0011 (7) / .863 |
| **gorilla (held-out)** | .250 / .0001 (1) / .962 | .250 / .0001 (1) / .962 |

Bar (FPR not higher, recall within .02 on both held-out): ⭐ **passes**. As predicted it is nearly identical, since
the reads that pass `Q` were already spliced out of a first exon of their own, but `Q1` is the definition to
state: *B has its own promoter iff ≥ 3-read clusters of transcript starts lie in the skipped region, and those
transcripts splice out of their own first exon into the downstream exons.* The remaining false removals (chr20:
TASP1, RRBP1, GINS1, RALGAPB, SYCP2, MTG2, TCFL5; gorilla: 1) are genes where an alternative last exon coincides
with an **internal alternative promoter** (an alternative first exon in the same intron). That is structurally a
second transcription unit, and RNA alone cannot call it anything else. External evidence (CAGE peaks, CpG islands
at the start cluster) could arbitrate; not tested.
