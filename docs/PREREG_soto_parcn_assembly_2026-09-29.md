# Pre-registration: exact assembly k-mer parCN for Soto 2025's SD98 genes, and their human-vs-ape calls

**Written 2026-09-29 (KEY=soto_parcn_asm) BEFORE any assembly k-mer count, parCN or famCN value was computed.**
Reference CHM13 v2.0 (human). Human (CHM13, HG002) and ape (chimpanzee, gorilla, orangutan) numbers are reported **per
genome and never pooled**. Nothing in `src/` is edited; nothing is committed. Scratch (code, k-mer tables, counts):
`/mnt/linuxdisk/tmp/rustle_figures_dev/soto_parcn_asm/`. Inputs: `bench/soto/soto_parCN_S1E.tsv` (sha1 0c4cb730,
Table S1E, 1,833 rows) and `bench/soto/soto_famCN_S1C.tsv` (sha1 d008a179, Table S1C). This file binds once
Amendment 1 records its sha1 and the instrument sha1s.

## 0. The question

Soto's paralog-specific copy number (parCN, S1E "Median parCN") is QuicK-mer2 run on 2,504 1KGP short-read genomes:
k-mers that occur once in the reference, counted in reads, GC-normalised. With complete assemblies the same k-mers can
be counted **exactly**. (Q1) Does an exact assembly count of QuicK-mer2's own k-mer set in HG002 reproduce Soto's
population parCN where the gene's CN is common? (Q2) Do the ape assemblies reproduce Soto's human-vs-great-ape calls
(S1C "famCN Status per Paralog", S1E "Family status")? (Q3) What do the named families (NPIP, TBC1D3, SRGAP2, ARHGAP11,
FAM72, NOTCH2NL) look like row by row?

## 1. What Soto and QuicK-mer2 did (established today from files, before this prereg)

- **Paper** (`ThesisVault/Fulltext/Soto_2025_brain_evolution_fulltext.md` l.1051): parCN = QuicK-mer2 on 1KGP 30×
  + 4 archaics, reference **T2T-CHM13 v1.0**; "genotyped across SD98 regions overlapping protein-encoding and
  unprocessed pseudogenes by calculating the **mean** parCN across the region of interest for each sample". S1E
  counts imply n = 2,504 (e.g. OR4F21 parCN>0.5 = 2,476 = 98.88%). **"Median parCN" = median over samples.**
- **Soto's QuicK-mer2 config** (their released `section_I&II/QM2_pipeline/readme.md` + `quickmer2.smk`, local copy in
  `rustle_figures_dev/soto_quickmer2/HSD_brain_evolution/`): `quicKmer2 search -c control.bed ref.fasta -e 1 -w 500`
  (`distance=1 window=500`), default `-k 30`, default `-d 100`.
- **QuicK-mer2 source** (`tools/QuicK-mer2/QuicKmer.c`): k = 30, **canonical** (min of k-mer and reverse complement);
  **case-insensitive** (the 2-bit code is `(ascii >> 1) & 3`, so soft-masked bases are counted like upper case; only
  `N` resets); the all-A k-mer (code 0) is skipped; occurrence counts cap at 255. A k-mer is **kept iff it occurs
  exactly once in the reference AND its "edit depth" < 100**, where edit depth = sum of the reference occurrence counts
  of its 90 Hamming-1 substitution variants (canonicalised) for `-e 1` (l.687-709, l.1222: dropped if
  `occr > 1 || edit_depth >= 100`). No RepeatMasker mask is used. Windows = 500 kept k-mers; window CN = GC-corrected
  mean k-mer depth / (control mean / 2).
- **Region** = S1E `SD98_v1.0` (gene ∩ SD98, v1.0) with its liftover `SD98_v2.0` (313 rows differ from `Gene Coords`).
  `SD98_v2.0` is BED half-open: PRAMEF2 `chr1:12401097-12405915` vs v2.0 RefSeq gene start 12,401,098 (1-based).
  1,831 rows have v2.0 coordinates (2 "Unmapped", both AP000522.1 on acrocentric short arms, excluded).
- **Human-vs-ape** (paper l.1047, not QuicK-mer2): WSSD famCN (read depth over the repeat-masked reference, all copies
  of a family). Per paralog: excluded if human median famCN > 10 ("Undetermined", 937 S1E rows); "Duplicated in
  humans" (121) / "Expanded in humans" (136) if human median (SGDP n = 269) > max over Clint, Mhudiblu, Kamilah,
  Susie, split at great-ape max 2.5; otherwise labelled **"Expanded in great apes"** (639 rows; the label covers every
  non-excluded paralog not called human; 39 of the 607 S1C rows so labelled have human > ape max — a Soto internal
  inconsistency, reported, not scored away). Family status "Human duplicated gene family" = ≥ 1 member
  duplicated/expanded OR ≥ 1 member non-syntenic with chimpanzee: 251 families, of which 124 have a CN-called member
  and 127 do not (synteny-only or excluded members).
- **parCN classes** (from the S1E columns): Fixed = parCN = 2 in ≥ 98% of 1KGP (340 rows, median of medians 2.01);
  Nearly-Fixed = parCN > 0.5 in ≥ 98% (978); Polymorphic = rest (515).

## 2. Instrument (fixed now)

**Reference for the k-mer set:** CHM13 v2.0 `chm13v2.0.fa` **without chrY** (= v1.1 content; Soto's v1.0 has no Y;
v1.0 → v1.1 changed the rDNA models on the acrocentric short arms and a few patches — disclosed, see §6).

**Query set Q:** every 30-mer (ACGT only, case-insensitive) fully inside a region: the 1,831 S1E `SD98_v2.0` regions +
the control regions (§2.4). Stored canonical (A0 C1 G2 T3, min with reverse complement), sorted unique u64.

**Counter:** `kc30` = the gorilla-trio `kc.c` hash counter with K = 30 (open addressing, canonical, case-insensitive,
non-ACGT resets; one FASTA stream on stdin; output int32 per Q entry). **Neighbour counter** `kn30`: for a batch of
candidate k-mers, builds the hash of their 90 canonical Hamming-1 variants, streams the reference once, and returns
Σ over variants of min(occurrence, 255) (QuicK-mer2's edit depth for `-e 1`; exact, not capped at 101 — equivalent for
the < 100 test). Batches keep RSS ≤ 8 GB.

**2.1 SPEC(g)** (paralog-specific k-mers, QuicK-mer2 rule): region k-mers with CHM13-noY count == 1, code ≠ 0, and edit
depth ≤ 99. **Sensitivity arm NOED:** without the edit-depth filter.

**2.2 parCN_G(g) = s_G × median over SPEC(g) of c_G(k)**, c_G = exact count in genome G's FASTA (every record). s = 1
for diploid sums: **HG002 v1.1** (`hg002v1.1.fasta.gz`, S3 human-pangenomics, 1,769,801,343 bytes, 47 records,
`chrN_MATERNAL` + `chrN_PATERNAL` + `chrX_MATERNAL` + `chrY_PATERNAL` + `chrM`, one pass = MAT + PAT sum) and
**gorilla mGorGor1** (`gorilla_haps/mat.fa` + `pat.fa`, two passes, summed; mat and pat also reported). s = 2 for
haploid primaries: **chimpanzee mPanTro3 v2.0 pri** and **orangutan mPonPyg2 v2.0 pri** (both carry one copy of each
autosome + X + Y; ×2 assumes homozygosity — a heterozygous CNV reads as 0 or 2 extra copies; caveat attached to every
ape number) and **CHM13** (sanity: 2 for every resolved gene by construction). Secondary statistic: s × mean (closer to
QuicK-mer2's window mean).

**2.3 RESOLVED iff |SPEC(g)| ≥ 100** — fixed a priori: one fifth of Soto's 500-k-mer window; below it QuicK-mer2's own
per-gene value is dominated by k-mers outside the gene, so there is nothing gene-specific to compare. Genes below are
**UNRESOLVED** and reported as a count, never imputed.

**2.4 Controls:** 300 autosomal RefSeq protein-coding genes (`Reference/chm13v2.0_RefSeq_full.gff.gz`, `gene`
features with `gene_biotype=protein_coding`), drawn with `random.Random(20260929)`, first 20 kb of the gene span,
excluding genes whose name is in S1E/S1C and genes within 1 Mb of any S1E region. Expected CN 2 everywhere.

**2.5 famCN analogue** (only for the human-vs-ape calls; WSSD in spirit = depth over repeat-masked sequence, all copies):
FAM(g) = region k-mers with **no soft-masked (lower-case) base** in `chm13v2.0.fa` (RepeatMasker/TRF soft mask ≈ WSSD's
masked reference), any CHM13 count. p_G(g) = share of FAM(g) positions with c_G ≥ 1. **famCN_G(g) = s_G × median of
c_G over the present positions** if |FAM(g)| ≥ 100 and p_G ≥ 0.10; **0 (absent)** if p_G < 0.10; unresolved if
|FAM(g)| < 100. Why condition on presence: an exact 30-mer survives divergence d with probability (1 − d)^30 — ≈ 0.70
chimp (d ≈ 1.2%), 0.62 gorilla (1.6%), 0.39 orangutan (3.1%) — so an unconditioned median would call most orangutan
orthologs absent. p ≥ 0.10 ⇔ d ≲ 7%, about what WSSD's read mapping tolerates.

**2.6 Calls.** SOTO-H = S1C status ∈ {Duplicated in humans, Expanded in humans}; SOTO-N = "Expanded in great apes";
Undetermined = excluded (famCN > 10), described, not scored. **ASM-H(g) iff famCN_HG002(g) > max(famCN_chimp,
famCN_gorilla, famCN_orangutan)**; sub-class "duplicated" iff that ape max < 2.5, else "expanded". Evaluable iff
FAM(g) resolved (the famCN of every genome is then defined). Sensitivity: human = CHM13 × 2. Family level (S1E
`Family ID`): SOTO-CN-H family = ≥ 1 SOTO-H member; ASM-H family = ≥ 1 evaluable ASM-H member.

**2.7 Paralog-level ape descriptor** (named families only, no bar): ape parCN on SPEC(g) and the ratio
r_G(g) = [share of SPEC(g) present in G] / p_G(g). A human-derived paralog should give r ≈ 0 in every ape (its
paralog-specific alleles did not exist before the duplication); its ancestral paralog r ≈ 1.

## 3. Clauses (evaluated in this order; each verdict string fixed now)

- **C0 instrument (gate; nothing else is claimed if it fails).** (a) `kc30` and `kn30` equal a Python brute force on a
  test (≥ 50 k query k-mers × 2 Mbp of CHM13, every value identical). (b) CHM13 parCN = 2 for every resolved gene.
  (c) HG002 parCN = 2 for ≥ 97% of resolved controls, and famCN = 2 for ≥ 90% of resolved controls in each ape.
  Fail → "INSTRUMENT NOT CALIBRATED", diagnose, stop.
- **C1 (primary; Fixed).** Among resolved Fixed genes, share with |parCN_HG002 − S1E Median parCN| ≤ 0.5.
  ≥ 0.90 → "ASSEMBLY parCN REPRODUCES QuicK-mer2 ON FIXED GENES"; [0.80, 0.90) → "PARTIAL"; < 0.80 → "DISAGREES"
  (each disagreement then classified: CHM13-private allele (SPEC present in no other genome), HG002 assembly gap/CNV,
  acrocentric short arm, SPEC < 300, other).
- **C2 (Nearly-Fixed).** Same statistic ≥ 0.70 → "REPRODUCES"; else "DOES NOT".
- **C3 (Polymorphic).** Same statistic reported, no bar (HG002 may legitimately differ from the population median).
- **C4 (rank).** Spearman(parCN_HG002, S1E Median parCN) over all resolved ≥ 0.40 → "RANKS AGREE".
- **C5 (resolution).** Share of the 1,831 resolved; reported with the NOED arm beside it.
- **H1.** SOTO "Duplicated in humans", evaluable: ASM-H share ≥ 0.80 → "HUMAN DUPLICATIONS CONFIRMED".
- **H2.** SOTO "Expanded in humans": ASM-H share ≥ 0.60 → "HUMAN EXPANSIONS CONFIRMED".
- **H3.** SOTO-N: ASM-H share ≤ 0.25 → "NON-CALLS CONFIRMED"; > 0.25 → "ASSEMBLIES CALL MORE HUMAN GAINS".
- **H4.** SOTO-CN-H families: share that are ASM-H families ≥ 0.80 → "FAMILY CALLS CONFIRMED".
- **H5 (descriptive).** ASM-H share among synteny-only human-duplicated families, Undetermined families, and SOTO-H vs
  SOTO-N stratified by Soto's own margin (human median − ape max).
- **F-cal (descriptive with bars).** Spearman(famCN_HG002, S1C Median famCN) ≥ 0.60 and Spearman(ape-max famCN_asm,
  S1C Max famCN Great Apes) ≥ 0.50 over evaluable SOTO-H ∪ SOTO-N rows.
- Every disagreement in H1/H2 and every ASM-H in H3 gets a one-line reason from a fixed menu: Soto margin < 1; ape
  individual differs (our ape max < Soto's by ≥ 1 or > by ≥ 1); FAM small (< 300); HG002 ≠ CHM13 × 2; p-floor hit in
  an ape; other.
- Named families (SRGAP2 ID_462, ARHGAP11 ID_145, FAM72 ID_354, NOTCH2NL ID_400, NPIP ID_149/151-155, TBC1D3
  ID_468/469): every S1E row printed with S1E parCN, HG002 / CHM13 / ape parCN, famCN per genome, r per ape, calls.

No parameter is tuned after counts exist; nothing is selected; if a clause's bar turns out to be ill-posed it is
reported as such, not re-cut.

## 4. Not blind (disclosed)

I have read S1E/S1C value distributions (class sizes, status counts, the named-family S1C rows: e.g. SRGAP2/SRGAP2B/C/D
"Duplicated in humans" with ape max ≈ 2.1-2.2; NPIP and TBC1D3 "Undetermined" with ape max 111-131 / 40-47; NOTCH2NLA/C/R
"Expanded in humans", NOTCH2NLB "Undetermined" at famCN 10.06). No assembly count of any k-mer has been made. The
08-02 memory records that parCN has little dynamic range (77% of families in [1.5, 2.5]) — which is why C1 (not C4) is
primary.

## 5. Predictions (this author, before any count)

- C0(c) passes: 0.85. C5 resolved share in [0.65, 0.85]: 0.55 (below 0.65: 0.35).
- C1 ≥ 0.90: 0.70; DISAGREES: 0.10. C2 ≥ 0.70: 0.65. C3 share in [0.3, 0.7]: 0.60. C4 ≥ 0.40: 0.60.
- H1 ≥ 0.80: 0.75. H2 ≥ 0.60: 0.50. H3 ≤ 0.25: 0.55. H4 ≥ 0.80: 0.70.
- F-cal human Spearman ≥ 0.60: 0.70; ape ≥ 0.50: 0.50.
- SRGAP2B/C/D, ARHGAP11B, FAM72B/C/D, NOTCH2NLA/B/C: r ≤ 0.2 in all three apes: 0.70; SRGAP2, ARHGAP11A: r ≥ 0.6 in
  chimp: 0.70. NPIP and TBC1D3 rows: not ASM-H (apes have more): 0.80.

## 6. Hostile self-review (applied above)

1. **Exact count ≠ read depth.** QuicK-mer2 values are GC-normalised means with ±0.2-0.3 noise; S1E is a population
   median; HG002 (Ashkenazi) is not in 1KGP. Only Fixed / Nearly-Fixed carry bars; polymorphic genes are descriptive.
2. **Reference version.** v2.0-noY (= v1.1) vs Soto's v1.0: uniqueness can differ on the acrocentric short arms
   (rDNA models changed). C1/C2 are also reported without genes on chr13/14/15/21/22 below 18 Mb.
3. **CHM13-private alleles.** A "unique" k-mer that is a CHM13 SNP allele, not a PSV, reads 0 in every other genome.
   The median over ≥ 100 k-mers tolerates < 50% such k-mers; SPEC present in no other genome is a named failure class.
4. **Integer counts vs continuous S1E.** The ±0.5 tolerance is the only fair comparison; the secondary mean is shown.
5. **Different individuals and 3 vs 4 apes.** Our apes are not Soto's (Clint / Mhudiblu / Kamilah / Susie; no bonobo;
   Bornean vs Sumatran orangutan). Max over 3 is ≤ max over 4 ⇒ ASM-H is biased toward calling human gains; H3 is the
   guard. Ape CN polymorphism is a named disagreement reason.
6. **Haploid primaries × 2** (chimp, orangutan) assume homozygosity — stated with every number.
7. **famCN analogue ≠ WSSD.** Presence-conditioning selects the ape copies closest to human; a family whose ape copies
   are all > 7% diverged reads as absent (0) — biases toward ASM-H. The p-floor is a named reason.
8. **HG002 is one person** vs Soto's SGDP median; CHM13 × 2 sensitivity shows how much the call depends on the person.
9. **No circularity:** S1E/S1C values never enter k-mer selection, thresholds or calls.
10. **The H clauses test Soto's WSSD-based calls, not QuicK-mer2**; they say whether assemblies reproduce the ape
    comparison, not whether parCN does.

## 7. Order and machine rules

(1) this file (sha1) → (2) build Q, `kc30`, `kn30`; C0(a) test → (3) **Amendment 1** (sha1s, |Q|, candidate count;
no parCN) → (4) counts, one genome per `tools/rlock.sh heavy` call (foreground, ≤ 600 s each; HG002 streamed from gz;
TMPDIR under `/mnt/linuxdisk`), `kn30` batches on CHM13-noY → (5) one analysis run under `tools/rlock.sh light` →
(6) Outcome + draft register rows (numbered 1161B…, **not appended** — a sibling task also drafts rows after 1160).
New data < 80 GB; delete only files this task created. Never `pkill -f`.

## Amendments

### Amendment 1 — the freeze (2026-09-29 14:05, written BEFORE any non-CHM13 genome was counted)

- The body above (§0-§7) is frozen at **sha1 2a76a930** (written 13:57).
- Instrument (scratch `soto_parcn_asm/`, 8-char sha1): `build_q.py` deb4ec95, `kc30.c` 44bb2ee8 (= `dna_cn/o3x/kc.c`
  with K = 30, header comment only otherwise), `kn30.c` c3aa5167, `test_c0a.py` a3a5307d, binaries `kc30` 1c1a6458 /
  `kn30` 38966225 (gcc -O3 -march=native).
- **C0(a) PASSES:** kc30 identical to brute force on 46,932 query k-mers (36,932 present) over 2 Mbp (CHM13
  chr16:28-29 Mb + chr1:121-122 Mb); kn30 identical on 2,998 candidates (372 with edit depth > 0, max 263).
- Q (`work/Q.u64` 5a0c73d5): 1,831 S1E regions + 300 controls (13,723 eligible controls; seed 20260929) → 33,928,439
  region positions, **17,538,361 distinct canonical 30-mers**; 62.4% of positions touch a soft-masked base.
  `regions.tsv` 01f1afaa.
- **Order deviation (disclosed):** the CHM13-noY count of Q (§7 step 4) was run before this amendment, because the
  kn30 candidate list needs it. It yields no parCN beyond CHM13's own (= 2 by construction): every Q k-mer is present
  (0 zero counts), 7,229,055 have count 1 (candidates, `cand_all.u64` b30d1611, split in 3 batches), 10,309,306 > 1.
- No HG002, gorilla, chimp or orangutan k-mer has been counted.

### Amendment 2 — the run (2026-09-29 14:40; no rule, threshold, clause or bar changed)

- Counts (`tools/rlock.sh heavy`, foreground, kc30 on each genome's full FASTA): HG002 v1.1 3:46 (47 records, 5.999
  Gbp), gorilla mat 2:04 (225 rec) / pat 1:55 (24 rec), chimp 1:56 (26 rec), orangutan 1:50 (26 rec); kn30 edit
  depth on CHM13-noY in 3 batches of 2,409,685 (7:32 / 6:15 / 6:13, 5.3 GB RSS each). The edit-depth filter drops
  31,861 of the 7,229,055 candidates (0.44%; 37.4% have edit depth > 0).
- `analyze.py` (10f0382e) was written after Amendment 1 and after the counts, **before any parCN / famCN value of a
  non-CHM13 genome was looked at**; it implements §2-§3 as written. One analysis run (`rlock light`, 10 s). Outputs:
  `results_s1e.tsv` 845aa80a, `results_ctrl.tsv` de54b561, `results_named.tsv` 84abe946, `summary.json` 02827af7.
  `describe.py` (8cb0537f) = descriptive tables only (per-ape rates, reason lists, size strata), no clause.
- Reporting detail (not a re-cut): for H3 the first reason on the menu ("Soto margin < 1") is true for every SOTO-N
  row by construction, so H3 reasons are reported as ALL matching reasons.

## Outcome (2026-09-29)

Per genome throughout; human (CHM13 / HG002) and ape numbers are never pooled. Ape famCN from chimpanzee and orangutan
primaries is haploid × 2.

**C0 — instrument CALIBRATED.** (a) passed (Amendment 1). (b) CHM13 parCN = 2 for all 1,163 resolved genes.
(c) Controls: HG002 parCN = 2 for **297/299 (99.3%)**; famCN = 2 for 98.3% (chimp), 98.3% (gorilla), 97.6%
(orangutan), 98.0% (HG002) of 297 resolved controls. Median share of control FAM k-mers present: HG002 0.999, chimp
0.731, gorilla 0.676, orangutan 0.442 — the (1 − d)^30 divergence loss §2.5 anticipated; without presence-conditioning
orangutan parCN of the controls is 2 in only 9.0%.

**C5 — resolution.** 1,163 / 1,831 (63.5%) S1E rows have ≥ 100 QuicK-mer2-rule k-mers (NOED arm 1,164): Fixed
322/340, Nearly-Fixed 629/978, Polymorphic 212/513. Unresolved regions are short (median 1.27 kb vs 11.4 kb resolved).

**C1 — "ASSEMBLY parCN REPRODUCES QuicK-mer2 ON FIXED GENES": 321/322 = 0.997** (|HG002 − S1E median| ≤ 0.5; mean
statistic 0.994; without acrocentric short arms 0.997; NOED 0.997). The one miss is NBPF3 (S1E 1.90, HG002 1; mean 1.62).

**C2 — "DOES NOT" (bar 0.70): Nearly-Fixed 408/629 = 0.649** (mean statistic 0.692; without acrocentric short arms
0.663; NOED 0.645). HG002 values: 430 × 2, 70 × 1, 48 × 0, 51 × 3, 27 ≥ 4. Disagreements (221; 127 HG002 below S1E,
94 above) by the pre-registered menu: HG002 famCN ≠ CHM13 famCN 73, SPEC < 300 52, CHM13-private 31, acrocentric 22,
other 43. Agreement rises with the number of paralog-specific k-mers: 0.53 (100-300), 0.63 (300-1k), 0.72 (1k-3k),
0.76 (≥ 3k). **Mechanism (descriptive):** in 22.4% of resolved Nearly-Fixed genes HG002 carries < 75% of the CHM13
paralog-specific k-mers (Fixed: 1.6%) — the "paralog-specific" alleles of one haploid reference are partly
haplotype-specific in these genes, so an exact per-k-mer median becomes bimodal (0 or 2) where QuicK-mer2's 500-k-mer
window mean, which also averages in flanking unique sequence, stays near 1-2. Example: 48 Nearly-Fixed genes (parCN >
0.5 in ≥ 98% of 1KGP) read 0 in HG002; 38 of them carry 15-49% of their SPEC k-mers (median 33%; NBPF4, CHRNA7, RGPD1,
NPIPB6/B9, CLEC18A …) and 6 carry < 5% (CHM13-private: NPIPB7, SPDYE17, MST1P2, GOLGA6L4 …).

**C3 (no bar):** Polymorphic 83/212 = 0.392 (mean 0.448). **C4 — "RANKS AGREE":** Spearman 0.414 (mean 0.488).

**H1 — "HUMAN DUPLICATIONS CONFIRMED": 105/109 = 0.963** (CHM13 × 2: 0.963); 108/109 have an assembly ape max < 2.5,
i.e. the assemblies also say "duplicated", not "expanded". The 4 misses all have < 300 FAM k-mers (DNM1P24,
AC244255.1 ×2, AC092666.2), 3 with Soto margin < 1.
**H2 — "HUMAN EXPANSIONS CONFIRMED": 105/132 = 0.795** (CHM13 × 2: 0.773; sub-class exp 62 / dup 43). Misses (27):
Soto margin < 1 in 15, our ape max differs from Soto's by ≥ 1 in 9 (AMYP1 ×3: ours 4 vs Soto 8.7-8.9; TRIM64DP/EP;
LRRC37BP1), HG002 ≠ CHM13 × 2 in 6, other 7 (SEC22B2P/3P, DDX12P, ANAPC1P5 …: human = chimp at 6 or 4).
**H3 — "NON-CALLS CONFIRMED": 131/631 = 0.208** ASM-H among "Expanded in great apes" (CHM13 × 2: 0.212). Of the 131,
102 have our ape max ≥ 1 below Soto's (our 3 apes are not Soto's 4; no bonobo), 93 have HG002 ≠ CHM13 × 2. Soto's own 41 internally inconsistent rows (human > ape max yet
labelled N; all margins ≤ 0.48; 40 evaluable) form the (0, 1] margin stratum, where ASM-H is 0.45 (vs 0.19 at margin ≤ 0).
**H4 — "FAMILY CALLS CONFIRMED": 105/118 = 0.890** families with a Soto CN-called member (CHM13 × 2: 0.864). Missed:
ID_12, 14, 43, 49, 50, 122, 128, 205, 225, 292, 325, 430, 448.
**H5 (descriptive):** synteny-only "Human duplicated" families 61/118 ASM-H (0.517); Soto-"Undetermined" families
72/335 (0.215); among Soto-excluded paralogs (famCN > 10) 234/699 are ASM-H — the subtelomeric DDX11L / WASH / OR4F
families (HG002 18-21 vs ape max 6) lead, which Soto's famCN > 10 exclusion removed from the comparison.
Per-ape (HG002 > that ape): Duplicated 0.972 chimp / 0.991 gorilla / 0.963 orangutan; Expanded 0.871 / 0.909 / 0.886;
"Expanded in great apes" 0.436 / 0.426 / 0.700. The ape max is chimp in 570, gorilla 210, orangutan 92 of 872 rows.
**F-cal:** Spearman(famCN_HG002, S1C Median famCN) **0.617** (CHM13 × 2 0.613; bar 0.60, passes narrowly);
Spearman(ape-max famCN, S1C Max famCN Great Apes) **0.886** (bar 0.50).

**Named families** (full rows `results_named.tsv`):
- **SRGAP2 (ID_462):** HG002 parCN 2/2/2/2 (S1E 2.06/2.05/1.70/1.60 for SRGAP2/B/C/D); famCN HG002 6/6/6/8, every ape 2
  → dup, = Soto. Paralog-specific k-mer presence ratio r in apes: SRGAP2 0.29 vs B/C/D 0.03-0.09 — the assemblies
  single out SRGAP2(A) as the ancestral copy.
- **ARHGAP11 (ID_145):** parCN 2/2 (S1E 2.05/2.08); famCN 4 vs apes 2 → dup, = Soto. r does NOT separate A (0.46-0.58)
  from B (0.39-0.50).
- **FAM72 (ID_354):** parCN 2 × 4; famCN 8 vs 2 → dup, = Soto; r FAM72A 0.27-0.31 vs B/C/D 0.01-0.12.
- **NOTCH2NL (ID_400):** famCN HG002 9 (CHM13 10) vs chimp 6 / gorilla 6 / orangutan 2 → exp, = Soto's "Expanded in
  humans" for NOTCH2 / NOTCH2NLA / C / R; NOTCH2NLB and NBPF26 (Soto excluded at famCN 10.06 / 10.23) also exp. HG002
  parCN NOTCH2NLA 1 (S1E 1.41), NOTCH2NLC **1 (S1E 2.01)**, NOTCH2NLR 2 (S1E 0.85). r 0.06-0.16 for every NOTCH2NL vs
  0.21-0.28 for NOTCH2.
- **NPIP:** core ID_154 (NPIPA2/3/5/7/8, NPIPB2/6-12/14P/15, NPIPP1) HG002 famCN 29-40 vs chimp 62-80 → not human, which
  agrees with Soto's numbers (human 45-48 < ape max 117-131) although Soto labels them Undetermined (famCN > 10); ID_151/
  152/155 likewise; ID_153 NPIPB5 exp (HG002 19 vs chimp 18). **ID_149 (NPIPA1, PKD1P1/P6, AC138932.1 …): ASM-H exp in
  all 9 rows (HG002 12-16 vs ape max 10-12) where Soto says "Expanded in great apes"** (Soto ape max 7.95-20.1): ape
  individual and HG002 ≠ CHM13 × 2 (HG002 +2 copies) both on the reason list. parCN NPIPB6/B7/B9 read 0 in HG002 (S1E
  1.1-1.4; 56-95% of their SPEC k-mers are CHM13-private).
- **TBC1D3:** ID_468 HG002 famCN 31 (CHM13 22) vs orangutan 28-30 / gorilla 20 / chimp 8-10 → ASM-H exp (barely),
  where Soto has human 18 < ape max 44-46 (Undetermined, famCN > 10) — our Bornean orangutan and Soto's Susie differ,
  and HG002's famCN (31) exceeds CHM13 × 2 (22) by 9. HG002 parCN 0-0.5 for TBC1D3B/I/G/H/D/E (SPEC 120-350). ID_469
  (TBC1D3P3/P4) HG002 4 vs orangutan 18 → not human.

**Predictions scored:** C0(c) pass ✓ (0.85); resolved in [0.65, 0.85] ✗ (0.635); C1 ≥ 0.90 ✓; C2 ≥ 0.70 ✗; C3 in
[0.3, 0.7] ✓; C4 ✓; H1 ✓; H2 ✓; H3 ✓; H4 ✓; F-cal human ✓ (0.617), ape ✓; derived paralogs r ≤ 0.2 in all apes ✗
(ARHGAP11B 0.39-0.50); SRGAP2 / ARHGAP11A r ≥ 0.6 ✗ (0.29 / 0.52 — r is diluted by divergence around each PSV and by
PSVs derived on the ancestral branch); NPIP and TBC1D3 not ASM-H ✗ (NPIP ID_149 and TBC1D3 ID_468 are ASM-H).

**Verdict.** On genes whose CN is fixed in humans, exact assembly counting of QuicK-mer2's own k-mer set reproduces
Soto's parCN almost perfectly (0.997). It does not reproduce the Nearly-Fixed class (0.649): there the CHM13-defined
"paralog-specific" k-mers are partly haplotype-specific, and a single individual's exact count disagrees with a
population median of window means. The assemblies reproduce Soto's human-vs-ape famCN calls for human duplications
(0.963) and expansions (0.795) and 89% of CN-called families; the 21% of Soto's non-calls that the assemblies call
human gains are dominated by ape-individual differences (102/131). Soto's famCN > 10 exclusion hides ~234 paralogs the
assemblies call human gains (subtelomeric families first).

### Draft register rows (NOT appended; suffix B avoids a sibling task's rows after 1160)

| 1161B | 2026-09-29 | CN instruments (human, CHM13 v2.0) | Exact assembly counting of QuicK-mer2's own k-mer set (k = 30 canonical, once in CHM13-noY, Hamming-1 edit depth < 100 as Soto's `-e 1`) in HG002 v1.1 reproduces Soto's S1E parCN | ✅ **On Fixed genes: 321/322 (0.997)** within 0.5 (prereg `PREREG_soto_parcn_assembly_2026-09-29.md`, sha1 2a76a930). Controls 297/299 = 2. Resolved (≥ 100 k-mers) 1,163/1,831. ⛔ **Nearly-Fixed 0.649 < 0.70 bar** (mean 0.692): in 22% of them HG002 carries < 75% of CHM13's "paralog-specific" k-mers — they are partly haplotype-specific; agreement rises with k-mer count (0.53 → 0.76). Polymorphic 0.392; Spearman 0.414. Edit-depth filter drops 0.44% of candidates (immaterial). |
| 1162B | 2026-09-29 | CN instruments (human vs ape assemblies) | Assembly famCN analogue (exact 30-mers over CHM13 soft-unmasked gene ∩ SD98, presence-conditioned median; HG002 vs max of chimp mPanTro3 ×2, gorilla mGorGor1 mat+pat, orangutan mPonPyg2 ×2) reproduces Soto's WSSD human-vs-ape calls | ✅ **Duplicated in humans 105/109 (0.963), Expanded in humans 105/132 (0.795), CN-called families 105/118 (0.890); non-calls 131/631 ASM-H (0.208, bar ≤ 0.25)**, 102/131 of those with our ape max ≥ 1 below Soto's (different apes, no bonobo). Calibrates: Spearman vs S1C human famCN 0.617, vs great-ape max 0.886. SRGAP2 / ARHGAP11 / FAM72 dup and NOTCH2NL exp = Soto; NPIP ID_149 and TBC1D3 ID_468 disagree (ape individual + HG002 ≠ CHM13). |
| 1163B | 2026-09-29 | Soto scope | Soto's famCN > 10 exclusion ("Undetermined") removes only families where apes are not below humans | ⚠ **No: 234/699 excluded paralogs are human > every ape in the assemblies** (DDX11L / WASH / OR4F subtelomeric families HG002 18-21 vs ape max 6; TBC1D3 31 vs 28; AMY1; CBWD; AGAP). Their NPIP core (ID_154) is correctly not human (HG002 29-40 vs chimp 62-80). Descriptive, one human, three apes. |
| 1164B | 2026-09-29 | O2 / paralog identity (descriptive) | Exact presence in ape assemblies of a human paralog's specific k-mers identifies the ancestral copy of a human-specific duplication | ⚠ **Partly.** Ratio r (specific-k-mer presence / family presence): SRGAP2 0.29 vs SRGAP2B/C/D 0.03-0.09; FAM72A 0.27-0.31 vs B/C/D 0.01-0.12; NOTCH2 0.21-0.28 vs NOTCH2NL 0.06-0.16 — but ARHGAP11A 0.46-0.58 vs ARHGAP11B 0.39-0.50 (not separated); predicted r ≥ 0.6 for ancestors failed (0.29/0.52). |
