# O3 — between-individual copy number: RNA cannot see it, DNA can (2026-09-25)

**Question (advisor):** copy numbers of duplicated genes differ between individuals; a sample that is not the
reference animal may carry copies the genome lacks. Can that be assessed from IsoSeq aligned to the gorilla
T2T genome, or at least *screened* from RNA and then confirmed with DNA?

**Answer:** No to both, on this data. Between-individual copy number is real and measurable with WGS
(17 / 516 autosomal families differ among three related gorillas). RNA does not measure it, and as a screen
it fires on 0 of those families while most of what it does flag is absent from DNA or allelic.

Advisor page: https://claude.ai/artifact/BsegZ9f3aCaMExP8HyRWKA
Scratch / code: `/mnt/linuxdisk/home/juanfraitu/dna_cn/o3x/` (`ks.c`, `analyze.py`, `results.json`,
`rna_copies.raw`, `page_data.json`).

⚠ **Not pre-registered.** The three cut-offs below (family difference; "present"; allele vs extra copy) were
chosen in this first pass. Treat the numbers as exploratory until a prereg reruns them.

---

## 1. Data

| animal | who | data | accession |
|---|---|---|---|
| Jim | KB3781, male, the mGorGor1 genome animal; fibroblast line | IsoSeq fibroblast (`fibroblasts/GCA_029281585.2_flnc_mm.bam`) + Illumina WGS 41.9 Gbp | SAMN04003007, SRR26039725 |
| Trib | mGorGor3, male | Illumina WGS 40.3 Gbp | SAMN35877945, SRR26044569 |
| Dolly | mGorGor2, female | Illumina WGS 44.2 Gbp | SAMN35877944, SRR26044623 |
| OR6737 | unrelated animal | IsoSeq testis (`winloci_data/GGO_mm.bam`), **no DNA** | — |

- WGS downloaded as `.sra` from NCBI ODP (AWS, ~11 MB/s; ENA gave ~0.5 MB/s) into
  `/mnt/linuxdisk/home/juanfraitu/gorilla_wgs/`, streamed with `fasterq-dump --fasta-unsorted -Z --split-spot`.
  Read and base counts match ENA exactly.
- ⚠ **The August downloads `gorilla_hifi/SRR2603972{5,6,7}_{1,2}.fastq.gz` are corrupt**: same byte size as
  ENA, md5 differs (`SRR26039725_1`: local `f7c5c7f2…` vs ENA `f10f22e8…`), gzip CRC errors, reads decode at
  a mean 232 bp instead of 151. Never use them.
- **Trio confirmed from data.** Using the mGorGor1 pat/mat haplotype assemblies: Jim's pat-only k-mers
  (188,220 seen in Jim) are in Trib 96.9% vs Dolly 39.1%; mat-only (1,157,911) in Dolly 99.7% vs Trib 89.3%.
  Trib = father, Dolly = mother.

## 2. Method (alignment-free k-mer dosage)

- Query set: canonical 21-mers from catalog exon interiors (`o1_gw/ggo_gw.copies.tsv`, 627 families,
  2,018 copies, built from the testis library) + autosomal and chrX single-copy CDS controls + every 21-mer of
  every `o3_rna_flag` consensus (fib + tes). 14.4 M distinct canonical k-mers (`qc.u64`).
- Counter: `ks.c`, sort-merge (LSD radix, 6×7-bit passes) against the sorted query, ~40 Mbp/s/core; verified
  identical to the Python `dna_cn/scan.py` counts. Counted in the primary assembly, pat, mat, and each WGS run.
- λ_hap = mean count of autosomal single-copy controls / 2: Jim 4.98, Trib 4.85, Dolly 5.19.
  chrX/auto: 0.503 / 0.496 / 1.002 (male, male, female) → a one-copy step is resolved.
- Family dosage (copies per haploid genome, assembly scale) = median over the family's k-mers of
  `w / (2 λ a)` × median `a`. Autosomal families with ≥ 200 k-mers: 516.
- **Difference rule:** |ΔD| ≥ 0.5 haploid copies **and** ≥ 20% of the smaller value.

## 3. Results

### 3a. DNA: copy number differs between individuals

| comparison | families differing (of 516) |
|---|---|
| Jim vs his own assembly (method floor) | **4** — GWFAM576 rDNA (~245 vs 17), GWFAM530 NPIP-like (4.8 vs 6), GWFAM106, GWFAM277 |
| Trib vs assembly / Dolly vs assembly | 12 / 12 |
| Trib vs Jim / Dolly vs Jim / Trib vs Dolly | 6 / 11 / 9 |
| **any pair of the three animals** | **17** |

Example: GWFAM535 (KDM8, JMJD5-like) and GWFAM536 (NSMCE1-like), chr18 (`NC_073242.2`) tandem SD, copies
~290 kb apart (orthologous to human 16p12). Dosage Jim 2.01 / 1.91, **Trib 2.78 / 2.88**, Dolly 1.99 / 1.93:
Trib carries the extra copy of the segment, and Jim inherited his other haplotype.
Full table of the 17 (+ GWFAM530) on the advisor page.

### 3b. RNA does not measure it

- **Absorption (r1089, arm A):** reads of a copy absent from the reference align 100% as MAPQ-60 primaries at
  the template for 0.5–5% divergence; 0 unmapped, 0 tied; only `de` rises (0.0015 → ≈ d). Same signature as a
  heterozygous allele.
- **Unmapped reads are negligible:** fibroblast 959 / 34.9 M alignments (0.003%), testis 5,506 / 10.7 M (0.05%).
- **Expressed copies ≠ copies, in the genome animal itself.** Fibroblast RNA sees every copy (≥ 2 primary
  reads over its exons) in 70% of the 516 families; Jim's DNA agrees with the assembly in 99%.

  | copies in family | n | RNA sees all | DNA agrees |
  |---|---|---|---|
  | 2 | 364 | 74.5% | 99.5% |
  | 3–4 | 100 | 60.0% | 100% |
  | 5–9 | 35 | 65.7% | 97.1% |
  | ≥ 10 | 17 | 35.3% | 94.1% |

  Within-family expression spread (fibroblast, expressed copies): median 7×, 90th percentile 389×, max 35,030×.
  An absorbed extra copy adds ~1.27× (excision panel), two orders of magnitude below that.

### 3c. RNA as a screen for DNA: it fails both ways

**Recall — does the RNA flagger fire where DNA shows a difference?** Of the 17 between-animal families:

| library | families with an expressed locus | fired (any `o3_rna_flag` status) | `reference_absent_candidate` |
|---|---|---|---|
| fibroblast (Jim) | 14 | 1 (GWFAM504, IGHA-like) | **0** |
| testis (OR6737) | 17 | 3 (GWFAM302, 504, 535) | **0** |

On the same-animal comparison (Jim's RNA vs Jim's 4 DNA-vs-assembly differences) the flagger fires on **0 of 4**.
⚠ Recall for the testis animal cannot be measured: it has no DNA, so its true differences are unknown.

**Precision — what the RNA flags are, once DNA is added.** "Present" = ≥ 10% of the candidate's
assembly-absent k-mers seen ≥ 3 times in the WGS. "Extra-copy-like" = the host k-mers shared with the consensus
sit at ≥ 2.6 λ (an added copy raises shared sequence to ~3; an allele leaves it at 2).

| RNA flags (`reference_absent_candidate`) | DNA | absent | allele-like | extra-copy-like | < 21 novel k-mers |
|---|---|---|---|---|---|
| fibroblast, 117 | Jim | **82** | 30 | 3 | 2 |
| fibroblast, 117 | Trib | 89 | 23 | 3 | 2 |
| fibroblast, 117 | Dolly | 86 | 26 | 3 | 2 |
| testis, 68 | Jim | 41 | 21 | 5 | 1 |
| testis, 68 | Trib | 32 | 32 | 3 | 1 |
| testis, 68 | Dolly | 32 | 23 | 12 | 1 |

- Fibroblast: 82 / 117 are not in the DNA of the animal the RNA came from → RNA-level artefacts or RNA
  processes. The 3 extra-copy-like are one locus, **CDK11B / SLC35E2B / LOC115933506**: host dosage 3.2 / 2.9 /
  3.9 λ in all three animals, carried by the maternal haplotype (0.991) → a **primary-assembly omission**, not
  between-individual variation. Extra-copy precision for between-individual CN: **0 / 117**.
- Controls: all 45 immunoglobulin-hypermutation loci have zero novel k-mers in Jim's DNA; 83% of `scattered`
  loci likewise.
- Testis: of 67 testable, 41 occur in ≥ 1 of the three genomes, 15 in Trib or Dolly but not Jim (population
  variation the reference animal lacks). RNA found expressed non-reference sequence; whether it is a copy in
  OR6737 needs OR6737's DNA.

### Why the screen is blind

1. Recent duplications (the CNV class) are > 99% identical: reads are indistinguishable from the reference copy,
   below the flagger's ~1% divergence floor.
2. A diverged extra copy and a heterozygous allele give the same `de` pile.
3. Dosage (~1.27×) is buried in expression variance (7× median, up to 35,000×).
4. A copy not expressed in the tissue produces no reads.

## 4. Consequence for O3

- The between-individual copy-number call is a DNA measurement. RNA's role is downstream: once DNA establishes
  the extra copy, RNA says whether it is expressed and which transcripts come from it (the O2 assignment
  problem).
- Updates `docs/O3_STATUS.md` §5 ("what would actually settle it"): the between-individual comparison is now
  executable on disk (Jim + Trib + Dolly WGS), but those parents have no RNA.

## 5. Next steps

1. Pre-register the difference rule and the allele/copy cut-off; rerun with all Illumina runs per animal
   (Jim 8 runs, Trib 5, Dolly 5 in PRJNA986879; ~3× the depth used here).
2. Obtain WGS for OR6737 (testis) to resolve the 15 parent-only testis sequences.
3. Record CDK11B / SLC35E2B as a primary-assembly omission (check against the HiFi reads).

## 6. Proposed register rows (not yet added — the register had uncommitted edits in another session)

| claim | verdict |
|---|---|
| Between-individual copy number can be assessed from IsoSeq aligned to the reference genome | ⛔ Not on this data. Absorbed at the paralogue (r1089); fibroblast RNA sees all copies in 70% of families vs DNA 99%; RNA flags do not follow genotype (genome animal 117, unrelated testis 68). |
| RNA can at least screen for copy-number differences that DNA then confirms | ⛔ Recall 0/17 candidates (1 fired, fib; 3 fired, testis) on DNA-variable families, 0/4 on Jim's own DNA-vs-assembly differences; of 117 fibroblast flags 82 absent from Jim's DNA, 30 allele-like, 3 extra-copy-like = one primary-assembly omission (CDK11B). |
| Copy number of gorilla gene families differs between individuals | ⭐ Yes, with WGS k-mer dosage: 17/516 autosomal families differ among Jim/Trib/Dolly, method floor 4 (Jim vs own assembly). Not pre-registered. |
