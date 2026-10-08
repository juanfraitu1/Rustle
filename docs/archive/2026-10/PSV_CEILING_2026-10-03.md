# PSV ceiling — how many copies are too similar for any PSV to place a read (2026-10-03)

Artifact: https://claude.ai/artifact/9BZzTF98t4sZmYggHDXpiq (private). Script `bench/psv_ceiling/psv_ceiling.py` (untracked); data
`/mnt/linuxdisk/tmp/psv_ceiling/{gorilla,human,chimp}.nearest.tsv`, `o2_bands.tsv`, `data.json`. Context: advisor doubts PSV-based
read-to-copy assignment is more than lucky / easy cases.

Method: per catalog copy, identity to its most similar same-family copy on the SPLICED copy sequence (minimap2 -cx asm20 all-vs-all,
best hit covering >= 50 % of the shorter copy, identity = matches / block length). "Beyond reach" = identity >= 1 - 1/2,162
(median KB3781 read length) = 99.954 %: fewer than one expected distinguishing site per read.

| species (catalog) | copies | families | no sibling alignment | identical | >= 99.9 % | beyond reach (>= 99.954 %) | median identity |
|---|---|---|---|---|---|---|---|
| gorilla KB3781 fibroblast (GWFAM) | 3,981 | 666 | 862 | 639 | 784 | 670 (21.5 % of 3,119 aligned) | 0.9936 |
| human testis (legacy catalog) | 2,986 | 663 | 1,989 | 156 | 163 | 156 (15.6 % of 997) | 0.9732 |
| chimpanzee (legacy catalog) | 8,833 | 1,507 | 6,027 | 213 | 238 | 218 (7.8 % of 2,806) | 0.9525 |

Gorilla copies by bin: <95 % 220 · 95-98 % 539 · 98-99 % 512 · 99-99.5 % 473 · 99.5-99.9 % 591 · 99.9-100 % 145 · identical 639.

Real gorilla reads (rna_allele o2_cat0..3, 41,263 contested reads, 309 families), by the bin of the copy the read sits in:
assigned 633 total, concentrated at 95-99.5 % (345 + 53 + 129) and absent at >= 99.9 % (0 + 1); "no distinguishing column" dominates
99.5-99.9 % (5,936 of 10,062) and identical (356 + 1,958 origin-rejected of 2,320); "matches no catalog copy" (origin-rejected: the
read disagrees with every candidate beyond sequencing error, i.e. the unit lacks the read's content or the catalog lacks the copy) is
the largest class in every bin below 99.5 % (e.g. 99-99.5 %: 5,819 of 6,021). Simulation with truth (human chr16, fig5 table): 163/163
assigned correct; assigned fraction 0.80 / 0.77 / 0.25 / 0.001 at 98-99 / 99-99.5 / 99.5-100 / identical.

Reading: one in five gorilla catalog copies with an aligning sibling is beyond PSV reach for a single read, and the test assigns
nothing there (abstains). Where distinguishing sites exist and the read matches a catalog copy, it assigns and is right in simulation.
The dominant real-data abstention is not missing PSVs but reads matching NO catalog copy — the O3 phenomenon (reference-absent copies,
incomplete units), which no assignment rule fixes. Caveats: catalog copies are de novo loci with reads; human/chimp legacy catalogs
have many copies with no sibling alignment (looser membership); the o2 runs cover 309 of 667 gorilla families (the pre-registration's
categories), with rna_allele's parameters (gtf copy set, origin_drop_indels).

## Follow-up (same day): the 99.5-99.9 % "no distinguishing column" block is ONE rRNA locus
Of the 5,936 such reads in that bin, 5,851 come from SM5 copy 50 (NC_073246.2:2,704,928-2,723,519, 3 exons, 1,800 bp spliced), whose
closest sibling SM5 copy 52 (NC_073246.2:3,746,169-4,110,853, 30 "exons", 365 kb: an rDNA array) differs from it at exactly 2 positions
(spliced 1,545: +C insertion; 1,751: A->G), both inside the 3' 800 bp where the reads sit (median read 800 bp, 3'-biased). Both columns
are discarded BY DESIGN: indels (`origin_drop_indels`, isoform structure) and A->G (`rna_editing_filter`). Both copies overlap RefSeq
rRNA genes (LOC129531215 etc.). Three catalog families overlap rRNA genes: SM5 (54 copies), SM7 (51), SM577 (2) = 105 copies, 7,551
contested reads. Without them the 99.5-99.9 % bin holds 3,463 reads, 76 with no distinguishing column; gorilla copies beyond reach
652 of 3,083 aligned (21.1 %), identical 624. The artifact has a toggle (default: rRNA families excluded). Catalog-hygiene point: rRNA
biotype loci should not enter the family catalog (library rRNA carry-over, homogenized tandem repeats).

## Annotation-free rRNA screen (same day; `bench/psv_ceiling/rrna_screen.sh`, `rrna_mature.fa`)
megablast of all 3,981 catalog copies against the MATURE human rRNAs (18S/5.8S/28S from U13369.1, 5S NR_023363.1), E <= 1e-20: hits in
exactly 3 families — SM5 (49/54 copies; 33 at >= 90 % over >= 300 bp), SM7 (51/51 copies, 5S at 100 % over 119 bp), SM577 (2/2, 5.8S at
89 % over 148 bp) — and NO other family. These are the same three the RefSeq rRNA biotype marks, so the screen needs no annotation of
the gorilla genome: mature rRNA is a universal sequence class (status of GT-AG / the poly(A) signal, not of a gene annotation).
Trap measured first: against the WHOLE 45S unit (spacer included) >100 families hit at ~83 % over ~300 bp (Alu-like spacer repeats).
No data-intrinsic signature singles these families out: SM0 has 488 copies, SM3 116 copies at identity 1.0, SM6/SM9/SM13 identity 1.0;
short 3'-biased reads occur in other families too (SM319, SM276). Recommended placement: read-level pre-filter before locus building
(library rRNA carry-over), or a catalog screen; Infernal/Rfam covariance models (not installed) are the species-independent form.

## What "matches no catalog copy" (origin-rejected) is — 1,498 sampled reads (seed 1, in_copy) realigned to _pri, mat, pat
Files: `/mnt/linuxdisk/tmp/psv_ceiling/orj_sample.{tsv,fa}`, `orj.{pri,mat,pat}.paf`, `orj_classes*.out`. Classes (best hit = most
matches; "good" = identity >= .999 and >= 90 % of the read aligned):
| class | reads | % |
|---|---|---|
| best _pri hit at a locus in NO catalog copy, identity 97-99.9 % (median 5 mismatches / 1,043 bp read; pri = hap identity in 397/424) | 858 | 57.3 |
| good fit at a locus in NO catalog copy: family copy homologous to that locus (>= 90 % id over >= 50 %) = UNCATALOGUED PARALOG | 258 | 17.2 |
| good fit at a locus in NO catalog copy: no homology to the family copy = unrelated gene | 142 | 9.5 |
| good fit at a locus in NO catalog copy: partial homology | 105 | 7.0 |
| noise or fragment (< 97 % or < 80 % aligned everywhere) | 50 | 3.3 |
| fits a haplotype elsewhere, not _pri at >= 99.9 % | 46 | 3.1 |
| another catalog copy (good 9 + divergent 21) | 30 | 2.0 |
| own copy (4 divergent) / allele (4) / haplotype-only copy absent from _pri (1) | 9 | 0.6 |
AS per base at the family copy (copy_assign `as_per_base_best`): origin-rejected p50 0.67 (p10 0.30) vs assigned 0.90, no-column 0.99.
Reading: ~60 % of the class is the aligner's secondary policy (-N 50 -p 0.1: weak secondaries of noisy reads from other loci; CCS error
0.1-1 % exceeds the certificate's 0.3 % budget), ~17 % are genuine paralogs the annotation-derived catalog lacks (e.g. SM244 copy at
NC_073227.2:119.1 Mb, 96 % identity over 99 % of the copy, 161 sampled reads; SM101 at NC_073228.2:86.1 Mb, 98.6 %), the rest unrelated
genes or partial homology. Genome-absent copies and alleles are rare (reference animal). The test that separates them is annotation-free:
realign the rejected reads to the genome, cluster the off-catalog loci, test homology of the family copies to them.

## Panel B population fixed (user, 2026-10-03): AS-tied reads only
`copy_assign` emits a row only for reads with >= 2 placements; its `contested` flag marks the AS ties (as_second = as_best). Panel B now uses
contested = 1 only: 40,160 reads over 277 families (1,103 non-tied rows dropped); without the three rRNA
families 32,610 reads. Uniquely placed reads never enter O2 and are in no metric. The "no sibling alignment" column was
mislabeled "not multi-mapping": those reads ARE AS-tied; their copy's siblings align over < 50 % of its length, so the identity covariate is
undefined (relabeled "identity undefined").
Counts assigned / insufficient / no column / no copy — all families: no alignment 9/10/594/11805 · <95% 40/5/44/5788 · 95-98% 0/6/14/2817 · 98-99% 7/1/5/354 · 99-99.5% 31/24/41/5817 · 99.5-99.9% 3/23/5936/4084 · 99.9-100% 0/4/43/339 · identical 1/5/356/1954.
rRNA families excluded: no alignment 9/8/584/10889 · <95% 40/5/44/5784 · 95-98% 0/6/14/2817 · 98-99% 7/1/5/354 · 99-99.5% 31/23/40/5817 · 99.5-99.9% 3/14/76/3354 · 99.9-100% 0/4/43/331 · identical 1/5/353/1948.

## Panel B rebuilt per COPY (user, 2026-10-03): does the typical copy's tied read see a distinguishing column?
Unit = catalog copy with >= 10 AS-tied reads (MAPQ < 60; primary + equal-AS secondary) and a defined sibling identity (siblings aligning
< 50 % of the copy = lost causes, excluded). y = share of the copy's tied reads with n_decisive >= 1 (a PSV or junction column where the
candidate copies differ, covered by the read), whether or not the certificate assigns. rRNA families excluded (toggle): 112 copies,
20,045 reads. Bin medians of the share with >= 1 column: <95 % 0.66 (10 copies) · 95-98 % 0.82 (13) · 98-99 % 0.92 (9) · 99-99.5 % 0.92 (23)
· 99.5-99.9 % 0.77 (30) · 99.9-100 % 0.42 (6) · identical 0.07 (21). Copies where the MAJORITY of tied reads cover a column: 5/10, 7/13,
6/9, 14/23, 17/30, 3/6, 7/21. Median share with >= 2 columns: .64/.30/.92/.92/.53/.23/.03. Copies with ANY assigned read: 1/0/0/2/1/0/0.
Reading: the advisor's "great majority of reads at similar copies have no resolvable PSV" holds for identical copies (and partly >= 99.9 %),
not below 99.9 %. Column existence != assignment: at most copies the tied reads match no catalog copy (noisy secondaries / uncatalogued
paralogs, see above), so assignment stays rare. Data `/mnt/linuxdisk/tmp/psv_ceiling/percopy_psv.out`, page version 5.

## Panel D (user, 2026-10-04): the NAMED families on real reads — per copy, share of AS-tied reads covering >= 1 distinguishing column
Data: human NPIP = bakeoff `human/ours_final2.assignments.tsv` (26 MCL copies); TBC1D3 = new `copy_assign --families` run on A119b over the
11 RefSeq TBC1D3* copies on chr17 (`/mnt/linuxdisk/tmp/psv_ceiling/named/tbc.*`; TBC1D30/31/32 are unrelated genes and were removed);
Y families = new run over RefSeq Y copies (DAZ1-4, CDY, RBMY1*, HSFY1/2, TSPY*, BPY2*, PRY, VCY; `named/yag.*`; one transcript per gene);
gorilla NPIP testis = bakeoff `mcl1_final2`. Tie rate = MAPQ-0 share of primaries in each copy's span (A119b.t2t.bam).
- NPIP human: 16,033 primaries, 24 % MAPQ 0; 23 copies with >= 10 tied reads; median share >= 1 column 0.55; 15 copies with assigned reads.
  Dominant tied copy = EIF3C/NPIPB9 fusion locus (3,348 primaries, 81 % tied; 5,326 tied reads, 1 % with a column): EIF3C has an identical
  duplicate (EIF3CL); not an NPIP PSV failure. Other large copies: NPIPA7 398 tied (0.62), NPIPB5 346 (0.55), NPIPB3 265 (0.16).
- TBC1D3 human: 9,325 primaries at the 14-gene span check, 84 % MAPQ 60, 3 % MAPQ 0; AS-tied molecules 959 (origin-rejected 667,
  contested 291: assigned 17, tied 17, ambiguous 257); copies 99.45-99.86 % to nearest sibling; per-copy share >= 1 column mostly 0.86-1.0
  (TBC1D3G 0.10 with 629 tied reads, 95 % match no copy). The aligner resolves TBC1D3 on CHM13; the tied residue mostly matches no RefSeq model.
- Y families (A119b IS testis, user 2026-10-04 — the samples registry's 'tissue unknown / notes disagree' is resolved; TSPY and VCY are simply low in this library, as the YAG prereg recorded): palindrome-arm copies HSFY1/2, RBMY1A1/B/D/J, BPY2/B/C, CDY2B, PRY are
  100 % identical to a sibling, 90-98 % of primaries tied, share >= 1 column 0.02-0.21 (RBMY1D 0.62, BPY2C 1.0 on 11 reads): the advisor's
  claim is TRUE here. DAZ2 (93.4 %) 0.98, DAZ4 (95.6 %) 0.49, TSPY4 (99.74 %) 0.96, TSPY9 0.99 — divergent Y copies do have columns.
  DAZ1 (97.3 %): 6,004 tied reads, 0.01 with a column, 91 % match no copy (DAZ repeat/isoform structure vs the single RefSeq model).
- Gorilla NPIP testis: 634 tied rows, 10 copies >= 10 reads, 9 of them 100 % "match no catalog copy": catalog mismatch (the gorilla NPIP
  loss), not a PSV question.
Message: for autosomal SD families (NPIP, TBC1D3, and the catalog families below 99.9 %) distinguishing columns exist for most copies'
tied reads; the Y palindrome families are the genuine exception; the bottleneck everywhere is reads matching no catalog copy, not PSV absence.
