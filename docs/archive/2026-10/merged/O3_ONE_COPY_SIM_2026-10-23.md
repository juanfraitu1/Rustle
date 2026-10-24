# O3 on one erased copy, simulated (2026-10-23)

**Why:**
- The user asked for one case instead of a 60-family panel: "can we somehow only focus on one erased copy? We could even simulate it
  so it's clean".
- They then named the base case: "absent copy that is hidden in multi mapping reads ... [the advisor] believes multi mapping reads
  are everywhere there is similarity so pseudogenes should be trackable".
- They chose to run both cases: the copy missing from the genome first, then the copy present but unannotated.

**Code:** `bench/o3_graph/o3sim_one.py` (pick / sim / report). Work dir: `/mnt/linuxdisk/tmp/o3sim1`.

## Setup
- **Family:** gorilla KB3781 GWFAM248, 5 copies, with GWFAM248:3 erased. The selection rule was fixed before any port run:
  - the simplest family (fewest copies, >= 3) whose erased copy's transcript is >= 97% identical over >= 80% of its length to
    >= 2 survivors;
  - >= half of the erased copy's reads must be multimappers, i.e. have a record within 10 AS of their best on >= 2 survivors.
  - The erased copy is 97.4% identical to survivors :1 and :2. 85% of its pick reads tie on >= 2 survivors.
- **Truth transcripts:** each copy's dominant fibroblast isoform (the most frequent intron chain of the fibroblast BAM's primary
  reads), taken on both haplotypes (A = `_pri`; B = lifted to the other haplotype, `isocon_sim.lift_blocks`).
  - In KB3781, :3 is homozygous (A = B).
  - :1 and :2 differ at 1 and 0 positions.
- **Reads:** 100 per allele per copy. The model is `isocon_sim`'s: 5' / 3' truncation, 0.5% substitutions, 0.25% insertions,
  0.25% deletions, seed 1. Every read's origin is known.
- **Reference:** a mini reference, i.e. the family's loci +- 50 kb, aligned with the Amendment 7 command (`-N 50 -p 0.1`).
- **Port:** v5 defaults, on the copy table without :3.

## Results
| case | erased copy's reads | the port | where the family's reads land outside the copy table |
|---|---|---|---|
| **missing** (masked in the genome, expressed) | 200 of 200 on survivors; 163 are multimappers (>= 2 piles) | one haplotype of 184 reads, **all from :3**; it takes the free second-allele slot of survivor :1 and is **not certified** | an unlisted similar locus, 473 secondary-only reads |
| missing, survivors heterozygous (**synthetic**: 6 substitutions on :1 / :2's B alleles) | the same | the same 184-read haplotype is **certified** (35 private alleles; 0 survivor reads in it); contig 1,306 bp = :3's transcript (1,623 bp) at **99.92%** (0 mismatches, 1 indel) | the same |
| **unannotated**, expressed (in the genome, not in the copy table) | its own locus (200 reads within 10 AS of their best) | nothing to certify on the survivors | **on :3: 715 reads** (200 near-best, 515 secondary) |
| unannotated, **silent** (a pseudogene: no reads of its own) | none | nothing | **on :3: 516 secondary-only reads** from the survivors' reads |

## What it shows
1. **The print of the missing copy is there and the port finds it.**
   - The missing copy's reads hide as multimappers across the survivors, as the advisor predicts.
   - Clustering by allele pattern pulls all of them into one haplotype that contains nothing else.
   - On survivor :1 that haplotype carries a non-reference allele at 26 of 27 informative columns; :1's own haplotype carries 0.
2. **Proving it is extra depends on the survivors (Theorem A).**
   - Where the survivors' alleles are distinguishable, both of their slots are full and the haplotype is certified, with a 99.9% contig.
   - KB3781's real survivors are homozygous here. A homozygous survivor has a free slot, and RNA cannot rule out that the haplotype
     is that survivor's second allele, so it is not certified.
   - That case is a putative copy: a slot-filler that disagrees with the copy's reference almost everywhere.
3. **Pseudogenes in the genome are trackable from multimappers alone.** The advisor's claim holds here:
   - Silent :3 collects 516 secondary alignments from its paralogs' reads.
   - So does an unlisted similar locus near the family (280-473 reads).
   - The signal is a pile of secondary alignments at a locus outside the copy table; none of them is within 10 AS of the read's best.
4. **A copy missing from the genome and silent leaves no reads,** so RNA has nothing to find (not simulated).

## Reproduce
```
O=/mnt/linuxdisk/tmp/o3sim1
python3 bench/o3_graph/o3sim_one.py pick --out $O
python3 bench/o3_graph/o3sim_one.py sim --out $O --name missing --case missing
python3 bench/o3_graph/o3sim_one.py sim --out $O --name missing_shet --case missing --s-het 6
python3 bench/o3_graph/o3sim_one.py sim --out $O --name unannot_expressed --case unannotated
python3 bench/o3_graph/o3sim_one.py sim --out $O --name unannot_silent --case unannotated --e-reads-per-allele 0
python3 bench/o3_graph/o3sim_one.py report --out $O --name <name>
```

## Sweep (port v6; `o3sim_one.py sweep`, 40 runs, `/mnt/linuxdisk/tmp/o3sim1/sweep.tsv`)
**Setup.**
- **Another individual** (the user: "the method should somehow account for individual differences of not having used the same
  individual for the genome and for the RNA samples"): the sample's own variants, which the reference lacks. Per copy, 0.1% of
  sites change on both alleles and 0.1% on one allele, a human-like rate. This is the CHM13-vs-HG002 situation.
- **Divergence:** E's 34 differences from its closest survivor :1 are kept with probability f, the rest reverted to :1's base.
  f = 0.5 / 0.25 / 0.1 / 0.05 gives about 1.3 / 0.65 / 0.26 / 0.13% divergence.
- **Heterozygous survivors** (synthetic): survivors whose real alleles differ at < 2 positions get a B allele with 6 substitutions.

| E (the erased copy) | same individual, survivors as in KB3781 | same individual, survivors heterozygous | other individual, survivors as in KB3781 | other individual, survivors heterozygous |
|---|---|---|---|---|
| 100 / 30 / 10 / 5 reads per allele | **putative** (pure print; contig 99.9-100%) | **certified** | certified at 100 / 30 / 10; missed at 5 | certified at 100 / 30 / 10; missed at 5 |
| 3 reads per allele | missed (below the 6-read floor) | missed | missed | missed |
| f 0.5 / 0.25 (1.3 / 0.65%) | silent in :1's free slot | certified | certified | certified |
| f 0.1 (0.26%) | silent | count right, but :1's own allele is named extra | certified | count right, wrong one named |
| f 0.05 (0.13%) | silent | count right, wrong one named | silent | count right, wrong one named |
| control (nothing erased) | 0 false calls | 0 false calls | 0 false calls | 0 false calls |

**What the sweep shows.**
1. **The reference being another individual does not create false calls (0 in every control).** The certificate never compares
   with the reference sequence to decide:
   - a column is informative only when the sample's own reads carry >= 2 alleles, so a variant the sample carries on both
     alleles is invisible to it;
   - haplotypes, the count against slots, and the putative-copy flag (tie patterns) all compare the sample with itself.
   Only two parts read the reference:
   - the slot order (closest to the reference first) decides which haplotype is named extra when the count exceeds the slots;
   - the DADA2 rule's reference parents.
2. **Another individual helps more than it hurts here.** The sample's own heterozygous sites make KB3781's homozygous survivors
   heterozygous, which fills their slots, so the erased copy becomes certifiable. It costs depth, though: the erased copy's two
   alleles now differ, so each needs 6 reads on its own (missed at 5 per allele).
3. **The limits are depth and divergence.**
   - Depth: about 5-10 reads per allele.
   - Divergence: below about 0.3% from its closest survivor (<= 4 differences), the extra copy is still counted, but RNA cannot
     say which of the three haplotypes is the extra one, so a survivor's allele can be named instead.
4. **The sweep found a real bug (v6).** Before the fix, the conflicting set counted 2-3-read artifact clusters. Three of them beat
   a real 93-read allele, freeing its slot for the erased copy (`expr100_6_other`: putative, not certified). The 6-read floor
   now applies before the count. On real development data this fix moved P1 from 5,762 to 6,129 D right (passes) with P2 still
   0 / 53 (`docs/PREREG_o3_port_2026-10-23.md`, development log).

## Correlation clustering in place of label propagation (`--clusterer cc`; sweep in `/mnt/linuxdisk/tmp/o3sim1_cc`)
Same 40 runs; the cells that change:

| run | label propagation | correlation clustering |
|---|---|---|
| other individual, 100 reads per allele: the erased copy's pure haplotype | 110 reads | **179** |
| other individual, 30 / 10 per allele | 22 / 12 reads | 49 / 16 |
| other individual, 5 per allele | missed | **certified** (7 reads) |
| same individual, 10 / 5 per allele | 15 / 6 reads | 19 / 9 |
| controls (nothing erased) | 0 false calls | 0 false calls |
| erased copy <= 0.26% from a survivor | count right, wrong one named | the same |

**Why it helps.** When the sample is another individual, the erased copy's two alleles differ at a site or two.
- Label propagation links only reads that disagree nowhere, so it splits the copy into two half-size haplotypes.
- Correlation clustering keeps them together (one differing column is one observation).

The divergence limit is the slot logic's, not the clustering's.
