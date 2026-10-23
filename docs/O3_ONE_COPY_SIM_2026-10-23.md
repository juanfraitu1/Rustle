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
   - That case is a "look harder" flag: a slot-filler that disagrees with the copy's reference almost everywhere.
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
