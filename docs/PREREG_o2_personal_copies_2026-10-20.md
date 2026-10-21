# Pre-registration: personal copies against the corrector, on gorilla KB3781's own reads with exact truth (2026-10-20)

**User (2026-10-20):**
- "lets go back to O2, and assume we have a family with members with different degrees of similarity (90 to slighly less
  than 100%); what is still missing to be able to disambibuate (ambiguous AS) to know the real provenance and use that for
  assembly"
- On the proposal to test per-copy sample sequences against the corrector's pileup: "yes please, also maybe we can use
  gorilla for this since we have the fibroblasts and T2T genome?"
- On the design below: "yes please proceed".

## The question

The corrector (`o2correct.py donor`, registers 1343-1346) learns each copy's sample alleles as a **pileup**: per copy, how
many of its confident reads carry each substitution or indel. A tied read's events are then checked one by one against
that pileup.

A **personal copy** is the copy's own sequence in the sample: the reference copy with the events its confident reads
carry written in, split into at most two haplotypes by the reads that carry them. A tied read is realigned to the personal
copies of its tied copies.

The personal copy keeps what the pileup loses: which events sit together on one haplotype, and the read's alignment to the
sample's sequence rather than to the reference.

The experiment asks:
1. Do personal copies place tied reads better than the pileup?
2. How close does either get to the animal's real copies?

## Substrate: one animal, his reads and his own genome

- **The reads.** KB3781 ("Jim", the animal of the mGorGor1 T2T assembly): fibroblast Iso-Seq, SRR27178662, 13,516,443
  molecules and 34,762,361 records. They were aligned to the RefSeq primary GCF_029281585.2 (`GGO.fasta`) by minimap2 2.31
  `-ax splice:hq -uf --eqx -Y -N 50 -p 0.1 --secondary=yes`
  (`/mnt/linuxdisk/home/juanfraitu/fibroblasts/GCA_029281585.2_flnc_mm.bam`).
- **The reference is a mosaic of his two haplotypes** (`rna_allele/chrmap.tsv`, sha1 e9182791, each chromosome checked
  identical in sequence):
  - his father's: chromosomes 1, 2 (NC_086017.1), 3, 4, 7, 8, 9, 12, 13, 14, 15, 16, 18, 19, 21 and Y;
  - his mother's: 5, 6, 10, 11, 17, 20, 22, 23 (NC_086018.1) and X.
  - This corrects the design note in chat, which said the autosomes were all his father's.
- **The ALT haplotype** is the other parent's copy of each autosome, plus every unplaced contig of either assembly. His
  assemblies are GCA_028885475.2 (paternal) and GCA_028885495.2 (maternal), in `gorilla_haps/{pat,mat}.fa`.
- **Two built-in groups of reads.** A read from the haplotype the reference was taken from has a sample copy identical to
  the reference: the **one-genome control**. A read from the ALT haplotype carries the other parent's alleles: the **test**.

## The families

- **Catalog.** The KB3781 fibroblast catalog GWFAM: `rna_allele/fibro_cat.copies.tsv`, sha1 7521851f, 666 families, 3,981
  copies.
- **Rule.** A family qualifies when every pair of its members' main transcripts aligns at identity in [0.90, 1.00), covering
  at least 50% of the shorter. The test is `minimap2 -X -c -x asm20` with `o2perfamily.qualifies`, the human list's rule.
  - It is run one family at a time (`o2personal.py ava`), since only same-family pairs count.
  - A family of more than 20 members is first tested on 20 of them; a pair failing there fails the family.
- **Result.** 283 families qualify, with 716 copies (`personal/qual.json`, sha1 f4bab888; `qual.copies.tsv`, sha1 78a3afa8).
  NPIP does not qualify.
- **The hold-out is by family.** The parity of sha1(family name), fixed before any result was seen (`o2personal.half`):
  - development: 140 families;
  - held-out: 143 families.

  | minimum pair identity | development | held-out |
  |---|---|---|
  | [0.90, 0.95) | 15 | 11 |
  | [0.95, 0.98) | 50 | 64 |
  | [0.98, 0.99) | 37 | 33 |
  | [0.99, 1.00) | 38 | 35 |

## The domain (the corrector's, unchanged)

A read is in the domain when:
- its alignments within 10 AS of its best lie on two or more copies of one qualifying family;
- none of them lies on a locus off those copies;
- minimap2's primary is among them.

The reads are every alignment of the 450,984 reads with any alignment on a qualifying copy (`personal/sub.bam`, 4,055,491
records, `names.txt` sha1 26dd6d3c).

## The truth: the read's best placement in the animal's own genome

**Alignment.** Each domain read is aligned to his maternal and his paternal assembly (`mGorGor1.{mat,pat}.splice.mmi`):
minimap2 2.30, `-ax splice:hq -uf --eqx -Y -N 50 -p 0.1 --secondary=yes`.

**Each alignment maps to a reference copy:**
- On a haplotype contig the reference is identical to, it maps to the copy it overlaps there, in the reference's
  coordinates (`summary.attribute`, the corrector's own attribution).
- On an ALT contig, it maps to the copy whose ALT ortholog it overlaps. A copy's ALT ortholog:
  - its reference span ± 20 kb, aligned with `minimap2 -x asm20 -c` to the other haplotype's whole assembly;
  - the longest primary alignment with MAPQ >= 20 through which its exons lift (`donorsim.pick_ortholog`);
  - on the homologous chromosome or an unplaced contig;
  - two copies whose orthologs overlap are both dropped (`donorsim.drop_ambiguous`).

**The truth** (`o2personal.truth_of`) is the one copy that every alignment within 5 AS of the read's best lies on. 5 AS is
one substitution under splice:hq (match 1 + mismatch 4).

| status | meaning |
|---|---|
| determined, **refcons** | The reference's own copy is among them: the sample's copy equals the reference over the read. |
| determined, **altonly** | Only the ALT haplotype's copy: the read carries the other parent's alleles. |
| undetermined | Two copies, or a copy and an off-copy locus, lie within 5 AS. No method can place the read. |
| off-copy | The best placement is on no reference copy: a copy the reference lacks, an unannotated locus, or an ortholog that did not lift. |

**Scored reads (S):** determined, with the truth among the read's tied copies.

## The arms

Each arm gives one copy per read. An arm that abstains keeps minimap2's primary.

1. **minimap2**: its primary.
2. **Corrector**: `o2correct.py donor`, version 2 without `--guard`, unchanged since bfa15f12.
3. **Personal copies, phased** (`o2personal.py`). Per copy:
   - **Reads:** the corrector's own confident reads (primary on the copy, every alignment on another copy of its family more
     than 9 AS lower).
   - **Events:** those the corrector calls explained (>= 2 reads, >= 10% of a depth >= 3). Events at >= 80% of the depth are
     **carried**.
   - **Haplotypes:** the other events split into at most two haplotypes by a 2-means over the reads (`o2personal.phase`).
   - **Sequence:** each haplotype is the reference span ± 5 kb with its events written in.
   - **Decision:** each tied read is realigned to the personal copies (the same minimap2 options, `-N 200`). The tied copy
     with the unique best AS takes it. A shared best, or a tied copy the read does not align to, keeps minimap2's pick.
4. **Oracle**: the same decision over the copies' true sequences: the reference span ± 5 kb and its ALT ortholog ± 5 kb.
5. **Exploratory, no bar**: personal copies unphased, one sequence per copy (the events more than half the depth carries).

## Scoring (`o2personal.py score`)

**Per arm, on S:**
- right and accuracy;
- against minimap2: **fixes** (minimap2 wrong, arm right) and **breaks** (minimap2 right, arm wrong);
- reads decided, and decided right.

**Strata:**
- refcons and altonly;
- exact ties (AS margin 0 between the two best tied copies) and near ties;
- the family's minimum pair identity, in the four bins above;
- each family.

**Also reported:** undetermined, off-copy, and truth outside the tie set.

## Bars, held-out half

- **G0 (gate):** oracle accuracy on S >= 0.99. Failing it means a pipeline error: INVALID.
- **P (power):** |S| >= 300, from >= 11 families, with minimap2 wrong on >= 50 of them.

**V1. Personal copies beat the corrector:**
- **G1:** personal accuracy >= the corrector's + 1.0 point.
- **S1:** personal fixes >= 5 × personal breaks.
- **F1:** >= 11 families where personal is right more often than the corrector, and no more than a fifth as many (rounded
  up) where it is right less often.
- **H1:** on refcons reads, personal accuracy >= minimap2's − 0.2 points.

**V2. The corrector beats minimap2** (secondary; its first real-read test with exact truth):
- **G2:** corrector accuracy >= minimap2's + 2.0 points.
- **S2:** fixes >= 5 × breaks.
- **F2:** >= 11 families where it fixes more than it breaks, and no more than a fifth as many (rounded up) where it breaks
  more than it fixes.
- **H2:** on refcons reads, corrector accuracy >= minimap2's − 0.2 points.

## Verdict

- **INVALID:** G0 fails.
- **UNDERPOWERED** (both V1 and V2): G0 holds and P fails.
- **PASS:** all four of that verdict's bars hold. **FAIL:** any of them does not.

**A V1 PASS would claim:** on one gorilla's own fibroblast reads, with exact truth from his own assemblies, the sample's
phased copy sequences place tied reads better than the corrector's pileup. This holds in more than 10 near-identical
families and leaves reads whose copy equals the reference unharmed.

**It would not claim:**
- other tissues or animals;
- copies absent from the reference (off-copy: O3);
- reads tied in the sample itself (undetermined);
- a better assembly (no second round here).

## Development, then the freeze

The development half (140 families) may be run and looked at freely. On it, I may:
- fix bugs;
- change the personal copies' construction, the realignment's minimap2 options and FLANK_P;
- change the truth's mechanics (contig names, ortholog lift).

Fixed for the whole experiment:
- the truth margin (5 AS);
- the domain;
- the corrector (bfa15f12);
- the bars.

Every change goes into Amendment 1, with its reason, at a commit. The held-out half is scored once, with the code at that
commit.

## Rules

- Gorilla numbers are never pooled with human ones.
- Heavy jobs run in the foreground, one at a time: the truth alignments (about 14 GB of index each) and the ALT ortholog
  alignments.
- The metric-auditor checks the held-out numbers before any verdict is reported.
- Not run, and kept as a later replication: the mirror, with reference = the other haplotype.

## Amendments

### Amendment 1 (2026-10-20): allele-blind O2, and an edit-based truth. Written before any arm was scored against it

**Why.**

1. **The user stopped the phased arm:** "I think O2 should not care about alleles, phasing, paternal and maternal
   haplotypes. It should only be about finding the true place of multimappers based on PSVs, BAM tags and other info".
   - Dropped: the personal copies (arms 3 and 5) and the oracle (arm 4, his assembled copies as targets), with gate G0.
   - His two assemblies stay as the **answer key only**. No arm sees them.
2. **AS was the wrong yardstick for the truth.** AS charges every intron, so a copy with fewer introns (a retrocopy) can
   outscore the spliced gene the read matches perfectly. This was seen on the development half only:
   - Of its 3,484 domain reads, 2,592 have minimap2's primary below the highest-AS tied copy. In 2,539 of them the primary
     has fewer edits (NM); in 2,540 it has more introns.
   - Example: SRR27178662.117798 matches SM438.1 (minimap2's primary, AS 1,196) with NM 0, and SM438.0 (AS 1,204, fewer
     introns) with NM 24.
   - The AS truth (5-AS margin) therefore called minimap2 wrong on 740 development refcons reads like these.
   - minimap2's primary has the fewest NM of the tied copies on 3,405 of the 3,484.

**What had been seen before this amendment:**
- **Development:**
  - minimap2's primary against the AS truth, by origin and margin;
  - the AS / NM / intron comparison above;
  - the corrector's output (`arm2.tsv`), scored against no truth.
- **Held-out:** only the AS truth's status counts: determined altonly 404, determined refcons 215, off-copy 7,
  undetermined 270. No arm's calls were compared with any truth.
- **Personal copies:** built (`personal.*.fa`), never realigned or scored. `o2personal._depths` was sped up (bisect),
  without changing its output.

**The truth** (replaces the section above). The alignments, and their mapping to reference copies, are unchanged
(`truth.mat.bam`, `truth.pat.bam`, `orthologs.tsv`).
- **An alignment's edits:** its NM plus its clipped read bases (soft and hard), i.e. the read bases it gets wrong or does
  not explain. Introns are free.
- **The truth** (`o2personal.truth_of`, now on edits): the one copy whose alignments hold the fewest edits, strictly fewer
  than every alignment on another copy or off every copy. Edits are whole numbers, so that is at least one edit fewer.
- **Statuses:**
  - **determined, refcons:** a reference-haplotype alignment holds the fewest edits;
  - **determined, altonly:** only an ALT one does;
  - **undetermined:** two copies, or a copy and an off-copy locus, share the fewest edits;
  - **off-copy:** the fewest edits lie off every reference copy.
- **S** is unchanged: held-out, determined, the truth among the read's tied copies.

**The arms** (replace the section above). Each gives one copy per read; a shared best keeps minimap2's pick.
1. **minimap2:** its primary.
2. **O2, fewest edits (the main arm):**
   - each tied copy is represented by its highest-AS alignment (the corrector's choice);
   - its edits are its substituted, inserted and deleted bases (`o2correct.alignment_events`, indels shifted left), counted
     over the read bases every tied alignment aligns (`o2correct.common_window`);
   - the unique fewest takes the read.

   Errors and alleles every tied copy shares cancel, so the difference is the PSVs.
3. **Fewest substitutions:** `o2correct.py run` (X ops, or NM less its indel bases, over the whole alignment).
4. **Highest AS (negative control):** the unique highest AS among the tied copies.
5. **Corrector:** `o2correct.py donor`, bfa15f12. A comparison only, with no bar.

**Strata** (one added): **structure**, i.e. whether the tied copies' alignments (each copy's highest-AS one) differ in
intron count ("intron-differ") or not ("same structure").

**Bars** (replace the section above). On the held-out half, O2 (arm 2) against minimap2:
- **P:** |S| >= 300, from >= 11 families, with minimap2 wrong on >= 50.
- **G:** accuracy >= minimap2's + 2.0 points.
- **S:** fixes >= 5 × breaks.
- **F:** >= 11 families where it fixes more than it breaks, and no more than a fifth as many (rounded up) where it breaks
  more than it fixes.
- **H:** on refcons reads, accuracy >= minimap2's − 0.2 points.

**Verdict:** PASS = all five hold; FAIL = P holds and another bar fails; UNDERPOWERED = P fails. This one verdict replaces
V1 and V2. Arms 3-5 are reported with the same numbers, without verdicts.

**Development and freeze** are unchanged: changes go into Amendment 2, at the freeze commit. Fixed from here on:
- the edit rule of the truth;
- the domain;
- the arms' definitions;
- the bars.

### Amendment 2 (2026-10-20): the development half, one secondary arm, and the freeze. Written before the held-out half was scored

**The development half, scored under Amendment 1** (`personal/score_dev.json`, 140 families):
- Domain statuses: determined refcons 1,805; determined altonly 1,661; truth outside the tie 7; off-copy 1; undetermined 10.
- Scored reads S: 3,466 from 31 families.
- The truth is robust: the runner-up copy has >= 5 more edits for 3,270 of them, 2-4 for 106, exactly 1 for 90.

| arm | accuracy | fixes / breaks against minimap2 |
|---|---|---|
| minimap2 | 0.9833 (58 wrong) | – |
| O2, fewest edits | 0.9645 | 4 / 69 |
| fewest substitutions | 0.9786 | 0 / 16 |
| highest AS (control) | 0.2726 | 4 / 2,467 |
| corrector | 0.9815 | 0 / 6 |

**What it showed:**
- **AS is the wrong yardstick between copies.** 3,220 of the 3,466 scored reads are ties between copies whose alignments
  differ in intron count. minimap2 is right on 99.97% of them; the highest-AS control is right on 24.7%.
- **minimap2's 58 errors are almost all beyond any allele-blind rule.** For 48 of them the true copy has more edits than
  minimap2's pick against the reference, and for 6 it has as many. Those reads carry his own alleles. Only 4 have the true
  copy with fewer edits; O2 fixes those 4.
- **O2's 69 breaks are one family's reads, not a bug.** 65 are SM251 (chr1, 3 copies at [0.93, 0.98]):
  - The reads match his maternal SM251.0 with 1-2 edits, inside SM251.0's ALT ortholog (MAPQ 60).
  - That maternal copy differs from the reference's own SM251.0 by about 41 edits per 2.5-kb read, so against the reference
    the reads are closer to SM251.1 (34 edits against 43).
  - It is a copy the reference effectively lacks (O3's case), which no allele-blind rule can place. minimap2's pick of
    SM251.0 is right.

**Changes:**
1. **Fewest substitutions ran on a smaller BAM.** `o2correct.py run` on `sub.bam` exceeded the 10-minute limit (it
   attributes all 4.06M records). It was run unchanged on `domain.bam`, every alignment of the 4,380 domain reads cut from
   `sub.bam`. All 4,380 reads got the same tie sets and primaries as `arm2.tsv`.
2. **Truth output.** The printout shows the development half only. A descriptive `runner_up` column was added (the fewest
   edits on any other copy or locus).
3. **A secondary arm, designed after the development result: O2, consistent.**
   - It takes O2's pick only where the read's bases contradict it nowhere in the common window. A contradiction is a read
     base that O2's copy mismatches (an X in its `--eqx` alignment) and that a rival tied copy's alignment matches (an =).
   - Otherwise it keeps minimap2's pick.
   - This is the O2 regimen's own principle, "its own bases contradict the call nowhere" (register 1342). In SM251 the
     maternal copy is a mosaic of SM251.0 and SM251.1 bases, so it vetoes those moves.
   - It gets the same five bars and its own verdict, secondary to O2's.
   - Development: 2 overrules, 0 fixes / 1 break, accuracy 0.9830. It vetoes O2's 4 fixes too.

Unchanged: O2's main arm, the truth's edit rule, the domain, the bars, and the primary verdict (O2, fewest edits, against
minimap2).

**Freeze.** The held-out half is scored once with the code at the commit that adds this amendment. Before scoring, every
output is rebuilt with that code and compared with the development run: `truth_edits.tsv`, `arms.tsv` and `run.tsv`.

## Results: the held-out half, scored once at 00edc953 (`personal/score_heldout.json`)

**Domain:** 896 held-out reads tied within 10 AS. 887 are determined, all with the truth among the tied copies (refcons 479,
altonly 408); 3 are undetermined, 6 off-copy.

**S:** 887 reads from 35 families; minimap2 is wrong on 91, so P holds.

| arm | accuracy | fixes / breaks against minimap2 | families fixed > broken / broken > fixed |
|---|---|---|---|
| minimap2 | 0.8974 | – | – |
| **O2, fewest edits (primary)** | **0.8985** | **2 / 1** | **2 / 1** |
| O2, consistent (secondary) | 0.8985 | 1 / 0 | 1 / 0 |
| fewest substitutions | 0.8884 | 0 / 8 | 0 / 1 |
| highest AS (control) | 0.6077 | 68 / 325 | 1 / 22 |
| corrector (descriptive, no bar) | 0.9899 | 82 / 0 | 2 / 0 |

**O2's bars:**
- P holds;
- **G fails:** +0.11 points, against +2.0;
- **S fails:** 2 fixes against 5 × 1 break;
- **F fails:** 2 families, against 11;
- H holds: −0.21 points.

**Verdict: FAIL by the pre-registered bars.** O2 changed 3 of minimap2's 887 calls: SM556 and SM456 fixed, SM463 broken.
All three are noisy reads, with 1.8-4.6% edits against a median of 0.29%. Substantively that is **no measurable effect**.
O2, consistent: FAIL (G, F).

**Two caveats on the bars** (audit):
- **No arm could have passed F.** minimap2's 91 errors lie in 5 families (SM392 82, SM502 6, SM556, SM456, SM277 1 each), so
  even an oracle would be better in 5 families against the 11 F needs. On the development half an oracle would gain +1.67
  points, short of G. Amendment 1 dropped the oracle arm and its gate, which would have shown this before the freeze. It is a
  design flaw of the bars, as in register 950.
- **P holds only through SM392.** Without it, S is 805 reads with 9 minimap2 errors: UNDERPOWERED.

**Where minimap2 errs:**
- **SM392 (82 errors):** two copies 4 kb apart on chr5 (0.989 identity), on a chromosome where the reference is his mother's.
  - Every scored read comes from his paternal SM392.1, whose ortholog is sound: order and spacing are kept (4,179 vs
    4,172 bp), and 0.913 of the reference SM392.1's 31-mers are in it. It differs from the reference's SM392.1 by a median 21
    edits per read.
  - On the reference, every one of the 82 primaries is a **two-copy chimera**: about 2.0 kb on SM392.0, a 12,220-bp "intron",
    then 475 bp on SM392.1. `summary.attribute` credits it to SM392.0.
  - minimap2's own AS favours SM392.1 on 67 of the 82, and the highest-AS control fixes those 67.
  - O2 does not fix them. Its edit count treats the 12-kb gap as a free intron, and the chimera stands for SM392.0.
  - So these reads are not beyond an allele-blind rule: AS, which charges for the gap, places most of them.
- **SM502 (6 errors):** fixed by the corrector only.
- **The other 3 errors:** one read in each of 3 families.
- **By origin:** 89 of the 91 errors are altonly reads (minimap2 0.7819 on them), 2 are refcons (0.9958).

**Structure of the ties.** 841 of the 887 scored reads tie copies whose alignments differ in intron count. Ranking by AS is
wrong on 39.2% of the scored reads (development: 72.7%), yet right on SM392's chimeras. Exact ties (AS margin 0): 29 reads,
right for every arm.

**The refcons stratum is close to self-scoring** (audit). On all 479 held-out refcons reads (and 1,805/1,805 in
development), the truth equals the copy with the fewest edits among the library's own reference alignments, which is what O2
minimises. minimap2's 0.9958 there, and H, hold largely by construction. The two assemblies add information on about 96 of
the 887 scored reads: SM392's 82 and 14 ties.

**How many reads his genome cannot place.** 3 of the 896 domain reads tie exactly under the strict-minimum rule (development:
10 of 3,484):
- only one (SM491) is a tie between two catalogued copies;
- the SM357 read ties with SM357.0's own maternal ortholog, which the ortholog lift missed (the maternal assembly carries
  0.95 of SM357.0's 31-mers at CM054600.2 ~5.07-5.11 Mb); 5 of the 6 off-copy reads are the same gap;
- the SM339 read ties with an uncatalogued paternal locus.

Another 24 scored reads are decided by a single edit and 12 by two. In register 1338, one-edit margins were 46.5% wrong on
HG002. Requiring a margin of at least 2 or at least 5 leaves the verdict unchanged (O2 2 / 1; corrector 82 / 0 and 76 / 0).

**The corrector is descriptive only.** It has no bar after Amendment 1, and it learns the sample's alleles, which is outside
the user's allele-blind definition of O2. 76 of its 82 fixes are SM392. Its 0 breaks depend on clipped bases counting as
edits: with an NM-only truth it is 78 / 3.

**The reads.** The domain reads come from four runs of the same animal: SRR27178662 (147), SRR27178663 (306), SRR27438212
(3,804) and SRR27438213 (123), all SAMN04003007 = KB3781. The substrate section names only the first.

**Reading.** By the pre-registered bars, O2's fewest-edits rule FAILS on held-out reads, with no measurable effect: 3 calls
changed in 887. Three things follow:
- For reads whose copy matches the reference, minimap2's pick is already the fewest-edit copy almost always, so an
  allele-blind rule has little to change there.
- minimap2's held-out errors are mostly one family, where its primary is a chimera spanning two copies. AS (which charges for
  the spurious intron) places most of them; the edit count, which makes every intron free, does not.
- The bars could not be passed by any arm (F), and the power bar holds only through that family.

## Audit (metric-auditor, 2026-10-20): holds with caveats

- **Every number reproduces to the digit**, from an independent rebuild of the truth (4,380 / 4,380 statuses, truths and
  origins agree) and an independent scorer. FAIL is the pre-registered verdict.
- **The hold-out order is clean:** prereg 725aafb2 → Amendment 1 99595983 → development score 19:00 → freeze 00edc953 at
  19:01:11 → held-out score 19:01:57. The frozen code matches 00edc953, and the frozen TSVs are byte-identical to the
  development-run files.
- **The secondary arm came from development only.** One point the auditor could not settle: whether any run printed held-out
  statuses under the edit truth before the freeze. The session shows the edit truth was first run with the
  development-only printout, so they were first seen at the held-out scoring (counts only in any case).
- **Corrections taken into the Results above:**
  - F was out of reach for any arm, and P holds only through SM392;
  - SM392's minimap2 errors are two-copy chimeric primaries that AS places, so not "beyond any allele-blind rule" and not
    O3's case;
  - "No allele-blind arm fixes any altonly read" was false: highest AS is +67 / −93 on altonly reads;
  - refcons is close to self-scoring;
  - "3 truly ambiguous" became 3 exact ties under the strict-minimum rule, of which 1 is a two-copy tie, plus 24 one-edit
    decisions;
  - the corrector's 0 breaks depend on counting clips;
  - four sequencing runs, not one.
- **Wording not to use:**
  - "O2 is refuted on real reads";
  - "P holds, so the test was powered";
  - "minimap2's errors are beyond any allele-blind rule";
  - "only the corrector fixed them";
  - "exact truth" without qualification;
  - any pass wording for the corrector.
- **Unclear, low risk:** the splice indexes were built before `gorilla_haps` was downloaded. Contigs map by unique length (24
  / 24 paternal, 223 / 225 maternal, no alignments on the 2 ambiguous ones), but sequence identity is not proven.
- The auditor's scripts are in `/mnt/linuxdisk/tmp/audit_2026-10-20_o2personal/`; the frozen copy is
  `personal/audit_frozen_00edc953/`.
