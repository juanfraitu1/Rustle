V5 = V4 + RETIREMENT OF THE SUPERSEDED PLACEHOLDER — FULL REPORT (agent 1 of 2)
Working dir /mnt/linuxdisk/home/juanfraitu/v5_retire. Nothing in src/ modified, nothing committed, no
subagent spawned, no MAPQ gate lowered, nothing written under ggo_npip/ or v4_gorilla/ (verified by
mtime scan at the end). All captured PAFs, node tables and replay wrappers reused; only absent query
keys were mapped (chimp: 86 tx + 77 body, two foreground minimap2 calls, rc=0 both).

=====================================================================================
0. PRE-REGISTRATION TRAIL (every file written before the numbers it governs)
=====================================================================================
  DECLARATIONS.txt                      md5 17c96d35cc5386c08f407f2a07159a05   11:47:24
  DECLARATIONS_ADDENDUM_A.txt           md5 702e3e11efdbb5edfca2430b2a00d85b   11:53:59  (base recovery)
  CHIMP_NOISE_FLOOR_FROZEN.txt          md5 d5b5d376e56f8bceecf7260dc2919f3e   11:54:57
  DECLARATIONS_ADDENDUM_B.txt           md5 da8ed67a5624f184df7ae295e24ad2af   12:03:53  (V5b)
  F_CHOICE_STATED_BEFORE_CHIMP.txt      md5 3e504e69a14a26bffffdf6eeee049393   12:08:55
  v4_gorilla/NOISE_FLOOR_FROZEN.txt     md5 69af973a85f582f9759a9bcfa634feda   re-verified UNCHANGED
All six md5s re-verified at the end of the run. The gorilla thresholds were never re-derived.
The decision rule was written verbatim as given, in DECLARATIONS section 3, before any number.

=====================================================================================
1. ANSWER
=====================================================================================
V5 IS A NO-OP. At f = 1.00, 0.90 and 0.50, on gorilla AND on chimp, V5 retires ZERO nodes and is
identical to V4 key for key — same node set, same edge set, same families. Every clause it passes or
fails, it passes or fails as V4. It is NOT A CANDIDATE.

The reason is structural, and it was measured rather than assumed. V4 relabels every unmeasured base
node 'U', and 'U' blocks a read-locus candidate of EITHER strand; only a MEASURED blocker is narrowed
to its own strand. A node V4 installs can therefore never overlap an unmeasured base node, and
containment requires overlap.

    gorilla   86 V4-added nodes, 86 overlap a base node exonically.
              Overlapped base nodes: 98 MEASURED OPPOSITE-STRAND, 0 unmeasured, 0 measured same-strand.
    chimp     86 V4-added nodes, 86 overlap a base node exonically.
              Overlapped base nodes: 96 MEASURED OPPOSITE-STRAND, 0 unmeasured, 0 measured same-strand.

The declared placeholder set is unreachable by construction. The tie rule never fires (0 abstentions
everywhere), so the f sweep is entirely flat.

WHAT THE FRAGMENTS ACTUALLY ARE. Opening the flagged families shows the superseded fragments ARE
contained by the node V4 installs — gorilla containments 1.000, 1.000, 1.000, 1.000, 1.000, 0.937,
0.652, 0.563 — but every one of them is MEASURED under V4's own predicate, because V4 counts a
junction-bearing read of ANY strand as measuring a node's strand. At exactly these loci the junction
reads that mark the '+' fragment "measured" are the '-' reads that form the new '-' node. The
fragment's own strand was never measured by anything; it is marked measured by evidence for the
opposite strand. The task's English ("a single-exon node whose strand was never measured") describes
these fragments. V4's predicate does not. That mismatch is the whole result.

=====================================================================================
2. THE ONE-LINE REPAIR, PRE-REGISTERED AND THEN REFUTED
=====================================================================================
V5b (ADDENDUM B, written before any V5b number and before any chimp variant number) is V4 with its
node construction COMPLETELY unchanged plus retirement using one substituted predicate:
    V5   placeholder = single-exon base node with NO junction read overlapping it.
    V5b  placeholder = single-exon base node with no junction read OF ITS OWN STRAND overlapping it.
Everything else — retainers, containment, the f grid, the tie rule, the read attribution, every
clause and threshold — identical. V5b is EXPLORATORY: gorilla is development for it, chimp is its
first and only held-out test, and Addendum B says in advance that one held-out pass would not be
adoption evidence.

  gorilla   placeholders 297 -> 675; retires 41 (f=1.00) / 43 (0.90) / 47 (0.50); 0 abstentions.
            Repairs what it can reach: V0 family 35's scatter disappears, four of V0 family 2's NPIP
            fragments are retired into retainers inside family 2. Still FAILS (b) and (c):
            lost 2, scatter 3, membership 73/76, FAMILY R 0.64. NOT A CANDIDATE.
  chimp     retires 31; STRICTLY WORSE THAN V4: scatter 1 -> 3, membership 41/43 -> 39/43,
            FAMILY R unchanged at 0.4375. It creates two splits V4 did not have — V0 family 9 splits
            across V5b families 10 and 18, V0 family 15 across 10 and 28 — because the retired node
            was holding each family's triangle together.
V5b is refuted on the held-out substrate. Do not adopt it, and do not re-tune it on gorilla.

=====================================================================================
3. THE FULL TABLE (SCORER v3; both membership readings; both precision conventions)
=====================================================================================
GORILLA (development). 25 membership loci, 3,112 expressed records, 76 V0 families.
Frozen thresholds: lost<=0, scatter<=0, membership>=76/76, orphaned==0.

 variant   nodes    (-,+)   opp/25  lost scat orph  memb(aware/strict) FAM_R  junk   R      P_cand P_old  F_old
 V0        2664   (-0,+0)    9       0    0    0     76/76 / 76/76     0.68    464  0.9235 0.8291 0.6847 0.7864
 V4        2750   (-0,+86)   3       2    3    0     71/76 / 71/76     0.64    482  0.9337 0.8131 0.6705 0.7805
 V5_100    2750   (-0,+86)   3       2    3    0     71/76 / 71/76     0.64    482  0.9337 0.8131 0.6705 0.7805
 V5_090    2750   (-0,+86)   3       2    3    0     71/76 / 71/76     0.64    482  0.9337 0.8131 0.6705 0.7805
 V5_050    2750   (-0,+86)   3       2    3    0     71/76 / 71/76     0.64    482  0.9337 0.8131 0.6705 0.7805
 V5b_100   2709  (-41,+86)   3       2    3    0     73/76 / 69/76     0.64    470  0.9332 0.8231 0.6803 0.7869
 V5b_090   2707  (-43,+86)   3       2    3    0     73/76 / 68/76     0.64    469  0.9327 0.8231 0.6805 0.7868
 V5b_050   2703  (-47,+86)   3       2    3    0     73/76 / 68/76     0.64    466  0.9322 0.8230 0.6811 0.7871

CHIMP (decision substrate, untouched by any strand variant). 16 membership loci, 842 expressed
records, 43 V0 families. Frozen thresholds: lost<=0, scatter<=0, membership>=43/43, orphaned==0.

 variant   nodes    (-,+)   opp/16  lost scat orph  memb(aware/strict) FAM_R  junk   R      P_cand P_old  F_old
 V0        1589   (-0,+0)    2       0    0    0     43/43 / 43/43     0.5625  603  0.8836 0.7546 0.4682 0.6121
 V4        1675   (-0,+86)   1       3    1    0     41/43 / 41/43     0.4375  633  0.8931 0.7217 0.4490 0.5975
 V5_100    1675   (-0,+86)   1       3    1    0     41/43 / 41/43     0.4375  633  0.8931 0.7217 0.4490 0.5975
 V5b_100   1644  (-31,+86)   1       3    3    0     39/43 / 37/43     0.4375  615  0.8919 0.7298 0.4568 0.6042

CLAUSE VERDICTS (a)(b)(c)(d). Every variant on every substrate: (a) PASS, (d) PASS, (b) FAIL,
(c) FAIL. No arm produces a CANDIDATE, so no CONFIRM is reachable and none is claimed.

Two reporting notes fixed in advance, not after the numbers. P_cand is conditioned on the prediction
and is reported, never gated (metric trap register). P_old rises for V5b only because the denominator
shrank — the node count is printed beside it every time and no V5b claim rests on it. MEMBERSHIP-
STRICT cannot be passed by any node-removing variant whose removals touch a family; both readings are
printed together every time, and the threshold 76/76 (43/43 on chimp) was not moved.

=====================================================================================
4. CLAUSE (d) — EVERY REMOVED NODE ACCOUNTED FOR
=====================================================================================
 variant     removed  all single-exon  all contained>=f  overlapped an expressed  best node of a
                                                          annotated record        membership locus
 V5 (all f)      0         yes               yes                 0                 V0 0, V4 0
 gorilla V5b_100 41        yes               yes                29                 V0 5, V4 0
 gorilla V5b_090 43        yes               yes                30                 V0 5, V4 0
 gorilla V5b_050 47        yes               yes                31                 V0 5, V4 0
 chimp   V5b_100 31        yes               yes                13                 V0 1, V4 0
Abstentions (>= 2 containing retainers): 0 in every arm at every f. The declared tie rule never fired.
gorilla V5b_100: 94,965 exonic bp retired, 81 reads reattributed. chimp: 85,132 bp, 43 reads.

The five gorilla nodes that were the V0 best node of a membership locus (LOC129527636, LOC109023568,
LOC115932744, LOC129527692, LOC115933039) were all already superseded in V4 — 0 of them was V4's best
node for its locus — so retiring them removes nothing the fix had not already replaced. Same on chimp
for LOC112205831. f = 0.50 is the one setting that loses real sequence: it retires
NC_073244.2:23011904-23016922 at containment 0.5175, overlapping MED26. That is hazard H9 realised,
and it is why f = 1.00 was carried forward.

=====================================================================================
5. EVERY FLAGGED FAMILY, OPENED — FIX OR GENUINE SPLIT
=====================================================================================
GORILLA, V4 (8 of 76 families flagged):
  fam 52 (2 nodes, J 0.0705, LOST+SCATTER) — FIX DISPLACING A PLACEHOLDER. The '+' fragment
     NC_073244.2:20955338-20957550 at LOC115933039, contained 1.000 by the new 45-exon '-' node,
     ends in no family at all.
  fam 35 (2 nodes, J 0.5599, SCATTER) — FIX DISPLACING A PLACEHOLDER. NC_073241.2:26331092-26333296,
     contained 1.000 by the new 19-exon '-' node, ends in no family.
  fam 51 (2 nodes, J 0.0876, LOST) — FIX-DRIVEN MERGE. Both '+' fragments are absorbed into the
     28-node family 2. Nothing splits; the Jaccard collapses because a 2-node family was swallowed.
  fam 3 (12, J 0.6000) and fam 6 (9, J 0.7239) — FOOTPRINT DILUTION, not splits. Every node stays in
     one family; J falls only because the newly added nodes widen that family's footprint.
  fam 2 (21 nodes, J 0.5788, SCATTER) — GENUINE SPLIT. 19 nodes stay in family 2; two break off into
     family 55: NC_073242.2:99634383-99636626 (LOC129527693) and NC_073242.2:103962940-103964751
     (LOC101134557). Neither is contained by any added node — there is no new node at those loci at
     all — so no containment rule, at any f, can reach this split.
  fam 4 (11 nodes, J 0.6315, SCATTER) — GENUINE RE-PARTITION, 3 to family 3 and 8 to family 5, driven
     by two long new '-' nodes bridging previously separate '+' fragment clusters.
  fam 11 (4 nodes, J 0.7719, SCATTER) — GENUINE SPLIT, 2/2 across families 45 and 7. The pair that
     leaves is contained only 0.563 and 0.652, so only f = 0.50 could touch it.
GORILLA, V5b_100 (8 flagged): fam 35's scatter is REPAIRED (fragment and retainer now in one family)
  and four of fam 2's NPIP fragments are retired into retainers inside family 2. Fam 2 still scatters
  because the LOC129527692 fragment's retainer sits in family 4 and family 55 still holds the two
  uncontained fragments. Fams 51 and 52 stay LOST: retirement moves their fragments into family 2,
  which is a merge the Jaccard still charges for.
CHIMP, V4 (7 of 43 flagged):
  fam 32 (J 0.108), fam 7 (J 0.162), fam 2 (J 0.594) — FIX-DRIVEN MERGE into the enlarged family 1.
  fam 23 (J 0.518), fam 3 (J 0.621) — FOOTPRINT DILUTION, every node stays.
  fam 37 (J 0.414, LOST+SCATTER) — FIX DISPLACING A PLACEHOLDER: NC_072416.2:44223894-44226204 falls
     out of every family once the new 7-exon '-' node takes the locus.
  fam 0 (42 nodes, J 0.973, SCATTER) — GENUINE SPLIT: the 16-exon '-' node at LOC112207901 relocates
     from family 0 to family 3 alongside the new 34-exon node at LOC112205831. A real spliced-node
     move, not a fragment.
CHIMP, V5b_100: the two NEW splits are V0 family 9 (4 nodes -> families 10 and 18, after
  NC_072416.2:43899999-43902241 is retired) and V0 family 15 (3 nodes -> families 10 and 28). This is
  the held-out damage, and it is the finding that decides against V5b.

=====================================================================================
6. THE CHIMP NULL, FROZEN BEFORE ANY CHIMP VARIANT NUMBER
=====================================================================================
30 relabellings of the unchanged chimp V0 graph (1,589 nodes, 918 edges, 43 families, 175 loci,
largest 42), seed base 20260918 — the same base as gorilla, so the two nulls are comparable. All five
per-draw assertions (non-identity permutation, edge set maps back, edge count, degree multiset,
n_reads multiset) passed in all 30. Result: lost {0:30}, scatter {0:30}, orphaned {0:30}, membership
{43/43 : 30}, mean and min best Jaccard 1.0000 in all 30. ZERO WIDTH, exactly like gorilla.
Consequence, declared in the frozen file rather than after the verdict: a clause calibrated to a
zero-width null is the clause "change nothing", and it charges any variant that adds 86 nodes for
changing the partition at all. The clause was still applied exactly as frozen, and every family it
flagged was opened.

=====================================================================================
7. GATES AND PROVENANCE
=====================================================================================
 C1 chimp replay: QUERY_PARITY tx IDENTICAL, body IDENTICAL; NO minimap2 run (captured PAF replayed);
    the replayed binary's own log reproduces 5048 skeletons -> 1483 reps, 1455 gene-level + 134
    read-locus = 1589 nodes, 558 exon + 859 body (918 pairs), 43 families / 175 loci.   PASS
 C2 chimp split 134 / 1455 — required ADDENDUM A (see below).                            PASS
 C3 chimp V0 re-derived exactly as a multiset of (chrom,strand,exons), no widening.      PASS
 C4 chimp V0 graph 558 + 859 / 918 / 43 / 175, re-verified after the new PAFs merged in. PASS
 C5 16 of 16 chimp membership loci resolve in PTR_genomic.gff.                           PASS
 P4 gorilla V0 reproduces 1194 + 297 / 1248 / 76 / 323.                                  PASS
 Parent parity: gorilla V0 and V4 node key-multisets IDENTICAL to the v4_gorilla run;
    V4 2,750 nodes, 1,453 pairs, 77 families, adds 86, opposite 9 -> 3.                  PASS
ADDENDUM A was forced: the parent run's raw-rep-span criterion gives 135 / 1454 on chimp and fails
C2 by one node, because consolidation merges overlapping same-strand reps before nodes are emitted
(chimp 1483 reps -> 1455 loci). Grouping reps by span overlap first gives exactly 1,455 groups
matching the Rust log, 134 non-matching nodes, all 134 read-locus candidates, zero ambiguity.
Cross-checked on gorilla: 2,469 -> 2,431 groups, 233 / 2431, and the recovered base set is IDENTICAL
key for key to the one the parent V4 run used. The amendment changes no gorilla number.

=====================================================================================
8. RECOMMENDATION
=====================================================================================
Do not adopt V5: it changes nothing, anywhere, at any f. Do not adopt V5b: its single held-out test
is negative. Record both as closed, with the structural reason, so neither is re-proposed.

The V4 family failure is not a redundant-containment problem. Opening all fifteen flagged families
across the two substrates shows it is two different things, and node retirement addresses neither:
  (i) FIX-DRIVEN MERGES — small 2-to-4-node families of '+' fragments being absorbed into the large
      family once the correct spliced '-' node appears (gorilla fams 51, 52; chimp fams 32, 7, 2).
      These look CORRECT. A zero-width Jaccard clause scores a correct merge as a total loss.
  (ii) GENUINE SPLITS at loci where no new node exists at all (gorilla family 55's pair; chimp family
      0's LOC112207901). Nothing contains these, so nothing can retire them.
The defensible next experiment is therefore on the CLAUSE, not the nodes: pre-register a family
statistic that distinguishes a V0 family entirely contained in one variant family (a merge) from one
whose nodes land in two or more (a split), re-run V4 against it on gorilla, and hold chimp back. If
V4 passes that clause, the strand fix is adoptable on the evidence already in hand — on both
substrates it improves the membership axis (gorilla 9 -> 3 opposite, chimp 2 -> 1, one real GAIN
each) and removes no node.
Worth recording separately as a measured defect in its own right, independent of V5: V4's MEASURED
predicate counts a junction read of ANY strand as measuring a node's strand. It classifies 863 - 297
= 566 gorilla single-exon base nodes and 800 - 323 = 477 chimp ones as "measured" on the strength of
opposite-strand junctions alone. That is what makes V5 unreachable. Fixing the predicate and then
retiring (V5b) is refuted; whether the predicate should be fixed for the SUPPRESSION test itself is a
separate, untested question and must not be inferred from this run
> **Verifier corrections, applied by the orchestrator** (verifier ok = true; the negative verdict stands):
> 1. **"8 of 76 families flagged" overstates it** — under the clause's own failure conditions exactly **six** gorilla families are flagged by V4 (2, 4, 11, 35, 51, 52); families 3 and 6 fail only the Jaccard reading while keeping every node.
> 2. **V5b strict membership is not 69/76 at every f**: 69/76 at f = 1.00, 68/76 at 0.90 and 0.50.
> 3. **"f = 0.50 adds 6" is measured against f = 1.00**, not against the preceding grid point (41 → 43 → 47).
> 4. **The chimp clause-(a) gain was truncated**: LOC112205831 goes 0.000 → **0.9812** same-strand coverage — the single largest per-copy strand repair measured anywhere in this line.
.
## Verification (independent recompute) — agent 2 of 2

**Scope.** Everything below was recomputed from my own code under `/mnt/linuxdisk/home/juanfraitu/v5_retire/verify/` (`v1_retire.py` retirement, `v2_eval.py` family panel / FAMILY R, `v3_ptr_truth.py`, `v4_ptr_base.py`, `v4a_reads.py`, `v5_ptr_build.py`, `vmm2_ptr.sh`), re-using only my own earlier verifier modules (`v4_gorilla/verify/w*.py`, `v4_heldout/verify/x*.py`, `strand_fix/verify/vmirror.py`). **No builder script under `v5_retire/` was read or imported.** I did my own chimp replay of `gw_family_catalog.bin` and my own minimap2 runs. Verdict: **the builder's report reproduces, number for number, on both substrates, and the recommendation follows.**

### 1. Declaration ordering (mtimes)
`DECLARATIONS.txt` 11:47:24 (md5 `17c96d35…`, unchanged) precedes every result file. `ADDENDUM_A` 11:53:59 precedes the chimp V0 graph 11:54:22; `CHIMP_NOISE_FLOOR_FROZEN.txt` 11:54:57 (md5 `d5b5d376…`, unchanged at the end) precedes every chimp *variant* artifact (earliest 12:09:27). `F_CHOICE_STATED_BEFORE_CHIMP.txt` 12:08:55 precedes them all. `ADDENDUM_B` 12:03:53 precedes the chimp arm entirely. `v4_gorilla/NOISE_FLOOR_FROZEN.txt` still md5 `69af973a85f582f9759a9bcfa634feda`. **One audit limitation, not a violation:** `out/ggo_variants.json`, which carries the gorilla V5 numbers Addendum B quotes, was last written 12:06:09, i.e. *after* Addendum B — the file was evidently regenerated when V5b was added (`g1_variants.py` mtime 12:04:21). The ordering claim for those V5 numbers therefore rests on the declaration text, not on mtimes. It does not matter for the verdict because I rebuilt them myself.

### 2. Gorilla: V4 and V5(f) rebuilt from scratch
Parent parity (my own artifacts, independent of this run): V0 **2,664** nodes / 1,194 exon + 297 body = **1,248** pairs / **76** families / 323 loci; V4 **2,750** nodes / 1,370 + 457 = **1,453** pairs / **77** families, +86 / −0. Matches.

V5(f), implemented from DECLARATIONS §1 by me: **59 retainers, 297 placeholders, 0 RETIRED, 0 ABSTAIN at f = 1.00, 0.90 and 0.50**, and the V5 node set, pair set and family partition are **equal key-for-key to V4** (checked explicitly). The 297 placeholders are exactly the V4-'U' set intersected with single-exon (`U = 297`; the junction-only predicate gives 304, of which the 7 extra are non-single-exon — the builder's "dropped-nonsingle 7").

**The "why" is structural, and I verified it twice.** (i) Measured: of the 86 nodes V4 adds, **all 86** overlap a base node exonically, and the overlapped base nodes are **98 MEASURED-OPPOSITE, 0 UNMEASURED, 0 MEASURED-SAME** (chimp: 86 added, **96 / 0 / 0**). (ii) Proved from the mirror source: `with_read_locus_nodes` builds its suppression index over `base` only, and `_compat(a,b,'U')` returns True whenever either side is `'U'`, so a 'U' base node blocks a candidate of *either* strand. Any candidate V4 installs therefore overlaps no unmeasured base node; containment implies overlap; so **V5 ≡ V4 is a theorem given V4's definition, not a data accident.** The builder's claim is correct as stated.

### 3. Gorilla: V5b (exploratory, Addendum B)
Reproduced exactly: placeholder set 297 → **675**; retired **41 / 43 / 47** at f = 1.00 / 0.90 / 0.50; **0 abstentions at every f** (the tie rule never fires on either substrate); **0** retirements with a same-strand *label* retainer. Nodes 2,709 / 2,707 / 2,703.

**Clause (d), re-verified node by node, not taken on trust.** For every removed node I independently re-tested (i) `len(exons)==1 and len(rep_exons)==1`, (ii) membership in the V5b-unmeasured set recomputed from the reads, (iii) that the node is a *base* node, and (iv) containment ≥ f. **0 violations at every f on both substrates.** Accounting at f=1.00: 41 removed, **94,965 exonic bp**, **81 reads reattributed**, **29/41** overlap an expressed annotated record, **5** were the best node of a membership locus under V0 (LOC129527636, LOC109023568, LOC115932744, LOC129527692, LOC115933039) and **0** under V4 — i.e. all five had already been superseded. Confirmed.

### 4. Gorilla clauses and flagged families (opened myself)
Scorer v3, my recompute: V0 junk 464, M 1824, C 2200, R 0.9235, P_cand 0.8291, P_old 0.6847, F_old 0.7864 · V4/V5 junk 482, M 1844, R 0.9337, P_cand 0.8131, P_old 0.6705, F_old 0.7805 · V5b_100 junk 470, M 1843, C 2239, R 0.9332, P_cand 0.8231, P_old 0.6803 @2,709 nodes, F_old 0.7869. Every figure matches `VERDICT.json`.

Clauses: (a) PASS for V4/V5/V5b, opposite-strand membership loci **9/25 → 3/25**, one GAIN on the ss_frac>0.5 criterion (LOC115933039 0.0910 → 0.7331), 0 LOSS, 0 carried. (b) FAIL — lost 2, scatter 3, membership 71/76 (V5) and 73/76 retire-aware / 69/76 strict (V5b) against the frozen 0 / 0 / 76/76. (c) FAIL — membership-locus FAMILY R **0.68 → 0.64** for every variant. (d) PASS. **V5 and V5b are not candidates.** Confirmed.

Opening the families myself reproduces the builder's account: fam 52 (J 0.0705) and fam 35 (J 0.5599) each lose a '+' fragment that ends in no family and is contained **1.000** by the new '-' node (45-exon and 19-exon respectively); fam 51 (J 0.0876) is a clean absorption of both '+' fragments into the 28-node family 2; fam 4 re-partitions 3/8 and fam 11 splits 2/2. **The key negative reproduces:** V0 family 2's breakaway pair is NC_073242.2:99634383-99636626 and NC_073242.2:103962940-103964751, and I confirm **zero V4-added nodes overlap either locus** — no containment rule can reach them. Under V5b, fam 35's scatter is repaired and family 2 gets *worse*, scattering over [2, 4, 55] because the LOC129527692 fragment's retainer (NC_073242.2:99211161-99344539) sits in family 4. All as reported.

### 5. Chimp: independent provenance
I replayed `gw_family_catalog.bin` myself (`RUSTLE_SHARED_DEFINITION=1 RUSTLE_LOCUS_AUDIT=1`, my own no-minimap2 wrapper, rc 0, 4:13 wall, 1.1 GB peak). **GATE C1 PASS:** tx and body query FASTAs byte-identical to `capture/ptr_dn/`; the log reproduces `5048 skeletons -> 1483 reps`, `1483 reps -> 1455 gene-level loci + 134 read-locus nodes = 1589 nodes`, `558 exon + 859 gene-body (918 pairs); 43 … families holding 175 loci`; and my `vptr.copies.tsv` / `vptr.families.tsv` are **md5-identical to the shipped `ptr_sd.*`**. **GATE C2 PASS:** 1,483 reps → **1,455** overlap-merged groups → **134 / 1,455**, all 134 being read-locus candidates. Addendum A checks out in detail: the raw-span reading (node `exons`-column span vs raw rep spans) gives 1,454 / 135 of which only 134 are candidates — the off-by-one the addendum describes; the gorilla cross-check gives 2,469 → **2,431** groups, 233 / 2,431, **base set identical key-for-key** to the parent V4 run's. **GATE C3 PASS:** V0 reproduces the published 1,589-node set with 0 missing / 0 extra. **GATE C4 PASS:** 558 + 859 = 918 pairs, 43 families, 175 loci, and the V0 partition is **identical to `ptr_sd.copies.tsv`**. **GATE C5 PASS:** 16/16 membership loci resolve in `PTR_genomic.gff` (1,582 gene+pseudogene records on NC_072416.2).

### 6. Chimp noise floor, recomputed with my own permutation code
30 relabellings, `random.Random(20260918+k)`, both scatter and gather orientations, all five per-draw assertions passing: **lost {0:30}, scatter {0:30}, orphaned {0:30}, membership {43/43: 30}, mean and min best Jaccard 1.0000 in all 30**, n_families 43, loci 175 in every draw. **Zero width**, exactly the frozen file. The same holds for the V4 graph (42 families, 189 loci, 0/0/0/42-of-42 in all 30). The consequence the builder states is correct: at this floor the clause is literally "change nothing".

### 7. Chimp variants and the held-out result
V0 1,589 nodes / 918 pairs / 43 families, opposite **2/16** (H10 ceiling: LOC112205831, LOC112207763), junk 603, R 0.8836, P_cand 0.7546, P_old 0.4682, F_old 0.6121, FAMILY R 0.5625. V4 1,675 (+86 / −0) / 1,004 pairs / 42 families, opposite **1**, junk 633, R 0.8931, P_cand 0.7217, P_old 0.4490, F_old 0.5975, FAMILY R 0.4375, lost 3 / scatter 1 / membership 41/43. **V5(1.00) is identically V4** (48 retainers, 323 placeholders, 0 retired). V5b(1.00) 1,644 nodes (−31, 0 abstentions, 0 violations of clause (d), 85,132 bp, 43 reads, 13/31 overlapping an expressed record, 1 best-of-locus in V0 and 0 in V4), 946 pairs, 42 families, junk 615, R 0.8919, P_cand 0.7298, P_old 0.4568, F_old 0.6042, FAMILY R 0.4375, **lost 3 / scatter 3 / membership 39/43 retire-aware (37/43 strict)**. Every figure matches `VERDICT.json`.

**The decisive held-out finding replicates: on the substrate no strand variant has touched, V5b is strictly worse than V4 on the family clause — scatter 1 → 3, membership 41/43 → 39/43, mean best J 0.9161 → 0.9029 — with no gain on (a) (opposite 1 in both) and none on (c) (FAMILY R 0.4375 in both).** Opening the families myself shows why, and it is retirement's own doing: V0 fam 9 scatters because the retired fragment NC_072416.2:43899999-43902241 is attributed to a retainer sitting in family 10 while its three siblings go to family 18; V0 fam 15 scatters as a knock-on of family 10 changing. Neither scatter exists under V4.

### 8. Traps
- **A "loss" that is a correct merge.** Checked explicitly and it is the dominant failure mode. On chimp, all three V4 "lost" families are absorptions, not splits: fams 7 (J 0.1620) and 32 (J 0.1084) keep **every node, all landing in one family** (fam 1, 23 nodes); the single V4 scatter is one 16-exon '-' node moving to fam 3. On gorilla, fams 51 and 52 are likewise absorptions into family 2. The clause is charging both variants for merges.
- **A gain that is denominator shrinkage.** Verified explicitly by the restriction route: the V5b node set **is** the V4 node set minus the removals (confirmed as a multiset identity), so P_old and junk are directly comparable. Gorilla: matched pairs 1844 → 1843 (−1) while the denominator drops 2,750 → 2,709; junk 482 → 470 but the junk **rate** barely moves (0.1753 → 0.1735) because only **12 of the 41** removed nodes were junk and **29** overlapped an expressed record. Chimp: M 752 → 751, junk rate 0.3779 → 0.3741, 18 of 31 junk. **Every apparent P_old/junk improvement is denominator shrinkage plus the loss of real sequence.** The builder rests no claim on it; nor do I.
- **A copy counted as recovered that was already correct.** Gorilla: 0 carried, 1 GAIN, 0 LOSS — nothing pre-correct is counted. Chimp: 8 carried (excluded from the gain), 1 GAIN, 0 LOSS. Clean.
- **A clause that cannot fail.** Clause (d) is near-tautological by construction — declared as such, and I made it a real test by re-deriving the three conditions from the reads and the base set independently (0 violations). Clause (b) cannot be *passed* by any variant that perturbs the partition, both nulls being zero-width — declared as H4. Both readings of membership are reported together, as required.
- **Retire-aware membership must not be allowed to flatter V5b.** It does flatter it: V5b's gorilla membership 71/76 → 73/76 over V4 comes **entirely** from the declared substitution (fams 35 and 52 stop failing only because the retired fragment is mapped to its retainer); strict membership is 69/76, i.e. worse than V4. Both numbers are printed by the builder every time, which is the correct handling.

### 9. Does the recommendation follow?
Yes. **V5 as pre-registered is a no-op by construction** and fails (b) and (c) on both substrates — not a candidate. **V5b** is not a gorilla candidate and is refuted on its one held-out test. The closing diagnosis — that the V4 family failure is not a redundant-containment problem, since family 2's breakaway pair has no new node at either locus — is confirmed directly. The forward suggestion (target the clause, not the nodes, because a zero-width null makes the clause "change nothing" and is charging V4 for merges) is supported by what I found in §8, and is offered as a proposal, not a result.