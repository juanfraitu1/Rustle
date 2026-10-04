# o3_candidates A13 (net by alignment, structural template) Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Fix the two measured causes of Amendment 12's failure in the `o3_candidates` stage and re-run the acceptance unchanged; flip the stage to default-on iff it passes.

**Architecture:** Changes confined to `src/bin/o3_candidates.rs` (pass B: unmapped reads attributed by one minimap2 run against the copies; template choice from the all-vs-all PAF), `src/rustle/vg_family/o3_candidates.rs` (new preset `MM2_ATTRIB`, splice `MM2_MEMBERS`, the structural-template scorer, the insertion vote order, the refine re-template, the empty-merge fallback), and the acceptance harness re-run into `/mnt/linuxdisk/tmp/rna_allele/a13/`.

**Tech Stack:** Rust (rustle crate), minimap2 2.30, the Amendment 12 harness (`bench/rna_allele/accept_o3_candidates.{sh,py}`, `panel_to_copies.py`).

**Spec:** `docs/superpowers/specs/2026-10-02-o3-candidates-design.md` as amended by prereg Amendment 13 in `docs/PREREG_rna_allele_haplotype_count_2026-10-01.md` (the amendment is the binding text for every rule below).

## Global Constraints

- Build/test only with `CARGO_TARGET_DIR=/mnt/linuxdisk/home/juanfraitu/rustle_target_m2 bash tools/rlock.sh heavy cargo ... --release`, output captured to a file; every heavy run (binary, minimap2, samtools, the harness) foreground through `bash tools/rlock.sh heavy ...`, each call < 10 min; no background jobs, no waiter loops, no `pkill -f`.
- The chain's registered rules do not move: `--delta 0.00958`, component merge rule, `--min-support 6`, `--min-cluster 3`, 0.98 tie ratio, `--max-reads 1000`, R13.
- Amendment 13b's exact values: `MM2_ATTRIB = [-x map-ont -c -N 5 -p 0.5]`; attribution set = unmapped reads >= 300 bp PLUS poorly placed un-netted reads (primary `de > 0.02` or MAPQ 0, no record on any family copy); targets = the families' mapped net reads (tagged `<family>|<read>`) + `--copies-fa`; a read joins the family of its best hit (most matches) iff the hit covers >= 50% of the READ (`(qe-qs)/qlen`) and `de <= 0.20`; template = member with the LOWEST total bases of indels >= 20 bp over its all-vs-all alignments to the other members (ties: longest, then name); `MM2_MEMBERS = [-x splice:hq -uf -c --cs -N 5 -p 0.5]`; insertion vote: >= 20 bp with >= 3 carriers first, then the < 20 bp 50% rule; refine re-templates by the same rule when the template is split off; an absorbing cluster with an empty re-polished consensus keeps the absorbed clusters separate.
- Existing tests stay green (1081 + new); the fixture integration test and `run_e2e.sh` still pass (one flagged candidate `cand_MCL0_0`, 850-950 bp).
- Commits end with `Co-Authored-By: Claude Fable 5.1 <noreply@anthropic.com>` and `Claude-Session: https://claude.ai/code/session_01DAyQQ6R8drUxY5GsM5wNkb`; author from the canonical repo's git config; never touch `/mnt/c/Users/jfris/Desktop/Rustle`; do not push.

## Review Focus

1. An unmapped read that aligns to copies of two families (a shared exon): it joins the family of its best hit only when that hit covers >= 50% of the read — a 40%-coverage best hit attributes nothing (Task 1 test).
2. A cluster whose longest member retains an intron: the template is a shorter clean member (Task 2 test `retained_intron_read_is_not_the_template`).
3. A cluster of two members with identical structure: ties resolve by length then name, deterministically (Task 2 test).
4. The fixture (all reads mapped) is unaffected by Task 1: identical outputs before/after (Task 1 check).
5. Pass B must not read sequences of unmapped records twice or hold all of them in memory when only some are long enough (Task 1: stream to the FASTA).

---

### Task 1: Net attribution by alignment (Amendment 13b: unmapped + poorly placed reads, targets = net reads + copies, map-ont, read coverage >= 0.5, de <= 0.20)

**Files:**
- Modify: `src/rustle/vg_family/o3_candidates.rs` (add `pub const MM2_ATTRIB`, `pub fn attribute_by_hits(hits: &[PafHit], family_of_target: &HashMap<String, String>) -> HashMap<String, String>` (read -> family), remove `FamilyKmerIndex`, `ATTRIB_MAX_FAMILIES`, their tests)
- Modify: `src/bin/o3_candidates.rs` (pass B writes `unmapped.fa` under `<out>.tmp/`, runs `minimap2(MM2_ATTRIB, copies_fa, unmapped.fa)`, attributes via `attribute_by_hits`; `families.tsv` gains nothing; the log line reports aligned / attributed counts; `--copies-fa` headers map target name -> family via `{fid}|...` prefix)
- Modify: `docs/MODULE_STATUS.md` only if the registry test demands it

- [ ] **Step 1: Failing tests** — `attribute_by_hits`: (a) best hit with `shorter_cov` 0.9 and de 0.05 -> attributed; (b) best hit coverage 0.4 -> None; (c) de 0.2 -> None; (d) two hits, the one with more matches decides even if the other has lower de.
- [ ] **Step 2: Run** (RED). **Step 3: Implement** in the library (pure) and wire pass B in the binary (stream unmapped records >= 300 bp to the FASTA while sweeping; one minimap2 call after the sweep; attributed reads join their family's net BEFORE the cap). Remove the k-mer index and its constants/tests. **Step 4: Run** `--lib o3_candidates`, `--test o3_candidates`, then the full suite (captured). Re-run the binary on the fixture: outputs byte-identical to the pre-change run (no unmapped reads there). **Step 5: Commit** — `o3_candidates: unmapped reads attributed by alignment to the copies (A13); k-mer index retired`.

### Task 2: Structural template and the consensus details

**Files:**
- Modify: `src/rustle/vg_family/o3_candidates.rs` (`MM2_MEMBERS` -> splice preset; `pub fn structural_template(members: &[usize], names: &[String], ava: &[PafHit], lens: &[usize]) -> usize` scoring big-indel bases from each member's `cs` against the other members; the insertion vote order in `consensus_from_template`; update the constant-pinning test)
- Modify: `src/bin/o3_candidates.rs` (`Net::longest` replaced by the structural template wherever a template is chosen: initial clusters, refine re-template, merged clusters; empty-merge fallback keeps absorbed clusters)

- [ ] **Step 1: Failing tests** — `retained_intron_read_is_not_the_template` (3 members: A full clean 900 bp, B = A + 300 bp intron inside, C = A with 0.2% errors; ava cs strings written by hand: B's hits carry a 300-bp insertion/deletion; template = A); `skipping_read_is_not_the_template` (one member lacks a 200-bp exon); `tie_breaks_by_length_then_name`; insertion vote: a 24-bp insertion with 3 carriers at a position where 4 of 8 carry a 1-bp insertion -> the 24-bp one is inserted (and the 1-bp one is not, since both cannot precede the same column — document the choice); refine re-template test (template split off -> new template chosen among the kept set); empty-merge fallback test.
- [ ] **Step 2: Run** (RED). **Step 3: Implement.** **Step 4: Run** focused + full suite; binary on the fixture: still one flagged candidate, union 850-950 bp (its length may change by a few bp — report it). **Step 5: Commit** — `o3_candidates: structurally central template, splice preset for votes, insertion vote by size class, refine re-template, empty-merge fallback (A13)`.

### Task 3: Acceptance A13 and the default flip

**Files:**
- Modify: `bench/rna_allele/accept_o3_candidates.sh` (work dir and prefixes parameterised: `A13` into `/mnt/linuxdisk/tmp/rna_allele/a13/`; reuse `A12.copies.*`), `docs/O3_CANDIDATES_ACCEPTANCE_A13_2026-10-03.md` (new), `docs/NEGATIVE_RESULTS_REGISTER.md` (rows 1221+), and — iff A13-1/2/3 all hold — `tools/rustle_pipeline.sh` (`CANDIDATES=1` default, `all` runs the stage, header/usage/README/AGENTS/REPRODUCE/figures docs updated, `run_e2e.sh` adapted) with the spec's §9b note.

- [ ] **Step 1:** rebuild the binary; run the stage in batches (as A12: 5 groups via `--families`, each under 10 min, `/usr/bin/time -v`); concatenate.
- [ ] **Step 2:** arm M exactly as A12 (rename contigs `iso_*`, index, realign the three parts, label contigs from the unmasked genome, `merge_test.py score`); A13-2 with the arm-M preset; A13-3 = summed batch time.
- [ ] **Step 3:** compute C = IsoCon's right D reads (Amendment 8's per-read calls, `merge_test.py` semantics on `linktest/RIL.bam`) over the truth-free attainable D reads (any record on a survivor in `R.bam`, or attributable by the Task 1 rule); A13-1 = stage D right >= 0.80 x C and false moves <= 5% (A12-1's 10,230 reported beside); cause table; attribution counts (unmapped / poorly placed: aligned, attributed, right family by `labels.tsv`); delta/2 and 2 x delta reruns.
- [ ] **Step 4:** verdicts as registered; write the doc and register rows; IF all three hold, flip the default (and say so in the doc); commit.


### Task 4: The no-deletion control for the stage (Amendment 14), added by ruling R22

**Files:**
- Modify: `bench/rna_allele/accept_o3_candidates.sh` / `.py` (a `CTRL` mode: unmasked `_pri`, `control/R0.bam`, a copies table of all 201 copies, the stage batched as A13, candidates classified against mat/pat with Amendment 9's rule — reuse `bench/rna_allele/control_test.py classify` logic — arm C = `_pri` + flagged unions)
- Create: `docs/O3_CANDIDATES_CONTROL_A14_2026-10-03.md`; modify `docs/NEGATIVE_RESULTS_REGISTER.md` (rows 1226+); iff C1'/C2' FAIL: revert `1f49d0f0` (the flip) with a dated note.

- [ ] **Step 1:** `panel_to_copies.py --all` -> `A14.copies.{tsv,fa}` (201 copies: `mask` + `keep` of `linktest/panel.json`, sequences from the UNMASKED `_pri`); the stage in 5 batches on `control/R0.bam` with `GGO.splice.mmi`; wall time recorded.
- [ ] **Step 2:** classify every flagged candidate (unions aligned to mat/pat, `splice:hq -uf -c -N 20`; lift of the copy intervals = `control/copies_lift.tsv`) into a / b / c as Amendment 9; C1' = families with >= 1 class-b/c flag <= 8 (15.7%); also the a+b+c rate and A9's 16/53.
- [ ] **Step 3:** arm C: `_pri` + flagged unions (renamed `iso_*`), realign the three control read parts (`linktest/scored.part*.fa`, same reads) with the pipeline flags, `control_test.py score`-style classification: C2' = false moves <= 5%.
- [ ] **Step 4:** doc, register rows, attribution counts and counters; decision per Amendment 14 (keep or revert the flip).

## Self-review notes
Spec coverage: Amendment 13's four bullets map to Tasks 1 (net), 2 (template + details), 3 (acceptance + flip). Review Focus 1-5 each name their test or check. No placeholders.
