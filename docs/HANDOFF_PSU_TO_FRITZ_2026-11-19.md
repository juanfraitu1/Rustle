# Handoff: what was done on PSU Roar (2026-11-18 / 19) and how a Fritz session takes it over

For the Claude session on Fritz. Everything below was done on PSU Roar (Penn State) while Fritz was in maintenance, in a mirror of
the Fritz setup at `/storage/group/kdm16/default/jxi21/fritz` (`psu:makova_group/fritz`). The code and documents are in git; the run
folders, the session notes (memory) and the ThesisVault edits are NOT, and steps 3-5 bring them over.

## 0. Already done from PSU (2026-11-19, evening, once Fritz was back)
- `~/Rustle` on Fritz fast-forwarded to origin/main (it was at 59c9b80d, tree clean). The rebuild was NOT done: `cargo build` on fritz1
  hung in a Lustre sync wait (`cl_sync_io_wait`) on its lock files in `~/rustle_target/release/` (/anvme, right after the maintenance);
  process 582694 on fritz1 may still be stuck there. When /anvme answers, rebuild (section 1) and check the SIRV md5; the binaries in
  `~/rustle_target` are the pre-PSU ones until then.
- Run folders of section 4 copied (light copy) into `~/ws/`, plus `r1828_support_profile` and `r1829_near_tie_census`.
- The PSU session notes and the five ThesisVault files of sections 5 and 6 are STAGED in `~/psu_handoff_2026-11-19/{memory,ThesisVault}`:
  merge them as sections 5 and 6 say (read from there instead of rsync from psu); nothing in the Fritz memory or vault was touched.
Still to do on Fritz: the SIRV md5 check (section 1), the merges, then the support-filter task (section 8).

**Done on Fritz, 2026-11-20:**
- Rebuild: the in-place build could not run (the hung `cargo` 582694 keeps the lock and cannot be killed). All binaries were built
  from `4e9bb9be` into `~/ws/rustle_target_b` and copied into `~/rustle_target/release/`; the earlier binaries are in
  `~/ws/rustle_target_bins_before_2026-11-20/`. The `.fingerprint` files there are older than the binaries: the next in-place
  `cargo build` recompiles once.
- SIRV md5 check: DONE 2026-11-20 11:25, after the storage target came back: `final.gtf` md5 `5cb0e7e4e7eb4cccc8f9def75cce2f65`
  with the binaries of `3a1005d1`, both from `~/ws/rustle_target_b` and from `~/rustle_target` (rebuilt in place the same hour,
  so it now carries `--polish-junction-ratio` and `RUSTLE_DUP_COUNT`); `~/ws/handoff_check_2026-11-20/run{2,3}`. Until then it
  could not run. Cause of the hang: Lustre target `nvme-OST000b` of /anvme is unavailable (`lfs check osts`); every
  file with data on it blocks its reader for good (`lfs find --ost 11 DIR` lists them). `r1801_sirv_wrong/sirv.bam` and several
  libraries of the Python environment (`libmpfr`, `zlib`, `_json`, parts of numpy / scipy) are among them.
- Session notes and the five vault files merged (earlier Fritz copies: `~/psu_handoff_2026-11-19/fritz_{memory,vault}_before_merge`).
- Section 8: row 1830 pre-registered, run (around the storage target: fresh tool environment, zebrafish re-mapped) and judged:
  FAIL, `docs/PREREG_support_filter_2026-11-20.md`, Result section.

## 1. Pull and rebuild (code)
```bash
cd ~/Rustle && git fetch && git log --oneline HEAD..origin/main   # every commit since fe9e7c23 (rows 1820-1829, the audit, this guide)
git pull --rebase
cargo build --release --bins            # with the pixi env on PATH, as usual; needed: src/family.rs and src/copy_assign_bin.rs changed
```
Rust change (commit e93f8cf3): opt-in `RUSTLE_SEED_SPLICED=support` (a spliced chain seeds only when >= 1 primary record carries it;
secondaries still add support; refused off the streaming `--assemble-only` path; driver comment in `tools/rustle_pipeline.sh`; unit
test `spliced_support_mode_seeds_only_primary_chains_and_counts_secondary_support`). Default output byte-identical (checked on PSU:
SIRV final.gtf md5 5cb0e7e4e7eb4cccc8f9def75cce2f65 and a whole HG002 quarter at `--seed-pool good`). Check the SIRV md5 again on
Fritz after the rebuild. No other default changed.
Later on 2026-11-19: the audit fixes of `docs/AUDIT_ROUNDS_2026-11.md` ("Audit of 2026-11-19") also need the rebuild. They do
not change the SIRV md5. (a) Finding 1: an `all` run started from `2bad9290` on (2026-11-19, before 04:00 on Fritz) seeded
with PRIMARIES although `PREFIX.assemble.env` says good; check with `grep -l 'seeding pool: primary' ~/ws/*/*.assemble.log`
against runs made by `all`. (b) The resume records no longer hold the thread / read-budget knobs, so an unfinished `all` or
`merged` redoes its stages once. (c) `copy_assign` refuses a FASTA without `.fai`; a malformed `RUSTLE_*` number stops the run.

## 2. What the rows say (register 1820-1828; read the documents, not this summary, before quoting)
| row | question | answer | document |
| 1820 | do reads that cannot be placed depict the copies they land on (pseudogene idea)? | no: C1 FAILS (1 copy gains its annotated structure) | `PREREG_undecidable_copies_2026-11-19.md` |
| 1821 | resolution-aware output: copy sets where reads cannot separate copies | R1 fails narrowly (+1% over the famwise layer), R2 passes (+2.5% over per-copy output, precision 0.888 -> 0.909), R3 passes; `RUSTLE_SEED_SPLICED=support` cuts wrong-copy chains 4x | `PREREG_resolution_call_2026-11-19.md`, `GOAL_multicopy_families_2026-11-19.md` |
| 1822 | borrowed structure; isoform-k=2; shared isoforms (KB3781) | K1, B1 FAIL (KB3781 leaves little room); shared isoforms descriptive | `PREREG_borrow_isoformk_2026-11-19.md` |
| 1823 | correlation clustering (Strainy port, copy_phasing) against HG002's own-haplotype truth | no gain for assignment; for copies the reference lacks the flag hits 13 / 13 families, false flags rise with depth | `CORRELATION_CLUSTERING_2026-11-19.md` |
| 1824 | oracle ceiling of correcting wrongly chosen primaries | +0.2% resolved units, no family member: not worth building | `PRIMARY_CORRECTION_CEILING_2026-11-19.md` |
| 1825 | are places reached only by secondaries real missed loci or shadows? | HG002 51 places, 6 real (pseudogenes, each dominated by another copy's reads); KB3781 0 | `SHADOW_LOCI_2026-11-19.md` |
| 1826 | tandem arrays: StringTie's behaviour, artifact reads, member size | StringTie joins units through reads > 1 period (50 -> 4 when removed); read classifier right on simulations; member-size rule exact on simulations but FAILED both held-out tests on OR6737 | `TANDEM_ARRAYS_2026-11-19.md`, `PREDICTION_tandem_heldout_2026-11-19.md` |
| 1828 | what do whole reads really support (IRFinder-S question)? | exact junction support; false SIRV transcripts are 1-3 bp shifts flagged by SHIFTED_DOMINANT / THIN_EXACT; no CNN needed | `SUPPORT_PROFILE_2026-11-19.md` |
| 1827 | are SDs where recent multi-copy families live? | yes on CHM13: 97% of human-specific paralogues in SDs (9% baseline), 99.5% of their pairs are one SD's two copies, 99% DNA-based | `SD_FAMILIES_2026-11-19.md` |
| 1829 | how many near-tied multimappers exist, and what could separate them? | near ties 0.3-3.5% of molecules, in SDs; HG002: 68% cannot be judged by its assembly, 88% of the rest undecidable (70% of those a read-length limit), reference PSVs place 82% of the decidable | `PREREG_near_tie_census_2026-11-19.md`; page https://claude.ai/artifact/JsK6mub28ZLP9sLDo7D8e8 |
Overview of the whole session: `MULTICOPY_FAMILIES_2026-11-19.md` (+ development ledger `MULTICOPY_DEV_LEDGER_2026-11-19.md`).
Also not registered: the user's family-member completeness test (StringTie as reference; missing members are mostly silent, none lost
to wrong placement; Rustle misses 24-75 expressed members StringTie finds), in the overview, step 7.
Answers for the advisor: `ADVISOR_QUESTIONS.md` (dated 2026-11-19 entries under Q3 and the Strainy item).
Open items: `PARKED.md` items 3 (multi-copy families), 4 (correlation clustering for missing copies), 5 (tandem member size).

## 3. Data spent (important for any future row)
- HG002 Kinnex on CHM13, all four hash quarters are now SPENT: 0-24 seen (rows 1776-1791), 25-49 row 1796, 50-74 row 1820,
  75-99 row 1821. No unseen HG002 quarter remains; reuse is allowed (user ruling) but say so.
- KB3781 judged in row 1822. OR6737 and A119b are spent for the tandem member-size rule (do not re-tune it on them).

## 4. Run folders (rsync from Fritz)
On Fritz first `ssh -fN psu` and approve the MFA push (see the PSU memory note `psu-roar.md`). Sizes on PSU:
| folder | full | without BAMs, indices, FASTA/FASTQ, binaries |
| `r1820_undecidable_copies` | 30 GB (23 GB are rebuilt HG002 haplotype `.mmi`) | 450 MB |
| `r1821_resolution_call` | 1.8 GB | 126 MB |
| `r1822_borrow_isoformk` | 85 MB | 42 MB |
| `r1829_near_tie_census` (layer2/sel.bam 170 MB is rebuildable) | 722 MB | about 400 MB (truth/placements.fa is 330 MB) |
| `goal_multicopy_dev` (rows 1823-1827 link into it; dev ledger, simulations, member checks) | 33 GB | 4.3 GB |
| `r1823_correlation_clustering` .. `r1827_sd_families` | symlinks into `goal_multicopy_dev` | |
Light copy (recommended; the BAMs can be rebuilt from the job scripts in each folder):
```bash
for d in r1820_undecidable_copies r1821_resolution_call r1822_borrow_isoformk goal_multicopy_dev \
         r1823_correlation_clustering r1824_primary_correction r1825_shadow_loci r1826_tandem_arrays r1827_sd_families r1829_near_tie_census; do
  rsync -a --exclude='*.mmi' --exclude='*.bam' --exclude='*.bai' --exclude='target/' --exclude='snap/target/' \
        --exclude='*.fa' --exclude='*.fq' psu:makova_group/fritz/ws/$d ~/ws/
done
```
Job scripts inside these folders use PSU paths (`/storage/group/kdm16/default/jxi21/fritz/...`, `--account=open`); on Fritz replace
the root with `/anvme/workspace/v138da13-rustle` and the Slurm header with the Fritz one.

## 5. Session notes (memory) changed on PSU: merge into the Fritz memory directory
The PSU copy of the memory lives at `psu:makova_group/fritz/fritz_memory/` (copied from Fritz on 2026-11-19, then edited). Files
changed on PSU: `MEMORY.md` (index lines for family-wise-layer and PSU Roar), `family-wise-layer.md`, `o1-families.md`,
`o2-assignment.md` (rows 1824, 1829), `working-rules.md` (no new locks; the 2026-11-19 driver audit), `o3-missing-copies.md`, `assembly.md`, and the NEW `psu-roar.md` (written by the Fritz session itself and copied
over). Fritz sessions may have edited the same files since; merge by appending the dated PSU paragraphs, not by overwriting:
```bash
rsync -a psu:makova_group/fritz/fritz_memory/ /tmp/psu_memory/
for f in MEMORY.md family-wise-layer.md o1-families.md o2-assignment.md o3-missing-copies.md assembly.md working-rules.md; do
  diff <(cat ~/.claude/projects/*/memory/$f) /tmp/psu_memory/$f | grep '^>' | head    # the PSU additions to merge
done
```

## 6. ThesisVault edits on PSU (not in git)
`~/ThesisVault` on PSU: `Papers/Kazantseva_2024_Strainy.md` ("Own results" section), `Topics/O2.md` and `Topics/O3.md` (one dated line
each under the Strainy entry), `Topics/SegDups.md` ("Own result (added 2026-11-19)"), `CHANGELOG.md` (two 2026-11-19 entries). Copy:
`rsync -a psu:makova_group/fritz/ThesisVault/{Papers/Kazantseva_2024_Strainy.md,Topics/O2.md,Topics/O3.md,Topics/SegDups.md,CHANGELOG.md} <vault>/...`
(check the Fritz vault did not change those files first).

## 7. Rules that applied (and still apply)
Register rows 1820-1827 were claimed on PSU; check the register tail on origin/main before claiming 1829. Commits on PSU included
only their own paths. Pre-registrations were committed before judging data were opened (1820, 1821, 1822; predictions for the
row-1826 held-out tests). Development results are labelled as such in every document.

## 8. Next task for the Fritz session (user, 2026-11-19): does the support profile, used as a FILTER, remove wrong transcripts?
Row 1828 (`docs/SUPPORT_PROFILE_2026-11-19.md`) only DETECTS. On SIRV `final.gtf` its flags (SHIFTED_DOMINANT / THIN_EXACT, relative)
caught 23 / 28 false transcripts and 2 / 48 true ones, but SIRV and the HG002 seen quarter were used to build the rule, so those
numbers are development results and say nothing about a filter. The user asked to test the filter on Fritz:
1. Claim the next free register row on origin/main (1829 is taken on PSU by the near-tie census; check the tail) and write `docs/PREREG_support_filter_<date>.md`, committed
   BEFORE any judging data is opened. Fix in it:
   - the filter: drop multi-exon transcripts of `PREFIX.final.gtf` flagged SHIFTED_DOMINANT or THIN_EXACT, rule exactly as
     committed in `bench/famwise_eval/support_profile.py` (no re-tuning; NEAR_MISS_DOMINANT, WEAK_JUNCTION, UNCOVERED, NO_FULL_READ
     are NOT filtered unless the pre-registration says so and why);
   - the judging data: NOT SIRV, NOT the HG002 0-24 quarter. Candidates on Fritz: another spike-in library with known truth
     (Sequins / SIRV-Set 4 / Lexogen, if one exists there), a fresh simulation not used before, or an HG002 quarter reused WITH that
     stated (all four are spent; reuse allowed by the user ruling of section 3);
   - truth labels by exact intron chain (`label_by_chain.py`), and the bars: false transcripts removed (share), true transcripts
     lost (share), precision and sensitivity before / after; pass = precision up with sensitivity loss below a stated bound
     (suggest: >= 50% of false removed, <= 2% of true lost; the user decides);
   - a GffCompare / SQANTI3 view as secondary, StringTie on the same BAM as the reference point.
2. Run, judge against the bars, write the result document, register the row, commit and push.
Default products must not change: the filter is a post-step on a copy of `final.gtf` until a row decides otherwise.
