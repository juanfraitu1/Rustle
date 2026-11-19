# Handoff: what was done on PSU Roar (2026-11-18 / 19) and how a Fritz session takes it over

For the Claude session on Fritz. Everything below was done on PSU Roar (Penn State) while Fritz was in maintenance, in a mirror of
the Fritz setup at `/storage/group/kdm16/default/jxi21/fritz` (`psu:makova_group/fritz`). The code and documents are in git; the run
folders, the session notes (memory) and the ThesisVault edits are NOT, and steps 3-5 bring them over.

## 1. Pull and rebuild (code)
```bash
cd ~/Rustle && git fetch && git log --oneline HEAD..origin/main   # 21 commits, fe9e7c23 .. e0a194bf (and this guide)
git pull --rebase
cargo build --release --bins            # with the pixi env on PATH, as usual; needed: src/family.rs and src/copy_assign_bin.rs changed
```
Rust change (commit e93f8cf3): opt-in `RUSTLE_SEED_SPLICED=support` (a spliced chain seeds only when >= 1 primary record carries it;
secondaries still add support; refused off the streaming `--assemble-only` path; driver comment in `tools/rustle_pipeline.sh`; unit
test `spliced_support_mode_seeds_only_primary_chains_and_counts_secondary_support`). Default output byte-identical (checked on PSU:
SIRV final.gtf md5 5cb0e7e4e7eb4cccc8f9def75cce2f65 and a whole HG002 quarter at `--seed-pool good`). Check the SIRV md5 again on
Fritz after the rebuild. No other default changed.
Later on 2026-11-19: the audit fixes of `docs/AUDIT_ROUNDS_2026-11.md` ("Audit of 2026-11-19") also need the rebuild. They do
not change the SIRV md5, but `all` runs before them seeded the assembly with PRIMARIES although `PREFIX.assemble.env` says good
(finding 1): an `all` result quoted as the good pool from before this commit was in fact the primary pool. `copy_assign` now
refuses a FASTA without `.fai`, and a malformed `RUSTLE_*` number stops the run.

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
| `goal_multicopy_dev` (rows 1823-1827 link into it; dev ledger, simulations, member checks) | 33 GB | 4.3 GB |
| `r1823_correlation_clustering` .. `r1827_sd_families` | symlinks into `goal_multicopy_dev` | |
Light copy (recommended; the BAMs can be rebuilt from the job scripts in each folder):
```bash
for d in r1820_undecidable_copies r1821_resolution_call r1822_borrow_isoformk goal_multicopy_dev \
         r1823_correlation_clustering r1824_primary_correction r1825_shadow_loci r1826_tandem_arrays r1827_sd_families; do
  rsync -a --exclude='*.mmi' --exclude='*.bam' --exclude='*.bai' --exclude='target/' --exclude='snap/target/' \
        --exclude='*.fa' --exclude='*.fq' psu:makova_group/fritz/ws/$d ~/ws/
done
```
Job scripts inside these folders use PSU paths (`/storage/group/kdm16/default/jxi21/fritz/...`, `--account=open`); on Fritz replace
the root with `/anvme/workspace/v138da13-rustle` and the Slurm header with the Fritz one.

## 5. Session notes (memory) changed on PSU: merge into the Fritz memory directory
The PSU copy of the memory lives at `psu:makova_group/fritz/fritz_memory/` (copied from Fritz on 2026-11-19, then edited). Files
changed on PSU: `MEMORY.md` (index lines for family-wise-layer and PSU Roar), `family-wise-layer.md`, `o1-families.md`,
`o2-assignment.md`, `o3-missing-copies.md`, `assembly.md`, and the NEW `psu-roar.md` (written by the Fritz session itself and copied
over). Fritz sessions may have edited the same files since; merge by appending the dated PSU paragraphs, not by overwriting:
```bash
rsync -a psu:makova_group/fritz/fritz_memory/ /tmp/psu_memory/
for f in MEMORY.md family-wise-layer.md o1-families.md o2-assignment.md o3-missing-copies.md assembly.md; do
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
