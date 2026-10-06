# Second machine (small disk): reproduce and extend the Soto evidence, then merge (2026-09-30)

For a second WSL machine with little disk, running its own Claude. Goal there: reproduce and extend two claims from
`docs/SOTO_VS_OURS_MEETING_2026-09-30.md`: (A) Soto's families sit inside ours and, with the right conditions, we find what they
find; (B) Soto's families are narrower than sequence homology and include fragments and pseudogenes. Nothing there needs
a build, a BAM, a genome or the 35 GB of copy-number tracks.

## 1. What to copy

| what | size | how |
|---|---|---|
| the repository | ~60 MB working tree (+ ~250 MB history) | `git clone --depth 1 git@github.com:juanfraitu1/Rustle.git` (shallow is enough: a branch pushed from a shallow clone merges normally here) |
| `famcn_ours_allwssd.tsv` | 0.2 MB | from `/mnt/linuxdisk/home/juanfraitu/winloci_data/soto_replication/`; any folder, e.g. `~/soto/` |
| `famcn_ours_all.tsv` (the 10-sample table) | 0.2 MB | same folder; only for the 10-sample ladder row |
| optional, for the second Claude's context | ~75 KB | from `~/.claude/projects/-mnt-c-Users-jfris-Desktop/memory/`: `project_soto_refinement_evidence.md`, `project_soto_full_replication.md` (45 KB), `reference_soto_2025_hsd_brain.md`, `reference_advisor_canzar.md`, `feedback_metric_traps.md`, `project_soto_family_pseudogene_fragment_audit.md`, `project_cover_and_jn_refuted.md` |
| optional, task B2 only | 0.8 MB | `sd98_gene_exons.tsv` from the same folder (exon footprints per member) |

Not needed, and what they would cost: the A119b BAM (96 GB), the WSSD tracks (35 GB), `t2t-chm13-v1.0.fa.gz` (0.9 GB),
`final_v1_clean.bed` (65 MB), `cat_v4.bed` (81 MB), any PAF. The repo already holds the frozen inputs: `bench/soto/soto_famCN_S1C.tsv`
(Soto's table), `shared_exons_5154_exon_mapback.tsv` (the exon edges), `soto_split_2026-09-29.tsv` (the dev / held-out split).
Until the newest commit is pushed, copy `bench/soto/soto_replication.py` and these docs by hand (they are small).

Environment: python3 ≥ 3.9 (numpy; scikit-learn optional, a fallback ARI is built in). No Rust.

## 2. Reproduce (must match exactly before anything else)

```
python3 bench/soto/soto_replication.py genesets --out-eligible elig.tsv --out-full full.tsv
python3 bench/soto/soto_replication.py nesting --shared bench/soto/shared_exons_5154_exon_mapback.tsv \
    --geneset elig.tsv --full-geneset full.tsv --famcn-ours ~/soto/famcn_ours_allwssd.tsv
python3 bench/soto/soto_replication.py ladder  --shared bench/soto/shared_exons_5154_exon_mapback.tsv \
    --geneset elig.tsv --full-geneset full.tsv --famcn-ours ~/soto/famcn_ours_allwssd.tsv \
    --famcn-ours10 ~/soto/famcn_ours_all.tsv --split bench/soto/soto_split_2026-09-29.tsv --drop-family ID_356
```

Expected `nesting`: 444 Soto families with ≥ 2 clean genes; sequence only 440 / 444 (99.1%) inside one of our clusters and 394 / 398
(99.0%) of our clusters are unions of whole Soto families; 33 of our clusters hold 87 Soto families; exceptions ID_192, ID_347,
ID_482, ID_62. Expected `ladder` (ALL ARI / exact): sequence only 0.7307 / 345; our famCN 10 samples 0.9198 / 373; 268 samples
(exons) 0.8855 / 375; 268 samples, Soto's interval 0.9277 / 411; S1C 0.9698 / 479 (held-out 0.9681 / 263).

## 3. Suggested tasks there (each a small script in a NEW folder, plus a short doc)

- **B1, Soto's table audit (S1C only):** per family, biotype composition from the `Biotype` column. Checks to match: 2,572 members /
  605 families, 56.5% of members pseudogene-biotype, 217 families (35.9%) entirely pseudogene, 287 (47.4%) with no protein-coding member.
- **B2, fragments bundled with full-length members (needs `sd98_gene_exons.tsv`):** family members' exonic footprints versus the
  family's largest member; check: 147 of 420 size-comparable families hold a member < 20% next to one ≥ 80%.
- **B3, Soto finer than ours:** for each of our 33 multi-family sequence clusters, list its Soto families with sizes and the genes that
  separate them (the copy-number values), and the 83 missing paralog pairs of `docs/SOTO_REPLICATION_STATUS_2026-09-28.md` §3.
- **A2 (optional):** nesting for Soto's published edges (`shared_exons_2334_finalv1_native.tsv`) versus the exon edges: does the 99.1% depend on the edge set?

Report numbers with the source file and the command. Do not call our families "better": the objectives differ and Soto's
table is the truth (see §4 of the meeting sheet).

## 4. Merge protocol (avoids conflicts with this machine)

1. There: `git checkout -b machine2/soto-evidence`. Add NEW files only: `bench/soto_m2/*.py`, `docs/SOTO_M2_*.md`.
2. Do NOT edit: `src/`, `tools/`, `bench/soto/soto_replication.py`, `README.md`, `REPRODUCE.md`, `docs/MODULE_STATUS.md`,
   `docs/NEGATIVE_RESULTS_REGISTER.md`, any existing `docs/PREREG_*`. Draft register rows in `docs/REGISTER_DRAFTS_machine2.md`
   (renumbered here on merge; this machine owns the register numbers).
3. Commit with the usual attribution lines; `git push -u origin machine2/soto-evidence`. Never push to `main`, never force-push.
   GitHub login there: `gh auth login` (device flow) or an HTTPS token; do not copy private keys between machines.
4. Here: `git fetch origin && git merge --no-ff origin/machine2/soto-evidence`, re-run the two commands above, renumber the register drafts, commit.
