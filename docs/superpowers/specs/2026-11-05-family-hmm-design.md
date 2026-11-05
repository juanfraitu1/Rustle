# Family HMM layer (opt-in): design

Date 2026-11-05. Status: approved in chat, spec for the implementation plan.

## Purpose

Give every multi-copy family a standard profile-HMM representation (the way Pfam defines a family), as an
opt-in bench layer that touches no shipped output. Two uses: (1) portable model files reviewers and other
labs can search with; (2) a standing concordance table that says whether the profiles agree with the MCL
families. It is NOT a member-discovery or read-assignment emitter: that arm ran as unmapped-rescue
Amendment 5 (2026-10-08) and gained 243 reads at double the wrong-join rate, carried by two families.

Honest limit, to be stated wherever the result is shown: Rustle copies sit at 90-99% nucleotide identity,
so a DNA profile is close to a weighted consensus and adds no sensitivity over minimap2. The value is the
standard model file and the Pfam link. The variation graph remains the thesis model.

## Inputs

A run PREFIX (e.g. `~/runs/human_A119b/a1`): `PREFIX.fam.copies.tsv` (family_id, copy_idx, tid, ...),
`PREFIX.fam.copies.fa` (headers `FAM|idx|chrom:start-end|strand|nexon=N`), `PREFIX.fam.loci.fa` (all de
novo loci, headers `chrom:start-end`). Tools from the pixi `bench` env: mafft, hmmbuild, hmmpress, nhmmer,
hmmscan, hmmemit, esl-translate. Pfam-A.hmm: one-time download into the workspace
(`$WS/pfam/Pfam-A.hmm`, pressed); absent = the script prints the download command and exits unless
`--no-pfam`.

## Script

`bench/famhmm/famhmm.py build|scan|concord|all --prefix PREFIX --out DIR [--only-families FILE]
[--threads N] [--no-pfam] [--pfam PATH] [--min-copies 2] [--min-orf-aa 50]`

- **build.** Per family with >= min-copies copies: MAFFT (`--auto`) the copies -> `DIR/dna/FAM.hmm`
  (`hmmbuild --dna -n FAM`). Longest ORF per copy (`esl-translate -l min-orf-aa`, longest record) ->
  MAFFT -> `DIR/prot/FAM.hmm` (`hmmbuild --amino`). Copies without an ORF are left out of the protein
  profile; a family with < 2 ORF copies gets `prot_status=none`. `DIR/orf.tsv`: family, copy tid, frame,
  nt offset, aa length: this is the recorded DNA-protein equivalence. Concatenate + hmmpress into
  `DIR/families.dna.hmm` and `DIR/families.prot.hmm`.
- **scan.** `hmmemit -c` the protein profiles, `hmmscan --domtblout` against Pfam-A ->
  `DIR/pfam_labels.tsv`: family, Pfam domains in query order (E <= 1e-5), best E, prot_status.
- **concord.** `nhmmer --tblout --cpu N DIR/families.dna.hmm PREFIX.fam.loci.fa`, parsed with
  `bench/unmapped_rescue/profiles.parse_tblout`. A hit counts if E <= 1e-10 and the hit span covers >= 50%
  of the profile length. `DIR/concordance.tsv`: family, n_copies, n_recovered (own copies hit by own
  profile), n_extra (non-member loci hit), n_top_own (member loci whose best-scoring profile is their
  own family). Footer line with totals.

## Pre-registered decision rule (fixed before any run)

The layer is CONCORDANT on a sample if, over families with >= 2 copies,
sum(n_recovered)/sum(n_copies) >= 0.95 AND sum(n_top_own)/sum(n_copies) >= 0.95.
Otherwise FAIL, reported as is; thresholds are not loosened post hoc. Extra loci are reported, not scored.
First run: the 57-family panel `bench/o2levers/panel.families.txt` on A119b with a 1 h wall; genome-wide
(1,565 families vs 76,161 loci) only after the panel completes and the per-family cost is known.

## Outputs and cost

Everything under `--out` in the workspace (`/anvme/workspace/v138da13-rustle/famhmm/<sample>`), never
home. No shipped file changes; `PREFIX.gtf`/families outputs stay byte-identical.

## Testing

`bench/famhmm/test_famhmm.py`: a toy 3-copy family with a 60-aa ORF (synthetic, in the test) -> build
produces both profiles, `orf.tsv` has frame/offset as planted, concord recovers 3/3 with n_top_own 3,
rule evaluates CONCORDANT. Skips if hmmbuild is not on PATH.

## Out of scope

HMM member discovery / read assignment; HMM-vs-HMM comparison (HH-suite); a driver stage (one case block
later if wanted); Rust changes.

## Result (A119b panel, 2026-11-05)

Build: 57 panel families, all with a DNA and a protein profile, 7 min at 16 threads (Pfam scan skipped: Pfam-A not
downloaded). Concordance target: `a1.fam.copies.fa` (42 MB, spliced copies), not `fam.loci.fa` (1.86 Gb; job 4442327
timed out at the 1 h wall with no output). nhmmer cost: ~1,300 CPU-min for 50 families on 72 cores (one single-cpu
nhmmer per family; `--cpu 16` on one process did 10 families in 1 h); 6 more needed another 4 h; MCL1312 (14 kb profile)
is still unfinished after 6 CPU-h (job 4443368 runs it alone). Genome-wide (1,565 families) is therefore days: the
concordance layer stays a panel-scale check; the export layer is cheap.

**Rule: FAIL as registered** on the 56 finished families: recovered 0.7131, top_own 0.7089 (both < 0.95); 417 extra
copy hits. Cause (diagnostic, not a rescue): every one of the 237 copies is hit by its own profile at E <= 1e-10; the 68
misses are hits below the 50% cover bar, because a MAFFT profile spans the union of all copies' columns (profile /
median copy length up to 2.95x; MCL884 has a 21 kb single-exon copy beside 4.4 kb ones), median cover of the misses
0.29. Post hoc, cover measured against the copy's own length instead of the profile's gives recovered 0.890 / top_own 0.890: still FAIL,
not claimed. Table: `/anvme/workspace/v138da13-rustle/famhmm/A119b_panel/concordance.tsv`.
