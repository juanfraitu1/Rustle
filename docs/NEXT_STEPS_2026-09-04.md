# Next steps after 2026-09-04 (ledger §6ds–§6ee)

Status: **O1 is two measurement sessions from a conclusion with stated limits. O2 has one honest claim on
NPIP (an abstention certificate) and needs a node fix and a new truth before more can be claimed.**
Canonical catalogs: `mcl_ann/rna_bp1_p9.clusters.tsv` (3 contigs) and `mcl_ann/gw_bp1_p9.clusters.tsv`
(genome-wide), both at `--prune 1e-9` (§6ed). Every number quoted from a catalog must say which prune.

## O1 — to conclude
| # | step | why | where it stands |
|---|---|---|---|
| O1-1 | **Non-overlapping-locus NODE rule** in `mcl_families` (overlapping annotation records on a contig = one locus; representative = longest exon-union) | nested annotations are the same DNA twice: 13 NPIP copy pairs, 607/1,221 O2 ties (§6ee); the nested MCL4 lncRNA (§6dy); NPIP giant models are real transcription units (§6ea) | **SHIPPED OFF, MEASURED (§6ef)**: exon-union-overlap loci, edges attributed to the representative; 3-contig 142→116 clusters, NPIP 44 records → 33 loci; ⚠ **NPIP reunites with MCL7 (U1) into 48 loci — the anchored boundary moves; the default flip is the user's call.** O2 on it: resolved 17/43 vs 9/39, Containment 10→2 (the 2 share no exon), MAPQ-0 still 1/78 |
| O1-2 | Run the audit columns on the genome-wide catalog `gw_bp1_p9`: curated-repeat fraction, member union-exonic-coverage, size flags | the 46-cluster audit (§6eb, `docs/audit46_2026-09-04.tsv`) was on the 3-contig 1e-5 slice; the reported object is genome-wide | scripts exist: `mcl_ann/adj/audit46/audit46.py`, `adj/gwaudit/gwaudit_p9.py` |
| O1-3 | **Cut certificate**, pre-registered: "members share ≥ half their exonic ≥0.85 partners outside the cluster" (§6ed) — positive control = the tandem 60/58 split (must flag each other), negative = NPIP | the 13 cell-A cuts are uncertified (§6dy); this is the first instrument sharing no input with MCL's edge density | statistic implemented (`gwaudit_p9.py`), thresholds NOT pre-registered as a cut test yet |
| O1-4 | Library-silent element handle: a de novo-family column (`rnd-*` family share + young-mode contig spread) | MCL4's element is absent from curated Dfam and "Unknown" de novo (§6dz); the curated column misses it | scripts in `adj/rmsk/`; needs the full-genome de novo .out parsed per cluster |
| O1-5 | Two-sided exonic clause at COPY level (γ) | priced on edges only (§6dy cell B) | offline candidate; run `mcl_families` with a both-sides floor and diff clusters |
| O1-6 | Write the limitations: node boundaries are copy-inconsistent in the annotation and RNA cannot fix them (§6ea); cohesion is slice-conditioned; families >~1,600 near-identical copies exceed even prune 1e-9 at I=2.8 | | text only |
| O1-7 | Restate §6dd/§6de (inflation plateau, "zero clusters >100") under prune 1e-9 in `docs/ADVISOR_QUESTIONS.md` | they were prune artefacts (§6ec) | text only |

## O2 — to claim anything beyond NPIP's abstention certificate
| # | step | why | where it stands |
|---|---|---|---|
| O2-1 | Same as O1-1: O2 must receive non-overlapping loci | 9/39 → 16/27 resolved copies when nested copies are dropped (§6ee) | DONE via O1-1 (`o2scale/fam_NPIPloci_073242`); ⚠ O2 peaked at **20.7 GB of 25** on 43 loci — memory must be looked at before any sweep |
| O2-2 | **A truth for the correction leg's re-threaded reads** | `--vg-realign-correct`'s +31 assignments are reads whose primary lies OUTSIDE every copy; unique-mapper agreement cannot score them (§6ee) | design: simulate reads from copy sequences (the existing sim harness, `-N 50`), or use copy-specific junctions as truth |
| O2-3 | Report O2 per copy WITH O1's nearest-paralogue identity as the abstention forecast (ρ −0.55; ≥0.99 ⟹ share 0) | the certificate that ties O1 to O2 | column to add to the joint catalog |
| O2-4 | Decision with the advisor: is "abstention certificate + identity forecast" the O2 result? | on NPIP in fibroblasts the multimapping pool is 592/6,545 reads and 98.6% of it is unassignable | user/advisor decision |
| O2-5 | Optional: the 86-family sweep (`o2scale/run_armA.sh`, resumable) EXCLUDING cell-B artefacts; ⚠ the 39-copy family peaked at 18.8 GB of 25 | only if O2-4 says O2 is reported catalog-wide | ready |
| O2-6 | Dispersed families: scope the claim ("assignment for local families, declared abstention for dispersed") or build a genome-wide mode | the genome-wide catalog is mostly dispersed families | user decision |

## Housekeeping
- The whole MCL line is committed (0634a0e, dcd41d7); everything after §6ed is uncommitted.
- Repeat library on disk: `winloci_data/rmsk/` (UCSC GenArk; curated Dfam 3.8 + RepeatModeler).
- Never run O2 on the cell-B artefacts (MCL28, MCL15, NEW101/MCL0): they are L1/ERV1/L1, not families.
