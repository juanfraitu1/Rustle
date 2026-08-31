# Superpowers plans & specs — archive index

These 30 design specs and implementation plans (2026-06-01 … 2026-07-26) were the *process* record for
work that has since shipped, been retired, or been superseded by the measured record in
`docs/o1_ledger.md` and `docs/NEGATIVE_RESULTS_REGISTER.md`. **The files were removed from the working
tree to cut documentation sprawl; every one remains in git history** — recover any with:

```
git log --all --diff-filter=D -- docs/superpowers/specs/<name>.md   # find the deleting commit
git show <commit>^:docs/superpowers/specs/<name>.md                 # print it
```

⚠ A spec records what was INTENDED at the time. Where a spec and the ledger disagree, **the ledger is
the measured record and wins** — several of these specs describe approaches later refuted (the EM
copy-assignment line, the DNA family repeat gate, the copy-graph objects).

| kind | file | size | subject |
|---|---|---|---|
| spec | `2026-07-17-copy-graph-objects-design.md` | 12K | Copy-Graph Objects (v1) — Design Spec |
| spec | `2026-07-17-copy-graph-v2-design.md` | 11K | Copy-Graph v2 — Exon Presence/Absence Arms + MI Tag — Design Spec |
| spec | `2026-07-18-augment-linearize-design.md` | 11K | Augment-and-Linearize Validation for Reference-Absent Copies — Design Spec |
| spec | `2026-07-18-noncoding-promotion-design.md` | 13K | Non-Coding-Aware Promotion of Reference-Absent Copies — Design Spec |
| spec | `2026-07-18-soto-dna-vg-ceiling-design.md` | 11K | DNA Variation-Graph Ceiling vs Soto SDs — Design Spec |
| spec | `2026-07-19-vg-o2-substrate-design.md` | 8K | The Family Variation Graph as the O2 Assignment Substrate — Design Spec |
| spec | `2026-07-20-mechanism-transparency-design.md` | 18K | Mechanism Transparency — Design (Deliverable A) |
| spec | `2026-07-21-code-consolidation-design.md` | 9K | Deliverable B — Byte-Identical Code Consolidation — Design |
| spec | `2026-07-21-identifiability-merge-design.md` | 8K | Identifiability-Respecting Locus Merge — Design |
| spec | `2026-07-21-nonexon-rescue-poc-design.md` | 4K | Non-Exon-Signal Rescue POC — Design |
| spec | `2026-07-21-unspliced-seeding-design.md` | 7K | Position-Aware Seeding of Unspliced Reads — Design |
| spec | `2026-07-22-tied-secondary-seeding-design.md` | 6K | Multi-mapper (tied-secondary) seeding — Phase 1 (benchmark-validated) |
| spec | `2026-07-22-tied-seed-existence-only-design.md` | 7K | Tied-seed copies as existence-only appends — design |
| spec | `2026-07-23-mischain-read-salvage-design.md` | 8K | Mis-chain Read Salvage — Design |
| spec | `2026-07-24-dna-family-repeat-gate-design.md` | 4K | DNA-family-fallback Repeat/Complexity Gate — Design |
| spec | `2026-07-25-mischain-G-cap-realign-followup.md` | 6K | Mis-chain rescue via targeted `-G`-capped re-alignment — follow-up |
| spec | `2026-07-26-dna-from-genome-mode-design.md` | 11K | `--from-genome` mode: reproduce Soto's families from the genome alone — Design |
| plan | `2026-07-17-copy-graph-objects.md` | 45K | Make every copy of a family a tagged, corroborable path in one Bandage-loadable family variation graph emitted by `copy_assign --phase`, with a REFERE |
| plan | `2026-07-17-copy-graph-v2.md` | 35K | Add the `MI` genome-map-identity tag and a sibling `<out>.exon.gfa` exon presence/absence graph to `copy_assign --phase`, so a copy differing by exon  |
| plan | `2026-07-18-augment-linearize.md` | 28K | Emit a threshold-free "linearize certificate" for each candidate reference-absent copy — augment the local family reference with the candidate, re-ali |
| plan | `2026-07-18-noncoding-promotion.md` | 32K | Recover credible non-coding / lncRNA reference-absent collapses that the protein-BLASTx funnel drops, by promoting them on a collapse-quality bar into |
| plan | `2026-07-18-soto-dna-vg-ceiling.md` | 28K | A new bench script `bench/soto/soto_dna_vg.py` that builds, per flagship Soto SD family, a base-level variation-graph GFA from the family's member gen |
| plan | `2026-07-19-vg-o2-substrate.md` | 24K | Make the per-family variation graph the explicit substrate the O2 copy-assignment decision runs on — an auditable ad-hoc reference (copies = paths) ea |
| plan | `2026-07-20-mechanism-transparency.md` | 41K | Make Rustle's single method legible end-to-end and disclose every heuristic, so the advisor sees two mechanisms + two consequences (not four approache |
| plan | `2026-07-21-code-consolidation.md` | 26K | Remove the duplicate implementations deliverable A disclosed (two refiners, two ≥2-loci predicates, three EM paths) and fix config hygiene (CLI flags, |
| plan | `2026-07-21-identifiability-merge.md` | 19K | Make the locus merge respect read-separability — two co-located same-strand copies merge only when no read distinguishes them (χ(H)) — so the 28 "dist |
| plan | `2026-07-21-unspliced-seeding.md` | 20K | Stop `pass1_skeletons_robust` from pooling all unspliced reads on a chromosome into one over-length skeleton (rejected by `MAX_SPLICED`), by clusterin |
| plan | `2026-07-22-tied-secondary-seeding.md` | 19K | Seed a candidate locus from AS-tied secondary reads that share an intron chain when no primary anchors it, to recover covered-but-tied K=0 Soto member |
| plan | `2026-07-23-mischain-read-salvage.md` | 20K | Split reads mis-chained across spurious giant introns into their local segments *before* skeleton seeding, so a real duplication locus is not lost whe |
| plan | `2026-07-26-dna-from-genome-mode.md` | 26K | Add a read-free, annotation-free `--from-genome` front-end to `gw_family_catalog` that discovers segmental-duplication families from the CHM13v2.0 gen |
