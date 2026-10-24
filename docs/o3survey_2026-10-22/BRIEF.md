# Shared brief for the O3 fingerprint survey (2026-10-22, 21:10 PDT; user has a meeting at about 23:10)

## Goal
Find ways to strengthen O3 using ONLY the sample's long RNA reads (PacBio HiFi / Kinnex) aligned to ONE reference (CHM13). The
user's rules:
- no DNA reads, no WGS, no per-sample assembly;
- donor assemblies are answer keys only.

**The advisor's view (Canzar):** "we should be able to cluster by PSV, BAM tags and find a certain print of a missing copy". It
does not need to be definitive (that is what DNA WGS is for), but it should suggest when to look harder.

So O3 is a SCREENING / PRIORITISATION flag (recall-oriented, ranked), not a hard decision. Canzar likes:
- clean combinatorial structure;
- provable statements;
- no arbitrary thresholds;
- no 1/k splitting;
- no minimizers.

## Where we are (read these)
- **O2-VC (the sample's own longcallR calls on confident reads, assign-or-abstain) results:**
  `/mnt/linuxdisk/home/juanfraitu/rustle_soto_subset/docs/archive/2026-10/merged/PREREG_o2_variant_calls_2026-10-22.md`, section "Results".
  - The held-out donor HG00733 fails at 6.3% wrong, mostly PRAMEF.
  - The errors concentrate in copy-number-variable families: PRAMEF, GTF2H2C, SPDYE1, NBPF1, ANKRD36C.
- **Exploratory gates:** `/mnt/linuxdisk/home/juanfraitu/rustle_soto_subset/docs/archive/2026-10/merged/O2O3_EXPLORE_2026-10-22.md` (§5 especially).
  - **Key insight:** a reference copy the donor LACKS (PRAMEF7.1 in HG00733) is taught by its paralog's reads and looks, from RNA
    alone, like a present copy with 4 private hom-alt SNVs.
  - The swap-pair gate is inert; the taught-by-paralog gate is harmful; a hindsight gate (H5) catches PRAMEF7.1.
- **Pre-registered expected failures**, with mechanisms and the literature: F1-F5:
  `/mnt/linuxdisk/home/juanfraitu/rustle_soto_subset/docs/archive/2026-10/merged/PREREG_o2_variant_calls_expected_failures_2026-10-22.md`
- **Prior O3 attempts.** Memory notes in `/home/juanfra/.claude/projects/-mnt-c-Users-jfris-Desktop/memory/`:
  - `project_o3_npip_maternal_copy_2026-10-21.md`: reference-free IsoCon chain; detection works, read recovery does not.
  - `project_samplecopy_2026-10-19.md`: EM haplotype mixture per family.
  - `feedback_o3_rna_only.md`.
  - `grep -il o3 *.md` for more.
- **Dead ends:** `/mnt/linuxdisk/home/juanfraitu/rustle_soto_subset/docs/NEGATIVE_RESULTS_REGISTER.md`. You MUST check it before
  proposing; grep for O3, missing copy, IsoCon, haplotype, samplecopy, PSV, flag.
- **Frozen outputs you may inspect** (read-only, if useful): `/mnt/linuxdisk/tmp/o3_graph/o2vc/audit_frozen_eccbc4e0/{hg002,hg00733}/`.
  - `calls.tsv`
  - `<family>.columns.tsv`: each read's columns, with every copy's ref / genotype / depth
  - `<family>.reads.tsv`: AS, NM, de, ms, subs, indels per tied alignment
  - `<family>.flags.tsv`
  - `<family>.teachers.tsv`
  - `vcf/*.vcf`: GT / DP / AF / PS

## Truth available later for evaluation
Use it only to plan the evaluation; do not tune on it here.
- **HG00733 (simulated, exact):** `/mnt/linuxdisk/tmp/o3_graph/donorsim/HG00733/sim/build.json` lists "copies without an ortholog"
  and "silent copies", i.e. reference copies the donor lacks. donorsim never simulates donor-only copies.
- **HG002 (real):** about 42 tied reads come from real donor-side insertions with no CHM13 counterpart (SPDYE1 28, RASA4B 8, USP17L11
  3, POTEM 2, BMS1P1 1). NPIPB2's 16 "no_counterpart" reads are a chain artifact, the donor's NPIPB2.2 with a 57-bp insertion.

## The vault
`/mnt/c/Users/jfris/Desktop/ThesisVault/`:
- `Papers/`: about 29 notes, each with "Method in brief";
- `Fulltext/`;
- `Topics/`: O2, O3, Phasing, SegDups, Y_Ampliconic_Genes, NPIP, TBC1D3, Terminology.

The closest papers: DuploMap, Paraphase (2023, 2025), Parascopy, SDA (Vollger 2019), Chaisson 2017, IsoCon, devider, pHapCompass,
longcallR, Clair3-RNA, Pal 2026, Dishuck 2025, Guitart 2024, Canzar 2016, Deschamps 2020.

## Output rules
- Write ONE markdown file (path given in your prompt), 80-160 lines, dense, with tables.
- Cite as `Paper_note §section` or with a short quote. Cite only what you read.
- For each candidate "fingerprint", give:
  - (a) the signal, computable from RNA reads after O2-style assignment;
  - (b) the mechanism: which missing / extra / absent-copy scenario produces it;
  - (c) literature support;
  - (d) the confounders: alleles, allele-specific expression, RNA editing, gene conversion, chimeras, liftover;
  - (e) whether it is threshold-free or combinatorial;
  - (f) how to compute it from our outputs;
  - (g) whether the register already killed it.
- Read-only everywhere else. No heavy jobs. No subagents.
- Reply in under 12 lines: the file path, plus your top 3 fingerprints in one line each.
