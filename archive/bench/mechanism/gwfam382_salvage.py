"""Task 6c: salvage GWFAM382 (V4b) -- corrected scoring against copy_assign's OWN detected
PSV alleles (intact run), NOT the broken catalog annotation for the deleted copy.

Background (see .superpowers/sdd/task-6b-report.md, "Attempt 1"): the prior GWFAM382 run
picked host=copy2, hidden(deleted)=[copy3, copy1]. copy3's CATALOG entry claims a 2.3kb span
(122769395-122771700, same size class as its 3 siblings) but copy_assign's own fresh
read-based transcript detection on the intact reference finds copy3's actual footprint is
122769398-123261551 -- a 492kb readthrough/mis-chain artifact. The prior driver
(real_tandem.py::pick_deletion_design) used copy3's QUANT.TSV span (the 492kb artifact) as the
mask span, so the "deleted copy" in the degraded reference was not a clean single tandem
paralog -- it barely touched the intended 2.3kb body inside the padded extraction window.

Fix: mask copy3 using its CATALOG core span (122769395-122771700, from GGO_gwcat.copies.tsv)
-- the real ~2.3kb copy body -- not the 492kb quant.tsv footprint. Score the admitted degraded
candidates against copy_assign's OWN intact-run PSV alleles for copy3 (psv_copies.tsv), which
is already what real_tandem.py's scoring step does (genome-pos matched via psv_cols.tsv on
both sides) -- that part of the method was already correct; only the MASK SPAN was wrong.

Design per task brief: delete ONE copy (copy3) first -- the family has 4 copies, so a
single-copy deletion should leave 3 in the reference. Whether that alone clears the hardcoded
min_clusters=3 admission gate (absent_copy.rs) is checked empirically (this locus already
shows collapsed_copies=38 in the INTACT run's families.tsv -- i.e. real, already-latent
multi-copy structure at copy2's locus even before any deletion, so it may not behave like the
"bare 1-of-N deletion never reaches 3 clusters" case sim_tandem.py's docstring describes for a
clean simulated tandem array).
"""
import json
import pathlib
import subprocess
import sys

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))
from real_tandem import (  # noqa: E402
    GGO_FASTA, GGO_BAM, run_copy_assign, load_psv, _read_tsv, identity,
    build_degraded_ref, extract_and_realign,
)
from sim_tandem import set_match  # noqa: E402

CONTIG = "NC_073239.2"
LOCUS_START, LOCUS_END = 122615409, 122771700
PAD = 3000
PAD_START, PAD_END = LOCUS_START - PAD, LOCUS_END + PAD

# copy3's TRUE catalog body span (GGO_gwcat.copies.tsv row: GWFAM382 3 DN_NC_073239.2_122769395_7
# NC_073239.2 122769395 122771700 7 - 3) -- NOT copy_assign's own 492kb quant.tsv footprint for
# this locus (122769398-123261551, a readthrough/mis-chain artifact, confirmed independently by
# both the Task 6b intact run and the fresh re-run below).
CATALOG_COPY3_SPAN = (122769395, 122771700)

WORK = pathlib.Path("/home/juanfra/winloci_scratch/mech_real/gwfam382_salvage")
INTACT_PREFIX = WORK / "gwfam382_intact"
DEGRADED_PREFIX = WORK / "gwfam382_degraded"


def main():
    WORK.mkdir(parents=True, exist_ok=True)

    # Step 2 (already run once foreground; re-check it's present, else run it).
    if not pathlib.Path(f"{INTACT_PREFIX}.families.tsv").exists():
        run_copy_assign(GGO_BAM, GGO_FASTA, f"{CONTIG}:{LOCUS_START}-{LOCUS_END}", INTACT_PREFIX,
                         extra_args=["--dump-psv"])

    quant = _read_tsv(f"{INTACT_PREFIX}.quant.tsv")
    families = _read_tsv(f"{INTACT_PREFIX}.families.tsv")
    n_copies_supported = int(families[0]["n_copies"]) if families else 0
    collapsed_copies = int(families[0]["collapsed_copies"]) if families else 0
    print(f"intact: n_copies_supported={n_copies_supported} collapsed_copies={collapsed_copies}",
          file=sys.stderr)
    for r in quant:
        print(f"  copy{r['copy_index']}: {r['copy_chrom']}:{r['copy_start']}-{r['copy_end']} "
              f"n_reads_hard={r['n_reads_hard']} anchored={r['anchored_reads']}", file=sys.stderr)

    cols, copies, tids = load_psv(str(INTACT_PREFIX))

    # Step 3: mask ONLY copy3's TRUE catalog body (2.3kb), keep copies 0,1,2 untouched in the
    # degraded reference.
    hidden = [3]
    keep = [0, 1, 2]
    mask_spans = [CATALOG_COPY3_SPAN]
    print(f"Design: keep={keep} hidden(deleted)={hidden} mask_spans={mask_spans} "
          f"(catalog span, NOT quant.tsv's 492kb artifact span)", file=sys.stderr)

    # Step 4: degrade + realign.
    degraded_fa = WORK / "degraded_ref.fa"
    new_contig = f"{CONTIG}_degraded"
    build_degraded_ref(GGO_FASTA, CONTIG, PAD_START, PAD_END, mask_spans, degraded_fa, new_contig)
    sorted_bam, n_reads = extract_and_realign(GGO_BAM, CONTIG, PAD_START, PAD_END, degraded_fa,
                                               WORK, "gwfam382_degraded_reads")
    degraded_region = f"{new_contig}:1-{PAD_END - PAD_START + 1}"
    run_copy_assign(sorted_bam, degraded_fa, degraded_region, DEGRADED_PREFIX,
                     extra_args=["--absent-copies", "--linearize", "--dump-psv"])

    d_families = _read_tsv(f"{DEGRADED_PREFIX}.families.tsv")
    dna_needs = _read_tsv(f"{DEGRADED_PREFIX}.dna_needs.tsv") \
        if pathlib.Path(f"{DEGRADED_PREFIX}.dna_needs.tsv").exists() else []
    d_n_copies = int(d_families[0]["n_copies"]) if d_families else 0
    d_collapsed = int(d_families[0]["collapsed_copies"]) if d_families else 0
    print(f"degraded: n_copies={d_n_copies} collapsed_copies={d_collapsed} "
          f"n_dna_needs={len(dna_needs)}", file=sys.stderr)

    n_reference_kept = len(keep)
    min_clusters_satisfied = d_n_copies > n_reference_kept

    result = {
        "locus": f"{CONTIG}:{LOCUS_START}-{LOCUS_END}",
        "family_id": "GWFAM382",
        "substrate": "GGO.fasta + GGO_mm.bam, gorilla-to-gorilla",
        "truth_source": "intact-run detected PSV alleles (genome-pos-matched), NOT catalog annotation",
        "n_copies_supported_intact": n_copies_supported,
        "collapsed_copies_intact": collapsed_copies,
        "design": {"kept_reference_copy_indices": keep, "hidden_deleted_indices": hidden,
                   "mask_spans_catalog": mask_spans,
                   "note": "mask span = catalog core body (2.3kb), NOT copy_assign's own "
                           "quant.tsv footprint for copy3 (492kb, readthrough/mis-chain artifact)"},
        "n_reads_realigned": n_reads,
        "degraded_n_copies": d_n_copies,
        "degraded_collapsed_copies": d_collapsed,
        "n_dna_needs": len(dna_needs),
        "dna_needs_reasons": [r.get("reason", "") for r in dna_needs],
        "min_clusters_satisfied": min_clusters_satisfied,
    }

    recovered = False
    identity_to_intact = None
    n_psv_cols_matched = 0
    match_detail = []
    all_admitted_scores = []
    if pathlib.Path(f"{DEGRADED_PREFIX}.psv_copies.tsv").exists():
        d_cols, d_copies, d_tids = load_psv(str(DEGRADED_PREFIX))
        offset = PAD_START - 1
        d_abs_cols = {ci: (pos + offset) for ci, pos in d_cols.items()}
        admitted_idxs = sorted(i for i in d_copies if i >= n_reference_kept)

        intact_pos_to_col = {pos: ci for ci, pos in cols.items()}
        shared_positions = sorted(set(p for p in d_abs_cols.values() if p in intact_pos_to_col))
        n_psv_cols_matched = len(shared_positions)

        def truth_allele_str(copy_idx):
            return "".join(copies[copy_idx][intact_pos_to_col[p]] for p in shared_positions)

        d_pos_to_col = {v: k for k, v in d_abs_cols.items()}

        def degraded_allele_str(copy_idx):
            return "".join(d_copies[copy_idx][d_pos_to_col[p]] for p in shared_positions)

        truth_seqs = [truth_allele_str(h) for h in hidden]
        recovered_seqs = [degraded_allele_str(i) for i in admitted_idxs]

        # score EVERY admitted candidate against the true deleted copy (not just the greedy
        # best-match set) so spurious/non-matching admissions are visible, not hidden.
        for ai, aidx in zip(admitted_idxs, recovered_seqs):
            idv, n = identity(aidx, truth_seqs[0]) if truth_seqs else (0.0, 0)
            all_admitted_scores.append({"admitted_copy_index": ai, "identity_to_copy3": idv, "n_cols": n})

        matches = set_match(recovered_seqs, truth_seqs) if recovered_seqs and truth_seqs else []
        match_detail = matches
        if matches:
            identity_to_intact = min(m[2] for m in matches)
            recovered = len(admitted_idxs) > 0 and identity_to_intact is not None and identity_to_intact > 0.90

    n_real_matches = sum(1 for s in all_admitted_scores if s["identity_to_copy3"] > 0.90)
    n_spurious = len(all_admitted_scores) - n_real_matches

    result.update({
        "n_admitted": len(all_admitted_scores),
        "n_psv_cols_matched": n_psv_cols_matched,
        "identity_to_intact": identity_to_intact,
        "recovered": recovered,
        "set_match": match_detail,
        "all_admitted_scores": all_admitted_scores,
        "n_real_matches_gt_0.90": n_real_matches,
        "n_spurious_admissions": n_spurious,
        "identity_resolution": "PSV-column (matched by absolute genome_pos between intact and degraded runs)",
    })

    out_path = WORK / "gwfam382_salvage_result.json"
    out_path.write_text(json.dumps(result, indent=2))
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
