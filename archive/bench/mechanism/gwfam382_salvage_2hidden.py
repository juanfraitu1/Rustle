"""Task 6c, second design: directly re-run Attempt 1's exact design (host=copy2,
hidden=[copy3, copy1]) with the CORRECTED mask span for copy3 (catalog 2.3kb core body,
122769395-122771700, instead of copy_assign's own 492kb quant.tsv readthrough footprint for
that locus). copy1's span is unchanged (catalog and quant.tsv already agree exactly:
122634082-122635806), so only copy3's mask span differs from Attempt 1.

Purpose: this is the most direct, apples-to-apples correction of the original reported
identity_to_intact=0.394-0.408 (Attempt 1, task-6b-report.md) -- same design, same host, same
2 hidden copies, only the broken mask span fixed. Run AFTER gwfam382_salvage.py (single-copy
design) and only once that process has exited (serial execution, WSL2 crash rule).
"""
import json
import pathlib
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

CATALOG_COPY1_SPAN = (122634082, 122635806)  # matches quant.tsv exactly -- unchanged from Attempt 1
CATALOG_COPY3_SPAN = (122769395, 122771700)  # CORRECTED: catalog span, not the 492kb quant artifact

WORK = pathlib.Path("/home/juanfra/winloci_scratch/mech_real/gwfam382_salvage")
INTACT_PREFIX = WORK / "gwfam382_intact"
DEGRADED_PREFIX = WORK / "gwfam382_degraded_2hidden"


def main():
    quant = _read_tsv(f"{INTACT_PREFIX}.quant.tsv")
    families = _read_tsv(f"{INTACT_PREFIX}.families.tsv")
    n_copies_supported = int(families[0]["n_copies"]) if families else 0
    cols, copies, tids = load_psv(str(INTACT_PREFIX))

    # Attempt 1's exact design: keep=[0,2], host=2, hidden=[3,1].
    keep = [0, 2]
    host = 2
    hidden = [3, 1]
    mask_spans = [CATALOG_COPY3_SPAN, CATALOG_COPY1_SPAN]
    print(f"Design (Attempt-1-matched): keep={keep} host={host} hidden(deleted)={hidden} "
          f"mask_spans={mask_spans} (copy3 span CORRECTED to catalog; copy1 span unchanged)",
          file=sys.stderr)

    degraded_fa = WORK / "degraded_ref_2hidden.fa"
    new_contig = f"{CONTIG}_degraded2"
    build_degraded_ref(GGO_FASTA, CONTIG, PAD_START, PAD_END, mask_spans, degraded_fa, new_contig)
    sorted_bam, n_reads = extract_and_realign(GGO_BAM, CONTIG, PAD_START, PAD_END, degraded_fa,
                                               WORK, "gwfam382_degraded2_reads")
    degraded_region = f"{new_contig}:1-{PAD_END - PAD_START + 1}"
    run_copy_assign(sorted_bam, degraded_fa, degraded_region, DEGRADED_PREFIX,
                     extra_args=["--absent-copies", "--linearize", "--dump-psv"])

    d_families = _read_tsv(f"{DEGRADED_PREFIX}.families.tsv")
    dna_needs = _read_tsv(f"{DEGRADED_PREFIX}.dna_needs.tsv") \
        if pathlib.Path(f"{DEGRADED_PREFIX}.dna_needs.tsv").exists() else []
    d_n_copies = int(d_families[0]["n_copies"]) if d_families else 0
    d_collapsed = int(d_families[0]["collapsed_copies"]) if d_families else 0
    print(f"degraded(2hidden): n_copies={d_n_copies} collapsed_copies={d_collapsed} "
          f"n_dna_needs={len(dna_needs)}", file=sys.stderr)

    n_reference_kept = len(keep)
    min_clusters_satisfied = d_n_copies > n_reference_kept

    result = {
        "locus": f"{CONTIG}:{LOCUS_START}-{LOCUS_END}",
        "family_id": "GWFAM382",
        "design_label": "Attempt-1-matched (host=copy2, hidden=[copy3,copy1]), copy3 mask span corrected",
        "n_copies_supported_intact": n_copies_supported,
        "design": {"kept_reference_copy_indices": keep, "host_index": host,
                   "hidden_deleted_indices": hidden, "mask_spans_catalog": mask_spans},
        "n_reads_realigned": n_reads,
        "degraded_n_copies": d_n_copies,
        "degraded_collapsed_copies": d_collapsed,
        "n_dna_needs": len(dna_needs),
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

        truth_seqs = [truth_allele_str(h) for h in hidden]  # [copy3, copy1]
        recovered_seqs = [degraded_allele_str(i) for i in admitted_idxs]

        for aidx, aseq in zip(admitted_idxs, recovered_seqs):
            best_id, best_h = -1.0, None
            for h, tseq in zip(hidden, truth_seqs):
                idv, n = identity(aseq, tseq)
                if idv > best_id:
                    best_id, best_h = idv, h
            all_admitted_scores.append({"admitted_copy_index": aidx, "best_match_hidden_copy": best_h,
                                         "identity": best_id})

        matches = set_match(recovered_seqs, truth_seqs) if recovered_seqs and truth_seqs else []
        match_detail = matches
        if matches:
            identity_to_intact = min(m[2] for m in matches)
            recovered = len(admitted_idxs) > 0 and identity_to_intact is not None and identity_to_intact > 0.90

    n_real_matches = sum(1 for s in all_admitted_scores if s["identity"] > 0.90)
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
    })

    out_path = WORK / "gwfam382_salvage_2hidden_result.json"
    out_path.write_text(json.dumps(result, indent=2))
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
