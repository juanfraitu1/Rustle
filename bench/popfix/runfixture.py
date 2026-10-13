"""Synthetic run directories for the popfix scorers' tests (a copy_assign --dump-star run in miniature)."""
import os

STAR_HEAD = "read_name\tfamily_id\tstatus\tassigned_copy\tcatalog_copy_idx\tn_candidates\tcandidates\tn_cols\tcolumns"
ASG_HEAD = ("read_name\tfamily_id\tassigned_copy\tstatus\tn_decisive\tmargin\tp_value\tmin_p_value\tas_best\tas_second\t"
            "as_margin\tas_per_base_best\tas_per_base_2nd\tin_copy\tcatalog_copy_idx\torigin_rejected\tn_candidates\t"
            "sole_candidate\tcontested\treadthrough_into\tprimary_local\ttie_outside_catalog\taligner_disagreement")
TIER_HEAD = "read_name\tfamily_id\tstatus_all\tcopy_all\tstatus_fixed\tcopy_fixed\tn_dec_fixed\ttier"
JOIN_HEAD = "family_id\tcopy_index\tcopy_tid\tcatalog_family_id\tcatalog_copy_idx\tchrom\tstart\tend\tn_reads_hard"
POP_HEAD = "chrom\tpos\tfamily_id\tcol_index\tstatus\tmethod\tn_men\tn_cover\tn_carry\tn_mismatch"


def write_run(d, stars, tiers=None, joins=None, positions=True, lengths=None, prefix="run"):
    """stars: list of dicts read, fam, status, copy, cand (list), cols [(readpos, obs, alleles, [pos per cand or None])],
    optional origin_rejected / tie_outside / contested. Returns (run prefix, lengths path)."""
    p = os.path.join(d, prefix)
    with open(p + ".star_reads.tsv", "w") as fh, open(p + ".assignments.tsv", "w") as ah:
        fh.write(STAR_HEAD + ("\tpositions" if positions else "") + "\n")
        ah.write(ASG_HEAD + "\n")
        for s in stars:
            cols = ",".join(f"{rp}:{o}:{al}" for rp, o, al, _ in s["cols"])
            pos = ",".join("/".join("." if x is None else str(x) for x in gp) for _, _, _, gp in s["cols"])
            fields = [s["read"], s["fam"], s["status"], str(s["copy"]), str(s["copy"]), str(len(s["cand"])),
                      ",".join(map(str, s["cand"])), str(len(s["cols"])), cols]
            if positions:
                fields.append(pos)
            fh.write("\t".join(fields) + "\n")
            ah.write("\t".join([s["read"], s["fam"], str(s["copy"]), s["status"], "0", "0", "1", "1", "0", "0", "0", "0",
                                "0", "true", str(s["copy"]), str(s.get("origin_rejected", 0)), str(len(s["cand"])), "0",
                                str(s.get("contested", 1)), "-", "1", str(s.get("tie_outside", 0)), "0"]) + "\n")
    if joins is not None:
        with open(p + ".family_join.tsv", "w") as jh:
            jh.write(JOIN_HEAD + "\n")
            for fam, idx, chrom in joins:
                jh.write(f"{fam}\t{idx}\tT{fam}_{idx}\t{fam}\t{idx}\t{chrom}\t1\t2\t0\n")
    if tiers is not None:
        with open(p + ".population_tier.tsv", "w") as th:
            th.write(TIER_HEAD + "\n")
            for t in tiers:
                th.write("\t".join(map(str, t)) + "\n")
    lp = os.path.join(d, "lengths.tsv")
    with open(lp, "w") as lh:
        for name, l in (lengths or {}).items():
            lh.write(f"{name}\t{l}\n")
    return p, lp


def write_pop(path, rows):
    """rows: (chrom, pos, status)."""
    with open(path, "w") as fh:
        fh.write(POP_HEAD + "\n")
        for i, (c, pos, st) in enumerate(rows):
            fh.write(f"{c}\t{pos}\tF\t{i}\t{st}\t{'K' if st != 'NOT_ASSESSED' else '-'}\t48\t48\t48\t0\n")


def write_pop_rows(path, rows):
    """rows: (chrom, pos, family_id, status, method)."""
    with open(path, "w") as fh:
        fh.write(POP_HEAD + "\n")
        for i, (c, pos, fam, st, method) in enumerate(rows):
            fh.write(f"{c}\t{pos}\t{fam}\t{i}\t{st}\t{method}\t48\t48\t48\t0\n")
