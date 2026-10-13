#!/usr/bin/env python3
"""psv_population.tsv (spec §2): per (chrom, pos) site, the alignment verdict where it exists, else the k-mer verdict.

    table.py PREFIX CATALOG_COPIES_TSV GENOME_FA KMER_STATUS ALIGN_STATUS OUT.tsv
"""
import os, sys
sys.path.insert(0, os.path.dirname(__file__))
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "psv_sun"))
import kmer_layer as KL

# Stricter status wins when a site belongs to several columns.
RANK = {"POLYMORPHIC": 0, "UNCOVERED": 1, "NOT_ASSESSED": 2, "FIXED": 3}


def combine(krow, arow):
    if krow["status"] == "NOT_ASSESSED":
        return "NOT_ASSESSED", "-"
    if arow is not None and arow["status"] in ("FIXED", "POLYMORPHIC"):
        return arow["status"], "A"
    return krow["status"], "K"


def merge_sites(cols, km, al):
    """cols: {(fam, col): [Site]} as from kmer_layer.columns; km / al: {(fam, int col): TSV row dict}, al may lack keys.
    -> {(chrom, pos): (status, method, (fam, col), n_men, n_cover, n_carry, n_mismatch)}.
    A site in several columns takes the stricter status; ties break on (fam, col), so input order never matters."""
    best = {}
    for key, sites in cols.items():
        st, meth = combine(km[key], al.get(key))
        a = al.get(key, {})
        row = (st, meth, key, km[key]["n_men"], km[key]["n_cover"], km[key]["n_carry_all"],
               a.get("n_mismatch", "0"))
        for _, c, p, *_ in sites:
            rank = (RANK[st], key)
            if (c, p) not in best or rank < best[(c, p)][0]:
                best[(c, p)] = (rank, row)
    return {site: row for site, (_, row) in best.items()}


def main(prefix, catalog_tsv, genome_fa, kmer_tsv, align_tsv, out):
    import psv_sun as P
    F = P.Fasta(genome_fa)
    genome = lambda c, s, e: F.fetch(c, s, e)
    km = {(r["family_id"], int(r["col_index"])): r for r in KL.tsv(kmer_tsv)}
    if not os.path.exists(align_tsv):
        sys.exit(f"no alignment table {align_tsv}; refusing to build a k-mer-only population table")
    al = {(r["family_id"], int(r["col_index"])): r for r in KL.tsv(align_tsv)}
    cols = KL.columns(prefix, catalog_tsv, genome)
    site = merge_sites(cols, km, al)
    with open(out, "w") as fh:
        fh.write("chrom\tpos\tfamily_id\tcol_index\tstatus\tmethod\tn_men\tn_cover\tn_carry\tn_mismatch\n")
        for (c, p), (st, meth, (fam, col), n, nc, ncar, nm) in sorted(site.items()):
            fh.write(f"{c}\t{p}\t{fam}\t{col}\t{st}\t{meth}\t{n}\t{nc}\t{ncar}\t{nm}\n")
    print(f"{len(site)} sites -> {out}")


if __name__ == "__main__":
    main(*sys.argv[1:7])
