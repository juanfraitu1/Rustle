#!/usr/bin/env python3
"""hybrid_isoform_assembly.py — StringTie base + exact-chain supplement in multi-copy families.

Builds a hybrid transcript set for the canonical known multi-copy families:
  1. Start with StringTie isoforms in each family region (bundle-based, clean single-copy models).
  2. Add exact-intron-chain de-novo isoforms that StringTie missed (overlap/intron match).

This is the "finds what StringTie finds but finds more" design for multi-copy regions.

Run: PYTHONHASHSEED=0 /home/juanfra/miniforge3/bin/python bench/hybrid_isoform_assembly.py
"""
import csv
import os
import collections

HERE = os.path.dirname(os.path.abspath(__file__))
SCRATCH = "/home/juanfra/winloci_scratch"

META_TSV = os.path.join(SCRATCH, "denovo_transcripts.meta.tsv")
SKEL_TSV = os.path.join(SCRATCH, "denovo_skeletons.tsv")
STRINGTIE_GTF = os.path.join(SCRATCH, "genome_st.gtf")
REFINE_TSV = os.path.join(HERE, "family_rna_refine.tsv")
SHOWCASE_TSV = os.path.join(HERE, "known_family_showcase.tsv")
OUT_GTF = os.path.join(HERE, "hybrid_isoforms_8families.gtf")
OUT_TSV = os.path.join(HERE, "hybrid_isoform_assembly.tsv")

PAD = 5_000


def load_showcase_families():
    out = {}
    with open(SHOWCASE_TSV) as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            out[r["family"]] = int(r["catalog_fid"])
    return out


def load_refine_members(fam_ids):
    out = collections.defaultdict(list)
    with open(REFINE_TSV) as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            fi = int(r["family_id"])
            if fi in fam_ids:
                out[fi].append((r["member_dn"], r["chrom"], int(r["start"]), int(r["end"])))
    return out


def load_meta_with_introns():
    skels = {}
    with open(SKEL_TSV) as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            key = (r["chrom"], int(r["start"]), int(r["end"]), int(r["n_exon"]), int(r["n_reads"]))
            introns = []
            if r["introns"].strip():
                for tok in r["introns"].split(";"):
                    a, b = tok.split("-")
                    introns.append((int(a), int(b)))
            skels[key] = tuple(introns)

    out = {}
    with open(META_TSV) as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            tid = r["id"]
            key = (r["chrom"], int(r["start"]), int(r["end"]), int(r["n_exon"]), int(r["n_reads"]))
            out[tid] = (r["chrom"], int(r["start"]), int(r["end"]), int(r["n_exon"]), int(r["n_reads"]), skels.get(key, ()), r["strand"])
    return out


def parse_stringtie_isoforms():
    """Parse StringTie GTF into dict transcript_id -> {chrom,start,end,strand,exons}."""
    transcripts = {}
    with open(STRINGTIE_GTF) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            p = line.rstrip("\n").split("\t")
            if len(p) < 9:
                continue
            chrom, typ, start, end, strand = p[0], p[2], int(p[3]), int(p[4]), p[6]
            attrs = p[8]
            tid = None
            for tok in attrs.split(";"):
                tok = tok.strip()
                if tok.startswith("transcript_id"):
                    tid = tok.split('"')[1]
                    break
            if tid is None:
                continue
            if typ == "transcript":
                transcripts[tid] = {"chrom": chrom, "start": start, "end": end,
                                    "strand": strand, "exons": []}
            elif typ == "exon" and tid in transcripts:
                transcripts[tid]["exons"].append((start, end))
    return transcripts


def family_region(members):
    by_chrom = collections.defaultdict(lambda: [None, None])
    for _, chrom, s, e in members:
        if by_chrom[chrom][0] is None or s < by_chrom[chrom][0]:
            by_chrom[chrom][0] = s
        if by_chrom[chrom][1] is None or e > by_chrom[chrom][1]:
            by_chrom[chrom][1] = e
    return [(c, max(0, s - PAD), e + PAD) for c, (s, e) in by_chrom.items()]


def in_region(chrom, start, end, regions):
    for rc, rs, re in regions:
        if chrom == rc and start <= re and end >= rs:
            return True
    return False


def overlaps(a, b, frac=0.5):
    if a[0] != b[0]:
        return False
    s1, e1 = a[1], a[2]
    s2, e2 = b[1], b[2]
    ov = max(0, min(e1, e2) - max(s1, s2))
    len1 = e1 - s1
    len2 = e2 - s2
    if len1 == 0 or len2 == 0:
        return False
    return ov / len1 >= frac or ov / len2 >= frac


def shares_intron(a, b):
    if a[0] != b[0]:
        return False
    return bool(set(a[5]) & set(b[5]))


def matches_any(query, targets, frac=0.5):
    for t in targets:
        if overlaps(query, t, frac) or shares_intron(query, t):
            return True
    return False


def exons_from_introns(chrom, start, end, introns, strand):
    """Reconstruct exon intervals from transcript span and introns."""
    exons = []
    cur = start
    for d, a in sorted(introns):
        exons.append((cur, d - 1))
        cur = a + 1
    exons.append((cur, end))
    return exons


def emit_gtf_transcript(o, chrom, start, end, strand, gene_id, tx_id, exons, source_tag):
    o.write(f"{chrom}\t{source_tag}\ttranscript\t{start}\t{end}\t1000\t{strand}\t.\t"
            f"gene_id \"{gene_id}\"; transcript_id \"{tx_id}\";\n")
    for i, (s, e) in enumerate(exons, 1):
        o.write(f"{chrom}\t{source_tag}\texon\t{s}\t{e}\t1000\t{strand}\t.\t"
                f"gene_id \"{gene_id}\"; transcript_id \"{tx_id}\"; exon_number \"{i}\";\n")


def main():
    showcase = load_showcase_families()
    fam_ids = set(showcase.values())
    members = load_refine_members(fam_ids)
    meta = load_meta_with_introns()
    stringtie = parse_stringtie_isoforms()

    rows = []
    total_added = 0
    total_st = 0

    with open(OUT_GTF, "w") as o:
        o.write("##gff-version 2\n")
        o.write("# hybrid isoforms: StringTie base + exact-chain supplement in 8 known multi-copy families\n")

        for name, fid in sorted(showcase.items(), key=lambda x: x[1]):
            regs = family_region(members[fid])

            # StringTie isoforms in region
            st_in_region = []
            for tid, t in stringtie.items():
                if in_region(t["chrom"], t["start"], t["end"], regs):
                    st_in_region.append((tid, t))

            # exact-chain isoforms in region
            exact_in_region = []
            for tid, (chrom, s, e, nx, nr, introns, strand) in meta.items():
                if in_region(chrom, s, e, regs):
                    exact_in_region.append((tid, chrom, s, e, nx, nr, introns, strand))

            # Add exact-chain isoforms not already represented by StringTie
            added = []
            for tid, chrom, s, e, nx, nr, introns, strand in exact_in_region:
                q = (chrom, s, e, nx, nr, introns, strand)
                st_list = []
                for _, t in st_in_region:
                    ex = sorted(t["exons"])
                    introns = tuple((ex[i][1] + 1, ex[i + 1][0] - 1) for i in range(len(ex) - 1))
                    st_list.append((t["chrom"], t["start"], t["end"], 0, 0, introns, t["strand"]))
                if not matches_any(q, st_list, frac=0.5):
                    added.append((tid, chrom, s, e, nx, introns, strand))

            gene_base = f"HYBRID.{fid}"
            tx_idx = 1

            # write StringTie isoforms
            for tid, t in st_in_region:
                gene_id = f"{gene_base}.ST"
                tx_id = f"{gene_base}.ST.{tx_idx}"
                tx_idx += 1
                emit_gtf_transcript(o, t["chrom"], t["start"], t["end"], t["strand"],
                                    gene_id, tx_id, sorted(t["exons"]), "StringTie")

            # write added exact-chain isoforms
            for tid, chrom, s, e, nx, introns, strand in added:
                gene_id = f"{gene_base}.EX"
                tx_id = f"{gene_base}.EX.{tx_idx}"
                tx_idx += 1
                exons = exons_from_introns(chrom, s, e, introns, strand)
                emit_gtf_transcript(o, chrom, s, e, strand, gene_id, tx_id, exons, "ExactChain")

            rows.append((name, fid, len(st_in_region), len(added), len(st_in_region) + len(added)))
            total_st += len(st_in_region)
            total_added += len(added)

    with open(OUT_TSV, "w") as o:
        o.write("family\tcatalog_fid\tstringtie_tx\texact_chain_added\thybrid_total\n")
        for row in rows:
            o.write("\t".join(str(x) for x in row) + "\n")

    print(f"{'family':<10} {'fid':>5} {'st_tx':>8} {'added':>8} {'hybrid':>8}")
    print("-" * 45)
    for row in rows:
        print(f"{row[0]:<10} {row[1]:>5} {row[2]:>8} {row[3]:>8} {row[4]:>8}")
    print("-" * 45)
    print(f"{'TOTAL':<10} {'':>5} {total_st:>8} {total_added:>8} {total_st + total_added:>8}")
    print(f"\n[+] wrote {OUT_GTF}")
    print(f"[+] wrote {OUT_TSV}")


if __name__ == "__main__":
    main()
