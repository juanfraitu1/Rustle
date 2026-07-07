#!/usr/bin/env python3
"""bundle_isoform_probe.py — do bundle-style assemblers capture more isoforms in multi-copy families?

Compares the shipped exact-intron-chain de-novo catalog against StringTie's bundle-based
output for the 8 canonical known multi-copy families. For each family we count how many
transcripts / unique intron chains overlap the family's genomic footprint in each assembly.

Run: PYTHONHASHSEED=0 /home/juanfra/miniforge3/bin/python bench/bundle_isoform_probe.py
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
OUT_TSV = os.path.join(HERE, "bundle_isoform_probe.tsv")

PAD = 5_000  # family region padding


def load_showcase_families():
    """Map canonical family name -> catalog family_id."""
    out = {}
    with open(SHOWCASE_TSV) as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            out[r["family"]] = int(r["catalog_fid"])
    return out


def load_refine_members(fam_ids):
    """{family_id: [(member_dn, chrom, start, end), ...]}"""
    out = collections.defaultdict(list)
    with open(REFINE_TSV) as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            fi = int(r["family_id"])
            if fi in fam_ids:
                out[fi].append((r["member_dn"], r["chrom"], int(r["start"]), int(r["end"])))
    return out


def load_meta_isoforms():
    """Exact-chain isoforms from de-novo meta table: list of (chrom,start,end,n_exon)."""
    out = []
    with open(META_TSV) as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            out.append((r["chrom"], int(r["start"]), int(r["end"]), int(r["n_exon"])))
    return out


def load_meta_with_introns():
    """Return {id: (chrom,start,end,n_exon, intron_tuple)} for exact-chain isoforms."""
    # skeletons keyed by (chrom,start,end,n_exon,n_reads)
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
            out[tid] = (r["chrom"], int(r["start"]), int(r["end"]), int(r["n_exon"]), skels.get(key, ()))
    return out


def parse_stringtie_isoforms():
    """Parse StringTie GTF into list of (chrom,start,end,intron_tuple)."""
    transcripts = {}  # transcript_id -> {chrom, start, end, exons[]}
    with open(STRINGTIE_GTF) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            p = line.rstrip("\n").split("\t")
            if len(p) < 9:
                continue
            chrom, typ, start, end, strand = p[0], p[2], int(p[3]), int(p[4]), p[6]
            attrs = p[8]
            # extract transcript_id
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

    out = []
    for tid, t in transcripts.items():
        exons = sorted(t["exons"])
        introns = []
        for i in range(len(exons) - 1):
            introns.append((exons[i][1] + 1, exons[i + 1][0] - 1))
        out.append((t["chrom"], t["start"], t["end"], tuple(introns)))
    return out


def family_region(members):
    """Return per-chromosome (chrom, start-pad, end+pad) regions for a family."""
    by_chrom = collections.defaultdict(lambda: [None, None])
    for _, chrom, s, e in members:
        if by_chrom[chrom][0] is None or s < by_chrom[chrom][0]:
            by_chrom[chrom][0] = s
        if by_chrom[chrom][1] is None or e > by_chrom[chrom][1]:
            by_chrom[chrom][1] = e
    return [(c, max(0, s - PAD), e + PAD) for c, (s, e) in by_chrom.items()]


def count_in_regions(isoforms, regions):
    """Count isoforms overlapping any of the regions, plus unique intron chains."""
    # regions: list of (chrom, start, end)
    by_chrom = collections.defaultdict(list)
    for c, s, e in regions:
        by_chrom[c].append((s, e))

    matched = []
    for chrom, start, end, n_exon, introns in isoforms:
        for rs, re in by_chrom.get(chrom, []):
            if start <= re and end >= rs:  # overlap
                matched.append((chrom, start, end, n_exon, introns))
                break

    tx_count = len(matched)
    chains = set((chrom, introns) for chrom, _, _, _, introns in matched)
    chain_count = len(chains)
    return tx_count, chain_count


def isoforms_in_regions(isoforms, regions):
    """Return list of isoforms overlapping any region."""
    by_chrom = collections.defaultdict(list)
    for c, s, e in regions:
        by_chrom[c].append((s, e))
    out = []
    for chrom, start, end, n_exon, introns in isoforms:
        for rs, re in by_chrom.get(chrom, []):
            if start <= re and end >= rs:
                out.append((chrom, start, end, n_exon, introns))
                break
    return out


def overlaps(a, b, min_frac=0.5):
    """a,b are (chrom,s,e). Reciprocal overlap >= min_frac on either side."""
    if a[0] != b[0]:
        return False
    s1, e1 = a[1], a[2]
    s2, e2 = b[1], b[2]
    ov = max(0, min(e1, e2) - max(s1, s2))
    len1 = e1 - s1
    len2 = e2 - s2
    if len1 == 0 or len2 == 0:
        return False
    return ov / len1 >= min_frac or ov / len2 >= min_frac


def shares_intron(a, b):
    """a,b are (chrom,s,e,n_exon,introns). True if they share any intron."""
    if a[0] != b[0]:
        return False
    return bool(set(a[4]) & set(b[4]))


def matched_isoforms(query_list, target_list, frac=0.5):
    """For each query isoform, return True if matched to any target by overlap or shared intron."""
    matched = 0
    for q in query_list:
        for t in target_list:
            if overlaps(q, t, frac) or shares_intron(q, t):
                matched += 1
                break
    return matched


def missed_isoforms(query_list, target_list, frac=0.5):
    """Return isoforms in query_list not matched to any target."""
    out = []
    for q in query_list:
        found = False
        for t in target_list:
            if overlaps(q, t, frac) or shares_intron(q, t):
                found = True
                break
        if not found:
            out.append(q)
    return out


def main():
    showcase = load_showcase_families()
    fam_ids = set(showcase.values())
    members = load_refine_members(fam_ids)

    # exact-chain isoforms
    exact = []
    for tid, (chrom, s, e, nx, introns) in load_meta_with_introns().items():
        exact.append((chrom, s, e, nx, introns))

    # StringTie bundle isoforms
    stringtie = []
    for chrom, s, e, introns in parse_stringtie_isoforms():
        stringtie.append((chrom, s, e, len(introns) + 1, introns))

    rows = []
    print(f"{'family':<10} {'fid':>5} {'regions':>7} {'exact_tx':>9} {'exact_ch':>9} {'st_tx':>8} {'st_ch':>8} {'st_missed':>10} {'exact_missed':>13}")
    print("-" * 95)
    for name, fid in sorted(showcase.items(), key=lambda x: x[1]):
        regs = family_region(members[fid])
        exact_list = isoforms_in_regions(exact, regs)
        st_list = isoforms_in_regions(stringtie, regs)
        etx, ech = len(exact_list), len(set((c, introns) for c, _, _, _, introns in exact_list))
        stx, sch = len(st_list), len(set((c, introns) for c, _, _, _, introns in st_list))
        # isoforms missed by the other method (overlap- or intron-based match)
        st_missed = stx - matched_isoforms(st_list, exact_list, frac=0.5)
        exact_missed = etx - matched_isoforms(exact_list, st_list, frac=0.5)
        rows.append((name, fid, len(regs), etx, ech, stx, sch, st_missed, exact_missed))
        print(f"{name:<10} {fid:>5} {len(regs):>7} {etx:>9} {ech:>9} {stx:>8} {sch:>8} {st_missed:>10} {exact_missed:>13}")

    print("-" * 95)
    tot_etx = sum(r[3] for r in rows); tot_ech = sum(r[4] for r in rows)
    tot_stx = sum(r[5] for r in rows); tot_sch = sum(r[6] for r in rows)
    tot_st_missed = sum(r[7] for r in rows); tot_exact_missed = sum(r[8] for r in rows)
    print(f"{'TOTAL':<10} {'':>5} {sum(r[2] for r in rows):>7} {tot_etx:>9} {tot_ech:>9} {tot_stx:>8} {tot_sch:>8} {tot_st_missed:>10} {tot_exact_missed:>13}")

    # Detailed breakdown of isoforms missed by exact-chain (captured by StringTie bundles)
    print("\n--- exact-chain isoforms missed by StringTie bundles (by family, exon-count) ---")
    for name, fid in sorted(showcase.items(), key=lambda x: x[1]):
        regs = family_region(members[fid])
        exact_list = isoforms_in_regions(exact, regs)
        st_list = isoforms_in_regions(stringtie, regs)
        missed = missed_isoforms(exact_list, st_list, frac=0.5)
        if not missed:
            continue
        by_exon = collections.Counter(x[3] for x in missed)
        print(f"{name}: n={len(missed)}  " + "  ".join(f"{k}ex:{v}" for k, v in sorted(by_exon.items())))

    print("\n--- StringTie bundle isoforms missed by exact-chain (by family, exon-count) ---")
    for name, fid in sorted(showcase.items(), key=lambda x: x[1]):
        regs = family_region(members[fid])
        exact_list = isoforms_in_regions(exact, regs)
        st_list = isoforms_in_regions(stringtie, regs)
        missed = missed_isoforms(st_list, exact_list, frac=0.5)
        if not missed:
            continue
        by_exon = collections.Counter(x[3] for x in missed)
        print(f"{name}: n={len(missed)}  " + "  ".join(f"{k}ex:{v}" for k, v in sorted(by_exon.items())))

    # Fragmentation: how many exact-chain isoforms match each StringTie isoform
    print("\n--- exact-chain fragmentation of StringTie isoforms ---")
    for name, fid in sorted(showcase.items(), key=lambda x: x[1]):
        regs = family_region(members[fid])
        exact_list = isoforms_in_regions(exact, regs)
        st_list = isoforms_in_regions(stringtie, regs)
        counts = []
        for s in st_list:
            n = sum(1 for e in exact_list if overlaps(s, e, 0.5) or shares_intron(s, e))
            counts.append(n)
        if counts:
            avg = sum(counts) / len(counts)
            print(f"{name}: ST isoforms={len(counts)}  avg exact matches per ST={avg:.2f}  max={max(counts)}  unmatched={sum(1 for c in counts if c == 0)}")

    with open(OUT_TSV, "w") as o:
        o.write("family\tcatalog_fid\tn_regions\texact_tx\texact_chains\tstringtie_tx\tstringtie_chains\tstringtie_missed_by_exact\texact_missed_by_stringtie\n")
        for row in rows:
            o.write("\t".join(str(x) for x in row) + "\n")
    print(f"\n[+] wrote {OUT_TSV}")


if __name__ == "__main__":
    main()
