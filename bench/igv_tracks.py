#!/usr/bin/env python
"""IGV tracks for the copy-assignment result — load these next to the BAM in IGV.

Produces, from a `copy_assign` run + the BAM it swept:
  <out>.tagged.bam (+ .bai) — the reads in the swept regions, each tagged with its assigned copy:
       cp:Z:<family>_c<idx>   (the copy id; IGV: right-click -> "Group by tag" / "Color by tag" -> cp)
       YC:Z:R,G,B             (per-read colour; IGV colours each read by its copy AUTOMATICALLY)
       HP:i:<idx+1>           (so IGV's haplotype colouring also splits the copies)
     A read's SECONDARY placements get the SAME colour, so you SEE the multimapping reads fan out across
     the copies. Tied / ambiguous reads get their own grey, so the K-frontier is visible, not hidden.
     A molecule claimed `assigned` by TWO families at once is tagged xf:Z:<stratum> (from
     `<prefix>.xfam_conflicts.tsv`, written under RUSTLE_XFAM_RECONCILE=report/abstain) and, for the
     cross_family_contradiction stratum, drawn near-black — it is not a copy call.
  <out>.copies.bed — one line per copy: its genomic extent (from its assigned reads), name = family|copy,
     itemRgb matching the read colour. A loci track to sit above the alignments.

Run:
  python bench/igv_tracks.py --assignments OUT.assignments.tsv --bam reads.bam --regions regions.txt --out OUT
Then in IGV load OUT.tagged.bam and OUT.copies.bed (+ the genome FASTA).
"""
import argparse
import subprocess
from collections import defaultdict

import pysam

# distinct, IGV-friendly colours per copy index; tied/ambiguous get greys (the unresolved classes).
PALETTE = [
    (228, 26, 28), (55, 126, 184), (77, 175, 74), (255, 127, 0), (152, 78, 163),
    (166, 86, 40), (247, 129, 191), (0, 158, 115), (153, 153, 0), (86, 180, 233),
]
TIED_RGB = (150, 150, 150)
AMB_RGB = (200, 200, 200)
# a molecule claimed `assigned` by two families at once: not a copy call, and not the K=0 wall either.
CONTESTED_RGB = (30, 30, 30)


def copy_color(idx):
    return PALETTE[idx % len(PALETTE)]


def load_assignments(path):
    """(read_name, family_id) -> (copy_idx, status), plus read_name -> [(family, copy, status), ...].

    CORRECTION (2026-08-31, measured). This function used to key on `read_name` alone and keep the FIRST
    row, on the stated ground that "the call is identical, so first wins". That is false, and by a wide
    margin. On the 12-family `--families` batch (mec/run_psv.sh: 79,175 rows over 77,372 molecules,
    1,587 of them with >= 2 rows):

      * first-wins names a DIFFERENT (family, copy, status) than the largest-margin ASSIGNED row for
        484/1,587 = 30.50% of multi-row molecules;
      * 136/1,587 have a NON-assigned first row while another row of the same molecule IS assigned, so
        first-wins painted them grey when the run had in fact resolved them;
      * only 324/1,587 have all their assigned rows naming one (family, copy) — i.e. the premise held
        for 20.4% of the cases it was asserted over.

    A row is a (molecule, FAMILY) genotyping result, so the key is `(read_name, family_id)`. The caller
    resolves which family a given BAM RECORD belongs to by genomic containment in that family's copy
    intervals (`load_copy_spans`), which is the only thing that makes the tag a statement about the
    record rather than about whichever row happened to be first.
    """
    amap = {}
    by_read = defaultdict(list)
    with open(path) as fh:
        next(fh)
        for ln in fh:
            f = ln.rstrip("\n").split("\t")
            amap[(f[0], f[1])] = (int(f[2]), f[3])
            by_read[f[0]].append((f[1], int(f[2]), f[3]))
    return amap, dict(by_read)


def load_copy_spans(prefix):
    """family_id -> [(chrom, start, end), ...] from `<prefix>.quant.tsv`, for record -> family resolution."""
    import os
    path = f"{prefix}.quant.tsv"
    if not os.path.exists(path):
        return {}
    spans = defaultdict(list)
    with open(path) as fh:
        next(fh)
        for ln in fh:
            f = ln.rstrip("\n").split("\t")
            if len(f) >= 6 and f[4].isdigit():
                spans[f[0]].append((f[3], int(f[4]), int(f[5])))
    return dict(spans)


def load_xfam_conflicts(prefix):
    """read_name -> worst stratum, from `<prefix>.xfam_conflicts.tsv` (RUSTLE_XFAM_RECONCILE=report/abstain).

    Absent unless that flag was on; when present it names exactly the molecules whose rows disagree
    ACROSS families, so they can be shown as contested instead of silently taking one row's colour.
    """
    import os
    path = f"{prefix}.xfam_conflicts.tsv"
    if not os.path.exists(path):
        return {}
    rank = {"shared_locus": 1, "readthrough_span": 2, "cross_family_contradiction": 3}
    worst = {}
    with open(path) as fh:
        next(fh)
        for ln in fh:
            f = ln.rstrip("\n").split("\t")
            st = f[1]
            if rank.get(st, 0) > rank.get(worst.get(f[0], ""), 0):
                worst[f[0]] = st
    return worst


def family_of_record(rec, rows, spans):
    """Which of this molecule's (family, copy, status) rows describes THIS alignment record.

    A record belongs to the family whose copy intervals it overlaps. Exactly one match => that row;
    no match or several => None, and the caller must not pretend to know (it tags the record contested
    rather than borrowing another placement's call).
    """
    if len(rows) == 1:
        return rows[0]
    hits = []
    for fam, ci, status in rows:
        for chrom, st, en in spans.get(fam, ()):
            if rec.reference_name == chrom and rec.reference_start < en and (rec.reference_end or 0) > st:
                hits.append((fam, ci, status))
                break
    return hits[0] if len(hits) == 1 else None


def tag_for(family, idx, status):
    if status == "assigned":
        return f"{family}_c{idx}", copy_color(idx), idx + 1
    if status == "tied":
        return f"{family}_tied", TIED_RGB, None
    return f"{family}_amb", AMB_RGB, None


def parse_regions(path):
    regs = []
    with open(path) as fh:
        for ln in fh:
            tok = ln.split()
            if not tok:
                continue
            c, rng = tok[0].split(":")
            s, e = rng.split("-")
            regs.append((c, int(s), int(e)))
    return regs


def write_psv_vcf(out, ca_prefix, regions, contig_lens):
    """PSV VCF (copies as samples) from copy_assign's --dump-psv matrix: <ca_prefix>.psv_cols.tsv (col->genome
    position) + .psv_copies.tsv (each copy's allele per column). Load in IGV alongside the tagged BAM: each PSV
    is a variant row whose per-copy genotype (0=ref allele, N=alt) is the reference the reads are matched to."""
    import os
    colf, copf = f"{ca_prefix}.psv_cols.tsv", f"{ca_prefix}.psv_copies.tsv"
    if not (os.path.exists(colf) and os.path.exists(copf)):
        print(f"  (no {colf}/.psv_copies.tsv — re-run copy_assign with --dump-psv for the PSV VCF)")
        return None
    cols = defaultdict(dict)   # family -> {col: genome_pos}
    for ln in open(colf):
        f = ln.rstrip("\n").split("\t")
        if len(f) >= 3 and f[1].isdigit():
            cols[f[0]][int(f[1])] = int(f[2])
    cops = defaultdict(dict)   # family -> {copy_idx: allele_string}
    for ln in open(copf):
        f = ln.rstrip("\n").split("\t")
        if len(f) >= 4 and f[1].isdigit():
            cops[f[0]][int(f[1])] = f[3]
    def chrom_of(pos):
        for c, s, e in regions:
            if s <= pos <= e:
                return c
        return regions[0][0] if regions else "."
    samples = [f"{fam}_c{ci}" for fam in sorted(cols) for ci in sorted(cops.get(fam, {}))]
    sidx = {s: i for i, s in enumerate(samples)}
    rows = []
    for fam in sorted(cols):
        cps = cops.get(fam, {})
        for col, pos in sorted(cols[fam].items()):
            if pos < 0:
                continue
            bases = {ci: cps[ci][col] for ci in cps if col < len(cps[ci]) and cps[ci][col] != '.'}
            if len(set(bases.values())) < 2:            # not decisive (all copies agree)
                continue
            ref = bases.get(0) or sorted(bases.values())[0]
            alts = [b for b in sorted(set(bases.values())) if b != ref]
            aidx = {b: i + 1 for i, b in enumerate(alts)}
            gts = ["."] * len(samples)
            for ci, b in bases.items():
                gts[sidx[f"{fam}_c{ci}"]] = "0" if b == ref else str(aidx[b])
            rows.append((chrom_of(pos), pos + 1, f"{fam}_col{col}", ref, ",".join(alts), fam, gts))
    vcf = f"{out}.psv.vcf"
    with open(vcf, "w") as fh:
        fh.write("##fileformat=VCFv4.2\n")
        fh.write('##INFO=<ID=FAM,Number=1,Type=String,Description="copy-family id">\n')
        fh.write('##FORMAT=<ID=GT,Number=1,Type=String,Description="copy allele: 0=REF, N=Nth ALT, .=gapped">\n')
        for c, l in contig_lens:
            fh.write(f"##contig=<ID={c},length={l}>\n")
        fh.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t" + "\t".join(samples) + "\n")
        for chrom, pos, rid, ref, alt, fam, gts in sorted(rows):
            fh.write(f"{chrom}\t{pos}\t{rid}\t{ref}\t{alt}\t.\tPASS\tFAM={fam}\tGT\t" + "\t".join(gts) + "\n")
    return vcf, len(rows), len(samples)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--assignments", required=True)
    ap.add_argument("--bam", required=True)
    ap.add_argument("--regions", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--samtools", default="samtools")
    a = ap.parse_args()

    import re as _re
    prefix = _re.sub(r"\.assignments\.tsv$", "", a.assignments)
    amap, by_read = load_assignments(a.assignments)
    spans = load_copy_spans(prefix)
    contested = load_xfam_conflicts(prefix)
    if contested:
        print(f"  {len(contested)} molecule(s) contested across families (from {prefix}.xfam_conflicts.tsv)")
    regions = parse_regions(a.regions)
    inbam = pysam.AlignmentFile(a.bam, "rb")
    bam_contigs = list(zip(inbam.references, inbam.lengths))  # captured before close (for the VCF header)
    tmp = f"{a.out}.tagged.unsorted.bam"
    out = pysam.AlignmentFile(tmp, "wb", template=inbam)

    extent = defaultdict(lambda: [None, 1 << 62, 0, None, None])  # key (fam,idx) -> [chrom, start, end, status, idx]
    seen = set()
    n_written = n_tagged = 0
    for (chrom, s, e) in regions:
        for rec in inbam.fetch(chrom, max(0, s), e):
            if rec.is_unmapped:
                continue
            key = (rec.query_name, rec.reference_start, rec.flag)
            if key in seen:
                continue
            seen.add(key)
            rows = by_read.get(rec.query_name)
            info = family_of_record(rec, rows, spans) if rows else None
            if rows and info is None:
                # The record falls in copy intervals of SEVERAL families (or of none). There is no
                # per-record call to draw, and borrowing another placement's call is exactly the
                # first-wins defect. Mark it contested and move on.
                strat = contested.get(rec.query_name, "unresolved_family")
                rgb = CONTESTED_RGB if strat == "cross_family_contradiction" else AMB_RGB
                rec.set_tag("xf", strat, "Z")
                rec.set_tag("cp", f"xfam_{strat}", "Z")
                rec.set_tag("YC", f"{rgb[0]},{rgb[1]},{rgb[2]}", "Z")
                n_tagged += 1
            if info:
                fam, idx, status = info
                cp, (r, g, b), hp = tag_for(fam, idx, status)
                strat = contested.get(rec.query_name)
                if strat:
                    # The molecule is claimed by >= 2 families. Say so on the record instead of letting one
                    # family's colour stand for a call the run did not make.
                    rec.set_tag("xf", strat, "Z")
                    if strat == "cross_family_contradiction":
                        cp, (r, g, b), hp = f"{fam}_xfam", CONTESTED_RGB, None
                rec.set_tag("cp", cp, "Z")
                rec.set_tag("YC", f"{r},{g},{b}", "Z")
                if hp is not None:
                    rec.set_tag("HP", hp, "i")
                n_tagged += 1
                if status == "assigned" and not contested.get(rec.query_name) == "cross_family_contradiction" \
                        and not (rec.is_secondary or rec.is_supplementary):
                    k = (fam, idx)
                    ex = extent[k]
                    ex[0] = rec.reference_name
                    ex[1] = min(ex[1], rec.reference_start)
                    ex[2] = max(ex[2], rec.reference_end or rec.reference_start)
                    ex[3], ex[4] = status, idx
            out.write(rec)
            n_written += 1
    out.close()
    inbam.close()

    sorted_bam = f"{a.out}.tagged.bam"
    subprocess.run([a.samtools, "sort", "-o", sorted_bam, tmp], check=True)
    subprocess.run([a.samtools, "index", sorted_bam], check=True)
    subprocess.run(["rm", "-f", tmp])

    # copies BED (genomic extent of each copy's assigned reads), itemRgb matching the read colour
    with open(f"{a.out}.copies.bed", "w") as fh:
        fh.write('track name="copies" itemRgb="On"\n')
        for (fam, idx), (chrom, st, en, status, ci) in sorted(extent.items()):
            if chrom is None:
                continue
            r, g, b = copy_color(idx)
            fh.write(f"{chrom}\t{st}\t{en}\t{fam}_c{idx}\t0\t+\t{st}\t{en}\t{r},{g},{b}\n")

    print(f"wrote {sorted_bam} (+ .bai): {n_written} alignments, {n_tagged} tagged by copy")
    print(f"wrote {a.out}.copies.bed: {sum(1 for v in extent.values() if v[0])} copies")

    # PSV VCF (copies as samples) from the --dump-psv matrix, if present
    ca_prefix = prefix
    reg_chroms = {c for c, _, _ in regions}
    contig_lens = [(c, l) for c, l in bam_contigs if c in reg_chroms]
    v = write_psv_vcf(a.out, ca_prefix, regions, contig_lens)
    if v:
        print(f"wrote {v[0]}: {v[1]} PSV sites x {v[2]} copies (IGV variant track — per-copy alleles)")

    print(f"IGV: load {sorted_bam} + {a.out}.copies.bed + {a.out}.psv.vcf + the genome FASTA. Reads auto-colour by")
    print(f"     copy (YC tag); right-click -> 'Group by tag' -> cp to stack copies; the VCF marks each PSV.")


if __name__ == "__main__":
    main()
