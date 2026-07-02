#!/usr/bin/env python3
"""
diploid_cn_oracle.py

Use the PHASED DIPLOID mGorGor1 assembly (maternal + paternal haplotypes, built
from DNA WGS, fully independent of the IsoSeq RNA reads) as a copy-number oracle
to give a NON-CIRCULAR corroboration of the RNA-derived multi-copy family method.

For each co-located family in validated_families.tsv:
  1. Extract the genomic sequence of ALL its reference copy loci from GGO.fasta.
  2. Align every copy locus to EACH haplotype index (minimap2 -cx asm20 -N 50
     -p 0.1 --secondary=yes).
  3. Count DISTINCT haplotype loci that carry a FULL-LENGTH copy of the family's
     gene unit, at identity >= 0.90. A hit only counts as a copy if the aligned
     QUERY coverage of that target region reaches >= COV_FRAC of the query length
     (a length-PROPORTIONAL criterion, NOT an absolute 1 kb floor) -- this stops a
     shared sub-gene domain / repeat fragment from counting as a full gene copy.
     Alignments of one query to the SAME target that OVERLAP in query space are
     kept as SEPARATE copies (tandem array units); query-disjoint collinear blocks
     of one query are stitched into one copy (a single alignment split by an
     indel). That distinct-copy count = haploid copy number (hap_CN_mat/pat).
     Distinct sub-COV_FRAC target regions are reported separately as `frag_*`.

copy-vs-allele resolution (sex-aware; mGorGor1 is MALE = maternal chrX, paternal
chrY, so an X-linked gene reads hap_CN_pat==0 and a Y-linked gene hap_CN_mat==0):
  autosome, hap_CN>=2 in BOTH haplotypes  -> multi_copy  (real paralogy; if the
                                             two haplotype counts differ it is ALSO
                                             flagged cnv_autosomal=1)
  autosome, hap_CN==1 in BOTH             -> single_locus_allele (RNA variants =
                                             allele/het, NOT copies)
  autosome, counts differ w/ one < 2      -> cnv          (autosomal CNV/deletion)
  sex-linked: present haplotype count h   -> multi_copy (h>=2) / single_locus_allele
                                             (h==1) / absent (h==0)  [male hemizygous]

CAVEAT: copy number / paralogy is largely SPECIES-FIXED, so the assembly is a
valid copy-number oracle regardless of whether the IsoSeq donor is mGorGor1.
Individual-exact het/allele calls for a specific variant require RNA==mGorGor1
(unknown). Lead with the robust copy-number use.
"""
import argparse
import os
import subprocess
import sys
from collections import defaultdict
from statistics import median

import pysam

MINIMAP2 = "/home/juanfra/miniforge3/bin/minimap2"
SCRATCH = "/home/juanfra/winloci_scratch"
GGO = os.path.join(SCRATCH, "GGO.fasta")
MAT_MMI = os.path.join(SCRATCH, "mGorGor1.mat.asm20.rebuild.mmi")
PAT_MMI = os.path.join(SCRATCH, "mGorGor1.pat.asm20.mmi")

VALIDATED = os.path.join(SCRATCH, "validated_families.tsv")
A1 = "bench/a1_read_consensus_o1.tsv"

MIN_IDENT = 0.90       # >= 90% identity  (Bailey/WGAC segdup)
PREFLOOR = 300         # discard tiny (<300 bp) alignment blocks while parsing
COV_FRAC = 0.5         # a copy must aggregate >= 50% of the QUERY length
Q_OVERLAP = 0.20       # blocks re-covering > 20% of query start a new copy (tandem)

# chrX / chrY contigs of the GGO reference (GCF_029281585.2). The maternal
# haplotype carries chrX, the paternal carries chrY, so X-linked genes read as
# hap_CN_pat == 0 (present on X, absent on Y) and Y-linked as hap_CN_mat == 0.
# That is expected male hemizygosity, NOT an inter-haplotype CNV.
CHRX = "NC_073247.2"
CHRY = "NC_073248.2"


def read_validated(path):
    """family -> list of (locus, chrom, start, end, nreads). Preserve file order."""
    fam = defaultdict(list)
    loc2fam = {}
    with open(path) as fh:
        fh.readline()
        for line in fh:
            f = line.rstrip("\n").split("\t")
            family = f[0]
            locus, chrom = f[1], f[2]
            start, end, nreads = int(f[3]), int(f[4]), int(f[5])
            fam[family].append((locus, chrom, start, end, nreads))
            loc2fam[locus] = family
    return fam, loc2fam


def fam_sexlink(loci):
    """X / Y / mixed_sex / auto from the reference chroms of a family's loci."""
    chroms = {c for (_l, c, _s, _e, _n) in loci}
    if chroms == {CHRX}:
        return "X"
    if chroms == {CHRY}:
        return "Y"
    if CHRX in chroms or CHRY in chroms:
        return "mixed_sex"
    return "auto"


def read_a1(path):
    """family -> dict(gene, n_loci_a1, chi_H, K_read, sedef_partners)."""
    a1 = {}
    if not os.path.exists(path):
        return a1
    with open(path) as fh:
        cols = fh.readline().rstrip("\n").split("\t")
        idx = {c: i for i, c in enumerate(cols)}
        for line in fh:
            f = line.rstrip("\n").split("\t")
            fam = f[idx["family"]]
            a1[fam] = dict(
                gene=f[idx["gene"]],
                n_loci_a1=f[idx["n_loci"]],
                chi_H=f[idx["chi_H"]],
                K_read=f[idx["K_read"]],
                sedef_partners=f[idx["sedef_partner_regions"]],
            )
    return a1


def build_query_fasta(fam, out_fa):
    """Extract every locus sequence from GGO.fasta -> multi-FASTA. Header=locus."""
    ref = pysam.FastaFile(GGO)
    n = 0
    with open(out_fa, "w") as out:
        for family in sorted(fam, key=lambda x: int(x)):
            for (locus, chrom, start, end, _nr) in fam[family]:
                seq = ref.fetch(chrom, start, end)
                out.write(">%s\n" % locus)
                for i in range(0, len(seq), 80):
                    out.write(seq[i:i + 80] + "\n")
                n += 1
    ref.close()
    return n


def run_minimap2(query_fa, mmi, paf_out, threads):
    cmd = [
        MINIMAP2, "-cx", "asm20", "-N", "50", "-p", "0.1",
        "--secondary=yes", "-t", str(threads),
        mmi, query_fa,
    ]
    with open(paf_out, "w") as out:
        subprocess.run(cmd, stdout=out, stderr=subprocess.DEVNULL, check=True)


def parse_paf_ident(fields):
    """identity = 1 - de:f if present else nmatch/blen."""
    nmatch = int(fields[9])
    blen = int(fields[10])
    de = None
    for tag in fields[12:]:
        if tag.startswith("de:f:"):
            de = float(tag[5:])
            break
    ident = (1.0 - de) if de is not None else (nmatch / blen if blen else 0.0)
    return ident, blen


def read_paf_hits(paf_path, loc2fam):
    """
    family -> list of per-alignment records:
        (qname, qlen, qs, qe, strand, target, ts, te, blen)
    filtered to identity >= MIN_IDENT and block >= PREFLOOR.
    """
    hits = defaultdict(list)
    with open(paf_path) as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) < 12:
                continue
            fam = loc2fam.get(f[0])
            if fam is None:
                continue
            ident, blen = parse_paf_ident(f)
            if ident < MIN_IDENT or blen < PREFLOOR:
                continue
            hits[fam].append((
                f[0], int(f[1]), int(f[2]), int(f[3]), f[4],
                f[5], int(f[7]), int(f[8]), blen))
    return hits


def _ov(a0, a1, b0, b1):
    return max(0, min(a1, b1) - max(a0, b0))


def _segment_query(blocks):
    """
    blocks for ONE (query, target, strand), each (qs, qe, ts, te). Segment into
    copies: a block joins the current copy if it does NOT re-cover the copy's
    query span by > Q_OVERLAP (else it is a distinct tandem unit) and its target
    is within an intron-scale gap of the copy (else a separate locus). Returns a
    list of copies, each = dict(qcov, t_lo, t_hi).
    """
    if not blocks:
        return []
    Lq = max(qe for (_qs, qe, _ts, _te) in blocks) or 1
    intron_cap = max(5000, int(0.5 * Lq))
    blocks = sorted(blocks, key=lambda b: b[2])   # by target start
    copies = []
    for (qs, qe, ts, te) in blocks:
        placed = False
        for cp in copies:
            if ts - cp["t_hi"] > intron_cap:       # too far on target (sorted -> only grows)
                continue
            qov = _ov(qs, qe, cp["q_lo"], cp["q_hi"])
            if qov > Q_OVERLAP * (qe - qs):        # re-covers query -> tandem/paralog unit
                continue
            cp["qivs"].append((qs, qe))
            cp["q_lo"] = min(cp["q_lo"], qs)
            cp["q_hi"] = max(cp["q_hi"], qe)
            cp["t_hi"] = max(cp["t_hi"], te)
            cp["t_lo"] = min(cp["t_lo"], ts)
            placed = True
            break
        if not placed:
            copies.append(dict(qivs=[(qs, qe)], q_lo=qs, q_hi=qe,
                               t_lo=ts, t_hi=te))
    out = []
    for cp in copies:
        ivs = sorted(cp["qivs"])
        cov, ce = 0, -1
        for (s, e) in ivs:
            s = max(s, ce)
            if e > s:
                cov += e - s
                ce = e
        out.append(dict(qcov=cov, t_lo=cp["t_lo"], t_hi=cp["t_hi"]))
    return out


def _merge_intervals(ivs):
    """count distinct target regions after gap-0 overlap merge; ivs=[(t,lo,hi)]."""
    by_t = defaultdict(list)
    for (t, lo, hi) in ivs:
        by_t[t].append((lo, hi))
    n = 0
    for t, xs in by_t.items():
        xs.sort()
        clo, chi = xs[0]
        for (lo, hi) in xs[1:]:
            if lo <= chi:                # strict overlap only (gap 0)
                chi = max(chi, hi)
            else:
                n += 1
                clo, chi = lo, hi
        n += 1
    return n


def count_cn(hits, l_fam, cov_frac=COV_FRAC):
    """
    hits = family alignment records. Returns (hap_CN, n_frag):
      hap_CN = distinct target regions where SOME query aligns a FULL copy of the
               family gene unit: aggregate collinear query coverage of that target
               region >= cov_frac * L (L = median reference locus length of the
               family, a length-PROPORTIONAL gene-unit criterion). This rejects a
               shared sub-gene domain that would look "full" only relative to a
               short/truncated reference locus (the fam76 heterogeneous-size leak).
      n_frag = distinct target regions covered ONLY by sub-cov_frac fragments
               (shared domain / repeat) and not overlapping any full copy.
    """
    thr = cov_frac * l_fam
    by_qt = defaultdict(list)
    for (q, _ql, qs, qe, strand, t, ts, te, _bl) in hits:
        by_qt[(q, t, strand)].append((qs, qe, ts, te))

    full, frag = [], []
    for (q, t, strand), blocks in by_qt.items():
        for cp in _segment_query(blocks):
            iv = (t, cp["t_lo"], cp["t_hi"])
            if cp["qcov"] >= thr:
                full.append(iv)
            else:
                frag.append(iv)

    hap_cn = _merge_intervals(full)
    # a fragment region only counts if it does not overlap any full-copy region
    full_by_t = defaultdict(list)
    for (t, lo, hi) in full:
        full_by_t[t].append((lo, hi))
    frag_only = []
    for (t, lo, hi) in frag:
        if not any(_ov(lo, hi, flo, fhi) > 0 for (flo, fhi) in full_by_t.get(t, [])):
            frag_only.append((t, lo, hi))
    return hap_cn, _merge_intervals(frag_only)


def classify(cm, cp, sex):
    """sex-aware primary class + autosomal-CNV flag (both>=2 differing = also cnv)."""
    if sex in ("X", "Y"):
        h = cm if sex == "X" else cp
        cls = ("multi_copy" if h >= 2 else
               "single_locus_allele" if h == 1 else "absent")
        return cls, 0
    # autosome / mixed
    if cm >= 2 and cp >= 2:
        return "multi_copy", int(cm != cp)
    if cm == 1 and cp == 1:
        return "single_locus_allele", 0
    if cm != cp:
        return "cnv", 0
    return "absent", 0            # 0/0


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default="bench/diploid_cn_oracle.tsv")
    ap.add_argument("--pafdir", default=os.path.join(SCRATCH, "cn_oracle_paf"))
    ap.add_argument("--threads", type=int, default=8)
    ap.add_argument("--reuse-paf", action="store_true",
                    help="reuse existing PAFs / query fasta if present")
    args = ap.parse_args()

    os.makedirs(args.pafdir, exist_ok=True)
    query_fa = os.path.join(args.pafdir, "query_loci.fa")
    mat_paf = os.path.join(args.pafdir, "mat.paf")
    pat_paf = os.path.join(args.pafdir, "pat.paf")

    fam, loc2fam = read_validated(VALIDATED)
    a1 = read_a1(A1)

    if not (args.reuse_paf and os.path.exists(query_fa)):
        nq = build_query_fasta(fam, query_fa)
        sys.stderr.write("extracted %d query loci -> %s\n" % (nq, query_fa))
    if not (args.reuse_paf and os.path.exists(mat_paf)):
        sys.stderr.write("aligning to maternal ...\n")
        run_minimap2(query_fa, MAT_MMI, mat_paf, args.threads)
    if not (args.reuse_paf and os.path.exists(pat_paf)):
        sys.stderr.write("aligning to paternal ...\n")
        run_minimap2(query_fa, PAT_MMI, pat_paf, args.threads)

    mat_hits = read_paf_hits(mat_paf, loc2fam)
    pat_hits = read_paf_hits(pat_paf, loc2fam)

    rows = []
    for family in sorted(fam, key=lambda x: int(x)):
        loci = fam[family]
        sex = fam_sexlink(loci)
        l_fam = int(median([e - s for (_l, _c, s, e, _n) in loci]))
        cm, fm = count_cn(mat_hits.get(family, []), l_fam)
        cp, fp = count_cn(pat_hits.get(family, []), l_fam)
        cls, cnv_auto = classify(cm, cp, sex)
        asm = (cm if sex == "X" else cp if sex == "Y" else max(cm, cp))
        meta = a1.get(family, {})
        rows.append((
            family, meta.get("gene", "NA"), len(loci),
            meta.get("K_read", "NA"), meta.get("chi_H", "NA"),
            meta.get("sedef_partners", "NA"),
            sex, cm, cp, fm, fp, asm, cls, cnv_auto,
        ))

    header = ["family", "gene", "n_loci_ref", "K_read", "chi_H",
              "sedef_partners", "sex", "hap_CN_mat", "hap_CN_pat",
              "frag_mat", "frag_pat", "asm_hapCN", "class", "cnv_autosomal"]
    with open(args.out, "w") as out:
        out.write("\t".join(header) + "\n")
        for r in rows:
            out.write("\t".join(str(x) for x in r) + "\n")
    sys.stderr.write("wrote %d families -> %s\n" % (len(rows), args.out))


if __name__ == "__main__":
    main()
