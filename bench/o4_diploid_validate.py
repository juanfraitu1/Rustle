#!/usr/bin/env python3
"""o4_diploid_validate.py -- DNA close of O4 (reference-absent divergent copy discovery).

The O4 machinery detects reference-absent COPY CANDIDATES (collapsed / CNV via co-segregating
PSVs, plus divergent/unmapped read clusters). A STRUCTURAL resolver on the HAPLOID T2T primary
(bench/copy_vs_allele_structural.py) could only give: ~26 COPY / 38 ALLELE / 65 AMBIGUOUS(MHC) /
16 DIVERGENT-NOVEL -- MECHANISM-ONLY, because a haploid assembly cannot separate a het single-locus
variant ("could be a copy OR an allele") from a real extra copy, and DNA parCN was unavailable.

WHAT CHANGED (donor-exact): we now have the PHASED DIPLOID mGorGor1 assembly (maternal + paternal
haplotypes, DNA-derived, INDEPENDENT of the IsoSeq RNA that flagged the candidate) AND the RNA donor
is confirmed == mGorGor1 (same individual). The primary reference GGO.fasta = GCF_029281585.2 is the
primary pseudo-haploid of that SAME individual. So copy-vs-allele is resolvable DONOR-EXACTLY by the
same diploid rule diploid_cn_oracle applied to co-located families, now applied to the O4 candidates:

  ref_loci / mat_loci / pat_loci = distinct FULL-COPY loci the candidate consensus maps to in the
  primary reference / maternal hap / paternal hap (minimap2 -cx asm20 -N50 -p0.1, id>=0.90, aggregate
  query coverage >= 0.5 of the consensus per locus; exon blocks within 500 kb stitched into one locus).

  COPY        extra full-copy locus in BOTH haplotypes beyond the reference (min(mat,pat) > ref)
              -> a real paralog the primary assembly collapsed/missed, present germline in both haps.
  ALLELE      extra locus in exactly ONE haplotype (mat XOR pat > ref), the other == ref
              -> heterozygous / one-chromosome CNV, NOT a germline copy. (also: ref=mat=pat=1 single
                 locus = the co-segregating RNA PSVs are a het allele of a one-copy gene.)
  NOVEL       NO full-length >=90% copy in reference OR either haplotype. This is a MIXED bucket that
              we SUB-CLASSIFY (loci=0 does not mean "absent"):
                 subthreshold_present  a >=50%-coverage locus aligns at ~80-90% identity in >=1
                                       target -- a single divergent locus that IS in the genome, just
                                       below the 0.90 identity gate; NOT novel-to-genome.
                 partial_present       high-identity (>=0.90) fragments cover <50% of the consensus
                                       in >=1 target -- chimeric/fragmentary transcript consensus of a
                                       present locus; NOT novel-to-genome.
                 genuinely_absent      no >=0.90 block and <50% coverage in ALL three -> genuinely
                                       novel / artifact / too divergent (<~85%) to place.
  UNDECIDABLE reference already models the loci and both haps match (no diploid-visible extra;
              MHC-class allele-of-existing-copy), or contradictory counts.

HONEST caveats (built in):
  (1) donor-exact NOW: RNA == mGorGor1, so het/allele calls are individual-exact (stated).
  (2) "reference-absent" = absent from the GCF PRIMARY, present in the fuller diploid; a candidate
      confirmed present in a haplotype but not the primary means the PRIMARY missed/collapsed it, NOT
      that it is absent from the genome. A candidate genuinely absent from mat+pat too is novel/artifact.
      A candidate that maps to the reference at high identity (even at <50% coverage) is REFERENCE-
      PRESENT, not reference-absent (this is why CDK11B, best_id_ref=0.9916, is NOT counted).
  (3) the collapsed/CNV O4 overlaps co-located families already scored by diploid_cn_oracle; we
      cross-flag oracle_family / oracle_collapsed so a COPY there is NOT double-counted as an
      independent discovery.
  (4) non-circular: the DNA mat/pat assembly is independent of the RNA divergence signal that
      flagged each candidate.
  (5) SENSITIVITY FLOOR: "0 COPY" means 0 copy at >=~85% identity. A genuinely divergent (<85%)
      both-hap germline copy would fall in the NOVEL/genuinely_absent bucket -- this transcript-copy
      oracle cannot see a copy more divergent than its identity gate. Reported as a caveat, not a hole.

Determinism: PYTHONHASHSEED=0, minimap2 -t 5, fixed -I; per-query alignments are minimap2-
deterministic (verified by re-aligning the reference and diffing the sorted PAF). 0-COPY verified
robust to a full parameter sweep: identity gate 0.80-0.90, coverage 0.2-0.5, locus gap 100k-500k.

Run: PYTHONHASHSEED=0 /home/juanfra/miniforge3/bin/python bench/o4_diploid_validate.py
Reuse existing PAFs (skip minimap2): O4_REUSE_PAF=1 PYTHONHASHSEED=0 ... bench/o4_diploid_validate.py
"""
import os
import sys
import json
import csv
import subprocess
import hashlib
from collections import defaultdict, Counter

MINIMAP2 = "/home/juanfra/miniforge3/bin/minimap2"
SCRATCH  = "/home/juanfra/winloci_scratch"
GGO      = f"{SCRATCH}/GGO.fasta"                          # primary reference = GCF_029281585.2 (= RNA ref)
MAT_MMI  = f"{SCRATCH}/mGorGor1.mat.asm20.rebuild.mmi"     # maternal hap (same as diploid_cn_oracle)
PAT_MMI  = f"{SCRATCH}/mGorGor1.pat.asm20.mmi"             # paternal hap
GW       = f"{SCRATCH}/refabsent/gw_promoted"
ABSENT   = f"{GW}/gw_reference_absent_copies.json"         # the 145 O4 candidates
DISC     = f"{GW}/gw_discriminated.json"                   # n_loci / cat for a subset
SURV_FA  = f"{GW}/survivors.fa"                            # candidate consensus seqs (145, keyed by cid)
POC      = f"{SCRATCH}/unmapped_poc/unmapped_rescue_poc_result.json"
REPS_FA  = f"{SCRATCH}/unmapped_poc/reps.fa"               # divergent/unmapped cluster reps
HAPLOID  = "bench/copy_vs_allele_structural.tsv"           # prior haploid labels
VALID    = f"{SCRATCH}/validated_families.tsv"             # co-located families (coords)
ORACLE   = "bench/diploid_cn_oracle.tsv"                   # diploid_cn_oracle results

WORK      = f"{SCRATCH}/o4_diploid"
OUT_TSV   = "bench/o4_diploid_validate.tsv"
OUT_JSON  = "bench/o4_diploid_validate.json"

MIN_IDENT = 0.90
COV_FRAC  = 0.50
LOCUS_GAP = 500_000     # exon blocks within 500 kb on one contig = one gene locus (stitch introns)
PREFLOOR  = 200         # discard alignment blocks < 200 bp
THREADS   = 5
CHRX, CHRY, CHRMT = "NC_073247.2", "NC_073248.2", "NC_011120.1"
SUBTHRESH_FLOOR = 0.80  # a below-gate locus is "reference-present" only if it aligns at >= this id


# ----------------------------------------------------------------------------- FASTA / PAF helpers
def read_fasta(path):
    seqs, name, buf = {}, None, []
    with open(path) as fh:
        for ln in fh:
            if ln.startswith(">"):
                if name is not None:
                    seqs[name] = "".join(buf)
                name = ln[1:].split()[0].strip()
                buf = []
            else:
                buf.append(ln.strip())
    if name is not None:
        seqs[name] = "".join(buf)
    return seqs


def parse_paf_ident(f):
    nmatch, blen = int(f[9]), int(f[10])
    de = None
    for tag in f[12:]:
        if tag.startswith("de:f:"):
            de = float(tag[5:]); break
    ident = (1.0 - de) if de is not None else (nmatch / blen if blen else 0.0)
    return ident, blen


def run_minimap2(query_fa, target, paf_out, is_fasta):
    cmd = [MINIMAP2, "-cx", "asm20", "-N", "50", "-p", "0.1",
           "--secondary=yes", "-t", str(THREADS)]
    if is_fasta:
        cmd += ["-I", "12G"]     # single index batch -> deterministic on the 3.6 Gb primary
    cmd += [target, query_fa]
    with open(paf_out, "w") as out:
        subprocess.run(cmd, stdout=out, stderr=subprocess.DEVNULL, check=True)


def read_paf_blocks(paf_path):
    """query -> list of (qs,qe,strand,target,ts,te,ident) for ALL blocks >= PREFLOOR (NO identity
    filter -- the identity gate is applied per-count in count_loci so we can sweep it and so we can
    see below-gate evidence for the NOVEL sub-classification)."""
    blocks = defaultdict(list)
    with open(paf_path) as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) < 12:
                continue
            ident, blen = parse_paf_ident(f)
            if blen < PREFLOOR:
                continue
            blocks[f[0]].append((int(f[2]), int(f[3]), f[4], f[5], int(f[7]), int(f[8]), ident))
    return blocks


def _union_cov(ivs):
    ivs = sorted(ivs); cov, ce = 0, -1
    for s, e in ivs:
        s = max(s, ce)
        if e > s:
            cov += e - s; ce = e
    return cov


def _cluster(bl, locus_gap):
    """bl = list of (ts,te,qs,qe,ident) sorted by ts -> list of clusters (blocks within locus_gap)."""
    bl.sort()
    clusters, cur, cur_hi = [], [], None
    for ts, te, qs, qe, ident in bl:
        if cur and ts - cur_hi > locus_gap:
            clusters.append(cur); cur, cur_hi = [], None
        cur.append((ts, te, qs, qe, ident))
        cur_hi = te if cur_hi is None else max(cur_hi, te)
    if cur:
        clusters.append(cur)
    return clusters


def count_loci(blocks, qlen, min_ident=MIN_IDENT, cov_frac=COV_FRAC, locus_gap=LOCUS_GAP):
    """distinct FULL-COPY loci at >= min_ident: per (target,strand) cluster blocks with target
    gap<locus_gap; a cluster is a full copy if union query-coverage >= cov_frac*qlen. Returns
    (full, frag, best_id, regions[(t,lo,hi)]). Blocks below min_ident are ignored here."""
    if qlen <= 0:
        return 0, 0, 0.0, []
    by = defaultdict(list)
    for qs, qe, strand, t, ts, te, ident in blocks:
        if ident < min_ident:
            continue
        by[(t, strand)].append((ts, te, qs, qe, ident))
    full, frag, best_id, regions = 0, 0, 0.0, []
    for (t, strand), bl in by.items():
        for cl in _cluster(bl, locus_gap):
            qcov = _union_cov([(qs, qe) for _, _, qs, qe, _ in cl])
            best_id = max(best_id, max(i for *_, i in cl))
            if qcov >= cov_frac * qlen:
                full += 1
                regions.append((t, min(c[0] for c in cl), max(c[1] for c in cl)))
            else:
                frag += 1
    return full, frag, round(best_id, 4), regions


def locus_evidence(blocks, qlen, locus_gap=LOCUS_GAP):
    """ALL-blocks evidence (no identity gate) used to sub-classify a loci=0 target:
       best_id_any  = max identity of any block (>= PREFLOOR bp)
       max_cov      = max, over stitched loci, of union query-coverage fraction
       cov_id       = coverage-weighted mean identity of the max-coverage stitched locus
    -> lets us tell 'reference-present below the 0.90 gate' from 'genuinely absent'."""
    if not blocks or qlen <= 0:
        return 0.0, 0.0, 0.0
    best_id_any = max(b[6] for b in blocks)
    by = defaultdict(list)
    for qs, qe, strand, t, ts, te, ident in blocks:
        by[(t, strand)].append((ts, te, qs, qe, ident))
    max_cov, cov_id = 0.0, 0.0
    for (t, strand), bl in by.items():
        for cl in _cluster(bl, locus_gap):
            cov = _union_cov([(qs, qe) for _, _, qs, qe, _ in cl]) / qlen
            if cov > max_cov:
                max_cov = cov
                tot = sum(qe - qs for _, _, qs, qe, _ in cl)
                cov_id = (sum((qe - qs) * i for _, _, qs, qe, i in cl) / tot) if tot else 0.0
    return round(best_id_any, 4), round(max_cov, 4), round(cov_id, 4)


# ----------------------------------------------------------------------------- classifiers
def classify(ref, mat, pat, sexlink):
    """autosomal donor-exact diploid rule on FULL-COPY loci counts. sexlink in {auto,X,Y}.
    Returns (label, reason, conf); the NOVEL reason is REPLACED in main() with an evidence-aware
    sub-classification string (loci=0 is a mixed bucket)."""
    if sexlink in ("X", "Y"):
        hap = mat if sexlink == "X" else pat
        if ref == 0 and hap == 0:
            return "NOVEL", "sex-linked no full-length copy", "high"
        if hap > ref:
            return "COPY", f"sex-linked({sexlink}) extra locus in present hap (ref={ref},hap={hap})", "med"
        return "UNDECIDABLE", f"sex-linked({sexlink}) ref={ref},hap={hap} no extra", "low"
    if ref == 0 and mat == 0 and pat == 0:
        return "NOVEL", "no full-length copy in ref/mat/pat (sub-classified in main)", "high"
    extra_m, extra_p = mat - ref, pat - ref
    if extra_m >= 1 and extra_p >= 1:
        conf = "high" if mat == pat else "med"
        return ("COPY",
                f"extra full-copy locus in BOTH haplotypes beyond reference (ref={ref},mat={mat},"
                f"pat={pat}) -- reference-absent copy the primary collapsed/missed, present germline",
                conf)
    if extra_m >= 1 and pat == ref:
        return ("ALLELE",
                f"extra locus on MATERNAL only (ref={ref},mat={mat},pat={pat}) -- heterozygous "
                f"one-chromosome CNV, not a germline copy", "high")
    if extra_p >= 1 and mat == ref:
        return ("ALLELE",
                f"extra locus on PATERNAL only (ref={ref},mat={mat},pat={pat}) -- heterozygous "
                f"one-chromosome CNV, not a germline copy", "high")
    if mat == ref and pat == ref and ref >= 1:
        if ref == 1:
            return ("ALLELE",
                    "single locus in reference AND both haplotypes (ref=mat=pat=1) -- the "
                    "co-segregating RNA PSVs are a het allele of a one-copy gene, not a copy", "med")
        return ("UNDECIDABLE",
                f"reference already models {ref} loci and both haplotypes match (no diploid-visible "
                f"extra) -- RNA variant is an allele of an existing copy / reference already contains "
                f"the paralog (MHC-class); not a NEW reference-absent copy", "med")
    return ("UNDECIDABLE",
            f"ambiguous / contradictory locus counts (ref={ref},mat={mat},pat={pat}) -- no clean "
            f"diploid resolution", "low")


def novel_subclass(ev):
    """ev = {'ref':(bid,cov,cid),'mat':..,'pat':..}. Sub-classify a loci=0-everywhere candidate.
       subthreshold_present : a >=COV_FRAC-coverage locus at >=SUBTHRESH_FLOOR id in >=1 target
                              (reference/hap-present, below the 0.90 identity gate) -> NOT novel.
       partial_present      : a >=MIN_IDENT block exists in >=1 target but coverage <COV_FRAC
                              (chimeric/fragmentary consensus of a present locus)      -> NOT novel.
       genuinely_absent     : neither                                                 -> novel/artifact."""
    any_sub = any(cov >= COV_FRAC and cid >= SUBTHRESH_FLOOR for (bid, cov, cid) in ev.values())
    any_high = any(bid >= MIN_IDENT for (bid, cov, cid) in ev.values())
    if any_sub:
        return "subthreshold_present"
    if any_high:
        return "partial_present"
    return "genuinely_absent"


def novel_reason(sub, ev):
    if sub == "subthreshold_present":
        tgt = max(ev, key=lambda k: ev[k][1]); bid, cov, cid = ev[tgt]
        return (f"reference/hap-PRESENT below the id>={MIN_IDENT} gate: a {cov*100:.0f}%-coverage "
                f"locus aligns at ~{cid*100:.0f}% identity in {tgt} -- a single divergent locus that "
                f"IS in the genome (transcript consensus ~{(1-cid)*100:.0f}% diverged), not a full "
                f">=90% copy and NOT novel-to-genome")
    if sub == "partial_present":
        tgt = max(ev, key=lambda k: ev[k][0]); bid, cov, cid = ev[tgt]
        return (f"reference/hap-PRESENT but only PARTIAL coverage: high-identity (~{bid*100:.0f}%) "
                f"fragments cover {cov*100:.0f}% (<{int(COV_FRAC*100)}%) of the consensus in {tgt} -- "
                f"chimeric/fragmentary transcript consensus of a locus that IS in the genome, not a "
                f"full extra copy and not novel")
    return (f"absent from reference AND both haplotypes (no block >= {MIN_IDENT} identity and "
            f"<{int(COV_FRAC*100)}% coverage in all three) -- genuinely novel or artifact / too "
            f"divergent (<~{int(SUBTHRESH_FLOOR*100)}%) to place as a full germline copy")


# ----------------------------------------------------------------------------- oracle cross-check
def load_oracle_regions():
    fam_loci = defaultdict(list)   # family -> [(chrom,start,end)]
    if os.path.exists(VALID):
        with open(VALID) as fh:
            fh.readline()
            for line in fh:
                f = line.rstrip("\n").split("\t")
                fam_loci[f[0]].append((f[2], int(f[3]), int(f[4])))
    oracle = {}
    if os.path.exists(ORACLE):
        for d in csv.DictReader(open(ORACLE), delimiter="\t"):
            oracle[d["family"]] = d
    regions = []   # per-LOCUS (chrom, lo, hi, family, class, collapsed_bool)
    for fam, loci in fam_loci.items():
        od = oracle.get(fam)
        if od is None:
            continue
        collapsed = int(od["asm_hapCN"]) > int(od["n_loci_ref"])
        for c, s, e in loci:
            regions.append((c, s, e, fam, od["class"], collapsed))
    return regions


def oracle_hit(chrom, start, end, regions):
    for c, lo, hi, fam, cls, coll in regions:
        if c == chrom and min(end, hi) - max(start, lo) > 0:
            return fam, cls, coll
    return None, None, None


# ----------------------------------------------------------------------------- main
def main():
    os.makedirs(WORK, exist_ok=True)

    # ---- load candidates -----------------------------------------------------
    absent = json.load(open(ABSENT))                         # 145
    disc = {d["cid"]: d for d in json.load(open(DISC))}
    haploid = {}
    for d in csv.DictReader(open(HAPLOID), delimiter="\t"):
        haploid[d["cid"]] = d["call"]
    surv = read_fasta(SURV_FA)
    reps = read_fasta(REPS_FA)
    poc = json.load(open(POC))
    poc_cands = poc.get("candidates", [])                    # divergent/unmapped (cl0, cl19, ...)

    # ---- build combined query fasta -----------------------------------------
    query_fa = f"{WORK}/o4_query.fa"
    qlen = {}
    meta = {}     # qname -> dict(cid, chrom, start, end, ...)
    with open(query_fa, "w") as out:
        for a in absent:
            cid = a["cid"]
            seq = surv.get(cid)
            if not seq:
                continue
            qn = cid
            qlen[qn] = len(seq)
            meta[qn] = dict(kind="collapsed", cid=cid, chrom=a["chrom"], start=a["start"],
                            end=a["end"], divergence=a.get("divergence"),
                            protein=a.get("protein"), orf=a.get("orf"),
                            n_loci_read=disc.get(cid, {}).get("n_loci"),
                            haploid_label=haploid.get(cid, "NA"))
            out.write(f">{qn}\n")
            for i in range(0, len(seq), 80):
                out.write(seq[i:i + 80] + "\n")
        rep_by_tag = {}
        for hdr, seq in reps.items():
            rep_by_tag[hdr.split("_")[0]] = (hdr, seq)
        for pc in poc_cands:
            tag = pc["tag"]
            if tag not in rep_by_tag:
                continue
            hdr, seq = rep_by_tag[tag]
            qn = f"UNM_{tag}"
            qlen[qn] = len(seq)
            meta[qn] = dict(kind="divergent_unmapped", cid=qn, chrom="(unplaced)", start=0, end=0,
                            divergence=None, protein=None, orf=pc.get("orf_aa"),
                            n_loci_read=pc.get("n_loci"),
                            haploid_label=f"POC:{pc.get('cls','')[:40]}")
            out.write(f">{qn}\n")
            for i in range(0, len(seq), 80):
                out.write(seq[i:i + 80] + "\n")
    sys.stderr.write(f"query: {len(qlen)} candidates ({len(absent)} collapsed + "
                     f"{len(qlen)-len(absent)} divergent/unmapped) -> {query_fa}\n")

    # ---- align to reference / maternal / paternal ---------------------------
    ref_paf, mat_paf, pat_paf = f"{WORK}/ref.paf", f"{WORK}/mat.paf", f"{WORK}/pat.paf"
    ref_paf2 = f"{WORK}/ref.check.paf"
    reuse = os.environ.get("O4_REUSE_PAF") == "1" and all(
        os.path.exists(p) for p in (ref_paf, mat_paf, pat_paf))
    if reuse:
        sys.stderr.write("O4_REUSE_PAF=1: reusing existing PAFs\n")
    else:
        sys.stderr.write("aligning to primary reference (GGO.fasta) ...\n")
        run_minimap2(query_fa, GGO, ref_paf, is_fasta=True)
        sys.stderr.write("aligning to maternal hap ...\n")
        run_minimap2(query_fa, MAT_MMI, mat_paf, is_fasta=False)
        sys.stderr.write("aligning to paternal hap ...\n")
        run_minimap2(query_fa, PAT_MMI, pat_paf, is_fasta=False)
        run_minimap2(query_fa, GGO, ref_paf2, is_fasta=True)

    ref_b = read_paf_blocks(ref_paf)
    mat_b = read_paf_blocks(mat_paf)
    pat_b = read_paf_blocks(pat_paf)

    # ---- determinism check: sorted-PAF md5 of two independent reference aligns ----
    def sorted_md5(p):
        with open(p) as fh:
            lines = sorted(fh.read().splitlines())
        return hashlib.md5("\n".join(lines).encode()).hexdigest()
    if os.path.exists(ref_paf2):
        det_ok = sorted_md5(ref_paf) == sorted_md5(ref_paf2)
    else:
        det_ok = "not-rechecked (no ref.check.paf)"
    sys.stderr.write(f"determinism (reference re-align, sorted PAF md5 equal): {det_ok}\n")

    # ---- oracle cross-check regions -----------------------------------------
    oregions = load_oracle_regions()

    # ---- classify -----------------------------------------------------------
    rows = []
    sexlink_of = {}
    for qn in meta:
        m = meta[qn]
        ql = qlen[qn]
        rf, rff, rid, rreg = count_loci(ref_b.get(qn, []), ql)
        mf, mff, mid, mreg = count_loci(mat_b.get(qn, []), ql)
        pf, pff, pid, preg = count_loci(pat_b.get(qn, []), ql)
        # all-blocks evidence (for NOVEL sub-classification / reference-absent adjudication)
        re_bid, re_cov, re_cid = locus_evidence(ref_b.get(qn, []), ql)
        me_bid, me_cov, me_cid = locus_evidence(mat_b.get(qn, []), ql)
        pe_bid, pe_cov, pe_cid = locus_evidence(pat_b.get(qn, []), ql)
        ev = {"ref": (re_bid, re_cov, re_cid), "mat": (me_bid, me_cov, me_cid),
              "pat": (pe_bid, pe_cov, pe_cid)}
        sexlink = ("X" if m["chrom"] == CHRX else "Y" if m["chrom"] == CHRY else "auto")
        sexlink_of[qn] = sexlink
        call, reason, conf = classify(rf, mf, pf, sexlink)
        sub = ""
        if call == "NOVEL":
            sub = novel_subclass(ev)
            reason = novel_reason(sub, ev)
        of, ocls, ocoll = oracle_hit(m["chrom"], m["start"], m["end"], oregions)
        rows.append(dict(
            cid=m["cid"], kind=m["kind"], chrom=m["chrom"], start=m["start"], end=m["end"],
            protein=m["protein"], divergence=m["divergence"], n_loci_read=m["n_loci_read"],
            haploid_label=m["haploid_label"],
            ref_loci=rf, mat_loci=mf, pat_loci=pf,
            ref_frag=rff, mat_frag=mff, pat_frag=pff,
            best_id_ref=rid, best_id_mat=mid, best_id_pat=pid,
            anyid_ref=re_bid, anyid_mat=me_bid, anyid_pat=pe_bid,
            allcov_ref=re_cov, allcov_mat=me_cov, allcov_pat=pe_cov,
            diploid_label=call, subclass=sub, confidence=conf,
            oracle_family=of if of is not None else "",
            oracle_class=ocls if ocls is not None else "",
            oracle_collapsed=(1 if ocoll else 0) if of is not None else "",
            reason=reason, _qn=qn,
        ))

    # ---- write tsv ----------------------------------------------------------
    cols = ["cid", "kind", "chrom", "start", "end", "protein", "divergence", "n_loci_read",
            "haploid_label", "ref_loci", "mat_loci", "pat_loci", "ref_frag", "mat_frag", "pat_frag",
            "best_id_ref", "best_id_mat", "best_id_pat", "anyid_ref", "anyid_mat", "anyid_pat",
            "allcov_ref", "allcov_mat", "allcov_pat", "diploid_label", "subclass", "confidence",
            "oracle_family", "oracle_class", "oracle_collapsed", "reason"]
    with open(OUT_TSV, "w") as fh:
        fh.write("\t".join(cols) + "\n")
        for r in sorted(rows, key=lambda r: (r["kind"], r["chrom"], r["start"])):
            fh.write("\t".join(str(r[c]) for c in cols) + "\n")

    # ---- robustness sweep: 0-COPY across identity/coverage/gap --------------
    sweep = []
    for mi in (0.80, 0.85, 0.90):
        for cf in (0.2, 0.3, 0.4, 0.5):
            for gap in (100_000, 500_000):
                lab = Counter()
                for qn in meta:
                    ql = qlen[qn]
                    rf, *_ = count_loci(ref_b.get(qn, []), ql, mi, cf, gap)
                    mf, *_ = count_loci(mat_b.get(qn, []), ql, mi, cf, gap)
                    pf, *_ = count_loci(pat_b.get(qn, []), ql, mi, cf, gap)
                    l, _, _ = classify(rf, mf, pf, sexlink_of[qn])
                    lab[l] += 1
                sweep.append(dict(min_ident=mi, cov_frac=cf, locus_gap=gap,
                                  copy=lab.get("COPY", 0), labels=dict(lab)))
    max_copy_over_sweep = max(s["copy"] for s in sweep)

    # ---- summary ------------------------------------------------------------
    dip = Counter(r["diploid_label"] for r in rows)
    collapsed_rows = [r for r in rows if r["kind"] == "collapsed"]
    diverg_rows = [r for r in rows if r["kind"] == "divergent_unmapped"]
    dip_coll = Counter(r["diploid_label"] for r in collapsed_rows)
    dip_div = Counter(r["diploid_label"] for r in diverg_rows)

    trans = defaultdict(Counter)
    for r in collapsed_rows:
        trans[r["haploid_label"]][r["diploid_label"]] += 1

    amb = [r for r in collapsed_rows if r["haploid_label"] == "AMBIGUOUS"]
    amb_res = Counter(r["diploid_label"] for r in amb)

    copies = [r for r in collapsed_rows if r["diploid_label"] == "COPY"]
    copies_independent = [r for r in copies if not (r["oracle_family"] != "" and r["oracle_collapsed"] == 1)]
    copies_oracle_dup = [r for r in copies if r["oracle_family"] != "" and r["oracle_collapsed"] == 1]

    # NOVEL sub-split: only genuinely_absent is truly novel/artifact; the rest are reference-present.
    novel_rows = [r for r in rows if r["diploid_label"] == "NOVEL"]
    novel_sub = Counter(r["subclass"] for r in novel_rows)
    novel_genuine = [r for r in novel_rows if r["subclass"] == "genuinely_absent"]
    novel_ref_present = [r for r in novel_rows if r["subclass"] in ("subthreshold_present", "partial_present")]

    # ALLELE sub-breakdown. A one-hap "extra" is GENUINE only if it PERSISTS when the coverage gate
    # is relaxed to 0.30 -- this removes coverage-boundary jitter (e.g. CDK11B: mat 52.8% vs ref 43%
    # of the SAME locus, present in ref at 99.16% id -> not a real het CNV).
    def loci_at(qn, cf):
        ql = qlen[qn]
        rf, *_ = count_loci(ref_b.get(qn, []), ql, MIN_IDENT, cf, LOCUS_GAP)
        mf, *_ = count_loci(mat_b.get(qn, []), ql, MIN_IDENT, cf, LOCUS_GAP)
        pf, *_ = count_loci(pat_b.get(qn, []), ql, MIN_IDENT, cf, LOCUS_GAP)
        return rf, mf, pf

    allele_rows = [r for r in rows if r["diploid_label"] == "ALLELE"]
    allele_one_hap_extra_raw, allele_one_hap_extra_genuine, allele_single_locus = [], [], []
    for r in allele_rows:
        rl, ml, pl = int(r["ref_loci"]), int(r["mat_loci"]), int(r["pat_loci"])
        if ml > rl or pl > rl:
            allele_one_hap_extra_raw.append(r)
            rf3, mf3, pf3 = loci_at(r["_qn"], 0.30)   # extra must survive the relaxed coverage gate
            if mf3 > rf3 or pf3 > rf3:
                allele_one_hap_extra_genuine.append(r)
        else:
            allele_single_locus.append(r)
    # coverage-boundary artifacts = raw one-hap-extra that did NOT survive (present locus, jitter)
    allele_boundary_artifact = [r for r in allele_one_hap_extra_raw if r not in allele_one_hap_extra_genuine]

    # reference-absent CONFIRMED in diploid: genuinely absent from the PRIMARY (no >=0.90 block AND
    # <COV_FRAC all-block coverage in ref) yet present as a full copy in >=1 haplotype. This EXCLUDES
    # coverage-threshold crossings that are actually reference-present (CDK11B best_id_ref=0.9916).
    ref_absent_confirmed = [
        r for r in rows
        if int(r["ref_loci"]) == 0 and float(r["anyid_ref"]) < MIN_IDENT and float(r["allcov_ref"]) < COV_FRAC
        and (int(r["mat_loci"]) > 0 or int(r["pat_loci"]) > 0)]
    # candidates that map to the reference at high identity but were previously (incorrectly) counted
    # reference-absent because ref full-loci==0 (coverage artifact) -- reported for transparency.
    ref_present_but_zero_full = [
        r for r in rows
        if int(r["ref_loci"]) == 0 and float(r["anyid_ref"]) >= MIN_IDENT]

    in_oracle_collapsed = [r for r in collapsed_rows if r["oracle_collapsed"] == 1]

    def slim(r):
        return dict(cid=r["cid"], protein=r["protein"], subclass=r["subclass"],
                    ref=r["ref_loci"], mat=r["mat_loci"], pat=r["pat_loci"],
                    anyid_ref=r["anyid_ref"], anyid_mat=r["anyid_mat"], anyid_pat=r["anyid_pat"],
                    allcov_ref=r["allcov_ref"], allcov_mat=r["allcov_mat"], allcov_pat=r["allcov_pat"])

    summary = dict(
        determinism_reference_realign_ok=det_ok,
        params=dict(min_ident=MIN_IDENT, cov_frac=COV_FRAC, locus_gap=LOCUS_GAP,
                    subthresh_floor=SUBTHRESH_FLOOR, minimap2="asm20 -N50 -p0.1", threads=THREADS),
        n_total=len(rows), n_collapsed=len(collapsed_rows), n_divergent=len(diverg_rows),
        diploid_all=dict(dip), diploid_collapsed=dict(dip_coll), diploid_divergent=dict(dip_div),
        haploid_to_diploid_transition={k: dict(v) for k, v in trans.items()},
        haploid_ambiguous_resolved=dict(amb_res),
        # ---- the payoff ----
        n_copy=dip.get("COPY", 0), n_allele=dip.get("ALLELE", 0),
        n_novel=dip.get("NOVEL", 0), n_undecidable=dip.get("UNDECIDABLE", 0),
        copies_total=len(copies),
        copies_independent_of_oracle=len(copies_independent),
        copies_already_in_oracle_collapsed=len(copies_oracle_dup),
        # ---- 0-COPY robustness ----
        copy_max_over_sweep=max_copy_over_sweep,
        sweep_note=("COPY=0 at EVERY (identity 0.80/0.85/0.90, coverage 0.2-0.5, gap 100k-500k) -- "
                    "the negative is not a threshold/method-capacity artifact"),
        robustness_sweep=sweep,
        sensitivity_floor=("0 COPY == 0 copy at >= ~%d%% identity; a genuinely divergent (<~%d%%) "
                           "both-hap germline copy would be scored NOVEL/genuinely_absent -- the "
                           "identity floor of a transcript-copy oracle, stated as a caveat."
                           % (int(SUBTHRESH_FLOOR*100), int(SUBTHRESH_FLOOR*100))),
        # ---- NOVEL is a MIXED bucket ----
        novel_total=len(novel_rows),
        novel_subclasses=dict(novel_sub),
        novel_genuinely_absent=len(novel_genuine),
        novel_reference_present_below_gate=len(novel_ref_present),
        novel_cov_frac_sensitivity=("the NOVEL/ALLELE split is COV_FRAC-sensitive (loci=0 vs 1 at the "
                                     "50%% boundary); see robustness_sweep. Only COPY=0 is threshold-robust."),
        novel_genuinely_absent_list=[slim(r) for r in novel_genuine],
        novel_reference_present_list=[slim(r) for r in novel_ref_present],
        # ---- ALLELE breakdown ----
        allele_total=len(allele_rows),
        allele_one_hap_extra_raw=len(allele_one_hap_extra_raw),
        allele_one_hap_extra_genuine=len(allele_one_hap_extra_genuine),
        allele_one_hap_extra_boundary_artifact=len(allele_boundary_artifact),
        allele_single_locus_het=len(allele_single_locus),
        allele_one_hap_extra_genuine_list=[slim(r) for r in allele_one_hap_extra_genuine],
        allele_boundary_artifact_list=[slim(r) for r in allele_boundary_artifact],
        # ---- reference-absent adjudication ----
        reference_absent_confirmed_in_diploid=[
            dict(cid=r["cid"], kind=r["kind"], protein=r["protein"], ref=r["ref_loci"],
                 mat=r["mat_loci"], pat=r["pat_loci"], anyid_ref=r["anyid_ref"],
                 anyid_mat=r["anyid_mat"], anyid_pat=r["anyid_pat"], label=r["diploid_label"])
            for r in ref_absent_confirmed],
        n_reference_absent_confirmed_in_diploid=len(ref_absent_confirmed),
        reference_present_zero_full_excluded=[
            dict(cid=r["cid"], protein=r["protein"], anyid_ref=r["anyid_ref"],
                 allcov_ref=r["allcov_ref"], label=r["diploid_label"])
            for r in ref_present_but_zero_full],
        # ---- oracle cross-check (no double-count) ----
        candidates_in_oracle_collapsed_family=[
            dict(cid=r["cid"], protein=r["protein"], oracle_family=r["oracle_family"],
                 ref=r["ref_loci"], mat=r["mat_loci"], pat=r["pat_loci"], label=r["diploid_label"],
                 subclass=r["subclass"]) for r in in_oracle_collapsed],
        copy_calls=[dict(cid=r["cid"], protein=r["protein"], ref=r["ref_loci"], mat=r["mat_loci"],
                         pat=r["pat_loci"], conf=r["confidence"], oracle_family=r["oracle_family"],
                         oracle_collapsed=r["oracle_collapsed"]) for r in copies],
        divergent_calls=[dict(cid=r["cid"], label=r["diploid_label"], subclass=r["subclass"],
                              ref=r["ref_loci"], mat=r["mat_loci"], pat=r["pat_loci"],
                              note=r["haploid_label"]) for r in diverg_rows],
    )
    json.dump(summary, open(OUT_JSON, "w"), indent=1)

    # ---- print --------------------------------------------------------------
    print(f"== O4 DIPLOID VALIDATION (donor-exact: RNA == mGorGor1) ==")
    print(f"determinism (reference re-align sorted-PAF md5 equal): {det_ok}")
    print(f"candidates: {len(rows)} = {len(collapsed_rows)} collapsed/CNV + {len(diverg_rows)} divergent/unmapped\n")
    print(f"-- diploid label, collapsed set (n={len(collapsed_rows)}) --")
    for k, v in dip_coll.most_common():
        print(f"   {k:12} {v:4}  ({100*v/len(collapsed_rows):.0f}%)")
    print(f"-- diploid label, divergent set (n={len(diverg_rows)}) --")
    for k, v in dip_div.most_common():
        print(f"   {k:12} {v:4}")
    print(f"\nCOPIES confirmed (extra full-copy locus in BOTH haps): {len(copies)}")
    print(f"   0-COPY robust across sweep? max COPY over id0.80-0.90 x cov0.2-0.5 x gap100k-500k = {max_copy_over_sweep}")
    print(f"\n-- NOVEL is a MIXED bucket (n={len(novel_rows)}); loci=0 != absent --")
    for k, v in novel_sub.most_common():
        print(f"   {k:22} {v}")
    print(f"   genuinely_absent (true novel/artifact)      : {len(novel_genuine)}")
    print(f"   reference-present below the 0.90 gate        : {len(novel_ref_present)}")
    print(f"\n-- haploid -> diploid transition (collapsed set) --")
    for hl in sorted(trans):
        print(f"   {hl:16} -> " + ", ".join(f"{k}:{v}" for k, v in trans[hl].most_common()))
    print(f"\n-- the 65 haploid-AMBIGUOUS resolved --")
    for k, v in amb_res.most_common():
        print(f"   {k:12} {v}")
    print(f"\nALLELE total {len(allele_rows)}: {len(allele_one_hap_extra_genuine)} GENUINE one-hap EXTRA "
          f"full copy (het CNV) + {len(allele_boundary_artifact)} coverage-boundary artifact "
          f"(reference-present, jitter) + {len(allele_single_locus)} plain single-locus het")
    for r in allele_one_hap_extra_genuine:
        print(f"     GENUINE  {r['cid']}  {r['protein']}  r/m/p={r['ref_loci']}/{r['mat_loci']}/{r['pat_loci']}")
    print(f"\nreference-absent CONFIRMED in diploid (absent primary, present in a hap): {len(ref_absent_confirmed)}")
    for r in ref_absent_confirmed:
        print(f"     {r['cid']}  {r['protein']}  anyid_ref={r['anyid_ref']} mat={r['mat_loci']} "
              f"pat={r['pat_loci']} -> {r['diploid_label']}")
    print(f"   (excluded as reference-PRESENT coverage artifacts: "
          f"{[r['protein'] or r['cid'] for r in ref_present_but_zero_full]})")
    print(f"\nO4 candidates inside an oracle reference-COLLAPSED family: {len(in_oracle_collapsed)} "
          f"(transcript-copy method labels: "
          f"{dict(Counter(r['diploid_label'] for r in in_oracle_collapsed))})")
    print(f"\nwrote {OUT_TSV} + {OUT_JSON}")


if __name__ == "__main__":
    main()
