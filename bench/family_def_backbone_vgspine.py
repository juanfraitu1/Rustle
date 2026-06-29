#!/usr/bin/env python3
"""family_def_backbone_vgspine.py — A3 (sequence-anchored): build the family VG
SPINE properly and measure junction-backbone CO-THREADING in the ALIGNED frame.

The exon-length-signature version (family_def_backbone_cothread.py) is too brittle:
it requires exon *lengths* to match, so real paralogs with length-drifted UTRs
(RFPL1/2/3) score low even though they share junction TOPOLOGY. The VARIATION-GRAPH
view fixes this: copies are parallel PATHS through one bubble chain anchored by a
SEQUENCE alignment, so a junction is "co-threaded" if it lands at a HOMOLOGOUS
sequence position of the reference copy's junction, regardless of length drift.

CONSTRUCTION (annotation-derived, no `rustle --vg` genome-wide):
  1. For each member gene take its representative transcript (most exons; tie->span).
     Build its SPLICED cDNA = concatenated exon genome-sequence (strand-correct).
     Record junction positions as cDNA offsets (cumulative exon lengths).
  2. Reference spine = member with the most junctions. Align every other member's
     cDNA to the reference cDNA with minimap2 (asm20, base-level CIGAR).
  3. PROJECT each member junction (cDNA offset) through the alignment to a reference
     cDNA coordinate. A member junction CO-THREADS the backbone if its projected
     position is within TOL bp of a reference junction AND the alignment actually
     covers that position (the two copies share that bubble boundary).
  4. shared_backbone_fraction = (# reference junctions co-threaded by ALL members)
                                / (# reference junctions).

This is the real VG criterion: "members are parallel paths sharing a junction
backbone (bubble chain) in one variation graph." Retrocopies (0 junctions) and
domain-sharers (junctions only over a sub-path) fail it; whole-gene-dup families
pass it even when lengths drift.

Run: /home/juanfra/miniforge3/bin/python bench/family_def_backbone_vgspine.py
Output: bench/family_def_backbone_vgspine.json
"""
import collections
import json
import os
import re
import subprocess
import tempfile

import pysam

HERE = os.path.dirname(os.path.abspath(__file__))
GFF = "/mnt/c/Users/jfris/Desktop/GGO_genomic.gff"
GENOME = "/home/juanfra/winloci_scratch/GGO.fasta"
REGION_CACHE = "/tmp/panel_gene_regions.txt"
FEATURE_CACHE = "/tmp/panel_features.gff"
TMP = tempfile.mkdtemp(prefix="vgspine_")

JUNC_TOL = 12          # bp tolerance for a projected junction to match a ref junction
CIG = re.compile(r"(\d+)([MIDNSHP=X])")

PANEL = [
    ("RABL2",          "real",   ["RABL2A", "RABL2B"]),
    ("APOBEC3",        "real",   ["APOBEC3D", "APOBEC3F"]),
    ("RFPL",           "real",   ["RFPL1", "RFPL2", "RFPL3"]),
    ("CASP8~FLACC1",   "domain", ["CASP8", "FLACC1"]),
    ("ASDURF~ASNSD1",  "domain", ["ASDURF", "ASNSD1"]),
    ("CREB1~METTL21A", "domain", ["CREB1", "METTL21A"]),
    ("GPR39~LYPD1",    "domain", ["GPR39", "LYPD1"]),
    ("EEF1A1_retro",   "retro",  ["EEF1A1", "LOC109023808"]),
    ("CNN2_retro",     "retro",  ["CNN2", "LOC129524764"]),
]

_REGIONS = None
_FEATURES = None  # (mrna_parent, mrna_strand, exons_by_mrna)


def _attr(attr, key):
    for kv in attr.split(";"):
        if kv.startswith(key + "="):
            return kv[len(key) + 1:]
    return None


def _load_caches():
    global _REGIONS, _FEATURES
    if _REGIONS is not None:
        return
    _REGIONS = {}
    for line in open(REGION_CACHE):
        name, chrom, s, e = line.rstrip("\n").split("\t")
        _REGIONS[name] = (chrom, int(s), int(e))
    mrna_parent, mrna_strand = {}, {}
    exons = collections.defaultdict(list)
    for line in open(FEATURE_CACHE):
        f = line.rstrip("\n").split("\t")
        typ, s, e, strand, attr = f[2], int(f[3]), int(f[4]), f[6], f[8]
        if typ == "mRNA":
            mid = _attr(attr, "ID")
            mrna_parent[mid] = _attr(attr, "Parent")
            mrna_strand[mid] = strand
        elif typ == "exon":
            exons[_attr(attr, "Parent")].append((s, e))
    _FEATURES = (mrna_parent, mrna_strand, exons)


def representative(name):
    _load_caches()
    mrna_parent, mrna_strand, exons = _FEATURES
    want = f"gene-{name}"
    best = None
    for mid, exs in exons.items():
        if mrna_parent.get(mid) != want:
            continue
        exs = sorted(exs)
        span = exs[-1][1] - exs[0][0]
        key = (len(exs), span)
        if best is None or key > best[0]:
            best = (key, mid, exs, mrna_strand.get(mid, "+"))
    if best is None:
        return None
    return best[1], best[2], best[3]   # mid, exons(sorted asc), strand


def revcomp(s):
    return s.translate(str.maketrans("ACGTNacgtn", "TGCANtgcan"))[::-1]


def spliced_cdna(name, fa):
    """Return (cdna_seq, [junction offsets in cDNA]) for the representative tx."""
    rep = representative(name)
    if rep is None:
        return None
    _load_caches()
    chrom = _REGIONS[name][0]
    _, exons, strand = rep
    blocks = [fa.fetch(chrom, a - 1, b).upper() for a, b in exons]  # GFF is 1-based
    if strand == "-":
        blocks = [revcomp(b) for b in reversed(blocks)]
    cdna = "".join(blocks)
    junc = []
    cum = 0
    for b in blocks[:-1]:
        cum += len(b)
        junc.append(cum)
    return cdna, junc, strand


def project_junctions(ref_seq, ref_junc, mem_seq, mem_junc):
    """Align member cDNA to reference cDNA; project member junction offsets to ref
    coordinates via the CIGAR. Return (co_threaded_ref_juncs:set, aln_cov_frac)."""
    refp = os.path.join(TMP, "ref.fa")
    memp = os.path.join(TMP, "mem.fa")
    open(refp, "w").write(f">ref\n{ref_seq}\n")
    open(memp, "w").write(f">mem\n{mem_seq}\n")
    # map member onto reference; secondary off, base CIGAR
    res = subprocess.run(
        ["minimap2", "-a", "--eqx", "-x", "asm20", "--secondary=no", refp, memp],
        capture_output=True, text=True)
    # parse the primary alignment SAM
    aln = None
    for line in res.stdout.splitlines():
        if line.startswith("@"):
            continue
        f = line.split("\t")
        flag = int(f[1])
        if flag & 0x4 or flag & 0x100 or flag & 0x800:
            continue
        aln = f
        break
    if aln is None:
        return set(), 0.0
    rstart = int(aln[3]) - 1     # 0-based ref start
    cigar = aln[5]
    is_rev = bool(int(aln[1]) & 0x10)
    # member query coords in the SAM are in the oriented read; if reverse, junctions
    # must be flipped to oriented coords.
    mlen = len(mem_seq)
    mj = [mlen - j for j in reversed(mem_junc)] if is_rev else list(mem_junc)
    # build query->ref position map by walking CIGAR
    qpos = 0
    rpos = rstart
    q2r = {}   # query offset (oriented) -> ref offset
    covered_q = 0
    for n, op in CIG.findall(cigar):
        n = int(n)
        if op in "=XM":
            for k in range(n):
                q2r[qpos + k] = rpos + k
            qpos += n
            rpos += n
            covered_q += n
        elif op in "ID":
            if op == "I":
                qpos += n
            else:
                rpos += n
        elif op in "SH":
            qpos += n
    # project member junctions
    co = set()
    ref_set = sorted(ref_junc)
    for j in mj:
        # nearest mapped query position
        rj = None
        for d in range(0, JUNC_TOL + 1):
            if j - d in q2r:
                rj = q2r[j - d]; break
            if j + d in q2r:
                rj = q2r[j + d]; break
        if rj is None:
            continue
        for rk in ref_set:
            if abs(rj - rk) <= JUNC_TOL:
                co.add(rk)
                break
    cov = covered_q / max(mlen, 1)
    return co, cov


def run_case(case, kind, members):
    fa = pysam.FastaFile(GENOME)
    recs = []
    for m in members:
        sc = spliced_cdna(m, fa)
        if sc is None:
            recs.append(dict(gene=m, seq="", junc=[], strand="?"))
            continue
        cdna, junc, strand = sc
        recs.append(dict(gene=m, seq=cdna, junc=junc, strand=strand, clen=len(cdna)))
    fa.close()
    ref = max(recs, key=lambda r: len(r["junc"]))
    ref_junc = ref["junc"]
    nref = len(ref_junc)
    per = []
    co_all = set(ref_junc)
    for r in recs:
        if r is ref:
            per.append(dict(gene=r["gene"], njunc=len(r["junc"]), co=nref,
                            frac=1.0 if nref else 0.0, aln_cov=1.0, is_ref=True))
            continue
        if not r["junc"] or nref == 0:
            co, cov = set(), 0.0
        else:
            co, cov = project_junctions(ref["seq"], ref_junc, r["seq"], r["junc"])
        per.append(dict(gene=r["gene"], njunc=len(r["junc"]), co=len(co),
                        frac=len(co) / nref if nref else 0.0,
                        aln_cov=round(cov, 3), is_ref=False))
        co_all &= co
    shared = len(co_all) / nref if nref else 0.0
    return dict(case=case, kind=kind, ref_gene=ref["gene"], ref_junctions=nref,
                shared_backbone_junctions=len(co_all),
                shared_backbone_fraction=round(shared, 3),
                min_member_fraction=round(min(p["frac"] for p in per), 3),
                any_intronless=any(len(r["junc"]) == 0 for r in recs),
                members=per,
                junc_counts={r["gene"]: len(r["junc"]) for r in recs})


def main():
    results = [run_case(*c) for c in PANEL]
    print("=== A3 (sequence-anchored VG spine): junction-backbone co-threading ===\n")
    print(f"{'case':16} {'kind':6} {'refJ':>4} {'shared':>6} {'minM':>5} {'intronless':>10}"
          f"  members(njunc,co,frac,alnCov)")
    for r in results:
        ms = "  ".join(f"{p['gene']}({p['njunc']},{p['co']},{p['frac']:.2f},{p['aln_cov']:.2f})"
                       for p in r["members"])
        print(f"{r['case']:16} {r['kind']:6} {r['ref_junctions']:>4} "
              f"{r['shared_backbone_fraction']:>6.2f} {r['min_member_fraction']:>5.2f} "
              f"{str(r['any_intronless']):>10}  {ms}")

    real = sorted(r["shared_backbone_fraction"] for r in results if r["kind"] == "real")
    domain = sorted(r["shared_backbone_fraction"] for r in results if r["kind"] == "domain")
    retro = sorted(r["shared_backbone_fraction"] for r in results if r["kind"] == "retro")
    print("\n--- shared_backbone_fraction by class ---")
    print(f"  REAL   : {real}   range [{min(real):.2f},{max(real):.2f}]")
    print(f"  DOMAIN : {domain}   range [{min(domain):.2f},{max(domain):.2f}]")
    print(f"  RETRO  : {retro}   range [{min(retro):.2f},{max(retro):.2f}]")
    m_rd = min(real) - max(domain)
    m_rr = min(real) - max(retro)
    print(f"\n  margin REAL_min - DOMAIN_max = {m_rd:+.3f}  (>0 => clean separation)")
    print(f"  margin REAL_min - RETRO_max  = {m_rr:+.3f}")

    json.dump(dict(panel=results, junc_tol=JUNC_TOL,
                   summary=dict(real=real, domain=domain, retro=retro,
                                margin_real_minus_domain=round(m_rd, 3),
                                margin_real_minus_retro=round(m_rr, 3))),
              open(os.path.join(HERE, "family_def_backbone_vgspine.json"), "w"), indent=2)
    print(f"\n[+] wrote {os.path.join(HERE, 'family_def_backbone_vgspine.json')}")


if __name__ == "__main__":
    main()
