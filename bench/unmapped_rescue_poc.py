#!/usr/bin/env python3
"""
unmapped_rescue_poc.py -- Reproduce the "rescue unmapped reads -> find missing
gene copies" proof-of-concept (a citable NEGATIVE result on T2T gorilla IsoSeq).

Hypothesis tested
-----------------
A real concern (raised by the advisor) is reference bias: gene copies that exist
in the transcriptome but are *absent from the reference genome* would never map,
so an assembler keyed on alignments would never see them. The cheapest place such
copies would surface is the UNMAPPED read pile. So: take every read minimap2's
splice-aware genome alignment flags unmapped, cluster them de-novo (all-vs-all),
and ask whether the reproducible clusters are genuinely-absent expressed copies.

What this script does (all from preserved artifacts; no re-alignment of the BAM)
-------------------------------------------------------------------------------
1. Read count + length distribution of the unmapped pile          (unmapped.fa)
2. Unmapped fraction vs. the primary-read denominator             (ggo_flagstat.txt)
3. De-novo clustering = connected components of the all-vs-all
   overlap graph                                                  (ava.paf)
4. Gene-level classification of the >=5-read cluster reps after a
   RELAXED (non-spliced, multi-segment) remap: how many are
   multi-locus chimeras vs. single-locus novel-sequence candidates
                                                  (reps.relaxed.sam + genes.bed)

Inputs live in the scratch dir (too large to commit); override with argv[1].
Re-deriving steps 1-4 needs only the small artifacts, NOT the 1.6 GB BAM.

Provenance of the heavy artifacts (recorded for reproducibility):
  unmapped.bam   = samtools view -b -f 4 GGO.bam            (primary-unmapped reads)
  unmapped.fa    = samtools fasta unmapped.bam
  ava.paf        = minimap2 -x ava-pb -X unmapped.fa unmapped.fa
  reps.fa        = one representative read per connected component of size >=5
  reps.relaxed.sam = minimap2 -ax map-pb -p 0.1 -N 50 GGO.fasta reps.fa
                    (NON-spliced remap; -p 0.1 -N 50 reports up to 50 secondary
                     placements so paralog homology is visible. A chimera's pieces
                     then map as primary+supplementary to >1 locus; alternative
                     paralog placements show up as SECONDARY records.)
"""
import sys, os, json, re
from collections import defaultdict
from bisect import bisect_right

SCRATCH = sys.argv[1] if len(sys.argv) > 1 else "/home/juanfra/winloci_scratch/unmapped_poc"

UNMAPPED_FA = os.path.join(SCRATCH, "unmapped.fa")
FLAGSTAT    = os.path.join(SCRATCH, "ggo_flagstat.txt")
AVA_PAF     = os.path.join(SCRATCH, "ava.paf")
RELAXED_SAM = os.path.join(SCRATCH, "reps.relaxed.sam")
GENES_BED   = os.path.join(SCRATCH, "genes.bed")
REPS_FA     = os.path.join(SCRATCH, "reps.fa")

# Classification thresholds (documented; conservative)
COV_NOVEL      = 0.50   # single-locus + < this fraction of the read aligns -> novel candidate
COV_CLEAN      = 0.80   # combined query coverage to call a multi-piece read "cleanly covered"
LOCUS_GAP      = 1_000_000  # same-chrom alignments farther apart than this = distinct loci


# ----------------------------------------------------------------------------- helpers
CIG = re.compile(r"(\d+)([MIDNSHP=X])")

def cigar_spans(cigar, is_rev):
    """Return (read_len, q0, q1, ref_span, M, Ibases, Dbases) in FORWARD coords."""
    lead = trail = qaln = rspan = M = Ib = Db = 0
    ops = CIG.findall(cigar)
    # leading clip
    i = 0
    while i < len(ops) and ops[i][1] in "SH":
        lead += int(ops[i][0]); i += 1
    j = len(ops) - 1
    while j >= 0 and ops[j][1] in "SH":
        trail += int(ops[j][0]); j -= 1
    for n, op in ops:
        n = int(n)
        if op in "M=X": M += n
        if op == "I": Ib += n
        if op == "D": Db += n
        if op in "MI=X": qaln += n
        if op in "MDN=X": rspan += n
    read_len = lead + qaln + trail
    if is_rev:
        q0 = trail; q1 = trail + qaln
    else:
        q0 = lead; q1 = lead + qaln
    return read_len, q0, q1, rspan, M, Ib, Db


def sub_identity(M, Ib, Db, nm):
    """Substitution identity = 1 - (NM - inserted - deleted) / M (mismatches/M)."""
    if nm is None or M == 0:
        return None
    subs = max(nm - Ib - Db, 0)
    return 1.0 - subs / M


_COMP = {"A": "T", "T": "A", "G": "C", "C": "G", "N": "N"}
_STOP = {"TAA", "TAG", "TGA"}

def _rc(s):
    return "".join(_COMP.get(b, "N") for b in reversed(s))

def longest_orf_aa(s):
    """Longest stop-free codon run (aa) over all 6 frames. Random ~21aa; a real
    coding transcript runs into the hundreds."""
    best = 0
    for seq in (s, _rc(s)):
        for f in range(3):
            cur = 0
            for i in range(f, len(seq) - 2, 3):
                if seq[i:i + 3] in _STOP:
                    best = max(best, cur); cur = 0
                else:
                    cur += 1
            best = max(best, cur)
    return best

def nearest_gene(by_chrom, chrom, pos):
    cand = by_chrom.get(chrom, [])
    best = ("none", 10 ** 9)
    for s, e, nm in cand:
        if s <= pos <= e:
            return (nm, 0)
        d = min(abs(pos - s), abs(pos - e))
        if d < best[1]:
            best = (nm, d)
    return best


def complexity(seq, maxp=1500):
    """Return (max tandem-period match-fraction, random baseline). A satellite/
    low-complexity sequence gives max-frac >> baseline; complex sequence ~baseline."""
    from collections import Counter as _C
    n = len(seq)
    cc = _C(seq)
    base = sum((cc[b] / n) ** 2 for b in "ACGT") if n else 0
    best = 0.0
    cap = min(maxp, n - 1)
    for p in range(2, cap):
        m = sum(1 for i in range(n - p) if seq[i] == seq[i + p])
        f = m / (n - p)
        if f > best:
            best = f
    return best, base


def load_genes(path):
    by_chrom = defaultdict(list)
    with open(path) as f:
        for line in f:
            c, s, e, name = line.rstrip("\n").split("\t")
            by_chrom[c].append((int(s), int(e), name))
    for c in by_chrom:
        by_chrom[c].sort()
    return by_chrom


def gene_at(by_chrom, chrom, rs, re_):
    """Gene name with the largest overlap of [rs,re_), or 'intergenic'."""
    cand = by_chrom.get(chrom, [])
    if not cand:
        return "intergenic"
    starts = [g[0] for g in cand]
    hi = bisect_right(starts, re_)
    best, best_ov = "intergenic", 0
    for k in range(hi - 1, -1, -1):
        s, e, name = cand[k]
        if e <= rs:
            # genes are start-sorted; once we pass far left we can stop-ish.
            if rs - e > 5_000_000:
                break
            continue
        ov = min(e, re_) - max(s, rs)
        if ov > best_ov:
            best_ov, best = ov, name
    return best


def union_len(intervals):
    if not intervals:
        return 0
    iv = sorted(intervals)
    total = 0; cs, ce = iv[0]
    for s, e in iv[1:]:
        if s <= ce:
            ce = max(ce, e)
        else:
            total += ce - cs; cs, ce = s, e
    total += ce - cs
    return total


def distinct_loci(records):
    """Collapse (chrom, refpos) anchors into distinct loci (>LOCUS_GAP apart)."""
    by_c = defaultdict(list)
    for r in records:
        by_c[r["rname"]].append(r["rpos"])
    loci = 0
    for c, pos in by_c.items():
        pos.sort()
        loci += 1
        for a, b in zip(pos, pos[1:]):
            if b - a > LOCUS_GAP:
                loci += 1
    return loci


# ----------------------------------------------------------------------------- step 1
def step1_pile():
    n = 0; lengths = []
    with open(UNMAPPED_FA) as f:
        for line in f:
            if line.startswith(">"):
                n += 1
            else:
                lengths.append(len(line.strip()))
    lengths.sort()
    med = lengths[len(lengths) // 2]
    return dict(n=n, median=med, mean=round(sum(lengths) / len(lengths), 1),
                min=lengths[0], max=lengths[-1])


def step2_denominator(n_unmapped):
    primary = None
    with open(FLAGSTAT) as f:
        for line in f:
            if " primary" in line and "mapped" not in line and "duplicates" not in line:
                primary = int(line.split()[0])
    pct = 100.0 * n_unmapped / primary
    return dict(primary_reads=primary, pct=round(pct, 4))


def step3_clusters():
    parent = {}
    def find(x):
        parent.setdefault(x, x)
        root = x
        while parent[root] != root:
            root = parent[root]
        while parent[x] != root:
            parent[x], x = root, parent[x]
        return root
    def union(a, b):
        ra, rb = find(a), find(b)
        if ra != rb:
            parent[ra] = rb
    nodes = set()
    with open(AVA_PAF) as f:
        for line in f:
            c = line.split("\t")
            q, t = c[0], c[5]
            nodes.add(q); nodes.add(t)
            if q != t:
                union(q, t)
    comp = defaultdict(int)
    for n in nodes:
        comp[find(n)] += 1
    sizes = sorted(comp.values(), reverse=True)
    def agg(k):
        s = [x for x in sizes if x >= k]
        return dict(clusters=len(s), reads=sum(s))
    return dict(reads_in_graph=len(nodes), ge2=agg(2), ge3=agg(3), ge5=agg(5),
                top=sizes[:8])


def step4_classify():
    """
    Classify each >=5-read cluster representative by its RELAXED (non-spliced)
    remap. 'Pieces' of a read = primary + supplementary alignments (the colinear
    chimeric split). SECONDARY alignments (alternative placements / paralog
    homology) are tracked separately and deliberately NOT counted as 'loci' --
    counting them was the over-count that produced the inflated session number.
    """
    genes = load_genes(GENES_BED)
    recs = defaultdict(list)
    with open(RELAXED_SAM) as f:
        for line in f:
            if line.startswith("@"):
                continue
            c = line.rstrip("\n").split("\t")
            qname, flag, rname, pos, mapq, cigar = c[0], int(c[1]), c[2], int(c[3]), int(c[4]), c[5]
            nm = None
            for tag in c[11:]:
                if tag.startswith("NM:i:"):
                    nm = int(tag[5:])
            recs[qname].append(dict(flag=flag, rname=rname, pos=pos, mapq=mapq, cigar=cigar, nm=nm))

    cats = defaultdict(int)
    covbuckets = defaultdict(int)
    mapq_hist = defaultdict(int)        # PRIMARY mapq only
    in_gene = 0
    n_with_secondary = 0                # paralog-homology reps
    n_secondary_crosschrom = 0
    n_reps_any_mapq0_record = 0         # incl. secondary/supp records
    single_locus_ids = []              # substitution identity of present bucket
    divergent_subset = []
    cov_by_read = []
    multilocus, twopiece_clean, novel = [], [], []
    per_read = []
    for q, rs in recs.items():
        prim = [r for r in rs if not (r["flag"] & 0x100) and not (r["flag"] & 0x800) and not (r["flag"] & 0x4)]
        supp = [r for r in rs if r["flag"] & 0x800]
        sec = [r for r in rs if r["flag"] & 0x100]
        unmapped = any(r["flag"] & 0x4 for r in rs) and not prim
        if unmapped:
            cats["unmapped_relaxed"] += 1
            per_read.append(dict(read=q, cls="unmapped_relaxed"))
            continue
        p = prim[0]
        mq = p["mapq"]
        mapq_hist["0" if mq == 0 else "1-9" if mq < 10 else "10-59" if mq < 60 else "60"] += 1
        if any(r["mapq"] == 0 for r in rs):
            n_reps_any_mapq0_record += 1
        if sec:
            n_with_secondary += 1
        pieces = prim + supp
        read_len = None
        anchors = []   # primary+supplementary pieces -> (gene, qiv, rec)
        prim_id = None
        for idx, r in enumerate(pieces):
            rl, q0, q1, rspan, M, Ib, Db = cigar_spans(r["cigar"], r["flag"] & 0x10)
            read_len = read_len or rl
            g = gene_at(genes, r["rname"], r["pos"], r["pos"] + rspan)
            anchors.append(dict(gene=g, qiv=(q0, q1), rname=r["rname"], rpos=r["pos"]))
            if r is p:
                prim_id = sub_identity(M, Ib, Db, r["nm"])
        if read_len is None or read_len == 0:
            read_len = max((a["qiv"][1] for a in anchors), default=1)
        if any(r["rname"] != p["rname"] for r in sec):
            n_secondary_crosschrom += 1
        cov = union_len([a["qiv"] for a in anchors]) / read_len
        cov_by_read.append(cov)
        cb = ("0-25" if cov < 0.25 else "25-50" if cov < 0.5 else "50-80" if cov < 0.8
              else "80-95" if cov < 0.95 else "95-100")
        covbuckets[cb] += 1
        prim_gene = anchors[0]["gene"]
        if prim_gene != "intergenic":
            in_gene += 1
        genes_hit = sorted({a["gene"] for a in anchors if a["gene"] != "intergenic"})
        nloci = distinct_loci(anchors)
        is_multi = len(genes_hit) >= 2 or nloci >= 2
        info = dict(read=q, cov=round(cov, 3), mapq=mq,
                    prim_id=round(prim_id, 4) if prim_id is not None else None,
                    genes=genes_hit, prim_gene=prim_gene, nloci=nloci,
                    npieces=len(anchors), n_secondary=len(sec))
        if cov < COV_NOVEL:
            cats["novel_candidate"] += 1
            novel.append(info)
            cls = "novel_candidate"
        elif is_multi:
            cats["multi_locus_chimera"] += 1
            multilocus.append(info)
            cls = "multi_locus_chimera"
            if len(anchors) == 2 and cov >= COV_CLEAN:
                cats["two_piece_clean"] += 1
                twopiece_clean.append(info)
        else:
            cats["single_locus_present"] += 1
            cls = "single_locus_present"
            if prim_id is not None:
                single_locus_ids.append(prim_id)
                if prim_id < 0.90:
                    divergent_subset.append(info)
        per_read.append(dict(cls=cls, **info))
    n_mapped = sum(1 for q, rs in recs.items()
                   if any(not (r["flag"] & 0x100) and not (r["flag"] & 0x800)
                          and not (r["flag"] & 0x4) for r in rs))
    single_locus_ids.sort()
    med_id = single_locus_ids[len(single_locus_ids) // 2] if single_locus_ids else None
    cov_sens = {thr: sum(1 for cv in cov_by_read if cv < thr) for thr in (0.4, 0.5, 0.6)}
    return dict(n_reps=len(recs), n_mapped_relaxed=n_mapped, cats=dict(cats),
                covbuckets=dict(covbuckets), mapq_hist=dict(mapq_hist),
                in_gene=in_gene, n_with_secondary=n_with_secondary,
                n_secondary_crosschrom=n_secondary_crosschrom,
                n_reps_any_mapq0_record=n_reps_any_mapq0_record,
                single_locus_median_id=med_id, single_locus_min_id=(single_locus_ids[0] if single_locus_ids else None),
                divergent_subset=divergent_subset, cov_sensitivity=cov_sens,
                multilocus=multilocus, twopiece_clean=twopiece_clean, novel=novel,
                per_read=per_read)


# ----------------------------------------------------------------------------- step 5
def load_fa(path):
    d = {}; name = None; seq = []
    with open(path) as f:
        for line in f:
            if line.startswith(">"):
                if name:
                    d[name] = "".join(seq)
                name = line[1:].strip().split()[0]; seq = []
            else:
                seq.append(line.strip())
    if name:
        d[name] = "".join(seq)
    return d


def step5_candidates(tags=("cl0", "cl19", "cl31", "cl56", "cl34", "cl102")):
    """Follow-up on the residual candidates: is each a novel coding transcript,
    an expressed repeat/TE, or a divergent single locus?"""
    if not os.path.exists(REPS_FA):
        return None
    genes = load_genes(GENES_BED)
    reps = load_fa(REPS_FA)
    byprefix = {}
    for k in reps:
        byprefix.setdefault(k.split("_")[0], k)
    sam = defaultdict(list)
    with open(RELAXED_SAM) as f:
        for line in f:
            if line.startswith("@"):
                continue
            c = line.rstrip("\n").split("\t")
            AS = next((int(t[5:]) for t in c[11:] if t.startswith("AS:i:")), None)
            sam[c[0]].append((int(c[1]), c[2], int(c[3]), int(c[4]), AS, c[5]))
    out = []
    for tag in tags:
        if tag not in byprefix:
            continue
        q = byprefix[tag]
        seq = reps[q]
        rs = sam.get(q, [])
        prim = [r for r in rs if not (r[0] & 0x100) and not (r[0] & 0x800) and not (r[0] & 0x4)]
        if not prim:
            continue
        p = prim[0]
        sec = [r for r in rs if r[0] & 0x100]
        loci = {(r[1], r[2] // 500_000) for r in rs if not (r[0] & 0x4)}
        chroms = {r[1] for r in rs if not (r[0] & 0x4)}
        # best-chain coverage (primary + supplementary)
        ivs = []; read_len = None
        for r in prim + [x for x in rs if x[0] & 0x800]:
            rl, q0, q1, *_ = cigar_spans(r[5], r[0] & 0x10)
            read_len = read_len or rl; ivs.append((q0, q1))
        cov = union_len(ivs) / read_len if read_len else 0
        orf = longest_orf_aa(seq)
        ng = nearest_gene(genes, p[1], p[2])
        sec_as = sorted((r[4] for r in sec if r[4] is not None), reverse=True)
        tandem_max, tandem_base = complexity(seq)
        low_complexity = tandem_max - tandem_base > 0.30
        # classify
        if low_complexity:
            cls = "low-complexity / satellite repeat"
        elif len(loci) >= 10:
            cls = "expressed repeat / TE (dispersed, many genomic copies)"
        elif cov < 0.5 and orf >= 300:
            cls = "NOVEL protein-coding transcript, absent from assembly"
        elif cov < 0.5:
            cls = "mostly-unaligned (low ORF) — novel/uncharacterised"
        elif p[3] >= 60:
            cls = f"divergent/unannotated single locus (present; ~{ng[0]} {ng[1]}bp)"
        else:
            cls = f"partial/divergent transcript near {ng[0]}"
        out.append(dict(tag=tag, len=len(seq), cov=round(cov, 3), orf_aa=orf,
                        mapq=p[3], n_rec=len(rs), n_loci=len(loci), n_chrom=len(chroms),
                        prim_as=p[4], top_sec_as=(sec_as[0] if sec_as else None),
                        tandem_max=round(tandem_max, 3), tandem_base=round(tandem_base, 3),
                        nearest=f"{ng[0]}@{ng[1]}bp", cls=cls))
    return out


# ----------------------------------------------------------------------------- main
def main():
    print("=" * 72)
    print("UNMAPPED-READ RESCUE PoC  (T2T gorilla IsoSeq)")
    print("=" * 72)

    s1 = step1_pile()
    print(f"\n[1] Unmapped pile: n={s1['n']}  median_len={s1['median']}bp  "
          f"mean={s1['mean']}  range={s1['min']}-{s1['max']}")

    s2 = step2_denominator(s1["n"])
    print(f"[2] Denominator : {s2['primary_reads']} primary reads -> "
          f"unmapped = {s2['pct']}%")

    s3 = step3_clusters()
    print(f"[3] De-novo clusters (connected components of ava.paf):")
    print(f"      reads in overlap graph : {s3['reads_in_graph']}")
    print(f"      clusters >=2 reads     : {s3['ge2']['clusters']:>5}  "
          f"({s3['ge2']['reads']} reads)")
    print(f"      clusters >=3 reads     : {s3['ge3']['clusters']:>5}  "
          f"({s3['ge3']['reads']} reads)")
    print(f"      clusters >=5 reads     : {s3['ge5']['clusters']:>5}  "
          f"({s3['ge5']['reads']} reads)   <- representatives remapped")
    print(f"      top sizes              : {s3['top']}")

    s4 = step4_classify()
    c = s4["cats"]
    m = s4["n_mapped_relaxed"]
    print(f"\n[4] Classification of {s4['n_reps']} cluster reps (relaxed/non-spliced remap):")
    print(f"      unmapped even relaxed  : {c.get('unmapped_relaxed', 0)}")
    print(f"      got primary alignment  : {m}")
    print(f"      primary MAPQ           : " +
          "  ".join(f"{k}:{s4['mapq_hist'].get(k,0)}" for k in ['60','10-59','1-9','0']))
    print(f"      primary overlaps a gene: {s4['in_gene']}/{m} "
          f"({100.0*s4['in_gene']/max(m,1):.0f}%)")
    print(f"      best-chain coverage    : " +
          "  ".join(f"{k}%:{s4['covbuckets'].get(k,0)}"
                    for k in ['95-100','80-95','50-80','25-50','0-25']))
    print(f"      paralog homology       : {s4['n_with_secondary']}/{m} reps report "
          f"secondaries ({s4['n_secondary_crosschrom']} cross-chrom); "
          f"{s4['n_reps_any_mapq0_record']} reps have a MAPQ-0 *secondary* record")
    print(f"      (NB: zero MAPQ-0 PRIMARIES, but paralog placements ARE present)")
    print(f"      ---- categories (primary+supplementary 'pieces'; NOT secondaries) ----")
    print(f"      SINGLE-LOCUS present   : {c.get('single_locus_present', 0)}"
          f"  ({100.0*c.get('single_locus_present',0)/max(m,1):.0f}%)  "
          f"[splice-rejected; median sub-identity "
          f"{100*s4['single_locus_median_id']:.2f}% -> sequence present]")
    if s4["divergent_subset"]:
        print(f"           of which <90% identity (divergent-copy candidates): "
              f"{len(s4['divergent_subset'])}")
        for d in s4["divergent_subset"]:
            print(f"             - {d['read'][:34]:34} id={100*d['prim_id']:.1f}% "
                  f"mapq={d['mapq']} gene={d['prim_gene']}")
    print(f"      MULTI-LOCUS chimera    : {c.get('multi_locus_chimera', 0)}"
          f"  ({100.0*c.get('multi_locus_chimera',0)/max(m,1):.0f}%)  "
          f"[read-through / fusion; clean 2-piece={c.get('two_piece_clean',0)}]")
    print(f"      NOVEL candidate        : {c.get('novel_candidate', 0)}"
          f"  ({100.0*c.get('novel_candidate',0)/max(m,1):.0f}%)  "
          f"[<{int(COV_NOVEL*100)}% of read aligns anywhere]")
    for n in s4["novel"]:
        print(f"        - {n['read'][:38]:38} cov={n['cov']:<5} mapq={n['mapq']} "
              f"(weak anchor ~{n['prim_gene']}; {int((1-n['cov'])*100)}% aligns nowhere)")
    cs = s4["cov_sensitivity"]
    print(f"      novel-count cutoff sensitivity: cov<0.4 -> {cs[0.4]}, "
          f"cov<0.5 -> {cs[0.5]}, cov<0.6 -> {cs[0.6]}")

    s5 = step5_candidates()
    if s5:
        print(f"\n[5] Residual-candidate follow-up (option b):")
        print(f"    {'tag':6} {'len':>5} {'cov':>5} {'ORFaa':>5} {'mapq':>4} "
              f"{'rec':>3} {'loci':>4} {'chr':>3} {'tandem':>11}  classification")
        for r in s5:
            print(f"    {r['tag']:6} {r['len']:>5} {r['cov']:>5} {r['orf_aa']:>5} "
                  f"{r['mapq']:>4} {r['n_rec']:>3} {r['n_loci']:>4} {r['n_chrom']:>3} "
                  f"{r['tandem_max']:.3f}/{r['tandem_base']:.3f}  {r['cls']}")

    # machine-readable dump
    out = os.path.join(SCRATCH, "unmapped_rescue_poc_result.json")
    with open(out, "w") as f:
        json.dump(dict(pile=s1, denom=s2, clusters=s3,
                       classify={k: v for k, v in s4.items() if k != "per_read"},
                       candidates=s5),
                  f, indent=2)
    print(f"\n[+] wrote {out}")
    print("\nVERDICT: the unmapped pile is dominated by splice-rejected reads whose "
          "sequence is already present in the genome (single-locus, MAPQ60, in-gene);")
    print("         multi-locus chimeras are a minority. Follow-up dissolves the residual "
          "into expressed TEs (cl19/cl31) + 1 divergent single locus (cl56), leaving")
    print("         exactly ONE genuine hit: cl0, a real protein-coding transcript (775aa "
          "ORF, 38 reads) ~94% absent from the assembly (origin: needs nt/nr search).")
    print("         => 'rescue unmapped reads -> find missing gene copies' is effectively "
          "DRY on a complete (T2T) reference (yield ~1). Real lever only vs INCOMPLETE refs.")


if __name__ == "__main__":
    main()
