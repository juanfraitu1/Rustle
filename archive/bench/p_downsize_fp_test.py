#!/usr/bin/env python3
"""
p_downsize_fp_test.py — is minimap2's permissive -p 0.1 (used to build GGO_mm.bam) a source of
FALSE POSITIVES (spurious multi-mappings -> spurious read-conflict / de-tie edges), and would
"downsizing" to the -p 0.8 default remove noise — or does the shipped de-tie gate already do -p's
job (better, on `de` not raw score), making -p 0.1 a safe wide net?

CLEAN TEST (unlike -G): -p ONLY filters secondary alignments; it never changes a primary or
re-aligns. So -p 0.8 is faithfully SIMULATED by post-filtering GGO_mm.bam: keep primary + only
secondaries whose score ratio to the global primary >= 0.8. Proxy = AS_sec/AS_prim (DP score);
the FAITHFUL predictor is s1_sec/s1_prim (chaining score = exactly what -p thresholds). We report
BOTH and VALIDATE against a real minimap2 -p 0.8 re-map of a deterministic read sample.

Prediction under test: the shipped gate (family_def_genomewide: DELTA=0.005, DE_MAX=0.05) requires
a read to fit BOTH loci COMPARABLY (small |de_i-de_j|, both de<=DE_MAX). A low-AS-ratio secondary
has high `de`, so the gate already rejects exactly the secondaries -p 0.8 would drop => -p 0.8 cleans
the RAW conflict graph but barely moves the GATED output. Risk on the other side: -p 0.8 can drop
GENUINE ties to divergent copies (~70-80% identity), the hard collapsed / K-frontier cases.

Parts:
  A. genome-wide AS- and s1-ratio distribution of secondaries (deterministic 1/32 crc32 subsample):
     % dropped at -p 0.8 and -p 0.5 (the noise volume).
  B. proxy validation: deterministic ~200-read sample -> minimap2 -p 0.1 AND -p 0.8 re-map;
     how well AS-ratio>=0.8 (and s1-ratio>=0.8) post-filter reproduces the real -p 0.8 secondary set,
     PLUS a primary-invariance check (the whole simulation rests on -p never touching primaries).
  C. co-located families (validated_families.tsv): raw multi-mapping degree + de-tie edges under
     -p 0.1 vs the -p 0.8 filter, using the FAITHFUL s1-ratio (chain score) for the drop test with an
     AS-ratio fallback; the pure AS-ratio run is kept for a faithfulness delta. Reports the GATE-vs-p
     correspondence (of secondaries -p 0.8 drops, gate-rejected junk vs genuine gate-kept ties), the
     divergence 2x2 (where truly divergent de>DE_MAX multimappers land), and the family concentration
     of the genuine edge loss.

Deterministic: PYTHONHASHSEED=0, crc32 sampling, fixed seeds, sorted before every write.
Run: PYTHONHASHSEED=0 /home/juanfra/miniforge3/bin/python bench/p_downsize_fp_test.py
"""
import collections
import csv
import json
import os
import re
import subprocess
import sys
import zlib
from bisect import bisect_right

# ---------------------------------------------------------------- config
BAM = "/home/juanfra/winloci_scratch/GGO_mm.bam"
MMI = "/home/juanfra/winloci_scratch/GGO.splice.mmi"
VALIDATED = "/home/juanfra/winloci_scratch/validated_families.tsv"
MINIMAP2 = "/home/juanfra/miniforge3/bin/minimap2"
SAMTOOLS = "samtools"
HERE = os.path.dirname(os.path.abspath(__file__))
SCRATCH = ("/tmp/claude-1000/-mnt-c-Users-jfris-Desktop/"
           "03f0156d-63b6-4d0f-bf09-862bddee491c/scratchpad")
OUT_TSV = os.path.join(HERE, "p_downsize_fp_test.tsv")
OUT_JSON = os.path.join(HERE, "p_downsize_fp_test.json")

# shipped gate (family_def_genomewide.py)
DELTA = 0.005
DE_MAX = 0.05
MIN_READS = 3
# -p thresholds to test
P_THRESH = [0.8, 0.5]
P_MAIN = 0.8
# genome-wide subsample: keep read iff crc32(qname) & MASK == 0  (1/32)
SUB_MASK = 0x1F
# proxy-validation sample size
N_SAMPLE = 200
POS_TOL = 100  # bp tolerance when matching two alignments to the "same locus"
CIG = re.compile(r"(\d+)([MIDNSHP=X])")


def log(*a):
    print(*a, flush=True)


def ref_span(cigar):
    s = 0
    for n, op in CIG.findall(cigar):
        if op in "MDN=X":
            s += int(n)
    return s


def tags(fields):
    """AS(int), s1(int), de(float) from optional SAM fields; None if absent."""
    a = s1 = de = None
    for t in fields:
        if t.startswith("AS:i:"):
            a = int(t[5:])
        elif t.startswith("s1:i:"):
            s1 = int(t[5:])
        elif t.startswith("de:f:"):
            de = float(t[5:])
    return a, s1, de


# ---------------------------------------------------------------- families
def load_families():
    """family -> list of copies (locus, chrom, start, end); per-chrom sorted copy index."""
    fam_copies = collections.defaultdict(list)
    with open(VALIDATED) as f:
        for row in csv.DictReader(f, delimiter="\t"):
            fam_copies[row["family"]].append(
                (row["locus"], row["chrom"], int(row["start"]), int(row["end"])))
    by_chrom = collections.defaultdict(list)  # chrom -> sorted [(start,end,fam,locus)]
    for fam, cps in fam_copies.items():
        for locus, c, s, e in cps:
            by_chrom[c].append((s, e, fam, locus))
    for c in by_chrom:
        by_chrom[c].sort()
    return fam_copies, by_chrom


def best_copy(by_chrom, chrom, rs, re_):
    """(fam, locus) of the family copy whose interval maximally overlaps [rs,re_), else None."""
    cand = by_chrom.get(chrom)
    if not cand:
        return None
    starts = [g[0] for g in cand]
    hi = bisect_right(starts, re_)
    best, best_ov = None, 0
    for k in range(hi - 1, -1, -1):
        s, e, fam, locus = cand[k]
        if e <= rs:
            if rs - e > 5_000_000:
                break
            continue
        ov = min(e, re_) - max(s, rs)
        if ov > best_ov:
            best_ov, best = ov, (fam, locus)
    return best


# ---------------------------------------------------------------- phase 1: region scan
def scan_family_regions(by_chrom):
    """One multi-region samtools pass over all family copies.
    Returns read -> {locus: placement} where placement = dict(fam,de,AS,s1,is_primary,chrom,pos).
    Keeps, per (read,locus), the record with the best AS."""
    bed = os.path.join(SCRATCH, "fam_regions.bed")
    ivs = []
    for c, lst in by_chrom.items():
        for s, e, fam, locus in lst:
            ivs.append((c, s, e))
    ivs.sort()
    with open(bed, "w") as o:
        for c, s, e in ivs:
            o.write(f"{c}\t{s}\t{e}\n")
    reads = collections.defaultdict(dict)
    p = subprocess.Popen([SAMTOOLS, "view", "-@", "4", "-M", "-L", bed, BAM],
                         stdout=subprocess.PIPE, text=True, bufsize=1 << 20)
    n = 0
    for line in p.stdout:
        f = line.rstrip("\n").split("\t")
        if len(f) < 11:
            continue
        flag = int(f[1])
        if flag & 0x804:  # unmapped or supplementary
            continue
        n += 1
        qname, chrom, pos, cigar = f[0], f[2], int(f[3]) - 1, f[5]
        hit = best_copy(by_chrom, chrom, pos, pos + ref_span(cigar))
        if hit is None:
            continue
        fam, locus = hit
        AS, s1, de = tags(f[11:])
        is_prim = (flag & 0x100) == 0
        d = reads[qname]
        prev = d.get(locus)
        if prev is None or (AS is not None and (prev["AS"] is None or AS > prev["AS"])):
            d[locus] = dict(fam=fam, de=de, AS=AS, s1=s1, is_primary=is_prim,
                            chrom=chrom, pos=pos)
        elif is_prim:
            d[locus]["is_primary"] = True
    p.wait()
    log(f"[phase1] family-region records parsed: {n:,}; reads touching families: {len(reads):,}")
    return reads


# ---------------------------------------------------------------- phase 2: proxy validation
comp_tab = str.maketrans("ACGTNacgtn", "TGCANtgcan")


def revcomp(s):
    return s.translate(comp_tab)[::-1]


def pick_sample(reads):
    """Deterministic ~N_SAMPLE reads that (a) have an in-family PRIMARY (SEQ available, co-located
    read) and (b) multimap onto >=2 family copies. Evenly spaced over sorted qnames."""
    elig = []
    for q, plc in reads.items():
        if len(plc) < 2:
            continue
        if any(v["is_primary"] for v in plc.values()):
            elig.append(q)
    elig.sort()
    if not elig:
        return []
    if len(elig) <= N_SAMPLE:
        return elig
    step = len(elig) / N_SAMPLE
    idx = sorted({int(i * step) for i in range(N_SAMPLE)})
    return [elig[i] for i in idx]


def extract_sample_fasta(sample, reads):
    """Pull each sample read's SEQ from its in-family primary record (restore orientation)."""
    want = {}  # qname -> (chrom, region) to fetch primary
    for q in sample:
        prim_loc = None
        for locus, v in reads[q].items():
            if v["is_primary"]:
                prim_loc = v
                break
        if prim_loc:
            want[q] = prim_loc
    fa = os.path.join(SCRATCH, "p_sample.fa")
    got = {}
    # group fetch by chrom window to limit samtools calls
    by_region = collections.defaultdict(list)
    for q, v in want.items():
        by_region[(v["chrom"], v["pos"] // 1_000_000)].append(q)
    qset = set(want)
    for (chrom, mb), qs in sorted(by_region.items()):
        reg = f"{chrom}:{mb*1_000_000+1}-{(mb+1)*1_000_000}"
        p = subprocess.run([SAMTOOLS, "view", "-F", "0x900", BAM, reg],
                           capture_output=True, text=True)
        for line in p.stdout.splitlines():
            i = line.find("\t")
            qn = line[:i]
            if qn not in qset or qn in got:
                continue
            f = line.split("\t")
            flag = int(f[1])
            seq = f[9]
            if seq == "*" or not seq:
                continue
            got[qn] = revcomp(seq) if (flag & 16) else seq
    with open(fa, "w") as o:
        for q in sorted(got):
            o.write(f">{q}\n{got[q]}\n")
    log(f"[phase2] sample reads eligible written to FASTA: {len(got)}")
    return fa, sorted(got)


def remap(fa, p_val):
    """minimap2 with the ORIGINAL params + given -p. Returns read -> list of alignments
    (is_secondary, chrom, pos, AS, s1, de, cigar, mapq). CIGAR/MAPQ are kept so we can
    prove -p never touches the PRIMARY alignment (only filters secondaries)."""
    cmd = [MINIMAP2, "-ax", "splice:hq", "-uf", "--eqx", "-Y", "-N50",
           "-p", str(p_val), "-t", "3", MMI, fa]
    p = subprocess.run(cmd, capture_output=True, text=True)
    out = collections.defaultdict(list)
    for line in p.stdout.splitlines():
        if line.startswith("@"):
            continue
        f = line.split("\t")
        if len(f) < 11:
            continue
        flag = int(f[1])
        if flag & 0x804:  # unmapped/supplementary
            continue
        AS, s1, de = tags(f[11:])
        out[f[0]].append(dict(sec=bool(flag & 0x100), chrom=f[2], pos=int(f[3]) - 1,
                              AS=AS, s1=s1, de=de, cigar=f[5], mapq=int(f[4])))
    return out


def loci_set(aligns, only_secondary=True):
    r = set()
    for a in aligns:
        if only_secondary and not a["sec"]:
            continue
        r.add((a["chrom"], a["pos"]))
    return r


def match_count(setA, setB):
    """count loci in setA that have a partner within POS_TOL in setB (same chrom)."""
    byc = collections.defaultdict(list)
    for c, p in setB:
        byc[c].append(p)
    for v in byc.values():
        v.sort()
    m = 0
    for c, p in setA:
        arr = byc.get(c)
        if not arr:
            continue
        # linear (small sets)
        if any(abs(p - q) <= POS_TOL for q in arr):
            m += 1
    return m


def proxy_filter(aligns, use_s1, thr):
    """Simulate -p `thr` on a -p0.1 remap: keep primary + secondaries with score-ratio>=thr.
    Returns the kept-SECONDARY loci set."""
    prim = None
    for a in aligns:
        if not a["sec"]:
            prim = a
            break
    if prim is None:
        return set()
    base = prim["s1"] if use_s1 else prim["AS"]
    if not base:
        return set()
    kept = set()
    for a in aligns:
        if not a["sec"]:
            continue
        sc = a["s1"] if use_s1 else a["AS"]
        if sc is None:
            continue
        if sc / base >= thr:
            kept.add((a["chrom"], a["pos"]))
    return kept


def primary_invariance(sample, r01, r08):
    """Load-bearing check: -p ONLY filters secondaries, it never re-aligns or changes a PRIMARY.
    For each sample read compare the primary alignment between the -p0.1 and -p0.8 remaps."""
    def primary_of(aligns):
        for a in aligns or []:
            if not a["sec"]:
                return a
        return None
    inv = dict(reads=0, both_mapped=0,
               identical_chrom_pos_as_cigar=0, mapq_diff_only=0,
               pos_or_as_or_cigar_diff=0, locus_moved=0,
               missing_in_p01=0, missing_in_p08=0,
               secondary_promoted_to_primary=0)
    for q in sample:
        inv["reads"] += 1
        p1 = primary_of(r01.get(q))
        p8 = primary_of(r08.get(q))
        if p1 is None and p8 is not None:
            inv["missing_in_p01"] += 1
            continue
        if p8 is None and p1 is not None:
            inv["missing_in_p08"] += 1
            continue
        if p1 is None and p8 is None:
            continue
        inv["both_mapped"] += 1
        same_locus = (p1["chrom"] == p8["chrom"] and abs(p1["pos"] - p8["pos"]) <= POS_TOL)
        if not same_locus:
            inv["locus_moved"] += 1
            # was the p0.8 primary a secondary at p0.1 (a genuine promotion)?
            if any(a["sec"] and a["chrom"] == p8["chrom"]
                   and abs(a["pos"] - p8["pos"]) <= POS_TOL for a in r01.get(q, [])):
                inv["secondary_promoted_to_primary"] += 1
            continue
        core_same = (p1["chrom"] == p8["chrom"] and p1["pos"] == p8["pos"]
                     and p1["AS"] == p8["AS"] and p1["cigar"] == p8["cigar"])
        if core_same and p1["mapq"] == p8["mapq"]:
            inv["identical_chrom_pos_as_cigar"] += 1
        elif core_same:
            inv["mapq_diff_only"] += 1
            inv["identical_chrom_pos_as_cigar"] += 1  # pos/AS/CIGAR still identical
        else:
            inv["pos_or_as_or_cigar_diff"] += 1
    return inv


def validate_proxy(fa, sample):
    """Compare AS-ratio>=0.8 and s1-ratio>=0.8 post-filters (on -p0.1 remap) vs real -p0.8 remap.
    Also verify the primary alignment is invariant under -p (the whole simulation rests on this)."""
    log("[phase2] re-mapping sample at -p 0.1 (reproduces original) and -p 0.8 (real drop)...")
    r01 = remap(fa, 0.1)
    r08 = remap(fa, 0.8)
    inv = primary_invariance(sample, r01, r08)
    log(f"[phase2] primary-invariance: {json.dumps(inv)}")
    agg = dict(reads=0, real_sec=0,
               as_recall_num=0, as_prec_den=0, as_prec_num=0,
               s1_recall_num=0, s1_prec_den=0, s1_prec_num=0,
               exact_as=0, exact_s1=0,
               p01_vs_orig_reads=0)
    per = []
    for q in sample:
        a01 = r01.get(q)
        a08 = r08.get(q)
        if not a01 or not a08:
            continue
        real = loci_set(a08, only_secondary=True)
        as_keep = proxy_filter(a01, use_s1=False, thr=P_MAIN)
        s1_keep = proxy_filter(a01, use_s1=True, thr=P_MAIN)
        agg["reads"] += 1
        agg["real_sec"] += len(real)
        # recall: real loci matched by proxy; precision: proxy loci matched by real
        as_rec = match_count(real, as_keep)
        as_prc = match_count(as_keep, real)
        s1_rec = match_count(real, s1_keep)
        s1_prc = match_count(s1_keep, real)
        agg["as_recall_num"] += as_rec
        agg["as_prec_num"] += as_prc
        agg["as_prec_den"] += len(as_keep)
        agg["s1_recall_num"] += s1_rec
        agg["s1_prec_num"] += s1_prc
        agg["s1_prec_den"] += len(s1_keep)
        if as_rec == len(real) and as_prc == len(as_keep) and len(real) == len(as_keep):
            agg["exact_as"] += 1
        if s1_rec == len(real) and s1_prc == len(s1_keep) and len(real) == len(s1_keep):
            agg["exact_s1"] += 1
        per.append((q, len(real), len(as_keep), len(s1_keep), as_rec, s1_rec))
    n = max(agg["reads"], 1)
    rd = agg["real_sec"] or 1
    res = dict(
        sample_reads=agg["reads"],
        real_p08_secondaries=agg["real_sec"],
        as_ratio_recall=round(agg["as_recall_num"] / rd, 4),
        as_ratio_precision=round(agg["as_prec_num"] / max(agg["as_prec_den"], 1), 4),
        s1_ratio_recall=round(agg["s1_recall_num"] / rd, 4),
        s1_ratio_precision=round(agg["s1_prec_num"] / max(agg["s1_prec_den"], 1), 4),
        exact_setmatch_as=round(agg["exact_as"] / n, 4),
        exact_setmatch_s1=round(agg["exact_s1"] / n, 4),
        primary_invariance=inv,
    )
    return res, per


# ---------------------------------------------------------------- phase 3: whole-BAM pass
def whole_bam_pass(famset):
    """One streaming pass over primary+secondary (exclude supplementary).
    (A) genome-wide 1/32 crc32 subsample: read -> {'prim':(AS,s1), 'sec':[(AS,s1,de)]}.
    (B) FAMSET global primary AS/s1 (exact -p reference score for family reads)."""
    sub = collections.defaultdict(lambda: {"prim": None, "sec": []})
    famprim = {}  # qname -> (AS, s1)  from the global primary record
    p = subprocess.Popen([SAMTOOLS, "view", "-@", "4", "-F", "0x800", BAM],
                         stdout=subprocess.PIPE, text=True, bufsize=1 << 20)
    n = 0
    for line in p.stdout:
        i = line.find("\t")
        qname = line[:i]
        in_sub = (zlib.crc32(qname.encode()) & SUB_MASK) == 0
        in_fam = qname in famset
        if not (in_sub or in_fam):
            continue
        f = line.rstrip("\n").split("\t")
        flag = int(f[1])
        if flag & 0x4:
            continue
        AS, s1, de = tags(f[11:])
        is_sec = bool(flag & 0x100)
        if in_sub:
            d = sub[qname]
            if is_sec:
                d["sec"].append((AS, s1, de))
            else:
                d["prim"] = (AS, s1)
        if in_fam and not is_sec:
            # keep best-AS primary (there should be exactly one)
            cur = famprim.get(qname)
            if cur is None or (AS is not None and (cur[0] is None or AS > cur[0])):
                famprim[qname] = (AS, s1)
        n += 1
    p.wait()
    log(f"[phase3] whole-BAM kept records: {n:,}; subsample reads: {len(sub):,}; "
        f"famset primaries found: {len(famprim):,}")
    return sub, famprim


def subsample_distribution(sub):
    """AS-ratio (and s1-ratio) of secondaries genome-wide; % dropped at each -p threshold."""
    as_ratios, s1_ratios = [], []
    n_sec = 0
    for q, d in sub.items():
        prim = d["prim"]
        if not prim or not prim[0]:
            continue
        pa, ps1 = prim
        for AS, s1, de in d["sec"]:
            n_sec += 1
            if AS is not None and pa:
                as_ratios.append(AS / pa)
            if s1 is not None and ps1:
                s1_ratios.append(s1 / ps1)
    as_ratios.sort()
    s1_ratios.sort()

    def frac_below(sorted_vals, thr):
        return round(sum(1 for v in sorted_vals if v < thr) / max(len(sorted_vals), 1), 4)

    def hist(vals):
        edges = [0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.0001]
        h = [0] * (len(edges) - 1)
        for v in vals:
            for k in range(len(edges) - 1):
                if edges[k] <= v < edges[k + 1]:
                    h[k] += 1
                    break
        return {f"[{edges[k]:.1f},{edges[k+1]:.1f})": h[k] for k in range(len(h))}

    return dict(
        secondaries_sampled=n_sec,
        as_ratio_available=len(as_ratios),
        s1_ratio_available=len(s1_ratios),
        pct_dropped_p08_as=frac_below(as_ratios, 0.8),
        pct_dropped_p05_as=frac_below(as_ratios, 0.5),
        pct_dropped_p08_s1=frac_below(s1_ratios, 0.8),
        pct_dropped_p05_s1=frac_below(s1_ratios, 0.5),
        as_ratio_hist=hist(as_ratios),
        s1_ratio_hist=hist(s1_ratios),
    )


# ---------------------------------------------------------------- phase 4: family copy structure
# de-histogram edges shared by the per-category divergence breakdown
DE_EDGES = [0.0, 0.005, 0.01, 0.02, 0.05, 0.1, 0.2, 0.3, 0.5, 1.0001]


def _de_bucket(de):
    for k in range(len(DE_EDGES) - 1):
        if DE_EDGES[k] <= de < DE_EDGES[k + 1]:
            return k
    return len(DE_EDGES) - 2


def family_structure(fam_copies, reads, famprim, use_s1=True):
    """Per family: degree + de-tie edges under -p0.1 vs the -p0.8 filter, and the gate-vs-p
    correspondence. The -p0.8 drop test uses the FAITHFUL chain-score ratio s1_sec/s1_prim
    (use_s1=True; the exact quantity minimap2 -p thresholds) with a per-placement fallback to
    AS_sec/AS_prim only when s1 is absent. use_s1=False reproduces the pure AS-ratio proxy so
    the two can be compared. Also splits secondaries by divergence to show where the truly
    divergent (de>DE_MAX) multimappers land in the drop/gate 2x2."""
    rows = []
    tot = dict(
        fam=0, reads_mm=0,
        deg_p01_sum=0, deg_p08_sum=0,
        edges_p01=0, edges_p08=0, edges_gate=0,
        edges_p08_lost=0, edges_gate_kept_p08_lost=0,
        sec=0, p08_drop=0, gate_reject=0,
        drop_and_reject=0, drop_but_gatekept=0, reject_but_p08kept=0, both_keep=0,
        ext_primary_reads=0,
        # divergent (de>DE_MAX) secondaries: where do they fall in the 2x2?
        div_sec=0, div_drop_and_reject=0, div_drop_but_gatekept=0,
        div_reject_but_p08kept=0, div_both_keep=0,
        nearid_sec=0,  # de<=DE_MAX secondaries
        sec_de_hist=[0] * (len(DE_EDGES) - 1),
        genuine_loss_de_hist=[0] * (len(DE_EDGES) - 1),
    )
    # index reads by family
    fam_reads = collections.defaultdict(list)  # fam -> list of (qname, {locus:placement})
    for q, plc in reads.items():
        fams = {v["fam"] for v in plc.values()}
        for fam in fams:
            sub = {loc: v for loc, v in plc.items() if v["fam"] == fam}
            if len(sub) >= 1:
                fam_reads[fam].append((q, sub))

    for fam in sorted(fam_copies, key=lambda x: (len(x), x)):
        rlist = fam_reads.get(fam, [])
        mm = [(q, s) for q, s in rlist if len(s) >= 2]
        if not mm:
            continue
        tot["fam"] += 1
        # per-read structures
        deg01 = deg08 = 0
        pair_p01 = collections.Counter()   # (loc_i,loc_j) -> reads supporting
        pair_p08 = collections.Counter()
        pair_gate = collections.Counter()
        f_sec = f_drop = f_reject = f_dr = f_dg = f_rp = f_bk = 0
        f_ext = 0
        for q, sub in mm:
            gp = famprim.get(q)
            as_prim_global = gp[0] if (gp and gp[0]) else None
            s1_prim_global = gp[1] if (gp and gp[1]) else None
            # local anchor = best-AS in-family placement (fallback if global primary AS missing)
            has_local_prim = any(v["is_primary"] for v in sub.values())
            if not has_local_prim:
                f_ext += 1
            anchor_loc = None
            best_as = -1
            for loc, v in sub.items():
                if v["is_primary"]:
                    anchor_loc = loc
                    break
            if anchor_loc is None:
                for loc, v in sub.items():
                    if v["AS"] is not None and v["AS"] > best_as:
                        best_as, anchor_loc = v["AS"], loc
            ref_as = as_prim_global
            if ref_as is None:  # fallback: best in-family AS
                ref_as = max((v["AS"] for v in sub.values() if v["AS"] is not None), default=None)
            ref_s1 = s1_prim_global
            if ref_s1 is None:
                ref_s1 = max((v["s1"] for v in sub.values() if v["s1"] is not None), default=None)
            anchor_de = sub[anchor_loc]["de"] if anchor_loc else None

            def score_ratio(v):
                """Faithful -p ratio: chain score s1_sec/s1_prim; fall back to AS only if s1 gone."""
                if use_s1 and v["s1"] is not None and ref_s1:
                    return v["s1"] / ref_s1
                if v["AS"] is not None and ref_as:
                    return v["AS"] / ref_as
                return None

            # degree p0.1 = all distinct copies; p0.8 = primary + score-ratio>=0.8
            kept08 = set()
            for loc, v in sub.items():
                is_anchor = (loc == anchor_loc)
                if is_anchor:
                    keep = True
                else:
                    r = score_ratio(v)
                    keep = (r is not None and r >= P_MAIN)
                if keep:
                    kept08.add(loc)
                # per-secondary (non-anchor) classification
                if not is_anchor:
                    f_sec += 1
                    de_v = v["de"]
                    if de_v is not None:
                        tot["sec_de_hist"][_de_bucket(de_v)] += 1
                    is_div = (de_v is not None and de_v > DE_MAX)
                    if is_div:
                        tot["div_sec"] += 1
                    elif de_v is not None:
                        tot["nearid_sec"] += 1
                    dropped = not keep
                    if dropped:
                        f_drop += 1
                    # gate: valid de-tie of this copy with anchor copy?
                    gate_ok = (anchor_de is not None and de_v is not None
                               and abs(de_v - anchor_de) <= DELTA
                               and max(de_v, anchor_de) <= DE_MAX)
                    rejected = not gate_ok
                    if rejected:
                        f_reject += 1
                    if dropped and rejected:
                        f_dr += 1
                        if is_div:
                            tot["div_drop_and_reject"] += 1
                    elif dropped and not rejected:
                        f_dg += 1   # GENUINE tie the gate keeps but -p0.8 drops
                        if de_v is not None:
                            tot["genuine_loss_de_hist"][_de_bucket(de_v)] += 1
                        if is_div:
                            tot["div_drop_but_gatekept"] += 1  # impossible by gate def (de<=DE_MAX)
                    elif rejected and not dropped:
                        f_rp += 1   # junk -p0.8 keeps but gate kills
                        if is_div:
                            tot["div_reject_but_p08kept"] += 1
                    else:
                        f_bk += 1   # both keep
                        if is_div:
                            tot["div_both_keep"] += 1
            deg01 += len(sub)
            deg08 += len(kept08)
            # pairwise edges
            locs = sorted(sub)
            for a in range(len(locs)):
                for b in range(a + 1, len(locs)):
                    la, lb = locs[a], locs[b]
                    va, vb = sub[la], sub[lb]
                    pair_p01[(la, lb)] += 1
                    if la in kept08 and lb in kept08:
                        pair_p08[(la, lb)] += 1
                    if (va["de"] is not None and vb["de"] is not None
                            and abs(va["de"] - vb["de"]) <= DELTA
                            and max(va["de"], vb["de"]) <= DE_MAX):
                        pair_gate[(la, lb)] += 1
        e01 = {k for k, c in pair_p01.items() if c >= MIN_READS}
        e08 = {k for k, c in pair_p08.items() if c >= MIN_READS}
        eg = {k for k, c in pair_gate.items() if c >= MIN_READS}
        lost = e01 - e08                       # edges -p0.8 removes vs raw
        genuine_lost = eg & lost               # gate-KEPT edges that -p0.8 removes = genuine loss
        rows.append(dict(
            family=fam, n_copies=len(fam_copies[fam]), reads_mm=len(mm),
            deg_p01_mean=round(deg01 / len(mm), 3), deg_p08_mean=round(deg08 / len(mm), 3),
            edges_p01=len(e01), edges_p08=len(e08), edges_gate=len(eg),
            edges_lost_p08=len(lost), edges_gatekept_lost_p08=len(genuine_lost),
            sec=f_sec, p08_drop=f_drop, gate_reject=f_reject,
            drop_and_reject=f_dr, drop_gatekept=f_dg, reject_p08kept=f_rp, both_keep=f_bk,
            ext_primary_reads=f_ext,
        ))
        tot["reads_mm"] += len(mm)
        tot["deg_p01_sum"] += deg01
        tot["deg_p08_sum"] += deg08
        tot["edges_p01"] += len(e01)
        tot["edges_p08"] += len(e08)
        tot["edges_gate"] += len(eg)
        tot["edges_p08_lost"] += len(lost)
        tot["edges_gate_kept_p08_lost"] += len(genuine_lost)
        tot["sec"] += f_sec
        tot["p08_drop"] += f_drop
        tot["gate_reject"] += f_reject
        tot["drop_and_reject"] += f_dr
        tot["drop_but_gatekept"] += f_dg
        tot["reject_but_p08kept"] += f_rp
        tot["both_keep"] += f_bk
        tot["ext_primary_reads"] += f_ext
    rows.sort(key=lambda r: (-r["reads_mm"], r["family"]))
    return rows, tot


# ---------------------------------------------------------------- main
def main():
    os.makedirs(SCRATCH, exist_ok=True)
    fam_copies, by_chrom = load_families()
    log(f"[load] families: {len(fam_copies)}; copies: {sum(len(v) for v in fam_copies.values())}")

    reads = scan_family_regions(by_chrom)
    famset = set(reads)

    sample = pick_sample(reads)
    log(f"[phase2] sample multimapping co-located reads: {len(sample)}")
    fa, sample_ok = extract_sample_fasta(sample, reads)
    proxy_res, proxy_per = validate_proxy(fa, sample_ok)
    log(f"[phase2] proxy validation: {json.dumps(proxy_res)}")

    sub, famprim = whole_bam_pass(famset)
    dist = subsample_distribution(sub)
    log(f"[phaseA] genome-wide noise: dropped@p0.8(AS)={dist['pct_dropped_p08_as']*100:.1f}% "
        f"dropped@p0.5(AS)={dist['pct_dropped_p05_as']*100:.1f}% "
        f"(s1: {dist['pct_dropped_p08_s1']*100:.1f}% / {dist['pct_dropped_p05_s1']*100:.1f}%)")

    # HEADLINE = faithful chain-score s1-ratio (the exact quantity minimap2 -p thresholds).
    rows, tot = family_structure(fam_copies, reads, famprim, use_s1=True)
    # COMPARISON = pure AS-ratio proxy (the weaker DP-score proxy) to quantify faithfulness delta.
    _, tot_as = family_structure(fam_copies, reads, famprim, use_s1=False)

    # ------- derived correspondence summary (on the faithful s1-ratio filter)
    sec = max(tot["sec"], 1)
    drop = max(tot["p08_drop"], 1)
    de_hist = {f"[{DE_EDGES[k]:.3f},{DE_EDGES[k+1]:.3f})": tot["sec_de_hist"][k]
               for k in range(len(tot["sec_de_hist"]))}
    genuine_de_hist = {f"[{DE_EDGES[k]:.3f},{DE_EDGES[k+1]:.3f})": tot["genuine_loss_de_hist"][k]
                       for k in range(len(tot["genuine_loss_de_hist"]))}
    corr = dict(
        score_used_for_p08_drop="s1_ratio (chain score; AS-ratio fallback if s1 absent)",
        families_with_multimappers=tot["fam"],
        multimap_reads=tot["reads_mm"],
        mean_degree_p01=round(tot["deg_p01_sum"] / max(tot["reads_mm"], 1), 3),
        mean_degree_p08=round(tot["deg_p08_sum"] / max(tot["reads_mm"], 1), 3),
        edges_p01=tot["edges_p01"], edges_p08=tot["edges_p08"], edges_gate=tot["edges_gate"],
        edges_lost_by_p08=tot["edges_p08_lost"],
        edges_gatekept_but_lost_by_p08=tot["edges_gate_kept_p08_lost"],
        frac_gate_edges_lost_by_p08=round(
            tot["edges_gate_kept_p08_lost"] / max(tot["edges_gate"], 1), 4),
        secondary_placements=tot["sec"],
        p08_dropped=tot["p08_drop"],
        gate_rejected=tot["gate_reject"],
        both_keep=tot["both_keep"],
        p08drop_AND_gatereject=tot["drop_and_reject"],
        p08drop_BUT_gatekept_GENUINE_LOSS=tot["drop_but_gatekept"],
        gatereject_BUT_p08kept_JUNK_GATE_CATCHES=tot["reject_but_p08kept"],
        frac_p08drops_that_are_gate_junk=round(tot["drop_and_reject"] / drop, 4),
        frac_p08drops_that_are_genuine_gatekept_ties=round(tot["drop_but_gatekept"] / drop, 4),
        frac_secondaries_dropped_by_p08=round(tot["p08_drop"] / sec, 4),
        frac_secondaries_rejected_by_gate=round(tot["gate_reject"] / sec, 4),
        ext_primary_reads=tot["ext_primary_reads"],
        secondary_de_hist=de_hist,
        genuine_loss_de_hist=genuine_de_hist,
    )

    # ------- where the truly DIVERGENT (de>DE_MAX) multimappers land in the drop/gate 2x2
    divs = max(tot["div_sec"], 1)
    divergent = dict(
        divergent_secondaries=tot["div_sec"],
        near_identical_secondaries=tot["nearid_sec"],
        div_killed_by_both=tot["div_drop_and_reject"],
        div_gate_kept_p08_dropped_GENUINE_LOSS=tot["div_drop_but_gatekept"],
        div_p08_kept_gate_rejected=tot["div_reject_but_p08kept"],
        div_both_keep=tot["div_both_keep"],
        frac_divergent_killed_by_gate=round(
            (tot["div_drop_and_reject"] + tot["div_reject_but_p08kept"]) / divs, 4),
        frac_divergent_killed_by_p08=round(
            (tot["div_drop_and_reject"] + tot["div_drop_but_gatekept"]) / divs, 4),
        note=("genuine-loss (p08drop & gate-kept) is <=DE_MAX by construction, so div_gate_kept "
              "is 0 by definition: truly divergent multimappers are rejected by the GATE itself, "
              "not uniquely by -p."),
    )

    # ------- concentration of the genuine edge loss across families
    loss_rows = sorted((r for r in rows if r["edges_gatekept_lost_p08"] > 0),
                       key=lambda r: (-r["edges_gatekept_lost_p08"], r["family"]))
    tot_gl = sum(r["edges_gatekept_lost_p08"] for r in rows) or 1
    concentration = dict(
        families_with_genuine_edge_loss=len(loss_rows),
        total_genuine_edges_lost=tot["edges_gate_kept_p08_lost"],
        top_family=(loss_rows[0]["family"] if loss_rows else None),
        top_family_genuine_edges_lost=(loss_rows[0]["edges_gatekept_lost_p08"] if loss_rows else 0),
        top_family_share=round((loss_rows[0]["edges_gatekept_lost_p08"] / tot_gl), 4)
        if loss_rows else 0.0,
        top3_share=round(sum(r["edges_gatekept_lost_p08"] for r in loss_rows[:3]) / tot_gl, 4)
        if loss_rows else 0.0,
    )

    # ------- faithfulness of the AS-ratio proxy vs the s1-ratio headline (aggregate delta)
    def _pct(a, b):
        return round(abs(a - b) / max(b, 1), 4)
    faithfulness = dict(
        p08_dropped_s1=tot["p08_drop"], p08_dropped_as=tot_as["p08_drop"],
        p08_dropped_rel_delta=_pct(tot_as["p08_drop"], tot["p08_drop"]),
        edges_lost_s1=tot["edges_p08_lost"], edges_lost_as=tot_as["edges_p08_lost"],
        genuine_edge_loss_s1=tot["edges_gate_kept_p08_lost"],
        genuine_edge_loss_as=tot_as["edges_gate_kept_p08_lost"],
        genuine_loss_placements_s1=tot["drop_but_gatekept"],
        genuine_loss_placements_as=tot_as["drop_but_gatekept"],
        note=("s1-ratio is the faithful chain-score proxy (sample recall 0.984/precision 1.00); "
              "AS-ratio is the weaker DP proxy (0.936/0.958). Headline uses s1; delta is small."),
    )

    log(f"[phaseC/s1] {json.dumps(corr)}")
    log(f"[phaseC/faithfulness] {json.dumps(faithfulness)}")
    log(f"[phaseC/divergent] {json.dumps(divergent)}")
    log(f"[phaseC/concentration] {json.dumps(concentration)}")

    result = dict(
        params=dict(BAM=BAM, DELTA=DELTA, DE_MAX=DE_MAX, MIN_READS=MIN_READS,
                    P_MAIN=P_MAIN, SUB_MASK=SUB_MASK, POS_TOL=POS_TOL, N_SAMPLE=N_SAMPLE),
        A_genomewide_noise=dist,
        B_proxy_validation=proxy_res,
        C_family_correspondence=corr,
        C_divergent_multimappers=divergent,
        C_genuine_loss_concentration=concentration,
        C_AS_vs_s1_faithfulness=faithfulness,
    )
    with open(OUT_JSON, "w") as o:
        json.dump(result, o, indent=2, sort_keys=True)

    cols = ["family", "n_copies", "reads_mm", "deg_p01_mean", "deg_p08_mean",
            "edges_p01", "edges_p08", "edges_gate", "edges_lost_p08",
            "edges_gatekept_lost_p08", "sec", "p08_drop", "gate_reject",
            "drop_and_reject", "drop_gatekept", "reject_p08kept", "both_keep",
            "ext_primary_reads"]
    with open(OUT_TSV, "w") as o:
        o.write("\t".join(cols) + "\n")
        for r in rows:
            o.write("\t".join(str(r[c]) for c in cols) + "\n")

    log(f"[write] {OUT_TSV}")
    log(f"[write] {OUT_JSON}")
    log("[done]")


if __name__ == "__main__":
    main()
