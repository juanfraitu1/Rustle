#!/usr/bin/env python3
"""Amendment 12 acceptance helpers (docs/PREREG_rna_allele_haplotype_count_2026-10-01.md): `o3_candidates` on Amendment 7's 53-family
held-out. `bench/rna_allele/accept_o3_candidates.sh` runs the subcommands in order; the stage and every minimap2 call run there (or, for
`keep`, here) under `tools/rlock.sh heavy`.

  nets     the stage's read nets recomputed from the BAM (pass A: primary / secondary records on a surviving copy; pass B: unmapped reads
           >= 300 bp given to a family by the 31-mer attribution of `FamilyKmerIndex`), the 1,000-read cap (`sample_net`: splitmix64
           Fisher-Yates, seed 1) and their D / S composition -> nets.tsv, net_reads.tsv (checked against the stage's families.tsv by `report`)
  plan     families -> batches of the stage (each batch one foreground call < 10 min), balanced on an estimated cost -> batches.txt
  concat   the batches' tables and FASTAs -> one set (`<prefix>.*`), one header, families in the stage's order
  contigs  the flagged candidates' unions renamed `iso_<family>_<k>` (the scorer counts only `iso_*` references as contigs) -> iso.contigs.fa,
           iso_names.tsv (iso -> cand)
  label    each contig's best UNMASKED hit (identity x coverage, `link_test.best_hits`) -> source D (the family's masked interval) / S:<copy> /
           elsewhere / none -> contigs.tsv (`link_test.py` layout, linked = 0), empty merge/paf/<family>.paf, links to R.bam, labels, panel
  keep     A12-2: each flagged candidate's reads (reads.tsv via clusters.tsv) against the unions and the cluster consensus sequences
           (minimap2 -c -x splice:hq -uf -N 10); kept = best AS on the read's own union >= 0.98 x its best AS over its candidate's
           cluster consensus sequences (`rep_choice.py`: no record on the union = not kept; no record on any of its consensus sequences =
           not measured). Pooled (all unions / all of clusters.fa as targets, the registered verdict) and isolated (per candidate: its union
           alone, its consensus sequences alone; the check)
  report   the stage's counts, the D / S labels, detection, candidates per family, clusters per candidate, the >= 2-cluster floor, the
           cause of each deleted copy without a D-derived flagged candidate -> report.json + stdout
  decompose  (post hoc, not a registered rule) arm M read by read for the stage and for IsoCon's Amendment 8 run (checked equal to
           merge_test.py score's totals), the deleted copies' reads split by where the R arm put them (on a surviving copy / mapped only
           elsewhere / unmapped) -> decompose.json + stdout

    accept_o3_candidates.py nets --w /mnt/linuxdisk/tmp/rna_allele/a12 --linktest /mnt/linuxdisk/tmp/rna_allele/linktest
"""
import argparse
import collections
import csv
import glob
import json
import os
import re
import subprocess

import pysam

K = 31
MIN_UNMAPPED_LEN = 300
ATTRIB_MAX_FAMILIES = 8
ATTRIB_MIN_FRAC, ATTRIB_MIN_RATIO = 0.30, 2.0
MAX_READS = 1000
TIE = 0.98
M64 = (1 << 64) - 1
GOLDEN = 0x9E3779B97F4A7C15
COMP = str.maketrans("ACGTNacgtn", "TGCANtgcan")
CODE = {"A": 0, "C": 1, "G": 2, "T": 3, "a": 0, "c": 1, "g": 2, "t": 3}
MM2_KEEP = ["-c", "-x", "splice:hq", "-uf", "-N", "10", "-t", "4"]


def name_key(s):
    """`name_order` of o3_candidates.rs: the text before the trailing digits, the number (digit count, digits), the whole name"""
    stem = s.rstrip("0123456789")
    digits = s[len(stem):].lstrip("0")
    return (stem, len(digits), digits, s)


def tsv(path):
    return list(csv.DictReader(open(path), delimiter="\t"))


def fasta(path):
    seqs, cur = {}, None
    for ln in open(path):
        if ln[0] == ">":
            cur = ln[1:].strip().split()[0]; seqs[cur] = []
        elif cur is not None:
            seqs[cur].append(ln.strip())
    return {k: "".join(v) for k, v in seqs.items()}


def oriented(rd):
    s = rd.query_sequence or ""
    return s.translate(COMP)[::-1] if rd.is_reverse else s


# ---- nets: the stage's pass A / pass B and the cap, recomputed ------------------------------------------------------------------------

def canonical_kmers(seq, k=K):
    mask, shift = (1 << (2 * k)) - 1, 2 * (k - 1)
    fw = rv = valid = 0
    out = []
    for ch in seq:
        c = CODE.get(ch)
        if c is None:
            fw = rv = valid = 0
            continue
        fw = ((fw << 2) | c) & mask
        rv = (rv >> 2) | ((3 - c) << shift)
        valid += 1
        if valid >= k:
            out.append(fw if fw < rv else rv)
    return out


def kmer_index(copies, fams):
    """FamilyKmerIndex::build: k-mer -> bitmask of the families carrying it; k-mers of > 8 families dropped"""
    fi = {f: i for i, f in enumerate(fams)}
    idx = {}
    for fam, seq in copies:
        bit = 1 << fi[fam]
        for km in set(canonical_kmers(seq)):
            idx[km] = idx.get(km, 0) | bit
    return {km: b for km, b in idx.items() if bin(b).count("1") <= ATTRIB_MAX_FAMILIES}


def attribute(idx, fams, read):
    """FamilyKmerIndex::attribute: the family with the most k-mer hits (multiplicity counted) when >= 30% of the read's k-mers hit it
    and it leads the runner-up >= 2x"""
    if len(read) < MIN_UNMAPPED_LEN:
        return None
    kms = canonical_kmers(read)
    hits = collections.Counter()
    for km in kms:
        b = idx.get(km, 0)
        while b:
            low = b & -b
            hits[low.bit_length() - 1] += 1
            b ^= low
    if not hits:
        return None
    v = sorted(hits.items(), key=lambda x: (-x[1], fams[x[0]]))
    best, n = v[0]
    second = v[1][1] if len(v) > 1 else 0
    if n < ATTRIB_MIN_FRAC * len(kms):
        return None
    if second > 0 and n < ATTRIB_MIN_RATIO * second:
        return None
    return fams[best]


def splitmix(x):
    z = (x + GOLDEN) & M64
    z = ((z ^ (z >> 30)) * 0xBF58476D1CE4E5B9) & M64
    z = ((z ^ (z >> 27)) * 0x94D049BB133111EB) & M64
    return z ^ (z >> 31)


def sample_net(names, cap=MAX_READS):
    v = sorted(set(names))
    if len(v) > cap:
        state = 1
        for i in range(len(v) - 1, 0, -1):
            r = splitmix(state); state = (state + GOLDEN) & M64
            j = r % (i + 1)
            v[i], v[j] = v[j], v[i]
        v = sorted(v[:cap])
    return v


def nets(a):
    rows = tsv(f"{a.w}/A12.copies.tsv")
    fams = list(dict.fromkeys(r["family_id"] for r in rows))
    cseq = {}
    cur = None
    for ln in open(f"{a.w}/A12.copies.fa"):
        if ln[0] == ">":
            h = ln[1:].strip().split("|"); cur = (h[0], int(h[1])); cseq[cur] = []
        else:
            cseq[cur].append(ln.strip())
    copies = [(r["family_id"], "".join(cseq[(r["family_id"], int(r["copy_idx"]))])) for r in rows]
    idx = kmer_index(copies, fams)
    lab = {r["read"]: r for r in tsv(f"{a.linktest}/labels.tsv")}
    bam = pysam.AlignmentFile(f"{a.linktest}/R.bam")
    names = {f: set() for f in fams}
    via = collections.defaultdict(dict)
    seqs = {}
    for r in rows:
        f, chrom = r["family_id"], r["chrom"]
        lo, hi = int(r["locus_start"]), int(r["locus_end"])
        for rd in bam.fetch(chrom, lo, max(hi, lo + 1)):
            if rd.is_unmapped or rd.is_supplementary:
                continue
            n = rd.query_name
            if not rd.is_secondary and n not in seqs:
                s = oriented(rd)
                if s:
                    seqs[n] = s
            names[f].add(n); via[f].setdefault(n, "A")
    need = {n for f in fams for n in names[f] if n not in seqs}
    n_unm = n_attr = 0
    for rd in pysam.AlignmentFile(f"{a.linktest}/R.bam").fetch(until_eof=True):
        if rd.is_unmapped:
            if rd.query_length < MIN_UNMAPPED_LEN:
                continue
            n_unm += 1
            s = oriented(rd)
            f = attribute(idx, fams, s)
            if f is None:
                continue
            n_attr += 1
            names[f].add(rd.query_name); via[f].setdefault(rd.query_name, "B")
            seqs.setdefault(rd.query_name, s)
            continue
        if not need or rd.is_secondary or rd.is_supplementary:
            continue
        if rd.query_name in need and rd.query_name not in seqs:
            s = oriented(rd)
            if s:
                seqs[rd.query_name] = s
    with open(f"{a.w}/nets.tsv", "w") as o, open(f"{a.w}/net_reads.tsv", "w") as nr:
        o.write("family\tn_net\tn_used\tD_net\tD_used\tS_net\tS_used\tother_net\tvia_B\tD_via_B\n")
        nr.write("family\tread\tused\trole\tvia\n")
        for f in fams:
            net = sorted(n for n in names[f] if n in seqs)
            used = set(sample_net(net))
            role = lambda n: lab[n]["role"] if n in lab and lab[n]["family"] == f else ("other" if n in lab else "unlabelled")
            c = collections.Counter((role(n), n in used) for n in net)
            vb = sum(1 for n in net if via[f][n] == "B")
            dvb = sum(1 for n in net if via[f][n] == "B" and role(n) == "D")
            o.write(f"{f}\t{len(net)}\t{len(used)}\t{c[('D', True)] + c[('D', False)]}\t{c[('D', True)]}\t{c[('S', True)] + c[('S', False)]}\t"
                    f"{c[('S', True)]}\t{sum(v for (r, _), v in c.items() if r not in ('D', 'S'))}\t{vb}\t{dvb}\n")
            for n in net:
                nr.write(f"{f}\t{n}\t{int(n in used)}\t{role(n)}\t{via[f][n]}\n")
    print(f"families {len(fams)}; k-mers indexed {len(idx)}; unmapped reads >= {MIN_UNMAPPED_LEN} bp {n_unm}, attributed {n_attr}; "
          f"nets total {sum(len([n for n in names[f] if n in seqs]) for f in fams)}")


def plan(a):
    """LPT over an estimated cost (s) = 25 x n_used / 1000 + 2 (calibrated on the registered run's first batch: 14-31 s of phase 1 per
    1,000-read family) of the families not already in the first --keep lines of batches.txt (batches that already ran stay as they are)"""
    rows = tsv(f"{a.w}/nets.tsv")
    bf = a.batch_file or f"{a.w}/batches.txt"
    kept = [ln.strip() for ln in open(bf)][: a.keep] if a.keep else []
    done = {f for ln in kept for f in ln.split(",")}
    cost = {r["family"]: 25 * int(r["n_used"]) / 1000 + 2 for r in rows if r["family"] not in done}
    bins = [[0.0, []] for _ in range(a.batches)]
    for f in sorted(cost, key=lambda f: (-cost[f], name_key(f))):
        b = min(bins, key=lambda b: b[0])
        b[0] += cost[f]; b[1].append(f)
    os.makedirs(os.path.dirname(os.path.abspath(bf)), exist_ok=True)
    with open(bf, "w") as o:
        for ln in kept:
            o.write(ln + "\n")
        for c, fs in bins:
            o.write(",".join(sorted(fs, key=name_key)) + "\n")
    for i, ln in enumerate(kept):
        print(f"batch {i}: {len(ln.split(','))} families (kept)")
    for i, (c, fs) in enumerate(bins):
        print(f"batch {len(kept) + i}: {len(fs)} families, estimated phase 1 {c:.0f} s")


# ---- concat / contigs / label -------------------------------------------------------------------------------------------------------

def concat(a):
    parts = sorted(glob.glob(f"{a.w}/{a.prefix}_g*.families.tsv"), key=lambda p: int(re.search(r"_g(\d+)\.families", p).group(1)))
    stems = [p[: -len(".families.tsv")] for p in parts]
    fam_of = {}
    seen = collections.Counter()
    for ext in ("families.tsv", "candidates.tsv", "clusters.tsv", "reads.tsv"):
        header, rows = None, []
        for s in stems:
            with open(f"{s}.{ext}") as fh:
                h = fh.readline()
                assert header is None or h == header, f"{s}.{ext}: header differs"
                header = h
                for ln in fh:
                    if ln.strip():
                        rows.append(ln)
        rows.sort(key=lambda ln: name_key(ln.split("\t", 2)[0 if ext != "reads.tsv" else 1]))
        with open(f"{a.w}/{a.prefix}.{ext}", "w") as o:
            o.write(header); o.writelines(rows)
        if ext == "families.tsv":
            seen.update(ln.split("\t", 1)[0] for ln in rows)
        if ext == "candidates.tsv":
            fam_of.update({ln.split("\t")[1]: ln.split("\t")[0] for ln in rows})
        if ext == "clusters.tsv":
            fam_of.update({ln.split("\t")[1]: ln.split("\t")[0] for ln in rows})
    dup = [f for f, n in seen.items() if n > 1]
    assert not dup, f"families in more than one batch: {dup}"
    for ext in ("contigs.fa", "clusters.fa"):
        recs = {}
        for s in stems:
            for k, v in fasta(f"{s}.{ext}").items():
                assert k not in recs, f"{ext}: {k} twice"
                recs[k] = v
        with open(f"{a.w}/{a.prefix}.{ext}", "w") as o:
            for k in sorted(recs, key=lambda k: (name_key(fam_of[k]), name_key(k))):
                o.write(f">{k}\n{recs[k]}\n")
    print(f"batches {len(stems)}; families {len(seen)}; candidates {sum(1 for k in fam_of if k.startswith('cand_'))}")


def contigs(a):
    cand = {r["candidate"]: r for r in tsv(f"{a.w}/{a.prefix}.candidates.tsv")}
    seqs = fasta(f"{a.w}/{a.prefix}.contigs.fa")
    with open(f"{a.w}/iso.contigs.fa", "w") as o, open(f"{a.w}/iso_names.tsv", "w") as t:
        t.write("iso\tcandidate\tfamily\tlength\tn_clusters\tn_reads\td\n")
        for c, s in seqs.items():
            r = cand[c]
            fam, k = r["family"], c.rsplit("_", 1)[1]
            assert c == f"cand_{fam}_{k}" and r["flagged"] == "1" and int(r["union_len"]) == len(s), c
            iso = f"iso_{fam}_{k}"
            o.write(f">{iso}\n{s}\n")
            t.write(f"{iso}\t{c}\t{fam}\t{len(s)}\t{r['n_clusters']}\t{r['n_reads']}\t{r['d']}\n")
    print(f"flagged contigs {len(seqs)} (candidates {len(cand)})")


def best_hits(path):
    """link_test.best_hits: per query, the hit with the highest identity x coverage"""
    b = {}
    for ln in open(path):
        f = ln.split("\t")
        s = int(f[9]) / max(1, int(f[10])) * (int(f[3]) - int(f[2])) / max(1, int(f[1]))
        if f[0] not in b or s > b[f[0]][0]:
            b[f[0]] = (s, f[5], int(f[7]), int(f[8]), int(f[9]), int(f[1]))
    return b


def label(a):
    P = {p["fam"]: p for p in json.load(open(f"{a.linktest}/panel.json"))}
    iso = tsv(f"{a.w}/iso_names.tsv")
    B = best_hits(f"{a.w}/iso.base.paf")
    os.makedirs(f"{a.w}/merge/paf", exist_ok=True)
    src_count = collections.Counter()
    with open(f"{a.w}/contigs.tsv", "w") as t:
        t.write("contig\tfamily\toutput\tlength\tbest_masked\td\tlinked\tsource\n")
        for r in iso:
            fam, b = r["family"], B.get(r["iso"])
            src = "none"
            if b:
                for lab_, (c, s0, e, g) in [("D", P[fam]["mask"])] + [("S:" + kk[3], kk) for kk in P[fam]["keep"]]:
                    if b[1] == c and b[2] < e and s0 < b[3]:
                        src = lab_; break
                else:
                    src = "elsewhere"
            src_count["S" if src.startswith("S:") else src] += 1
            t.write(f"{r['iso']}\t{fam}\t{r['candidate']}\t{r['length']}\tNA\t{r['d']}\t0\t{src}\n")
    for fam in sorted({r["family"] for r in iso}):
        open(f"{a.w}/merge/paf/{fam}.paf", "w").close()
    for f in ("R.bam", "R.bam.bai", "labels.tsv", "panel.json"):
        dst = f"{a.w}/{f}"
        if not os.path.lexists(dst):
            os.symlink(f"{a.linktest}/{f}", dst)
    print(f"contigs {len(iso)} by source: {dict(src_count)}; families with contigs {len({r['family'] for r in iso})}")


# ---- keep (A12-2) -------------------------------------------------------------------------------------------------------------------

def max_as(paf):
    """(query, target) -> best AS over the PAF lines"""
    best = {}
    for ln in open(paf):
        f = ln.rstrip("\n").split("\t")
        s = next((int(x[5:]) for x in f[12:] if x.startswith("AS:i:")), None)
        if s is None:
            continue
        key = (f[0], f[5])
        if s > best.get(key, -10 ** 9):
            best[key] = s
    return best


def mm2(target, query, out):
    with open(out, "w") as o, open(out + ".log", "w") as e:
        subprocess.run(["minimap2", *MM2_KEEP, target, query], stdout=o, stderr=e, check=True)


def keep(a):
    cand = {r["candidate"]: r for r in tsv(f"{a.w}/{a.prefix}.candidates.tsv") if r["flagged"] == "1"}
    cl_of = collections.defaultdict(list)                       # candidate -> its clusters
    for r in tsv(f"{a.w}/{a.prefix}.clusters.tsv"):
        if r["candidate"] in cand:
            cl_of[r["candidate"]].append(r["cluster"])
    cand_of_cl = {c: k for k, cs in cl_of.items() for c in cs}
    pairs = collections.defaultdict(set)                        # candidate -> reads
    for r in tsv(f"{a.w}/{a.prefix}.reads.tsv"):
        k = cand_of_cl.get(r["cluster"])
        if k:
            pairs[k].add(r["read"])
    reads_needed = set().union(*pairs.values()) if pairs else set()
    rseq = {}
    for p in sorted(glob.glob(f"{a.linktest}/scored.part*.fa")):
        for n, s in fasta(p).items():
            if n in reads_needed:
                rseq[n] = s
    missing = reads_needed - set(rseq)
    assert not missing, f"{len(missing)} reads without a sequence, e.g. {sorted(missing)[:3]}"
    unions = fasta(f"{a.w}/{a.prefix}.contigs.fa")
    clseq = fasta(f"{a.w}/{a.prefix}.clusters.fa")
    kd = f"{a.w}/keep"
    os.makedirs(f"{kd}/iso", exist_ok=True)
    with open(f"{kd}/reads.fa", "w") as o:
        for n in sorted(reads_needed):
            o.write(f">{n}\n{rseq[n]}\n")
    # pooled: every union / every cluster consensus as the targets (the registered verdict)
    mm2(f"{a.w}/{a.prefix}.contigs.fa", f"{kd}/reads.fa", f"{kd}/pooled.union.paf")
    mm2(f"{a.w}/{a.prefix}.clusters.fa", f"{kd}/reads.fa", f"{kd}/pooled.clusters.paf")
    pooled = (max_as(f"{kd}/pooled.union.paf"), max_as(f"{kd}/pooled.clusters.paf"))
    # isolated: per candidate, its union alone and its own consensus sequences alone (the check)
    iso_u, iso_c = {}, {}
    for k in sorted(pairs, key=name_key):
        d = f"{kd}/iso/{k}"
        with open(f"{d}.reads.fa", "w") as o:
            for n in sorted(pairs[k]):
                o.write(f">{n}\n{rseq[n]}\n")
        with open(f"{d}.union.fa", "w") as o:
            o.write(f">{k}\n{unions[k]}\n")
        with open(f"{d}.clusters.fa", "w") as o:
            for c in cl_of[k]:
                o.write(f">{c}\n{clseq[c]}\n")
        mm2(f"{d}.union.fa", f"{d}.reads.fa", f"{d}.union.paf")
        mm2(f"{d}.clusters.fa", f"{d}.reads.fa", f"{d}.clusters.paf")
        iso_u.update(max_as(f"{d}.union.paf")); iso_c.update(max_as(f"{d}.clusters.paf"))
    res = {}
    for name, (U, C) in (("pooled", pooled), ("isolated", (iso_u, iso_c))):
        per = {}
        tot = collections.Counter()
        for k in sorted(pairs, key=name_key):
            c = collections.Counter()
            for n in pairs[k]:
                best = max((C[(n, cl)] for cl in cl_of[k] if (n, cl) in C), default=None)
                if best is None:
                    c["no_record_on_clusters"] += 1; continue
                u = U.get((n, k))
                c["kept" if u is not None and u >= TIE * best else "no_record_on_union" if u is None else "lost"] += 1
            per[k] = dict(c)
            tot.update(c)
        measured = tot["kept"] + tot["lost"] + tot["no_record_on_union"]
        frac = tot["kept"] / max(1, measured)
        fr = sorted(v.get("kept", 0) / max(1, v.get("kept", 0) + v.get("lost", 0) + v.get("no_record_on_union", 0)) for v in per.values())
        below = [k for k, v in per.items() if v.get("kept", 0) < 0.95 * (v.get("kept", 0) + v.get("lost", 0) + v.get("no_record_on_union", 0))]
        res[name] = dict(total=dict(tot), measured=measured, kept_frac=frac, per_candidate=per, below95=sorted(below, key=name_key),
                         dist=dict(min=fr[0] if fr else None, p10=fr[len(fr) // 10] if fr else None, median=fr[len(fr) // 2] if fr else None))
        print(f"[{name}] (read, candidate) pairs {sum(tot.values())}: kept {tot['kept']}, lost {tot['lost']}, no record on the union "
              f"{tot['no_record_on_union']}, not measured (no record on its consensus sequences) {tot['no_record_on_clusters']} -> kept "
              f"{frac:.2%} of {measured} measured (bar 95%: {'PASSES' if frac >= 0.95 else 'FAILS'}); per candidate kept fraction min "
              f"{res[name]['dist']['min']}, p10 {res[name]['dist']['p10']}, median {res[name]['dist']['median']}; candidates below 95%: "
              f"{len(below)}/{len(per)} {res[name]['below95']}")
    json.dump(res, open(f"{a.w}/keep.json", "w"), indent=0)


# ---- report -------------------------------------------------------------------------------------------------------------------------

def report(a):
    P = {p["fam"]: p for p in json.load(open(f"{a.linktest}/panel.json"))}
    lab = {r["read"]: r for r in tsv(f"{a.linktest}/labels.tsv")}
    fams = {r["family"]: r for r in tsv(f"{a.w}/{a.prefix}.families.tsv")}
    cands = tsv(f"{a.w}/{a.prefix}.candidates.tsv")
    clus = tsv(f"{a.w}/{a.prefix}.clusters.tsv")
    reads = tsv(f"{a.w}/{a.prefix}.reads.tsv")
    ctg = {r["output"]: r for r in tsv(f"{a.w}/contigs.tsv")}       # cand id -> contigs.tsv row (source)
    nets_rows = {r["family"]: r for r in tsv(f"{a.nets}/nets.tsv")}
    netr = collections.defaultdict(dict)
    for r in tsv(f"{a.nets}/net_reads.tsv"):
        netr[r["family"]][r["read"]] = r
    out = {}
    # 1. the replicated nets against the stage's own counts
    mism = [f for f in fams if (int(fams[f]["n_net"]), int(fams[f]["n_used"])) != (int(nets_rows[f]["n_net"]), int(nets_rows[f]["n_used"]))]
    out["net_check"] = dict(families=len(fams), mismatches=mism)
    print(f"net check: {len(fams) - len(mism)}/{len(fams)} families with n_net and n_used equal to the stage's" + (f"; MISMATCH {mism}" if mism else ""))
    # 2. stage counts
    tot = collections.Counter()
    for r in fams.values():
        for k in ("n_net", "n_used", "n_clusters", "n_in_reference", "n_linked", "n_new", "n_candidates", "n_flagged"):
            tot[k] += int(r[k])
    out["stage"] = dict(tot)
    print("stage totals:", dict(tot), f"| families with >= 1 flagged candidate {sum(1 for r in fams.values() if int(r['n_flagged']) > 0)}/{len(fams)}")
    # 3. flagged candidates by label; detection
    flagged = [c for c in cands if c["flagged"] == "1"]
    src = {c["candidate"]: ctg[c["candidate"]]["source"] for c in flagged}
    kind = lambda s: "S" if s.startswith("S:") else s
    by = collections.Counter(kind(s) for s in src.values())
    dfam = collections.Counter(c["family"] for c in flagged if src[c["candidate"]] == "D")
    sfam = {c["family"] for c in flagged if kind(src[c["candidate"]]) != "D"}
    out["flagged"] = dict(n=len(flagged), by_source=dict(by), families_with_D=len(dfam), families_exactly_one_D=sum(1 for v in dfam.values() if v == 1),
                          families_with_non_D=len(sfam))
    print(f"flagged candidates {len(flagged)} by source {dict(by)}; deleted copies with >= 1 D-derived flagged candidate {len(dfam)}/{len(P)}; "
          f"exactly one D-derived candidate in {out['flagged']['families_exactly_one_D']}/{len(dfam)}; D-derived candidates per such family "
          f"{dict(sorted(collections.Counter(dfam.values()).items()))}; families with a non-D flagged candidate {len(sfam)}")
    # 4. candidates per family, clusters per candidate
    cpf = collections.Counter(int(r["n_candidates"]) for r in fams.values())
    fpf = collections.Counter(int(r["n_flagged"]) for r in fams.values())
    cpc = collections.Counter(int(c["n_clusters"]) for c in cands)
    cpc_f = collections.Counter(int(c["n_clusters"]) for c in flagged)
    rpc = sorted(int(c["n_reads"]) for c in flagged)
    out["per_family"] = dict(candidates=dict(sorted(cpf.items())), flagged=dict(sorted(fpf.items())))
    out["per_candidate"] = dict(clusters_all=dict(sorted(cpc.items())), clusters_flagged=dict(sorted(cpc_f.items())),
                                reads_flagged_median=rpc[len(rpc) // 2] if rpc else None)
    print(f"candidates per family {dict(sorted(cpf.items()))}; flagged per family {dict(sorted(fpf.items()))}; clusters per candidate (all) "
          f"{dict(sorted(cpc.items()))}, (flagged) {dict(sorted(cpc_f.items()))}; reads per flagged candidate median "
          f"{out['per_candidate']['reads_flagged_median']}, min {rpc[0] if rpc else None}, max {rpc[-1] if rpc else None}")
    # 5. the >= 2-cluster floor (every cluster holds >= --min-cluster 3 reads, so >= 2 clusters implies >= 6 reads: a subset of the flagged)
    two = [c for c in cands if int(c["n_clusters"]) >= 2]
    assert all(c["flagged"] == "1" for c in two), "a >= 2-cluster candidate below the 6-read floor"
    by2 = collections.Counter(kind(src[c["candidate"]]) for c in two)
    d2 = {c["family"] for c in two if src[c["candidate"]] == "D"}
    out["floor2"] = dict(n=len(two), by_source=dict(by2), families_with_D=len(d2))
    print(f">= 2-cluster floor: {len(two)} candidates (of {len(flagged)} flagged at >= 6 reads) by source {dict(by2)}; deleted copies with a "
          f"D-derived one {len(d2)}/{len(P)}")
    # 6. causes for the deleted copies without a D-derived flagged candidate
    cl = {r["cluster"]: r for r in clus}
    candrow = {c["candidate"]: c for c in cands}
    keep_iv = lambda fam: [(k[0], k[1], k[2]) for k in P[fam]["keep"]]
    causes, table = collections.Counter(), []
    for fam in sorted(P, key=name_key):
        if fam in dfam:
            continue
        nr = netr.get(fam, {})
        d_net = [n for n, r in nr.items() if r["role"] == "D"]
        d_used = [n for n in d_net if nr[n]["used"] == "1"]
        fate = collections.Counter()
        linked_to = collections.Counter()
        in_cl = {r["read"]: r["cluster"] for r in reads if r["family"] == fam}
        for n in d_used:
            c = in_cl.get(n)
            if c is None:
                fate["no_cluster"] += 1; continue
            row = cl[c]
            if row["fate"] == "linked":
                loc = row["linked_to"]
                ch, span = loc.rsplit(":", 1); s, e = map(int, span.split("-"))
                on = any(ch == kc and s < ke and ks < e for kc, ks, ke in keep_iv(fam))
                fate["linked"] += 1; linked_to["survivor" if on else "elsewhere"] += 1
            else:
                k = candrow[row["candidate"]]
                fate["flagged_not_D" if k["flagged"] == "1" else "below_support"] += 1
        if not d_net:
            cause = "i: no D read in the net"
        elif not d_used:
            cause = "i': D reads in the net, none under the 1,000-read cap"
        else:
            order = ["flagged_not_D", "below_support", "linked", "no_cluster"]
            top = max(order, key=lambda k: (fate[k], -order.index(k)))
            cause = {"flagged_not_D": "iv: D reads in a flagged candidate labelled S / elsewhere",
                     "below_support": "ii: D reads' component below --min-support",
                     "linked": "iii: D reads' clusters linked (allele of a reference locus)",
                     "no_cluster": "ii: D reads in no cluster reaching the chain (< --min-cluster, split off, or an in-reference cluster)"}[top]
        causes[cause] += 1
        f = fams[fam]
        table.append(dict(family=fam, D_net=len(d_net), D_used=len(d_used), n_net=int(f["n_net"]), n_used=int(f["n_used"]), clusters=int(f["n_clusters"]),
                          in_ref=int(f["n_in_reference"]), linked=int(f["n_linked"]), new=int(f["n_new"]), candidates=int(f["n_candidates"]),
                          flagged=int(f["n_flagged"]), D_fates=dict(fate), D_linked_to=dict(linked_to), cause=cause))
    out["causes"] = dict(counts=dict(causes), table=table)
    print(f"deleted copies without a D-derived flagged candidate: {len(table)}")
    for c, n in sorted(causes.items()):
        print(f"  {n:3d}  {c}")
    print("family\tD_net\tD_used\tn_net\tn_used\tclusters\tin_ref\tlinked\tnew\tcands\tflagged\tD reads by fate\tcause")
    for t in table:
        print(f"{t['family']}\t{t['D_net']}\t{t['D_used']}\t{t['n_net']}\t{t['n_used']}\t{t['clusters']}\t{t['in_ref']}\t{t['linked']}\t{t['new']}\t"
              f"{t['candidates']}\t{t['flagged']}\t{t['D_fates']} {t['D_linked_to'] or ''}\t{t['cause']}")
    # 7. D reads in the nets overall, and where the R arm (masked genome) puts them: unmapped / a record on a surviving copy of the family
    # (= pass A) / mapped only elsewhere
    dn = sum(int(r["D_net"]) for r in nets_rows.values()); du = sum(int(r["D_used"]) for r in nets_rows.values())
    nD = sum(1 for r in lab.values() if r["role"] == "D")
    rec = collections.defaultdict(list)
    for rd in pysam.AlignmentFile(f"{a.linktest}/R.bam").fetch(until_eof=True):
        r = lab.get(rd.query_name)
        if rd.is_supplementary or not r or r["role"] != "D":
            continue
        rec[rd.query_name].append(None if rd.is_unmapped else (rd.reference_name, rd.reference_start, rd.reference_end,
                                                                rd.get_tag("de") if rd.has_tag("de") else None))
    status = collections.defaultdict(collections.Counter)
    de_surv = collections.defaultdict(list)                     # family -> per D read, the smallest de over its records on a surviving copy
    for n, rs in rec.items():
        fam = lab[n]["family"]
        rs = [x for x in rs if x]
        on = [x[3] for x in rs if x[3] is not None and any(x[0] == kc and x[1] < ke and ks < x[2] for kc, ks, ke in keep_iv(fam))]
        st = "unmapped" if not rs else "on_survivor" if any(x[0] == kc and x[1] < ke and ks < x[2] for x in rs for kc, ks, ke in keep_iv(fam)) else "elsewhere"
        status[fam][st] += 1
        if on:
            de_surv[fam].append(min(on))
    med = lambda v: sorted(v)[len(v) // 2] if v else None
    for t in table:
        t["D_de_on_survivor_median"] = med(de_surv[t["family"]])
    print("cause-table families with D reads on a surviving copy, median de of those reads (the D-to-survivor divergence; delta = 0.00958):")
    print("  " + "; ".join(f"{t['family']} {t['D_de_on_survivor_median']:.4f} (n {len(de_surv[t['family']])}, {t['cause'][:3].strip()})"
                         for t in table if t["D_de_on_survivor_median"] is not None))
    allst = collections.Counter()
    for v in status.values():
        allst.update(v)
    c1 = [t["family"] for t in table if t["cause"].startswith("i:")]
    st1 = collections.Counter()
    for f in c1:
        st1.update(status[f])
    out["D_reads"] = dict(total=nD, in_nets=dn, used=du, R_status=dict(allst), R_status_cause_i=dict(st1),
                          per_family={f: dict(status[f]) for f in sorted(status, key=name_key)})
    print(f"D reads: {nD} scored; in a net of their family {dn}; under the cap {du}; R arm: {dict(allst)}; in the {len(c1)} cause-i families: "
          f"{dict(st1)} (families whose D reads are all unmapped: {sum(1 for f in c1 if set(status[f]) == {'unmapped'})})")
    json.dump(out, open(f"{a.w}/report.json", "w"), indent=0)


# ---- decompose (post hoc): the read-level result of arm M by where the R arm put each read, stage vs IsoCon -----------------------------

def arm_m(w, lab, P):
    """merge_test.score's arm M, read by read: (role, class) per scored read; contigs.tsv (linked = 0) + merge/components.tsv as loci"""
    rows = [r for r in tsv(f"{w}/contigs.tsv") if r["linked"] == "0"]
    fam_of = {r["contig"]: r["family"] for r in rows}
    srcs = {r["contig"]: r["source"] for r in rows}
    mp = {r["contig"]: r["component"] for r in tsv(f"{w}/merge/components.tsv")}
    assert set(mp) == set(fam_of), f"{w}: components.tsv and contigs.tsv disagree"
    holds, lfam = collections.defaultdict(set), {}
    for c, L in mp.items():
        holds[L].add(srcs[c]); lfam[L] = fam_of[c]
    recs = collections.defaultdict(list)
    for rd in pysam.AlignmentFile(f"{w}/RIL.bam").fetch(until_eof=True):
        if rd.query_name not in lab or rd.is_supplementary:
            continue
        recs[rd.query_name].append(None if rd.is_unmapped else (rd.is_secondary, rd.reference_name, rd.reference_start, rd.reference_end,
                                                                 rd.get_tag("AS") if rd.has_tag("AS") else 0))

    def locus(chrom, s, e, fam):
        if chrom.startswith("iso_"):
            return ("ctg", mp[chrom])
        for kk in P[fam]["keep"]:
            if chrom == kk[0] and s < kk[2] and kk[1] < e:
                return ("copy", kk[3])
        return ("other", f"{chrom}:{s // 100000}")
    out, fm_to = {}, collections.Counter()
    for n, r in lab.items():
        rs = [x for x in recs.get(n, []) if x]
        fam = r["family"]
        call = None
        if rs:
            prim = next((x for x in rs if not x[0]), rs[0])
            srt = sorted(rs, key=lambda x: -x[4])
            if len(srt) > 1 and srt[1][4] > 0 and srt[1][4] >= TIE * srt[0][4] and \
                    len({locus(x[1], x[2], x[3], fam) for x in srt if x[4] >= TIE * srt[0][4]}) > 1:
                call = None
            else:
                call = locus(prim[1], prim[2], prim[3], fam)
        if call is None:
            out[n] = (r["role"], "unplaced"); continue
        kind, key = call
        own = kind == "ctg" and lfam[key] == fam
        if r["role"] == "D":
            out[n] = ("D", "right" if own and "D" in holds[key] else "wrong")
        elif kind == "copy" and key == r["copy"]:
            out[n] = ("S", "stay")
        elif kind == "ctg":
            out[n] = ("S", "stay" if own and ("S:" + r["copy"]) in holds[key] else "false_move")
            if out[n][1] == "false_move":                     # where the false move went: the target locus' sources, own family or not
                tgt = "+".join(sorted("S:other" if x.startswith("S:") else x for x in holds[key]))
                fm_to[(fam, "own family" if own else "other family", tgt)] += 1
        else:
            out[n] = ("S", "elsewhere")
    out["__false_moves__"] = fm_to
    return out


def decompose(a):
    P = {p["fam"]: p for p in json.load(open(f"{a.linktest}/panel.json"))}
    lab = {r["read"]: r for r in tsv(f"{a.linktest}/labels.tsv")}
    rec = collections.defaultdict(list)
    for rd in pysam.AlignmentFile(f"{a.linktest}/R.bam").fetch(until_eof=True):
        if rd.query_name in lab and not rd.is_supplementary:
            rec[rd.query_name].append(None if rd.is_unmapped else (rd.reference_name, rd.reference_start, rd.reference_end))
    status = {}
    for n, r in lab.items():
        rs = [x for x in rec.get(n, []) if x]
        keep_iv = [(k[0], k[1], k[2]) for k in P[r["family"]]["keep"]]
        status[n] = "unmapped" if not rs else "on_survivor" if any(c == kc and s < ke and ks < e for c, s, e in rs for kc, ks, ke in keep_iv) else "elsewhere"
    res = {}
    for name, w in (("IsoCon (Amendment 8)", a.linktest), ("o3_candidates", a.w)):
        calls = arm_m(w, lab, P)
        fm_to = calls.pop("__false_moves__")
        tot = collections.Counter(calls.values())
        ref = json.load(open(f"{w}/merge/score.json"))["M"]
        mism = {k: (v, ref.get("|".join(k), 0)) for k, v in tot.items() if ref.get("|".join(k), 0) != v}
        assert not mism, f"{name}: the replication differs from merge_test.py score: {mism}"
        x = collections.Counter((status[n], c) for n, (role, c) in calls.items() if role == "D")
        per = collections.defaultdict(collections.Counter)
        for n, (role, c) in calls.items():
            if role == "D":
                per[lab[n]["family"]][c] += 1
        res[name] = dict(totals={"|".join(k): v for k, v in tot.items()}, D_by_R_status={"|".join(k): v for k, v in x.items()},
                         D_right_per_family={f: per[f]["right"] for f in sorted(per, key=name_key)},
                         false_moves_to={"|".join(k): v for k, v in sorted(fm_to.items(), key=lambda kv: -kv[1])})
        print(f"[{name}] arm M equals merge_test.py score (D right {tot[('D', 'right')]}, S false moves {tot[('S', 'false_move')]}); D reads by "
              f"their R-arm status:")
        for st in ("on_survivor", "elsewhere", "unmapped"):
            n_st = sum(v for (s_, _), v in x.items() if s_ == st)
            print(f"    {st:12s} {n_st:6d}: right {x[(st, 'right')]:6d}  wrong {x[(st, 'wrong')]:6d}  unplaced {x[(st, 'unplaced')]:6d}")
        print("    S false moves by (family, target's family, target's sources): " +
              "; ".join(f"{f} {o} {t}: {v}" for (f, o, t), v in sorted(fm_to.items(), key=lambda kv: -kv[1])))
    i, o = res["IsoCon (Amendment 8)"]["D_right_per_family"], res["o3_candidates"]["D_right_per_family"]
    print("per family D right (IsoCon -> o3_candidates), families where either > 0:")
    print("  " + "; ".join(f"{f} {i.get(f, 0)}->{o.get(f, 0)}" for f in sorted(set(i) | set(o), key=name_key) if i.get(f, 0) or o.get(f, 0)))
    json.dump(res, open(f"{a.w}/decompose.json", "w"), indent=0)


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("cmd", choices=["nets", "plan", "concat", "contigs", "label", "keep", "report", "decompose"])
    ap.add_argument("--w", required=True, help="the run's work dir")
    ap.add_argument("--linktest", default="/mnt/linuxdisk/tmp/rna_allele/linktest", help="Amendment 7's work dir")
    ap.add_argument("--prefix", default="cand")
    ap.add_argument("--batches", type=int, default=7, help="plan: the number of NEW batches")
    ap.add_argument("--keep", type=int, default=0, help="plan: keep the first N lines of the batch file (batches already run)")
    ap.add_argument("--batch-file", default=None, help="plan: the batch file (default <w>/batches.txt)")
    ap.add_argument("--nets", default=None, help="dir with nets.tsv / net_reads.tsv (default --w)")
    a = ap.parse_args(argv)
    a.nets = a.nets or a.w
    globals()[a.cmd](a)


if __name__ == "__main__":
    main()
