#!/usr/bin/env python3
"""The `o3_candidates` acceptance helpers on Amendment 7's 53-family held-out (docs/archive/2026-10/PREREG_rna_allele_haplotype_count_2026-10-01.md:
Amendment 12, and Amendment 13 + 13b-13e). `bench/rna_allele/accept_o3_candidates.sh` runs the subcommands in order; the stage and every
minimap2 call run there (or, for `nets` and `keep`, here) under `tools/rlock.sh heavy`.

  nets     (A13) the nets as the stage built them, replicated and checked against the stage's own products: pass A from the BAM (primary /
           secondary records on a surviving copy), pass B re-run per batch exactly as the stage runs it (Amendment 13b/13c: unmapped reads
           and the poorly placed reads in no net of THIS batch (R18), >= 300 bp, aligned with map-ont to the batch's net reads + every
           copy; a read joins the family of its best hit iff it covers >= 50% of the read and de <= 0.20), the 1,000-read cap
           (`sample_net`). Checked: the per-class counts of every batch log's pass-B line, n_net / n_used of families.tsv, every reads.tsv
           read under the cap of its family, nets.fa = the whole nets of the flagged families -> nets.tsv, net_reads.tsv, attrib.json
           (+ stdout: the attribution counts, right family by labels.tsv for every joined read and for the reads.tsv subset, the reads in
           nets of two batches). (A12's k-mer replica of the retired FamilyKmerIndex is at commit fde90c0a.)
  plan     families -> batches of the stage (each batch one foreground call < 10 min), balanced on an estimated cost -> batches.txt
  concat   the batches' tables and FASTAs -> one set (`<prefix>.*`), one header, families in the stage's order
  contigs  the flagged candidates' unions renamed `iso_<family>_<k>` (the scorer counts only `iso_*` references as contigs) -> iso.contigs.fa,
           iso_names.tsv (iso -> cand)
  label    each contig's best UNMASKED hit (identity x coverage, `link_test.best_hits`) -> source D (the family's masked interval) / S:<copy> /
           elsewhere / none -> contigs.tsv (`link_test.py` layout, linked = 0), empty merge/paf/<family>.paf, links to R.bam, labels, panel
  keep     A12-2 / A13-2: each flagged candidate's reads (reads.tsv via clusters.tsv) against the unions and the cluster consensus sequences
           (minimap2 -c -x splice:hq -uf -N 10); kept = best AS on the read's own union >= 0.98 x its best AS over its candidate's
           cluster consensus sequences (`rep_choice.py`: no record on the union = not kept; no record on any of its consensus sequences =
           not measured). Pooled (all unions / all of clusters.fa as targets, the registered verdict) and isolated (per candidate: its union
           alone, its consensus sequences alone; the check)
  report   the stage's counts, the stage log's per-family counters (A13: templates by kind, kept sets re-templated, absorptions undone),
           the D / S labels, detection, candidates per family, clusters per candidate, the >= 2-cluster floor, the cause of each deleted copy
           without a D-derived flagged candidate -> report.json + stdout
  comparator  (A13) Amendment 13b's comparator C = IsoCon's right D reads (Amendment 8's arm M, read by read) over the truth-free
           attainable D reads: (1) a record (primary or secondary) on a surviving copy of their family in R.bam, (2) attributed into their
           own family's net by this run (`nets`), else; A13-1 = stage D right >= 0.80 x C and false moves <= 5% of S reads (A12-1's
           10,230 beside); also C with part 2 read from reads.tsv alone (the reads that reached a reported cluster) -> comparator.json
  decompose  (post hoc, not a registered rule) arm M read by read for the stage and for IsoCon's Amendment 8 run (checked equal to
           merge_test.py score's totals), the deleted copies' reads split by where the R arm put them (on a surviving copy / mapped only
           elsewhere / unmapped) -> decompose.json + stdout
  Amendment 14 (the no-deletion control, ACC=a14; `nets` with --bam control/R0.bam, `label` with --all-copies):
  classify (C1') every flagged union against KB3781's mat / pat by Amendment 9's rule (`control_test.classify`: a haplotype-only locus /
           b allele / c unmatched / pri), self-checked on A9's own 28 candidates; overlaps with A9's candidates and with the A13 run's
           survivor-derived flags -> candidates_classified.tsv, classify.json
  cscore   (C2') arms R0 and C scored as `control_test.score` (each flagged union its own locus), self-checked on A9's R0 -> score.json, calls.tsv
  creport  wall time, stage totals, the per-family counters, pass B per batch, what the joined reads become, a per-family table -> report.json
  wplan / wreport  ruling R23's whole-BAM cost: batches of the 378-family copies table; per batch Elapsed, peak RSS, attribution set, flagged

    accept_o3_candidates.py nets --w /mnt/linuxdisk/tmp/rna_allele/a13 --linktest /mnt/linuxdisk/tmp/rna_allele/linktest
"""
import argparse
import collections
import csv
import glob
import json
import os
import re
import shutil
import struct
import subprocess
import sys
import tempfile

import pysam

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import control_test  # noqa: E402  (Amendment 9's haplotype names and truth paths; Amendment 14 reuses its classify rule)

MAX_READS = 1000
TIE = 0.98
M64 = (1 << 64) - 1
GOLDEN = 0x9E3779B97F4A7C15
COMP = str.maketrans("ACGTNacgtn", "TGCANtgcan")
MM2_KEEP = ["-c", "-x", "splice:hq", "-uf", "-N", "10", "-t", "4"]
# o3_candidates.rs (A13), pass B: MIN_UNMAPPED_LEN (both classes, Amendment 13c / R19), POORLY_PLACED_DE (an f32 compared with the f32 `de`),
# ATTRIB_MIN_READ_COV, ATTRIB_MAX_DE, MM2_ATTRIB (the stage appends `-t <threads>`; its `--threads 4`)
MIN_ATTRIB_LEN = 300
POORLY_PLACED_DE = struct.unpack("f", struct.pack("f", 0.02))[0]
ATTRIB_MIN_READ_COV, ATTRIB_MAX_DE = 0.5, 0.20
MM2_ATTRIB = ["-x", "map-ont", "-c", "-N", "5", "-p", "0.5", "-t", "4"]
PASS_B = re.compile(r"pass B: unmapped >= \d+ bp (\d+); poorly placed >= \d+ bp (\d+) \(.*?; (\d+) below the floor\); aligned (\d+): unmapped (\d+), "
                    r"poorly placed (\d+); attributed (\d+) \(.*?\): unmapped (\d+), poorly placed (\d+); joined this run's families (\d+)")
PASS_B_KEYS = ("unmapped", "poorly_placed", "poorly_placed_short", "aligned", "aligned_unmapped", "aligned_poorly_placed", "attributed",
               "attributed_unmapped", "attributed_poorly_placed", "joined")
# the per-family phase-1 line of the A13 stage log (`ClusterLog`'s Display in src/bin/o3_candidates.rs)
CLUSTER_LOG = re.compile(r"\] (\S+): net \d+ reads \(\d+ used\) \| (\d+) read clusters >= --min-cluster; (\d+) empty consensus dropped; refinement "
                         r"split off (\d+) reads \((\d+) new clusters, (\d+) clusters fell under --min-cluster, (\d+) kept sets re-templated\); "
                         r"significance merge absorbed (\d+) clusters in (\d+) rounds \((\d+) absorptions undone[^)]*\); templates without an eligible "
                         r"member [^:]*: (\d+) the longest member with an aligned partner, (\d+) the longest member \(no member aligned\); (\d+) clusters")
CLUSTER_LOG_KEYS = ("read_clusters", "empty", "split_off", "split_clusters", "fell_under", "retemplated", "absorbed", "rounds", "undone",
                    "longest_aligned", "longest_unaligned", "final")


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


# ---- nets (A13): the stage's nets replicated (pass A from the BAM, pass B re-run per batch) and checked against its own products -------

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


def stage_batches(w):
    """[(g, the batch's families, its pass-B counts)] from the batch logs (`logs/stage_g<g>.log`): the families of the `--families` argument
    in time -v's `Command being timed` line, the counts of the stage's pass-B line"""
    out = []
    for p in sorted(glob.glob(f"{w}/logs/stage_g*.log"), key=lambda p: int(re.search(r"stage_g(\d+)\.log$", p).group(1))):
        g = int(re.search(r"stage_g(\d+)\.log$", p).group(1))
        text = open(p).read()
        fams = re.search(r'Command being timed: ".*?--families ([^\s"]+)', text).group(1).split(",")
        m = PASS_B.search(text)
        assert m, f"{p}: no pass-B line (not an A13 stage log?)"
        out.append((g, fams, dict(zip(PASS_B_KEYS, map(int, m.groups())))))
    return out


def copy_rows(path):
    rows = tsv(path)
    assert all(r.get("member_status", "member") != "partner" for r in rows), f"{path}: partner rows (not handled here)"
    return rows


def nets(a):
    rows = copy_rows(a.copies)
    order = list(dict.fromkeys(r["family_id"] for r in rows))           # group_families: first appearance in --copies = the run order
    lab = {r["read"]: r for r in tsv(f"{a.linktest}/labels.tsv")}
    fams_tsv = {r["family"]: r for r in tsv(f"{a.w}/{a.prefix}.families.tsv")}
    batches = stage_batches(a.w)
    batch_of = {f: g for g, fs, _ in batches for f in fs}
    assert sorted(batch_of) == sorted(fams_tsv), "the batch logs and families.tsv name different families"
    # pass A, for every family: names with a primary / secondary record overlapping a copy's locus extent (no supplementary)
    bam = pysam.AlignmentFile(a.bam)
    names = collections.defaultdict(set)
    for r in rows:
        lo, hi = (int(r["locus_start"]), int(r["locus_end"])) if r.get("locus_start", "NA") not in ("", "NA") else (int(r["start"]), int(r["end"]))
        for rd in bam.fetch(r["chrom"], lo, max(hi, lo + 1)):
            if not (rd.is_unmapped or rd.is_supplementary):
                names[r["family_id"]].add(rd.query_name)
    # pass B's sweep, once (file order): every unmapped record and every primary record, with what the attribution set is chosen on
    recs, seqs = [], {}
    for rd in pysam.AlignmentFile(a.bam).fetch(until_eof=True):
        if rd.is_unmapped:
            recs.append((True, rd.query_name, oriented(rd), None, None))
        elif not (rd.is_secondary or rd.is_supplementary):
            s = oriented(rd)
            if s:
                seqs.setdefault(rd.query_name, s)
            recs.append((False, rd.query_name, s, rd.get_tag("de") if rd.has_tag("de") else None, rd.mapping_quality))
    copies_fa = open(a.copies_fa).read()
    checked = {(r["family_id"], int(r["copy_idx"])) for r in rows}
    family_of_copy = {}
    for h in (ln[1:].split()[0] for ln in copies_fa.splitlines() if ln.startswith(">")):
        f, idx = h.split("|")[:2]
        if (f, int(idx)) in checked:
            family_of_copy[h] = f
    net, via, cls_of, res = {}, {}, {}, {}
    tmp = tempfile.mkdtemp(prefix="nets.", dir=a.w)
    for g, fams_g, logged in batches:
        run = [f for f in order if f in set(fams_g)]
        netted = set().union(*(names[f] for f in run))
        rep = collections.Counter()
        att_class = {}
        with open(f"{tmp}/attrib.fa", "w") as o:
            for unmapped, n, s, de, mapq in recs:
                if unmapped:
                    if len(s) < MIN_ATTRIB_LEN:
                        continue
                    c = "unmapped"
                elif n in netted or not ((de if de is not None else 0.0) > POORLY_PLACED_DE or mapq == 0):
                    continue
                elif len(s) < MIN_ATTRIB_LEN:
                    rep["poorly_placed_short"] += 1; continue
                else:
                    c = "poorly_placed"
                rep[c] += 1
                att_class.setdefault(n, c)
                o.write(f">{n} {c}\n{s}\n")
        fot = dict(family_of_copy)
        with open(f"{tmp}/targets.fa", "w") as o:
            for f in run:
                for n in sorted(names[f]):
                    if n in seqs:
                        o.write(f">{f}|{n}\n{seqs[n]}\n"); fot[f"{f}|{n}"] = f
            o.write(copies_fa)
        with open(f"{tmp}/attrib.paf", "w") as out, open(f"{tmp}/attrib.log", "w") as err:
            subprocess.run(["minimap2", *MM2_ATTRIB, f"{tmp}/targets.fa", f"{tmp}/attrib.fa"], stdout=out, stderr=err, check=True)
        best = {}                                                          # attribute_by_hits: best hit by matches, the first on a tie
        for ln in open(f"{tmp}/attrib.paf"):
            x = ln.rstrip("\n").split("\t")
            if len(x) < 12 or x[5] not in fot:
                continue
            m = int(x[9])
            if x[0] not in best or m > best[x[0]][0]:
                de = next((float(t[5:]) for t in x[12:] if t.startswith("de:f:")), 1.0)
                best[x[0]] = (m, max(0, int(x[3]) - int(x[2])) / max(1, int(x[1])), de, fot[x[5]])
        attributed = {q: v[3] for q, v in best.items() if v[1] >= ATTRIB_MIN_READ_COV and v[2] <= ATTRIB_MAX_DE}
        for q in best:
            rep["aligned_" + att_class[q]] += 1
        for q in attributed:
            rep["attributed_" + att_class[q]] += 1
        rep["aligned"], rep["attributed"] = len(best), len(attributed)
        joined = {q: f for q, f in attributed.items() if f in set(run)}
        rep["joined"] = len(joined)
        assert all(rep[k] == logged[k] for k in PASS_B_KEYS), f"batch {g}: replicated pass B {dict(rep)} != the stage log {logged}"
        right = collections.Counter()
        for q, f in attributed.items():
            r = lab.get(q)
            right[(att_class[q], "joined" if q in joined else "elsewhere", "right" if r and r["family"] == f else "wrong")] += 1
        for f in run:
            net[f] = sorted({n for n in names[f] if n in seqs} | {q for q, jf in joined.items() if jf == f})
            via[f] = {n: ("B" if n in joined else "A") for n in net[f]}
            for q in (q for q, jf in joined.items() if jf == f):
                cls_of[(f, q)] = att_class[q]
        # nets.fa: each read of the flagged families' whole nets once (R9)
        nfa = set(fasta(f"{a.w}/{a.prefix}_g{g}.nets.fa"))
        flagged = [f for f in run if int(fams_tsv[f]["n_flagged"]) > 0]
        want = set().union(*(set(net[f]) for f in flagged)) if flagged else set()
        assert nfa == want, f"batch {g}: nets.fa ({len(nfa)} reads) != the replicated nets of the flagged families ({len(want)})"
        res[g] = dict(families=run, logged=logged, replicated=dict(rep), attributed_by={"|".join(k): v for k, v in sorted(right.items())})
    shutil.rmtree(tmp)
    # the cap and the stage's own products
    in_rt = collections.defaultdict(set)
    for r in tsv(f"{a.w}/{a.prefix}.reads.tsv"):
        in_rt[r["family"]].add(r["read"])
    used = {}
    for f in fams_tsv:
        used[f] = set(sample_net(net[f]))
        st = fams_tsv[f]
        assert (len(net[f]), len(used[f])) == (int(st["n_net"]), int(st["n_used"])), \
            f"{f}: replicated n_net / n_used {len(net[f])} / {len(used[f])} != the stage's {st['n_net']} / {st['n_used']}"
        assert in_rt[f] <= used[f], f"{f}: {len(in_rt[f] - used[f])} reads.tsv reads outside the replicated used net"
    role = lambda f, n: lab[n]["role"] if n in lab and lab[n]["family"] == f else ("other" if n in lab else "unlabelled")
    with open(f"{a.w}/nets.tsv", "w") as o, open(f"{a.w}/net_reads.tsv", "w") as nr:
        o.write("family\tbatch\tn_net\tn_used\tD_net\tD_used\tS_net\tS_used\tother_net\tvia_B\tD_via_B\tD_via_B_used\tD_via_B_in_reads_tsv\t"
                "other_via_B\n")
        nr.write("family\tread\tused\trole\tvia\tclass\tin_reads_tsv\n")
        for f in sorted(fams_tsv, key=name_key):
            c = collections.Counter((role(f, n), n in used[f]) for n in net[f])
            vb = [n for n in net[f] if via[f][n] == "B"]
            dvb = [n for n in vb if role(f, n) == "D"]
            o.write(f"{f}\t{batch_of[f]}\t{len(net[f])}\t{len(used[f])}\t{c[('D', True)] + c[('D', False)]}\t{c[('D', True)]}\t"
                    f"{c[('S', True)] + c[('S', False)]}\t{c[('S', True)]}\t{sum(v for (r, _), v in c.items() if r not in ('D', 'S'))}\t{len(vb)}\t"
                    f"{len(dvb)}\t{sum(1 for n in dvb if n in used[f])}\t{sum(1 for n in dvb if n in in_rt[f])}\t{len(vb) - len(dvb)}\n")
            for n in net[f]:
                nr.write(f"{f}\t{n}\t{int(n in used[f])}\t{role(f, n)}\t{via[f][n]}\t{cls_of.get((f, n), 'pass_A')}\t{int(n in in_rt[f])}\n")
    # reads in nets of more than one batch (ruling R18: a read netted in one batch may be attributed in another)
    batches_of = collections.defaultdict(set)
    viab = collections.defaultdict(set)
    for f in fams_tsv:
        for n in net[f]:
            batches_of[n].add(batch_of[f])
            if via[f][n] == "B":
                viab[n].add(batch_of[f])
    multi = [n for n, bs in batches_of.items() if len(bs) > 1]
    multi_b = [n for n in multi if viab[n]]
    # report
    tot_log = collections.Counter()
    for g, _, logged in batches:
        tot_log.update(logged)
    print("pass B per batch (the stage log's counts; the replication re-runs the attribution and equals them in every batch):")
    print("batch\tfamilies\tunmapped>=300\tpoorly placed>=300 (short)\taligned u / p\tattributed u / p\tjoined")
    for g, fs, lg in batches:
        print(f"{g}\t{len(fs)}\t{lg['unmapped']}\t{lg['poorly_placed']} ({lg['poorly_placed_short']})\t{lg['aligned_unmapped']} / "
              f"{lg['aligned_poorly_placed']}\t{lg['attributed_unmapped']} / {lg['attributed_poorly_placed']}\t{lg['joined']}")
    print(f"sum\t{len(batch_of)}\t{tot_log['unmapped']}\t{tot_log['poorly_placed']} ({tot_log['poorly_placed_short']})\t"
          f"{tot_log['aligned_unmapped']} / {tot_log['aligned_poorly_placed']}\t{tot_log['attributed_unmapped']} / "
          f"{tot_log['attributed_poorly_placed']}\t{tot_log['joined']}")
    rt = collections.Counter()
    for k in res.values():
        for key, v in k["attributed_by"].items():
            rt[key] += v
    print("attributed reads by class, joined (a family of the batch) or not, right family by labels.tsv (summed over batches; a read may be "
          "attributed in several batches):")
    for c in ("unmapped", "poorly_placed"):
        for j in ("joined", "elsewhere"):
            r_, w_ = rt[f"{c}|{j}|right"], rt[f"{c}|{j}|wrong"]
            print(f"  {c:13s} {j:9s}: {r_ + w_:5d}, right family {r_:5d} ({r_ / max(1, r_ + w_):.1%})")
    jn = [(f, n) for f in fams_tsv for n in net[f] if via[f][n] == "B"]
    jr = sum(1 for f, n in jn if n in lab and lab[n]["family"] == f)
    jv = [(f, n) for f, n in jn if n in in_rt[f]]
    jvr = sum(1 for f, n in jv if n in lab and lab[n]["family"] == f)
    jd = sum(1 for f, n in jn if role(f, n) == "D")
    print(f"joined reads (in a net of their batch): {len(jn)} (unmapped {sum(1 for f, n in jn if cls_of[(f, n)] == 'unmapped')}, poorly placed "
          f"{sum(1 for f, n in jn if cls_of[(f, n)] == 'poorly_placed')}); right family {jr} ({jr / max(1, len(jn)):.1%}; D reads of the family "
          f"{jd}); under the cap {sum(1 for f, n in jn if n in used[f])}; in reads.tsv (a reported cluster) {len(jv)}, right family {jvr} "
          f"({jvr / max(1, len(jv)):.1%}) — reads.tsv lists only the new-copy and linked clusters' reads, so it sees {len(jv)} of {len(jn)}")
    print(f"reads in nets of more than one batch: {len(multi)} (of {len(batches_of)} netted reads); attributed (pass B) in at least one of "
          f"them: {len(multi_b)} (ruling R18)")
    print(f"checks: n_net and n_used of families.tsv equal in {len(fams_tsv)}/{len(fams_tsv)} families; every reads.tsv read under the cap of "
          f"its family; nets.fa = the flagged families' whole nets in every batch")
    json.dump(dict(batches=res, multi_batch=len(multi), multi_batch_attributed=len(multi_b), joined=len(jn), joined_right=jr, joined_D=jd,
                   joined_in_reads_tsv=len(jv), joined_in_reads_tsv_right=jvr, totals_logged=dict(tot_log)),
              open(f"{a.w}/attrib.json", "w"), indent=0)


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
                first = ("S:" + P[fam]["mask"][3]) if a.all_copies else "D"
                for lab_, (c, s0, e, g) in [(first, P[fam]["mask"])] + [("S:" + kk[3], kk) for kk in P[fam]["keep"]]:
                    if b[1] == c and b[2] < e and s0 < b[3]:
                        src = lab_; break
                else:
                    src = "elsewhere"
            src_count["S" if src.startswith("S:") else src] += 1
            t.write(f"{r['iso']}\t{fam}\t{r['candidate']}\t{r['length']}\tNA\t{r['d']}\t0\t{src}\n")
    for fam in sorted({r["family"] for r in iso}):
        open(f"{a.w}/merge/paf/{fam}.paf", "w").close()
    # merge_test.py score's inputs (arm M); Amendment 14 scores its own arm C (`cscore`) and never reads the masked run's R.bam
    for f in ("R.bam", "R.bam.bai", "labels.tsv", "panel.json") if not a.all_copies else ():
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
    # 2b. the stage log's per-family phase-1 counters (A13 logs: kept sets re-templated, absorptions undone, fallback templates by kind)
    counters, logs = {}, sorted(glob.glob(f"{a.w}/logs/stage_g*.log"))
    for p in logs:
        for m in CLUSTER_LOG.finditer(open(p).read()):
            counters[m.group(1)] = dict(zip(CLUSTER_LOG_KEYS, map(int, m.groups()[1:])))
    if counters:
        ct = collections.Counter()
        for v in counters.values():
            ct.update(v)
        fam_with = lambda k: sorted((f for f, v in counters.items() if v[k] > 0), key=name_key)
        out["counters"] = dict(total=dict(ct), per_family=counters)
        print(f"stage log counters over {len(counters)} families (phase 1): read clusters >= --min-cluster {ct['read_clusters']}; empty consensus "
              f"dropped {ct['empty']}; refinement split off {ct['split_off']} reads into {ct['split_clusters']} new clusters, {ct['fell_under']} "
              f"clusters fell under --min-cluster, kept sets re-templated {ct['retemplated']} (families {fam_with('retemplated')}); significance "
              f"merge absorbed {ct['absorbed']} clusters in {ct['rounds']} rounds, absorptions undone {ct['undone']}; fallback templates (no eligible "
              f"member): longest with an aligned partner {ct['longest_aligned']}, longest overall {ct['longest_unaligned']} (every other template "
              f"choice is the structural medoid, or the kept template); final clusters {ct['final']}")
    elif logs:
        print("stage log counters: none (not an A13 stage log)")
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


def comparator(a):
    """Amendment 13b's comparator and the A13-1 verdict (see the module doc)"""
    P = {p["fam"]: p for p in json.load(open(f"{a.linktest}/panel.json"))}
    lab = {r["read"]: r for r in tsv(f"{a.linktest}/labels.tsv")}
    nD = sum(1 for r in lab.values() if r["role"] == "D"); nS = len(lab) - nD
    iso = arm_m(a.linktest, lab, P)
    iso.pop("__false_moves__")
    tot = collections.Counter(iso.values())
    ref = json.load(open(f"{a.linktest}/merge/score.json"))["M"]
    assert {"|".join(k): v for k, v in tot.items()} == ref, f"IsoCon's arm M replication differs from merge_test.py score: {dict(tot)} vs {ref}"
    # (1) a record, primary or secondary, on a surviving copy of the read's family in R.bam
    on_surv = set()
    for rd in pysam.AlignmentFile(f"{a.linktest}/R.bam").fetch(until_eof=True):
        r = lab.get(rd.query_name)
        if rd.is_unmapped or rd.is_supplementary or not r or r["role"] != "D":
            continue
        if any(rd.reference_name == k[0] and rd.reference_start < k[2] and k[1] < rd.reference_end for k in P[r["family"]]["keep"]):
            on_surv.add(rd.query_name)
    # (2) attributed (pass B) into the read's OWN family's net by this run, without a record on a survivor
    netr = tsv(f"{a.w}/net_reads.tsv")
    own = lambda r: r["read"] in lab and lab[r["read"]]["role"] == "D" and lab[r["read"]]["family"] == r["family"]
    att = {r["read"] for r in netr if r["via"] == "B" and own(r)} - on_surv
    att_cls = {r["read"]: r["class"] for r in netr if r["via"] == "B" and own(r)}
    # (2') the same read from reads.tsv alone: the D reads of the family in a reported cluster (new-copy or linked), without a record on a survivor
    rt = {r["read"] for r in tsv(f"{a.w}/{a.prefix}.reads.tsv") if own(r)} - on_surv
    right = lambda S: sum(1 for n in S if iso[n] == ("D", "right"))
    C1, C2, C2rt = right(on_surv), right(att), right(rt)
    C, Crt = C1 + C2, C1 + C2rt
    by_cls = collections.Counter((att_cls[n], iso[n] == ("D", "right")) for n in att)
    out = dict(attainable_on_survivor=len(on_surv), attainable_attributed=len(att), attainable_attributed_reads_tsv=len(rt),
               C1=C1, C2=C2, C=C, bar=0.8 * C, C2_reads_tsv=C2rt, C_reads_tsv=Crt, bar_reads_tsv=0.8 * Crt,
               isocon_right_total=tot[("D", "right")], attributed_by_class={f"{c}|{'right' if k else 'not_right'}": v for (c, k), v in by_cls.items()})
    print(f"IsoCon (Amendment 8, arm M read by read; equals merge_test.py score: D right {tot[('D', 'right')]}, S false moves "
          f"{tot[('S', 'false_move')]})")
    print(f"attainable D reads (truth-free, Amendment 13b): (1) a record on a surviving copy of their family in R.bam {len(on_surv)}; (2) else "
          f"attributed into their own family's net by this run {len(att)} (unmapped {sum(v for (c, _), v in by_cls.items() if c == 'unmapped')}, "
          f"poorly placed {sum(v for (c, _), v in by_cls.items() if c == 'poorly_placed')}); total {len(on_surv) + len(att)} of {nD}")
    print(f"C = IsoCon's right D reads over them = C1 {C1} + C2 {C2} = {C}; bar 0.80 x C = {0.8 * C:.1f}")
    print(f"  (part 2 read from reads.tsv alone — the attributed D reads that reached a reported cluster: {len(rt)} reads, C2 {C2rt}, C {Crt}, "
          f"bar {0.8 * Crt:.1f})")
    sj = f"{a.w}/merge/score.json"
    if os.path.exists(sj):
        st = json.load(open(sj))["M"]
        dr, fm = st.get("D|right", 0), st.get("S|false_move", 0)
        ok = dr >= 0.8 * C and fm <= 0.05 * nS
        okrt = dr >= 0.8 * Crt and fm <= 0.05 * nS
        stage = arm_m(a.w, lab, P)
        stage.pop("__false_moves__")
        dr_att = sum(1 for n in on_surv | att if stage[n] == ("D", "right"))
        out.update(stage_D_right=dr, stage_false_moves=fm, n_S=nS, A13_1=ok, A13_1_reads_tsv=okrt, stage_D_right_on_attainable=dr_att,
                   A12_1_bar=10230, A12_1=dr >= 10230 and fm <= 0.05 * nS)
        print(f"stage (arm M, merge_test.py score): D right {dr} ({dr / max(1, C):.1%} of C; on the attainable reads {dr_att}), false moves {fm}/{nS} "
              f"= {fm / nS:.2%}")
        print(f"A13-1: D right {dr} >= 0.80 x C = {0.8 * C:.1f} {'yes' if dr >= 0.8 * C else 'NO'}; false moves {fm / nS:.2%} <= 5% "
              f"{'yes' if fm <= 0.05 * nS else 'NO'} -> {'PASSES' if ok else 'FAILS'} (with part 2 from reads.tsv alone: "
              f"{'PASSES' if okrt else 'FAILS'})")
        print(f"beside (not decided on): A12-1's bar 10,230 (80% of IsoCon's 12,787): {'met' if dr >= 10230 else 'not met'}")
    else:
        print(f"no {sj}: arm M not scored for this run (its flagged contig set equals the registered run's, or the step has not run)")
    json.dump(out, open(f"{a.w}/comparator.json", "w"), indent=0)


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


# ---- Amendment 14: the no-deletion control of the stage (classify = Amendment 9's rule, cscore = arm C, creport) and R23's plan -------

def all_copies(linktest):
    """`control_test.panel`: family -> every copy [chrom, start, end, gene], the copy Amendment 7 masked first"""
    return {p["fam"]: [p["mask"]] + p["keep"] for p in json.load(open(f"{linktest}/panel.json"))}


def hap_best(paths):
    """`control_test.classify`'s best haplotype hit per query over the mat, then the pat PAF: (identity x coverage, hap, accession, start,
    end), the first on a tie; a target named chr<N>_<hap>_hsa* (the haplotype splice indexes) -> its accession by <hap>.len.tsv, any other
    target skipped"""
    alias = control_test.hap_alias()
    hb = {}
    for h, path in paths:
        for ln in open(path):
            f = ln.rstrip("\n").split("\t")
            m = re.fullmatch(r"chr(\w+?)_(mat|pat)_hsa[^_]*", f[5])
            if not m:
                continue
            acc = alias[(m.group(2), m.group(1))]
            s = int(f[9]) / max(1, int(f[10])) * (int(f[3]) - int(f[2])) / max(1, int(f[1]))
            if f[0] not in hb or s > hb[f[0]][0]:
                hb[f[0]] = (s, h, acc, int(f[7]), int(f[8]))
    return hb


def hap_class(c, fam, hb, P, bside, liftd):
    """`control_test.classify`'s per-contig class (its `cls_contig`): c_unmatched = no haplotype hit at identity x coverage >= 0.999; pri =
    the hit lies on the haplotype `_pri` took that chromosome from (in `_pri` after all: an artefact); b_allele = on the other (B)
    haplotype, inside the lifted B interval of a copy of the family (lift_frac >= 0.5); a_haplotype_only = on B, outside every such
    interval (a genuine reference-absent locus: the only TRUE flag)"""
    b = hb.get(c)
    if not b or b[0] < 0.999:
        return "c_unmatched"
    _, h, acc, s, e = b
    if (h, acc) not in bside:
        return "pri"
    for cc, s0, e0, g in P[fam]:
        L = liftd.get(g)
        if L and float(L["lift_frac"]) >= 0.5 and L["B_chrom"] == acc and int(L["B_start"]) < e and s < int(L["B_end"]):
            return "b_allele"
    return "a_haplotype_only"


def same_transcript(path):
    """`control_test.classify`'s overlap rule on an asm20 PAF: the (query, target) pairs of the same family (iso_<family>_<k> names) at
    identity >= 0.999 over >= 90% of the shorter sequence"""
    pairs = set()
    for ln in open(path):
        f = ln.split("\t")
        ident = int(f[9]) / max(1, int(f[10]))
        short = min(int(f[1]), int(f[6]))
        cov = (int(f[3]) - int(f[2])) / short if int(f[1]) <= int(f[6]) else (int(f[8]) - int(f[7])) / short
        if ident >= 0.999 and cov >= 0.9 and f[5].split("_")[1] == f[0].split("_")[1]:
            pairs.add((f[0], f[5]))
    return pairs


A13_DETECTED = 25          # A13's family-level detection rate 25/53 (docs/archive/2026-10/O3_CANDIDATES_ACCEPTANCE_A13_2026-10-03.md)
C1_MAX_FAMILIES = 8        # Amendment 14: <= 1/3 of 25/53 = 0.157 -> <= 8 of 53 families
A9_FALSE_FAMILIES = 16     # Amendment 9: the IsoCon chain, 16/53 families with a false candidate at any support


def classify(a):
    """Amendment 14's C1': every flagged union (iso.contigs.fa, iso_<family>_<k>) classified against KB3781's mat / pat assemblies by
    Amendment 9's rule (iso.mat.paf / iso.pat.paf: minimap2 -c -x splice:hq -uf -N 20; chrmap.tsv; control/copies_lift.tsv). Self-check
    first: Amendment 9's 28 candidates re-classified from its own PAFs must equal control/candidates.tsv. Beside: the a + b + c rate, A9's
    16/53, the overlap with A9's 28 candidates (control/contigs_L.fa, components of control/merge/components.tsv) and with the A13 run's
    survivor-derived flags (A13 contigs.tsv source S:*), both by `same_transcript` -> candidates_classified.tsv, classify.json + stdout"""
    P = all_copies(a.linktest)
    n_fam = len(P)
    bside = {(r["B_hap"], r["B_name"]) for r in tsv(f"{control_test.TRUTH}/chrmap.tsv") if r["B_name"]}
    liftd = {r["gene_id"]: r for r in tsv(f"{a.control}/copies_lift.tsv")}
    hb9 = hap_best([("mat", f"{a.control}/contigs_L.mat.paf"), ("pat", f"{a.control}/contigs_L.pat.paf")])
    a9 = tsv(f"{a.control}/candidates.tsv")
    mism = [r["candidate"] for r in a9 if hap_class(r["best_contig"], r["family"], hb9, P, bside, liftd) != r["class"]]
    assert not mism, f"Amendment 9's classes not reproduced from its own PAFs: {mism}"
    print(f"self-check: Amendment 9's {len(a9)} candidates re-classified from control/contigs_L.{{mat,pat}}.paf equal control/candidates.tsv "
          f"({dict(collections.Counter(r['class'] for r in a9))})")
    iso = tsv(f"{a.w}/iso_names.tsv")
    src = {r["contig"]: r["source"] for r in tsv(f"{a.w}/contigs.tsv")}
    hb = hap_best([("mat", f"{a.w}/iso.mat.paf"), ("pat", f"{a.w}/iso.pat.paf")])
    comp9 = {r["contig"]: r["component"] for r in tsv(f"{a.control}/merge/components.tsv")}
    ov9 = same_transcript(f"{a.w}/overlap_a9.paf")
    a13_s = {r["contig"]: r for r in tsv(f"{a.a13}/contigs.tsv") if r["source"].startswith("S:")}
    ov13 = {(q, t) for q, t in same_transcript(f"{a.w}/overlap_a13.paf") if t in a13_s}
    rows = []
    for r in iso:
        c, fam = r["iso"], r["family"]
        b = hb.get(c)
        rows.append(dict(iso=c, candidate=r["candidate"], family=fam, n_clusters=int(r["n_clusters"]), n_reads=int(r["n_reads"]), d=r["d"],
                         length=int(r["length"]), source=src[c], cls=hap_class(c, fam, hb, P, bside, liftd),
                         hap_hit=f"{b[1]}:{b[2]}:{b[3]}-{b[4]} {b[0]:.4f}" if b else "-",
                         a9=sorted({comp9[t] for q, t in ov9 if q == c}), a13=sorted({t for q, t in ov13 if q == c})))
    with open(f"{a.w}/candidates_classified.tsv", "w") as o:
        o.write("iso\tcandidate\tfamily\tn_clusters\tn_reads\tlength\td\tsource\tclass\thap_hit\tmatches_A9_candidates\tmatches_A13_survivor_flags\n")
        for x in rows:
            o.write(f"{x['iso']}\t{x['candidate']}\t{x['family']}\t{x['n_clusters']}\t{x['n_reads']}\t{x['length']}\t{x['d']}\t{x['source']}\t{x['cls']}\t"
                    f"{x['hap_hit']}\t{','.join(x['a9']) or '-'}\t{','.join(x['a13']) or '-'}\n")
    classes = collections.Counter(x["cls"] for x in rows)
    fam_any = {x["family"] for x in rows}
    fam_false = {x["family"] for x in rows if x["cls"] != "a_haplotype_only"}
    fam_true = {x["family"] for x in rows if x["cls"] == "a_haplotype_only"}
    ok = len(fam_false) <= C1_MAX_FAMILIES
    det = A13_DETECTED / n_fam
    lr = det / (len(fam_false) / n_fam) if fam_false else float("inf")
    per_fam = collections.defaultdict(collections.Counter)
    for x in rows:
        per_fam[x["family"]][x["cls"]] += 1
    by_src = collections.Counter((("S" if x["source"].startswith("S:") else x["source"]), x["cls"]) for x in rows)
    m9 = [x for x in rows if x["a9"]]
    a9_hit = sorted({cpn for x in rows for cpn in x["a9"]})
    a9_fams = {r["family"] for r in a9}
    m13 = [x for x in rows if x["a13"]]
    a13_hit = sorted({t for x in rows for t in x["a13"]})
    out = dict(flagged=len(rows), classes=dict(classes), families_any=sorted(fam_any, key=name_key), families_false=sorted(fam_false, key=name_key),
               families_true=sorted(fam_true, key=name_key), n_families=n_fam, C1_prime=ok, C1_bar_families=C1_MAX_FAMILIES, LR=lr,
               by_source_class={"|".join(k): v for k, v in sorted(by_src.items())},
               per_family={f: dict(v) for f, v in sorted(per_fam.items(), key=lambda kv: name_key(kv[0]))},
               match_A9=dict(flags=len(m9), A9_candidates_reproduced=len(a9_hit), A9_candidates=len(a9), families_both=sorted(fam_any & a9_fams, key=name_key),
                             A9_families=len(a9_fams)),
               match_A13_survivor=dict(flags=len(m13), A13_survivor_flags_reproduced=len(a13_hit), A13_survivor_flags=len(a13_s)))
    print(f"flagged candidates {len(rows)} in {len(fam_any)} families; classes {dict(classes)}")
    print(f"by stage label x class: {out['by_source_class']}")
    print(f"C1': families with >= 1 FALSE flag (b + c + pri) {len(fam_false)}/{n_fam} = {len(fam_false) / n_fam:.1%} (bar <= {C1_MAX_FAMILIES} = "
          f"{C1_MAX_FAMILIES / n_fam:.1%}, i.e. <= 1/3 of A13's detection {A13_DETECTED}/{n_fam} = {det:.3f}) -> {'HOLDS' if ok else 'FAILS'}; "
          f"family-level LR of a flag {det:.3f} / {len(fam_false) / n_fam:.3f} = {lr:.2f} (bar 3)")
    print(f"beside: families with any flag (a + b + c) {len(fam_any)}/{n_fam} = {len(fam_any) / n_fam:.1%}; with a TRUE flag (a) {len(fam_true)} "
          f"{sorted(fam_true, key=name_key)}; A9's IsoCon chain (any support): {A9_FALSE_FAMILIES}/{n_fam} = {A9_FALSE_FAMILIES / n_fam:.1%} with a false candidate")
    print(f"overlap with A9's {len(a9)} candidates (same family, identity >= 0.999 over >= 90% of the shorter): {len(m9)} of the {len(rows)} flags match "
          f"one; A9 candidates reproduced {len(a9_hit)}/{len(a9)}; families flagged by both {len(fam_any & a9_fams)} (A9 {len(a9_fams)}, A14 {len(fam_any)})")
    print(f"overlap with the A13 run's {len(a13_s)} survivor-derived flags: {len(m13)} of the {len(rows)} flags match one; A13 survivor-derived flags "
          f"reproduced {len(a13_hit)}/{len(a13_s)}")
    def describe(pop):
        hh = sorted(float(x["hap_hit"].split()[-1]) for x in pop if x["hap_hit"] != "-")
        ds = sorted(float(x["d"]) for x in pop)
        rd = sorted(x["n_reads"] for x in pop)
        med = lambda v: (v[(len(v) - 1) // 2] + v[len(v) // 2]) / 2
        return (f"n {len(pop)}; best haplotype hit (identity x coverage, {len(hh)} with a hit) median {med(hh):.4f} ({hh[0]:.4f}-{hh[-1]:.4f}), "
                f">= 0.99 {sum(v >= 0.99 for v in hh)}, >= 0.95 {sum(v >= 0.95 for v in hh)}, < 0.90 {sum(v < 0.90 for v in hh)}; whole-length d "
                f"median {med(ds):.5f} ({ds[0]:.5f}-{ds[-1]:.5f}), <= 0.02 {sum(v <= 0.02 for v in ds)}, > 0.10 {sum(v > 0.10 for v in ds)}; "
                f"reads median {med(rd):g} ({rd[0]}-{rd[-1]}), >= 100 {sum(v >= 100 for v in rd)}, <= 10 {sum(v <= 10 for v in rd)}") if pop else "n 0"
    print(f"the unmatched (c) flags: {describe([x for x in rows if x['cls'] == 'c_unmatched'])}")
    print(f"all false flags (b + c + pri): {describe([x for x in rows if x['cls'] != 'a_haplotype_only'])}")
    print("per family (flagged candidates by class; a = haplotype-only, b = allele, c = unmatched):")
    for f in sorted(per_fam, key=name_key):
        v = per_fam[f]
        print(f"  {f}\t" + " ".join(f"{k.split('_')[0]}:{n}" for k, n in sorted(v.items())) + ("\tFALSE" if f in fam_false else "\ttrue only"))
    print("candidate\tfamily\treads\tclusters\tlength\td\tlabel\tclass\tbest haplotype hit\tA9\tA13-S")
    for x in rows:
        print(f"{x['iso']}\t{x['family']}\t{x['n_reads']}\t{x['n_clusters']}\t{x['length']}\t{x['d']}\t{x['source']}\t{x['cls']}\t{x['hap_hit']}\t"
              f"{','.join(x['a9']) or '-'}\t{','.join(x['a13']) or '-'}")
    json.dump(out, open(f"{a.w}/classify.json", "w"), indent=0)


def cscore(a):
    """Amendment 14's C2': arms R0 (control/R0.bam) and C (C.bam = the scored reads realigned to `_pri` + the flagged unions) scored as
    `control_test.score` (each flagged candidate its own locus: the stage's candidate is already its component): stay / other_copy /
    own_candidate (a candidate labelled with the read's own copy) / false_move (any other candidate) / elsewhere / unplaced (no record, or
    a 0.98 AS tie across loci). Self-check: arm R0 equals Amendment 9's control/score.json -> score.json, calls.tsv + stdout"""
    P = all_copies(a.linktest)
    lab = {r["read"]: r for r in tsv(f"{a.linktest}/labels.tsv")}
    src = {r["contig"]: r["source"] for r in tsv(f"{a.w}/contigs.tsv")}
    want = set(lab)

    def locus(chrom, s, e, fam):
        if chrom.startswith("iso_"):
            return ("ctg", chrom)
        for c, s0, e0, g in P[fam]:
            if chrom == c and s < e0 and s0 < e:
                return ("copy", g)
        return ("other", f"{chrom}:{s // 100000}")

    def calls(bam, ties=None):
        recs = collections.defaultdict(list)
        for rd in pysam.AlignmentFile(bam).fetch(until_eof=True):
            if rd.query_name not in want or rd.is_supplementary:
                continue
            recs[rd.query_name].append(None if rd.is_unmapped else
                                       (rd.is_secondary, rd.reference_name, rd.reference_start, rd.reference_end,
                                        rd.get_tag("AS") if rd.has_tag("AS") else 0))
        out = {}
        for n, rs in recs.items():
            rs = [r for r in rs if r]
            if not rs:
                out[n] = ("unplaced", None); continue
            fam = lab[n]["family"]
            prim = next((r for r in rs if not r[0]), rs[0])
            srt = sorted(rs, key=lambda r: -r[4])
            if ties is not None:                      # the loci at >= 0.98 x the best AS (what a tie is between; beside, no rule)
                ties[n] = {locus(r[1], r[2], r[3], fam) for r in srt if r[4] >= TIE * srt[0][4]}
            if len(srt) > 1 and srt[1][4] > 0 and srt[1][4] >= TIE * srt[0][4]:
                if len({locus(r[1], r[2], r[3], fam) for r in srt if r[4] >= TIE * srt[0][4]}) > 1:
                    out[n] = ("unplaced", None); continue
            out[n] = ("placed", locus(prim[1], prim[2], prim[3], fam))
        return out

    def cls(n, call):
        r = lab[n]
        st, L = call
        if st == "unplaced":
            return "unplaced"
        kind, key = L
        if kind == "copy":
            return "stay" if key == r["copy"] else "other_copy"
        if kind == "ctg":
            return "own_candidate" if src[key] == "S:" + r["copy"] else "false_move"
        return "elsewhere"
    res, per, ties_c = {}, {}, {}
    for arm, bam in (("R0", f"{a.control}/R0.bam"), ("C", f"{a.w}/C.bam")):
        cc = calls(bam, ties_c if arm == "C" else None)
        per[arm] = {n: (cls(n, cc.get(n, ("unplaced", None))), cc.get(n, ("unplaced", None))[1]) for n in lab}
        res[arm] = collections.Counter(v[0] for v in per[arm].values())
        print(f"[{arm}]", dict(sorted(res[arm].items())))
    ref = json.load(open(f"{a.control}/score.json"))["R0"]
    assert dict(res["R0"]) == ref, f"arm R0 differs from Amendment 9's control/score.json: {dict(res['R0'])} vs {ref}"
    print(f"self-check: arm R0 equals Amendment 9's control/score.json {ref}")
    n = len(lab)
    fm = res["C"]["false_move"]
    ok = fm <= 0.05 * n
    print(f"C2': false moves {fm}/{n} = {fm / n:.2%} (bar <= 5%) -> {'HOLDS' if ok else 'FAILS'}; own-candidate placements {res['C']['own_candidate']} "
          f"({res['C']['own_candidate'] / n:.1%}); unplaced {res['R0']['unplaced']} -> {res['C']['unplaced']}; stay {res['R0']['stay']} -> {res['C']['stay']}")
    to = collections.Counter()
    own = collections.Counter()
    unpl = collections.Counter()
    for r, (k, L) in per["C"].items():
        fam = lab[r]["family"]
        if k == "false_move":
            tgt = L[1]
            to[(fam, tgt, src[tgt], "own family" if tgt.split("_")[1] == fam else "other family")] += 1
        elif k == "own_candidate":
            own[(fam, L[1])] += 1
        if k == "unplaced" and per["R0"][r][0] != "unplaced":
            unpl[fam] += 1
    print("false moves by (family, candidate, its label, own family or not): " +
          ("; ".join(f"{f} -> {t} ({s_}, {o}): {v}" for (f, t, s_, o), v in sorted(to.items(), key=lambda kv: -kv[1])) or "none"))
    print("own-candidate placements by (family, candidate): " + "; ".join(f"{f} {t}: {v}" for (f, t), v in sorted(own.items(), key=lambda kv: -kv[1])))
    print("reads placed in R0 and unplaced in C, by family: " + "; ".join(f"{f} {v}" for f, v in sorted(unpl.items(), key=lambda kv: -kv[1])))
    # what those reads tie between in C (the loci at >= 0.98 x their best AS): their own copy and a flag derived from it (false or
    # true by `classify`'s candidates_classified.tsv, when it exists), or not
    kpath = f"{a.w}/candidates_classified.tsv"
    klass = {r["iso"]: r["class"] for r in tsv(kpath)} if os.path.exists(kpath) else {}
    tie_kind = collections.Counter()
    for r, (k, L) in per["C"].items():
        if k != "unplaced" or per["R0"][r][0] == "unplaced":
            continue
        own_copy = lab[r]["copy"]
        t = ties_c.get(r)
        flags = [key for kind, key in (t or ()) if kind == "ctg"]
        if not t:
            tie_kind["no mapped record in C"] += 1
        elif ("copy", own_copy) not in t:
            tie_kind["do not tie with their own copy"] += 1
        elif any(src[f] != "S:" + own_copy for f in flags):
            tie_kind["tie their own copy with a flag derived from another copy (or labelled elsewhere)"] += 1
        elif flags:
            k = "" if not klass else "the TRUE " if all(klass.get(f) == "a_haplotype_only" for f in flags) else "a false "
            tie_kind[f"tie their own copy with {k or 'a '}flag derived from it"] += 1
        else:
            tie_kind["tie their own copy with other loci only"] += 1
    print(f"those {sum(unpl.values())} reads by what they tie between in C: " + "; ".join(f"{v} {k}" for k, v in tie_kind.most_common()))
    with open(f"{a.w}/calls.tsv", "w") as o:
        o.write("read\tfamily\tcopy\trole\tR0\tC\tC_target\n")
        for r in sorted(lab):
            k, L = per["C"][r]
            o.write(f"{r}\t{lab[r]['family']}\t{lab[r]['copy']}\t{lab[r]['role']}\t{per['R0'][r][0]}\t{k}\t{L[1] if L else '-'}\n")
    json.dump(dict(arms={k: dict(v) for k, v in res.items()}, n_reads=n, false_moves=fm, C2_prime=ok,
                   false_moves_to={"|".join(k): v for k, v in to.items()}, own_candidate_by={"|".join(k): v for k, v in own.items()},
                   newly_unplaced_by_family=dict(unpl), newly_unplaced_ties=dict(tie_kind)), open(f"{a.w}/score.json", "w"), indent=0)


ELAPSED = re.compile(r"Elapsed \(wall clock\) time \(h:mm:ss or m:ss\): (\S+)")
MAXRSS = re.compile(r"Maximum resident set size \(kbytes\): (\d+)")
DONE = re.compile(r"done: (\d+) families, (\d+) clusters, (\d+) candidates, (\d+) flagged in (\d+) families; (\d+) minimap2 calls; ([\d.]+) s "
                  r"\(nets ([\d.]+) s, clusters ([\d.]+) s, genome ([\d.]+) s\)")
PASS_A = re.compile(r"BAM: pass A (\d+) reads by a primary record and (\d+) secondary records on the copies; pass B (\d+) of (\d+) secondary-only")


def hms(t):
    s = 0.0
    for x in t.split(":"):
        s = s * 60 + float(x)
    return s


def batch_logs(w, stem):
    """per batch log (logs/<stem>_g<g>.log): /usr/bin/time's Elapsed (s) and peak RSS (GB = 10^6 kB, as the A13 doc), the stage's done line,
    its pass A / pass B lines, and whether the run was stopped (a signal: the wstage step's time limit)"""
    out = {}
    for p in sorted(glob.glob(f"{w}/logs/{stem}_g*.log"), key=lambda p: int(re.search(r"_g(\d+)\.log$", p).group(1))):
        g = int(re.search(r"_g(\d+)\.log$", p).group(1))
        text = open(p).read()
        e, m, d, pb, pa = ELAPSED.search(text), MAXRSS.search(text), DONE.search(text), PASS_B.search(text), PASS_A.search(text)
        out[g] = dict(elapsed=hms(e.group(1)) if e else None, rss_gb=int(m.group(1)) / 1e6 if m else None,
                      done=dict(zip(("families", "clusters", "candidates", "flagged", "flagged_families", "mm2_calls", "total", "nets", "clusters_s",
                                     "genome"), map(float, d.groups()))) if d else None,
                      pass_b=dict(zip(PASS_B_KEYS, map(int, pb.groups()))) if pb else None,
                      pass_a=dict(zip(("primary", "secondary", "found", "need"), map(int, pa.groups()))) if pa else None,
                      stopped="Command terminated by signal" in text or "Command exited with non-zero status 124" in text)
    return out


def creport(a):
    """Amendment 14's reported numbers beside C1' / C2' (prereg: candidates per family by class, the attribution counts without a deletion
    and what the joined reads become, the per-family counters): wall time per batch, stage totals, the stage log's phase-1 counters, pass B
    per batch, the joined reads (nets replication: right family, under the cap, in a reported cluster / a flagged candidate, R0 and arm C
    class), and a per-family table -> report.json + stdout"""
    lab = {r["read"]: r for r in tsv(f"{a.linktest}/labels.tsv")}
    fams = {r["family"]: r for r in tsv(f"{a.w}/{a.prefix}.families.tsv")}
    cands = {r["candidate"]: r for r in tsv(f"{a.w}/{a.prefix}.candidates.tsv")}
    clus = {r["cluster"]: r for r in tsv(f"{a.w}/{a.prefix}.clusters.tsv")}
    cl_of = collections.defaultdict(dict)
    for r in tsv(f"{a.w}/{a.prefix}.reads.tsv"):
        cl_of[r["family"]][r["read"]] = r["cluster"]
    klass = {r["candidate"]: r for r in tsv(f"{a.w}/candidates_classified.tsv")} if os.path.exists(f"{a.w}/candidates_classified.tsv") else {}
    calls = {r["read"]: r for r in tsv(f"{a.w}/calls.tsv")} if os.path.exists(f"{a.w}/calls.tsv") else {}
    # ACC=a15h classifies and scores with control_test.py (A9's layouts: candidates.tsv / score.json, no candidates_classified.tsv /
    # calls.tsv): creport then reports the stage / nets / pass-B numbers without the per-flag class and per-read arm-C columns
    netr = tsv(f"{a.w}/net_reads.tsv")
    nets_rows = {r["family"]: r for r in tsv(f"{a.w}/nets.tsv")}
    out = {}
    bl = batch_logs(a.w, "stage")
    print("batch\tfamilies\tElapsed s\tstage clock s (nets / clusters / genome)\tpeak RSS GB\tflagged (families)")
    for g, v in bl.items():
        d = v["done"]
        print(f"{g}\t{int(d['families'])}\t{v['elapsed']:.1f}\t{d['total']:.1f} ({d['nets']:.1f} / {d['clusters_s']:.1f} / {d['genome']:.1f})\t{v['rss_gb']:.1f}\t"
              f"{int(d['flagged'])} ({int(d['flagged_families'])})")
    tot_e = sum(v["elapsed"] for v in bl.values())
    print(f"sum\t{sum(int(v['done']['families']) for v in bl.values())}\t{tot_e:.1f} = {tot_e / 60:.1f} min\t"
          f"{sum(v['done']['total'] for v in bl.values()):.1f}\t{max(v['rss_gb'] for v in bl.values()):.1f} (max)")
    out["batches"] = bl
    tot = collections.Counter()
    for r in fams.values():
        for k in ("n_net", "n_used", "n_clusters", "n_in_reference", "n_linked", "n_new", "n_candidates", "n_flagged"):
            tot[k] += int(r[k])
    out["stage"] = dict(tot)
    print("stage totals:", dict(tot), f"| families with >= 1 flagged candidate {sum(1 for r in fams.values() if int(r['n_flagged']) > 0)}/{len(fams)}")
    counters = {}
    for p in sorted(glob.glob(f"{a.w}/logs/stage_g*.log")):
        for m in CLUSTER_LOG.finditer(open(p).read()):
            counters[m.group(1)] = dict(zip(CLUSTER_LOG_KEYS, map(int, m.groups()[1:])))
    ct = collections.Counter()
    for v in counters.values():
        ct.update(v)
    fam_with = lambda k: sorted((f for f, v in counters.items() if v[k] > 0), key=name_key)
    out["counters"] = dict(total=dict(ct), per_family=counters)
    print(f"stage log counters over {len(counters)} families (phase 1): read clusters >= --min-cluster {ct['read_clusters']}; empty consensus dropped "
          f"{ct['empty']}; refinement split off {ct['split_off']} reads into {ct['split_clusters']} new clusters, {ct['fell_under']} clusters fell under "
          f"--min-cluster, kept sets re-templated {ct['retemplated']} (families {fam_with('retemplated')}); significance merge absorbed {ct['absorbed']} "
          f"clusters in {ct['rounds']} rounds, absorptions undone {ct['undone']}; fallback templates (no eligible member): longest with an aligned "
          f"partner {ct['longest_aligned']}, longest overall {ct['longest_unaligned']}; final clusters {ct['final']}")
    print("pass B per batch (stage logs): batch\tunmapped>=300\tpoorly placed>=300 (short)\taligned u / p\tattributed u / p\tjoined")
    pbt = collections.Counter()
    for g, v in bl.items():
        b = v["pass_b"]
        pbt.update(b)
        print(f"  {g}\t{b['unmapped']}\t{b['poorly_placed']} ({b['poorly_placed_short']})\t{b['aligned_unmapped']} / {b['aligned_poorly_placed']}\t"
              f"{b['attributed_unmapped']} / {b['attributed_poorly_placed']}\t{b['joined']}")
    print(f"  sum\t{pbt['unmapped']}\t{pbt['poorly_placed']} ({pbt['poorly_placed_short']})\t{pbt['aligned_unmapped']} / {pbt['aligned_poorly_placed']}\t"
          f"{pbt['attributed_unmapped']} / {pbt['attributed_poorly_placed']}\t{pbt['joined']}")
    out["pass_b_total"] = dict(pbt)
    # the joined reads (a read entering a net of its batch by pass B) and what they become
    jn = [r for r in netr if r["via"] == "B"]
    fate = collections.Counter()
    for r in jn:
        f, n = r["family"], r["read"]
        right = n in lab and lab[n]["family"] == f
        cl = cl_of[f].get(n)
        if r["used"] != "1":
            where = "over the cap"
        elif cl is None:
            where = "no reported cluster (in-reference / unreported)"
        elif clus[cl]["fate"] == "linked":
            where = "linked cluster"
        else:
            k = cands.get(clus[cl]["candidate"])
            where = "flagged candidate" if k and k["flagged"] == "1" else "candidate below support"
        fate[(r["class"], "right" if right else "wrong", where, calls[n]["C"] if n in calls else "-")] += 1
    out["joined"] = {"|".join(k): v for k, v in sorted(fate.items())}
    def by(i):
        c = collections.Counter()
        for k, v in fate.items():
            c[k[i]] += v
        return c
    print(f"joined reads (in a net of their batch via pass B): {len(jn)} joins of {len({r['read'] for r in jn})} reads; by class "
          f"{dict(collections.Counter(r['class'] for r in jn))}; right family {sum(v for k, v in fate.items() if k[1] == 'right')}; where they end "
          f"{dict(by(2))}; their arm-C class {dict(by(3))}")
    for k, v in sorted(fate.items(), key=lambda kv: -kv[1]):
        print(f"    {v:5d}  {k[0]:13s} {k[1]:5s}  {k[2]:48s} arm C: {k[3]}")
    multi = collections.Counter(r["read"] for r in netr)
    print(f"reads in nets of more than one family/batch: {sum(1 for v in multi.values() if v > 1)} (ruling R18)")
    # per family
    cpf = collections.defaultdict(list)
    for r in klass.values():
        cpf[r["family"]].append(r)
    fm = collections.Counter(c["family"] for c in calls.values() if c["C"] == "false_move")
    oc = collections.Counter(c["family"] for c in calls.values() if c["C"] == "own_candidate")
    print("family\tbatch\tnet (used)\tjoined\tclusters ref / linked / new\tcandidates (flagged)\tflagged: class/label\tC: own-candidate / false moves")
    table = []
    for f in sorted(fams, key=name_key):
        r = fams[f]
        fl = sorted(cpf.get(f, []), key=lambda x: name_key(x["candidate"]))
        t = dict(family=f, batch=nets_rows[f]["batch"], n_net=int(r["n_net"]), n_used=int(r["n_used"]), joined=int(nets_rows[f]["via_B"]),
                 in_ref=int(r["n_in_reference"]), linked=int(r["n_linked"]), new=int(r["n_new"]), candidates=int(r["n_candidates"]),
                 flagged=int(r["n_flagged"]), flags=[(x["class"], x["source"], int(x["n_reads"])) for x in fl], own_candidate=oc[f], false_moves=fm[f])
        table.append(t)
        print(f"{f}\t{t['batch']}\t{t['n_net']} ({t['n_used']})\t{t['joined']}\t{t['in_ref']} / {t['linked']} / {t['new']}\t{t['candidates']} ({t['flagged']})\t"
              f"{' '.join(c.split('_')[0] + '/' + s_ + '/' + str(nr) for c, s_, nr in t['flags']) or '-'}\t{t['own_candidate']} / {t['false_moves']}")
    out["per_family"] = table
    json.dump(out, open(f"{a.w}/report.json", "w"), indent=0)


def wplan(a):
    """ruling R23's whole-BAM batches: the families of --copies in table order, --batch-size per line -> --batch-file"""
    order = list(dict.fromkeys(r["family_id"] for r in tsv(a.copies)))
    with open(a.batch_file, "w") as o:
        for i in range(0, len(order), a.batch_size):
            o.write(",".join(order[i:i + a.batch_size]) + "\n")
    print(f"families {len(order)} -> {(len(order) + a.batch_size - 1) // a.batch_size} batches of <= {a.batch_size} ({a.batch_file})")


def wreport(a):
    """ruling R23: the whole-BAM batches' cost from their logs (logs/wstage_g<g>.log): Elapsed, peak RSS, the pass A / pass B counts (the
    attribution set: unmapped / poorly placed written; the targets = this batch's net reads by pass A + every copy record), flagged"""
    bl = batch_logs(a.w, "wstage")
    n_copy_records = sum(1 for ln in open(a.copies_fa) if ln.startswith(">"))
    print("batch\tElapsed s\tpeak RSS GB\tpass A reads (secondary records)\tattribution set: unmapped / poorly placed (short)\taligned\tattributed\t"
          "joined\ttargets (net reads + copies)\tflagged (families)\tstage clock (nets / clusters / genome)")
    for g, v in bl.items():
        b, d, pa = v["pass_b"], v["done"], v["pass_a"]
        fams = {r["family"]: r for r in tsv(f"{a.w}/cand_g{g}.families.tsv")} if os.path.exists(f"{a.w}/cand_g{g}.families.tsv") else {}
        n_targets = sum(int(r["n_net"]) for r in fams.values()) - b["joined"] if (fams and b) else None
        print(f"{g}\t{v['elapsed']}\t{v['rss_gb']}\t{pa and pa['primary']} ({pa and pa['secondary']})\t"
              f"{b and b['unmapped']} / {b and b['poorly_placed']} ({b and b['poorly_placed_short']})\t{b and b['aligned']}\t{b and b['attributed']}\t"
              f"{b and b['joined']}\t{n_targets} + {n_copy_records}\t{d and int(d['flagged'])} ({d and int(d['flagged_families'])})\t"
              f"{d and d['total']} ({d and d['nets']} / {d and d['clusters_s']} / {d and d['genome']})" + ("\tSTOPPED (did not finish)" if v["stopped"] else ""))
    done = [v for v in bl.values() if v["done"]]
    if done:
        print(f"sum over {len(done)} finished batches: Elapsed {sum(v['elapsed'] for v in done):.1f} s = {sum(v['elapsed'] for v in done) / 60:.1f} min; "
              f"peak RSS max {max(v['rss_gb'] for v in done):.1f} GB; flagged {sum(int(v['done']['flagged']) for v in done)}")
    json.dump(bl, open(f"{a.w}/wreport.json", "w"), indent=0)


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("cmd", choices=["nets", "plan", "concat", "contigs", "label", "keep", "report", "comparator", "decompose", "classify", "cscore",
                                    "creport", "wplan", "wreport"])
    ap.add_argument("--w", required=True, help="the run's work dir")
    ap.add_argument("--linktest", default="/mnt/linuxdisk/tmp/rna_allele/linktest", help="Amendment 7's work dir")
    ap.add_argument("--prefix", default="cand")
    ap.add_argument("--batches", type=int, default=7, help="plan: the number of NEW batches")
    ap.add_argument("--keep", type=int, default=0, help="plan: keep the first N lines of the batch file (batches already run)")
    ap.add_argument("--batch-file", default=None, help="plan: the batch file (default <w>/batches.txt)")
    ap.add_argument("--nets", default=None, help="dir with nets.tsv / net_reads.tsv (default --w)")
    ap.add_argument("--copies", default=None, help="nets: the stage's --copies table (default <w>/A12.copies.tsv)")
    ap.add_argument("--copies-fa", default=None, help="nets: the stage's --copies-fa (default <w>/A12.copies.fa)")
    ap.add_argument("--bam", default=None, help="nets: the stage's --bam (default <linktest>/R.bam; Amendment 14: control/R0.bam)")
    ap.add_argument("--all-copies", action="store_true", help="label: every copy is in the reference (Amendment 14): the masked copy is S:<copy> too")
    ap.add_argument("--control", default="/mnt/linuxdisk/tmp/rna_allele/control", help="Amendment 9's work dir (R0.bam, copies_lift.tsv, its candidates)")
    ap.add_argument("--a13", default="/mnt/linuxdisk/tmp/rna_allele/a13", help="classify: the A13 run's work dir (its survivor-derived flags)")
    ap.add_argument("--batch-size", type=int, default=50, help="wplan: families per whole-BAM batch")
    a = ap.parse_args(argv)
    a.nets = a.nets or a.w
    a.bam = a.bam or f"{a.linktest}/R.bam"
    a.copies = a.copies or f"{a.w}/A12.copies.tsv"
    a.copies_fa = a.copies_fa or f"{a.w}/A12.copies.fa"
    globals()[a.cmd](a)


if __name__ == "__main__":
    main()
