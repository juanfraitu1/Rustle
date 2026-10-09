#!/usr/bin/env python3
"""Amendment 41: maternal-assembly truth for the real-read recall test (docs/PREREG_unmapped_rescue_2026-10-08.md). Miniforge python.

Pure helpers (tested in test_mat_truth.py): cs_divergence, mask, assembly_version, group_loci, verdict.
Steps (resumable, each one call): see main()."""
import re

CS = re.compile(r"(:\d+|\*[a-z][a-z]|[+-][a-z]+|~[a-z]{2}\d+[a-z]{2})")
CIG = re.compile(r"(\d+)([MIDNSHP=X])")
DELTA = 0.00958


def cs_divergence(tstart, cs):
    """-> ([(target position, divergent bases)], target end): substitutions (1 base), deletions from the target (their length, at their start), insertions into the
    query (their length, at the target position where they occur); the cs walks the target forward from tstart"""
    pos, ev = tstart, []
    for tok in CS.findall(cs):
        c = tok[0]
        if c == ":":
            pos += int(tok[1:])
        elif c == "*":
            ev.append((pos, 1))
            pos += 1
        elif c == "-":
            ev.append((pos, len(tok) - 1))
            pos += len(tok) - 1
        elif c == "+":
            ev.append((pos, len(tok) - 1))
        else:  # ~ intron (not produced by asm5)
            pos += int(re.findall(r"\d+", tok)[0])
    return ev, pos


def _merge(iv):
    out = []
    for a, b in sorted(iv):
        if out and a <= out[-1][1]:
            out[-1][1] = max(out[-1][1], b)
        else:
            out.append([a, b])
    return [(a, b) for a, b in out]


def mask(tlen, covered, events, window=500, max_div=DELTA / 2):
    """target intervals with no counterpart within max_div: not covered by any record, or in a window whose divergent bases exceed max_div x window"""
    cov = _merge([(max(0, a), min(tlen, b)) for a, b in covered if a < b])
    out, prev = [], 0
    for a, b in cov:
        if a > prev:
            out.append((prev, a))
        prev = max(prev, b)
    if prev < tlen:
        out.append((prev, tlen))
    per = {}
    for p, n in events:
        per[p // window] = per.get(p // window, 0) + n
    for w, n in per.items():
        if n > max_div * window:
            out.append((w * window, min(tlen, (w + 1) * window)))
    return _merge(out)


def assembly_version(fetch, pos0, cigar):
    """the reference bases under an alignment's blocks (M = X D), introns (N) skipped, read insertions and clips ignored; fetch(a, b) -> reference[a:b]"""
    out, pos = [], pos0
    for n, op in CIG.findall(cigar):
        n = int(n)
        if op in "M=XD":
            out.append(fetch(pos, pos + n))
            pos += n
        elif op == "N":
            pos += n
    return "".join(out)


def group_loci(reads, gap=1000):
    """reads = [(contig, start, end, name)] -> [dict(contig, start, end, names)], joined when a read starts within `gap` of the current locus end"""
    loci = []
    for c, s, e, n in sorted(reads):
        if loci and loci[-1]["contig"] == c and s <= loci[-1]["end"] + gap:
            loci[-1]["end"] = max(loci[-1]["end"], e)
            loci[-1]["names"].append(n)
        else:
            loci.append(dict(contig=c, start=s, end=e, names=[n]))
    return loci


def verdict(statuses):
    """TRUE when most placed reads (mat / pat / shared) are mat- or pat-specific; WRONG otherwise (including no placed read)"""
    placed = [s for s in statuses if s in ("mat", "pat", "shared")]
    spec = sum(s in ("mat", "pat") for s in placed)
    return "TRUE" if placed and spec > len(placed) / 2 else "WRONG"


# ---------------------------------------------------------------- steps
import json  # noqa: E402
import os  # noqa: E402
import subprocess  # noqa: E402
import sys  # noqa: E402

RA = "/mnt/linuxdisk/tmp/rna_allele"
HAPS = "/mnt/linuxdisk/home/juanfraitu/gorilla_haps"
PRI = "/mnt/linuxdisk/home/juanfraitu/_from_wsl/winloci_scratch/GGO.fasta"
FIB = "/mnt/linuxdisk/home/juanfraitu/fibroblasts/GCA_029281585.2_flnc_mm.bam"
OUT = "/mnt/linuxdisk/tmp/o3_rescue/mattruth"


def chrmap():
    rows = [ln.rstrip("\n").split("\t") for ln in open(f"{RA}/chrmap.tsv")][1:]
    return [dict(pri=r[0], chrom=r[1], same=r[2], bhap=r[5] if len(r) > 5 else "", bname=r[6] if len(r) > 6 else "") for r in rows]


def step_mask(hap="mat"):
    """mask.<hap>.bed: <hap> bases with no primary counterpart within delta/2 (chromosomes whose primary is the other haplotype) + every unplaced <hap> contig"""
    os.makedirs(OUT, exist_ok=True)
    lens = {ln.split("\t")[0]: int(ln.split("\t")[1]) for ln in open(f"{HAPS}/{hap}.fa.fai")}
    placed = {r["bname"] for r in chrmap() if r["bhap"] == hap} | {r_ for r in chrmap() for r_ in [r.get("same_name")] if r_}
    rows, tot = [], 0
    for r in chrmap():
        if r["bhap"] != hap:
            continue
        covered, events = [], []
        for ln in open(f"{RA}/out/chr{r['chrom']}.paf"):
            f = ln.rstrip("\n").split("\t")
            if "tp:A:P" not in f[12:]:
                continue
            ts, te = int(f[7]), int(f[8])
            covered.append((ts, te))
            cs = next(t[5:] for t in f[12:] if t.startswith("cs:Z:"))
            ev, _ = cs_divergence(ts, cs)
            events += ev
        for a, b in mask(lens[r["bname"]], covered, events):
            rows.append((r["bname"], a, b, f"chr{r['chrom']}"))
            tot += b - a
    cm = {ln.rstrip("\n").split("\t")[6] for ln in open(f"{RA}/chrmap.tsv") if ln.count("\t") >= 6} | {ln.split("\t")[3] for ln in open(f"{RA}/chrmap.tsv")}
    un = [c for c in lens if c not in cm]
    for c in un:
        rows.append((c, 0, lens[c], "unplaced"))
    with open(f"{OUT}/mask.{hap}.bed", "w") as o:
        for c, a, b, lab in rows:
            o.write(f"{c}\t{a}\t{b}\t{lab}\n")
    print(f"mask.{hap}: {len(rows) - len(un)} intervals, {tot:,} bp on the placed chromosomes; {len(un)} unplaced contigs, {sum(lens[c] for c in un):,} bp")


RCT = str.maketrans("ACGTNacgtn", "TGCANtgcan")


def load_best(paf_files):
    """read -> PAF fields of its primary-chain record with the most matching bases"""
    best = {}
    for fn in paf_files:
        for ln in open(fn):
            f = ln.rstrip("\n").split("\t")
            if "tp:A:P" not in f[12:]:
                continue
            if f[0] not in best or int(f[9]) > int(best[f[0]][9]):
                best[f[0]] = f
    return best


def tag(f, key):
    return next((t[5:] for t in f[12:] if t.startswith(key + ":")), None)


def testable_contigs(hap="mat"):
    """contigs of <hap> that are NOT identical to a primary chromosome: the partner chromosomes of the other haplotype's primary ones + unplaced contigs"""
    lens = [ln.split("\t")[0] for ln in open(f"{HAPS}/{hap}.fa.fai")]
    same = {r[3] for r in (ln.rstrip("\n").split("\t") for ln in list(open(f"{RA}/chrmap.tsv"))[1:]) if r[2] == hap}
    return {c for c in lens if c not in same}


def step_versions(which, hap="mat"):
    """versions.<which>.<hap>.fa: the <hap> genomic bases under each clean read's alignment blocks, in the read's orientation, for reads on testable contigs"""
    import glob
    import pysam
    sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
    from o3_maternal import common as C
    al = C.alias()
    best = load_best(sorted(glob.glob(f"{OUT}/{which}_{hap}/chunk*.paf")))
    for f in best.values():
        f[5] = C.accession(f[5], al)       # index names chrN_<hap>_hsaX -> the FASTA's accessions
    fa = pysam.FastaFile(f"{HAPS}/{hap}.fa")
    test = testable_contigs(hap)
    n_clean = n_test = 0
    with open(f"{OUT}/versions.{which}.{hap}.fa", "w") as o, open(f"{OUT}/versions.{which}.{hap}.tsv", "w") as t:
        t.write("read\tcontig\tstart\tend\tstrand\tde\tmapq\n")
        for r, f in best.items():
            de = float(tag(f, "de"))
            if de > DELTA:
                continue
            n_clean += 1
            if f[5] not in test:
                continue
            n_test += 1
            seq = assembly_version(lambda x, y, c=f[5]: fa.fetch(c, x, y).upper(), int(f[7]), tag(f, "cg"))
            if f[4] == "-":
                seq = seq.translate(RCT)[::-1]
            o.write(f">{r}\n{seq}\n")
            t.write(f"{r}\t{f[5]}\t{f[7]}\t{f[8]}\t{f[4]}\t{de}\t{f[11]}\n")
    json.dump(dict(mapped=len(best), clean=n_clean, testable=n_test), open(f"{OUT}/versions.{which}.{hap}.json", "w"))
    print(f"{which}/{hap}: mapped {len(best)}, clean (de <= {DELTA}) {n_clean}, on testable contigs {n_test}")


def step_score(which, hap="mat"):
    """truth.<which>.<hap>.tsv: identity x coverage (with terminal-exon rescue, gate off) of every assembly version on the primary; specific = score < 1 - delta"""
    import glob
    import pysam
    sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
    import discover as D
    import seeds as SD
    vers = SD.read_fa(f"{OUT}/versions.{which}.{hap}.fa")
    rec = load_best(sorted(glob.glob(f"{OUT}/vpri_{which}_{hap}/chunk*.paf")))
    pri = pysam.FastaFile(D.PRIMARY_FA)
    meta = {ln.split("\t")[0]: ln.rstrip("\n").split("\t") for ln in list(open(f"{OUT}/versions.{which}.{hap}.tsv"))[1:]}
    n_spec = 0
    with open(f"{OUT}/truth.{which}.{hap}.tsv", "w") as o:
        o.write("read\tcontig\tstart\tend\tscore\tspecific\n")
        for r, seq in vers.items():
            if not seq:
                continue
            sc = D.idcov_with_rescue(seq, rec[r], pri, D.ident, False, f"{OUT}/tmp_rescue") if r in rec else 0.0
            spec = sc < 1 - DELTA
            n_spec += spec
            m = meta[r]
            o.write(f"{r}\t{m[1]}\t{m[2]}\t{m[3]}\t{sc:.5f}\t{int(spec)}\n")
    print(f"{which}/{hap}: {len(vers)} assembly versions scored; {hap}-specific (score < {1 - DELTA:.5f}) {n_spec}")


def read_status(r, best_mat, spec_mat, best_pat, spec_pat):
    """status of a read from its best alignments on the two haplotype assemblies: the haplotype where it aligns with the lower de decides (ties: not specific);
    mat / pat = clean there and its assembly version has no primary counterpart within delta; shared = clean there with a counterpart; unplaced = clean on neither"""
    cand = []
    for hap, best, spec in (("mat", best_mat, spec_mat), ("pat", best_pat, spec_pat)):
        f = best.get(r)
        if f is None:
            continue
        de = float(tag(f, "de"))
        if de <= DELTA:
            cand.append((de, 0 if spec.get(r) else 1, hap, bool(spec.get(r))))
    if not cand:
        return "unplaced"
    cand.sort(key=lambda x: (x[0], -x[1]))           # lowest de; on a tie the non-specific placement first
    de, _, hap, sp = cand[0]
    return hap if sp else "shared"


def load_spec(fn):
    out = {}
    if os.path.exists(fn):
        for ln in list(open(fn))[1:]:
            f = ln.rstrip("\n").split("\t")
            out[f[0]] = f[5] == "1"
    return out


def step_loci():
    """loci.json: the net's mat-specific reads grouped by maternal position"""
    reads = []
    for ln in list(open(f"{OUT}/truth.net.mat.tsv"))[1:]:
        f = ln.rstrip("\n").split("\t")
        if f[5] == "1":
            reads.append((f[1], int(f[2]), int(f[3]), f[0]))
    loci = group_loci(reads)
    json.dump(loci, open(f"{OUT}/loci.json", "w"))
    n3 = [l for l in loci if len(l["names"]) >= 3]
    print(f"truth loci (net-reachable): {len(loci)} with >= 1 mat-specific read ({len(reads)} reads); {len(n3)} with >= 3 ({sum(len(l['names']) for l in n3)} reads)")
    ston = [l for l in loci if l["contig"] == "CM054594.2" and l["start"] < 96_100_000 and l["end"] > 95_900_000]
    print("STON1-GTF2A1L positive control:", [(l["contig"], l["start"], l["end"], len(l["names"])) for l in ston])
    big = sorted(loci, key=lambda l: -len(l["names"]))[:10]
    print("largest loci:", [(l["contig"], l["start"], l["end"], len(l["names"])) for l in big])


def step_outside():
    """the 1% sample of reads outside the net: how many mat-specific, how many inside a net truth locus, how many distinct loci outside"""
    loci = json.load(open(f"{OUT}/loci.json"))
    spec = [ln.rstrip("\n").split("\t") for ln in list(open(f"{OUT}/truth.sample.mat.tsv"))[1:]]
    sp = [(f[1], int(f[2]), int(f[3]), f[0]) for f in spec if f[5] == "1"]
    inside = [r for r in sp if any(l["contig"] == r[0] and r[1] <= l["end"] + 1000 and r[2] >= l["start"] - 1000 for l in loci)]
    out = [r for r in sp if r not in inside]
    out_loci = group_loci(out)
    meta = json.load(open(f"{OUT}/versions.sample.mat.json"))
    print(f"sample outside the net: {meta['clean']} clean reads, {meta['testable']} on testable contigs, {len(sp)} mat-specific; inside a net truth locus {len(inside)}; "
          f"outside every net locus {len(out)} in {len(out_loci)} distinct loci (each sampled read stands for ~100)")
    json.dump(dict(specific=len(sp), inside=len(inside), outside=len(out), outside_loci=len(out_loci)), open(f"{OUT}/outside.json", "w"))


GOOD = ("PRESENT", "ALLELE-LIKE")
TAG = os.environ.get("MATTRUTH_TAG", "")          # "" = Amendment 41's seed-round run; "_binned" = Amendment 42
RUN = f"{OUT}/net_run{TAG}"


def run_clusters():
    import csv
    import collections
    cl = collections.defaultdict(list)
    for r in csv.DictReader(open(f"{RUN}/clusters.tsv"), delimiter="\t"):
        cl[r["cluster"]].append(r["read"])
    rows = json.load(open(f"{RUN}/classes.json"))
    return cl, rows


def flagged(rows):
    return [r for r in rows if r["cls"] not in GOOD and r["cls"] != "UNSUPPORTED"]


def step_flag_reads():
    """flag_reads.fa: every read of a flagged cluster (to place them on the paternal assembly)"""
    sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
    import seeds as SD
    cl, rows = run_clusters()
    names = sorted({x for r in flagged(rows) for x in cl[r["key"]]})
    seqs = SD.read_fa(f"{OUT}/net_all.fa")
    SD.write_fa(f"{OUT}/flag{TAG}_reads.fa", {n: seqs[n] for n in names}, names)
    print(f"{len(flagged(rows))} flagged clusters, {len(names)} reads")


def step_eval():
    import collections
    import glob
    sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
    cl, rows = run_clusters()
    loci = json.load(open(f"{OUT}/loci.json"))
    locus_of = {n: i for i, l in enumerate(loci) for n in l["names"]}
    best_mat = load_best(sorted(glob.glob(f"{OUT}/net_mat/chunk*.paf")))
    best_pat = load_best(sorted(glob.glob(f"{OUT}/flag{TAG}_pat/chunk*.paf")))
    spec_mat = load_spec(f"{OUT}/truth.net.mat.tsv")
    spec_pat = load_spec(f"{OUT}/truth.flag{TAG}.pat.tsv")
    flags = flagged(rows)
    verdicts, recovered = [], set()
    for r in flags:
        st = [read_status(x, best_mat, spec_mat, best_pat, spec_pat) for x in cl[r["key"]]]
        v = verdict(st)
        c = collections.Counter(st)
        mats = [locus_of[x] for x in cl[r["key"]] if x in locus_of]
        loc = collections.Counter(mats).most_common(1)[0][0] if mats else None
        if v == "TRUE" and loc is not None:
            recovered.add(loc)
        verdicts.append(dict(k=r["k"], cls=r["cls"], reads=r["reads"], length=r["length"], R=r["R"], verdict=v, status=dict(c), locus=loc))
    n_true = sum(v["verdict"] == "TRUE" for v in verdicts)
    by_cls = collections.Counter((v["cls"], v["verdict"]) for v in verdicts)
    all_cls = collections.Counter(r["cls"] for r in rows)
    print(f"pipeline: {len(rows)} clusters ({sum(r['reads'] for r in rows)} reads) of the 526,772-read net; classes {dict(all_cls)}")
    print(f"FLAGGED clusters: {len(flags)}; TRUE {n_true}; WRONG {len(flags) - n_true} ({(len(flags) - n_true) / max(1, len(flags)):.1%}); by class and verdict {dict(by_cls)}")
    # recall
    in_cluster = {x: k for k, rs in cl.items() for x in rs}
    cls_of = {r["key"]: r["cls"] for r in rows}
    flagkeys = {r["key"] for r in flags}
    for floor in (1, 3):
        sel = [i for i, l in enumerate(loci) if len(l["names"]) >= floor]
        rec = [i for i in sel if i in recovered]
        why = collections.Counter()
        for i in sel:
            if i in recovered:
                continue
            names = loci[i]["names"]
            ks = [in_cluster.get(x) for x in names]
            if sum(k is None for k in ks) > len(ks) / 2:
                why["not clustered"] += 1
                continue
            cs = collections.Counter(cls_of[k] for k in ks if k is not None)
            top = cs.most_common(1)[0][0]
            if top in GOOD:
                why["clustered, not flagged (consensus held by the primary)"] += 1
            elif top == "UNSUPPORTED":
                why["clustered, consensus unsupported"] += 1
            else:
                why["flagged, but the cluster is WRONG or dominated by another locus"] += 1
        print(f"recall over truth loci with >= {floor} mat-specific reads: {len(rec)} of {len(sel)} = {len(rec) / max(1, len(sel)):.3f}; misses by first failing step: {dict(why)}")
    ston = [i for i, l in enumerate(loci) if l["contig"] == "CM054594.2" and l["start"] < 96_100_000 and l["end"] > 95_900_000]
    print("STON1-GTF2A1L positive control recovered:", [(i, i in recovered) for i in ston])
    json.dump(dict(verdicts=verdicts, recovered=sorted(recovered)), open(f"{OUT}/eval{TAG}.json", "w"))


if __name__ == "__main__":
    cmd = sys.argv[1]
    if cmd == "mask":
        step_mask(sys.argv[2] if len(sys.argv) > 2 else "mat")
    elif cmd == "versions":
        step_versions(sys.argv[2], sys.argv[3] if len(sys.argv) > 3 else "mat")
    elif cmd == "score":
        step_score(sys.argv[2], sys.argv[3] if len(sys.argv) > 3 else "mat")
    elif cmd == "loci":
        step_loci()
    elif cmd == "outside":
        step_outside()
    elif cmd == "flagreads":
        step_flag_reads()
    elif cmd == "eval":
        step_eval()
