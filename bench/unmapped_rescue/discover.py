#!/usr/bin/env python3
"""Discovery of reference-absent sequence from reads that do not align to the reference (docs/PREREG_unmapped_rescue_2026-10-08.md Amendment 28). Miniforge python.

    discover.py <name> <reads.fa>      resumable (exit 75 = run again); writes W/discover_<name>/{clusters.tsv, cons.fa, classes.json}
Cluster (adopted rule: --proper edge, star step), abPOA consensus, align each consensus to the primary (R) and to the mother's (Tm) and the father's (Tp) assembly,
classify CONFIRMED / NOVEL / ALLELE-LIKE / PRESENT (flagmetric.discovery_class)."""
import collections
import csv
import functools
import itertools
import json
import os
import random
import subprocess
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.dirname(HERE))
import chain as CH  # noqa: E402
import flagmetric as FM  # noqa: E402
import graph as GR  # noqa: E402
import libsig  # noqa: E402
import consensus_support as CS  # noqa: E402
import partition as PT  # noqa: E402
import rescue as RS  # noqa: E402
import run_polish as RP  # noqa: E402
import seeds as SD  # noqa: E402
from o3_maternal import common as C  # noqa: E402

W = "/mnt/linuxdisk/tmp/o3_rescue"
M = "/mnt/linuxdisk/tmp/o3_mat"
DELTA = 0.00958
RCT = str.maketrans("ACGTacgt", "TGCAtgca")
PRIMARY_FA = "/mnt/linuxdisk/tmp/rustle_figures/liftoff/gorilla/genome.fa"
HAP_FA = "/mnt/linuxdisk/home/juanfraitu/gorilla_haps/{}.fa"
WINDOW = 5000
WD = "/mnt/linuxdisk/home/juanfraitu/winloci_data"
CHM13_IDX = "/mnt/linuxdisk/home/juanfraitu/npip_ladder/idx/target.splice.mmi"
CHM13_FA = "/mnt/linuxdisk/home/juanfraitu/npip_ladder/idx/target.fa"
PTR_FA = f"{WD}/GCF_028858775.2_NHGRI_mPanTro3-v2.0_pri_genomic.fna"
PPY_FA = f"{WD}/GCF_028885625.2_NHGRI_mPonPyg2-v2.0_pri_genomic.fna"
AD = f"{W}/animals"
TIME_BUDGET = 360   # seconds of alignment before the run pauses (exit 75) so one Bash call stays under ten minutes


def ident(n):
    return n


@functools.lru_cache(maxsize=1)
def _alias():
    return C.alias()


def hap_name(n):
    return C.accession(n, _alias())


GGO_R = (RP.GENOME_IDX, PRIMARY_FA, ident)
# Amendment 33. R = the animal's own reference (index or FASTA, FASTA for the end rescue, contig-name function); others = every other assembly on disk.
ANIMALS = {
    "ggo_testis": dict(bam=f"{WD}/GGO_mm.bam", R=GGO_R,
                       others={"Tm": (C.HAP_IDX.format("mat"), HAP_FA.format("mat"), hap_name), "Tp": (C.HAP_IDX.format("pat"), HAP_FA.format("pat"), hap_name)}),
    "a119b": dict(bam=f"{WD}/A119b.t2t.bam", R=(CHM13_IDX, CHM13_FA, ident), others={"GGO": GGO_R}),
    "ptr": dict(bam=f"{WD}/PTR_mm.bam", R=(PTR_FA, PTR_FA, ident), others={"HSA": (CHM13_IDX, CHM13_FA, ident), "GGO": GGO_R}),
    "ppy": dict(bam=f"{WD}/PPY_mm.bam", R=(PPY_FA, PPY_FA, ident), others={"HSA": (CHM13_IDX, CHM13_FA, ident), "GGO": GGO_R}),
}


def dense_reads(recs, bin_size=5000, min_reads=5):
    """names of the reads that share a (contig, position // bin_size) bin with at least min_reads reads; recs = (name, contig, pos)"""
    n = collections.Counter((c, p // bin_size) for _, c, p in recs)
    return [name for name, c, p in recs if n[(c, p // bin_size)] >= min_reads]


def sample_names(names, k, seed):
    names = sorted(names)
    return names if len(names) <= k else sorted(random.Random(seed).sample(names, k))


def class_counts(rows):
    by = collections.defaultdict(lambda: [0, 0])
    for r in rows:
        by[r["cls"]][0] += 1
        by[r["cls"]][1] += r["reads"]
    return by


def specific_fraction(rows):
    """share of the clustered reads that lie in PRESENT or ALLELE-LIKE clusters (the control bar); None with no clusters"""
    tot = sum(r["reads"] for r in rows)
    return sum(r["reads"] for r in rows if r["cls"] in ("PRESENT", "ALLELE-LIKE")) / tot if tot else None


def prepare(animal, n_windows=12, width=5_000_000, limit=200_000, control_size=959, seed=5):
    """the library gate (5' G signature of its clean mapped primaries) and the specificity control (959 clean reads from dense loci); one pass over windows of the longest contigs"""
    import pysam
    a = ANIMALS[animal]
    d = f"{AD}/{animal}"
    os.makedirs(d, exist_ok=True)
    h = pysam.AlignmentFile(a["bam"], "rb")
    ctgs = sorted(zip(h.references, h.lengths), key=lambda x: -x[1])[:n_windows]
    regs = [f"{c}:{max(0, L // 2 - width // 2) + 1}-{min(L, L // 2 + width // 2)}" for c, L in ctgs]
    proc = subprocess.Popen(["samtools", "view", "-F", "2308", "-q", "10", a["bam"], *regs], stdout=subprocess.PIPE, text=True)
    keep, n = [], [0]

    def lines():
        for ln in itertools.islice(proc.stdout, limit):
            n[0] += 1
            f = ln.rstrip("\n").split("\t")
            de = next((float(t[5:]) for t in f[11:] if t.startswith("de:f:")), 1.0)
            if de <= DELTA:
                keep.append((f[0], f[2], int(f[3]), f[9].translate(RCT)[::-1] if int(f[1]) & 16 else f[9]))
            yield ln
    sig = libsig.signature(lines())
    proc.terminate()
    ok, p = libsig.gate(sig)
    json.dump(dict(gate=ok, p=p, signature=sig, windows=regs, primaries_read=n[0]), open(f"{d}/gate.json", "w"), indent=1)
    dense = set(dense_reads([(x[0], x[1], x[2]) for x in keep]))
    pick = set(sample_names(dense, control_size, seed))
    seqs = {x[0]: x[3] for x in keep if x[0] in pick}
    SD.write_fa(f"{d}/control.fa", seqs, sorted(seqs))
    print(f"{animal}: {n[0]} primaries read in {len(regs)} windows, {len(keep)} clean (de <= {DELTA}), {len(dense)} in dense bins, control {len(seqs)} reads; "
          f"gate {'OPEN' if ok else 'closed'} (p {p:.3g}; pure clips {sig['pure']})")


def cluster(name, reads_fa):
    d = f"{W}/discover_{name}"
    os.makedirs(f"{d}/rounds", exist_ok=True)
    open(f"{d}/reads_path.txt", "w").write(reads_fa)
    seqs = SD.read_fa(reads_fa)
    lens = {n: len(s) for n, s in seqs.items()}
    try:
        comp = SD.run_rounds(lens, SD.minimap_map_fn(seqs, f"{d}/rounds", DELTA, 0.5, proper=True), n_seeds=2000, min_size=3, max_rounds=6)
    except SD.Pause as e:
        print(f"round {e} finished; run again")
        sys.exit(75)
    cl = {c: rs for c, rs in GR.clusters(comp, 3).items()}
    new = CH.refine(cl, seqs, CH.minimap_allvsall(f"{d}/chain_tmp"), DELTA)
    with open(f"{d}/clusters.tsv", "w") as o:
        o.write("read\tcluster\tsize\n")
        for c, rs in new.items():
            for r in rs:
                o.write(f"{r}\t{c}\t{len(rs)}\n")
    cons = {f"cl{c}|n={len(rs)}": PT.abpoa_consensus([seqs[r] for r in rs]) for c, rs in new.items()}
    SD.write_fa(f"{d}/cons.fa", cons, list(cons))
    print(f"{name}: {len(seqs)} reads, {len(cl)} components, {len(new)} clusters after the star step, {sum(len(v) for v in new.values())} reads clustered")
    return d, seqs, new


def best_records(paf):
    b = {}
    for ln in open(paf):
        f = ln.rstrip("\n").split("\t")
        k = f[0].split("|")[0]
        if k not in b or int(f[9]) > int(b[k][9]):
            b[k] = f
    return b


def idcov_with_rescue(cons, f, fasta, to_fasta_name, gate, tmp):
    """gate-aware identity x coverage of a consensus on one assembly, with the unaligned ends of the best record rescued inside the locus window (Amendment 29)"""
    qlen, qs, qe = int(f[1]), int(f[2]), int(f[3])
    matches, blk = int(f[9]), int(f[10])
    segs = RS.end_segments(qs, qe, qlen)
    iv, pm, pb = [], 0, 0
    if segs:
        chrom, ts, te = to_fasta_name(f[5]), int(f[7]), int(f[8])
        lo, hi = max(0, ts - WINDOW), min(fasta.get_reference_length(chrom), te + WINDOW)
        iv, pm, pb = RS.rescued(segs, RS.map_pieces(fasta.fetch(chrom, lo, hi).upper(), segs, cons, tmp))
    aligned, ident = RS.combine(qlen, (qs, qe, matches, blk), (iv, pm, pb))
    lead = min([qs] + [a for a, _ in iv])
    return FM.gate_aware_total(ident, aligned, qlen, lead, cons[:lead], gate)


def classify(name, confirmable=(), animal=None):
    """animal=None: the KB3781 setting of Amendment 28 (R = primary, Tm and Tp = the haplotype assemblies; CONFIRMED / NOVEL / ...). animal=<key of ANIMALS>: Amendment 33
    (R = the animal's own reference, E = the best score over the other assemblies; ELSEWHERE / NOVEL / ...)."""
    import pysam
    d = f"{W}/discover_{name}"
    t0 = time.time()
    if animal is None:
        gate = json.load(open(f"{W}/o3hap_mat/gate.json"))["gate"]
        al = C.alias()
        assemblies = {"R": (RP.GENOME_IDX, PRIMARY_FA, lambda n: n), "Tm": (C.HAP_IDX.format("mat"), HAP_FA.format("mat"), lambda n: C.accession(n, al)),
                      "Tp": (C.HAP_IDX.format("pat"), HAP_FA.format("pat"), lambda n: C.accession(n, al))}
    else:
        gate = json.load(open(f"{AD}/{animal}/gate.json"))["gate"]
        assemblies = {"R": ANIMALS[animal]["R"], **ANIMALS[animal]["others"]}
    cons = {n.split("|")[0]: s for n, s in RP.read_cons(f"{d}/cons.fa").items()}
    recs = {}
    for which, (idx, _fa, _nm) in assemblies.items():
        paf = f"{d}/cons.{which}.paf"
        if not os.path.exists(paf + ".done"):
            if time.time() - t0 > TIME_BUDGET:
                print(f"{which}: paused after {time.time() - t0:.0f} s; run again")
                sys.exit(75)
            subprocess.run(f"minimap2 -c -x splice:hq -uf -N 5 -t 4 {idx} {d}/cons.fa > {paf}", shell=True, check=True, stderr=subprocess.DEVNULL)
            open(paf + ".done", "w").write("ok")
        recs[which] = best_records(paf)
    fas = {w: pysam.FastaFile(a[1]) for w, a in assemblies.items()}
    cl = collections.defaultdict(list)
    for r in csv.DictReader(open(f"{d}/clusters.tsv"), delimiter="\t"):
        cl["cl" + r["cluster"]].append(r["read"])
    reads_fa = SD.read_fa(open(f"{d}/reads_path.txt").read().strip())
    rows = []
    for k, rs in cl.items():
        sc = {w: (idcov_with_rescue(cons[k], recs[w][k], fas[w], assemblies[w][2], gate, f"{d}/tmp_{w}") if k in recs[w] else None) for w in assemblies}
        med = CS.median_read_divergence(cons[k], {n: reads_fa[n] for n in sorted(rs)[:100]}, f"{d}/tmp_support")
        unsupported = med is None or med > DELTA
        if animal is None:
            cls = "UNSUPPORTED" if unsupported else FM.discovery_class(sc["R"], sc["Tm"], sc["Tp"])
        else:
            cls = "UNSUPPORTED" if unsupported else FM.elsewhere_class(sc["R"], [v for w, v in sc.items() if w != "R"])
        hits = {w: (r[5], int(r[7]), int(r[8])) for w, rr in recs.items() for r in [rr.get(k)] if r}
        rows.append(dict(k=k, reads=len(rs), length=len(cons[k]), scores=sc, R=sc["R"], Tm=sc.get("Tm"), Tp=sc.get("Tp"), median_read_divergence=med, cls=cls, hits=hits,
                         mat_hit=hits.get("Tm")))
    json.dump(rows, open(f"{d}/classes.json", "w"), indent=1)
    by = class_counts(rows)
    print(f"{name}: gate {'OPEN' if gate else 'closed'}; {len(rows)} clusters, {sum(r['reads'] for r in rows)} reads; PRESENT or ALLELE-LIKE share of clustered reads "
          f"{specific_fraction(rows) if specific_fraction(rows) is None else round(specific_fraction(rows), 4)}")
    order = ("CONFIRMED", "NOVEL", "DIVERGED", "UNSUPPORTED", "ALLELE-LIKE", "PRESENT") if animal is None else ("ELSEWHERE", "NOVEL", "DIVERGED", "UNSUPPORTED", "ALLELE-LIKE", "PRESENT")
    for c in order:
        print(f"  {c:12s} clusters {by[c][0]:4d}  reads {by[c][1]:5d}")
    if confirmable:
        cls_of = {r_: next(x["cls"] for x in rows if x["k"] == k) for k, rs in cl.items() for r_ in rs}
        inn = sum(cls_of.get(r) == "CONFIRMED" for r in confirmable)
        print(f"  confirmable reads (map on the mother's or the father's assembly): {len(confirmable)}; in CONFIRMED clusters {inn} ({inn / len(confirmable):.0%}); in any cluster "
              f"{sum(r in cls_of for r in confirmable)}")
    big = sorted(rows, key=lambda x: -x["reads"])[:5]
    f3 = lambda v: None if v is None else round(v, 4)
    print("  largest clusters (reads, length, class, scores by assembly, median read divergence):", [(x["reads"], x["length"], x["cls"], {w: f3(v) for w, v in x["scores"].items()}, f3(x["median_read_divergence"])) for x in big])
    for c in order[:4]:
        sel = sorted((x for x in rows if x["cls"] == c), key=lambda x: -x["reads"])
        print(f"  {c} clusters: {len(sel)}; reads {[x['reads'] for x in sel][:15]}; lengths {[x['length'] for x in sel][:15]}")


if __name__ == "__main__":
    if sys.argv[1] == "prepare":
        prepare(sys.argv[2])
        sys.exit(0)
    name, reads_fa = sys.argv[1], sys.argv[2]
    animal = sys.argv[3] if len(sys.argv) > 3 else None
    conf = ()
    if name == "unm":
        import pysam
        s = set()
        for h in ("mat", "pat"):
            for rd in pysam.AlignmentFile(f"{M}/map/R_unm.{h}.bam", "rb").fetch(until_eof=True):
                if not rd.is_unmapped and not rd.is_secondary and not rd.is_supplementary:
                    s.add(rd.query_name)
        conf = sorted(s)
    if not os.path.exists(f"{W}/discover_{name}/clusters.tsv"):
        os.makedirs(f"{W}/discover_{name}", exist_ok=True)
        cluster(name, reads_fa)
    classify(name, conf, animal)
