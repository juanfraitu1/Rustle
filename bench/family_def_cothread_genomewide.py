#!/usr/bin/env python3
"""family_def_cothread_genomewide.py — deploy junction-POSITION co-threading at genome scale.

The panel prototype (family_def_backbone_vgspine.py) showed junction-position co-threading cleanly
separates real families (RABL2 0.89, APOBEC3 0.67, RFPL 0.50) from domain-sharers/retrocopies (0.0)
on n=9. The open question (verifier challenge): does a POSITION-based metric cut the multi-exon
DOMAIN BRIDGES of the over-merge while KEEPING length/count-drifted real paralogs (e.g. the held-out
GGT1~GGTLC2, 15 vs 4 introns) connected? Count/length concordance could not (it destroys real
paralogs). Position might -- a real partial-dup shares a CONTIGUOUS, HOMOLOGOUS-POSITION run of
junctions; a domain bridge shares at most a scattered sub-path.

This harness does NOT hard-code one metric. For every candidate edge it aligns the two spliced cDNAs
(minimap2, batched per reference) and emits rich position features so the discriminating form can be
chosen empirically against the held-out positives:
  ref_j / mem_j        : junction counts (ref = richer backbone)
  co                   : # ref junctions co-threaded by the member (projected within TOL bp)
  frac_ref = co/ref_j  : share of the RICHER backbone (domain-sharers low; whole-dup high)
  frac_mem = co/mem_j  : share of the SMALLER backbone (partial-dup of a small paralog can be high)
  aln_cov_small        : fraction of the SHORTER cDNA the alignment covers
  contig               : longest contiguous run of co-threaded ref junctions / co  (co-linearity)
  seq_id               : cDNA identity from rep_ava (context)

Edge sets: BRIDGE = internal both-spliced edges of the large post-retrocopy components (the over-merge
to break); REAL = panel families + Compara checkable multi-exon pairs + GGT1~GGTLC2 (must stay).

Run: /home/juanfra/miniforge3/bin/python bench/family_def_cothread_genomewide.py [--min-comp N] [--all]
Out: bench/family_def_cothread_genomewide.{tsv,json}
"""
import collections
import json
import os
import re
import subprocess
import sys
import tempfile
from multiprocessing import Pool

import pysam

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from family_def_read_filters import dna_homology
from family_def_retrocopy_filter import build_family_graph, load_intron_index

GENOME = "/home/juanfra/winloci_scratch/GGO.fasta"
EXON_INDEX = "/home/juanfra/winloci_scratch/gene_exon_index.json"
COMPARA = os.path.join(HERE, "compara_paralog_relation.json")
JUNC_TOL = 12
CIG = re.compile(r"(\d+)([MIDNSHP=X])")

MIN_COMP = 20            # default: only de-bridge components with >= this many members
EXON_IDX = None          # globals for forked workers
CACHE = None             # gene -> (cdna_seq, [junc offsets], strand)


# ---------------- cDNA construction ----------------
def revcomp(s):
    return s.translate(str.maketrans("ACGTNacgtn", "TGCANtgcan"))[::-1]


def build_cdna(name, fa):
    v = EXON_IDX.get(name)
    if v is None or len(v["exons"]) == 0:
        return None
    blocks = [fa.fetch(v["chrom"], a - 1, b).upper() for a, b in v["exons"]]
    if v["strand"] == "-":
        blocks = [revcomp(b) for b in reversed(blocks)]
    cdna = "".join(blocks)
    junc, cum = [], 0
    for b in blocks[:-1]:
        cum += len(b)
        junc.append(cum)
    return cdna, junc, v["strand"]


# ---------------- junction projection through one primary alignment ----------------
def project(ref_junc, mem_seq, mem_junc, sam_fields):
    rstart = int(sam_fields[3]) - 1
    cigar = sam_fields[5]
    is_rev = bool(int(sam_fields[1]) & 0x10)
    mlen = len(mem_seq)
    mj = [mlen - j for j in reversed(mem_junc)] if is_rev else list(mem_junc)
    qpos = rpos = covered = 0
    q2r = {}
    rpos = rstart
    for n, op in CIG.findall(cigar):
        n = int(n)
        if op in "=XM":
            for k in range(n):
                q2r[qpos + k] = rpos + k
            qpos += n; rpos += n; covered += n
        elif op == "I":
            qpos += n
        elif op == "D":
            rpos += n
        elif op in "SH":
            qpos += n
    co_positions = []
    ref_set = sorted(ref_junc)
    for j in mj:
        rj = None
        for d in range(0, JUNC_TOL + 1):
            if j - d in q2r:
                rj = q2r[j - d]; break
            if j + d in q2r:
                rj = q2r[j + d]; break
        if rj is None:
            continue
        for ri, rk in enumerate(ref_set):
            if abs(rj - rk) <= JUNC_TOL:
                co_positions.append(ri); break
    co_positions = sorted(set(co_positions))
    # longest contiguous run of co-threaded ref-junction indices
    longest = run = 0
    for i, ri in enumerate(co_positions):
        if i and ri == co_positions[i - 1] + 1:
            run += 1
        else:
            run = 1
        longest = max(longest, run)
    return co_positions, covered / max(mlen, 1), longest


# ---------------- per-reference worker: batch-map all partners onto one reference ----------------
def work_reference(args):
    ref, partners, seq_id = args      # ref gene; partners list; seq_id dict {(a,b):id}
    rc = CACHE.get(ref)
    if rc is None:
        return []
    ref_seq, ref_junc, _ = rc
    tmp = tempfile.mkdtemp(prefix=f"ct_{os.getpid()}_")
    refp = os.path.join(tmp, "r.fa")
    memp = os.path.join(tmp, "m.fa")
    with open(refp, "w") as fh:
        fh.write(f">ref\n{ref_seq}\n")
    name_of = {}
    with open(memp, "w") as fh:
        for i, p in enumerate(partners):
            mc = CACHE.get(p)
            if mc is None:
                continue
            qn = f"m{i}"
            name_of[qn] = p
            fh.write(f">{qn}\n{mc[0]}\n")
    if not name_of:
        return []
    res = subprocess.run(["minimap2", "-a", "--eqx", "-x", "asm20", "--secondary=no", refp, memp],
                         capture_output=True, text=True)
    prim = {}
    for line in res.stdout.splitlines():
        if line.startswith("@"):
            continue
        f = line.split("\t")
        flag = int(f[1])
        if flag & 0x4 or flag & 0x100 or flag & 0x800:
            continue
        prim.setdefault(f[0], f)   # first primary per query
    out = []
    nref = len(ref_junc)
    for qn, mem in name_of.items():
        ms, mj, _ = CACHE[mem]
        edge = tuple(sorted((ref, mem)))
        rec = dict(ref=ref, mem=mem, ref_j=nref, mem_j=len(mj),
                   seq_id=seq_id.get(edge, 0.0))
        f = prim.get(qn)
        if f is None or nref == 0:
            rec.update(co=0, frac_ref=0.0, frac_mem=0.0, aln_cov_small=0.0, contig=0)
        else:
            co, cov, longest = project(ref_junc, ms, mj, f)
            small_len = min(len(ref_seq), len(ms))
            # aln_cov of the smaller gene: covered query bases / smaller cDNA length
            rec.update(co=len(co),
                       frac_ref=round(len(co) / nref, 3),
                       frac_mem=round(len(co) / max(len(mj), 1), 3),
                       aln_cov_small=round((cov * len(ms)) / max(small_len, 1), 3),
                       contig=longest)
        out.append(rec)
    for fn in (refp, memp):
        try: os.remove(fn)
        except OSError: pass
    try: os.rmdir(tmp)
    except OSError: pass
    return out


def compara_real_pairs(idx):
    """Held-out / external real paralog pairs: Compara checkable, both multi-exon, both in exon idx."""
    pairs = set()
    try:
        cj = json.load(open(COMPARA))
    except (OSError, ValueError):
        return pairs
    for fam, info in cj.get("families", {}).items():
        if not info.get("checkable_pair"):
            continue
        gs = [g for g in info.get("mapped_genes", []) if g in EXON_IDX
              and (idx.get(g, {}).get("n_intron", 0) or 0) >= 1]
        for i in range(len(gs)):
            for j in range(i + 1, len(gs)):
                pairs.add(tuple(sorted((gs[i], gs[j]))))
    return pairs


def main():
    global EXON_IDX, CACHE, MIN_COMP
    do_all = "--all" in sys.argv
    if "--min-comp" in sys.argv:
        MIN_COMP = int(sys.argv[sys.argv.index("--min-comp") + 1])
    Hd, _ = dna_homology()
    idx = load_intron_index()
    EXON_IDX = json.load(open(EXON_INDEX))
    seq_id = {tuple(sorted((a, b))): r.get("id", 0.0) for (a, b), r in Hd.items()}

    # post-retrocopy graph -> both-spliced edges, partitioned by component size
    import networkx as nx
    G, prov, part = build_family_graph(Hd, idx, filter_retrocopy=True)
    both_spliced = set(map(lambda e: tuple(sorted(e)), part["both_spliced_kept"]))
    comp_of = {}
    big_comps = []
    for ci, c in enumerate(nx.connected_components(G)):
        for g in c:
            comp_of[g] = ci
        if len(c) >= MIN_COMP:
            big_comps.append((ci, len(c)))
    big_ids = {ci for ci, _ in big_comps}

    # BRIDGE candidate edges: both-spliced, inside a big component
    bridge_edges = {e for e in both_spliced if comp_of.get(e[0]) in big_ids and comp_of.get(e[0]) == comp_of.get(e[1])}
    if do_all:
        bridge_edges = both_spliced

    # REAL control edges (must survive): panel + GGT + Compara checkable multi-exon
    panel_real = [("RABL2A", "RABL2B"), ("APOBEC3D", "APOBEC3F"),
                  ("RFPL1", "RFPL2"), ("RFPL2", "RFPL3"), ("RFPL1", "RFPL3"),
                  ("GGT1", "GGTLC2")]
    real_edges = {tuple(sorted(e)) for e in panel_real} | compara_real_pairs(idx)
    # keep only real edges that are actually homology edges we can align (both in exon idx)
    real_edges = {e for e in real_edges if e[0] in EXON_IDX and e[1] in EXON_IDX}

    label = {}
    for e in bridge_edges:
        label[e] = "bridge"
    for e in real_edges:
        label[e] = "real"          # real overrides if overlap
    all_edges = list(label)
    print(f"[setup] big components (>= {MIN_COMP}): {sorted([s for _,s in big_comps], reverse=True)[:8]}")
    print(f"[setup] bridge-candidate edges: {len(bridge_edges)} ; real-control edges: {len(real_edges)} "
          f"; total to align: {len(all_edges)}")

    # build cDNA cache for all genes involved
    genes = {g for e in all_edges for g in e}
    fa = pysam.FastaFile(GENOME)
    cache = {}
    miss = 0
    for g in genes:
        c = build_cdna(g, fa)
        if c is None:
            miss += 1
        else:
            cache[g] = c
    fa.close()
    CACHE = cache
    print(f"[setup] cDNA built for {len(cache)}/{len(genes)} genes ({miss} missing exon coords)")

    # group edges by reference = endpoint with MORE junctions (richer backbone)
    def njunc(g):
        c = cache.get(g)
        return len(c[1]) if c else -1
    by_ref = collections.defaultdict(list)
    for a, b in all_edges:
        if a not in cache or b not in cache:
            continue
        ref, mem = (a, b) if njunc(a) >= njunc(b) else (b, a)
        by_ref[ref].append(mem)
    jobs = [(ref, mems, seq_id) for ref, mems in by_ref.items()]
    print(f"[run] {len(jobs)} reference groups, aligning {sum(len(m) for _,m in by_ref.items())} edges ...")

    records = []
    with Pool(5) as pool:
        for batch in pool.imap_unordered(work_reference, jobs, chunksize=8):
            records.extend(batch)

    # attach label
    for r in records:
        r["label"] = label.get(tuple(sorted((r["ref"], r["mem"]))), "bridge")

    # write TSV
    cols = ["label", "ref", "mem", "ref_j", "mem_j", "co", "frac_ref", "frac_mem",
            "aln_cov_small", "contig", "seq_id"]
    with open(os.path.join(HERE, "family_def_cothread_genomewide.tsv"), "w") as out:
        out.write("\t".join(cols) + "\n")
        for r in records:
            out.write("\t".join(str(r.get(c, "")) for c in cols) + "\n")

    # summary: distributions per label + where the real controls land
    def summ(recs, key):
        vals = sorted(r[key] for r in recs)
        if not vals:
            return "n=0"
        n = len(vals)
        q = lambda p: vals[min(n - 1, int(p * n))]
        return f"n={n} med={q(0.5):.2f} p10={q(0.1):.2f} p90={q(0.9):.2f} min={vals[0]:.2f} max={vals[-1]:.2f}"
    reals = [r for r in records if r["label"] == "real"]
    bridges = [r for r in records if r["label"] == "bridge"]
    print("\n=== position-feature distributions (REAL controls vs BRIDGE candidates) ===")
    for key in ("frac_ref", "frac_mem", "aln_cov_small", "contig"):
        print(f"  {key:14} REAL  : {summ(reals, key)}")
        print(f"  {key:14} BRIDGE: {summ(bridges, key)}")
    print("\n=== the named real controls (must stay high on the right feature) ===")
    named = {("GGT1", "GGTLC2"), ("RABL2A", "RABL2B"), ("APOBEC3D", "APOBEC3F"),
             ("RFPL2", "RFPL3"), ("RFPL1", "RFPL2")}
    for r in records:
        if tuple(sorted((r["ref"], r["mem"]))) in named:
            print(f"  {r['ref']:10}~{r['mem']:12} ref_j={r['ref_j']:2} mem_j={r['mem_j']:2} "
                  f"frac_ref={r['frac_ref']:.2f} frac_mem={r['frac_mem']:.2f} "
                  f"alnCovSmall={r['aln_cov_small']:.2f} contig={r['contig']} id={r['seq_id']:.2f}")

    json.dump(dict(min_comp=MIN_COMP, do_all=do_all,
                   n_bridge=len(bridges), n_real=len(reals),
                   prov=prov, records=records[:5000]),
              open(os.path.join(HERE, "family_def_cothread_genomewide.json"), "w"), indent=2)
    print(f"\n[+] wrote family_def_cothread_genomewide.{{tsv,json}} ({len(records)} edges)")


if __name__ == "__main__":
    main()
