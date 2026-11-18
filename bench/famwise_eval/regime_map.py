"""PREREG r1795 (docs/PREREG_famwise_origin_2026-11-18.md), row c: where, in differences between copies and in their
expression, perfect provenance changes an assembly several-fold (bench/famsim, planted truth).

  regime_map.py --out DIR --bin RELEASE_DIR [--seeds 1,2,3,4,5] [--only CELL,...] [--jobs N] [--threads T] [--promoted-mapq Q]

Cells: k{0,1,2,4,8,20}_{1to1,10to1,silent}_{full,frag} (two copies A and B of one 6-exon template, B with k exonic
SNPs; reads A / B 200 / 200, 200 / 20, 200 / 0; full-length reads or 3' fragments) and array_{full,frag} (8 identical
copies, 25 reads each). Each copy has the full chain and one exon-skipping isoform (weight 0.3).
Per cell and seed, one alignment and three assemblies: P = primary alignments only, D = the default, O = each read kept
only at its planted copy (the record there becomes the primary); then the stage `famwise` on D.
DIR/map.tsv   one row per cell, seed and input: planted (isoform, copy) pairs with reads, how many are emitted with their
              exact intron chain, multi-exon chains emitted on a copy without reads, share of the cell's reads whose two best
              records tie exactly on ms / AS.
DIR/cells.tsv per cell the medians over the seeds and `several_fold` (median O / median P >= 2, or P none and O some)."""
import argparse, collections, csv, json, multiprocessing, os, shutil, statistics, subprocess, sys
HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.dirname(os.path.dirname(HERE))
DRIVER = os.path.join(REPO, "tools", "rustle_pipeline.sh")

EXONS, INTRONS = [300, 200, 400, 250, 150, 700], [300] * 5       # 2,000 bp of exons: k differences = rate k / 2000
KS, EXPR, READS = (0, 1, 2, 4, 8, 20), {"1to1": (200, 200), "10to1": (200, 20), "silent": (200, 0)}, ("full", "frag")
ISO = [{"skip": [3], "weight": 0.3}]
CELLS = [f"k{k}_{e}_{r}" for k in KS for e in EXPR for r in READS] + [f"array_{r}" for r in READS]

def spec_for(cell, seed, template):
    """the famsim spec of a cell; template = path of the saved template, or None for the synthetic one"""
    *what, reads = cell.split("_")
    if what == ["array"]:
        copies = [{"id": f"A{i}", "expression": 25, "isoforms": ISO, "ops": []} for i in range(1, 9)]
    else:
        k, (ea, eb) = int(what[0][1:]), EXPR[what[1]]
        copies = [{"id": "A", "expression": ea, "isoforms": ISO, "ops": []},
                  {"id": "B", "expression": eb, "isoforms": ISO, "ops": [{"op": "snp", "rate": k / sum(EXONS), "region": "exons"}] if k else []}]
    return {"name": cell, "seed": seed, "background": {"source": "random", "length": 200000},
            "template": {"source": "file", "path": template} if template else {"source": "synthetic", "exons": EXONS, "introns": INTRONS},
            "copies": copies, "reads": dict({"per_copy": 200, "jitter": 30}, **({"trunc5_frac": 1.0, "trunc5_max": 0.7} if reads == "frag" else {}))}

def several_fold(pairs):
    """pairs = [(recovered P, recovered O)] per seed -> (median P, median O, yes / no)"""
    p, o = statistics.median(x for x, _ in pairs), statistics.median(y for _, y in pairs)
    return p, o, "yes" if (o > 0 and p == 0) or (p > 0 and o / p >= 2) else "no"

def sh(cmd, log):
    with open(log, "a") as fh:
        fh.write(f"$ {cmd}\n"); fh.flush()
        if subprocess.run(cmd, shell=True, stdout=fh, stderr=fh, cwd=REPO).returncode:
            raise RuntimeError(f"FAILED: {cmd} (see {log})")

def origin_bam(d, mapq):
    """d/O.bam: every read kept only at its planted copy (read name copy|isoform|i) -> counts"""
    import pysam
    span = {r["copy_id"]: (r["contig"], int(r["start"]) - 1, int(r["end"])) for r in csv.DictReader(open(f"{d}/copies.tsv"), delimiter="\t")}
    by, n = collections.defaultdict(list), collections.Counter()
    with pysam.AlignmentFile(f"{d}/reads.bam") as fh:
        for r in fh:
            if not (r.is_unmapped or r.is_supplementary):
                by[r.query_name].append(r)
        with pysam.AlignmentFile(f"{d}/O.unsorted.bam", "wb", template=fh) as w:
            for q, rs in by.items():
                n["reads"] += 1
                for tag in ("ms", "AS"):     # the read's two best records tie exactly
                    v = sorted((r.get_tag(tag) for r in rs), reverse=True); n[f"tied_{tag}"] += len(v) > 1 and v[0] == v[1]
                c, s, e = span[q.split("|")[0]]
                home = [r for r in rs if r.reference_name == c and r.reference_start < e and s < r.reference_end]
                if not home:
                    n["no_home"] += 1; home = [r for r in rs if not r.is_secondary][:1] or rs[:1]
                elif all(r.is_secondary for r in home):
                    n["moved"] += 1; r = max(home, key=lambda r: r.get_tag("ms")); r.flag &= ~0x100
                    if mapq is not None:
                        r.mapping_quality = mapq
                n["dropped_records"] += len(rs) - len(home)
                for r in home:
                    w.write(r)
    pysam.sort("-o", f"{d}/O.bam", f"{d}/O.unsorted.bam"); pysam.index(f"{d}/O.bam"); os.remove(f"{d}/O.unsorted.bam")
    return n

def emitted(gtf):
    """[(contig, start, end, frozenset of introns)] of the multi-exon transcripts of a GTF (0-based half-open)"""
    sys.path.insert(0, os.path.join(REPO, "bench"))
    from famsim import evaluate
    out = []
    for L in evaluate.gtf_loci(gtf).values():
        for ex in L["transcripts"].values():
            ex = sorted(ex)
            if len(ex) > 1:
                out.append((L["contig"], ex[0][0], ex[-1][1], frozenset((ex[i][1], ex[i + 1][0]) for i in range(len(ex) - 1))))
    return out

def score(d, gtf):
    """-> (planted pairs with reads, recovered, chains on copies without reads)"""
    sys.path.insert(0, os.path.join(REPO, "bench"))
    from famsim import chromosome
    em, pairs, rec, silent = emitted(gtf), 0, 0, set()
    for p in chromosome.load_planted(d)[0]:
        if p.kind != "copy" or not p.in_reference:
            continue
        if p.expression == 0:
            silent |= {(p.id, ch) for c, s, e, ch in em if c == p.contig and s < p.end and p.pos < e}
            continue
        for _, skip, weight in p.isoform_list():
            if weight > 0:
                pairs += 1; rec += any(c == p.contig and ch == frozenset(p.junctions(skip)) for c, _, _, ch in em)
    return pairs, rec, len(silent)

def unit(job):
    """one cell and seed -> (rows of map.tsv, layer line)"""
    cell, seed, out, bindir, threads, mapq = job
    d = f"{out}/{cell}.s{seed}"
    shutil.rmtree(d, ignore_errors=True); os.makedirs(d)
    log = f"{d}.log"; open(log, "w").close()
    *what, reads = cell.split("_")
    head = [cell, "array" if what == ["array"] else what[0][1:], "8x25" if what == ["array"] else what[1], reads, seed]
    try:
        json.dump(spec_for(cell, seed, f"{out}/template.json"), open(f"{d}.json", "w"))
        sh(f"{sys.executable} bench/famsim make {d}.json --out {d}", log)
        subprocess.run(f"{sys.executable} bench/famsim verify {d}", shell=True, stdout=open(log, "a"), stderr=subprocess.STDOUT, cwd=REPO)
        fails = sorted({l.split("\t")[0] for l in open(f"{d}/verify.tsv") if "\tFAIL\t" in l})
        if set(fails) - {"pair_identity_genomic"}:
            raise RuntimeError(f"verify FAIL: {','.join(fails)}")
        sh(f"{sys.executable} bench/famsim align {d} --threads {threads}", log)
        n = origin_bam(d, mapq)
        drv = lambda stage, bam, x, flags="": f"bash {DRIVER} {stage} --bam {d}/{bam} --fasta {d}/genome.ref.fa --out {d}/{x}/run --bin {bindir} --threads {threads} {flags}"
        rows = []
        for x, bam, flags in (("P", "reads.bam", "--no-cache --seed-pool primary"), ("D", "reads.bam", "--no-cache"), ("O", "O.bam", "--no-cache --seed-pool all")):
            os.makedirs(f"{d}/{x}")
            sh(drv("assemble", bam, x, flags), log)
            note = ";".join(filter(None, [f"verify_fail={','.join(fails)}" if fails else "", f"moved={n['moved']},no_home={n['no_home']}" if x == "O" else ""]))
            rows.append([*head, x, *score(d, f"{d}/{x}/run.gtf"), *(f"{n[t] / max(1, n['reads']):.3f}" for t in ("tied_ms", "tied_AS")), note])
        try:
            sh(drv("famwise", "reads.bam", "D"), log)
            layer = open(f"{d}/D/run.famwise.log").read().strip().replace("\n", " | ")
        except (RuntimeError, OSError) as e:
            layer = f"famwise did not run: {e}"
        return rows, [*head, layer]
    except Exception as e:      # a failed cell is reported, the map goes on
        return [[*head, "FAILED", "", "", "", "", "", str(e).replace("\n", " ")[:300]]], [*head, "FAILED"]

def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--out", required=True); ap.add_argument("--bin", required=True); ap.add_argument("--seeds", default="1,2,3,4,5")
    ap.add_argument("--only"); ap.add_argument("--jobs", type=int, default=1); ap.add_argument("--threads", type=int, default=2)
    ap.add_argument("--promoted-mapq", type=int, help="MAPQ written on the record that becomes a read's primary in O (default: unchanged)")
    a = ap.parse_args()
    out = os.path.abspath(a.out); os.makedirs(out, exist_ok=True)
    cells = a.only.split(",") if a.only else CELLS
    assert not set(cells) - set(CELLS), f"unknown cell (cells: {', '.join(CELLS)})"
    if not os.path.exists(f"{out}/template.json"):      # one template for every cell and seed: seeds change reads and SNP positions only
        t = f"{out}/_template"; json.dump(spec_for("k0_1to1_full", 1, None), open(f"{t}.json", "w"))
        sh(f"{sys.executable} bench/famsim make {t}.json --out {t}", f"{t}.log"); shutil.copy(f"{t}/template.json", f"{out}/template.json")
    jobs = [(c, int(s), out, a.bin, a.threads, a.promoted_mapq) for c in cells for s in a.seeds.split(",")]
    with multiprocessing.Pool(a.jobs) as pool:
        res = pool.map(unit, jobs, chunksize=1)
    cols = ("cell", "differences", "expression", "reads", "seed", "input", "planted_pairs", "recovered", "silent_chains", "tied_ms", "tied_AS", "note")
    with open(f"{out}/map.tsv", "w") as fh:
        fh.write("\t".join(cols) + "\n"); fh.writelines("\t".join(map(str, r)) + "\n" for rows, _ in res for r in rows)
    with open(f"{out}/layer.tsv", "w") as fh:
        fh.write("cell\tdifferences\texpression\treads\tseed\tfamwise\n"); fh.writelines("\t".join(map(str, l)) + "\n" for _, l in res)
    by = collections.defaultdict(dict)
    for rows, _ in res:
        for r in rows:
            if r[5] != "FAILED":
                by[tuple(r[:4])].setdefault(r[4], {})[r[5]] = r[7]
    with open(f"{out}/cells.tsv", "w") as fh:
        fh.write("cell\tdifferences\texpression\treads\tseeds\tP_median\tD_median\tO_median\tseveral_fold\n")
        for k, seeds in by.items():
            p, o, sf = several_fold([(v["P"], v["O"]) for v in seeds.values()])
            fh.write("\t".join(map(str, (*k, len(seeds), p, statistics.median(v["D"] for v in seeds.values()), o, sf))) + "\n")
    print(open(f"{out}/cells.tsv").read(), end="")
    failed = [r for rows, _ in res for r in rows if r[5] == "FAILED"]
    if failed:
        print(f"{len(failed)} cell x seed FAILED (map.tsv, column note)", file=sys.stderr)

if __name__ == "__main__":
    main()
