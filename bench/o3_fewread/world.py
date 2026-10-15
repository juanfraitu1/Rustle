#!/usr/bin/env python3
"""Few-read false-flag worlds for `candidate_copies` (docs/PREREG_o3_fewread_flags_2026-10-14.md, Phase 1).

    world.py build OUT --seed S          # synthetic world (bench/famsim via synth_world.build_world) + the stage's inputs
    world.py run OUT --variant E0|E3     # map the variant's reads (sim.MM2), run candidate_copies -> OUT/E0/P.cand.*

Family types (per seed, 48 families): NULL3, NULLclose, NULLlow (nothing missing), ALLELE@d (a second allele absent from the
reference), ERASED (one copy absent, at a range of expression levels). The copies table holds every IN-REFERENCE copy with its
truth exons (`mcl_families`' P.fam.copies.tsv / .fa layout), so the stage sees exactly the reference's copies.
"""
import argparse
import json
import os
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
BENCH = os.path.dirname(HERE)
sys.path.insert(0, BENCH)
sys.path.insert(0, os.path.join(BENCH, "unmapped_rescue"))

LOW_C = (4, 6, 8, 12, 16, 24, 32, 48)           # NULLlow: copy C's reads, one level per family
ERASED_E = (4, 4, 6, 6, 8, 8, 12, 12, 16, 16, 40, 40, 80, 80)
ALLELE_D = (0.0025, 0.005, 0.01, 0.02, 0.04)
COPIES_HEADER = ("family_id\tcopy_idx\ttid\tchrom\tstart\tend\tn_exon\tstrand\tn_reads\texons\tmax_family_identity\tsource\tgene_id\t"
                 "core_hull\tsd_depth\tcore_bp\trep_frac\tmember_status\tlocus_start\tlocus_end")


def specs(seed):
    from famsim.model import seeded
    snp = lambda r: [{"op": "snp", "rate": r}] if r else []
    plan = ([("NULL3", i) for i in range(8)] + [("NULLclose", i) for i in range(8)] + [("NULLlow", i) for i in range(8)] +
            [(f"ALLELE@{d}", d) for d in ALLELE_D for _ in range(2)] + [("ERASED", i) for i in range(14)])
    out = []
    for k, (sc, x) in enumerate(plan):
        name = f"f{k:02d}"
        rng = seeded(seed, "family", name)
        n = rng.randint(5, 10)
        exons = [rng.randint(100, 250)] + [rng.randint(90, 350) for _ in range(n - 2)] + [rng.randint(250, 800)]
        introns = [rng.randint(300, 2500) for _ in range(n - 1)]
        if sc == "NULL3":
            cp = [("A", 0, True, None), ("B", 0.02, True, None), ("C", 0.03, True, None)]
        elif sc == "NULLclose":
            cp = [("A", 0, True, None), ("B", 0.01, True, None), ("C", 0.015, True, None)]
        elif sc == "NULLlow":
            cp = [("A", 0, True, None), ("B", 0.02, True, None), ("C", 0.03, True, LOW_C[x])]
        elif sc == "ERASED":
            cp = [("A", 0.02, True, None), ("B", 0.03, True, None), ("E", 0, False, ERASED_E[x])]
        else:
            cp = [("A", 0, True, None), ("B", 0.03, True, None), ("H", x, False, None)]
        copies, kinds = [], {}
        for tag, r, inref, expr in cp:
            c = {"id": name + tag, "ops": snp(r)}
            if not inref:
                c["in_reference"] = False
            if expr is not None:
                c["expression"] = expr
            copies.append(c)
            kinds[name + tag] = "present" if inref else ("copy" if tag == "E" else "allele")
        out.append({"name": name, "seed": seed * 1000 + k, "scenario": sc, "kinds": kinds,
                    "divergence_class": x if sc.startswith("ALLELE") else 0.02,
                    "background": {"source": "random", "length": 60000}, "template": {"source": "synthetic", "exons": exons, "introns": introns},
                    "copies": copies, "reads": {"per_copy": 80, "err": 0.001, "indel": 0.0003, "jitter": 30}})
    return out


def read_fa(path):
    from synth_world import read_fa as rf
    return rf(path)


def truth_exons(gtf, gene):
    """iso0 exons of `gene` in a famsim truth.gtf, 0-based half-open, and the transcript's 1-based start"""
    ex, tstart = [], None
    for ln in open(gtf):
        r = ln.rstrip("\n").split("\t")
        if len(r) < 9 or f'gene_id "{gene}"' not in r[8] or f'transcript_id "{gene}.iso0"' not in r[8]:
            continue
        if r[2] == "transcript":
            tstart = int(r[3])
        elif r[2] == "exon":
            ex.append((int(r[3]) - 1, int(r[4])))
    return sorted(ex), tstart


def write_copies(out):
    ref = read_fa(f"{out}/genome.ref.fa")
    truth = read_fa(f"{out}/genome.truth.fa")
    rows, fa = [], []
    idx = {}
    for ln in list(open(f"{out}/copies.tsv"))[1:]:
        copy, fam, role, D, ctg, pos0, end, strand, n_ex, inref = ln.rstrip("\n").split("\t")
        if inref != "True":
            continue
        pos0, end = int(pos0), int(end)
        ex, tstart = truth_exons(f"{out}/fam/{fam}/truth.gtf", copy)
        off = pos0 - (tstart - 1)
        ex = [(s + off, e + off) for s, e in ex]
        tctg = ctg                                    # world contigs are named <family>_<famsim contig> in both FASTAs
        assert ref[ctg][pos0:end] == truth[tctg][tstart - 1:tstart - 1 + (end - pos0)], f"{copy}: reference span differs from truth"
        assert ex[0][0] == pos0 and ex[-1][1] == end, f"{copy}: exon hull {ex[0][0]}-{ex[-1][1]} vs span {pos0}-{end}"
        seq = "".join(ref[ctg][s:e] for s, e in ex)
        i = idx.get(fam, 0)
        idx[fam] = i + 1
        rows.append([fam, str(i), f"SIM_{copy}", ctg, str(pos0), str(end), str(len(ex)), strand, "0", ",".join(f"{s}-{e}" for s, e in ex),
                     "1.000", "famsim", copy, "NA", "1", str(len(seq)), "0.000", "kept", str(pos0), str(end)])
        fa.append(f">{fam}|{i}|{ctg}:{pos0}-{end}|{strand}|nexon={len(ex)}\n{seq}\n")
    with open(f"{out}/P.fam.copies.tsv", "w") as f:
        f.write(COPIES_HEADER + "\n" + "".join("\t".join(r) + "\n" for r in rows))
    open(f"{out}/P.fam.copies.fa", "w").write("".join(fa))
    return len(rows)


def cmd_build(a):
    from synth_world import build_world
    build_world(a.out, a.seed, log=lambda *x, **k: None, specs=specs(a.seed))
    n = write_copies(a.out)
    sh(f"samtools faidx {a.out}/genome.ref.fa && minimap2 -x splice:hq -d {a.out}/genome.ref.splice.mmi {a.out}/genome.ref.fa 2>/dev/null")
    print(f"[o3_fewread] world {a.out} seed {a.seed}: {n} in-reference copies", flush=True)


def sh(cmd):
    subprocess.run(cmd, shell=True, check=True, executable="/bin/bash")


def cmd_run(a):
    import sim
    w, v = a.out, a.variant
    os.makedirs(f"{w}/{v}", exist_ok=True)
    bam = f"{w}/{v}/reads.bam"
    if not os.path.exists(bam + ".bai"):
        sh(f"{sim.MM2} -t {a.threads} {w}/genome.ref.fa {w}/reads.{v}.fq 2>/dev/null | samtools sort -@2 -o {bam} - && samtools index {bam}")
    binary = a.bin or os.path.join(os.environ.get("CARGO_TARGET_DIR", os.path.expanduser("~/rustle_target")), "release", "candidate_copies")
    sh(f"{binary} --bam {bam} --fasta {w}/genome.ref.fa --copies {w}/P.fam.copies.tsv --copies-fa {w}/P.fam.copies.fa "
       f"--index {w}/genome.ref.splice.mmi --threads {a.threads} --out {w}/{v}/P.cand 2> {w}/{v}/P.cand.log")
    print(open(f"{w}/{v}/P.cand.log").read().splitlines()[-1], flush=True)


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd", required=True)
    p = sub.add_parser("build"); p.add_argument("out"); p.add_argument("--seed", type=int, required=True); p.set_defaults(func=cmd_build)
    p = sub.add_parser("run"); p.add_argument("out"); p.add_argument("--variant", choices=("E0", "E3"), required=True)
    p.add_argument("--threads", type=int, default=4); p.add_argument("--bin"); p.set_defaults(func=cmd_run)
    a = ap.parse_args()
    a.func(a)
