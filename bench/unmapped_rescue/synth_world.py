#!/usr/bin/env python3
"""Supervised synthetic world (docs/PREREG_unmapped_rescue_2026-10-08.md Amendment 10): random genome, synthetic gene families, one copy per family
erased from the reference, reads simulated from every copy. Built with bench/famsim. Exact truth for every read and every copy.

    synth_world.py build <out_dir> [--seed 20261008]    # writes <out_dir>/{genome.ref.fa, genome.truth.fa, copies.tsv, transcripts.fa, targets.fa, reads.E0|E1|E2.fq, reads.truth.tsv}
"""
import argparse
import json
import os
import random
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
BENCH = os.path.dirname(HERE)
sys.path.insert(0, BENCH)

SEED = 20261008
CLASSES = (0.005, 0.01, 0.02, 0.04, 0.08)


def _seeded(seed, *parts):
    from famsim.model import seeded
    return seeded(seed, *parts)


def family_specs(seed=SEED, n_per_class=6, classes=CLASSES, per_copy=80, lead_g=False):
    """famsim specs, one per family: A (template, present), B (present, SNP rate 1.5 D), E (erased, SNP rate D); D = the divergence class."""
    out, k = [], 0
    for D in classes:
        for _ in range(n_per_class):
            name = f"f{k:02d}"
            rng = _seeded(seed, "family", name)
            n = rng.randint(5, 10)
            exons = [rng.randint(100, 250)] + [rng.randint(90, 350) for _ in range(n - 2)] + [rng.randint(250, 800)]
            introns = [rng.randint(300, 2500) for _ in range(n - 1)]
            out.append({"name": name, "seed": seed * 1000 + k, "divergence_class": D,
                        "background": {"source": "random", "length": 60000},
                        "template": {"source": "synthetic", "exons": exons, "introns": introns},
                        "copies": [{"id": name + "A"},
                                   {"id": name + "B", "ops": [{"op": "snp", "rate": 1.5 * D}]},
                                   {"id": name + "E", "in_reference": False, "ops": [{"op": "snp", "rate": D}]}],
                        "reads": {"per_copy": per_copy, "err": 0.001, "indel": 0.0003, "jitter": 30}})
            if lead_g:
                out[-1]["lead_g"] = k % 4          # the template's first exon begins with this many G (GC-rich 5' ends, Amendment 16)
            k += 1
    return out


def partition_specs(seed=SEED, n_each=12, per_copy=100):
    """Amendment 18: n_each families of each scenario. CTL: one erased copy, one isoform. ISO: the erased copy also makes an isoform skipping internal exon 3
    (reads 2 : 1). SIB: two erased copies, E1 the template and E2 0.5% SNPs from it (inside the edge cutoff). Survivors A (2% SNPs) and B (3%)."""
    out, k = [], 0
    for sc in ("CTL", "ISO", "SIB"):
        for _ in range(n_each):
            name = f"f{k:02d}"
            rng = _seeded(seed, "family", name)
            n = rng.randint(5, 10)
            exons = [rng.randint(100, 250)] + [rng.randint(90, 350) for _ in range(n - 2)] + [rng.randint(250, 800)]
            introns = [rng.randint(300, 2500) for _ in range(n - 1)]
            copies = [{"id": name + "A", "ops": [{"op": "snp", "rate": 0.02}]}, {"id": name + "B", "ops": [{"op": "snp", "rate": 0.03}]}]
            if sc == "SIB":
                copies += [{"id": name + "E1", "in_reference": False}, {"id": name + "E2", "in_reference": False, "ops": [{"op": "snp", "rate": 0.005}]}]
            else:
                e = {"id": name + "E", "in_reference": False}
                if sc == "ISO":
                    e["isoforms"] = [{"skip": [3], "weight": 0.5}]
                copies.append(e)
            out.append({"name": name, "seed": seed * 1000 + k, "scenario": sc, "divergence_class": 0.02, "background": {"source": "random", "length": 60000},
                        "template": {"source": "synthetic", "exons": exons, "introns": introns}, "copies": copies,
                        "reads": {"per_copy": per_copy, "err": 0.001, "indel": 0.0003, "jitter": 30}})
            k += 1
    return out


def add_polya(fq_lines, seed):
    """a 3' polyA tail of 20 to 30 A on every read of a 4-line FASTQ stream (the E1 read-end variant)"""
    rng = random.Random(seed)
    out = []
    for i in range(0, len(fq_lines), 4):
        n = rng.randint(20, 30)
        out += [fq_lines[i], fq_lines[i + 1] + "A" * n, fq_lines[i + 2], fq_lines[i + 3] + "I" * n]
    return out


G5_P = 0.95
G5_LENGTHS = ((1, 0.22), (2, 0.43), (3, 0.24), (4, 0.07), (5, 0.03), (6, 0.01))   # measured on the real dev-A reads (Amendment 13)


def add_g5(fq_lines, seed, p=G5_P, lengths=G5_LENGTHS):
    """an untemplated 5' run of G on a fraction p of the reads of a 4-line FASTQ stream, run length from `lengths` (the E3 read-end variant)"""
    rng = random.Random(seed)
    out = []
    for i in range(0, len(fq_lines), 4):
        n = 0
        if rng.random() < p:
            u, acc = rng.random() * sum(w for _, w in lengths), 0.0
            for k, w in lengths:
                acc += w
                n = k
                if u <= acc:
                    break
        out += [fq_lines[i], "G" * n + fq_lines[i + 1], fq_lines[i + 2], "I" * n + fq_lines[i + 3]]
    return out


def read_fa(path):
    out, name = {}, None
    for ln in open(path):
        if ln[0] == ">":
            name = ln[1:].strip().split()[0]
            out[name] = []
        else:
            out[name].append(ln.strip())
    return {k: "".join(v) for k, v in out.items()}


def write_fa(path, seqs):
    with open(path, "w") as o:
        for k, v in seqs.items():
            o.write(f">{k}\n{v}\n")


def build_world(out, seed=SEED, log=print, n_per_class=6, classes=CLASSES, lead_g=False, specs=None):
    from famsim import chromosome, reads as FR
    os.makedirs(f"{out}/fam", exist_ok=True)
    ref, truth, copies, trans, targets = {}, {}, [], {}, {}
    truth_hdr, truth_rows, fq = None, [], {"E0": [], "E2": [], "E4": [], "E6": []}
    vtruth = {"E0": [], "E2": [], "E4": [], "E6": []}          # per-variant read truth tables (the read starts differ between variants)
    for spec in (specs if specs is not None else family_specs(seed, n_per_class, classes, lead_g=lead_g)):
        name, D = spec["name"], spec["divergence_class"]
        fdir = f"{out}/fam/{name}"
        fspec = {k: v for k, v in spec.items() if k not in ("divergence_class", "lead_g", "scenario")}
        if "lead_g" in spec:
            from famsim.template import resolve_template
            os.makedirs(fdir, exist_ok=True)
            m = resolve_template(spec["template"], _seeded(spec["seed"], "template"))
            j = spec["lead_g"]
            m.seq = "G" * j + m.seq[j:]
            m.save(f"{fdir}/template.in.json")
            fspec["template"] = {"source": "file", "path": f"{fdir}/template.in.json"}
        planted, man = chromosome.build(fspec, fdir, lambda *a, **k: None)
        FR.simulate(planted, spec["reads"], fdir, spec["seed"], lambda *a, **k: None)
        e2dir = f"{fdir}/E2"
        os.makedirs(e2dir, exist_ok=True)
        FR.simulate(planted, dict(spec["reads"], trunc5_frac=0.3, trunc5_max=0.3), e2dir, spec["seed"] + 1, lambda *a, **k: None)
        for ctg, s in read_fa(f"{fdir}/genome.truth.fa").items():
            truth[f"{name}_{ctg}"] = s
        for ctg, s in read_fa(f"{fdir}/genome.ref.fa").items():
            ref[f"{name}_{ctg}"] = s
        for p in planted:
            ctg = f"{name}_{p.contig}"
            copies.append((p.id, name, "E" if not p.in_reference else p.id[-1], D, ctg, p.pos, p.end, p.strand, len(p.model.exons), p.in_reference))
            trans[p.id] = p.model.chain_seq()
            for iso_name, skip, _w in p.isoform_list()[1:]:
                trans[f"{p.id}.{iso_name}"] = p.model.chain_seq(skip)
            if p.in_reference:
                targets[f"{name}:{p.id[-1]}"] = ref[ctg][p.pos:p.end]
        e4dir = f"{fdir}/E4"
        os.makedirs(e4dir, exist_ok=True)
        FR.simulate(planted, dict(spec["reads"], jitter=3), e4dir, spec["seed"] + 2, lambda *a, **k: None)
        e6dir = f"{fdir}/E6"
        os.makedirs(e6dir, exist_ok=True)
        FR.simulate(planted, dict(spec["reads"], jitter=1), e6dir, spec["seed"] + 3, lambda *a, **k: None)
        for tag, d in (("E0", fdir), ("E2", e2dir), ("E4", e4dir), ("E6", e6dir)):
            fq[tag] += open(f"{d}/reads.fq").read().splitlines()
            vtruth[tag] += open(f"{d}/reads.truth.tsv").read().splitlines()[1:]
        lines = open(f"{fdir}/reads.truth.tsv").read().splitlines()
        truth_hdr = lines[0]
        truth_rows += lines[1:]
        log(f"{name} D={D}: {len(planted)} copies")
    json.dump({sp["name"]: sp.get("scenario") for sp in (specs or [])}, open(f"{out}/world.json", "w"))
    write_fa(f"{out}/genome.ref.fa", ref)
    write_fa(f"{out}/genome.truth.fa", truth)
    write_fa(f"{out}/transcripts.fa", trans)
    write_fa(f"{out}/targets.fa", targets)
    with open(f"{out}/copies.tsv", "w") as o:
        o.write("copy\tfamily\trole\tD\tcontig\tpos0\tend\tstrand\tn_exons\tin_reference\n")
        for c in copies:
            o.write("\t".join(str(x) for x in c) + "\n")
    with open(f"{out}/reads.truth.tsv", "w") as o:
        o.write(truth_hdr + "\n" + "\n".join(truth_rows) + "\n")
    open(f"{out}/reads.E0.fq", "w").write("\n".join(fq["E0"]) + "\n")
    open(f"{out}/reads.E2.fq", "w").write("\n".join(fq["E2"]) + "\n")
    for tag, rows in vtruth.items():
        open(f"{out}/reads.truth.{tag}.tsv", "w").write(truth_hdr + "\n" + "\n".join(rows) + "\n")
    open(f"{out}/reads.E4.fq", "w").write("\n".join(fq["E4"]) + "\n")
    open(f"{out}/reads.E6.fq", "w").write("\n".join(fq["E6"]) + "\n")
    open(f"{out}/reads.E5.fq", "w").write("\n".join(add_g5(fq["E6"], seed)) + "\n")
    open(f"{out}/reads.E3.fq", "w").write("\n".join(add_g5(fq["E4"], seed)) + "\n")
    open(f"{out}/reads.E1.fq", "w").write("\n".join(add_polya(fq["E0"], seed)) + "\n")
    log(f"world: {len(ref)} reference contigs, {len(truth) - len(ref)} omitted (erased copies), {len(copies)} copies, {len(truth_rows)} reads")


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("cmd", choices=["build"])
    ap.add_argument("out")
    ap.add_argument("--seed", type=int, default=SEED)
    ap.add_argument("--n-per-class", type=int, default=6)
    ap.add_argument("--partition-world", action="store_true", help="Amendment 18 scenarios (CTL / ISO / SIB)")
    ap.add_argument("--lead-g", action="store_true", help="force 0-3 leading G on each template's first exon (Amendment 16)")
    ap.add_argument("--classes", default=",".join(str(c) for c in CLASSES), help="comma-separated divergence classes")
    a = ap.parse_args()
    build_world(a.out, a.seed, n_per_class=a.n_per_class, classes=tuple(float(x) for x in a.classes.split(",")), lead_g=a.lead_g, specs=partition_specs(a.seed) if a.partition_world else None)
