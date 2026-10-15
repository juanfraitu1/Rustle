#!/usr/bin/env python3
"""Structure worlds for O3 (docs/superpowers/specs/2026-10-14-o3-structure-projection-design.md): every family has two copies in the
reference (A, B) and one copy E absent from it, whose reads `candidate_copies` must turn into a candidate with E's sequence AND E's exon
structure. E differs from its nearest present copy by SNPs at a divergence ladder, or carries a structural change (exon loss, gain,
duplication, shuffle, 5' truncation, a resized intron) at 2%, or is lowly expressed. NULL families (nothing missing) are the control.

    world.py build OUT --seed S          # synthetic world + the stage's inputs (same layout as bench/o3_fewread/world.py)
    world.py run OUT --variant E0|E2|E3  # map the variant's reads, run candidate_copies -> OUT/<variant>/P.cand.*

Seeds 11 and 12 are development worlds; the held-out seeds are fixed in the pre-registration.
"""
import argparse
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
BENCH = os.path.dirname(HERE)
sys.path.insert(0, BENCH)
sys.path.insert(0, os.path.join(BENCH, "unmapped_rescue"))
sys.path.insert(0, os.path.join(BENCH, "o3_fewread"))
import world as FW  # noqa: E402  (bench/o3_fewread/world.py: write_copies, run)

DIV = (0.005, 0.01, 0.02, 0.04, 0.08)


def plan():
    """(scenario, E's ops, E's expression or None, divergence of A from the template) per family, in family order"""
    out = [(f"DIV@{d}", [], None, d) for d in DIV for _ in range(2)]
    out += [("LOSS", [{"op": "exon_delete", "exon": k}], None, 0.02) for k in (2, 3, 4)]
    out += [("GAIN", [{"op": "exon_insert", "after": k, "length": 120}], None, 0.02) for k in (1, 2, 3)]
    out += [("DUP", [{"op": "exon_insert", "after": k, "source_exon": k}], None, 0.02) for k in (2, 3)]
    out += [("SHUFFLE", [{"op": "exon_shuffle", "a": 2, "b": 3}], None, 0.02), ("SHUFFLE", [{"op": "exon_shuffle", "a": 2, "b": 4}], None, 0.02)]
    out += [("TRUNC5", [{"op": "truncate", "side": 5, "exons": k}], None, 0.02) for k in (1, 2)]
    out += [("INTRON", [{"op": "intron_resize", "intron": 2, "length": L}], None, 0.02) for L in (150, 4000)]
    out += [("LOWEXPR", [], n, 0.02) for n in (8, 12, 16)]
    out += [("NULL", None, None, 0.02) for _ in range(4)]
    return out


def strand_of(k, tag, mixed):
    """--strands: family k % 4 == 1 has every copy on -, k % 4 == 3 has A on + and B, E (or C) on -; the rest +"""
    if not mixed:
        return "+"
    return "-" if k % 4 == 1 or (k % 4 == 3 and tag != "A") else "+"


def specs(seed, mixed=False):
    from famsim.model import seeded
    out = []
    for k, (sc, e_ops, e_expr, d) in enumerate(plan()):
        name = f"f{k:02d}"
        rng = seeded(seed, "family", name)
        n = rng.randint(6, 10)                      # >= 6 exons: every op above has an internal exon to act on
        exons = [rng.randint(100, 250)] + [rng.randint(90, 350) for _ in range(n - 2)] + [rng.randint(250, 800)]
        introns = [rng.randint(300, 2500) for _ in range(n - 1)]
        copies = [{"id": name + "A", "ops": [{"op": "snp", "rate": d}], "strand": strand_of(k, "A", mixed)},
                  {"id": name + "B", "ops": [{"op": "snp", "rate": d + 0.01}], "strand": strand_of(k, "B", mixed)}]
        kinds = {name + "A": "present", name + "B": "present"}
        if e_ops is None:                           # NULL: a third present copy instead of an absent one
            copies.append({"id": name + "C", "ops": [{"op": "snp", "rate": 0.03}], "strand": strand_of(k, "C", mixed)})
            kinds[name + "C"] = "present"
        else:
            e = {"id": name + "E", "in_reference": False, "ops": list(e_ops), "strand": strand_of(k, "E", mixed)}
            if e_expr is not None:
                e["expression"] = e_expr
            copies.append(e)
            kinds[name + "E"] = "copy"
        out.append({"name": name, "seed": seed * 1000 + k, "scenario": sc, "kinds": kinds, "divergence_class": d,
                    "background": {"source": "random", "length": 80000}, "template": {"source": "synthetic", "exons": exons, "introns": introns},
                    "copies": copies, "reads": {"per_copy": 40, "err": 0.001, "indel": 0.0003, "jitter": 30}})
    return out


def cmd_build(a):
    from synth_world import build_world
    build_world(a.out, a.seed, log=lambda *x, **k: None, specs=specs(a.seed, a.strands))
    n = FW.write_copies(a.out)
    FW.sh(f"samtools faidx {a.out}/genome.ref.fa && minimap2 -x splice:hq -d {a.out}/genome.ref.splice.mmi {a.out}/genome.ref.fa 2>/dev/null")
    print(f"[o3_structure] world {a.out} seed {a.seed}: {n} in-reference copies", flush=True)


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd", required=True)
    p = sub.add_parser("build"); p.add_argument("out"); p.add_argument("--seed", type=int, required=True)
    p.add_argument("--strands", action="store_true", help="mixed copy orientations (strand_of); the dev worlds 11 and 12 were built without")
    p.set_defaults(func=cmd_build)
    p = sub.add_parser("run"); p.add_argument("out"); p.add_argument("--variant", choices=("E0", "E2", "E3"), required=True)
    p.add_argument("--threads", type=int, default=4); p.add_argument("--bin"); p.set_defaults(func=FW.cmd_run)
    a = ap.parse_args()
    a.func(a)
