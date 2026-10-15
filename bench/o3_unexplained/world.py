#!/usr/bin/env python3
"""Worlds for O3 as the advisor defines it (docs/superpowers/specs/2026-10-15-o3-unexplained-multimappers-design.md): reads of a copy
absent from the GENOME land on its paralogues; does O2's origin certificate find no adequate copy for them, and does
`copy_assign --flag-missing-copies` cluster them?

    world.py build OUT --seed S
    world.py run OUT --variant E0|E2 --arm T|P      # map reads (sim.MM2), copy_assign --families --flag-missing-copies -> OUT/<variant>.<arm>/

Scenarios: see `specs` (reads 40 per expressed copy unless noted). The copy E (or allele H) is absent from the GENOME, so its reads
land on its paralogues: as exact ties on identical twins (REMOVED_TWIN, arm T sees them) or uniquely on its nearest copy (REMOVED_UNIQ,
REMOVED_MULTI, LOWEXPR: arm P only).
Arm T = the shipped contested set (exact AS ties); arm P = every molecule placed on a family copy (`--no-as-tied-only`).
Seeds 31-33 are development worlds; held-out seeds are fixed in the pre-registration.
"""
import argparse
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
BENCH = os.path.dirname(HERE)
sys.path.insert(0, BENCH)
sys.path.insert(0, os.path.join(BENCH, "unmapped_rescue"))
sys.path.insert(0, os.path.join(BENCH, "o3_fewread"))
import world as FW  # noqa: E402  (bench/o3_fewread/world.py: write_copies, sh)

DIV = (0.003, 0.006, 0.01, 0.02, 0.04)
ALLELE_D = (0.0025, 0.005, 0.008)


def plan():
    """(scenario, d, E's expression or None) per family, in family order"""
    out = [("NULL_TWIN", None, None)] * 3 + [("NULL_MULTI", None, None)] * 3 + [("NULL_CLOSE", None, None)] * 2
    out += [(f"REMOVED_TWIN@{d}", d, None) for d in DIV for _ in range(2)]
    out += [(f"REMOVED_UNIQ@{d}", d, None) for d in DIV for _ in range(2)]
    out += [(f"REMOVED_MULTI@{d}", d, None) for d in (0.006, 0.01, 0.02) for _ in range(2)]
    out += [("LOWEXPR@0.02", 0.02, n) for n in (8, 12, 16)]
    out += [(f"ALLELE@{d}", d, None) for d in ALLELE_D for _ in range(2)]
    out += [(f"ADDBACK@{d}", d, None) for d in (0.01, 0.02)] + [("DECOY", None, None)] * 2
    out += [("COMBO_EH", 0.02, None)] * 2 + [("DOUBLE", 0.01, None)] * 2
    return out


def strand_of(k, mixed):
    return "-" if mixed and k % 4 == 1 else "+"


def specs(seed, mixed=False):
    """Copies per scenario (A = the template unless noted; E / H absent from the genome):
      NULL_TWIN       A1 = A2, B 3%                 NULL_MULTI   A, B 1.5%, C 3%, D 5%        NULL_CLOSE  A, B 0.3%
      REMOVED_TWIN@d  A1 = A2, B 3%, E = A + d      REMOVED_UNIQ@d  A, B 3%, E = A + d
      REMOVED_MULTI@d A, B 1.5%, C 3%, D 5%, E = A + d                LOWEXPR@0.02  as REMOVED_UNIQ, E at 8-16 reads
      ALLELE@d        A, B 3%, H = A + d at A's expression (a heterozygous allele the reference lacks)
      ADDBACK@d       REMOVED_TWIN@d with E in the reference          DECOY  NULL_TWIN + C 5% in the reference, not expressed
      COMBO_EH        A, B 3%; absent copy E = A + 2% AND allele H = A + 0.5%: two groups absorbed onto A
      DOUBLE          A, B 4%; two absent copies E = A + 1%, F = A + 2.5%"""
    from famsim.model import seeded
    snp = lambda r: [{"op": "snp", "rate": r}] if r else []
    out = []
    for k, (sc, d, e_expr) in enumerate(plan()):
        name = f"f{k:02d}"
        rng = seeded(seed, "family", name)
        n = rng.randint(6, 10)
        exons = [rng.randint(100, 250)] + [rng.randint(90, 350) for _ in range(n - 2)] + [rng.randint(250, 800)]
        introns = [rng.randint(300, 2500) for _ in range(n - 1)]
        kind = sc.split("@")[0]
        twin = [("A1", 0, True, None), ("A2", 0, True, None), ("B", 0.03, True, None)]
        multi = [("A", 0, True, None), ("B", 0.015, True, None), ("C", 0.03, True, None), ("D", 0.05, True, None)]
        cp = {"NULL_TWIN": twin, "NULL_MULTI": multi, "NULL_CLOSE": [("A", 0, True, None), ("B", 0.003, True, None)],
              "REMOVED_TWIN": twin + [("E", d, False, None)], "REMOVED_UNIQ": [("A", 0, True, None), ("B", 0.03, True, None), ("E", d, False, None)],
              "REMOVED_MULTI": multi + [("E", d, False, None)],
              "LOWEXPR": [("A", 0, True, None), ("B", 0.03, True, None), ("E", d, False, e_expr)],
              "ALLELE": [("A", 0, True, None), ("B", 0.03, True, None), ("H", d, False, None)],
              "ADDBACK": twin + [("E", d, True, None)], "DECOY": twin + [("C", 0.05, True, 0)],
              "COMBO_EH": [("A", 0, True, None), ("B", 0.03, True, None), ("E", 0.02, False, None), ("H", 0.005, False, None)],
              "DOUBLE": [("A", 0, True, None), ("B", 0.04, True, None), ("E", 0.01, False, None), ("F", 0.025, False, None)]}[kind]
        copies, kinds = [], {}
        for tag, r, inref, expr in cp:
            c = {"id": name + tag, "ops": snp(r), "strand": strand_of(k, mixed)}
            if not inref:
                c["in_reference"] = False
            if expr is not None:
                c["expression"] = expr
            copies.append(c)
            kinds[name + tag] = "present" if inref else ("allele" if tag == "H" else "absent")
        out.append({"name": name, "seed": seed * 1000 + k, "scenario": sc, "kinds": kinds, "divergence_class": d or 0.0,
                    "background": {"source": "random", "length": 120000}, "template": {"source": "synthetic", "exons": exons, "introns": introns},
                    "copies": copies, "reads": {"per_copy": 40, "err": 0.001, "indel": 0.0003, "jitter": 30}})
    return out


def cmd_build(a):
    from synth_world import build_world
    build_world(a.out, a.seed, log=lambda *x, **k: None, specs=specs(a.seed, a.strands))
    n = FW.write_copies(a.out)
    FW.sh(f"samtools faidx {a.out}/genome.ref.fa")
    print(f"[o3_unexplained] world {a.out} seed {a.seed}: {n} in-reference copies", flush=True)


def cmd_run(a):
    import sim
    w, v = a.out, a.variant
    d = f"{w}/{v}.{a.arm}"
    os.makedirs(d, exist_ok=True)
    bam = f"{w}/{v}.reads.bam"
    if not os.path.exists(bam + ".bai"):
        FW.sh(f"{sim.MM2} -t {a.threads} {w}/genome.ref.fa {w}/reads.{v}.fq 2>/dev/null | samtools sort -@2 -o {bam} - 2>/dev/null && samtools index {bam}")
    regions = f"{w}/regions.txt"
    FW.sh(f"samtools view -H {bam} | awk '$1==\"@SQ\"{{sub(\"SN:\",\"\",$2); sub(\"LN:\",\"\",$3); print $2\":1-\"$3}}' > {regions}")
    binary = a.bin or os.path.join(os.environ.get("CARGO_TARGET_DIR", os.path.expanduser("~/rustle_target")), "release", "copy_assign")
    extra = " --no-as-tied-only" if a.arm == "P" else ""
    FW.sh(f"{binary} --bam {bam} --fasta {w}/genome.ref.fa --regions {regions} --families {w}/P.fam.copies.tsv --copies-fa {w}/P.fam.copies.fa "
          f"--threads {a.threads} --flag-missing-copies{extra} --out {d}/A > {d}/A.log 2>&1")
    print(f"[o3_unexplained] {d}: done", flush=True)


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd", required=True)
    p = sub.add_parser("build"); p.add_argument("out"); p.add_argument("--seed", type=int, required=True)
    p.add_argument("--strands", action="store_true", help="families k % 4 == 1 on the - strand")
    p.set_defaults(func=cmd_build)
    p = sub.add_parser("run"); p.add_argument("out"); p.add_argument("--variant", choices=("E0", "E2"), required=True)
    p.add_argument("--arm", choices=("T", "P"), required=True); p.add_argument("--threads", type=int, default=4); p.add_argument("--bin")
    p.set_defaults(func=cmd_run)
    a = ap.parse_args()
    a.func(a)
