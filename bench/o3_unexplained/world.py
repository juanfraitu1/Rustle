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
import json
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


HARD = (
    # (scenario, sources' truth note) -- see hard_specs
    "HET_BG@0.001", "HET_BG@0.002", "HET_BG@0.004", "HET_COPY", "MANY_NULL", "MANY_REMOVED", "CONV_ABSENT", "NULL_CONV",
    "MOSAIC", "HIGHERR_NULL", "HIGHERR@0.006", "HIGHERR@0.01", "IMBALANCE", "ISOFORMS", "STRUCT_ABSENT",
    "EDIT_NULL", "EDIT_COPY", "LOW4", "SPLIT",
)


def hard_specs(seed, mixed=True):
    """The harder world (`--design hard`), 2 families per scenario (HET_BG 3). A = the template; absent sequences are E / F (copies)
    or H (an allele: the individual's second haplotype of a present copy, at that copy's expression).
      HET_BG@d       A, B 3%; H = A + d                      -- the heterozygous allele background a diploid carries everywhere
      HET_COPY       A, B 3%; H = A + 0.2%, E = A + 2%       -- three groups on one host
      MANY_NULL      8 present copies at 0, 0.5, 1, 1.5, 2, 3, 4, 5%           MANY_REMOVED  MANY_NULL + E = A + 1% absent
      CONV_ABSENT    A, B 3%; E = A + 2% with exons 2-3 converted back from A (identical tracts)
      NULL_CONV      A, B 2% with exons 2-3 converted from A; C 4%: all present
      MOSAIC         A, B 2%; E = B with exons 1..n/2 converted from A, + 0.5%: its reads split between A and B
      HIGHERR_NULL   A1 = A2, B 3%, reads err 0.005 / indel 0.002 (5x; a pessimistic stress case, above Kinnex HiFi's ~0.1%)
      HIGHERR@d      A, B 3%, E = A + d absent, same reads
      IMBALANCE      A at 200 reads, B 3%, E = A + 2% at 6 reads
      ISOFORMS       A and E (= A + 2%, absent) each with an exon-3-skipping isoform at weight 0.5; B 3%
      STRUCT_ABSENT  A, B 3%; E = A + 2% with exon 3 deleted
      EDIT_NULL      A, B 3% present; A-to-I editing (A->G in the reads) at 8 sites of A, each edited in 40% of A's reads
      EDIT_COPY      EDIT_NULL + E = A + 2% absent
      LOW4           A, B 3%; E = A + 2% at 4 reads
      SPLIT          A, B = A + 0.4% present; E = A + 1% absent: its reads split between A and B"""
    from famsim.model import seeded
    snp = lambda r: [{"op": "snp", "rate": r}] if r else []
    plan = [x for x in HARD for _ in range(3 if x.startswith("HET_BG") else 2)]
    out = []
    for k, sc in enumerate(plan):
        name = f"f{k:02d}"
        rng = seeded(seed, "family", name)
        n = rng.randint(6, 10)
        exons = [rng.randint(100, 250)] + [rng.randint(90, 350) for _ in range(n - 2)] + [rng.randint(250, 800)]
        introns = [rng.randint(300, 2500) for _ in range(n - 1)]
        kind, d = (sc.split("@")[0], float(sc.split("@")[1])) if "@" in sc else (sc, None)
        reads = {"per_copy": 40, "err": 0.001, "indel": 0.0003, "jitter": 30}
        A = ("A", snp(0), True, None, None)
        B3 = ("B", snp(0.03), True, None, None)
        if kind == "HET_BG":
            cp = [A, B3, ("H", snp(d), False, None, None)]
        elif kind == "HET_COPY":
            cp = [A, B3, ("H", snp(0.002), False, None, None), ("E", snp(0.02), False, None, None)]
        elif kind in ("MANY_NULL", "MANY_REMOVED"):
            cp = [A] + [(t, snp(r), True, None, None) for t, r in (("B", 0.005), ("C", 0.01), ("D", 0.015), ("G", 0.02), ("I", 0.03),
                                                                   ("J", 0.04), ("K", 0.05))]
            if kind == "MANY_REMOVED":
                cp.append(("E", snp(0.01), False, None, None))
        elif kind == "CONV_ABSENT":
            cp = [A, B3, ("E", snp(0.02) + [{"op": "convert", "from": name + "A", "exon": 2}, {"op": "convert", "from": name + "A", "exon": 3}],
                          False, None, None)]
        elif kind == "NULL_CONV":
            cp = [A, ("B", snp(0.02) + [{"op": "convert", "from": name + "A", "exon": 2}, {"op": "convert", "from": name + "A", "exon": 3}],
                      True, None, None), ("C", snp(0.04), True, None, None)]
        elif kind == "MOSAIC":
            conv = [{"op": "convert", "from": name + "A", "exon": i} for i in range(1, n // 2 + 1)]
            cp = [A, ("B", snp(0.02), True, None, None), ("E", snp(0.02) + conv + snp(0.005), False, None, None)]
        elif kind == "HIGHERR_NULL":
            cp = [("A1", snp(0), True, None, None), ("A2", snp(0), True, None, None), B3]
            reads = dict(reads, err=0.005, indel=0.002)
        elif kind == "HIGHERR":
            cp = [A, B3, ("E", snp(d), False, None, None)]
            reads = dict(reads, err=0.005, indel=0.002)
        elif kind == "IMBALANCE":
            cp = [("A", snp(0), True, 200, None), B3, ("E", snp(0.02), False, 6, None)]
        elif kind == "ISOFORMS":
            iso = [{"skip": [3], "weight": 0.5}]
            cp = [("A", snp(0), True, None, iso), B3, ("E", snp(0.02), False, None, iso)]
        elif kind == "STRUCT_ABSENT":
            cp = [A, B3, ("E", snp(0.02) + [{"op": "exon_delete", "exon": 3}], False, None, None)]
        elif kind == "EDIT_NULL":
            cp = [A, B3]
        elif kind == "EDIT_COPY":
            cp = [A, B3, ("E", snp(0.02), False, None, None)]
        elif kind == "LOW4":
            cp = [A, B3, ("E", snp(0.02), False, 4, None)]
        elif kind == "SPLIT":
            cp = [A, ("B", snp(0.004), True, None, None), ("E", snp(0.01), False, None, None)]
        copies, kinds = [], {}
        for tag, ops, inref, expr, iso in cp:
            c = {"id": name + tag, "ops": ops, "strand": strand_of(k, mixed)}
            if not inref:
                c["in_reference"] = False
            if expr is not None:
                c["expression"] = expr
            if iso:
                c["isoforms"] = iso
            copies.append(c)
            kinds[name + tag] = "present" if inref else ("allele" if tag == "H" else "absent")
        out.append({"name": name, "seed": seed * 1000 + k, "scenario": sc, "kinds": kinds, "divergence_class": d or 0.0,
                    "background": {"source": "random", "length": 200000}, "template": {"source": "synthetic", "exons": exons, "introns": introns},
                    "copies": copies, "reads": reads})
    return out


EDIT_SITES, EDIT_P, EDIT_CTX = 8, 0.4, 15


def apply_editing(out, seed):
    """A-to-I editing in the EDIT_* families of a hard world: 8 sites of copy A's transcript (an A with a unique 15-bp left context),
    each read of A edited A->G at each site it carries with probability 0.4, independently (seeded per read and site). Rewrites every
    reads.*.fq in place and writes edits.tsv (family, site, context)."""
    import random
    import zlib
    scen = json.load(open(f"{out}/world.json"))
    tx = read_fa_any(f"{out}/transcripts.fa")
    ctx = {}
    for fam, sc in scen.items():
        if not sc.startswith("EDIT_"):
            continue
        t = tx[fam + "A"]
        rng = random.Random(zlib.crc32(f"{seed}|{fam}".encode()))
        cand = [i for i in range(EDIT_CTX, len(t)) if t[i] == "A" and t.count(t[i - EDIT_CTX:i]) == 1]
        ctx[fam] = sorted(rng.sample(cand, EDIT_SITES))
        ctx[fam] = [t[i - EDIT_CTX:i] for i in ctx[fam]]
    with open(f"{out}/edits.tsv", "w") as fh:
        for fam, cs in ctx.items():
            for k, c in enumerate(cs):
                fh.write(f"{fam}\t{k}\t{c}\n")
    for fq in sorted(f for f in os.listdir(out) if f.startswith("reads.") and f.endswith(".fq")):
        lines = open(f"{out}/{fq}").read().splitlines()
        for i in range(0, len(lines), 4):
            name = lines[i][1:]
            fam = name[:3]
            if fam not in ctx or not name.startswith(fam + "A|"):
                continue
            seq = list(lines[i + 1])
            sq = lines[i + 1]
            for k, c in enumerate(ctx[fam]):
                j = sq.find(c)
                if j >= 0 and j + EDIT_CTX < len(seq) and seq[j + EDIT_CTX] == "A":
                    if random.Random(zlib.crc32(f"{seed}|{name}|{k}".encode())).random() < EDIT_P:
                        seq[j + EDIT_CTX] = "G"
            lines[i + 1] = "".join(seq)
        open(f"{out}/{fq}", "w").write("\n".join(lines) + "\n")


def read_fa_any(path):
    out, name = {}, None
    for ln in open(path):
        if ln.startswith(">"):
            name = ln[1:].split()[0]
            out[name] = []
        else:
            out[name].append(ln.strip())
    return {k: "".join(v) for k, v in out.items()}


def cmd_build(a):
    from synth_world import build_world
    sp = hard_specs(a.seed, a.strands) if a.design == "hard" else specs(a.seed, a.strands)
    build_world(a.out, a.seed, log=lambda *x, **k: None, specs=sp)
    if a.design == "hard":
        apply_editing(a.out, a.seed)
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
    p.add_argument("--design", choices=("base", "hard"), default="base", help="hard: hard_specs (allele background, 8 copies, conversion, ...)")
    p.set_defaults(func=cmd_build)
    p = sub.add_parser("run"); p.add_argument("out"); p.add_argument("--variant", choices=("E0", "E2"), required=True)
    p.add_argument("--arm", choices=("T", "P"), required=True); p.add_argument("--threads", type=int, default=4); p.add_argument("--bin")
    p.set_defaults(func=cmd_run)
    a = ap.parse_args()
    a.func(a)
