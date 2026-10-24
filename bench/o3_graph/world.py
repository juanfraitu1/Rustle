#!/usr/bin/env python3
"""Ideal-scenario worlds for the graph-O3 detector (prereg docs/archive/2026-10/merged/PREREG_o3_graph_sim_2026-10-15.md): one family per (curated
template, scenario), every copy expressing the template's curated isoforms (RefSeq / CAT; templates.py), reads by famsim.

    world.py build OUT TEMPLATES.json --seed S --scenario SC [--first 0 --n 12] [--per-copy 60]

One world per scenario, one family per template: every family's copy A IS its template, so two families of one template in a
world would be copies of each other (reads of one family's A align perfectly to the other's).

Scenarios (copies; A = the template, B = A + 3% SNPs, both present in the genome and the catalog):
  NULL          A, B
  ALLELE        A, B, H = A + 0.4% (absent: an allele the reference lacks)
  ABSENT@d      A, B, E = A + d (absent)
  BRANCH        A, B, E = A + a novel 150-bp exon in the long intron nearest the middle + 0.2% SNPs (absent): structure, not
                sequence, differs
  OUTSIDE@d     A, B, D = A + d in the GENOME but not in the catalog (an uncatalogued paralog)
  TWIN_ALT      A, T = A + 15% SNPs in the first cassette exon only, B; the representative skips that exon, so A and T look
                identical to a representative-based method
Isoforms: the canonical (iso0, famsim weight 1.0), the MAIN isoform (skips the first cassette exon; weight MAIN_WEIGHT, so it is
each copy's most expressed isoform and the catalog's representative) and up to two more curated ones (OTHER_WEIGHT).
Products besides synth_world's: P.fam.copies.tsv / .fa (the catalog: present, catalogued copies; representative = most
expressed isoform), isoforms.json (every catalogued copy's curated isoforms in reference coordinates), kinds.json.
"""
import argparse
import collections
import json
import os
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
BENCH = os.path.dirname(HERE)
for p in (BENCH, os.path.join(BENCH, "unmapped_rescue"), os.path.join(BENCH, "o3_fewread")):
    sys.path.insert(0, p)

SCENARIOS = ("NULL", "ALLELE", "ABSENT@0.01", "ABSENT@0.02", "BRANCH", "OUTSIDE@0.01", "OUTSIDE@0.03", "TWIN_ALT")
MAIN_WEIGHT, OTHER_WEIGHT, ISO_CAP = 2.0, 0.5, 4
TWIN_EXON_RATE, ALLELE_RATE, BRANCH_SNP, BRANCH_LEN, INTRON_CAP = 0.15, 0.004, 0.002, 150, 2000


def curated_isoforms(tpl, cap=ISO_CAP - 1):
    """The famsim isoforms of a template beside iso0: the main one first (skips the first cassette exon, fewest skips), then
    the template's others in order, `cap` in all."""
    k = tpl["cassettes"][0]
    isos = [x for x in tpl["isoforms"]]
    main = min((x for x in isos if k in x["skip"]), key=lambda x: (len(x["skip"]), isos.index(x)))
    out = [{"skip": list(main["skip"]), "weight": MAIN_WEIGHT}]
    for x in isos:
        if x is not main and len(out) < cap:
            out.append({"skip": list(x["skip"]), "weight": OTHER_WEIGHT})
    return out


def shift_skips(isos, after):
    """Skip indices after a newly inserted exon (1-based position `after` + 1) move up by one."""
    return [dict(x, skip=[k + 1 if k > after else k for k in x["skip"]]) for x in isos]


def rep_isoform(isos):
    """The most expressed isoform's famsim name (iso0 weighs 1.0; a tie keeps iso0)."""
    best, name = 1.0, "iso0"
    for i, x in enumerate(isos):
        if x["weight"] > best:
            best, name = x["weight"], f"iso{i + 1}"
    return name


def branch_after(tpl, min_len=300):
    """The exon (1-based, transcript orientation) after which BRANCH's novel exon goes: the intron nearest the middle of the
    intron list among those >= min_len bp (famsim needs >= 100 bp around it), else the longest intron."""
    introns = tpl["introns"]
    centre = (len(introns) + 1) / 2
    ok = [i for i, n in enumerate(introns, 1) if n >= min_len]
    if not ok:
        return max(range(1, len(introns) + 1), key=lambda i: (introns[i - 1], -i))
    return min(ok, key=lambda i: (abs(i - centre), i))


def specs(templates, seed, genome, annotation, per_copy=60):
    snp = lambda r: [{"op": "snp", "rate": r}]  # noqa: E731
    out = []
    for ti, t in enumerate(templates):
        iso = curated_isoforms(t)
        strand = "-" if ti % 2 else "+"
        after = branch_after(t) if "introns" in t else t["n_exons"] // 2 + 1
        k = t["cassettes"][0]
        for si, sc in enumerate(SCENARIOS):
            kind, _, d = sc.partition("@")
            d = float(d) if d else None
            name = f"{kind}{'' if d is None else str(d).replace('0.', '')}_{t['gene']}".replace("-", "")
            A = ("A", [], True, "present", iso)
            B = ("B", snp(0.03), True, "present", iso)
            cp = {"NULL": [A, B],
                  "ALLELE": [A, B, ("H", snp(ALLELE_RATE), False, "allele", iso)],
                  "ABSENT": [A, B, ("E", snp(d or 0), False, "absent", iso)],
                  "BRANCH": [A, B, ("E", [{"op": "exon_insert", "after": after, "length": BRANCH_LEN}] + snp(BRANCH_SNP), False,
                                    "absent", shift_skips(iso, after))],
                  "OUTSIDE": [A, B, ("D", snp(d or 0), True, "outside", iso)],
                  "TWIN_ALT": [A, ("T", [{"op": "snp", "rate": TWIN_EXON_RATE, "region": "exon", "exon": k}], True, "present", iso), B],
                  }[kind]
            copies, kinds = [], {}
            for tag, ops, inref, kd, isos in cp:
                c = {"id": name + tag, "ops": ops, "strand": strand, "isoforms": isos}
                if not inref:
                    c["in_reference"] = False
                copies.append(c)
                kinds[name + tag] = kd
            out.append({"name": name, "seed": seed * 1000 + ti * len(SCENARIOS) + si, "scenario": sc, "kinds": kinds,
                        "divergence_class": d or 0.0, "background": {"source": "random", "length": 60000},
                        "template": {"source": "annotation", "genome": genome, "annotation": annotation, "gene": t["gene"],
                                     "chrom": t["chrom"], "transcript": t["canonical"], "intron_cap": INTRON_CAP},
                        "copies": copies, "reads": {"per_copy": per_copy, "err": 0.001, "indel": 0.0003, "jitter": 30}})
    return out


def truth_isoforms(gtf, copy):
    """{isoN: exons (0-based half-open, truth coordinates)} of a copy in a famsim truth.gtf, and iso0's 1-based start."""
    ex, tstart = collections.defaultdict(list), None
    for ln in open(gtf):
        r = ln.rstrip("\n").split("\t")
        if len(r) < 9 or f'gene_id "{copy}"' not in r[8]:
            continue
        tid = r[8].split('transcript_id "')[1].split('"')[0]
        iso = tid.rsplit(".", 1)[1]
        if r[2] == "transcript" and iso == "iso0":
            tstart = int(r[3])
        elif r[2] == "exon":
            ex[iso].append((int(r[3]) - 1, int(r[4])))
    return {k: sorted(v) for k, v in ex.items()}, tstart


def write_catalog(out, sp):
    """P.fam.copies.tsv / .fa (present, catalogued copies; exons = the most expressed isoform) and isoforms.json."""
    import world as FW  # bench/o3_fewread/world.py
    from synth_world import read_fa
    ref = read_fa(f"{out}/genome.ref.fa")
    by_copy = {c["id"]: (s, c) for s in sp for c in s["copies"]}
    rows, fa, isoj = [], [], collections.defaultdict(list)
    idx = collections.Counter()
    for ln in list(open(f"{out}/copies.tsv"))[1:]:
        copy, fam, _role, _D, ctg, pos0, end, strand, _n, inref = ln.rstrip("\n").split("\t")
        s, c = by_copy[copy]
        if inref != "True" or s["kinds"][copy] != "present":
            continue
        isos, tstart = truth_isoforms(f"{out}/fam/{fam}/truth.gtf", copy)
        off = int(pos0) - (tstart - 1)
        isos = {k: [(a + off, b + off) for a, b in v] for k, v in isos.items()}
        rep = rep_isoform(c["isoforms"])
        ex = isos[rep]
        seq = "".join(ref[ctg][a:b] for a, b in ex)
        i = idx[fam]
        idx[fam] += 1
        rows.append([fam, str(i), f"SIM_{copy}", ctg, str(ex[0][0]), str(ex[-1][1]), str(len(ex)), strand, "0",
                     ",".join(f"{a}-{b}" for a, b in ex), "1.000", "famsim", copy, "NA", "1", str(len(seq)), "0.000", "kept",
                     pos0, end])
        fa.append(f">{fam}|{i}|{ctg}:{ex[0][0]}-{ex[-1][1]}|{strand}|nexon={len(ex)}\n{seq}\n")
        isoj[fam].append(dict(copy=copy, idx=i, contig=ctg, start=int(pos0), end=int(end), strand=strand, rep=rep,
                              isoforms={k: v for k, v in sorted(isos.items())}))
    with open(f"{out}/P.fam.copies.tsv", "w") as f:
        f.write(FW.COPIES_HEADER + "\n" + "".join("\t".join(r) + "\n" for r in rows))
    open(f"{out}/P.fam.copies.fa", "w").write("".join(fa))
    json.dump(isoj, open(f"{out}/isoforms.json", "w"), indent=1)
    return len(rows)


def cmd_build(a):
    from synth_world import build_world
    tj = json.load(open(a.templates))
    tpls = tj["templates"][a.first:a.first + a.n]
    sp = [s for s in specs(tpls, a.seed, tj["genome"], tj["annotation"], a.per_copy) if s["scenario"] == a.scenario]
    os.makedirs(a.out, exist_ok=True)
    json.dump(sp, open(f"{a.out}/specs.json", "w"), indent=1)
    build_world(a.out, a.seed, log=lambda *x, **k: None, specs=sp)
    n = write_catalog(a.out, sp)
    subprocess.run(f"samtools faidx {a.out}/genome.ref.fa", shell=True, check=True)
    print(f"[o3_graph] world {a.out}: {len(sp)} families from {len(tpls)} templates, {n} catalogued copies", flush=True)


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd", required=True)
    p = sub.add_parser("build")
    p.add_argument("out")
    p.add_argument("templates")
    p.add_argument("--seed", type=int, required=True)
    p.add_argument("--scenario", required=True, choices=SCENARIOS)
    p.add_argument("--first", type=int, default=0)
    p.add_argument("--n", type=int, default=12)
    p.add_argument("--per-copy", type=int, default=60)
    p.set_defaults(func=cmd_build)
    a = ap.parse_args(argv)
    a.func(a)
    return 0


if __name__ == "__main__":
    sys.exit(main())
