#!/usr/bin/env python3
"""Phased-data toy — docs/archive/2026-10/merged/PREREG_phased_toy_2026-10-12.md.

One real gene (SNRPB, CHM13 chr20, CAT) on a 30-kb real slice, two haplotypes per scenario:
  S0 control  hap2 = hap1
  SA allele   hap2's only copy carries V (heterozygous SNVs)
  SB copy     hap2 also carries an extra copy B = unit + V, inserted downstream; hap2's own copy = hap1's
V and the per-sequence read counts are the same in SA and SB, so a collapsed (hap1-only) reference receives
identically distributed reads in both. The phased arm maps to each haplotype separately (as with mat/pat).

    python3 bench/phased_toy/phased_toy.py [--out DIR] [--seed 7]
"""
import argparse, json, os, random, re, subprocess

CHR20 = "/mnt/linuxdisk/home/juanfraitu/bakeoff/human_chr20/chr20.fa"
GTF = "/mnt/linuxdisk/home/juanfraitu/bakeoff/human_chr20/chr20_ref_cat.gtf"
TX = "CHM13_T0135794"
SLICE = ("chr20", 2481001, 2511000)  # 1-based inclusive
UNIT_PAD, INSERT_AT, FLANK = 500, 2506000, 2000
N_EXONIC_V, V_RATE = 5, 0.005
ERR, INDEL, JITTER = 0.001, 0.0003, 30
PLAN = {"S0": [(1, "A", 50), (2, "A", 50)],
        "SA": [(1, "A", 50), (2, "A", 50)],
        "SB": [(1, "A", 25), (2, "A", 25), (2, "B", 50)]}
COMP = str.maketrans("ACGTN", "TGCAN")
CS = re.compile(r"(:\d+|\*[a-z]{2}|\+[a-z]+|-[a-z]+|~[a-z]{2}\d+[a-z]{2})")


def rc(s):
    return s.translate(COMP)[::-1]


def sh(cmd, out=None):
    with open(out, "w") if out else open(os.devnull, "w") as fh:
        r = subprocess.run(cmd, stdout=subprocess.PIPE if out is None else fh, stderr=subprocess.DEVNULL, text=True, check=True)
    return r.stdout


def slice_seq():
    c, s, e = SLICE
    return "".join(sh(["samtools", "faidx", CHR20, f"{c}:{s}-{e}"]).split("\n")[1:]).upper()


def exons():
    ex, strand = [], None
    for line in open(GTF):
        f = line.split("\t")
        if len(f) > 8 and f[2] == "exon" and f'transcript_id "{TX}"' in f[8]:
            ex.append((int(f[3]), int(f[4])))
            strand = f[6]
    return sorted(ex), strand


def mutate(seq, var):
    s = list(seq)
    for p, alt in var.items():
        s[p] = alt
    return "".join(s)


def build(rng):
    """Haplotypes per scenario, copy offsets, the unit, exons relative to the unit, V (unit coords -> alt)."""
    hap1 = slice_seq()
    ex, strand = exons()
    us = ex[0][0] - UNIT_PAD - SLICE[1]          # unit start, 0-based in slice
    ue = ex[-1][1] + UNIT_PAD - SLICE[1] + 1
    unit = hap1[us:ue]
    rex = [(s - SLICE[1] - us, e - SLICE[1] - us + 1) for s, e in ex]   # half-open, unit coords
    exonic = [p for s, e in rex for p in range(s, e)]
    ex_set = set(exonic)
    other = [p for p in range(len(unit)) if p not in ex_set]
    pos = rng.sample(exonic, N_EXONIC_V) + [p for p in other if rng.random() < V_RATE]
    var = {p: rng.choice([b for b in "ACGT" if b != unit[p]]) for p in pos}
    ins = INSERT_AT - SLICE[1]
    bseq = mutate(unit, var)
    haps = {"S0": (hap1, hap1),
            "SA": (hap1, hap1[:us] + mutate(unit, var) + hap1[ue:]),
            "SB": (hap1, hap1[:ins] + bseq + hap1[ins:])}
    offs = {"A": us, "B": ins}                    # unit start of each copy in the haplotype that carries it
    return haps, offs, unit, rex, strand, var, us


def transcript(hap, off, rex, strand):
    t = "".join(hap[off + s: off + e] for s, e in rex)
    return rc(t) if strand == "-" else t


def tx_positions(rex, strand, var):
    """Transcript coordinates of the exonic V positions."""
    out, k = [], 0
    for s, e in rex:
        out += [k + p - s for p in var if s <= p < e]
        k += e - s
    n = k
    return sorted(n - 1 - p for p in out) if strand == "-" else sorted(out)


def noisy(t, rng):
    out = []
    for b in t:
        r = rng.random()
        if r < ERR:
            out.append(rng.choice([x for x in "ACGT" if x != b]))
        elif r < ERR + INDEL / 2:
            continue
        elif r < ERR + INDEL:
            out += [b, rng.choice("ACGT")]
        else:
            out.append(b)
    return "".join(out)


def simulate(scen, haps, offs, rex, strand, vtx, rng, fa):
    reads = {}
    with open(fa, "w") as fh:
        for hap, copy, n in PLAN[scen]:
            t = transcript(haps[hap - 1], offs[copy], rex, strand)
            for i in range(n):
                a, b = rng.randint(0, JITTER), rng.randint(0, JITTER)
                name = f"{scen}|h{hap}|{copy}|{i}"
                reads[name] = {"hap": hap, "copy": copy, "informative": any(a <= p < len(t) - b for p in vtx)}
                fh.write(f">{name}\n{noisy(t[a:len(t) - b], rng)}\n")
    return reads


def write_fa(path, name, seq):
    with open(path, "w") as fh:
        fh.write(f">{name}\n{seq}\n")


def map_reads(ref, fa):
    """Primary alignments: name -> {pos, end, nm, subs: {pos: alt}, blocks}."""
    out = {}
    for line in sh(["minimap2", "-ax", "splice:hq", "-uf", "--secondary=no", "--cs", "-t", "1", ref, fa]).splitlines():
        if line.startswith("@"):
            continue
        f = line.split("\t")
        if int(f[1]) & 0x904:
            continue
        tags = dict((t[:2], t[5:]) for t in f[11:])
        pos, subs, blocks = int(f[3]) - 1, {}, []
        for tok in CS.findall(tags["cs"]):
            if tok[0] == ":":
                blocks.append((pos, pos + int(tok[1:]))); pos += int(tok[1:])
            elif tok[0] == "*":
                subs[pos] = tok[2].upper(); blocks.append((pos, pos + 1)); pos += 1
            elif tok[0] == "-":
                pos += len(tok) - 1
            elif tok[0] == "~":
                pos += int(tok[3:-2])
        out[f[0]] = {"pos": int(f[3]) - 1, "end": pos, "nm": int(tags["NM"]), "subs": subs, "blocks": blocks}
    return out


def collapsed_summary(recs):
    cover, alt = {}, {}
    for r in recs.values():
        for s, e in r["blocks"]:
            for c in range(s, e):
                cover[c] = cover.get(c, 0) + 1
        for c, a in r["subs"].items():
            alt[(c, a)] = alt.get((c, a), 0) + 1
    cols = sorted({c for (c, a), n in alt.items() if cover[c] >= 10 and n / cover[c] >= 0.2})
    carriers = [n for n, r in recs.items() if any(c in r["subs"] for c in cols)]
    return {"mapped": len(recs), "variable_columns": cols, "alt_fraction": round(len(carriers) / max(len(recs), 1), 3)}


def has_partner(home_seq, other_ref, s, e, work):
    qs, qe = max(0, s - FLANK), min(len(home_seq), e + FLANK)
    q = os.path.join(work, "locus.fa")
    write_fa(q, "locus", home_seq[qs:qe])
    qlen = qe - qs
    for line in sh(["minimap2", "-cx", "asm10", "-t", "1", other_ref, q]).splitlines():
        f = line.split("\t")
        a, b, ts, te = int(f[2]), int(f[3]), int(f[7]), int(f[8])
        if a <= 100 and b >= qlen - 100 and abs((te - ts) - (b - a)) <= 0.05 * qlen:
            return True
    return False


def merge(spans):
    out = []
    for s, e in sorted(spans):
        if out and s <= out[-1][1]:
            out[-1][1] = max(out[-1][1], e)
        else:
            out.append([s, e])
    return out


def phased(rec1, rec2, haps, refs, work):
    names = set(rec1) | set(rec2)
    call, home = {}, {}
    for n in names:
        if n not in rec1 and n not in rec2:
            continue
        n1 = rec1[n]["nm"] if n in rec1 else 10 ** 9
        n2 = rec2[n]["nm"] if n in rec2 else 10 ** 9
        if n1 == n2:
            call[n] = "SHARED"
        else:
            home[n] = 1 if n1 < n2 else 2
    loci = []
    for h, rec in ((1, rec1), (2, rec2)):
        for s, e in merge([(rec[n]["pos"], rec[n]["end"]) for n in home if home[n] == h]):
            loci.append({"hap": h, "start": s, "end": e,
                         "partner": has_partner(haps[h - 1], refs[2 - h], s, e, work)})
    for n, h in home.items():
        r = (rec1 if h == 1 else rec2)[n]
        L = next(l for l in loci if l["hap"] == h and l["start"] <= r["pos"] < l["end"])
        call[n] = "ALLELE" if L["partner"] else "EXTRA_COPY"
    return call, home, loci


def verify(haps, unit, var, us, rex, strand, offs):
    h1, a2 = haps["SA"]
    assert haps["S0"][0] == haps["S0"][1]
    assert len(h1) == len(a2) and sum(x != y for x, y in zip(h1, a2)) == len(var)
    b2 = haps["SB"][1]
    assert len(b2) == len(h1) + len(unit) and b2[offs["B"]:offs["B"] + len(unit)] == mutate(unit, var)
    t1, t2 = transcript(h1, us, rex, strand), transcript(a2, us, rex, strand)
    assert sum(x != y for x, y in zip(t1, t2)) == N_EXONIC_V
    return {"slice_bp": len(h1), "unit_bp": len(unit), "V": len(var), "V_exonic": N_EXONIC_V,
            "unit_divergence": round(len(var) / len(unit), 4), "mrna_bp": len(t1)}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default="/mnt/linuxdisk/tmp/phased_toy")
    ap.add_argument("--seed", type=int, default=7)
    a = ap.parse_args()
    os.makedirs(a.out, exist_ok=True)
    rng = random.Random(a.seed)
    haps_all, offs, unit, rex, strand, var, us = build(rng)
    res = {"seed": a.seed, "verify": verify(haps_all, unit, var, us, rex, strand, offs)}
    vtx = tx_positions(rex, strand, var)
    for scen in ("S0", "SA", "SB"):
        haps = haps_all[scen]
        refs = []
        for h in (1, 2):
            p = os.path.join(a.out, f"{scen}.hap{h}.fa"); write_fa(p, f"hap{h}", haps[h - 1]); refs.append(p)
        fa = os.path.join(a.out, f"{scen}.reads.fa")
        truth = simulate(scen, haps, offs, rex, strand, vtx, rng, fa)
        rec1, rec2 = map_reads(refs[0], fa), map_reads(refs[1], fa)
        call, home, loci = phased(rec1, rec2, haps, refs, a.out)
        rows = [{"read": n, **t, "call": call.get(n, "UNMAPPED"), "home": home.get(n, 0)} for n, t in truth.items()]
        with open(os.path.join(a.out, f"{scen}.calls.tsv"), "w") as fh:
            fh.write("read\thap\tcopy\tinformative\tcall\thome\n")
            for r in rows:
                fh.write("\t".join(str(r[k]) for k in ("read", "hap", "copy", "informative", "call", "home")) + "\n")
        res[scen] = {"collapsed": collapsed_summary(rec1), "loci": loci, "rows": rows}

    def frac(rows, pred):
        return round(sum(map(pred, rows)) / max(len(rows), 1), 3)
    S0, SA, SB = (res[s]["rows"] for s in ("S0", "SA", "SB"))
    inf_A = [r for r in SA if r["informative"]]
    inf_B = [r for r in SB if r["informative"] and r["copy"] == "B"]
    called = [r for r in SA + SB if r["call"] in ("ALLELE", "EXTRA_COPY")]
    cA, cB = res["SA"]["collapsed"], res["SB"]["collapsed"]
    v1 = sorted(us + p for s, e in rex for p in var if s <= p < e)
    m = {
        "P0_S0_nonshared": frac(S0, lambda r: r["call"] != "SHARED"),
        "P1_same_columns": cA["variable_columns"] == cB["variable_columns"],
        "P1_columns_are_exonic_V": cA["variable_columns"] == v1,
        "P1_alt_fraction_diff": round(abs(cA["alt_fraction"] - cB["alt_fraction"]), 3),
        "P2_SA_informative_allele": frac(inf_A, lambda r: r["call"] == "ALLELE"),
        "P2_SA_informative_extra": frac(inf_A, lambda r: r["call"] == "EXTRA_COPY"),
        "P2_SB_informativeB_extra": frac(inf_B, lambda r: r["call"] == "EXTRA_COPY"),
        "P2_SB_all_allele": frac(SB, lambda r: r["call"] == "ALLELE"),
        "P3_hap_correct": frac(called, lambda r: r["home"] == r["hap"]),
        "n": {"S0": len(S0), "SA_informative": len(inf_A), "SB_informative_B": len(inf_B), "called": len(called)},
    }
    m["PASS"] = {"P0": m["P0_S0_nonshared"] <= 0.02,
                 "P1": m["P1_same_columns"] and m["P1_alt_fraction_diff"] <= 0.10,
                 "P2": m["P2_SA_informative_allele"] >= 0.95 and m["P2_SA_informative_extra"] <= 0.02
                       and m["P2_SB_informativeB_extra"] >= 0.95 and m["P2_SB_all_allele"] <= 0.05,
                 "P3": m["P3_hap_correct"] >= 0.98}
    res["metrics"] = m
    for s in ("S0", "SA", "SB"):
        res[s]["calls"] = {k: sum(r["call"] == k for r in res[s]["rows"]) for k in ("SHARED", "ALLELE", "EXTRA_COPY", "UNMAPPED")}
        del res[s]["rows"]
    with open(os.path.join(a.out, "results.json"), "w") as fh:
        json.dump(res, fh, indent=1)
    print(json.dumps({k: res[k] for k in ("verify", "metrics")}, indent=1))
    for s in ("S0", "SA", "SB"):
        c = res[s]["collapsed"]
        print(s, "collapsed:", len(c["variable_columns"]), "variable cols, alt fraction", c["alt_fraction"],
              "| phased:", res[s]["calls"], "| loci:", [(l["hap"], l["start"], l["end"], l["partner"]) for l in res[s]["loci"]])


if __name__ == "__main__":
    main()
