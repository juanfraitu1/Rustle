#!/usr/bin/env python3
"""r1642 (docs/PREREG_exon_space_placement_2026-11-12.md): place reads by an unspliced alignment to the catalog copies'
spliced sequences, score against strict_truth dumps (STRICT_DUMP_OUTSIDE=1). Counts molecules (read names).
Usage: exon_place.py PAF COPIES.tsv PANEL_DUMP HELDOUT_DUMP [MOLECULES.tsv: primary splice shape, column 11]"""
import collections, csv, os, sys

LOCUS = os.environ.get("EXON_PLACE_RIGHT") == "locus"  # r1644 post hoc: a pick overlapping a truth copy is right

def load_paf(path):
    """read -> {(family, copy): (AS, query cover)} from the best-AS record per copy."""
    best = collections.defaultdict(dict)
    if path.endswith(".picks.tsv"):  # r1647: a tool's own picks (read, family, copy); absent reads abstain
        for line in open(path):
            q, fam, idx = line.rstrip("\n").split("\t"); best[q][(fam, idx)] = (1, 1.0)
        return best
    with open(path) as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            a = next(int(t[5:]) for t in f[12:] if t.startswith("AS:i:"))
            fam, idx = f[5].split("|")[:2]
            cov = (int(f[3]) - int(f[2])) / int(f[1])
            k = (fam, idx); q = f[0][2:] if f[0].startswith("R|") else f[0]
            if k not in best[q] or a > best[q][k][0]: best[q][k] = (a, cov)
    return best

def place(recs, span, margin=1, cover=0.90):
    """(copy, margin to the best other LOCUS) under the registered rule, or (None, why)."""
    if not recs: return None, "no_record"
    top = max(recs, key=lambda k: recs[k][0])
    if top[0] == "OTHER": return None, "other"  # r1647: best hit is a transcript outside the catalog
    def same(a, b): return span[a][0] == span[b][0] and span[a][1] < span[b][2] and span[b][1] < span[a][2]
    rival = max((recs[k][0] for k in recs if not same(k, top)), default=None)
    m = 10 ** 6 if rival is None else recs[top][0] - rival
    if m < margin: return None, "tie"
    if recs[top][1] < cover: return None, "cover"
    return top, m

def main(paf, copies, panel, heldout, molecules=None):
    span = {}
    with open(copies) as fh:
        for r in csv.DictReader(fh, delimiter="\t"): span[(r["family_id"], r["copy_idx"])] = (r["chrom"], int(r["start"]), int(r["end"]))
    best = load_paf(paf); span[("OTHER", "0")] = ("OTHER", 0, 1)
    shape = {}
    if molecules:
        with open(molecules) as fh:
            for line in fh:
                if line[0] != "#": f = line.split("\t", 12); shape[f[0]] = f[10] if len(f) > 10 else "NA"  # as_table column 11 = primary spliced
    mol = {}  # read -> dict(set, fam, truth, al, o2)
    for name, path in (("panel", panel), ("heldout", heldout)):
        with open(path) as fh:
            for line in fh:
                q, fam, st, v, al, _, t = line.rstrip("\n").split("\t")
                m = mol.setdefault(q, {"set": name, "fam": fam, "al": al, "o2": "abstain",
                                       "truth": None if t == "OUTSIDE" else {tuple(x.split(":")) for x in t.split(";")}})
                if st == "assigned": m["o2"] = v if m["o2"] != "wrong" else "wrong"
    if os.environ.get("EXON_PLACE_PRIMARY"):  # r1648: primary right = its RECORD (name chrom start end) overlaps a truth copy
        rec = {}
        for line in open(os.environ["EXON_PLACE_PRIMARY"]):
            f = line.split("\t")
            if f[0] in mol: rec[f[0]] = (f[1], int(f[2]), int(f[3]))
        for q, m in mol.items():
            if m["truth"] and q in rec: m["al"] = "al_right" if any(span[t][0] == rec[q][0] and rec[q][1] < span[t][2] and span[t][1] < rec[q][2] for t in m["truth"]) else "al_wrong"
    only = os.environ.get("EXON_PLACE_SHAPE")  # r1644: restrict every table to one primary shape (0 unspliced, 1 spliced)
    if only is not None: mol = {q: m for q, m in mol.items() if shape.get(q) == only}
    dec = {q: m for q, m in mol.items() if m["truth"]}; neg = {q: m for q, m in mol.items() if not m["truth"]}
    big = collections.Counter(m["fam"] for m in dec.values()).most_common(1)[0][0]

    dump = open(os.environ["EXON_PLACE_DUMP"], "w") if os.environ.get("EXON_PLACE_DUMP") else None  # per molecule, pooled table only
    def table(label, qs, margin=1, cover=0.90):
        c = collections.Counter()
        for q in qs:
            m = dec[q]; pick, why = place(best.get(q, {}), span, margin, cover)
            if dump and label == "REGISTERED pooled":
                r = pick is not None and (pick in m["truth"] or any(span[pick][0] == span[t][0] and span[pick][1] < span[t][2] and span[t][1] < span[pick][2] for t in m["truth"]))
                dump.write(f"{q}\t{m['fam']}\t{'assigned' if pick else why}\t{'NA' if pick is None else int(r)}\t{int(m['al'] == 'al_right')}\t{m['o2']}\n")
            if pick is None: c["abstain_" + str(why)] += 1; continue
            arm = pick in m["truth"]; al = m["al"] == "al_right"
            if LOCUS: arm = arm or any(span[pick][0] == span[t][0] and span[pick][1] < span[t][2] and span[t][1] < span[pick][2] for t in m["truth"])  # as the primary is scored
            c["assigned"] += 1; c["right" if arm else "wrong"] += 1; c["prim_right" if al else "prim_wrong"] += 1
            c["arm+prim-" if arm and not al else "arm-prim+" if al and not arm else "agree"] += 1
            if m["o2"] != "abstain": c["both_o2"] += 1; c["o2_right_on_both"] += m["o2"] == "right"; c["arm_right_on_both"] += arm
        n = len(qs)
        print(f"{label}\tn={n}\tassigned={c['assigned']} ({c['assigned'] / max(n, 1):.3f})\tarm {c['right']}/{c['wrong']}\tprimary {c['prim_right']}/{c['prim_wrong']}"
              f"\tarm+prim-={c['arm+prim-']}\tarm-prim+={c['arm-prim+']}\tO2-both n={c['both_o2']} arm {c['arm_right_on_both']} o2 {c['o2_right_on_both']}"
              f"\tabstain tie/cover/norec/other={c['abstain_tie']}/{c['abstain_cover']}/{c['abstain_no_record']}/{c['abstain_other']}")
        return c

    print(f"molecules: decidable {len(dec)} (panel {sum(m['set'] == 'panel' for m in dec.values())}), negative {len(neg)}; largest family {big}")
    table("REGISTERED pooled", list(dec))
    table("REGISTERED heldout", [q for q in dec if dec[q]["set"] == "heldout"])
    table("REGISTERED panel", [q for q in dec if dec[q]["set"] == "panel"])
    table(f"REGISTERED minus {big}", [q for q in dec if dec[q]["fam"] != big])
    na = sum(place(best.get(q, {}), span)[0] is not None for q in neg)
    print(f"NEGATIVE assigned {na}/{len(neg)} ({na / max(len(neg), 1):.3f})")
    print("-- readout, no bar")
    for mg in (5, 10): table(f"margin>={mg}", list(dec), mg)
    table("cover gate off", list(dec), 1, 0.0)
    print("truth copy has no record:", sum(not (dec[q]["truth"] & set(best.get(q, {}))) for q in dec), "of", len(dec))
    for s in sorted({shape.get(q, "NA") for q in dec}): table(f"primary spliced={s}", [q for q in dec if shape.get(q, "NA") == s])
    fams = collections.Counter(m["fam"] for m in dec.values()).most_common(8)
    for f, _ in fams: table(f"family {f}", [q for q in dec if dec[q]["fam"] == f])

def _selfcheck():
    span = {("F", "0"): ("c", 0, 10), ("F", "1"): ("c", 5, 15), ("F", "2"): ("c", 100, 110)}
    assert place({("F", "0"): (100, 1.0), ("F", "1"): (100, 1.0), ("F", "2"): (99, 1.0)}, span)[0] == ("F", "0")  # overlapping copies = one locus
    assert place({("F", "0"): (100, 1.0), ("F", "2"): (100, 1.0)}, span) == (None, "tie")
    assert place({("F", "0"): (100, 0.5), ("F", "2"): (90, 1.0)}, span) == (None, "cover")
    assert place({}, span) == (None, "no_record")
    span[("OTHER", "0")] = ("OTHER", 0, 1)
    assert place({("OTHER", "0"): (100, 1.0), ("F", "0"): (99, 1.0)}, span) == (None, "other")
    assert place({("F", "0"): (100, 1.0), ("OTHER", "0"): (100, 1.0)}, span)[0] is None  # a tie with OTHER abstains

if __name__ == "__main__":
    _selfcheck()
    main(*sys.argv[1:6])
