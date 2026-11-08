#!/usr/bin/env python3
"""Register 1582 scorer (docs/PREREG_crossind_copyphasing_2026-11-08.md): copy_phasing hits at the crossind strong loci.
usage: copyphasing_lookup.py STRONG.tsv OR_PREFIX OR_COPIES KB_PREFIX KB_COPIES"""
import csv, sys

def copies(path):
    by_key, rows = {}, []
    for r in csv.DictReader(open(path), delimiter="\t"):
        c = (r["chrom"], int(r["start"]), int(r["end"]), r["family_id"])
        rows.append(c)
        for k in (r["tid"], f'{r["family_id"]}:{r["copy_idx"]}', (r["family_id"], r["copy_idx"])):
            by_key[k] = c
    return by_key, rows

def hits(prefix, by_key):
    out = []  # (copy, kind)
    def res(fam, home):
        for k in (home, f"{fam}:{home}", (fam, home)):
            if k in by_key:
                return by_key[k]
        sys.exit(f"unresolved home {home!r} in family {fam} ({prefix})")
    for r in csv.DictReader(open(f"{prefix}.candidates.tsv"), delimiter="\t"):
        if r["flagged"] == "1":
            out.append((res(r["family"], r["nearest_locus"]), "certified"))
    for r in csv.DictReader(open(f"{prefix}.putative.tsv"), delimiter="\t"):
        out.append((res(r["family"], r["copy"]), "putative"))
    return out

ov = lambda c, w: c[0] == w[0] and c[1] < w[2] and w[1] < c[2]

def main(strong, orp, orc, kbp, kbc):
    loci = [(r[0], int(r[1]), int(r[2])) for r in csv.reader(open(strong), delimiter="\t")]
    S = {}
    for name, p, cp in (("OR", orp, orc), ("KB", kbp, kbc)):
        by_key, rows = copies(cp)
        h = hits(p, by_key)
        hit_copies = {c for c, _ in h}
        S[name] = (rows, h, hit_copies)
        print(f"{name}: copies {len(rows)}, hit copies {len(hit_copies)}, background {len(hit_copies)/len(rows):.4f} "
              f"(certified {sum(k == 'certified' for _, k in h)}, putative {sum(k == 'putative' for _, k in h)})")
    print("locus\tOR_copies\tOR_fams\tOR_cert\tOR_put\tKB_copies\tKB_fams\tKB_cert\tKB_put")
    tally = dict(or_test=0, or_hit=0, or_cert=0, kb_test=0, kb_hit=0, win_copies=0, win_hit=0)
    for w in loci:
        row = [f"{w[0]}:{w[1]}-{w[2]}"]
        for name in ("OR", "KB"):
            rows, h, hc = S[name]
            cs = [c for c in rows if ov(c, w)]
            cert = sum(k == "certified" and ov(c, w) for c, k in h)
            put = sum(k == "putative" and ov(c, w) for c, k in h)
            row += [len(cs), ",".join(sorted({c[3] for c in cs})) or "-", cert, put]
            testable = bool(cs) and w[0] != "NC_073247.2"  # chrX excluded from the verdict (prereg; amended: the prereg named chr2 by mistake)
            if name == "OR" and testable:
                tally["or_test"] += 1; tally["or_hit"] += (cert + put) > 0; tally["or_cert"] += cert > 0
                tally["win_copies"] += len(cs); tally["win_hit"] += sum(c in hc for c in cs)
            if name == "KB" and testable and S["OR"] and any(ov(c, w) for c in S["OR"][0]):
                tally["kb_test"] += 1; tally["kb_hit"] += (cert + put) > 0
        print("\t".join(map(str, row)))
    t = tally
    bg = len(S["OR"][2]) / len(S["OR"][0])
    win = t["win_hit"] / t["win_copies"] if t["win_copies"] else 0.0
    need = 4 if t["or_test"] >= 8 else -(-t["or_test"] // 2)
    b1, b2, b3 = t["or_hit"] >= need and t["or_test"] > 0, t["kb_hit"] <= 1, win >= 2 * bg
    print(f"OR testable {t['or_test']}, OR hit {t['or_hit']} (certified-only {t['or_cert']}), need {need}: bar1 {'PASS' if b1 else 'FAIL'}")
    print(f"KB testable {t['kb_test']}, KB hit {t['kb_hit']}: bar2 {'PASS' if b2 else 'FAIL'}")
    print(f"OR per-copy in windows {t['win_hit']}/{t['win_copies']} = {win:.4f} vs background {bg:.4f}: bar3 {'PASS' if b3 else 'FAIL'}")
    v = "PASS" if b1 and b2 and b3 else "FAMILY COMPLEXITY" if b1 and not b2 else "NOT CORROBORATED" if not b1 else "FAIL (bar 3)"
    print("VERDICT", v)

if __name__ == "__main__":
    main(*sys.argv[1:6])
