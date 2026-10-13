#!/usr/bin/env python3
"""Arm B of docs/PREREG_yag_o2_2026-10-12.md: O2's read-star certificate on the full read and on each half.

    split_half.py RUN_PREFIX READ_LENGTHS.tsv OUT.json
RUN_PREFIX = a `copy_assign --dump-star` run (needs .star_reads.tsv and .assignments.tsv); READ_LENGTHS.tsv =
`read_name<TAB>length`. Certificate (PSV columns only, `read_star_junctions false`): for candidates a, b the columns
where both have a base, the bases differ and the read has a base; LLR(a, b) = lr (k - m) with k / m = columns matching
a / b; best = maximin LLR (ties: PSV score, then lower copy index); Tied if some rival shares no such column with the
best, Assigned if every rival's p = P(Binomial(n, e/3) >= k) < alpha / #rivals and its LLR > 0, else Ambiguous.
"""
import collections, json, math, sys

ERR, ALPHA = 0.003, 1e-3
EPS = ERR / 3
LR = math.log((1 - ERR) / EPS)


def tail(k, n):
    return sum(math.comb(n, j) * EPS ** j * (1 - EPS) ** (n - j) for j in range(k, n + 1))


def parse_cols(field):
    if not field:
        return []
    out = []
    for tok in field.split(","):
        pos, obs, al = tok.split(":", 2)
        out.append((int(pos), obs, al))
    return out


def certify(cand, cols, keep=lambda pos: True):
    """(status, copy) by the read-star certificate on the kept columns."""
    cols = [c for c in cols if keep(c[0]) and c[1] not in ".-"]
    nc = len(cand)
    if nc < 2:
        return "tied", None
    pair = {}
    for a in range(nc):
        for b in range(nc):
            if a == b:
                continue
            n = k = m = 0
            for _, o, al in cols:
                x, y = al[a], al[b]
                if x in ".-" or y in ".-" or x == y:
                    continue
                n += 1
                k += o == x
                m += o == y
            pair[(a, b)] = (n, k, m)
    def psv_score(a):
        return sum((1 if o == al[a] else -1) for _, o, al in cols if al[a] not in ".-")
    def key(a):
        mins = min(LR * (pair[(a, b)][1] - pair[(a, b)][2]) for b in range(nc) if b != a)
        return (mins, psv_score(a), -a)
    bk = max(range(nc), key=key)
    rivals = [b for b in range(nc) if b != bk]
    if any(pair[(bk, b)][0] == 0 for b in rivals):
        return "tied", cand[bk]
    thr = ALPHA / len(rivals)
    p_read = max(tail(pair[(bk, b)][1], pair[(bk, b)][0]) for b in rivals)
    margin = min(LR * (pair[(bk, b)][1] - pair[(bk, b)][2]) for b in rivals)
    return ("assigned" if p_read < thr and margin > 0 else "ambiguous"), cand[bk]


def tsv(path):
    with open(path) as fh:
        head = fh.readline().rstrip("\n").split("\t")
        return [dict(zip(head, line.rstrip("\n").split("\t"))) for line in fh]


def main(prefix, lengths, out):
    L = {}
    if lengths != "-":
        for line in open(lengths):
            n, l = line.split()
            L[n] = int(l)
    asg = {(r["read_name"], r["family_id"]): r for r in tsv(prefix + ".assignments.tsv")}
    c = collections.Counter()
    per_fam = collections.defaultdict(collections.Counter)
    import os
    only = set(json.load(open(os.environ["YAG_FAMILIES"]))["families"]) if os.environ.get("YAG_FAMILIES") else None
    for r in tsv(prefix + ".star_reads.tsv"):
        if only is not None and r["family_id"] not in only:  # YAG families only (user, 2026-10-12)
            continue
        a = asg.get((r["read_name"], r["family_id"]))
        cand = [int(x) for x in r["candidates"].split(",")] if r["candidates"] else []
        cols = parse_cols(r.get("columns", ""))
        st, cp = certify(cand, cols)
        o2_st, o2_cp = r["status"], r["assigned_copy"]
        gated = a is not None and (a.get("origin_rejected") == "1" or a.get("tie_outside_catalog") == "1")
        if not gated:
            c["p0_n"] += 1
            same = st == o2_st and (st != "assigned" or str(cp) == o2_cp)
            c["p0_same"] += same
            if not same:
                c[f"p0_diff|{o2_st}->{st}"] += 1
        if r["read_name"] not in L:
            c["no_length"] += 1
            continue
        half = L[r["read_name"]] / 2
        s5, c5 = certify(cand, cols, lambda p: p < half)
        s3, c3 = certify(cand, cols, lambda p: p >= half)
        fam = r["family_id"]
        per_fam[fam]["n"] += 1
        if s5 == "assigned" or s3 == "assigned":
            c["some_half_assigns"] += 1
            per_fam[fam]["some_half_assigns"] += 1
        if s5 == "assigned" and s3 == "assigned":
            c["both_halves_assign"] += 1
            c["both_halves_agree"] += c5 == c3
            per_fam[fam]["both"] += 1
            per_fam[fam]["agree"] += c5 == c3
        if o2_st == "assigned":
            for s, cc in ((s5, c5), (s3, c3)):
                if s == "assigned":
                    c["o2_and_half_assign"] += 1
                    c["o2_and_half_agree"] += str(cc) == o2_cp
    res = {"counts": dict(c)}
    res["B_P0_reproduction"] = round(c["p0_same"] / c["p0_n"], 4) if c["p0_n"] else None
    res["B_P0_pass"] = res["B_P0_reproduction"] is not None and res["B_P0_reproduction"] >= 0.99
    nb = c["both_halves_assign"]
    res["B_P1_split_half_concordance"] = {"n": nb, "agree": c["both_halves_agree"],
                                          "rate": round(c["both_halves_agree"] / nb, 4) if nb else None,
                                          "verdict": "UNDERPOWERED" if nb < 30 else ("PASS" if c["both_halves_agree"] / nb >= 0.98 else "FAIL")}
    no = c["o2_and_half_assign"]
    res["B_P2_full_vs_half"] = {"n": no, "agree": c["o2_and_half_agree"],
                                "rate": round(c["o2_and_half_agree"] / no, 4) if no else None,
                                "verdict": "UNDERPOWERED" if no < 30 else ("PASS" if c["o2_and_half_agree"] / no >= 0.98 else "FAIL")}
    res["per_family"] = {f: dict(v) for f, v in sorted(per_fam.items(), key=lambda x: -x[1]["n"])[:15]}
    json.dump(res, open(out, "w"), indent=1)
    print(json.dumps({k: v for k, v in res.items() if k != "per_family"}, indent=1))


if __name__ == "__main__":
    main(*sys.argv[1:4])
