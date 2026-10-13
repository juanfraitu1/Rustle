#!/usr/bin/env python3
"""Arm B of docs/PREREG_yag_o2_2026-10-12.md: O2's read-star certificate on the full read and on each half.

    split_half.py RUN_PREFIX READ_LENGTHS.tsv OUT.json
RUN_PREFIX = a `copy_assign --dump-star` run (needs .star_reads.tsv and .assignments.tsv); READ_LENGTHS.tsv =
`read_name<TAB>length`. Certificate (PSV columns only, `read_star_junctions false`): for candidates a, b the columns
where both have a base, the bases differ and the read has a base; LLR(a, b) = lr (k - m) with k / m = columns matching
a / b; best = maximin LLR (ties: PSV score, then lower copy index); Tied if some rival shares no such column with the
best, Assigned if every rival's p = P(Binomial(n, e/3) >= k) < alpha / #rivals and its LLR > 0, else Ambiguous.

Population mode (docs/PREREG_popfix_2026-10-12.md, ruling R27): env POPULATION=psv_population.tsv re-runs the same
certificate on the columns that are FIXED across the population, PER PAIR as Rust `fixed_only_certificate` does: a
column counts for the pair (a, b) iff both candidates have a genome position there and both (chrom, pos) are in the
site set S (chrom of a candidate from <RUN_PREFIX>.family_join.tsv); a candidate's PSV-score tie key counts the columns
where its own position is in S. S = the FIXED sites, or every site of the table with any status when POPULATION_ANY=1
(arm c'). Contested rows only. Under POPULATION (FIXED mode) the full-read masked certificate must also reproduce
<RUN_PREFIX>.population_tier.tsv (`status_fixed`, and `copy_fixed` when assigned): gate P0-fixed, >= 0.99.
"""
import collections, json, math, os, sys

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


def parse_positions(field, n_cols, nc):
    """`positions` field of a star dump -> per column a list of nc genome positions (None = no position)."""
    if not field:
        return [[None] * nc for _ in range(n_cols)]
    out = []
    for entry in field.split(","):
        parts = entry.split("/")
        if len(parts) != nc:
            raise ValueError(f"positions entry {entry!r} has {len(parts)} candidates, expected {nc}")
        out.append([None if x == "." else int(x) for x in parts])
    if len(out) != n_cols:
        raise ValueError(f"positions has {len(out)} columns, expected {n_cols}")
    return out


def load_sites(path, any_status=False):
    """The site set S of a psv_population.tsv: {(chrom, pos)} of the FIXED rows, or of every row with any_status."""
    sites = set()
    with open(path) as fh:
        head = fh.readline().rstrip("\n").split("\t")
        ic, ip, ist = head.index("chrom"), head.index("pos"), head.index("status")
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if any_status or f[ist] == "FIXED":
                sites.add((f[ic], int(f[ip])))
    return sites


def read_chroms(prefix):
    """{(family_id, copy_index): chrom} from <prefix>.family_join.tsv."""
    path = prefix + ".family_join.tsv"
    if not os.path.exists(path):
        raise SystemExit(f"population mode needs {path} (chrom of each candidate copy)")
    return {(r["family_id"], int(r["copy_index"])): r["chrom"] for r in tsv(path)}


def mask_columns(cols, gpos, chroms, sites):
    """Add to each column (readpos, obs, alleles) the per-candidate flags `(chrom_k, gpos[j][k]) in sites`.
    gpos = None (no positions) -> no candidate is in S anywhere (R20: no fixed column, the certificate ties)."""
    out = []
    for j, (pos, obs, al) in enumerate(cols):
        gp = gpos[j] if gpos is not None else [None] * len(chroms)
        ok = [g is not None and ch is not None and (ch, g) in sites for ch, g in zip(chroms, gp)]
        out.append((pos, obs, al, ok))
    return out


def certify(cand, cols, keep=lambda pos: True):
    """(status, copy) by the read-star certificate on the kept columns. Columns are (readpos, obs, alleles) or, from
    mask_columns, (readpos, obs, alleles, in_s) with in_s[k] = candidate k's position is in S: a column then counts for
    the pair (a, b) only if in_s[a] and in_s[b], and for a's PSV score only if in_s[a] (R27)."""
    cols = [c for c in cols if keep(c[0]) and c[1] not in ".-"]
    masked = bool(cols) and len(cols[0]) > 3
    nc = len(cand)
    if nc < 2:
        return "tied", None
    pair = {}
    for a in range(nc):
        for b in range(nc):
            if a == b:
                continue
            n = k = m = 0
            for c in cols:
                o, al = c[1], c[2]
                x, y = al[a], al[b]
                if x in ".-" or y in ".-" or x == y:
                    continue
                if masked and not (c[3][a] and c[3][b]):
                    continue
                n += 1
                k += o == x
                m += o == y
            pair[(a, b)] = (n, k, m)
    def psv_score(a):
        return sum((1 if c[1] == c[2][a] else -1) for c in cols if c[2][a] not in ".-" and (not masked or c[3][a]))
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


def load_population(prefix, header):
    """(sites, any_mode, chroms, tiers) for population mode, None when POPULATION is unset."""
    path = os.environ.get("POPULATION")
    any_mode = os.environ.get("POPULATION_ANY") == "1"
    if not path:
        if any_mode:
            print("split_half: POPULATION_ANY=1 is ignored without POPULATION", file=sys.stderr)
        return None
    if "positions" not in header:
        raise SystemExit("population mode needs a star dump with a `positions` column (copy_assign --dump-star, "
                         "built after the popfix change)")
    sites = load_sites(path, any_status=any_mode)
    chroms = read_chroms(prefix)
    tier_path = prefix + ".population_tier.tsv"
    tiers = {(r["read_name"], r["family_id"]): r for r in tsv(tier_path)} if os.path.exists(tier_path) else None
    return sites, any_mode, chroms, tiers


def main(prefix, lengths, out):
    L = {}
    if lengths != "-":
        for line in open(lengths):
            n, l = line.split()
            L[n] = int(l)
    asg = {(r["read_name"], r["family_id"]): r for r in tsv(prefix + ".assignments.tsv")}
    c = collections.Counter()
    per_fam = collections.defaultdict(collections.Counter)
    only = set(json.load(open(os.environ["YAG_FAMILIES"]))["families"]) if os.environ.get("YAG_FAMILIES") else None
    with open(prefix + ".star_reads.tsv") as fh:
        star_header = fh.readline().rstrip("\n").split("\t")
    pop = load_population(prefix, star_header)
    if pop:
        sites, any_mode, chroms, tiers = pop
        p0_enabled = not any_mode and tiers is not None  # the tier file is the FIXED-only certificate's
        p0_seen, p0_mismatch = set(), []
        # tier rows the gate is about: YAG families, a tier (status_fixed != NA)
        p0_rows = {k: v for k, v in (tiers or {}).items()
                   if v["status_fixed"] != "NA" and (only is None or k[1] in only)}
    for r in tsv(prefix + ".star_reads.tsv"):
        if only is not None and r["family_id"] not in only:  # YAG families only (user, 2026-10-12)
            continue
        a = asg.get((r["read_name"], r["family_id"]))
        cand = [int(x) for x in r["candidates"].split(",")] if r["candidates"] else []
        cols = parse_cols(r.get("columns", ""))
        fam = r["family_id"]
        if pop:
            key = (r["read_name"], fam)
            contested = a is not None and a.get("contested") == "1"
            tier = p0_rows.get(key) if p0_enabled else None
            if not contested and tier is None:  # arms (c), (c'), (e) read contested rows only (Amendment 1)
                c["uncontested_skipped"] += 1
                continue
            cols = mask_columns(cols, parse_positions(r.get("positions", ""), len(cols), len(cand)),
                                [chroms.get((fam, ci)) for ci in cand], sites)
            c["cand_chrom_unknown"] += sum(chroms.get((fam, ci)) is None for ci in cand)
        st, cp = certify(cand, cols)
        if pop:
            if tier is not None:
                p0_seen.add(key)
                partners = r.get("n_candidates", "") not in ("", str(len(cand)))  # Python knows no §6ft partners
                same = st == tier["status_fixed"] and (st != "assigned" or str(cp) == tier["copy_fixed"])
                c["p0fixed_same"] += same
                c["p0fixed_with_partners"] += partners
                per_fam[fam]["p0fixed_n"] += 1
                per_fam[fam]["p0fixed_same"] += same
                if not same:
                    c["p0fixed_mismatch_with_partners"] += partners
                    p0_mismatch.append({"read_name": r["read_name"], "family_id": fam,
                                        "tier": [tier["status_fixed"], tier["copy_fixed"]],
                                        "python": [st, str(cp)], "n_dec_fixed": tier["n_dec_fixed"],
                                        "candidates_with_partners": partners})
            if not contested:
                continue
            c["pop_full|" + st] += 1
            per_fam[fam]["full|" + st] += 1
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
    if pop:
        # the gate P0-fixed (Amendment 1): the masked full-read certificate vs the tier file's FIXED-only call
        res["population"] = {"mode": "any" if any_mode else "fixed", "path": os.environ["POPULATION"],
                             "n_sites": len(sites),
                             "note": "B_P0_* compares the MASKED certificate to the run's all-column status: informational "
                                     "here; the gate is P0_fixed_*. B_P1/B_P2 and the halves are on masked columns, "
                                     "contested rows only.",
                             "full_read_status": {k[len("pop_full|"):]: v for k, v in sorted(c.items()) if k.startswith("pop_full|")}}
        if p0_enabled:
            missing = [k for k in p0_rows if k not in p0_seen]
            n_rows = len(p0_rows)
            same = c["p0fixed_same"]
            miss_rows = [{"read_name": k[0], "family_id": k[1], "reason": "no_star_row",
                          "tier": [p0_rows[k]["status_fixed"], p0_rows[k]["copy_fixed"]]} for k in missing]
            res["P0_fixed_reproduction"] = round(same / n_rows, 4) if n_rows else None
            res["P0_fixed_pass"] = res["P0_fixed_reproduction"] is not None and res["P0_fixed_reproduction"] >= 0.99
            res["population"]["P0_fixed"] = {
                "n_tier_rows": n_rows, "n_compared": len(p0_seen), "n_no_star_row": len(missing), "n_same": same,
                "n_rows_with_partners": c["p0fixed_with_partners"],
                "n_mismatch_with_partners": c["p0fixed_mismatch_with_partners"],
                "first_mismatches": (p0_mismatch + miss_rows)[:5]}
        else:
            res["P0_fixed_reproduction"] = None
            res["P0_fixed_pass"] = None
            res["population"]["P0_fixed"] = ("not applicable: POPULATION_ANY is not the FIXED-only certificate" if any_mode
                                             else "no population_tier.tsv beside the run")
    res["per_family"] = {f: dict(v) for f, v in sorted(per_fam.items(), key=lambda x: -x[1]["n"])[:15]}
    json.dump(res, open(out, "w"), indent=1)
    print(json.dumps({k: v for k, v in res.items() if k != "per_family"}, indent=1))


if __name__ == "__main__":
    main(*sys.argv[1:4])
