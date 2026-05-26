#!/usr/bin/env python3
"""Diff two parity_decisions JSONL logs (rustle + stringtie).

Usage:
    diff.py rustle.jsonl stringtie.jsonl [--step junction_accept] [--range 22149137-22155355]

Joins events by (step, start, end, strand). Reports:
- Events present in both, payloads identical
- Events present in both, payloads differ → diff fields
- Events only in rustle
- Events only in stringtie

Special handling for bundlenode_list and graphnode_list: parses the compact "s-e:cov"
node list and reports structural differences (nodes only in rustle, only in ST, or with
different cov). Both steps share the same n_nodes/nodes payload schema.
"""
import json, sys, argparse, collections
from typing import Dict, List, Tuple, Any, Optional

EventKey = Tuple[str, int, int, str]   # (step, start, end, strand)


def load(path: str, step_filter: str = None, range_filter: Tuple[int, int] = None):
    """Returns dict keyed by EventKey → list of payloads (multi-emit allowed)."""
    by_key: Dict[EventKey, List[dict]] = collections.defaultdict(list)
    n_total = 0
    with open(path) as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            try:
                ev = json.loads(line)
            except json.JSONDecodeError:
                continue
            step = ev.get("step", "")
            if step.startswith("_"):
                continue
            if step_filter and step != step_filter:
                continue
            s, e = int(ev.get("start", 0)), int(ev.get("end", 0))
            if range_filter:
                rs, re_ = range_filter
                if not (s == 0 and e == 0):
                    if e < rs or s > re_:
                        continue
            strand = ev.get("strand", "?")
            key = (step, s, e, strand)
            by_key[key].append(ev.get("payload", {}))
            n_total += 1
    return by_key, n_total


def diff_payloads(a: dict, b: dict) -> List[str]:
    """Return list of "key: a_val ≠ b_val" strings for differing fields."""
    out = []
    keys = sorted(set(a.keys()) | set(b.keys()))
    for k in keys:
        av, bv = a.get(k), b.get(k)
        if av != bv:
            out.append(f"{k}: rustle={av} stringtie={bv}")
    return out


# ---------------------------------------------------------------------------
# bundlenode_list specialised diff
# ---------------------------------------------------------------------------

def parse_nodes(nodes_str: str) -> List[Tuple[int, int, float]]:
    """Parse "s1-e1:cov1,s2-e2:cov2,..." into [(start, end, cov), ...]."""
    if not nodes_str:
        return []
    result = []
    for part in nodes_str.split(","):
        part = part.strip()
        if not part:
            continue
        try:
            se, cov = part.rsplit(":", 1)
            s, e = se.split("-")
            result.append((int(s), int(e), float(cov)))
        except ValueError:
            continue
    return result


def diff_bundlenode_payloads(rp: dict, sp: dict, cov_tol: float = 0.5) -> List[str]:
    """Structural diff of two bundlenode_list payloads.

    Matches nodes by (start, end). Reports:
    - n_nodes mismatch summary
    - nodes only in rustle (by start-end)
    - nodes only in stringtie (by start-end)
    - nodes in both with |cov_r - cov_s| > cov_tol
    """
    out = []
    rn = rp.get("n_nodes", "?")
    sn = sp.get("n_nodes", "?")
    if rn != sn:
        out.append(f"n_nodes: rustle={rn} stringtie={sn}")

    r_nodes = parse_nodes(rp.get("nodes", ""))
    s_nodes = parse_nodes(sp.get("nodes", ""))

    r_by_se = {(n[0], n[1]): n[2] for n in r_nodes}
    s_by_se = {(n[0], n[1]): n[2] for n in s_nodes}

    only_r = sorted(set(r_by_se) - set(s_by_se))
    only_s = sorted(set(s_by_se) - set(r_by_se))
    both   = sorted(set(r_by_se) & set(s_by_se))

    cov_diffs = [(se, r_by_se[se], s_by_se[se]) for se in both
                 if abs(r_by_se[se] - s_by_se[se]) > cov_tol]

    if only_r:
        out.append(f"  nodes only in rustle ({len(only_r)}):")
        for se in only_r[:10]:
            out.append(f"    {se[0]}-{se[1]} cov={r_by_se[se]:.1f}")
        if len(only_r) > 10:
            out.append(f"    ... and {len(only_r)-10} more")
    if only_s:
        out.append(f"  nodes only in stringtie ({len(only_s)}):")
        for se in only_s[:10]:
            out.append(f"    {se[0]}-{se[1]} cov={s_by_se[se]:.1f}")
        if len(only_s) > 10:
            out.append(f"    ... and {len(only_s)-10} more")
    if cov_diffs:
        out.append(f"  nodes with cov diff > {cov_tol} ({len(cov_diffs)}):")
        for (se, rc, sc) in cov_diffs[:5]:
            out.append(f"    {se[0]}-{se[1]}: rustle={rc:.1f} stringtie={sc:.1f}")
        if len(cov_diffs) > 5:
            out.append(f"    ... and {len(cov_diffs)-5} more")

    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("rustle_log")
    ap.add_argument("stringtie_log")
    ap.add_argument("--step", default=None, help="filter to one step name")
    ap.add_argument("--range", default=None, help="filter to coord range LO-HI")
    ap.add_argument("--show-equal", action="store_true", help="also print events where payloads agree")
    ap.add_argument("--max", type=int, default=20, help="max events of each kind to print")
    ap.add_argument("--cov-tol", type=float, default=0.5,
                    help="coverage tolerance for bundlenode_list node matching (default 0.5)")
    args = ap.parse_args()

    rng = None
    if args.range:
        lo, hi = args.range.split("-")
        rng = (int(lo), int(hi))

    rustle, n_r = load(args.rustle_log, args.step, rng)
    string, n_s = load(args.stringtie_log, args.step, rng)
    print(f"loaded rustle:    {len(rustle):5d} unique keys ({n_r} events)")
    print(f"loaded stringtie: {len(string):5d} unique keys ({n_s} events)")

    common = set(rustle) & set(string)
    only_r = set(rustle) - set(string)
    only_s = set(string) - set(rustle)
    print(f"\ncommon keys: {len(common)}")
    print(f"only rustle: {len(only_r)}")
    print(f"only stringtie: {len(only_s)}")

    # Steps that use the compact n_nodes/nodes payload schema
    _NODE_LIST_STEPS = {"bundlenode_list", "graphnode_list"}

    def _is_node_list(step: str) -> bool:
        return step in _NODE_LIST_STEPS

    # Within common, count payload divergences
    eq, ne = 0, 0
    diffs = []
    for k in sorted(common):
        rps = rustle[k]
        sps = string[k]
        rp, sp = rps[0], sps[0]
        if _is_node_list(k[0]):
            structural = diff_bundlenode_payloads(rp, sp, cov_tol=args.cov_tol)
            if not structural:
                eq += 1
            else:
                ne += 1
                diffs.append((k, rp, sp, structural))
        else:
            if rp == sp:
                eq += 1
            else:
                ne += 1
                diffs.append((k, rp, sp, None))

    print(f"common w/ equal payload:    {eq}")
    print(f"common w/ DIFFERING payload: {ne}")

    if args.show_equal and eq > 0:
        print(f"\n--- {min(args.max, eq)} equal-payload events ---")
        n = 0
        for k in sorted(common):
            rp, sp = rustle[k][0], string[k][0]
            is_eq = (_is_node_list(k[0]) and not diff_bundlenode_payloads(rp, sp, args.cov_tol)) \
                    or (not _is_node_list(k[0]) and rp == sp)
            if is_eq:
                if _is_node_list(k[0]):
                    print(f"  {k} :: n_nodes={rp.get('n_nodes')}")
                else:
                    print(f"  {k} :: {rp}")
                n += 1
                if n >= args.max:
                    break

    if ne > 0:
        print(f"\n--- {min(args.max, ne)} DIFFERING-payload events ---")
        for entry in diffs[:args.max]:
            k, rp, sp = entry[0], entry[1], entry[2]
            structural = entry[3]
            print(f"  {k}")
            if structural:
                for line in structural:
                    print(f"    {line}")
            else:
                for d in diff_payloads(rp, sp):
                    print(f"    {d}")

    if only_r:
        print(f"\n--- {min(args.max, len(only_r))} ONLY in rustle ---")
        for k in sorted(only_r)[:args.max]:
            p = rustle[k][0]
            if _is_node_list(k[0]):
                print(f"  {k} :: n_nodes={p.get('n_nodes')}")
            else:
                print(f"  {k} :: {p}")

    if only_s:
        print(f"\n--- {min(args.max, len(only_s))} ONLY in stringtie ---")
        for k in sorted(only_s)[:args.max]:
            p = string[k][0]
            if _is_node_list(k[0]):
                print(f"  {k} :: n_nodes={p.get('n_nodes')}")
            else:
                print(f"  {k} :: {p}")


if __name__ == "__main__":
    main()
