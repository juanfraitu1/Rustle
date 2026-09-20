#!/usr/bin/env python3
"""Compare Rustle vs StringTie `pred_kill` parity events.

Loads two JSONL parity logs (one per tool), tallies kills per reason for each
tool, and — for the predcluster pairwise reasons we now reproduce
(`retained_intron`, `included_drop`) — reports how many of StringTie's kills
Rustle also reproduces.

Match key: (start, end, strand, nexons). StringTie's pred_kill payload does NOT
carry the intron chain, so we match on the kill record's genomic span + strand +
exon count. This is a coordinate match, not an intron-chain match: two distinct
predictions sharing span+strand+nexons would collide, but in practice predcluster
kills at the same (start,end,strand,nexons) are the same prediction. Documented as
an approximation.

Usage:
    pred_kill_diff.py RUSTLE_PK.jsonl STRINGTIE_PK.jsonl

Defaults to /tmp/ru_pk.jsonl and /tmp/st_pk.jsonl if no args given.
"""
import json
import sys
from collections import Counter, defaultdict

# Reasons emitted by the predcluster pairwise stage that Task 3 reproduces.
PAIRWISE_REASONS = ("retained_intron", "included_drop")


def load(path):
    """Return list of (reason, key) for every pred_kill record in `path`."""
    out = []
    with open(path) as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            rec = json.loads(line)
            if rec.get("step") != "pred_kill":
                continue
            payload = rec.get("payload", {})
            reason = payload.get("reason", "?")
            key = (
                int(rec["start"]),
                int(rec["end"]),
                rec.get("strand", "."),
                int(payload.get("nexons", -1)),
            )
            out.append((reason, key))
    return out


def main():
    ru_path = sys.argv[1] if len(sys.argv) > 1 else "/tmp/ru_pk.jsonl"
    st_path = sys.argv[2] if len(sys.argv) > 2 else "/tmp/st_pk.jsonl"

    ru = load(ru_path)
    st = load(st_path)

    ru_by_reason = Counter(r for r, _ in ru)
    st_by_reason = Counter(r for r, _ in st)

    print("=== pred_kill counts per reason ===")
    all_reasons = sorted(set(ru_by_reason) | set(st_by_reason))
    print(f"{'reason':<24}{'rustle':>8}{'stringtie':>11}")
    for reason in all_reasons:
        print(f"{reason:<24}{ru_by_reason.get(reason, 0):>8}{st_by_reason.get(reason, 0):>11}")
    print(f"{'TOTAL':<24}{len(ru):>8}{len(st):>11}")

    # Rustle keys grouped by reason, for matching.
    ru_keys_by_reason = defaultdict(set)
    for reason, key in ru:
        ru_keys_by_reason[reason].add(key)
    # Also a reason-agnostic Rustle key set (Rustle may kill via a different reason
    # label than ST for the same prediction; report both strict and any-reason).
    ru_keys_any = set(k for _, k in ru)

    print()
    print("=== ST retained_intron + included_drop kills reproduced by Rustle ===")
    print("(match on start,end,strand,nexons; ST payload lacks introns)")
    total_st_pairwise = 0
    total_repro_same_reason = 0
    total_repro_any_reason = 0
    for reason in PAIRWISE_REASONS:
        st_keys = [k for r, k in st if r == reason]
        total_st_pairwise += len(st_keys)
        same = sum(1 for k in st_keys if k in ru_keys_by_reason.get(reason, set()))
        anyr = sum(1 for k in st_keys if k in ru_keys_any)
        total_repro_same_reason += same
        total_repro_any_reason += anyr
        print(
            f"  {reason:<18} ST={len(st_keys):>4}  "
            f"reproduced(same reason)={same:>4}  "
            f"reproduced(any reason)={anyr:>4}"
        )
    print(
        f"  {'TOTAL':<18} ST={total_st_pairwise:>4}  "
        f"reproduced(same reason)={total_repro_same_reason:>4}  "
        f"reproduced(any reason)={total_repro_any_reason:>4}"
    )


if __name__ == "__main__":
    main()
