#!/usr/bin/env python3
"""Diff transfrag_collapse provenance between Rustle and StringTie.

Both tools emit one `transfrag_collapse` event per kept representative of the
keeptrf containment-collapse, recording which member transfrags were folded in
(intron-chain + abundance) and the consolidated group coverage. This tool joins
those events by (strand, rep_introns) and reports provenance divergence:
  - rep chains present in both / Rustle-only / ST-only
  - on the common set: how many match on group_cov (within 5%) and on n_members
  - sample of the largest group_cov mismatches (Rustle folds in much more/less)

Inputs: /tmp/ru_tc.jsonl  /tmp/st_tc.jsonl
Capture with PARITY_FILTER_STEPS=transfrag_collapse on a single chrom.
"""
import json
import sys


def load(path):
    """(strand, rep_introns) -> dict(group_cov, n_members, start, end)."""
    out = {}
    try:
        fh = open(path)
    except OSError as e:
        print(f"ERROR: cannot open {path}: {e}", file=sys.stderr)
        return out
    for line in fh:
        line = line.strip()
        if not line:
            continue
        try:
            e = json.loads(line)
        except Exception:
            continue
        if e.get("step") != "transfrag_collapse":
            continue
        p = e.get("payload", e)
        key = (str(e.get("strand", "?")), str(p.get("rep_introns", "")))
        rec = {
            "group_cov": float(p.get("group_cov", 0.0)),
            "n_members": int(p.get("n_members", 0)),
            "start": e.get("start"),
            "end": e.get("end"),
        }
        # If a rep chain appears more than once (multiple bundles), keep the
        # one with the largest group_cov so the join is deterministic.
        prev = out.get(key)
        if prev is None or rec["group_cov"] > prev["group_cov"]:
            out[key] = rec
    return out


def main():
    ru = load("/tmp/ru_tc.jsonl")
    st = load("/tmp/st_tc.jsonl")
    print(f"Rustle transfrag_collapse reps: {len(ru)}")
    print(f"ST     transfrag_collapse reps: {len(st)}")

    ru_keys = set(ru)
    st_keys = set(st)
    both = ru_keys & st_keys
    ru_only = ru_keys - st_keys
    st_only = st_keys - ru_keys

    print()
    print(f"rep chains in BOTH:      {len(both)}")
    print(f"rep chains Rustle-only:  {len(ru_only)}")
    print(f"rep chains ST-only:      {len(st_only)}")

    if not both:
        print("\n(no common rep chains -- nothing to compare)")
        return

    cov_match = 0
    nmem_match = 0
    mismatches = []  # (abs_rel_diff, key, ru_cov, st_cov, ru_n, st_n)
    for k in both:
        r = ru[k]
        s = st[k]
        rc, sc = r["group_cov"], s["group_cov"]
        denom = max(abs(rc), abs(sc), 1e-9)
        rel = abs(rc - sc) / denom
        if rel <= 0.05:
            cov_match += 1
        if r["n_members"] == s["n_members"]:
            nmem_match += 1
        mismatches.append((rel, k, rc, sc, r["n_members"], s["n_members"]))

    n = len(both)
    print()
    print(f"group_cov match (within 5%): {cov_match}/{n} "
          f"({100.0 * cov_match / n:.1f}%)")
    print(f"n_members match:             {nmem_match}/{n} "
          f"({100.0 * nmem_match / n:.1f}%)")

    mismatches.sort(reverse=True)
    print()
    print("=== Largest group_cov mismatches (Rustle vs ST) ===")
    print(f"{'rel_diff':>9}  {'strand':>6}  {'ru_cov':>10}  {'st_cov':>10}  "
          f"{'ru_n':>5}  {'st_n':>5}  rep_introns")
    for rel, k, rc, sc, rn, sn in mismatches[:20]:
        if rel <= 0.05:
            break
        strand, introns = k
        disp = introns if len(introns) <= 60 else introns[:57] + "..."
        print(f"{rel:>9.3f}  {strand:>6}  {rc:>10.3f}  {sc:>10.3f}  "
              f"{rn:>5}  {sn:>5}  {disp}")


if __name__ == "__main__":
    main()
