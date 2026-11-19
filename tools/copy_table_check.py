#!/usr/bin/env python3
"""The rules a copy table (PREFIX.fam.copies.{tsv,fa}) must keep for `assign` (audit 2026-11-19, finding 13):
every family's copy_idx runs 0..n-1, (family_id, copy_idx) is unique, and the FASTA holds exactly one record per row, headed
`family|copy_idx|...`.

  copy_table_check.py PREFIX          exit 0 when PREFIX.tsv / PREFIX.fa keep the rules, else 1 with the first problems
  copy_table_check.py --selftest"""
import collections, os, sys, tempfile

def problems(prefix, limit=10):
    out = []
    rows = collections.defaultdict(list)
    with open(prefix + ".tsv") as fh:
        head = fh.readline().rstrip("\n").split("\t")
        fi, ci = head.index("family_id"), head.index("copy_idx")
        for n, l in enumerate(fh, 2):
            f = l.rstrip("\n").split("\t")
            try:
                rows[f[fi]].append(int(f[ci]))
            except (IndexError, ValueError):
                out.append(f"{prefix}.tsv line {n}: no integer copy_idx")
    keys = set()
    for fam, idx in rows.items():
        if sorted(idx) != list(range(len(idx))):
            out.append(f"{fam}: copy_idx {sorted(idx)[:8]} is not 0..{len(idx) - 1}")
        keys.update((fam, i) for i in idx)
    seen = collections.Counter()
    with open(prefix + ".fa") as fh:
        for l in fh:
            if l.startswith(">"):
                p = l[1:].split("|")
                seen[(p[0], int(p[1]) if len(p) > 1 and p[1].isdigit() else -1)] += 1
    dup = [k for k, c in seen.items() if c > 1]
    if dup:
        out.append(f"{len(dup)} FASTA record(s) repeated, first {dup[0]}")
    if set(seen) - keys:
        out.append(f"{len(set(seen) - keys)} FASTA record(s) with no table row, first {sorted(set(seen) - keys)[0]}")
    if keys - set(seen):
        out.append(f"{len(keys - set(seen))} table row(s) with no FASTA record, first {sorted(keys - set(seen))[0]}")
    return out[:limit]

def _selftest():
    d = tempfile.mkdtemp()
    p = os.path.join(d, "t")
    hdr = "family_id\tcopy_idx\ttid\n"
    open(p + ".tsv", "w").write(hdr + "F1\t0\ta\nF1\t1\tb\nF2\t0\tc\n")
    open(p + ".fa", "w").write(">F1|0|c:1-2|+|nexon=1\nA\n>F1|1|x\nC\n>F2|0|x\nG\n")
    assert problems(p) == [], problems(p)
    open(p + ".tsv", "w").write(hdr + "F1\t0\ta\nF1\t2\tb\nF2\t0\tc\n")
    assert any("not 0..1" in x for x in problems(p)), problems(p)
    open(p + ".tsv", "w").write(hdr + "F1\t0\ta\nF1\t1\tb\n")
    assert any("no table row" in x for x in problems(p)), problems(p)
    open(p + ".tsv", "w").write(hdr + "F1\t0\ta\nF1\t1\tb\nF2\t0\tc\nF3\t0\td\n")
    assert any("no FASTA record" in x for x in problems(p)), problems(p)

if __name__ == "__main__":
    if sys.argv[1:] == ["--selftest"]:
        _selftest(); print("ok"); sys.exit(0)
    if len(sys.argv) != 2:
        sys.exit(__doc__)
    bad = problems(sys.argv[1])
    for x in bad:
        print(f"[copy_table_check] {x}", file=sys.stderr)
    sys.exit(1 if bad else 0)
