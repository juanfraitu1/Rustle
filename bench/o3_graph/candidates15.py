#!/usr/bin/env python3
"""Round 15's family pool (2026-10-17), fixed before any build: every CAT protein-coding gene-name STEM (the name less a final
number and optional letter) with >= 3 genes, none of them matched by any earlier group's pattern (families.GROUPS), built with
round 12's rule; a stem forms a family only where families.py links >= 3 copies (identity >= 0.90 over >= 0.8 of a locus).

    candidates15.py OUT_DIR      -> OUT_DIR/candidates15.json (the stems), OUT_DIR/families15_all.json (the build)
"""
import collections
import gzip
import importlib.util
import json
import os
import re
import sys

HERE = os.path.dirname(os.path.abspath(__file__))


def stems(cat_gff, used_patterns, min_genes=3):
    """{stem: regex} over protein-coding gene names, stems with >= min_genes genes and no gene an earlier pattern matches."""
    names = set()
    with gzip.open(cat_gff, "rt") as fh:
        for ln in fh:
            if "\tgene\t" not in ln or "gene_biotype=protein_coding" not in ln:
                continue
            m = re.search(r"gene_name=([^;\n]+)", ln)
            if m:
                names.add(m.group(1))
    by = collections.defaultdict(set)
    for n in names:
        s = re.sub(r"[0-9]+[A-Z]?$", "", n)
        if s and s != n and "-" not in n and not re.match(r"^(ENSG|MSTRG|LOC|LINC|AC[0-9]|AL[0-9]|AP[0-9])", n):
            by[s].add(n)
    pats = [re.compile(p) for p in used_patterns]
    return {s: rf"^{re.escape(s)}[0-9]+[A-Z]?$" for s, v in sorted(by.items())
            if len(v) >= min_genes and not any(p.match(x) for x in v for p in pats)}


def main(argv):
    out = argv[0]
    os.makedirs(out, exist_ok=True)
    spec = importlib.util.spec_from_file_location("families", os.path.join(HERE, "families.py"))
    F = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(F)
    used = [p for g in F.GROUPS.values() for p in g.values()]
    cand = stems(f"{F.CAT}/chm13v2.0_CAT_Liftoff.slim.gff3.gz", used)
    json.dump(cand, open(f"{out}/candidates15.json", "w"), indent=0)
    print(f"[candidates15] {len(cand)} stems -> {out}/candidates15.json", flush=True)
    if "--build" in argv:
        F.GROUPS = {"auto15": cand}
        return F.main(["build", f"{out}/families15_all.json", "--cap", "1000", "--link", "locus", "--min-cover", "0.8"])
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
