#!/usr/bin/env python3
"""Do the simulation's families hold domain sharers? (user, 2026-10-16: "lets test if the links that cover >= 80-90% of the
shorter transcript eliminate the domain sharers")

    domains.py test FAMILIES_ALL.json OUT.json [--covers 0.5,0.8,0.9]

FAMILIES_ALL.json: a `families.py build --cap 1000` output (every linked gene, link rule >= half the shorter transcript).
Per family, per link-coverage threshold: the largest component under the transcript rule (families.largest_component) is
KEPT, the rest DROPPED. An independent yardstick grades every gene against the kept genes: each isoform aligned, spliced
(minimap2 -cx splice:hq), to the other gene's genomic LOCUS (+-1 kb), blocks at >= 0.90 identity, both directions (grade):
`whole` (each locus holds >= 80% of the other's best isoform), `fragment` (one direction only: a truncated copy), `partial`
(50-80%), `segment` (< 50% both ways: they share only a part, a DOMAIN SHARER); a gene takes its best grade over the kept
genes. Aligning to the partner's genome, not its annotation, keeps intron divergence (old duplicates: SSX, AMY) and the
canonical-transcript choice (NPIPB6-9) out of the yardstick. A threshold eliminates domain sharers when the kept genes are
whole copies (or fragments) and the dropped ones are segments.
"""
import argparse
import collections
import importlib.util
import json
import os
import subprocess
import sys
import tempfile

HERE = os.path.dirname(os.path.abspath(__file__))
WHOLE, PART, LOCUS_IDENT, LOCUS_PAD = 0.8, 0.5, 0.90, 1000


def _load(name):
    spec = importlib.util.spec_from_file_location(f"o3_graph_{name}", os.path.join(HERE, f"{name}.py"))
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


_F = _load("families")
union_cover, tx_locus_cover = _F.union_cover, _F.tx_locus_cover    # the yardstick is the build's `--link locus` measure


RANK = {"whole": 3, "fragment": 2, "partial": 1, "segment": 0}


def grade(cab, cba, whole=WHOLE, part=PART):
    """A gene pair by how much of each gene's transcripts the other gene's LOCUS holds (cab: A's best isoform on B's locus):
    whole copies (both >= whole), a fragment (one only), partial, or a segment (a domain sharer: both < part)."""
    hi, lo = max(cab, cba), min(cab, cba)
    return "whole" if lo >= whole else "fragment" if hi >= whole else "partial" if hi >= part else "segment"


def _asm20(seqs):
    with tempfile.TemporaryDirectory() as d:
        with open(f"{d}/s.fa", "w") as fh:
            fh.writelines(f">{k}\n{v}\n" for k, v in seqs.items())
        return subprocess.run(["minimap2", "-cx", "asm20", "-X", f"{d}/s.fa", f"{d}/s.fa"], capture_output=True, text=True,
                              check=True).stdout.splitlines()


def family_test(f, fa, F, covers):
    names = [c["gene_id"] for c in f["copies"]]
    by = {c["gene_id"]: c for c in f["copies"]}
    tx = {g: F.spliced(fa, by[g]["chrom"], by[g]["strand"], [tuple(x) for x in by[g]["transcripts"][by[g]["main"]]]) for g in names}
    hits = [(x[0], x[5], int(x[9]) / int(x[10]), int(x[10]), int(x[1]), int(x[6])) for x in (ln.split("\t") for ln in _asm20(tx))]
    best_tx = collections.defaultdict(float)
    for a, b, ident, blk, la, lb in hits:
        if a != b and ident >= F.MIN_IDENT:
            for g in (a, b):
                best_tx[g] = max(best_tx[g], blk / min(la, lb))
    with tempfile.TemporaryDirectory() as d:                        # every isoform, spliced, onto every gene's locus
        with open(f"{d}/q.fa", "w") as fh:
            for g in names:
                for t, ex in by[g]["transcripts"].items():
                    fh.write(f">{g}|{t}\n{F.spliced(fa, by[g]['chrom'], by[g]['strand'], [tuple(x) for x in ex])}\n")
        with open(f"{d}/t.fa", "w") as fh:
            for g in names:
                fh.write(f">{g}\n{fa.fetch(by[g]['chrom'], max(0, by[g]['start'] - LOCUS_PAD), by[g]['end'] + LOCUS_PAD).upper()}\n")
        paf = subprocess.run(["minimap2", "-cx", "splice:hq", "-N", "1000", "-p", "0.01", "--secondary=yes", f"{d}/t.fa", f"{d}/q.fa"],
                             capture_output=True, text=True, check=True).stdout.splitlines()
    rows = [(x[0], int(x[1]), int(x[2]), int(x[3]), x[5], int(x[9]) / int(x[10])) for x in (ln.split("\t") for ln in paf)]
    cov = tx_locus_cover(rows)
    out = {}
    for t in covers:
        kept = F.largest_component(names, hits, min_cover=t)
        dropped = [n for n in names if n not in kept]

        def best(n):
            gs = [grade(cov.get((n, m), 0.0), cov.get((m, n), 0.0)) for m in kept if m != n]
            return max(gs, key=RANK.get) if gs else "segment"
        out[str(t)] = dict(kept=len(kept), dropped=[dict(gene=by[n]["gene_name"], gene_id=n, cls=best(n), tx_cover=round(best_tx[n], 3))
                                                    for n in dropped],
                           kept_class=dict(collections.Counter(best(n) for n in kept)) if len(kept) > 1 else {},
                           kept_not_whole=[dict(gene=by[n]["gene_name"], cls=best(n)) for n in kept
                                           if len(kept) > 1 and best(n) != "whole"])
    return out


def cmd_test(a):
    import pysam
    F = _load("families")
    fa = pysam.FastaFile(F.GENOME)
    covers = [float(x) for x in a.covers.split(",")]
    res, tot = {}, collections.defaultdict(collections.Counter)
    for st, fs in json.load(open(a.families))["sets"].items():
        for f in fs:
            if len(f["copies"]) < 2:
                continue
            r = family_test(f, fa, F, covers)
            res[f"{st}/{f['family']}"] = r
            for t, v in r.items():
                tot[t]["families"] += 1
                tot[t]["kept"] += v["kept"]
                tot[t]["dropped"] += len(v["dropped"])
                tot[t]["families_changed"] += bool(v["dropped"])
                for k, n in v["kept_class"].items():
                    tot[t][f"kept_{k}"] += n
                for x in v["dropped"]:
                    tot[t][f"dropped_{x['cls']}"] += 1
    json.dump(dict(totals={t: dict(v) for t, v in tot.items()}, families=res), open(a.out, "w"), indent=1)
    for t, v in tot.items():
        print(f"link cover >= {t}: " + ", ".join(f"{k} {n}" for k, n in sorted(v.items())))
    return 0


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd", required=True)
    p = sub.add_parser("test")
    p.add_argument("families")
    p.add_argument("out")
    p.add_argument("--covers", default="0.5,0.8,0.9")
    p.set_defaults(func=cmd_test)
    a = ap.parse_args(argv)
    return a.func(a)


if __name__ == "__main__":
    sys.exit(main())
