"""Usage: python3 bench/g1/misses_g1.py   (light; post hoc, descriptive)

G1 human misses (prereg §5 'if it fails': the missed copies listed by stratum), with fig_loci._write_guided's sets."""
import collections, sys
from pathlib import Path
sys.path.insert(0, "/mnt/linuxdisk/home/juanfraitu/rustle_soto_subset/figures")
import figlib, _liftoff as L
cfg = figlib.load_inputs(None)
t = L.species_target(cfg, "human")
loci, recs = L.read_loci(t.loci), L.read_records(t)
spans = L.Index([{"contig": r["contig"], "iv": [(r["start0"], r["end"])]} for r in recs])
cands = [c for c in L.read_candidates(L.finder_dir(t) / "candidates.tsv") if c["identity"] >= 0.80 - 1e-9]
cidx, allc = L.Index(cands), L.Index(L.read_candidates(L.finder_dir(t) / "candidates.tsv"))
copies = [r for r in loci if r["cls"] == "extra_copy"]
for r in copies:
    r["inside"] = any(L.ov(spans.rows[i]["iv"][0][0], spans.rows[i]["iv"][0][1], r["iv"][0][0], r["iv"][-1][1]) > 0
                      for i in spans.near(r["contig"], [(r["iv"][0][0], r["iv"][-1][1])]))
ref = [r for r in copies if not r["short"] and r["seqid_f"] >= 0.95 - 1e-9 and not r["inside"]]
per_src = collections.Counter(r["source_id"] for r in copies)
miss, found = [], []
for r in ref:
    c, i = cidx.best_cov(r["contig"], r["iv"])
    (found if c >= L.MATCH_COV else miss).append((r, c))
print("ref", len(ref), "found", len(found), "missed", len(miss))
def tab(key):
    a, b = collections.Counter(key(r) for r, _ in found), collections.Counter(key(r) for r, _ in miss)
    for k in sorted(set(a) | set(b), key=lambda k: -(a[k] + b[k])):
        print(f"   {k}: missed {b[k]} of {a[k] + b[k]}")
print("by stratum"); tab(lambda r: r["stratum"])
print("by rtype/biotype"); tab(lambda r: f"{r['rtype']}/{r['biotype']}")
print("by Liftoff copies of the source"); tab(lambda r: "1" if per_src[r["source_id"]] == 1 else ("2-5" if per_src[r["source_id"]] <= 5 else ("6-50" if per_src[r["source_id"]] <= 50 else ">50")))
print("by exon union"); tab(lambda r: "<1 kb" if r["exonic_len"] < 1000 else ("1-5 kb" if r["exonic_len"] < 5000 else ">=5 kb"))
part = [c for r, c in miss if c > 0]
print("misses with some candidate overlap (<50%):", len(part), "| covered >=50% by a candidate below identity 0.80:",
      sum(allc.best_cov(r["contig"], r["iv"])[0] >= L.MATCH_COV for r, _ in miss))
src = collections.Counter(r["source_id"] for r, _ in miss)
print("missed copies come from", len(src), "source records; top:", src.most_common(6))
