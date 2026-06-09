"""Analysis (c): can any feature separate real-missing from spurious-missing ST chains?

AUC computed via the Mann-Whitney U statistic (ties counted as 0.5) -> no numpy needed.
"""
import json, sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import lib

FEATURES = ["longcov", "cov", "entry_abund", "nexons"]

def auc_for_feature(rows, feat):
    pos = [r[feat] for r in rows if r["label"] == "real"]
    neg = [r[feat] for r in rows if r["label"] == "spurious"]
    if not pos or not neg:
        return float("nan")
    wins = 0.0
    for p in pos:
        for n in neg:
            wins += 1.0 if p > n else (0.5 if p == n else 0.0)
    return wins / (len(pos) * len(neg))

def build_rows():
    misses = [json.loads(l) for l in open(f"{lib.RG}/st_only_classified.jsonl")]
    ann_by_chrom = {}
    rows = []
    for m in misses:
        chrom = m["chrom"]
        wd = f"{lib.CACHE}/{chrom}/{m['tid']}"
        if not os.path.exists(f"{wd}/.done"):
            continue
        if chrom not in ann_by_chrom:
            ref = lib.parse_gtf(f"{lib.PERCHROM}/{chrom}/ref.gtf")
            ann_by_chrom[chrom] = {d["introns"] for d in ref.values() if d["introns"]}
        ann = ann_by_chrom[chrom]
        rustle = {lib.parse_chain_str(e["payload"]["introns"])
                  for e in lib.read_parity(f"{wd}/r.jsonl", "path_extracted")
                  if e["payload"].get("introns")}
        seen = set()
        for e in lib.read_parity(f"{wd}/s.jsonl", "path_extracted"):
            s = e["payload"].get("introns")
            if not s:
                continue
            ch = lib.parse_chain_str(s)
            if ch in rustle or ch in seen:
                continue
            seen.add(ch)
            p = e["payload"]
            rows.append({
                "label": "real" if ch in ann else "spurious",
                "longcov": float(p.get("longcov", 0.0)),
                "cov": float(p.get("cov", 0.0)),
                "entry_abund": float(p.get("entry_abund", 0.0)),
                "nexons": float(p.get("nexons", 0)),
            })
    return rows

def main():
    rows = build_rows()
    nr = sum(1 for r in rows if r["label"] == "real")
    ns = sum(1 for r in rows if r["label"] == "spurious")
    print(f"=== separability: ST chains rustle misses (real={nr} spurious={ns}) ===")
    for feat in FEATURES:
        auc = auc_for_feature(rows, feat)
        print(f"  {feat:12s} AUC(real>spurious) = {auc:.3f}")
    print("\nGate input: any AUC >> 0.5 (or << 0.5) means a discriminator exists.")
    with open(f"{lib.ROOT}/bench/flow_recall_phase0/separability_rows.jsonl", "w") as fh:
        for r in rows:
            fh.write(json.dumps(r) + "\n")

if __name__ == "__main__":
    main()
