"""Filter-rescue feasibility: among rustle's GENERATED-then-KILLED chains, can a cov
threshold separate the real (annotated, want to rescue) from the spurious (correctly killed)?

Uses the existing cache: per locus, rustle path_extracted = generated; rustle r.gtf = kept;
killed = generated - kept. Label each killed chain real if its intron chain is annotated on
that chrom, else spurious. Sweep a cov threshold and report rescued-real vs added-FP.

This is the CORRECT population for the rescue decision (rustle's own killed chains), unlike
the Phase 0 separability which measured ST's missing chains.
"""
import json, os, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import lib

def gtf_chain_set(path):
    chains = set()
    for d in lib.parse_gtf(path).values():
        if d["introns"]:
            chains.add(d["introns"])
    return chains

def auc(real_vals, spur_vals):
    if not real_vals or not spur_vals:
        return float("nan")
    wins = 0.0
    for p in real_vals:
        for n in spur_vals:
            wins += 1.0 if p > n else (0.5 if p == n else 0.0)
    return wins / (len(real_vals) * len(spur_vals))

def main():
    misses = [json.loads(l) for l in open(f"{lib.RG}/st_only_classified.jsonl")]
    ann_by_chrom = {}
    real_cov, spur_cov = [], []   # cov of killed chains, by label
    seen = set()
    for m in misses:
        chrom = m["chrom"]
        wd = f"{lib.CACHE}/{chrom}/{m['tid']}"
        if not os.path.exists(f"{wd}/.done"):
            continue
        if chrom not in ann_by_chrom:
            ann_by_chrom[chrom] = gtf_chain_set(f"{lib.PERCHROM}/{chrom}/ref.gtf")
        ann = ann_by_chrom[chrom]
        kept = gtf_chain_set(f"{wd}/r.gtf")
        # generated chains with their max cov
        gen_cov = {}
        for e in lib.read_parity(f"{wd}/r.jsonl", "path_extracted"):
            s = e["payload"].get("introns")
            if not s:
                continue
            ch = lib.parse_chain_str(s)
            gen_cov[ch] = max(gen_cov.get(ch, 0.0), float(e["payload"].get("cov", 0.0)))
        for ch, cov in gen_cov.items():
            if ch in kept:
                continue                      # survived -> not killed
            key = (chrom, ch)
            if key in seen:
                continue
            seen.add(key)
            if ch in ann:
                real_cov.append(cov)          # annotated but killed -> want to rescue
            else:
                spur_cov.append(cov)          # unannotated and killed -> correctly killed

    nr, ns = len(real_cov), len(spur_cov)
    print(f"=== rustle GENERATED-then-KILLED chains (correct rescue population) ===")
    print(f"  real (annotated, want to rescue): {nr}")
    print(f"  spurious (unannotated, correctly killed): {ns}")
    print(f"  cov AUC(real > spurious): {auc(real_cov, spur_cov):.3f}")
    print(f"\n=== cov-threshold sweep (rescue chains with cov >= t) ===")
    print(f"  {'t':>6} {'rescued_real':>12} {'added_FP':>9} {'real:FP':>8}")
    for t in [0, 1, 2, 3, 5, 8, 12, 20, 50]:
        rr = sum(1 for c in real_cov if c >= t)
        fp = sum(1 for c in spur_cov if c >= t)
        ratio = (rr / fp) if fp else float("inf")
        print(f"  {t:>6} {rr:>12} {fp:>9} {ratio:>8.2f}")

if __name__ == "__main__":
    main()
