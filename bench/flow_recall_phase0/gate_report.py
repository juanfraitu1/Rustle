"""Apply the pre-registered go/no-go gate (spec: recall ceiling + discriminator)."""
import json, sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import lib

RECALL_MIN = 150
AUC_HI, AUC_LO = 0.70, 0.30

def decide(recall_ceiling, best_auc):
    if recall_ceiling < RECALL_MIN:
        return {"verdict": "STOP", "reason": f"recall ceiling {recall_ceiling} < {RECALL_MIN} (small ceiling)"}
    if not (best_auc >= AUC_HI or best_auc <= AUC_LO):
        return {"verdict": "STOP", "reason": f"no discriminator (best AUC {best_auc:.2f} in [{AUC_LO},{AUC_HI}])"}
    return {"verdict": "PROCEED", "reason": f"ceiling {recall_ceiling} >= {RECALL_MIN}, discriminator AUC {best_auc:.2f}"}

def main():
    import recall_oracle, separability
    agg = recall_oracle.fsm_union_per_chrom()
    ceiling = sum(a["st_only"] for a in agg.values())
    rows = []
    sp = f"{lib.ROOT}/bench/flow_recall_phase0/separability_rows.jsonl"
    if os.path.exists(sp):
        rows = [json.loads(l) for l in open(sp)]
    aucs = {f: separability.auc_for_feature(rows, f) for f in separability.FEATURES} if rows else {}
    best_auc = max((max(a, 1 - a) for a in aucs.values() if a == a), default=0.5)
    v = decide(ceiling, best_auc)
    print("=== PHASE 0 GATE ===")
    print(f"  recall ceiling (st_only genome-wide): {ceiling}")
    print(f"  feature AUCs: " + ", ".join(f"{f}={a:.2f}" for f, a in aucs.items()))
    print(f"  best |discrimination| AUC: {best_auc:.2f}")
    print(f"\n  VERDICT: {v['verdict']} — {v['reason']}")

if __name__ == "__main__":
    main()
