import os, sys
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import recall_oracle

def test_fsm_union_known_chrom():
    # NC_073227.2 results.jsonl: both=970, rustle_only=53, st_only=20.
    agg = recall_oracle.fsm_union_per_chrom()
    rec = agg["NC_073227.2"]
    assert rec["both"] == 970 and rec["rustle_only"] == 53 and rec["st_only"] == 20
    assert rec["union_fsm"] == 970 + 53 + 20      # 1043
    assert rec["rustle_fsm"] == 970 + 53          # 1023 = current rustle recall

def test_delta_recall_positive():
    agg = recall_oracle.fsm_union_per_chrom()
    assert agg["NC_073227.2"]["union_fsm"] - agg["NC_073227.2"]["rustle_fsm"] == 20
