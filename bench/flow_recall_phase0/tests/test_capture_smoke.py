import os, sys
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import lib

def test_ref_chain_lookup():
    # XM_004027824.5 on NC_073224.2 is a known container miss; rustle GENERATES it.
    ch = lib.ref_chain("NC_073224.2", "rna-XM_004027824.5")
    assert ch is not None
    span_lo, span_hi, introns = ch
    assert len(introns) == 10            # 11 exons -> 10 introns
    assert span_lo < span_hi
