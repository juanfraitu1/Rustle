import os, sys
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import gen_census

def test_classify_generated():
    ref = ((201, 300), (501, 800))
    rustle_chains = {((201, 300), (501, 800))}
    assert gen_census.classify(ref, rustle_chains) == "generated"

def test_classify_non_generated():
    ref = ((201, 300), (501, 800))
    rustle_chains = {((201, 300), (777, 888))}   # different second intron
    assert gen_census.classify(ref, rustle_chains) == "non_generated"
