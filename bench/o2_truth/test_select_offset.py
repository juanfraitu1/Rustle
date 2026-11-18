"""NT_SELECT_OFFSET (PREREG r1796): the hash sample of neutral_truth.py select can be shifted to reads not yet looked at."""
import hashlib, os, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import neutral_truth as NT

def test_offset_selects_a_disjoint_quarter():
    names = [f"m84036_230422/{i}/ccs/{i % 97}_{i * 7}" for i in range(2000)]
    old = {q for q in names if not int(hashlib.sha256(q.encode()).hexdigest()[:8], 16) % 100 >= 25}   # the expression before the offset
    a, b = {q for q in names if NT.selected(q, 25)}, {q for q in names if NT.selected(q, 25, 25)}
    assert a == old and not a & b
    assert 400 < len(a) < 600 and 400 < len(b) < 600
    assert all(NT.selected(q, 100, 25) for q in names)
