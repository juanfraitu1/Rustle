import os, sys
sys.path.insert(0, os.path.dirname(__file__))
from sequence_families import covered

def test_union_of_overlapping_alignments():
    paf = ['q1\t100\t0\t40\t+\tr\t500\t0\t40\t40\t40\t60', 'q1\t100\t30\t60\t+\tr\t500\t90\t120\t30\t30\t60',
           'q1\t100\t80\t90\t-\tr\t500\t0\t10\t10\t10\t60', 'q2\t50\t0\t50\t+\tr\t500\t0\t50\t50\t50\t60', 'short line']
    assert covered(paf) == {'q1': 0.7, 'q2': 1.0}
