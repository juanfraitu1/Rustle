import os, sys
sys.path.insert(0, os.path.dirname(__file__))
from trunc5 import truncated

A, B, C = (100, 200), (300, 400), (500, 600)  # introns left to right

def test_strand_decides_which_side_is_five_prime():
    cnt = {(B, C): 2, (A, B, C): 3}                       # the longer reads continue on the left
    assert truncated(cnt, {(B, C): '+'}) == {(B, C)}      # left = 5' on plus: truncation
    assert truncated(cnt, {(B, C): '-'}) == set()         # left = 3' on minus: kept
    cnt = {(A, B): 2, (A, B, C): 2}                       # continue on the right
    assert truncated(cnt, {(A, B): '-'}) == {(A, B)} and truncated(cnt, {(A, B): '+'}) == set()

def test_needs_a_full_read_and_at_least_as_many_container_reads():
    assert truncated({(A, B, C): 5}, {(B, C): '+'}) == set()                 # no read spells the chain: kept
    assert truncated({(B, C): 3, (A, B, C): 2}, {(B, C): '+'}) == set()      # fewer container reads than full reads
    assert truncated({(B, C): 1, (A, C): 9}, {(B, C): '+'}) == set()         # not a contiguous run of the longer chain
