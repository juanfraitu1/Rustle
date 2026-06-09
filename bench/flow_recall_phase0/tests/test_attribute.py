import os, sys
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import attribute

def test_flow_enumeration_when_all_junctions_present():
    ref = ((201, 300), (501, 800))
    junctions = {(201, 300), (501, 800), (900, 950)}
    assert attribute.attribute_cause(ref, junctions) == "flow_enumeration"

def test_graph_missing_when_a_junction_absent():
    ref = ((201, 300), (501, 800))
    junctions = {(201, 300)}              # second intron not accepted
    assert attribute.attribute_cause(ref, junctions) == "graph_missing"
