import importlib.util
import os

HERE = os.path.dirname(os.path.abspath(__file__))
spec = importlib.util.spec_from_file_location("o3_graph_candidates17", os.path.join(HERE, "candidates17.py"))
C = importlib.util.module_from_spec(spec)
spec.loader.exec_module(C)


def test_round_17_takes_the_stems_of_exactly_two_genes():
    two_plus = {"ABC": "^ABC[0-9]+[A-Z]?$", "DEF": "^DEF[0-9]+[A-Z]?$", "XYZ": "^XYZ[0-9]+[A-Z]?$"}
    assert C.pairs(two_plus, {"DEF": "^DEF[0-9]+[A-Z]?$"}) == {"ABC": "^ABC[0-9]+[A-Z]?$", "XYZ": "^XYZ[0-9]+[A-Z]?$"}
