import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import tolerant as T  # noqa: E402


def test_a_cluster_needs_its_reads_beyond_delta_from_every_locus():
    de = {"a": 0.02, "b": 0.03, "c": 0.001, "d": 0.0}
    assert T.median_de({"a", "b", "c"}, de) == 0.02
    assert T.beyond({"a", "b", "c"}, de)               # median 0.02 > delta
    assert not T.beyond({"a", "c", "d"}, de)           # median 0.001: the genome explains the reads
    assert T.median_de({"x"}, de) is None and not T.beyond({"x"}, de)
