import importlib.util
import os
import random

HERE = os.path.dirname(os.path.abspath(__file__))


def _load(name):
    spec = importlib.util.spec_from_file_location(f"o3_graph_{name}", os.path.join(HERE, f"{name}.py"))
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


T = _load("trim")
A5 = "AGATGTGTATAAGAGACAGGGGGGGGGGTGGTATCAACGCAGAGTACG"
A3 = "AAAAAAAAAAAAAAAGCATCGCTGTCTCTTATACACATCT"


def _mut(s, n, rng):
    s = list(s)
    for p in rng.sample(range(len(s)), n):
        s[p] = rng.choice([b for b in "ACGT" if b != s[p]])
    return "".join(s)


def test_adapters_with_a_few_errors_are_trimmed_from_both_ends():
    rng = random.Random(2)
    insert = "".join(rng.choice("ACGT") for _ in range(600))
    read = _mut(A5, 2, rng) + insert + _mut(A3, 2, rng)
    assert T.trim_read(read, A5, A3) == insert


def test_a_partial_adapter_is_trimmed_and_a_read_without_one_is_left_alone():
    rng = random.Random(3)
    insert = "".join(rng.choice("ACGT") for _ in range(500))
    assert T.trim_read(A5[-30:] + insert + A3[:25], A5, A3) == insert
    assert T.trim_read(insert, A5, A3) == insert
