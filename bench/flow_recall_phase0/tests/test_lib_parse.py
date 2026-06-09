import os, sys
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import lib
FX = os.path.join(os.path.dirname(__file__), "fixtures")

def test_gtf_intron_chain():
    tx = lib.parse_gtf(os.path.join(FX, "mini.gtf"))
    assert tx["t1"]["strand"] == "+"
    assert tx["t1"]["introns"] == ((201, 300), (501, 800))
    assert tx["t1"]["span"] == (100, 900)

def test_chain_str_roundtrip():
    assert lib.chain_str(((201, 300), (501, 800))) == "201-300,501-800"
    assert lib.parse_chain_str("201-300,501-800") == ((201, 300), (501, 800))

def test_parse_path_extracted():
    evs = lib.read_parity(os.path.join(FX, "pe.jsonl"), "path_extracted")
    assert len(evs) == 1
    assert evs[0]["payload"]["introns"] == "201-300,501-800"

def test_parse_junctions():
    js = lib.read_parity(os.path.join(FX, "pe.jsonl"), "junction_accept")
    assert (201, 300) in {(j["start"], j["end"]) for j in js}
