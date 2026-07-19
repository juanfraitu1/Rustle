#!/usr/bin/env python3
"""Real-abpoa smoke for soto_dna_vg. Run via the interpreter that has pyabpoa:
  /home/juanfra/miniforge3/bin/python bench/soto/smoke_soto_dna_vg.py
Builds PMS2P (ID_8) end-to-end and asserts the DNA-ceiling presence claim. Exit 0 = pass.
(A standalone script, not a pytest test, because the pytest interpreter lacks pyabpoa.)"""
import sys
sys.path.insert(0, "bench/soto")
import soto_dna_vg as V

bed = open("bench/soto/80_fams.chr.bed").read().splitlines()
fa = V.read_fasta(V.MEMFA)
detection = V.load_detection(open("bench/soto/soto_member_detection.tsv").read().splitlines())
members = V.parse_family_members(bed, "ID_8")            # PMS2P, moderate ~5-15 kb members
assert len(members) >= 5, members
r = V.build_family("ID_8", members, fa, detection)
assert r["missing"] == [], f"members missing from graph: {r['missing']}"   # the DNA-ceiling claim
assert r["n_present"] == r["n_members"], (r["n_present"], r["n_members"])
kinds = {l[0] for l in r["gfa"].splitlines() if l}
assert {"H", "S", "P"} <= kinds, kinds
assert r["colours"].splitlines()[0] == "Node,Colour" and len(r["colours"].splitlines()) > 1
plines = {l.split("\t")[1] for l in r["gfa"].splitlines() if l.startswith("P\t")}
assert plines == {g for g, _, _, _ in members}, plines
n_nodes = sum(1 for l in r["gfa"].splitlines() if l.startswith("S"))
print(f"[smoke] PMS2P ID_8: {r['n_present']}/{r['n_members']} copies as paths, {n_nodes} nodes  OK")
