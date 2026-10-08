#!/usr/bin/env python3
"""IsoCon arm -> the common candidate tables.   adapt_isocon.py <workdir>   (workdir = W/isoc, after control_test.py contigs)

cands.tsv: candidate, family, n_transcripts (= contigs in the merge component, Amendment 8's rule), contigs; cands.fa: the new-copy contigs."""
import collections
import os
import shutil
import sys
import types

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.join(HERE, "..", "rna_allele"))
import common as C  # noqa: E402
import control_test  # noqa: E402


def main(w):
    comp, fam_of = control_test.comps_at(types.SimpleNamespace(w=w), C.DELTA)
    members = collections.defaultdict(list)
    for c, cid in comp.items():
        members[cid].append(c)
    with open(f"{w}/cands.tsv", "w") as o:
        o.write("candidate\tfamily\tn_transcripts\tcontigs\n")
        for cid, cs in sorted(members.items()):
            o.write(f"{cid}\t{fam_of[cs[0]]}\t{len(cs)}\t{','.join(sorted(cs))}\n")
    if os.path.exists(f"{w}/contigs_L.fa"):
        shutil.copy(f"{w}/contigs_L.fa", f"{w}/cands.fa")
    else:
        open(f"{w}/cands.fa", "w").close()
    flagged = sum(1 for cs in members.values() if len(cs) >= 2)
    print(f"candidates {len(members)} in {len({fam_of[cs[0]] for cs in members.values()})} families; flagged (>= 2 transcripts) {flagged}")


if __name__ == "__main__":
    main(sys.argv[1])
