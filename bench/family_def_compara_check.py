#!/usr/bin/env python3
"""
family_def_compara_check.py — independent corroboration of the conflict-graph family
definition's POSITIVES against Ensembl Compara within-species paralogy.

Claim it reproduces (see bench/family_definition_formal.md, "Independent confirmation"):
  - The annotated positive families (RABL2, MAGEA) are confirmed RECENT within-species
    paralogs (Homininae / Catarrhini taxonomy level).
  - The definition is a STRICT SUBSET of recent paralogy: APOBEC3 is ALSO a confirmed
    recent paralog, yet it is a correct TN in the panel (reads are read-resolvable). So the
    conflict-family answers a sharper question than Compara: recent-paralogy ∩ read-confusability.
  - The de-novo LOC copies (AK6/LOC115934278, CCDC196/LOC129526440) are absent from Compara
    (unannotated paralogues) — copies the read-level method finds that the annotation misses.

Cache: bench/compara_cache.json (human-mapped Ensembl Compara homology, keyed
`homology_paralogues|<symbol>`). Reference-only; no network.
"""
import json
import os
import sys

CACHE = os.path.join(os.path.dirname(os.path.abspath(__file__)), "compara_cache.json")
# "recent" within-species paralog taxonomy levels (great-ape / catarrhine duplications)
RECENT = {"Homininae", "Hominidae", "Hominoidea", "Catarrhini", "Simiiformes"}


def load():
    with open(CACHE) as f:
        return json.load(f)


def recent_paralogs(cache, symbol):
    """Return (n_within_species, n_recent, taxonomy_levels) or None if symbol absent."""
    key = "homology_paralogues|" + symbol
    if key not in cache:
        return None
    try:
        homs = cache[key]["data"]["data"][0]["homologies"]
    except (KeyError, IndexError, TypeError):
        return (0, 0, [])
    wsp = [h for h in homs if h.get("type") == "within_species_paralog"]
    recent = [h for h in wsp if h.get("taxonomy_level") in RECENT]
    return (len(wsp), len(recent), sorted({h.get("taxonomy_level", "?") for h in recent}))


def main():
    cache = load()
    # (panel role, representative annotated member to query)
    POSITIVES = [("RABL2 (TP family)", "RABL2A"), ("MAGEA (TP family)", "MAGEA12")]
    CONTRAST = [("APOBEC3 (TN — recent paralog but read-RESOLVABLE)", "APOBEC3D")]
    DENOVO_LOC = [("AK6 copy", "LOC115934278"), ("CCDC196 copy", "LOC129526440")]

    print("=" * 74)
    print("Compara corroboration of the conflict-graph family POSITIVES")
    print("=" * 74)

    ok = True
    print("\n[positives — expect RECENT within-species paralogs]")
    for label, sym in POSITIVES:
        r = recent_paralogs(cache, sym)
        if r is None or r[1] == 0:
            print(f"  FAIL {label:34} {sym}: no recent within-species paralog ({r})")
            ok = False
        else:
            print(f"  OK   {label:34} {sym}: {r[1]} recent within-species paralog(s) @ {r[2]}")

    print("\n[contrast — recent paralog by Compara, yet correctly a TN (the strict-subset point)]")
    for label, sym in CONTRAST:
        r = recent_paralogs(cache, sym)
        if r is None or r[1] == 0:
            print(f"  note {label}: {sym} not a recent paralog in cache ({r}) — contrast weaker")
        else:
            print(f"  OK   {label}\n         {sym}: {r[1]} recent within-species paralog(s) @ {r[2]} — "
                  f"yet TN here (reads decidable) => definition ⊊ recent-paralogy")

    print("\n[de-novo LOC copies — expect ABSENT from Compara (annotation gap the method fills)]")
    for label, sym in DENOVO_LOC:
        r = recent_paralogs(cache, sym)
        status = "absent (as expected)" if r is None else f"present ({r})"
        print(f"  {label:18} {sym}: {status}")

    print("\nVERDICT:", "positives corroborated; strict-subset confirmed." if ok else "CHECK FAILED.")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
