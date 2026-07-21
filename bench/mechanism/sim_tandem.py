"""V4a: simulate a 3-copy tandem array, delete copies from the reference, recover
via the certificate, and set-match recovered vs deleted copies (permutation-invariant).
Divergent arm proves recovery; identical arm proves the honest K=0 abstain floor.
Uses a seeded PRNG (no Math.random / wall-clock) so runs are reproducible."""
import random, pathlib, subprocess, json, argparse, shutil, re, sys

BASES = "ACGT"


def build_tandem(seed, n_copies, divergence, seed_rng):
    """Return (reference_with_n_tandem_copies, [copy_seqs]). Divergence = per-base
    substitution prob planted independently per copy (shared coordinate frame => PSVs)."""
    rng = random.Random(seed_rng)
    base = list(seed)
    copies = []
    for _ in range(n_copies):
        c = base[:]
        if divergence > 0:
            for i in range(len(c)):
                if rng.random() < divergence:
                    c[i] = rng.choice([b for b in BASES if b != c[i]])
        copies.append("".join(c))
    # simple tandem: copies concatenated with a short spacer
    spacer = "N" * 50
    ref = spacer.join(copies)
    return ref, copies


def _identity(a, b):
    n = min(len(a), len(b))
    if n == 0:
        return 0.0
    same = sum(1 for i in range(n) if a[i] == b[i])
    return same / max(len(a), len(b))


def set_match(recovered, truth):
    """Greedy permutation-invariant matching: each recovered copy -> best-identity truth copy."""
    matches = []
    used = set()
    for ri, r in enumerate(recovered):
        best, best_id = None, -1.0
        for ti, t in enumerate(truth):
            if ti in used:
                continue
            idv = _identity(r, t)
            if idv > best_id:
                best, best_id = ti, idv
        if best is not None:
            used.add(best)
            matches.append((ri, best, best_id))
    return matches
