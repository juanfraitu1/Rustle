"""The mutation operators. Each takes a GeneModel (mutated IN PLACE), a spec dict and a random.Random, and appends a
record with the realised coordinates (copy-local, AFTER the op) to `model.ops`. Exons are numbered 1..n in the spec
(transcript order, counting every DNA exon, in_rna or not).

Hygiene (memory: tandem-sim lessons): splice dinucleotides are never mutated by snp/indel; every exon inserted carries
AG/GT flanks so its junctions are canonical; an inverted exon leaves the RNA chain (its sites face the wrong way).
"""
import random

from .model import BASES, Exon, GeneModel, rc, random_seq


def _exon_index(model, spec, key="exon"):
    i = spec.get(key)
    if i is None:
        raise ValueError(f"op {spec.get('op')} needs '{key}'")
    if isinstance(i, str):
        if i == "middle":
            return len(model.exons) // 2
        if i == "first":
            return 0
        if i == "last":
            return len(model.exons) - 1
        for k, e in enumerate(model.exons):
            if e.label == i:
                return k
        raise ValueError(f"no exon labelled {i}")
    i = int(i)
    if not 1 <= i <= len(model.exons):
        raise ValueError(f"exon {i} out of range 1..{len(model.exons)}")
    return i - 1


def _rec(model, **kw):
    model.ops.append(kw)
    return model


# ---------------------------------------------------------------- point changes
def op_snp(model, spec, rng):
    """`rate` substitutions per base over the region (`all` | `exons` | `introns`), exactly round(rate * L) of them at
    distinct positions, splice dinucleotides protected."""
    rate = float(spec["rate"])
    region = spec.get("region", "all")
    prot = model.protected()
    if region == "exons":
        cand = [p for e in model.exons for p in range(e.start, e.end)]
    elif region == "introns":
        cand = [p for s, e in model.introns() for p in range(s, e) if p not in prot]
    else:
        cand = [p for p in range(len(model.seq)) if p not in prot]
    n = min(int(round(rate * len(cand))), len(cand))
    pos = sorted(rng.sample(cand, n))
    s = list(model.seq)
    for p in pos:
        s[p] = rng.choice([b for b in BASES if b != s[p]])
    model.seq = "".join(s)
    return _rec(model, op="snp", rate=rate, region=region, n=n, positions=pos)


def op_indel(model, spec, rng):
    """`rate` indels per base (round(rate * L) of them), each 1..`max_len` bp, deletion or insertion with equal
    probability, inside one exon or one intron and away from splice sites (redrawn otherwise)."""
    rate, max_len = float(spec["rate"]), int(spec.get("max_len", 3))
    n = int(round(rate * len(model.seq)))
    events = []
    tries = 0
    while len(events) < n and tries < 50 * n + 50:
        tries += 1
        p = rng.randrange(0, len(model.seq))
        k = rng.randint(1, max_len)
        if rng.random() < 0.5:
            try:
                model.edit(p, k, "")
                events.append({"pos": p, "del": k})
            except ValueError:
                continue
        else:
            ins = random_seq(k, rng)
            try:
                model.edit(p, 0, ins)
                events.append({"pos": p, "ins": ins})
            except ValueError:
                continue
    return _rec(model, op="indel", rate=rate, max_len=max_len, n=len(events), events=events)


# ---------------------------------------------------------------- exon count
def op_exon_delete(model, spec, rng):
    """Remove exon `exon` from the DNA. Internal exon: its bases go, the two introns merge (GT of the left one, AG of
    the right one kept, so the chain exon-1 -> exon+1 is canonical). Terminal exon: the exon AND its adjacent intron go
    (the copy still starts/ends at an exon)."""
    i = _exon_index(model, spec)
    e = model.exons[i]
    n = len(model.exons)
    if n < 2:
        raise ValueError("cannot delete the only exon")
    if i == 0:
        a, b = 0, model.exons[1].start
    elif i == n - 1:
        a, b = model.exons[n - 2].end, len(model.seq)
    else:
        a, b = e.start, e.end
    label = e.label
    removed = model.seq[a:b]
    # drop the exon from the list first so edit() sees the interval as intronic
    model.exons = [x for k, x in enumerate(model.exons) if k != i]
    if i == 0 or i == n - 1:
        # the removed interval includes an intron end; bypass the splice-site guard deliberately
        delta = -(b - a)
        model.seq = model.seq[:a] + model.seq[b:]
        model.exons = [Exon(x.start + (delta if x.start >= b else 0), x.end + (delta if x.start >= b else 0),
                            x.label, x.in_rna, x.inverted) for x in model.exons]
        model.check()
    else:
        model.edit(a, b - a, "")
    return _rec(model, op="exon_delete", exon=i + 1, label=label, removed_bp=len(removed), at=a)


def op_splice_kill(model, spec, rng):
    """Kill the donor of internal exon `exon` (GT -> CT, 2 bp of DNA): the chain skips the exon, the DNA keeps it."""
    i = _exon_index(model, spec)
    if i == 0 or i == len(model.exons) - 1:
        raise ValueError("splice_kill needs an internal exon (a terminal exon has no donor/acceptor pair to kill)")
    e = model.exons[i]
    d = e.end
    model.seq = model.seq[:d] + "CT" + model.seq[d + 2:]
    e.in_rna = False
    return _rec(model, op="splice_kill", exon=i + 1, label=e.label, donor_at=d, was="GT", now="CT")


def op_exon_insert(model, spec, rng):
    """Insert a new exon into the intron after exon `after` (1..n-1): AG + exon + GT at `offset` bp into the intron
    (default: its middle), at least 30 bp from either existing splice site. Sequence: random of `length` bp, or a copy
    of exon `source_exon` (an internal duplication). Label `ins<k>` or `dup<label>`."""
    i = _exon_index(model, spec, "after")
    if i >= len(model.exons) - 1:
        raise ValueError("'after' must be an exon followed by an intron")
    s, e = model.exons[i].end, model.exons[i + 1].start
    if e - s < 100:
        raise ValueError("intron too short to host an exon (< 100 bp)")
    if "source_exon" in spec:
        j = _exon_index(model, spec, "source_exon")
        src = model.exons[j]
        new = model.seq[src.start:src.end]
        label = f"dup{src.label}"
    else:
        L = int(spec.get("length", 120))
        new = random_seq(L, rng)
        label = f"ins{sum(1 for o in model.ops if o['op'] == 'exon_insert') + 1}"
    off = spec.get("offset")
    pos = s + (e - s) // 2 if off is None else s + int(off)
    pos = min(max(pos, s + 30), e - 30)
    ins = "AG" + new + "GT"
    # insertion inside an intron, away from its ends: edit() accepts it; then register the exon
    model.edit(pos, 0, ins)
    model.exons.append(Exon(pos + 2, pos + 2 + len(new), label))
    model.exons.sort(key=lambda x: x.start)
    model.check()
    return _rec(model, op="exon_insert", after=i + 1, label=label, start=pos + 2, end=pos + 2 + len(new), length=len(new),
                source=spec.get("source_exon"))


def op_exon_shuffle(model, spec, rng):
    """Swap the sequences (and labels) of exons `a` and `b`: the chain reads ... b ... a ... (the generalised
    `missing-copy shuffled`)."""
    ia, ib = sorted((_exon_index(model, spec, "a"), _exon_index(model, spec, "b")))
    if ia == ib:
        raise ValueError("a and b must differ")
    A, B = model.exons[ia], model.exons[ib]
    sa, sb = model.seq[A.start:A.end], model.seq[B.start:B.end]
    la, lb = A.label, B.label
    # replace B first (downstream) so A's coordinates stay valid
    model.exons[ib].label = la
    _replace_exon(model, ib, sa)
    model.exons[ia].label = lb
    _replace_exon(model, ia, sb)
    return _rec(model, op="exon_shuffle", a=ia + 1, b=ib + 1, labels=[la, lb])


def _replace_exon(model, i, new):
    e = model.exons[i]
    old_len = e.end - e.start
    # edit() refuses an interval that equals a whole exon only when it touches a splice site; an exon's own bases never do
    model.edit(e.start, old_len, new)


# ---------------------------------------------------------------- inversions
def op_invert(model, spec, rng):
    """Reverse-complement a window in place.
      exon: N      that exon (coordinates unchanged, marked inverted, dropped from the chain unless keep_in_rna)
      intron: N    the interior of intron N (first/last 6 bp kept, so splicing is unaffected: silent at RNA)
      span: [a,b]  copy-local window; exons fully inside are reversed (order and sequence) and dropped from the chain,
                   no exon or splice site may be cut
      whole: true  is handled by the planter as a strand flip (same transcript, other strand) and recorded there."""
    if spec.get("whole"):
        return _rec(model, op="invert", whole=True)
    if "exon" in spec:
        i = _exon_index(model, spec)
        e = model.exons[i]
        a, b = e.start, e.end
        kind = f"exon {i + 1}"
    elif "intron" in spec:
        k = int(spec["intron"])
        intr = model.introns()
        if not 1 <= k <= len(intr):
            raise ValueError(f"intron {k} out of range 1..{len(intr)}")
        s, e = intr[k - 1]
        if e - s < 40:
            raise ValueError("intron too short to invert its interior")
        a, b = s + 6, e - 6
        kind = f"intron {k}"
    else:
        a, b = int(spec["span"][0]), int(spec["span"][1])
        kind = "span"
    if not 0 <= a < b <= len(model.seq):
        raise ValueError("inversion window outside the sequence")
    prot = model.protected()
    for e in model.exons:
        inside = e.start >= a and e.end <= b
        outside = e.end <= a or e.start >= b
        if not inside and not outside:
            raise ValueError(f"inversion [{a},{b}) cuts exon {e}")
    if "intron" in spec or "span" in spec:
        if any(p in prot for p in (a, a - 1, b, b - 1)) and "span" in spec:
            raise ValueError("inversion window touches a splice site")
    model.seq = model.seq[:a] + rc(model.seq[a:b]) + model.seq[b:]
    flipped = []
    for e in model.exons:
        if e.start >= a and e.end <= b:
            ns, ne = a + (b - e.end), a + (b - e.start)
            e.start, e.end = ns, ne
            e.inverted = not e.inverted
            e.in_rna = bool(spec.get("keep_in_rna", False))
            flipped.append(e.label)
    model.exons.sort(key=lambda x: x.start)
    model.check()
    return _rec(model, op="invert", kind=kind, start=a, end=b, exons_inverted=flipped, keep_in_rna=bool(spec.get("keep_in_rna", False)))


# ---------------------------------------------------------------- partial copies
def op_truncate(model, spec, rng):
    """Drop the first (`side` 5) or last (`side` 3) `exons` exons with their introns, or `bp` bases from that end; a cut
    inside an exon shortens it, a cut inside an intron is moved to the next exon boundary so the copy still starts and
    ends with an exon."""
    side = str(spec.get("side", 5))
    n = len(model.exons)
    if "exons" in spec:
        k = int(spec["exons"])
        if not 1 <= k < n:
            raise ValueError(f"truncate exons must be 1..{n - 1}")
        cut = model.exons[k].start if side == "5" else model.exons[n - k - 1].end
    else:
        bp = int(spec["bp"])
        cut = bp if side == "5" else len(model.seq) - bp
    if side == "5":
        # move a cut that lands in an intron forward to the next exon start
        for e in model.exons:
            if e.start <= cut < e.end:
                break
            if cut < e.start:
                cut = e.start
                break
        removed = cut
        model.seq = model.seq[cut:]
        new = []
        for e in model.exons:
            if e.end <= cut:
                continue
            new.append(Exon(max(e.start, cut) - cut, e.end - cut, e.label, e.in_rna, e.inverted))
        model.exons = new
    else:
        for e in reversed(model.exons):
            if e.start < cut <= e.end:
                break
            if cut > e.end:
                cut = e.end
                break
        removed = len(model.seq) - cut
        model.seq = model.seq[:cut]
        model.exons = [Exon(e.start, min(e.end, cut), e.label, e.in_rna, e.inverted) for e in model.exons if e.start < cut]
    model.check()
    return _rec(model, op="truncate", side=side, removed_bp=removed, exons_left=len(model.exons),
                labels=[e.label for e in model.exons])


# ---------------------------------------------------------------- gene conversion
def op_convert(model, spec, rng, donors=None):
    """Replace exon `exon` (by label, matched in the donor by label) with the DONOR copy's current sequence of it;
    `introns: true` also converts the flanking intron halves? (no — exon only; a span conversion needs identical
    coordinate systems and is refused unless both models have the same exon table)."""
    donor = (donors or {}).get(spec["from"])
    if donor is None:
        raise ValueError(f"convert: donor copy {spec['from']!r} must be defined BEFORE this copy")
    if "exon" in spec:
        i = _exon_index(model, spec)
        label = model.exons[i].label
        src = next((e for e in donor.exons if e.label == label), None)
        if src is None:
            raise ValueError(f"donor has no exon labelled {label}")
        new = donor.seq[src.start:src.end]
        _replace_exon(model, i, new)
        return _rec(model, op="convert", source=spec["from"], exon=i + 1, label=label, bp=len(new))
    a, b = int(spec["span"][0]), int(spec["span"][1])
    if [(e.start, e.end) for e in model.exons] != [(e.start, e.end) for e in donor.exons] or len(model.seq) != len(donor.seq):
        raise ValueError("span conversion needs identical exon tables (apply it before indels / structural ops)")
    model.seq = model.seq[:a] + donor.seq[a:b] + model.seq[b:]
    return _rec(model, op="convert", source=spec["from"], start=a, end=b, bp=b - a)


def op_intron_resize(model, spec, rng):
    """Resize intron `intron` to `length` bp, keeping its first and last 20 bp (splice sites and flanks)."""
    k, L = int(spec["intron"]), int(spec["length"])
    intr = model.introns()
    if not 1 <= k <= len(intr):
        raise ValueError(f"intron {k} out of range")
    if L < 60:
        raise ValueError("resized intron must be >= 60 bp")
    s, e = intr[k - 1]
    keep = 20
    mid_old = e - s - 2 * keep
    mid_new = L - 2 * keep
    if mid_new <= mid_old:
        model.edit(s + keep, mid_old - mid_new, "")
    else:
        model.edit(s + keep, 0, random_seq(mid_new - mid_old, rng))
    return _rec(model, op="intron_resize", intron=k, old=e - s, new=L)


OPS = {
    "snp": op_snp, "indel": op_indel, "exon_delete": op_exon_delete, "splice_kill": op_splice_kill,
    "exon_insert": op_exon_insert, "exon_shuffle": op_exon_shuffle, "invert": op_invert, "truncate": op_truncate,
    "convert": op_convert, "intron_resize": op_intron_resize,
}


def apply_ops(model, ops, rng, donors=None):
    """Apply a list of op specs in order to `model` (in place). `donors`: {copy_id: GeneModel} already built."""
    for k, spec in enumerate(ops):
        name = spec.get("op")
        if name not in OPS:
            raise ValueError(f"unknown op {name!r} (known: {', '.join(OPS)})")
        sub = random.Random(rng.getrandbits(62))   # one stream per op, so an edit to op k leaves ops < k unchanged
        if name == "convert":
            op_convert(model, spec, sub, donors)
        else:
            OPS[name](model, spec, sub)
    return model
