#!/usr/bin/env python3
"""The decision rule of O2 from the sample's own variant calls (design 2026-10-22,
docs/superpowers/specs/2026-10-22-o2-variant-calls-design.md). Pure functions: no BAM, VCF or truth is read here (o2vc_data.py
reads the BAM and the VCF, o2vc.py runs the rule and scores it).

The design's premise: an exact AS tie means the read covers no difference between the copies that the REFERENCE shows (r1350,
r1362); where minimap2 is wrong on real reads, the sample's copies differ from the reference's (r1329). The sample's other
reads say how: the reads minimap2 places on a copy with an AS margin >= 10 (its teachers) go to a variant caller (longcallR)
against the one reference, and the calls give each copy a genotype at every position. A tied read is then held against those
genotypes, base by base:

- genotype(): a copy's alleles at a position -- the record's GT alleles whatever its FILTER, the reference base when there is
  no record, UNKNOWN (None) under MIN_DEPTH teacher reads, at an RnaEdit / Dense record or at a GT calling no allele.
- assign(): a column SEPARATES c and d when both are known and share no allele; the read's base (quality >= MIN_BASEQ) VOTES
  for the one whose genotype holds it; c BEATS d with a vote for c and none for d. The read is placed on the copy that beats
  every other copy, else it abstains with the copies no copy beats. An unknown copy is in no separating column, so it is never
  beaten and never beats. No threshold is tuned: MIN_DEPTH and MIN_BASEQ are the caller's own.
- foreign(): a base that no tied copy holds, every copy being known, is FOREIGN (an allele the reference's copies do not
  explain); A>G in the transcript frame is "editing" instead and never counts.
- flags(): foreign alleles that recur and co-occur are a possible MISSING COPY. Two recurrent alleles at different columns of
  a read are linked when >= 2 reads give the pair; a connected component of linked alleles is the missing copy's haplotype,
  assembled allele by allele, and its reads are those carrying its alleles at >= 2 columns.

Everything is compared in the read's own frame, as sequenced (FLNC = transcript sense): a copy aligned opposite to the read has
its alleles complemented (to_read_frame()), or a read would be placed on strand-mismatched evidence. Sets are only built,
counted and tested for membership and every other loop is sorted, so the output depends neither on hash order nor, for
flags(), on the order the reads arrive.
"""
import collections

MIN_DEPTH = 6                                                   # teacher reads for a copy to be known (the caller's HiFi min-depth)
MIN_BASEQ = 10                                                  # a read base votes from this quality (the caller's --min-baseq)
UNKNOWN_FILTERS = frozenset({"RnaEdit", "Dense"})
BASES = frozenset("ACGT")
COMP = str.maketrans("ACGTNacgtn", "TGCANtgcan")
# one VCF record: REF, the ALT alleles, GT as allele indices (0 = REF), FILTER, DP
Call = collections.namedtuple("Call", "ref alts gt filter dp")


def genotype(call, depth, ref_base):
    """A copy's alleles at a position as a frozenset in the genome's strand, or None when the copy is unknown there: fewer than
    MIN_DEPTH teacher reads (depth: the record's DP where the caller wrote one, else the teacher depth counted as the caller
    counts it, o2vc_data.teacher_depth), an RnaEdit or Dense record (a Dense window is what a collapsed extra copy leaves on
    its nearest reference copy), or a record whose GT calls no allele (gt == (), from "./."). Otherwise the record's GT
    alleles whatever its FILTER (a LowQual het is an unphased het: it can only make two copies overlap, never separate them
    falsely), or the reference base when the caller left no record."""
    if depth < MIN_DEPTH or (call is not None and (call.filter in UNKNOWN_FILTERS or not call.gt)):
        return None
    if call is None:
        return frozenset({ref_base})
    alleles = (call.ref, *call.alts)
    return frozenset(alleles[i] for i in call.gt)


def to_read_frame(g, flip):
    """A genotype in the read's frame: complemented when the copy is aligned opposite to the read (flip), as inverted paralogs
    are; None (unknown) stays None."""
    return frozenset(a.translate(COMP) for a in g) if flip and g is not None else g


def assign(copies, cols):
    """(the placed copy or None, the copies no copy beats, the separating columns voting for the placed copy -- 0 when none is
    placed). cols: [dict(base, qual, g)], g = {copy: genotype in the read's frame, None = unknown}, a copy missing from g being
    unknown. A column separates c and d when both genotypes are known and share no allele; it votes, at qual >= MIN_BASEQ, for
    the one whose genotype holds the read's base; c beats d when >= 1 column votes for c against d and none for d against c.
    The read is placed on the copy that beats every other copy; else it abstains with the copies nothing beats. That set can
    hold one copy with no placement (A beats B and B beats C, but no column separates A from C). Fewer than two copies leave
    nothing to separate: no placement."""
    cs = sorted(set(copies))
    if len(cs) < 2:
        return None, frozenset(cs), 0
    votes = collections.defaultdict(list)                       # (c, d) -> the columns voting for c against d
    for n, col in enumerate(cols):
        if col["qual"] < MIN_BASEQ:
            continue
        g, base = col["g"], col["base"]
        known = [c for c in cs if g.get(c) is not None]
        hold = [c for c in known if base in g[c]]
        lack = [c for c in known if base not in g[c]]
        for c in hold:
            for d in lack:
                if not g[c] & g[d]:                             # d is separated from c, and the read's base is c's
                    votes[c, d].append(n)
    beats = sorted((c, d) for c, d in votes if (d, c) not in votes)             # c beats d
    beaten = {d for _, d in beats}
    n_beats = collections.Counter(c for c, _ in beats)
    win = next((c for c in cs if n_beats[c] == len(cs) - 1), None)
    if win is None:
        return None, frozenset(cs) - beaten, 0
    return win, frozenset({win}), len({n for d in cs if d != win for n in votes[win, d]})


def foreign(col):
    """"foreign" when the read's base (one of ACGT, at qual >= MIN_BASEQ) is held by none of the tied copies and every copy is
    known (g lists every tied copy; an unknown copy never stands in with its reference base, or any ordinary SNP of a weakly
    expressed copy would look foreign); "editing" instead when the base is G and every copy holds A (A>G in the read's
    transcript frame); else ""."""
    base, g = col["base"], list(col["g"].values())
    if base not in BASES or col["qual"] < MIN_BASEQ or not g or any(x is None for x in g):
        return ""
    if any(base in x for x in g):
        return ""
    return "editing" if base == "G" and all("A" in x for x in g) else "foreign"


def flags(by_read):
    """The possible missing copies in foreign alleles: [dict(id, alleles, reads)]. by_read: {read: [(column id within the read,
    keys)]}, a key being (copy, chrom, pos, the genome strand's base) -- one foreign allele named once for each tied copy it
    falls on. A key is recurrent when >= 2 reads carry it; two recurrent keys of one read at different columns are a pair, linked
    when >= 2 reads give it (the two copies' keys of one column never link by themselves); a flag is a connected component of
    linked keys with the reads carrying its keys at >= 2 columns. Flags with fewer than two reads are dropped, then numbered
    from 0 in the order of their smallest key; alleles and reads are sorted."""
    carriers = collections.defaultdict(set)                     # key -> the reads carrying it
    for r in sorted(by_read):
        for _, keys in by_read[r]:
            for k in sorted(keys):
                carriers[k].add(r)
    recurrent = {k for k, rs in carriers.items() if len(rs) >= 2}
    given = collections.defaultdict(set)                        # (key, key) -> the reads giving the pair
    for r in sorted(by_read):
        at = sorted((k, i) for i, keys in by_read[r] for k in sorted(keys) if k in recurrent)
        for x, (k1, i) in enumerate(at):
            for k2, j in at[x + 1:]:
                if i != j and k1 != k2:
                    given[k1, k2].add(r)
    parent = {}

    def find(k):
        parent.setdefault(k, k)
        while parent[k] != k:
            parent[k] = parent[parent[k]]
            k = parent[k]
        return k

    for a, b in sorted(p for p, rs in given.items() if len(rs) >= 2):
        ra, rb = find(a), find(b)
        if ra != rb:
            parent[max(ra, rb)] = min(ra, rb)                   # a component's root is its smallest key
    members = collections.defaultdict(list)                     # root -> its keys, sorted
    for k in sorted(parent):
        members[find(k)].append(k)
    carried = collections.defaultdict(lambda: collections.defaultdict(set))    # root -> read -> the columns carrying its keys
    for r in sorted(by_read):
        for i, keys in by_read[r]:
            for k in sorted(keys):
                if k in parent:
                    carried[find(k)][r].add(i)
    out = []
    for root in sorted(members):
        reads = sorted(r for r, ids in carried[root].items() if len(ids) >= 2)
        if len(reads) >= 2:                                     # always so: a linked pair has two reads carrying both keys
            out.append(dict(id=len(out), alleles=members[root], reads=reads))
    return out
