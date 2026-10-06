"""GeneModel: one copy's genomic sequence in TRANSCRIPT orientation plus its exons.

Coordinates are 0-based half-open on `seq`; exons are sorted and non-overlapping, so exon i's donor is the two bases
at `exons[i].end` and exon i+1's acceptor the two bases before `exons[i+1].start`. Every exon carries a `label` (the
template's `e<k>`, or `ins<k>` / `dup<k>` for inserted ones) and `in_rna` (False after splice_kill / an exon inversion:
the DNA keeps the exon, the chain skips it). The model knows nothing about where it is planted; `chromosome` does.
"""
import json
import random

from lib import rc  # noqa: F401  (re-exported for the other modules)

BASES = "ACGT"


class Exon:
    __slots__ = ("start", "end", "label", "in_rna", "inverted")

    def __init__(self, start, end, label, in_rna=True, inverted=False):
        self.start, self.end, self.label, self.in_rna, self.inverted = int(start), int(end), label, bool(in_rna), bool(inverted)

    def __len__(self):
        return self.end - self.start

    def to_json(self):
        return {"start": self.start, "end": self.end, "label": self.label, "in_rna": self.in_rna, "inverted": self.inverted}

    @classmethod
    def from_json(cls, d):
        return cls(d["start"], d["end"], d["label"], d.get("in_rna", True), d.get("inverted", False))

    def __repr__(self):
        return f"Exon({self.start},{self.end},{self.label}{'' if self.in_rna else ',skip'}{',inv' if self.inverted else ''})"


class GeneModel:
    def __init__(self, seq, exons, ops=None, provenance=None):
        self.seq = seq.upper()
        self.exons = [e if isinstance(e, Exon) else Exon(*e) if isinstance(e, (tuple, list)) else Exon.from_json(e) for e in exons]
        self.ops = list(ops or [])
        self.provenance = dict(provenance or {})
        self.check()

    # ---------------------------------------------------------------- invariants
    def check(self):
        L = len(self.seq)
        prev = -1
        for e in self.exons:
            if not (0 <= e.start < e.end <= L):
                raise ValueError(f"exon {e} outside sequence of length {L}")
            if e.start < prev:
                raise ValueError(f"exons overlap or are unsorted at {e}")
            prev = e.end
        if self.exons and (self.exons[0].start != 0 or self.exons[-1].end != L):
            raise ValueError("a model spans exactly its first exon start to its last exon end")
        return self

    def copy(self):
        return GeneModel(self.seq, [Exon(e.start, e.end, e.label, e.in_rna, e.inverted) for e in self.exons],
                         [dict(o) for o in self.ops], dict(self.provenance))

    # ---------------------------------------------------------------- derived
    def introns(self):
        """[(start, end)] of the introns between consecutive exons (DNA, regardless of in_rna)."""
        return [(a.end, b.start) for a, b in zip(self.exons, self.exons[1:])]

    def protected(self):
        """Positions never touched by snp/indel: the donor GT and acceptor AG of every intron."""
        p = set()
        for s, e in self.introns():
            p.update((s, s + 1, e - 2, e - 1))
        return p

    def rna_exons(self, skip=()):
        """The exons of an isoform: every in_rna exon, minus the labels/indices in `skip`."""
        out = []
        for i, e in enumerate(self.exons):
            if e.in_rna and i + 1 not in skip and e.label not in skip:
                out.append(e)
        return out

    def chain_seq(self, skip=()):
        return "".join(self.seq[e.start:e.end] for e in self.rna_exons(skip))

    def chain_junctions(self, skip=()):
        """[(donor_pos, acceptor_pos)] of an isoform: local intron intervals (0-based half-open) between its exons."""
        ex = self.rna_exons(skip)
        return [(a.end, b.start) for a, b in zip(ex, ex[1:])]

    def exon_seqs(self):
        return {e.label: self.seq[e.start:e.end] for e in self.exons}

    def junction_motifs(self):
        """(donor dinucleotide, acceptor dinucleotide) per intron, on the model's own strand."""
        return [(self.seq[s:s + 2], self.seq[e - 2:e]) for s, e in self.introns()]

    def canonical(self):
        return all(d == "GT" and a == "AG" for d, a in self.junction_motifs())

    # ---------------------------------------------------------------- editing primitive
    def edit(self, pos, del_len, ins):
        """Replace seq[pos:pos+del_len] with `ins`, shifting every exon after it. The edited interval must lie entirely
        inside ONE exon or ONE intron and must not touch a splice dinucleotide (the callers guarantee that; this raises
        otherwise), so exon identity is never ambiguous."""
        ins = ins.upper()
        end = pos + del_len
        if not (0 <= pos <= end <= len(self.seq)):
            raise ValueError("edit outside the sequence")
        delta = len(ins) - del_len
        prot = self.protected()
        if any(p in prot for p in range(pos, end)) or (del_len == 0 and pos in prot and pos - 1 in prot):
            raise ValueError(f"edit at {pos} touches a splice site")
        new = []
        for e in self.exons:
            if e.end <= pos:                       # entirely before
                new.append(Exon(e.start, e.end, e.label, e.in_rna, e.inverted))
            elif e.start >= end:                   # entirely after
                new.append(Exon(e.start + delta, e.end + delta, e.label, e.in_rna, e.inverted))
            elif e.start <= pos and end <= e.end:  # inside this exon
                if del_len == 0 and (pos == e.start or pos == e.end):
                    raise ValueError(f"insertion at an exon boundary {pos}")
                new.append(Exon(e.start, e.end + delta, e.label, e.in_rna, e.inverted))
            else:
                raise ValueError(f"edit [{pos},{end}) straddles exon {e}")
        self.seq = self.seq[:pos] + ins + self.seq[end:]
        self.exons = new
        return self.check()

    # ---------------------------------------------------------------- io
    def to_json(self):
        return {"seq": self.seq, "exons": [e.to_json() for e in self.exons], "ops": self.ops, "provenance": self.provenance}

    @classmethod
    def from_json(cls, d):
        return cls(d["seq"], [Exon.from_json(e) for e in d["exons"]], d.get("ops"), d.get("provenance"))

    def save(self, path):
        with open(path, "w") as fh:
            json.dump(self.to_json(), fh)

    @classmethod
    def load(cls, path):
        with open(path) as fh:
            return cls.from_json(json.load(fh))

    def summary(self):
        return (f"{len(self.seq)} bp, {len(self.exons)} exons ({sum(1 for e in self.exons if e.in_rna)} in RNA), "
                f"chain {len(self.chain_seq())} bp, canonical={self.canonical()}, ops={[o['op'] for o in self.ops]}")


def random_seq(n, rng, gc=0.41):
    """Random sequence with the given GC fraction."""
    w = [(1 - gc) / 2, gc / 2, gc / 2, (1 - gc) / 2]
    return "".join(rng.choices(BASES, weights=w, k=n))


def synthetic_model(exon_lengths, intron_lengths, rng, gc=0.41):
    """A random gene with canonical GT…AG introns: len(intron_lengths) == len(exon_lengths) - 1."""
    if len(intron_lengths) != len(exon_lengths) - 1:
        raise ValueError("need one intron length per exon pair")
    seq, exons, pos = [], [], 0
    for i, L in enumerate(exon_lengths):
        exons.append(Exon(pos, pos + L, f"e{i + 1}"))
        seq.append(random_seq(L, rng, gc)); pos += L
        if i < len(intron_lengths):
            I = intron_lengths[i]
            if I < 30:
                raise ValueError("introns shorter than 30 bp are not simulated")
            seq.append("GT" + random_seq(I - 4, rng, gc) + "AG"); pos += I
    m = GeneModel("".join(seq), exons, provenance={"source": "synthetic", "exons": list(exon_lengths), "introns": list(intron_lengths)})
    return m


def seeded(seed, *parts):
    """A random.Random seeded stably from an integer seed and string parts (never hash())."""
    import zlib
    return random.Random((int(seed) * 1_000_003) ^ zlib.crc32("\x1f".join(str(p) for p in parts).encode()))
