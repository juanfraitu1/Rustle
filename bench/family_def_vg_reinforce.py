#!/usr/bin/env python3
"""
family_def_vg_reinforce.py — clean genome-wide multi-copy family definition:
de-novo loci (tight vertices) + per-family VG BACKBONE-REINFORCEMENT subgraph + guard.

  DISCOVER : read-conflict (de-tie) graph over de-novo loci -> candidate families.
  REINFORCE: for each family, build every member's all-isoform exon-union copy model,
             align them, and form the FAMILY BACKBONE GRAPH (edge = two members share a
             real backbone). Real copies mutually reinforce (a member backbone-linked to
             several mates is strongly real); a repeat-bridge member has NO backbone link
             and is rejected. Size-2 families survive (the pair reinforces itself iff it
             shares a backbone). A guard never rejects a member too sparse to validate.
  RESULT   : validated families = backbone-connected cores; bridge members rejected;
             over-merged families split into their real backbone cores.

OOM-safe: models built ONLY for family-member loci (~1.3k); reads capped; streamed.
Run (background): /home/juanfra/miniforge3/bin/python bench/family_def_vg_reinforce.py
"""
import collections
import json
import os
import subprocess
import sys

import pysam
import networkx as nx

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from family_def_genomewide import (scan as conflict_scan, pair_evidence, components,
                                    best_gene, DELTA, DE_MAX, MIN_READS)

BAM = "/home/juanfra/winloci_scratch/GGO.bam"
GENOME = "/home/juanfra/winloci_scratch/GGO.fasta"
META = "/home/juanfra/winloci_scratch/denovo_transcripts.meta.tsv"
GENES_BED = "/home/juanfra/winloci_scratch/unmapped_poc/genes.bed"
SEQS = "/home/juanfra/winloci_scratch/vg_reinforce_copies.fa"
PAF = "/home/juanfra/winloci_scratch/vg_reinforce_copies.paf"
PAD = 5000
MAX_READS = 3000
GUARD_READS = 20      # below this, a model is "unvalidatable" -> never rejected (guard)
MIN_LEN = 300
ID_MIN = 0.80
COV_MIN = 0.30        # backbone edge: id>=.80 AND max-coverage of the pair >= this


def load_denovo():
    by_chrom = collections.defaultdict(list)
    info = {}
    with open(META) as f:
        next(f)
        for line in f:
            p = line.rstrip("\n").split("\t")
            vid, c, s, e = p[0], p[1], int(p[2]), int(p[3])
            by_chrom[c].append((s, e, vid))
            info[vid] = (c, s, e, int(p[6]))
    for c in by_chrom:
        by_chrom[c].sort()
    return by_chrom, info


def merge(iv):
    iv.sort(); out = []
    for s, e in iv:
        if out and s <= out[-1][1]:
            out[-1][1] = max(out[-1][1], e)
        else:
            out.append([s, e])
    return out


def build_model(bam, genome, chrom, start, end):
    blocks = []; n = 0
    try:
        it = bam.fetch(chrom, max(0, start - PAD), end + PAD)
    except (ValueError, KeyError):
        return "", 0
    for r in it:
        if r.is_unmapped or r.is_supplementary or r.is_secondary:
            continue  # primary only: the new -N50 -p0.1 BAM has 65x more secondaries that would pollute S(v)
        if r.reference_end is None or r.reference_end < start or r.reference_start > end:
            continue
        b = r.get_blocks()
        if b:
            blocks.extend(b); n += 1
            if n >= MAX_READS:
                break
    if not blocks:
        return "", n
    union = merge([list(x) for x in blocks])
    return "".join(genome.fetch(chrom, s, e).upper() for s, e in union), n


def main():
    print("[1/4 discover] read-conflict over de-novo loci ...", flush=True)
    by_chrom, info = load_denovo()
    mm, _, _ = conflict_scan(by_chrom)
    ev = pair_evidence(mm)
    edges, fams = components(ev, DELTA, DE_MAX, MIN_READS)
    member_loci = sorted({g for c in fams for g in c})
    def chrom_of(v): return info[v][0] if v in info else None
    n_bridge0 = sum(1 for c in fams if len({chrom_of(g) for g in c}) >= 3)
    print(f"  candidate: {len(fams)} families, {len(edges)} edges, {len(member_loci)} loci; "
          f"cross-chrom bridges={n_bridge0}", flush=True)

    models_present = set()
    if os.path.exists(SEQS):
        print("[2/4 build] CACHED: reusing existing copy models (loci are deterministic)", flush=True)
        for line in open(SEQS):
            if line.startswith(">"):
                models_present.add(line[1:].strip())
    else:
        print(f"[2/4 build] exon-union copy models for {len(member_loci)} member loci ...", flush=True)
        bam = pysam.AlignmentFile(BAM, "rb"); genome = pysam.FastaFile(GENOME)
        with open(SEQS, "w") as out:
            for i, v in enumerate(member_loci):
                c, s, e, _ = info[v]
                seq, _ = build_model(bam, genome, c, s, e)
                if len(seq) >= MIN_LEN:
                    out.write(f">{v}\n{seq}\n"); models_present.add(v)
                if (i + 1) % 300 == 0:
                    print(f"    ... {i+1}/{len(member_loci)}", flush=True)
        bam.close(); genome.close()
    print(f"  {len(models_present)} copy models", flush=True)
    # guard: a de-novo locus is "validatable" (trustworthy enough to reject on) iff it has
    # a model AND >= GUARD_READS reads (from the meta n_reads, info[v][3]).
    nreads = {v: (info[v][3] if v in models_present else 0) for v in member_loci}

    print("[3/4 reinforce] backbone alignment + intersection graph (~R ∩ ~B) ...", flush=True)
    if not os.path.exists(PAF):
        with open(PAF, "w") as pf:
            subprocess.run(["minimap2", "-cx", "asm20", "-N", "30", "-p", "0.05", "-t", "6",
                            SEQS, SEQS], stdout=pf, stderr=subprocess.DEVNULL)
    H = {}
    with open(PAF) as f:
        for line in f:
            x = line.split("\t")
            if len(x) < 11:
                continue
            q, ql, qs, qe, t, m, bl = x[0], int(x[1]), int(x[2]), int(x[3]), x[5], int(x[9]), int(x[10])
            if q == t:
                continue
            ident = m / bl if bl else 0; cov = (qe - qs) / ql if ql else 0
            key = (q, t) if q < t else (t, q)
            r = H.setdefault(key, {"id": 0.0, "ca": 0.0, "cb": 0.0})
            r["id"] = max(r["id"], ident)
            if q == key[0]:
                r["ca"] = max(r["ca"], cov)
            else:
                r["cb"] = max(r["cb"], cov)

    def backbone(a, b):
        # ~B: the two copies share a backbone (reciprocal homology over a real fraction
        # of both). Repeat-bridges fail it (one side aligns only over a short element).
        r = H.get((a, b) if a < b else (b, a))
        return bool(r and r["id"] >= ID_MIN and min(r["ca"], r["cb"]) >= COV_MIN)

    # THE DEFINITION (two-stage): a candidate family is an ~R-component; the VALIDATED
    # family is its ~B-connected core. A member is kept iff it shares a backbone with
    # SOME family-mate (membership is defined w.r.t. the family, not a specific edge),
    # so a repeat-bridge member -- backbone-linked to no one -- falls out, while real
    # copies stay connected even where an individual conflict edge has a weak backbone.
    validated = []
    rejected_members = []
    for fam in fams:
        members = sorted(fam)
        BG = nx.Graph(); BG.add_nodes_from(members)
        for i in range(len(members)):
            for j in range(i + 1, len(members)):
                if backbone(members[i], members[j]):
                    BG.add_edge(members[i], members[j])
        for comp in nx.connected_components(BG):
            if len(comp) >= 2:
                validated.append(set(comp))
        for v in members:                    # isolated in ~B AND validatable -> bridge
            if BG.degree[v] == 0 and nreads.get(v, 0) >= GUARD_READS:
                rejected_members.append(v)

    n_bridge1 = sum(1 for c in validated if len({chrom_of(g) for g in c}) >= 3)

    # dump the VALIDATED families' member loci (id + coords) so the genome-wide PSV
    # resolution scan can rebuild each family's variation graph (psv_graph_genomewide.py).
    fam_dump = "/home/juanfra/winloci_scratch/validated_families.tsv"
    with open(fam_dump, "w") as o:
        o.write("family\tlocus\tchrom\tstart\tend\tnreads\n")
        for fi, fam in enumerate(sorted(validated, key=lambda c: -len(c))):
            for v in sorted(fam, key=lambda x: info[x][:3]):
                c, s, e, nr = info[v]
                o.write(f"{fi}\t{v}\t{c}\t{s}\t{e}\t{nr}\n")
    print(f"  [+] wrote {fam_dump} ({len(validated)} families)")

    print("\n[4/4 result]")
    print(f"  candidate families   : {len(fams)}   ({n_bridge0} cross-chrom bridges)")
    print(f"  VALIDATED families   : {len(validated)} (backbone-reinforced cores)   "
          f"({n_bridge1} cross-chrom bridges)")
    print(f"  rejected bridge members (validatable, 0 backbone links): {len(rejected_members)}")
    sizes0 = collections.Counter(len(c) for c in fams)
    sizes1 = collections.Counter(len(c) for c in validated)
    print(f"  size-2 families: {sizes0[2]} -> {sizes1[2]} (kept iff the pair shares a backbone)")

    # DNA-truth audit: map each rejected/kept ORIGINAL edge to genes, classify
    gby = collections.defaultdict(list)
    with open(GENES_BED) as f:
        for line in f:
            c, s, e, name = line.rstrip("\n").split("\t")
            gby[c].append((int(s), int(e), name))
    for c in gby:
        gby[c].sort()
    def gene_of(v):
        c, s, e, _ = info[v]
        return best_gene(gby, c, s, e)
    try:
        from family_def_read_filters import dna_homology
        Hd, dna = dna_homology()
        def cls_edge(a, b):
            ga, gb = gene_of(a), gene_of(b)
            if not ga or not gb or ga == gb:
                return "no_gene"
            r = Hd.get((ga, gb) if ga < gb else (gb, ga))
            if r is None or r["id"] == 0:
                return "no_homology"
            if r["id"] >= 0.90 and max(r["cov_a"], r["cov_b"]) >= 0.30:
                return "DNA_paralog"
            return "sub_bar" if r["id"] >= 0.80 else "weak"
        rej_set = set(rejected_members)
        rej_edges = [(a, b) for a, b, n in edges if a in rej_set or b in rej_set]
        rc = dict(collections.Counter(cls_edge(a, b) for a, b in rej_edges))
        print("\n  rejected-member edges by DNA class (want mostly no_homology):")
        for k, v in sorted(rc.items(), key=lambda kv: -kv[1]):
            print(f"     {v:5d}  {k}")
    except Exception as ex:
        rc = {"error": str(ex)}
        print(f"  (DNA audit skipped: {ex})")

    json.dump(dict(candidate_families=len(fams), validated_families=len(validated),
                   rejected_members=len(rejected_members),
                   bridges_before=n_bridge0, bridges_after=n_bridge1,
                   size2_before=sizes0[2], size2_after=sizes1[2],
                   rejected_edge_dna_classes=rc),
              open(os.path.join(HERE, "family_def_vg_reinforce.json"), "w"), indent=2)
    print(f"\n[+] wrote {os.path.join(HERE,'family_def_vg_reinforce.json')}")


if __name__ == "__main__":
    main()
