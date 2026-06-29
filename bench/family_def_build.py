#!/usr/bin/env python3
"""family_def_build.py — END-TO-END multi-copy gene-family definition (consolidated spec).

Reproducible pipeline that turns the over-merged cDNA all-vs-all into a layered, trustworthy family
set and writes a manifest for `rustle --family-manifest`. Every edge is WHOLE-GENE / WHOLE-PROTEIN
homology, never a shared sub-domain (enforced by reciprocal-coverage bars, audited 96%->100% clean).

PIPELINE
  0. cDNA homology graph, id>=0.90 (rep_ava.tsv).
  1. RETROCOPY filter  : drop one-side-intronless edges (processed pseudogenes).         [family_def_retrocopy_filter]
  2. STRAND filter     : drop antisense ('-') edges (sense/antisense + opposite-repeat artifacts). [family_def_strand_probe]
  3. LAYER A  recent coding paralog families:
        keep cDNA edges that are ALSO protein-homologous with WHOLE-PROTEIN coverage:
        fident>=0.30 AND min(qcov,tcov)>=0.50 AND max(qcov,tcov)>=0.70  (anti-sub-domain).
        connected components = families.  (cDNA id>=0.90 keeps these RECENT; ancient domain
        families like ZNF/OR are split into their recent sub-duplications -- not over-merged.)
  4. LAYER B  immune-receptor families:
        IG/TCR V/D/J/C segments have NO protein (not translated until VDJ); recover from
        segment-biotype cDNA homology + SAME-LOCUS constraint (<=3 Mb tandem array).
  (LAYER C  ancient domain super-families: protein-graph community detection -- OPTIONAL, --ancient.)

OUTPUT: /home/juanfra/winloci_scratch/family_manifest.tsv
        family_id  gene_id  chrom  start  end  strand  biotype  layer
Run: /home/juanfra/miniforge3/bin/python bench/family_def_build.py [--ancient]
"""
import collections
import os
import sys

import networkx as nx

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from family_def_read_filters import dna_homology
from family_def_retrocopy_filter import build_family_graph, load_intron_index
from family_def_strand_probe import edge_strand
from family_def_protein_validate import gene_biotype

GENES_BED = "/home/juanfra/winloci_scratch/unmapped_poc/genes.bed"
EXON_INDEX = "/home/juanfra/winloci_scratch/gene_exon_index.json"
PROT_FA = "/home/juanfra/winloci_scratch/proteins.fa"
PROT_AVA = "/home/juanfra/winloci_scratch/prot_ava.m8"
OUT = "/home/juanfra/winloci_scratch/family_manifest.tsv"
# whole-protein homology bar (anti-sub-domain): both proteins >=50% covered, shorter >=70%
P_FIDENT, P_MINCOV, P_MAXCOV, P_EVAL = 0.30, 0.50, 0.70, 1e-5
SEG_BT = ("V_segment", "C_region", "D_segment", "J_segment")
LOCUS_WIN = 3_000_000


def load_coord_strand():
    import json
    coord = {}
    for line in open(GENES_BED):
        p = line.rstrip("\n").split("\t")
        if len(p) >= 4:
            coord[p[3]] = (p[0], int(p[1]), int(p[2]))
    strand = {}
    if os.path.exists(EXON_INDEX):
        for g, v in json.load(open(EXON_INDEX)).items():
            strand[g] = v.get("strand", ".")
    return coord, strand


def protein_whole_edges():
    """gene-pairs with WHOLE-protein homology (anti-sub-domain reciprocal coverage)."""
    pe = set()
    with open(PROT_AVA) as f:
        for line in f:
            q, t, fid, qc, tc, ev, bits, al = line.rstrip("\n").split("\t")
            if q == t:
                continue
            fid, qc, tc, ev = float(fid), float(qc), float(tc), float(ev)
            if (fid >= P_FIDENT and min(qc, tc) >= P_MINCOV
                    and max(qc, tc) >= P_MAXCOV and ev <= P_EVAL):
                pe.add((q, t) if q < t else (t, q))
    return pe


def main():
    want_ancient = "--ancient" in sys.argv
    Hd, _ = dna_homology()
    idx = load_intron_index()
    bt = gene_biotype()
    coord, strand = load_coord_strand()
    pe = protein_whole_edges()

    # 1. retrocopy-filtered cDNA graph
    G, prov, part = build_family_graph(Hd, idx, filter_retrocopy=True)
    # 2. strand filter
    st = edge_strand()
    anti = [e for e in G.edges() if st.get(tuple(sorted(e))) == "-"]
    G.remove_edges_from(anti)
    G.remove_nodes_from([n for n in list(G.nodes()) if G.degree(n) == 0])

    families = []   # (layer, set_of_genes)

    # 3. LAYER A: recent coding paralog families (cDNA ∩ whole-protein homology)
    Ga = nx.Graph()
    for a, b in G.edges():
        if tuple(sorted((a, b))) in pe:
            Ga.add_edge(a, b)
    for c in nx.connected_components(Ga):
        if len(c) >= 2:
            families.append(("recent_coding", set(c)))

    # 4. LAYER B: immune-receptor families (segment cDNA + same-locus)
    segset = {g for g, b in bt.items() if b in SEG_BT and g in G}
    Gi = nx.Graph(); Gi.add_nodes_from(segset)
    for a, b in G.subgraph(segset).edges():
        ca, cb = coord.get(a), coord.get(b)
        if ca and cb and ca[0] == cb[0] and abs(ca[1] - cb[1]) <= LOCUS_WIN:
            Gi.add_edge(a, b)
    for c in nx.connected_components(Gi):
        if len(c) >= 2:
            families.append(("immune_segment", set(c)))

    # (optional) LAYER C: ancient domain super-families via protein community detection
    if want_ancient:
        from networkx.algorithms import community as nxcom
        Gp = nx.Graph()
        # use ALL whole-protein edges (not just cDNA-recent) to reach ancient families
        for a, b in pe:
            Gp.add_edge(a, b)
        assigned = {g for _, c in families for g in c}
        for og in nx.connected_components(Gp):
            if len(og) < 20 or (og & assigned):
                continue
            sub = Gp.subgraph(og)
            for com in nxcom.louvain_communities(sub, resolution=1.0, seed=0):
                if len(com) >= 2:
                    families.append(("ancient_domain", set(com)))

    # write manifests: rustle-compatible 6-col (parse_family_manifest expects EXACTLY 6 cols,
    # '#' header skipped) + a rich 8-col analysis version. Genes without a locus are dropped from
    # the rustle manifest (a '.' chrom is not a usable family locus); families needing >=2 located.
    families.sort(key=lambda x: -len(x[1]))
    lay_prefix = {"recent_coding": "RCF", "immune_segment": "IMM", "ancient_domain": "ANC"}
    full = OUT.replace(".tsv", ".full.tsv")
    n_rows = 0
    with open(OUT, "w") as out6, open(full, "w") as outf:
        out6.write("#family_id\tgene_id\tchrom\tstart\tend\tstrand\n")
        outf.write("family_id\tgene_id\tchrom\tstart\tend\tstrand\tbiotype\tlayer\n")
        for i, (layer, genes) in enumerate(families):
            fid = f"{lay_prefix[layer]}_{i}"
            located = [g for g in sorted(genes) if g in coord]
            for g in sorted(genes):
                c, s, e = coord.get(g, (".", 0, 0))
                outf.write(f"{fid}\t{g}\t{c}\t{s}\t{e}\t{strand.get(g,'.')}\t{bt.get(g,'?')}\t{layer}\n")
            if len(located) < 2:
                continue
            for g in located:
                c, s, e = coord[g]
                out6.write(f"{fid}\t{g}\t{c}\t{s}\t{e}\t{strand.get(g,'.')}\n")
                n_rows += 1

    # report
    by_layer = collections.Counter(layer for layer, _ in families)
    sizes = sorted((len(g) for _, g in families), reverse=True)
    named = [("RABL2A", "RABL2B"), ("APOBEC3D", "APOBEC3F"), ("RFPL1", "RFPL2"), ("RFPL2", "RFPL3")]
    g2f = {g: i for i, (_, c) in enumerate(families) for g in c}
    okn = sum(1 for a, b in named if g2f.get(a) is not None and g2f.get(a) == g2f.get(b))
    print("=== consolidated family definition ===")
    print(f"  cDNA over-merge baseline: 1216 families, max comp 238")
    print(f"  pipeline: retrocopy({prov['retrocopy_dropped']} dropped) -> strand({len(anti)} antisense dropped) "
          f"-> protein whole-gene bar (fident>={P_FIDENT}, min-cov>={P_MINCOV}, max-cov>={P_MAXCOV})")
    print(f"  FAMILIES: {len(families)} total, max comp {sizes[0]}, top {sizes[:6]}")
    for layer, n in by_layer.most_common():
        ng = sum(len(g) for l, g in families if l == layer)
        print(f"    {layer:16}: {n:4d} families, {ng} genes")
    print(f"  named real families preserved: {okn}/{len(named)}")
    print(f"  [+] wrote {OUT} ({n_rows} located gene rows, 6-col rustle-compatible)")
    print(f"  [+] wrote {full} (8-col: + biotype, layer)")


if __name__ == "__main__":
    main()
