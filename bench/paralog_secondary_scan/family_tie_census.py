#!/usr/bin/env python3
"""Family-resolved tie census: how common are TRULY-TIED multi-mapping reads, per paralog family?

A multimapper is FULLY TIED if its best-two placements are byte-identical on every discriminating
minimap2 tag (AS, NM, de, ms) -> impossible to assign (the DAZ-39 situation). Families are defined
DATA-DRIVEN: genes (from the RefSeq GFF) linked by shared multimapping reads = paralog copies
(connected components). Per family we report the multimapper count, the fully-tied fraction, and the
gene symbols at its loci. Answers: is the high DAZ tie rate a DAZ quirk or a property of its class?

Usage: python3 family_tie_census.py <contig> [<contig> ...]
"""
import sys, re, subprocess, collections

BAM = "/mnt/c/Users/jfris/Desktop/GGO.bam"
GFF = "/mnt/c/Users/jfris/Desktop/GGO_genomic.gff"
SAMTOOLS = "samtools"

def load_genes(contig):
    """Return sorted list of (start, end, symbol) genes on the contig."""
    genes = []
    for line in open(GFF):
        if line.startswith("#"):
            continue
        f = line.rstrip("\n").split("\t")
        if len(f) < 9 or f[2] != "gene" or f[0] != contig:
            continue
        m = re.search(r'gene=([^;]+)', f[8])
        genes.append((int(f[3]), int(f[4]), m.group(1) if m else "?"))
    genes.sort()
    return genes

def gene_at(genes, pos):
    """Symbol of the gene overlapping pos (1-based), or None. Linear scan ok for one contig."""
    for s, e, sym in genes:
        if s <= pos <= e:
            return sym
        if s > pos:
            break
    return None

def aln_tags(fields):
    """Extract (AS, NM, de, ms) from a SAM line's optional fields; missing -> sentinel."""
    AS = NM = ms = None
    de = None
    for t in fields[11:]:
        if t.startswith("AS:i:"): AS = int(t[5:])
        elif t.startswith("NM:i:"): NM = int(t[5:])
        elif t.startswith("ms:i:"): ms = int(t[5:])
        elif t.startswith("de:f:"): de = t[5:]   # keep string for byte-identity
    return AS, NM, de, ms

def census(contigs):
    for contig in contigs:
        genes = load_genes(contig)
        # Collect placements per read: read -> list of (pos, AS, NM, de, ms, gene)
        reads = collections.defaultdict(list)
        proc = subprocess.Popen([SAMTOOLS, "view", BAM, contig],
                                stdout=subprocess.PIPE, text=True, bufsize=1 << 20)
        n_aln = 0
        for line in proc.stdout:
            f = line.rstrip("\n").split("\t")
            flag = int(f[1])
            if flag & 0x800:   # skip supplementary (chimeric pieces, not separate placements)
                continue
            if flag & 0x4:     # unmapped
                continue
            n_aln += 1
            pos = int(f[3])
            AS, NM, de, ms = aln_tags(f)
            if AS is None:
                continue
            reads[f[0]].append((pos, AS, NM, de, ms))
        proc.wait()

        # Multimappers = reads with >=2 placements at DISTINCT loci (>50kb apart = different copy).
        # For each, best-two by AS; tied if identical on (AS,NM,de,ms). Link the genes of all placements.
        gene_mm = collections.Counter()    # gene-pair edges weighted by shared multimappers
        gene_link = collections.defaultdict(collections.Counter)
        per_read = []  # (genes_set, tied_bool)
        for rn, places in reads.items():
            if len(places) < 2:
                continue
            # distinct loci
            places_sorted = sorted(places, key=lambda p: -p[1])  # by AS desc
            # keep placements at distinct loci (>50kb from any kept)
            kept = []
            for p in places_sorted:
                if all(abs(p[0] - k[0]) > 50000 for k in kept):
                    kept.append(p)
            if len(kept) < 2:
                continue
            b0, b1 = kept[0], kept[1]
            tied = (b0[1] == b1[1] and b0[2] == b1[2] and b0[3] == b1[3] and b0[4] == b1[4])
            gset = set()
            for p in kept:
                g = gene_at(genes, p[0])
                gset.add(g if g else f"intergenic:{p[0]//100000*100000//1000}kb")
            per_read.append((frozenset(gset), tied))
            gl = sorted(gset)
            for i in range(len(gl)):
                for j in range(i + 1, len(gl)):
                    gene_link[gl[i]][gl[j]] += 1
                    gene_link[gl[j]][gl[i]] += 1

        # Connected components over the gene-link graph (>=2 shared multimappers = an edge).
        adj = collections.defaultdict(set)
        for a, nb in gene_link.items():
            for b, w in nb.items():
                if w >= 2:
                    adj[a].add(b); adj[b].add(a)
        seen = set(); comps = []
        for node in adj:
            if node in seen: continue
            stack = [node]; comp = set()
            while stack:
                x = stack.pop()
                if x in seen: continue
                seen.add(x); comp.add(x)
                stack.extend(adj[x] - seen)
            comps.append(comp)
        # Map each multimapper to a family (the component containing >=2 of its genes)
        node2comp = {}
        for ci, comp in enumerate(comps):
            for g in comp:
                node2comp[g] = ci
        fam_total = collections.Counter(); fam_tied = collections.Counter()
        fam_genes = collections.defaultdict(set)
        n_mm = 0; n_tied = 0
        for gset, tied in per_read:
            n_mm += 1; n_tied += tied
            cs = set(node2comp[g] for g in gset if g in node2comp)
            if len(cs) == 1:
                ci = next(iter(cs))
                fam_total[ci] += 1; fam_tied[ci] += tied
                fam_genes[ci].update(g for g in gset if not g.startswith("intergenic"))

        print("=" * 78)
        print(f"CONTIG {contig}: {n_aln:,} primary/secondary alns; "
              f"{n_mm:,} multimappers (>=2 distinct loci); "
              f"{n_tied:,} fully-tied ({100*n_tied/max(1,n_mm):.1f}%)")
        print(f"  {len(comps)} multimapper-linked families")
        print(f"  {'family (gene symbols)':42s} {'n_mm':>6s} {'tied':>6s} {'tie%':>6s}")
        for ci, tot in fam_total.most_common(25):
            syms = sorted(s for s in fam_genes[ci])[:4]
            label = ",".join(syms) if syms else "(intergenic)"
            print(f"  {label:42.42s} {tot:6d} {fam_tied[ci]:6d} {100*fam_tied[ci]/tot:5.1f}%")

if __name__ == "__main__":
    census(sys.argv[1:] if len(sys.argv) > 1 else ["NC_086018.1"])
