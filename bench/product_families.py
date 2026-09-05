#!/usr/bin/env python3
"""Multi-copy family TRUTH from annotation PRODUCT strings — independent of the sequence clustering.

⚠ Reads CDS products, pseudogenes AND lncRNAs. A truncated index holding only mRNA-level products once
produced a retracted answer key (project-gene-naming-traps): "symbol-only intervals = 0/26, not 7/26".

⭐ This is a TRUTH PROXY, not the definition. The definition is MCL over sequence. Products come from a
different source (curator descriptions) than the clustering (sequence), so recovery is a real test.
⚠ It is still ANNOTATION, so it inherits annotation's blind spots: a family whose members are all
"uncharacterized" is invisible here, which is exactly the NPIP situation and why NPIP needs its own
coordinate truth.
"""
import re, sys, collections

GFF = '/mnt/linuxdisk/home/juanfraitu/winloci_data/GGO_genomic.gff'
STOP = re.compile(r'\b(low quality protein|isoform|variant|partial|putative|uncharacterized|'
                  r'hypothetical|protein|like|pseudogene|transcript)\b', re.I)

def norm(p):
    p = p.split(' isoform')[0].split(', transcript variant')[0]
    p = re.sub(r'\bLOW QUALITY PROTEIN:\s*', '', p, flags=re.I)
    p = re.sub(r'[-_]?like\b', '', p, flags=re.I)
    # ⚠ DO NOT strip a trailing number: "zinc finger protein 700" and "...763" are DIFFERENT genes.
    # Stripping it collapsed 553 ZNFs into one "family" — the superfamily artefact MCL exists to avoid
    # (NEGATIVE_RESULTS_REGISTER.md:474: transitive closure gave 145- and 114-gene superfamilies).
    p = p.strip()
    return re.sub(r'\s+', ' ', p).strip().lower()

def load():
    gene = {}                      # gene_id -> (chrom,start,end,name,biotype)
    prod = collections.defaultdict(set)   # gene_id -> products seen anywhere under it
    tx2gene = {}
    for line in open(GFF):
        if line.startswith('#'): continue
        f = line.rstrip('\n').split('\t')
        if len(f) < 9: continue
        a = f[8]
        gid = re.search(r'\bgene=([^;]+)', a)
        if f[2] in ('gene', 'pseudogene'):
            i = re.search(r'ID=([^;]+)', a)
            n = re.search(r'Name=([^;]+)', a)
            bt = re.search(r'gene_biotype=([^;]+)', a)
            if i: gene[n.group(1) if n else i.group(1)] = (f[0], int(f[3])-1, int(f[4]),
                                                           n.group(1) if n else '?',
                                                           bt.group(1) if bt else '?')
        pm = re.search(r'\bproduct=([^;]+)', a)
        if pm and gid: prod[gid.group(1)].add(pm.group(1))
        dm = re.search(r'\bdescription=([^;]+)', a)
        if dm and gid: prod[gid.group(1)].add(dm.group(1))
    return gene, prod

if __name__ == '__main__':
    gene, prod = load()
    fam = collections.defaultdict(list)
    for g, (c, s, e, n, bt) in gene.items():
        ps = prod.get(g, set())
        if not ps: continue
        key = norm(sorted(ps, key=len)[0])
        if not key or len(key) < 8: continue
        if key in ('uncharacterized protein','uncharacterized','protein'): continue
        # ncRNA arrays (U6 snRNA n=1315, 5S rRNA, snoRNA...) are repeat arrays, not gene families.
        if bt not in ('protein_coding','pseudogene'): continue
        fam[key].append((c, s, e, n, bt))
    multi = {k: v for k, v in fam.items() if len(v) >= 3}
    print(f"genes with a product/description: {sum(1 for g in gene if prod.get(g))}/{len(gene)}")
    print(f"product families with >=3 members: {len(multi)}   members: {sum(len(v) for v in multi.values())}")
    disp = {k: v for k, v in multi.items() if len({x[0] for x in v}) > 1}
    print(f"  of which DISPERSED across >1 contig: {len(disp)}")
    print("\nlargest product families (size, %LOC-named, contigs):")
    for k, v in sorted(multi.items(), key=lambda x: -len(x[1]))[:15]:
        loc = sum(1 for x in v if x[3].startswith('LOC'))
        print(f"  {len(v):>4}  LOC {loc/len(v):>5.0%}  contigs={len({x[0] for x in v}):>2}  {k[:62]}")
    import json
    json.dump({k: [list(x) for x in v] for k, v in multi.items()},
              open('/mnt/linuxdisk/home/juanfraitu/mcl_ann/product_families.json','w'))
    print("\nwritten product_families.json")
