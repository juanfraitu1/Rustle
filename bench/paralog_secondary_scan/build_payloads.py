#!/usr/bin/env python3
"""Build per-hit verification payloads for the workflow agents.
Emits /tmp/verify_payloads.json : {gene_id -> payload} for the
expressed_real_copy + pure_spillover hits, and prints the gene-id lists.
"""
import json, sys

hits = json.load(open('/tmp/hits_nm.json'))
genes = {}
for line in open('/tmp/gene_introns.tsv'):
    if line.startswith('gene_id'):
        continue
    gid, chrom, s, e, strand, ntx, chain = line.rstrip('\n').split('\t')
    genes[gid] = dict(chrom=chrom, start=int(s), end=int(e), strand=strand,
                      chain=chain)

KEEP = {'expressed_real_copy', 'pure_spillover'}
payloads = {}
ids = {'expressed_real_copy': [], 'pure_spillover': []}
for h in hits:
    if h['nm_taxonomy'] not in KEEP:
        continue
    gid = h['gene']; g = genes[gid]
    sib = h.get('nm_dom_sibling') or h.get('dom_sibling')
    sg = genes.get(sib, {})
    payloads[gid] = dict(
        gene=gid, gene_name=h.get('gene_name', ''),
        gene_desc=h.get('gene_desc', ''), biotype=h.get('gene_biotype', ''),
        chrom=g['chrom'], start=g['start'], end=g['end'], strand=g['strand'],
        n_introns=h.get('n_introns'), K=h.get('K'), chain=g['chain'],
        klass=h['klass'], nm_taxonomy=h['nm_taxonomy'],
        nm_anchor_G=h['nm_anchor_G'], nm_anchor_sibling=h['nm_anchor_sibling'],
        nm_anchor_ambiguous=h['nm_anchor_ambiguous'],
        nm_margin_median=h.get('nm_margin_median'),
        n_chain_reads=h.get('n_chain_reads'),
        sibling=sib, sibling_name=h.get('sibling_name', ''),
        sibling_desc=h.get('sibling_desc', ''),
        sibling_chrom=sg.get('chrom'), sibling_start=sg.get('start'),
        sibling_end=sg.get('end'), sibling_strand=sg.get('strand'),
        sibling_chain=sg.get('chain', ''), locus_rel=h.get('locus_rel'))
    ids[h['nm_taxonomy']].append(gid)

json.dump(payloads, open('/tmp/verify_payloads.json', 'w'))
json.dump(ids, open('/tmp/verify_ids.json', 'w'))
print(f"payloads: {len(payloads)}  "
      f"expressed={len(ids['expressed_real_copy'])}  "
      f"spillover={len(ids['pure_spillover'])}", file=sys.stderr)
