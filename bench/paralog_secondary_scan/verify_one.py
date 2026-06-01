#!/usr/bin/env python3
"""Independent per-hit re-derivation for a verification agent.
Given a gene_id, query the gene locus AND the sibling locus directly from the
BAM, join reads by name, and recompute the raw-dNM copy anchor from scratch
(independent of the genome-wide pipeline). Also emit artifact signals.

Usage: verify_one.py GENE_ID [BAM]
"""
import sys, os, re, json, subprocess, collections, statistics as st

T = 2  # distinguishing-mismatch threshold

def parse_cig(c):
    return re.findall(r'(\d+)([MIDNSH=X])', c)

def junctions(pos, cig):
    cur = pos; js = []
    for ln, op in parse_cig(cig):
        ln = int(ln)
        if op == 'N':
            js.append((cur - 1, cur + ln)); cur += ln
        elif op in 'MD=X':
            cur += ln
    return js

def alen(c):
    return sum(int(l) for l, o in parse_cig(c) if o in 'M=XD') or 1

def softclip_frac(c):
    sc = sum(int(l) for l, o in parse_cig(c) if o in 'SH')
    q = sum(int(l) for l, o in parse_cig(c) if o in 'M=XIS H'.replace(' ', ''))
    return sc / max(1, q)

def query(bam, chrom, start, end, introns, strand, K):
    """qname -> best chain-matching alignment dict, for reads on `strand`."""
    if not chrom:
        return {}
    out = subprocess.run(['samtools', 'view', bam, f"{chrom}:{start}-{end}"],
                         capture_output=True, text=True, timeout=180).stdout
    best = {}
    for line in out.splitlines():
        f = line.split('\t')
        if len(f) < 11:
            continue
        flag = int(f[1])
        if flag & 0x800:
            continue
        st_ = '-' if flag & 16 else '+'
        if st_ != strand:
            continue
        js = junctions(int(f[3]), f[5])
        nmatch = sum(1 for j in js if j in introns)
        if nmatch < K:
            continue
        nm = next((int(t[5:]) for t in f[11:] if t.startswith('NM:i:')), None)
        if nm is None:
            continue
        d = dict(nm=nm, alen=alen(f[5]), mapq=int(f[4]),
                 sc=softclip_frac(f[5]), spliced=('N' in f[5]), nmatch=nmatch)
        b = best.get(f[0])
        if b is None or nm < b['nm']:
            best[f[0]] = d
    return best

def main():
    gid = sys.argv[1]
    bam = sys.argv[2] if len(sys.argv) > 2 else '../GGO.bam'
    P = json.load(open('/tmp/verify_payloads.json'))[gid]
    iG = set(tuple(map(int, p.split(':'))) for p in P['chain'].split(',') if p)
    iS = set(tuple(map(int, p.split(':'))) for p in P['sibling_chain'].split(',')
             if p) if P.get('sibling_chain') else set()
    nintr = int(P['n_introns'])
    K = max(2, min(8, (nintr + 2) // 3))
    G = query(bam, P['chrom'], P['start'], P['end'], iG, P['strand'], K)
    Ssib = {}
    if P.get('sibling_chrom') and iS:
        Ssib = query(bam, P['sibling_chrom'], P['sibling_start'],
                     P['sibling_end'], iS, P['sibling_strand'], max(2, K))
    ga = sa = amb = 0
    for qn, g in G.items():
        s = Ssib.get(qn)
        if s is None or s['alen'] < 0.7 * g['alen']:
            ga += 1; continue            # uniquely G (no comparable sibling aln)
        dnm = s['nm'] - g['nm']
        if dnm >= T:
            ga += 1
        elif dnm <= -T:
            sa += 1
        else:
            amb += 1
    mapqs = [g['mapq'] for g in G.values()]
    rep = dict(
        gene=gid, gene_name=P['gene_name'], gene_desc=P['gene_desc'][:60],
        nm_taxonomy=P['nm_taxonomy'], sibling=P['sibling'],
        sibling_name=P['sibling_name'], locus_rel=P['locus_rel'],
        n_chain_reads_at_gene=len(G), reanchor_G=ga, reanchor_sibling=sa,
        reanchor_ambiguous=amb,
        sibling_chain_reads=len(Ssib),
        median_mapq=st.median(mapqs) if mapqs else None,
        frac_softclipped_gt30=round(
            sum(1 for g in G.values() if g['sc'] > 0.3) / max(1, len(G)), 2),
        frac_spliced=round(
            sum(1 for g in G.values() if g['spliced']) / max(1, len(G)), 2),
        n_introns=nintr,
        pipeline_anchor=dict(G=P['nm_anchor_G'], sib=P['nm_anchor_sibling'],
                             amb=P['nm_anchor_ambiguous']))
    print(json.dumps(rep, indent=1))

if __name__ == '__main__':
    main()
