#!/usr/bin/env python3
"""PREREG ISOCON-GORILLA (2026-10-31): targets, read pools and IsoCon candidate GTFs for the IsoCon vs our-arm comparison.
  isocon_targets.py select  OUR.gtf BAM OUT_DIR          -> OUT_DIR/targets.tsv (target, class, gene, name, chrom, start0, end, strand, records)
  isocon_targets.py pool    BAM OUT_DIR                  -> OUT_DIR/pools/<target>.fa (every primary/secondary record on any span, by read name)
  isocon_targets.py gtf     GENOME OUT_DIR               -> OUT_DIR/isocon.gtf from each target's IsoCon final_candidates.fa (minimap2 splice:hq)
Targets: the 4 Yoo families (famtruth.GGO_MAINSTAYS, class multi) + 30 single-copy genes (class single): protein-coding, no Compara paralog
at any level, a de novo locus in OUR.gtf, 50-5,000 records on the span, random.Random(0) over the eligible ids sorted."""
import sys, os, re, collections, subprocess, random
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import famtruth as F
cmd = sys.argv[1]
def sh(c): return subprocess.run(c, shell=True, capture_output=True, text=True).stdout
RC = str.maketrans('ACGTacgtNn', 'TGCAtgcaNn')
if cmd == 'select':
    gtf, bam, out = sys.argv[2:5]; os.makedirs(out, exist_ok=True)
    G = F.load_ggo(F.load_cat()); idx = F.exon_index(G)
    rows = []
    for name, pat in F.GGO_MAINSTAYS.items():
        rx = re.compile(pat)
        for g, d in G.items():
            if rx.search(d['name']) or rx.search(d.get('desc', '')):
                rows.append((name, 'multi', g, d))
    levels = set()
    for v in F.LEVELS.values(): levels |= set(v)
    pairs, universe = F.truth_pairs(G, F.PARALOGS, levels)
    para = {g for p in pairs for g in p}
    expressed = F.expressed_from_gtf(gtf, idx)
    elig = sorted(g for g, d in G.items() if d['biotype'] == 'protein_coding' and g in universe and g not in para and g in expressed and d['ens'])
    rnd = random.Random(0); rnd.shuffle(elig)
    picked = []
    for g in elig:
        d = G[g]; s, e = d['exons'][0][0], d['exons'][-1][1]
        n = int(sh(f"samtools view -c {bam} {d['chrom']}:{s + 1}-{e}").strip() or 0)
        if 50 <= n <= 5000: picked.append((g, d, n))
        if len(picked) == 30: break
    with open(f"{out}/targets.tsv", 'w') as w:
        w.write("target\tclass\tgene\tname\tchrom\tstart0\tend\tstrand\trecords\n")
        for name, cls, g, d in rows:
            s, e = d['exons'][0][0], d['exons'][-1][1]
            n = sh(f"samtools view -c {bam} {d['chrom']}:{s + 1}-{e}").strip()
            w.write(f"{name}\t{cls}\t{g}\t{d['name']}\t{d['chrom']}\t{s}\t{e}\t{d['strand']}\t{n}\n")
        for g, d, n in picked:
            s, e = d['exons'][0][0], d['exons'][-1][1]
            w.write(f"{d['name']}\tsingle\t{g}\t{d['name']}\t{d['chrom']}\t{s}\t{e}\t{d['strand']}\t{n}\n")
    print(f"targets: {len(rows)} multi-copy members in {len(F.GGO_MAINSTAYS)} families, {len(picked)} single-copy genes (eligible {len(elig)})")
elif cmd == 'pool':
    bam, out = sys.argv[2:4]; os.makedirs(f"{out}/pools", exist_ok=True)
    spans = collections.defaultdict(list)
    for l in list(open(f"{out}/targets.tsv"))[1:]:
        t, cls, g, name, chrom, s, e, st, n = l.rstrip('\n').split('\t'); spans[t].append(f"{chrom}:{int(s) + 1}-{e}")
    for t, regs in spans.items():
        seen = {}
        for line in sh(f"samtools view -F 2052 {bam} {' '.join(regs)}").splitlines():
            f = line.split('\t', 11)
            if f[0] in seen or f[9] == '*': continue
            seq = f[9].translate(RC)[::-1] if int(f[1]) & 16 else f[9]
            seen[f[0]] = seq
        with open(f"{out}/pools/{t}.fa", 'w') as w:
            for k, v in seen.items(): w.write(f">{k}\n{v}\n")
        print(f"{t}\t{len(seen)} reads")
elif cmd == 'gtf':
    genome, out = sys.argv[2:4]
    with open(f"{out}/isocon.gtf", 'w') as w:
        for t in sorted(os.listdir(f"{out}/isocon")):
            fa = f"{out}/isocon/{t}/final_candidates.fa"
            if not os.path.exists(fa) or os.path.getsize(fa) == 0: print(f"{t}\tno candidates"); continue
            sam = sh(f"minimap2 -ax splice:hq -uf --eqx -N 50 -p 0.1 --secondary=yes -t 8 {genome} {fa} 2>/dev/null")
            n = 0
            for line in sam.splitlines():
                if line.startswith('@'): continue
                f = line.split('\t', 11)
                flag = int(f[1])
                if flag & 2308 or f[2] == '*': continue
                pos = int(f[3]); exons = []; cur = pos
                for ln, op in re.findall(r'(\d+)([MIDNSHP=X])', f[5]):
                    ln = int(ln)
                    if op in 'M=XD':
                        if not exons or exons[-1][1] != cur: exons.append([cur, cur + ln])
                        else: exons[-1][1] = cur + ln
                        cur += ln
                    elif op == 'N': cur += ln
                exons = [(a, b - 1) for a, b in exons]
                strand = '-' if flag & 16 else '+'
                tid = f"{t}|{f[0]}"; n += 1
                w.write(f"{f[2]}\tIsoCon\ttranscript\t{exons[0][0]}\t{exons[-1][1]}\t.\t{strand}\t.\tgene_id \"{t}\"; transcript_id \"{tid}\";\n")
                for i, (a, b) in enumerate(exons): w.write(f"{f[2]}\tIsoCon\texon\t{a}\t{b}\t.\t{strand}\t.\tgene_id \"{t}\"; transcript_id \"{tid}\"; exon_number \"{i + 1}\";\n")
            print(f"{t}\t{n} candidates placed")
