"""Breakpoint groups in real RNA per family (BPCLUST, LEDGER_O3 prereg 156ad0a7).
usage: bpclust.py SPLIT OUTDIR      SPLIT = dev | heldout
Events per primary record overlapping an autosomal catalog copy: SPLIT (SA tag), CLIP (>= 100 bp soft clip, no SA), JUNCTION
(unannotated, non-canonical motif). Groups = same type, junction within 10 bp (+ partner within 10 bp / same side), >= 2 reads;
kept if 200 bp of read across the junction does not align to the reference (cov >= 0.95, identity >= 0.95)."""
import collections, os, pickle, re, subprocess, sys
import pysam

BAM = '/mnt/linuxdisk/home/juanfraitu/fibroblasts/GCA_029281585.2_flnc_mm.bam'
GENOME = '/mnt/linuxdisk/home/juanfraitu/_from_wsl/winloci_scratch/GGO.fasta'
GFF = '/mnt/linuxdisk/home/juanfraitu/winloci_data/GGO_genomic.gff'
COPIES = '/mnt/linuxdisk/home/juanfraitu/winloci_data/GGO_gwcat.copies.tsv'
TRUTH = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'bpclust_truth.tsv')
SEX = {'NC_073247.2', 'NC_073248.2', 'NC_011120.1'}
WIN, MIN_CLIP, MIN_READS, FLANK = 10, 100, 2, 100
CANON = {('GT', 'AG'), ('GC', 'AG'), ('AT', 'AC'), ('CT', 'AC'), ('CT', 'GC'), ('GT', 'AT')}


def gff_introns(cache):
    if os.path.exists(cache):
        return pickle.load(open(cache, 'rb'))
    ex = collections.defaultdict(list)
    for l in open(GFF):
        if l.startswith('#'):
            continue
        f = l.split('\t')
        if len(f) < 9 or f[2] != 'exon':
            continue
        par = re.search(r'Parent=([^;\n]+)', f[8])
        if par:
            ex[(f[0], par.group(1))].append((int(f[3]) - 1, int(f[4])))
    introns = set()
    for (c, _), xs in ex.items():
        xs.sort()
        for (a, b), (s, e) in zip(xs, xs[1:]):
            introns.add((c, b, s))
    pickle.dump(introns, open(cache, 'wb'))
    return introns


def near_annotated(introns, c, s, e):
    """a non-canonical junction within WIN bp (both ends) of an annotated intron is aligner jitter (dev finding, LEDGER_O3)"""
    return any((c, s + a, e + b) in introns for a in range(-WIN, WIN + 1) for b in range(-WIN, WIN + 1))


def clips(rec):
    ct = rec.cigartuples
    return (ct[0][1] if ct[0][0] == 4 else 0), (ct[-1][1] if ct[-1][0] == 4 else 0)


def sa_entries(rec):
    out = []
    for e in rec.get_tag('SA').rstrip(';').split(';'):
        c, pos, strand, cig = e.split(',')[:4]
        ops = re.findall(r'(\d+)([MIDNSHP=X])', cig)
        lead = int(ops[0][0]) if ops[0][1] in 'SH' else 0
        trail = int(ops[-1][0]) if ops[-1][1] in 'SH' else 0
        alen = sum(int(n) for n, o in ops if o in 'MI=X')
        rlen = sum(int(n) for n, o in ops if o in 'MDN=X')
        qs = trail if strand == '-' else lead  # start on the original read
        out.append((c, int(pos) - 1, int(pos) - 1 + rlen, strand, qs, qs + alen))
    return out


def events(rec, genome, introns):
    """(type, junction_pos, partner, query_pos_in_alignment_orientation)"""
    q = rec.query_length
    lc, rc = clips(rec)
    if rec.has_tag('SA'):
        o_s, o_e = (rc, q - lc) if rec.is_reverse else (lc, q - rc)
        for c, s, e, strand, qs, qe in sa_entries(rec):
            after = qs >= (o_s + o_e) / 2  # supplementary part comes later in the read
            # later in the read = reference end for a forward primary, start for a reverse one
            at_end = after != rec.is_reverse
            jpos = rec.reference_end if at_end else rec.reference_start
            facing_start = after == (strand == '+')
            ppos = s if facing_start else e
            yield 'SPLIT', jpos, (c, ppos), (q - rc if at_end else lc)
    else:
        if lc >= MIN_CLIP:
            yield 'CLIP', rec.reference_start, 'L', lc
        if rc >= MIN_CLIP:
            yield 'CLIP', rec.reference_end, 'R', q - rc
    p, qp = rec.reference_start, 0
    for op, n in rec.cigartuples:
        if op == 3:
            motif = (genome.fetch(rec.reference_name, p, p + 2).upper(), genome.fetch(rec.reference_name, p + n - 2, p + n).upper())
            if motif not in CANON and not near_annotated(introns, rec.reference_name, p, p + n):
                yield 'JUNCTION', p, p + n, qp
        if op in (0, 2, 3, 7, 8):
            p += n
        if op in (0, 1, 4, 7, 8):
            qp += n


def group(evs):
    """single linkage on junction position, then on partner"""
    out = []
    by = collections.defaultdict(list)
    for ev in evs:
        t, j, partner, rname, chrom, qp, seq = ev
        side = partner[0] if t == 'SPLIT' else partner if t == 'CLIP' else None
        by[(t, chrom, side)].append(ev)
    for key, xs in by.items():
        xs.sort(key=lambda x: x[1])
        runs, cur = [], [xs[0]]
        for x in xs[1:]:
            if x[1] - cur[-1][1] <= WIN:
                cur.append(x)
            else:
                runs.append(cur); cur = [x]
        runs.append(cur)
        for run in runs:
            if key[0] == 'CLIP':
                out.append(run); continue
            pv = (lambda x: x[2][1]) if key[0] == 'SPLIT' else (lambda x: x[2])
            run.sort(key=pv)
            sub = [run[0]]
            for x in run[1:]:
                if pv(x) - pv(sub[-1]) <= WIN:
                    sub.append(x)
                else:
                    out.append(sub); sub = [x]
            out.append(sub)
    return [g for g in out if len({x[3] for x in g}) >= MIN_READS]


def main():
    split, outdir = sys.argv[1], sys.argv[2]
    os.makedirs(outdir, exist_ok=True)
    truth = {}
    for l in open(TRUTH):
        f = l.split()
        if f[0] != 'family' and f[2] == split:
            truth[f[0]] = f[1]
    spans = collections.defaultdict(list)
    for l in open(COPIES):
        f = l.split('\t')
        if f[0] in truth and f[3] not in SEX:
            spans[f[0]].append((f[3], int(f[4]), int(f[5])))
    introns = gff_introns(os.path.join(outdir, '..', 'gff_introns.pkl'))
    genome = pysam.FastaFile(GENOME)
    bam = pysam.AlignmentFile(BAM)
    groups = []
    nreads = collections.Counter()
    for fam in sorted(spans, key=lambda x: int(x[5:])):
        seen, evs = set(), []
        for c, s, e in spans[fam]:
            for rec in bam.fetch(c, s, e):
                if rec.is_secondary or rec.is_supplementary or rec.is_unmapped or rec.query_name in seen:
                    continue
                seen.add(rec.query_name)
                for t, j, partner, qp in events(rec, genome, introns):
                    if any(cc == c2 and s2 <= j <= e2 for cc in [rec.reference_name] for c2, s2, e2 in spans[fam]):
                        evs.append((t, j, partner, rec.query_name, rec.reference_name, qp, rec.query_sequence))
        nreads[fam] = len(seen)
        for g in group(evs) if evs else []:
            groups.append((fam, g))
    with open(os.path.join(outdir, 'junctions.fa'), 'w') as o:
        for k, (fam, g) in enumerate(groups):
            t, j, partner, rname, chrom, qp, seq = g[0]
            o.write(f'>g{k}\n{seq[max(0, qp - FLANK):qp + FLANK]}\n')
    paf = os.path.join(outdir, 'junctions.paf')
    subprocess.run(f'minimap2 -c -x sr -N 50 -t 4 {GENOME} {outdir}/junctions.fa > {paf} 2> {outdir}/junctions.log',
                   shell=True, check=True)
    present = set()
    for l in open(paf):
        f = l.split('\t')
        if (int(f[3]) - int(f[2])) / int(f[1]) >= 0.95 and int(f[9]) / int(f[10]) >= 0.95:
            present.add(f[0])
    with open(os.path.join(outdir, 'groups.tsv'), 'w') as o:
        o.write('family\ttruth\tgroup\ttype\tchrom\tjunction\tpartner\tn_reads\tref_present\n')
        for k, (fam, g) in enumerate(groups):
            t, j, partner, *_ = g[0]
            o.write(f'{fam}\t{truth[fam]}\tg{k}\t{t}\t{g[0][4]}\t{j}\t{partner}\t{len({x[3] for x in g})}\t{int(f"g{k}" in present)}\n')
    # stage-2 rule (frozen from the dev census, LEDGER_O3): SPLIT or CLIP, >= 3 reads, partner not the mitochondrion
    with open(os.path.join(outdir, 'flags.tsv'), 'w') as o:
        o.write('family\ttruth\tflag\n')
        flag = {fam for k, (fam, g) in enumerate(groups) if f'g{k}' not in present and g[0][0] in ('SPLIT', 'CLIP')
                and len({x[3] for x in g}) >= 3 and not (g[0][0] == 'SPLIT' and g[0][2][0] in SEX)}
        for fam in sorted(spans, key=lambda x: int(x[5:])):
            o.write(f'{fam}\t{truth[fam]}\t{int(fam in flag)}\n')
    with open(os.path.join(outdir, 'families.tsv'), 'w') as o:
        o.write('family\ttruth\tn_reads\t' + '\t'.join(f'{t}_kept' for t in ('SPLIT', 'CLIP', 'JUNCTION')) + '\n')
        kept = collections.Counter((fam, g[0][0]) for k, (fam, g) in enumerate(groups) if f'g{k}' not in present)
        for fam in sorted(spans, key=lambda x: int(x[5:])):
            o.write(f'{fam}\t{truth[fam]}\t{nreads[fam]}\t' + '\t'.join(str(kept[(fam, t)]) for t in ('SPLIT', 'CLIP', 'JUNCTION')) + '\n')


if __name__ == '__main__':
    main()
