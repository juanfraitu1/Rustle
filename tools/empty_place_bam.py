#!/usr/bin/env python3
"""The alignments that lie where an assembly has no transcript on their strand (register 1935).

    empty_place_bam.py ASSEMBLY.gtf IN.bam OUT.bam [--threads N]

The multi-copy stage assembles a second time at a one-read floor and keeps only transcripts at places ASSEMBLY.gtf leaves
empty on that strand (multicopy_rescue.py `candidates`: no overlap with the merged transcript spans of the strand). An
alignment that overlaps such a span cannot support a kept transcript, so the second assembly does not need it. This
writes the others: every record (primary, secondary, supplementary) whose span overlaps no transcript span of ASSEMBLY.gtf
on the record's transcript strand (minimap2 ts tag, as tools/eq_classes.py), coordinate-sorted and indexed.
Same definition of an occupied place as multicopy_rescue.py (its `spans` / `occupied`).
RESULT (register 1935, OR6737): NOT equivalent. Assembling at a one-read floor on these alignments alone gives 63% more
candidates and 3.4 x the rescued transcripts, and loses 4% / 25% of the full run's: the alignments on occupied places are
the context that makes the assembler collapse or reject chains next to them. Kept as a recorded arm; no stage calls it."""
import sys, os, argparse, multiprocessing, subprocess
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import multicopy_rescue as M, eq_classes as E
G = {}

def work(job):
    import pysam
    i, contig = job; part = f"{G['out']}.part.{i:06d}.bam"; kept = seen = 0
    B = pysam.AlignmentFile(G['bam'])
    with pysam.AlignmentFile(part, 'wb', template=B) as w:
        for r in B.fetch(contig):
            seen += 1
            if r.is_unmapped: continue
            if M.occupied(G['sp'], contig, E.tx_strand(r), r.reference_start + 1, r.reference_end): continue
            w.write(r); kept += 1
    return part, seen, kept

def run(gtf, bam, out, threads):
    import pysam
    G.update(sp=M.spans(gtf), bam=bam, out=out)
    refs = list(pysam.AlignmentFile(bam).references)
    with multiprocessing.get_context('fork').Pool(threads) as pool: res = pool.map(work, list(enumerate(refs)), chunksize=1)
    parts = [p for p, _, _ in res]
    pysam.cat('-o', out + '.tmp', *parts)
    for p in parts: os.remove(p)
    os.replace(out + '.tmp', out); pysam.index(out)
    seen = sum(s for _, s, _ in res); kept = sum(k for _, _, k in res)
    print(f'empty_place_bam: {kept} of {seen} records kept ({kept / max(1, seen):.3f}) -> {out}')
    return seen, kept

def selftest():
    import tempfile, pysam
    d = tempfile.mkdtemp(); gtf = os.path.join(d, 'x.gtf'); bam = os.path.join(d, 'x.bam')
    open(gtf, 'w').write('c\tx\ttranscript\t101\t600\t.\t+\t.\tgene_id "G"; transcript_id "T";\n')
    hdr = {'HD': {'VN': '1.6', 'SO': 'coordinate'}, 'SQ': [{'SN': 'c', 'LN': 5000}]}
    def rd(name, pos, flag):
        r = pysam.AlignedSegment(); r.query_name = name; r.reference_id = 0; r.reference_start = pos; r.cigarstring = '100M'; r.flag = flag; r.mapping_quality = 60; r.query_sequence = 'A' * 100
        return r
    with pysam.AlignmentFile(bam, 'wb', header=hdr) as w:
        for r in (rd('in', 150, 0), rd('antisense', 200, 16), rd('edge', 550, 0), rd('out', 1000, 0), rd('sec', 1200, 256)): w.write(r)
    pysam.index(bam); out = os.path.join(d, 'o.bam')
    assert run(gtf, bam, out, 1) == (5, 3)
    assert [r.query_name for r in pysam.AlignmentFile(out)] == ['antisense', 'out', 'sec']
    print('empty_place_bam selftest ok')

if __name__ == '__main__':
    if sys.argv[1:] == ['--selftest']: selftest()
    else:
        ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
        ap.add_argument('gtf'); ap.add_argument('bam'); ap.add_argument('out'); ap.add_argument('--threads', type=int, default=8)
        a = ap.parse_args(); run(a.gtf, a.bam, a.out, a.threads)
