#!/home/juanfra/miniforge3/envs/flair/bin/python
"""Convert an EXISTING coordinate-sorted BAM into the BED12 that `flair collapse -q` expects.

Why not just run `flair align`? Because that would re-align the reads with FLAIR's own
minimap2 settings (`-ax splice --secondary=no`), which would confound the comparison:
we want StringTie, FLAIR and isoseq collapse to all see the SAME alignments
(minimap2 2.31 -ax splice:hq -uf --eqx -Y -N 50 -p 0.1 --secondary=yes).

`flair align` = doalignment() then dofiltering(); this calls only dofiltering(), which is
the exact BAM->BED12 step, so the BED is what flair would have produced from this BAM.
"""
import argparse, sys, logging
from types import SimpleNamespace
from flair.flair_align import dofiltering

def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('-b', '--bam', required=True, help='coordinate-sorted input BAM')
    p.add_argument('-o', '--output', required=True, help='output base (writes <base>.bed)')
    p.add_argument('--quality', type=int, default=0, help='minimum MAPQ (flair align default: 0)')
    p.add_argument('--filtertype', default='removesup',
                   choices=['keepsup', 'removesup', 'separate'],
                   help='chimeric-alignment handling (flair align default: removesup)')
    p.add_argument('--remove_singleexon', action='store_true')
    p.add_argument('--remove_internal_priming', action='store_true')
    p.add_argument('-g', '--genome', default=None, help='only needed with --remove_internal_priming')
    p.add_argument('-f', '--gtf', default=None, help='only needed with --remove_internal_priming')
    p.add_argument('--intprimingthreshold', type=int, default=12)
    p.add_argument('--intprimingfracAs', type=float, default=0.75)
    a = p.parse_args()

    logging.basicConfig(level=logging.INFO, format='%(asctime)s %(message)s', stream=sys.stderr)
    args = SimpleNamespace(output=a.output, quality=a.quality, filtertype=a.filtertype,
                           remove_singleexon=a.remove_singleexon,
                           remove_internal_priming=a.remove_internal_priming,
                           genome=a.genome, gtf=a.gtf,
                           intprimingthreshold=a.intprimingthreshold,
                           intprimingfracAs=a.intprimingfracAs)
    dofiltering(args, a.bam)

if __name__ == '__main__':
    main()
