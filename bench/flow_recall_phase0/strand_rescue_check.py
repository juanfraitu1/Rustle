"""Does RUSTLE_STRAND_PURE_MINORITY recover the strand_swamp st_only misses?

For each strand_swamp case (antisense-dominated locus where rustle drops the minority-strand
gene), slice the locus and run rustle -L WITHOUT vs WITH the flag; check whether the target
annotated intron chain appears in rustle's output (FSM-equivalent: exact intron chain present).
OOM-safe: per-locus slices. The convergent antisense partner is co-located, so the slice
captures the bundle the flag acts on.
"""
import json, os, sys, subprocess
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import lib

def run_get_chains(slice_bam, chrom, out_gtf, flag):
    env = dict(os.environ)
    if flag:
        env["RUSTLE_STRAND_PURE_MINORITY"] = "1"
    subprocess.run([lib.RUSTLE, "-L", slice_bam, "-o", out_gtf], env=env,
                   stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    return {d["introns"] for d in lib.parse_gtf(out_gtf).values() if d["introns"]}

def main():
    probed = [json.loads(l) for l in open(f"{lib.ROOT}/bench/copy_recovery_eval/"
                                          f"results_genomewide/no_overlap_probed.jsonl")]
    strand = [p for p in probed if p["bucket"] == "strand_swamp"]
    print(f"strand_swamp cases: {len(strand)}\n")
    wd = "/tmp/strand_check"
    os.makedirs(wd, exist_ok=True)
    rec_off = rec_on = 0
    for p in strand:
        chrom, tid = p["chrom"], p["tid"]
        ch = lib.ref_chain(chrom, tid)
        if ch is None:
            continue
        lo, hi, introns = ch
        introns = tuple(introns)
        sb = f"{wd}/s.bam"
        reg = f"{chrom}:{max(1, lo-8000)}-{hi+8000}"
        with open(sb, "wb") as fh:
            subprocess.run(["samtools", "view", "-b",
                            f"{lib.PERCHROM}/{chrom}/c.bam", reg],
                           stdout=fh, stderr=subprocess.DEVNULL)
        subprocess.run(["samtools", "index", sb], stderr=subprocess.DEVNULL)
        off = introns in run_get_chains(sb, chrom, f"{wd}/off.gtf", False)
        on = introns in run_get_chains(sb, chrom, f"{wd}/on.gtf", True)
        rec_off += off
        rec_on += on
        flag_recovers = on and not off
        print(f"  {chrom} {tid:22s} sense={p['sense_frac']:.2f} "
              f"off={'Y' if off else '.'} on={'Y' if on else '.'}"
              f"{'  <- RECOVERED by flag' if flag_recovers else ''}")
    print(f"\n=== recovered (exact chain present) ===")
    print(f"  default (-L):                 {rec_off}/{len(strand)}")
    print(f"  +RUSTLE_STRAND_PURE_MINORITY: {rec_on}/{len(strand)}")
    print(f"  net gain from flag: +{rec_on - rec_off}")

if __name__ == "__main__":
    main()
