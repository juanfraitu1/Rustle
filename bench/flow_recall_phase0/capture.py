"""Driver: ensure cached cross-tool logs exist for every st_only miss. Resumable."""
import json, sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import lib

def main():
    with open(f"{lib.RG}/st_only_classified.jsonl") as fh:
        misses = [json.loads(l) for l in fh]
    done = 0
    for i, m in enumerate(misses):
        wd = lib.ensure_locus_logs(m["chrom"], m["tid"])
        if wd:
            done += 1
        if (i + 1) % 25 == 0:
            print(f"{i+1}/{len(misses)} captured ({done} ok)", flush=True)
    print(f"DONE {done}/{len(misses)}")

if __name__ == "__main__":
    main()
