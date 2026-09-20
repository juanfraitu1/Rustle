#!/usr/bin/env python3
"""Wave 3: archive committed experiment OUTPUTS, keeping code, docs and fixtures.

What a reviewer sees on first contact is the problem this solves: `bench/` holds ~1,040 committed
outputs (tsv/json/gfa/csv/png) of experiments whose scripts and write-ups are already in the tree, so
the repo reads as a scratch directory rather than a method.

KEEP: all source, all documentation, all test fixtures, and any data file the Rust build or tests
actually read. The last of those is DERIVED FROM THE SOURCES, not hardcoded -- `src/` builds some paths
with `concat!(env!("CARGO_MANIFEST_DIR"), "/bench/...")`, which a naive grep for `"bench/..."` misses.
Dropping one of those broke a test when this was first tried (§6r0).

DRY RUN BY DEFAULT; `--apply` performs `git mv` into archive/.
Validate afterwards the way §6r0 did: build the reduced tree in a scratch directory, run the full test
suite there, and re-run a REPRODUCE.md command -- its GTF must come back byte-identical.
"""
import argparse, os, re, subprocess, sys

ROOT = subprocess.run(['git','rev-parse','--show-toplevel'],capture_output=True,text=True).stdout.strip()
os.chdir(ROOT)
def git(*a): return [x for x in subprocess.run(['git',*a],capture_output=True,text=True).stdout.split('\n') if x]

def rust_referenced_paths():
    """Every repo-relative path the Rust sources or tests name, however the literal is built."""
    out=set()
    for d in ('src','tests'):
        for dp,_,fs in os.walk(d):
            for f in fs:
                if not f.endswith('.rs'): continue
                body=open(os.path.join(dp,f),errors='ignore').read()
                # "bench/x", "/bench/x" (concat! with CARGO_MANIFEST_DIR), "tests/fixtures/x"
                for m in re.findall(r'"/?((?:bench|tests|test_data|docs|analysis)/[A-Za-z0-9_./-]+)"', body):
                    out.add(m)
    return out

def main():
    ap=argparse.ArgumentParser(); ap.add_argument('--apply',action='store_true'); a=ap.parse_args()
    files=git('ls-files')
    needed=rust_referenced_paths()

    def keep(f):
        if f.startswith(('archive/','figures/')): return False
        if f in ('find_missed_loci.py','analyze_novel.py','analyze_missed_patterns.py'): return False
        if f in needed: return True
        if f.startswith('analysis/') and not f.endswith(('.R','.md')): return False
        if f.startswith('bench/'): return f.endswith(('.py','.sh','.md'))
        return True

    k=[f for f in files if keep(f)]; move=[f for f in files if not keep(f)]
    size=lambda L: sum(os.path.getsize(x) for x in L if os.path.exists(x))
    print(f"KEEP   {len(k):>5d} files {size(k)/1e6:>7.1f} MB")
    print(f"MOVE   {len(move):>5d} files {size(move)/1e6:>7.1f} MB")
    print(f"  data files kept because the Rust sources read them: "
          f"{sorted(f for f in needed if f.startswith(('bench/','analysis/')) and os.path.exists(f))}")
    if not a.apply:
        print("\n(dry run -- pass --apply to act)"); return
    n=0
    for p in move:
        if not os.path.exists(p): continue
        dest=os.path.join('archive',p); os.makedirs(os.path.dirname(dest),exist_ok=True)
        r=subprocess.run(['git','mv','-k','--',p,dest],capture_output=True,text=True)
        if r.returncode: print("  SKIP",p,r.stderr.strip(),file=sys.stderr); continue
        n+=1
    print(f"\nmoved {n}. Now validate: fresh build + full test suite + a REPRODUCE.md command (byte-identical GTF).")
main()
