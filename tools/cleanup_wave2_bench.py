#!/usr/bin/env python3
"""Wave 2: invert the default in bench/ -- archive every script that no authoritative document names.

`bench/` is the lab notebook: 781 tracked scripts, of which only ~20% are reachable from anything that
records a result. Rather than triage the rest one by one (the audit's transitive anchoring is unreliable
in both directions -- its own KEEP-CITED control was 4/12 dead), this inverts the default: KEEP exactly
what an authoritative document names, plus the closure of what those keepers import, and archive the rest.

Authoritative = the ledger, the negative-results register, the pre-registrations, the top-level docs
(README, REPRODUCE, AGENTS, docs/*.md incl. ACTIVE_WORKING_SET and DATA), the bench/*.md write-ups (a
script named by a report is that report's provenance), and the Rust sources and tests.

Only `.py`/`.sh` under bench/ are considered. Reports, fixtures and data stay where they are.
DRY RUN BY DEFAULT; `--apply` performs `git mv` into archive/ (tracked) or a plain move (untracked).
"""
import argparse, os, re, subprocess, sys, collections

ROOT = subprocess.run(['git','rev-parse','--show-toplevel'],capture_output=True,text=True).stdout.strip()
os.chdir(ROOT)

def git(*a):
    return subprocess.run(['git',*a],capture_output=True,text=True).stdout.split('\n')

def authoritative_texts():
    pats = ['README.md','REPRODUCE.md','AGENTS.md']
    pats += [f'docs/{f}' for f in os.listdir('docs') if f.endswith('.md')]
    for d in ('docs/experiments',):
        if os.path.isdir(d): pats += [f'{d}/{f}' for f in os.listdir(d) if f.endswith('.md')]
    pats += [p for p in git('ls-files','bench/*.md') if p]
    for d in ('src','tests'):
        for dp,_,fs in os.walk(d):
            pats += [os.path.join(dp,f) for f in fs if f.endswith('.rs')]
    out=[]
    for p in pats:
        if os.path.isfile(p):
            try: out.append(open(p,errors='ignore').read())
            except OSError: pass
    return '\n'.join(out)

def main():
    ap=argparse.ArgumentParser(); ap.add_argument('--apply',action='store_true'); a=ap.parse_args()

    scripts=[p for p in git('ls-files','bench/*.py','bench/*.sh') if p]
    untracked=[p for p in git('ls-files','--others','--exclude-standard','bench/*.py','bench/*.sh') if p]
    allscripts=sorted(set(scripts)|set(untracked))

    text=authoritative_texts()
    named={m for m in re.findall(r'bench/[A-Za-z0-9_./-]+\.(?:py|sh)', text)}
    keep={p for p in allscripts if p in named}

    # transitive closure over imports (python) and sourced/invoked siblings (shell)
    by_stem=collections.defaultdict(list)
    for p in allscripts: by_stem[os.path.splitext(os.path.basename(p))[0]].append(p)
    frontier=list(keep); closure=set()
    while frontier:
        p=frontier.pop()
        if not os.path.isfile(p): continue
        try: body=open(p,errors='ignore').read()
        except OSError: continue
        refs=set(re.findall(r'^\s*(?:from|import)\s+([A-Za-z_][A-Za-z0-9_]*)', body, re.M))
        refs|=set(re.findall(r'bench/([A-Za-z0-9_./-]+)\.(?:py|sh)', body))
        refs|={os.path.splitext(os.path.basename(x))[0] for x in re.findall(r'(?:source|bash|sh|python3?)\s+([A-Za-z0-9_./-]+\.(?:py|sh))', body)}
        for r in refs:
            for t in by_stem.get(os.path.basename(r), []):
                if t not in keep and t not in closure:
                    closure.add(t); frontier.append(t)
    keep_all = keep | closure
    move = [p for p in allscripts if p not in keep_all]

    print(f"bench scripts: {len(allscripts)}  (tracked {len(scripts)}, untracked {len(untracked)})")
    print(f"  named by an authoritative doc : {len(keep)}")
    print(f"  pulled in by import closure   : {len(closure)}")
    print(f"  KEEP total                    : {len(keep_all)}")
    print(f"  ARCHIVE                       : {len(move)}")
    if closure:
        print("  closure additions:", ', '.join(sorted(os.path.basename(c) for c in closure)[:12]))
    if not a.apply:
        print("\n(dry run -- pass --apply to act)  first 10 to archive:")
        for p in move[:10]: print("   ", p)
        return

    tracked=set(scripts); n=0
    for p in move:
        if not os.path.exists(p): continue
        dest=os.path.join('archive', p) if p in tracked else os.path.join('archive','untracked', p)
        os.makedirs(os.path.dirname(dest), exist_ok=True)
        if p in tracked:
            r=subprocess.run(['git','mv','-k','--',p,dest],capture_output=True,text=True)
            if r.returncode: print("  SKIP",p,r.stderr.strip(),file=sys.stderr); continue
        else:
            os.replace(p,dest)
        n+=1
    print(f"\nmoved {n} scripts into archive/. Next: cargo test, re-run the audit, commit.")

main()
