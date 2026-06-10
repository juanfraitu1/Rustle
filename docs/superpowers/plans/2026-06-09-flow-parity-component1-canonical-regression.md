# Flow-Parity Component 1: Characterize the Canonical Regression — Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Name exactly why `RUSTLE_PARSE_TRFLONG_ST_CANONICAL=1` regresses (223 rustle-only vs ST, worse than the 187 default), producing a §6p findings note + a first-divergence worklist, and a verdict on whether the existing canonical port is *close-but-broken* (debug it) or *fundamentally-off* (rebuild it).

**Architecture:** Pure characterization — NO source changes to the flow. Regenerate the canonical / StringTie / default `path_extracted` traces, do a 3-way chain-set diff, classify the chains canonical adds/removes (TP/FP vs RefSeq, in-ST-final), cluster the divergent loci into a worklist, then first-divergence node-flux trace 2-3 representative loci. The only committed artifacts are a reusable bench script and the §6p findings note.

**Tech Stack:** Python 3 (analysis), the release `target/release/rustle`, instrumented `tools/stringtie/stringtie`, gffcompare not needed (exact intron-chain matching in Python).

**Spec:** `docs/superpowers/specs/2026-06-09-flow-enumeration-parity-port-design.md` (Component 1).

**Working directory for all commands:** `/mnt/c/Users/jfris/Desktop/Rustle`. Branch: `vg/flow-capacity-apportionment`.

**Key framing (do not lose):** rustle's flow is *different-not-worse* (it finds +19 real isoforms ST misses). This sub-step is characterization only — it does NOT change flow behavior. The 223-vs-187 regression means canonical extracts chains that even ST does not; Component 1 names them.

---

### Task 1: Regenerate the three ground-truth traces

`/tmp` is cleared between sessions, so regenerate. The instrumented StringTie is `tools/stringtie/stringtie`; `GGO_19.bam` is a repo-root symlink.

**Files:** none (produces `/tmp` artifacts).

- [ ] **Step 1: Regenerate StringTie trace (the parity reference)**

Run:
```bash
STRINGTIE_PARITY_LOG=/tmp/stP.jsonl tools/stringtie/stringtie -L GGO_19.bam -o /tmp/stP.gtf 2>/dev/null
echo "ST path_extracted: $(grep -c '"step":"path_extracted"\|"step": "path_extracted"' /tmp/stP.jsonl); tx: $(grep -cP '\ttranscript\t' /tmp/stP.gtf)"
```
Expected: ST path_extracted ≈ 5199; tx ≈ 1845.

- [ ] **Step 2: Regenerate default rustle trace (the 187 baseline)**

Run:
```bash
RUSTLE_PARITY_LOG=/tmp/ru_def.jsonl RUSTLE_PARITY_FILTER_STEPS=path_extracted,parse_trflong_seed \
  ./target/release/rustle -L GGO_19.bam -o /tmp/ru_def.gtf 2>/dev/null
python3 bench/gtf_chain_diff.py /tmp/ru_def.gtf /tmp/stP.gtf | tail -2
```
Expected: `Rustle-only: 187  ST-only: 104`.

- [ ] **Step 3: Regenerate canonical rustle trace (the 223 regression)**

Run:
```bash
RUSTLE_PARSE_TRFLONG_ST_CANONICAL=1 RUSTLE_PARITY_LOG=/tmp/ru_canon.jsonl \
  RUSTLE_PARITY_FILTER_STEPS=path_extracted,parse_trflong_seed \
  ./target/release/rustle -L GGO_19.bam -o /tmp/ru_canon.gtf 2>/dev/null
python3 bench/gtf_chain_diff.py /tmp/ru_canon.gtf /tmp/stP.gtf | tail -2
```
Expected: `Rustle-only: 223  ST-only:` (some value ≥104). If canonical rustle-only is NOT ~223, STOP and report — the binary or env may be stale (canonical must be `RUSTLE_PARSE_TRFLONG_ST_CANONICAL=1`).

---

### Task 2: Build the 3-way canonical-divergence classifier

A reusable script that, given the default/canonical/ST GTFs + the canonical path_extracted jsonl + the RefSeq GFF, reports: which final chains canonical ADDS vs default and vs ST, which it REMOVES, and the TP/FP nature (vs RefSeq) of each bucket. This answers "are the +36 canonical-adds new FPs (bad) or ST-extra recoveries (good)?"

**Files:**
- Create: `bench/canonical_divergence.py`

- [ ] **Step 1: Write a sanity-asserting test invocation (TDD-lite)**

There is no unit-test harness for bench scripts; the "test" is a self-check the script prints. Write the script (Step 2) to END with an assertion block that fails loudly if the chain-set arithmetic is inconsistent (added/removed counts must reconcile default↔canonical). Step 3 runs it and confirms the assertion passes.

- [ ] **Step 2: Create `bench/canonical_divergence.py`**

```python
#!/usr/bin/env python3
"""3-way canonical-vs-default-vs-StringTie chain-set divergence classifier.
Usage: python3 bench/canonical_divergence.py <default.gtf> <canonical.gtf> <st.gtf> <annot.gff> <chrom>
Reports the final chains canonical adds/removes vs default and vs ST, and their TP/FP nature.
"""
import re, sys, collections

def chains_gtf(path):
    tx = collections.defaultdict(list); strand = {}
    for line in open(path):
        if line.startswith('#'): continue
        f = line.rstrip('\n').split('\t')
        if len(f) < 9 or f[2] != 'exon': continue
        m = re.search(r'transcript_id "([^"]+)"', f[8])
        if not m: continue
        tx[m.group(1)].append((int(f[3]), int(f[4]))); strand[m.group(1)] = f[6]
    out = set()
    for t, ex in tx.items():
        ex.sort()
        if len(ex) < 2: continue
        out.add((strand[t], tuple((ex[i][1], ex[i+1][0]) for i in range(len(ex)-1))))
    return out

def chains_gff(path, chrom):
    tx = collections.defaultdict(list); strand = {}
    for line in open(path):
        if line.startswith('#'): continue
        f = line.rstrip('\n').split('\t')
        if len(f) < 9 or f[0] != chrom or f[2] != 'exon': continue
        m = re.search(r'Parent=([^;]+)', f[8])
        if not m: continue
        tx[m.group(1)].append((int(f[3]), int(f[4]))); strand[m.group(1)] = f[6]
    out = set()
    for t, ex in tx.items():
        ex.sort()
        if len(ex) < 2: continue
        out.add((strand[t], tuple((ex[i][1], ex[i+1][0]) for i in range(len(ex)-1))))
    return out

dft, can, stp, annot, chrom = sys.argv[1:6]
D = chains_gtf(dft); C = chains_gtf(can); S = chains_gtf(stp); A = chains_gff(annot, chrom)

def tpfp(chs):
    tp = sum(1 for k in chs if k in A); return tp, len(chs) - tp

added = C - D          # chains canonical introduced (not in default)
removed = D - C        # chains canonical dropped (were in default)
can_only = C - S       # canonical's rustle-only (the 223)
def_only = D - S       # default's rustle-only (the 187)
print(f"default tx-chains={len(D)} canonical tx-chains={len(C)} ST tx-chains={len(S)} annot={len(A)}")
print(f"canonical rustle-only(vs ST)={len(can_only)}  default rustle-only={len(def_only)}  delta={len(can_only)-len(def_only)}")
print()
for name, chs in [("ADDED by canonical (in C, not D)", added),
                  ("REMOVED by canonical (in D, not C)", removed)]:
    tp, fp = tpfp(chs)
    in_st = sum(1 for k in chs if k in S)
    print(f"{name}: n={len(chs)}  TP(RefSeq)={tp} FP={fp}  in_ST_final={in_st}")
# The decisive split: of the chains canonical ADDS, how many are ST-shared (good: matching ST)
# vs canonical-only (bad: new FPs ST does not have)?
add_st_shared = sum(1 for k in added if k in S)
add_canon_only_fp = sum(1 for k in added if k not in S and k not in A)
add_canon_only_tp = sum(1 for k in added if k not in S and k in A)
print()
print(f"  ADDED breakdown: ST-shared(converging)={add_st_shared}  "
      f"canonical-only-FP(regressing)={add_canon_only_fp}  canonical-only-TP(real)={add_canon_only_tp}")
# Reconciliation assertion: canonical_only == default_only + added_not_in_S - removed_not_in_S
import sys as _s
lhs = len(can_only)
rhs = len(def_only) + len([k for k in added if k not in S]) - len([k for k in removed if k not in S])
assert lhs == rhs, f"chain-set arithmetic inconsistent: can_only={lhs} != reconstructed={rhs}"
print("\n[OK] chain-set arithmetic reconciles")
```

- [ ] **Step 3: Run it and confirm the assertion passes + read the verdict**

Run:
```bash
python3 bench/canonical_divergence.py /tmp/ru_def.gtf /tmp/ru_canon.gtf /tmp/stP.gtf ../GGO_genomic.gff NC_073243.2
```
Expected: prints the ADDED/REMOVED breakdown and ends with `[OK] chain-set arithmetic reconciles`. The decisive line is the `ADDED breakdown`:
- If `canonical-only-FP (regressing)` dominates the additions → canonical introduces NEW FPs ST does not have → the port is *fundamentally-off* in those mechanisms.
- If `ST-shared (converging)` dominates → canonical is moving toward ST but losing elsewhere → *close-but-broken*.
Record the actual numbers.

- [ ] **Step 4: Commit the classifier**

```bash
git add bench/canonical_divergence.py
git commit -m "test(parity): add 3-way canonical-divergence chain classifier

Co-Authored-By: Claude Fable 5 <noreply@anthropic.com>"
```

---

### Task 3: Cluster the canonical-only FP chains into a first-divergence worklist

The chains canonical adds that are NOT in ST (the regressing additions) are the loci where canonical's `long_max_flow_st` diverges from ST. Cluster them by genomic locus so Component 2+ can trace one at a time.

**Files:** none (prints the worklist; it goes into the §6p note in Task 5).

- [ ] **Step 1: Print the clustered worklist**

Run:
```bash
python3 - <<'PY'
import re, collections
def cg(p):
    tx=collections.defaultdict(list); s={}
    for l in open(p):
        if l.startswith('#'): continue
        f=l.rstrip().split('\t')
        if len(f)<9 or f[2]!='exon': continue
        m=re.search(r'transcript_id "([^"]+)"',f[8])
        if not m: continue
        tx[m.group(1)].append((int(f[3]),int(f[4]))); s[m.group(1)]=f[6]
    o={}
    for t,ex in tx.items():
        ex.sort()
        if len(ex)<2: continue
        o[(s[t],tuple((ex[i][1],ex[i+1][0]) for i in range(len(ex)-1)))]=(ex[0][0],ex[-1][1])
    return o
D=cg('/tmp/ru_def.gtf'); C=cg('/tmp/ru_canon.gtf'); S=cg('/tmp/stP.gtf')
# canonical-only-FP additions = added by canonical AND not in ST
regress=[(k,C[k]) for k in C if k not in D and k not in S]
regress.sort(key=lambda kv:(kv[0][0][0], kv[1][0]))
# cluster by overlapping span, same strand
clusters=[]
for k,(s_,e_) in regress:
    strand=k[0]
    if clusters and clusters[-1]['strand']==strand and s_<=clusters[-1]['end']+1000:
        clusters[-1]['end']=max(clusters[-1]['end'],e_); clusters[-1]['n']+=1
    else:
        clusters.append({'strand':strand,'start':s_,'end':e_,'n':1})
print(f"canonical-only-FP chains (the regression): {len(regress)} in {len(clusters)} loci")
for c in sorted(clusters,key=lambda c:-c['n'])[:15]:
    print(f"  {c['strand']} {c['start']}-{c['end']}  chains={c['n']}")
PY
```
Expected: a list of loci (the worklist). Record the top loci — these are the Component-2 trace targets.

---

### Task 4: First-divergence node-flux trace on 2 representative loci

For the top 2 worklist loci, capture per-node flux from BOTH tools and find the first node where canonical's `long_max_flow_st` flux/path diverges from ST's `long_max_flow`. This confirms the divergence is in the flow algorithm (vs upstream graph/seed) and pins the mechanism.

**Files:** none (trace output recorded in §6p).

- [ ] **Step 1: Pick the top locus from Task 3 and trace StringTie**

Substitute `<LOCUS>` with the top locus span (e.g. `31707587-31722117`) from Task 3:
```bash
ST_TRACE_COV_NODES=1 tools/stringtie/stringtie -L GGO_19.bam -o /dev/null 2>/tmp/st_covnode.log
grep "ST_COVNODE" /tmp/st_covnode.log | awk -F'nid=' '{print $2}' | head -40 > /tmp/st_flux_locus.txt
echo "ST node-flux lines captured: $(wc -l < /tmp/st_flux_locus.txt)"
```
(ST_TRACE_COV_NODES dumps all nodes; filter to the locus coordinates by the `nid=<id> <start>-<end>` field when reading.)

- [ ] **Step 2: Trace rustle canonical at the same locus**

```bash
RUSTLE_PARSE_TRFLONG_ST_CANONICAL=1 RUSTLE_COV_DEBUG=1 RUSTLE_TRACE_LOCUS=<LOCUS> \
  ./target/release/rustle -L GGO_19.bam -o /dev/null 2>/tmp/ru_covnode.log
grep "COV_DEBUG" /tmp/ru_covnode.log | head -60 > /tmp/ru_flux_locus.txt
echo "rustle node-flux lines captured: $(wc -l < /tmp/ru_flux_locus.txt)"
```
Expected: both files non-empty for the locus.

- [ ] **Step 3: Identify the first divergent node**

Read `/tmp/st_flux_locus.txt` and `/tmp/ru_flux_locus.txt` side by side. For the shared backbone nodes (matched by coordinate), find the FIRST node where `nodeflux` (or the chosen path) differs materially between ST and canonical. Record: the node coordinate, ST's nodeflux/noderate, canonical's nodeflux/noderate, and which direction canonical diverges (extracts more flux → extra path; less → missed path). This is the concrete first-divergence for Component 2.

- [ ] **Step 4: Repeat for the 2nd locus**

Repeat Steps 1-3 for the 2nd-ranked worklist locus to confirm whether the divergence mechanism is the SAME across loci (one fix helps many) or locus-specific (many fixes).

---

### Task 5: Write the §6p findings note + verdict

**Files:**
- Modify: `docs/STRINGTIE_PARITY_FINDINGS.md` (append §6p)

- [ ] **Step 1: Append §6p**

Append a `## §6p — Canonical flow port regression characterized` section recording: the 3-way numbers (Task 2/3 actuals), the ADDED breakdown (ST-shared vs canonical-only-FP vs canonical-only-TP), the worklist of top loci (Task 3), the 2 first-divergence node-flux traces (Task 4), and a VERDICT:
  - **close-but-broken** (additions are mostly ST-shared; the canonical-only-FP regression is concentrated in 1-2 reproducible `long_max_flow_st` mechanisms) → debug the existing port; the worklist is the Component-2 task list.
  - **fundamentally-off** (additions are mostly canonical-only FPs spread across many unrelated mechanisms) → recommend rebuild over debug; escalate the debug-vs-rebuild decision to the user.

- [ ] **Step 2: Commit**

```bash
git add docs/STRINGTIE_PARITY_FINDINGS.md
git commit -m "docs(parity): characterize canonical flow-port regression (§6p)

Co-Authored-By: Claude Fable 5 <noreply@anthropic.com>"
```

- [ ] **Step 3: Update memory**

Append the §6p verdict + worklist pointer to `project_flow_parity_scope.md` and refresh its `MEMORY.md` index line.

---

## Notes for the implementer

- **This is characterization, not a flow change.** Do NOT edit `parse_trflong_st.rs` or any flow code in Component 1. The only code artifact is `bench/canonical_divergence.py`.
- **`/tmp` is cleared between sessions** — Task 1 regenerates everything; later tasks depend on Task 1's artifacts existing in the same session.
- **Trace JSONL gotcha:** the StringTie log has ~0.7% malformed lines; parse with `try/except json.JSONDecodeError: continue` if you load `/tmp/stP.jsonl` directly (the GTF-based classifier in Task 2 avoids this).
- **Default-OFF guard:** confirm `./target/release/rustle -L GGO_19.bam` (no env) still gives 187/104 before trusting any canonical measurement — the canonical changes must never touch the default path.
- **The verdict is the deliverable.** Component 1 exists to answer debug-vs-rebuild; do not start porting (Component 2) until §6p reaches a verdict.
