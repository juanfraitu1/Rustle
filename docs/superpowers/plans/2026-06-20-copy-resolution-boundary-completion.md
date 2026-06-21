# Copy-Resolution Boundary Completion Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Complete the copy-resolution boundary into reproducible tooling — an exhaustive per-family K=0 census, a splice-divergence resolver (Tier-2 lever), and the Tier-3 co-quantification proposition.

**Architecture:** Three independent Python bench units sharing the K=0 panel: Task 1 extends `bench/copy_resolution_census.py` to per-family labels; Task 2 is a standalone splice-divergence resolver; Task 3 adds a proposition + executable check to the theory note. No Rust, no pipeline integration.

**Tech Stack:** Python 3, pysam, samtools/minimap2 (already present), markdown + LaTeX for the proposition.

## Global Constraints

- Python bench tooling only; no Rust / pipeline integration.
- Substrate (absolute paths — the scripts must NOT rely on the working directory):
  - BAM `/home/juanfra/winloci_scratch/GGO.bam`
  - reference `/home/juanfra/winloci_scratch/GGO.fasta`
  - annotation `/home/juanfra/winloci_scratch/GGO_genomic.gff`
  - families `/mnt/c/Users/jfris/Desktop/Rustle/bench/denovo_families.tsv` (1130 families; header `family_id\tn_copies\tmembers`; `members` = comma-separated `DN_<chrom>_<start>_<n>` ids; chrom contains dots, e.g. `NC_073247.2`).
- K=0 definition (established, do not change): a co-located copy pair's cross-mapping reads = reads with the same `query_name` aligned to BOTH loci with ≥1 MAPQ-0 alignment. `frac_same` = fraction of those reads with `NM_A == NM_B`. Per-pair class: `low_support` if `n_xmap < 3`; else `k0_strict` if `frac_same == 1.0`, `k0` if `frac_same >= 0.95`, `resolvable` otherwise. (`k0_strict` is a subset of `k0` for counting.)
- Established headline (the check targets): of the n_xmap≥3 assignment-relevant pairs, ≈77.5% resolvable, ≈83% of cross-mapping reads resolvable; K=0 confined to the MAGEA / X-linked cancer-testis-antigen inverted-dup tail.
- K=0 MAGEA panel (NC_073247.2): pair0 161251228-161257000 ~ 161458538-161464324 (Tier-2, ~33% of junction-reads splice-resolvable); pair2 164381222-164384848 ~ 164442447-164446101 (Tier-3, splice-identical); pair3 164397061-164401095 ~ 164426194-164430228 (Tier-3, splice-identical). All three loci are members of family `DNFAM1`.
- Every component ships a deterministic check runnable standalone (`python3 bench/<file>.py`); Task 3's check goes in `bench/copy_assignment_theory_checks.py` (suite must still exit 0).
- Commit messages end with the two trailers shown in each task's commit step (Co-Authored-By + Claude-Session). Work on the current branch; add only the new/modified files named in each task — do not touch the unrelated pre-existing modified files in the tree.

---

## Task 1: Exhaustive per-family census

**Files:**
- Modify: `bench/copy_resolution_census.py`
- Create (output): `bench/copy_resolution_census.tsv`

**Interfaces:**
- Consumes (existing in the file): `classify_pair(chrom, sA, sB) -> dict|None`. For `n_xmap >= 1` the dict has keys `chrom, sA, sB, spanA, spanB, n_xmap, nm_same, nm_diff, frac_same` (where `spanA`/`spanB` are `(start, end)` tuples). Requires same chrom + non-overlapping loci. Calls `locus_span` (cached) and `fetch_locus`.
- Produces: `per_family_census(win=2_000_000) -> (rows, all_pairs)`, writing `bench/copy_resolution_census.tsv`; `check_census()`.

- [ ] **Step 1: Read the existing script and confirm the interface**

Run: `sed -n '1,10p;75,104p' bench/copy_resolution_census.py`
Confirm: line 3 opens the BAM with a **relative** path `pysam.AlignmentFile("GGO.bam","rb")`; lines ~80-103 are a **top-level** block that does `json.load(open('/tmp/coloc_pairs.json'))` (this crashes on re-run — `/tmp/coloc_pairs.json` is gone). Confirm `classify_pair`'s return keys match the Interfaces block above.

- [ ] **Step 2: Edit the BAM path to absolute, and delete the broken top-level run block**

Edit 1 — make the BAM path absolute (line 3):
```python
bam=pysam.AlignmentFile("/home/juanfra/winloci_scratch/GGO.bam","rb")
```
Edit 2 — DELETE the legacy full-run block: everything from the comment line
`# ---- Full run over co-located pairs that actually cross-map ----`
through the final `json.dump(xmap_pairs, open('/tmp/xmap_results.json','w'))` line (the block that reads `/tmp/coloc_pairs.json`). KEEP the MAGEA sanity loop just above it (the `for (s1,s2) in [...]: ... print("MAGEA", ...)` lines — harmless and fast).

- [ ] **Step 3: Append the per-family census + verdict + TSV writer + check**

Append to `bench/copy_resolution_census.py`:

```python
import re as _re, os as _os, subprocess as _sub

FAM_TSV = "/mnt/c/Users/jfris/Desktop/Rustle/bench/denovo_families.tsv"
GFF     = "/home/juanfra/winloci_scratch/GGO_genomic.gff"
OUT_TSV = "/mnt/c/Users/jfris/Desktop/Rustle/bench/copy_resolution_census.tsv"

# Memoize the BAM-fetch hot path. classify_pair calls fetch_locus by global name (late binding),
# so rebinding the global here makes exhaustive pairing bounded by distinct loci, not pair count.
_fetch_cache = {}
_orig_fetch_locus = fetch_locus
def fetch_locus(chrom, span):
    key = (chrom, span)
    if key not in _fetch_cache:
        _fetch_cache[key] = _orig_fetch_locus(chrom, span)
    return _fetch_cache[key]

def _parse_member(m):
    mm = _re.match(r"DN_(N[CW]_\d+\.\d+)_(\d+)_(\d+)$", m)
    return (mm.group(1), int(mm.group(2))) if mm else None

def _pair_class(p):
    if p["n_xmap"] < 3: return "low_support"
    if p["frac_same"] >= 1.0: return "k0_strict"
    if p["frac_same"] >= 0.95: return "k0"
    return "resolvable"

def per_family_census(win=2_000_000):
    """Per-family K=0 labels over ALL co-located copy pairs (same chrom, within `win`, threshold n_xmap>=1).
    low_support pairs (n_xmap<3) are counted in n_colocated_pairs but excluded from the family verdict."""
    rows, all_pairs = [], []
    for line in open(FAM_TSV):
        if line.startswith("family_id"): continue
        parts = line.rstrip("\n").split("\t")
        if len(parts) < 3: continue
        fid, ncopies, members = parts[0], parts[1], parts[2]
        locs = [x for x in (_parse_member(m) for m in members.split(",")) if x]
        by_chrom = {}
        for c, s in locs: by_chrom.setdefault(c, []).append(s)
        classed = []
        for c, starts in by_chrom.items():
            starts = sorted(set(starts))
            for i in range(len(starts)):
                for j in range(i + 1, len(starts)):
                    if starts[j] - starts[i] > win: break  # sorted -> no closer pair beyond the window
                    try:
                        p = classify_pair(c, starts[i], starts[j])
                    except Exception:
                        continue
                    if p is None or p["n_xmap"] < 1: continue
                    p["fid"] = fid
                    classed.append(p); all_pairs.append(p)
        if not classed: continue
        cls   = [_pair_class(p) for p in classed]
        n_k0  = cls.count("k0") + cls.count("k0_strict")
        n_k0s = cls.count("k0_strict")
        n_res = cls.count("resolvable")
        n_ar  = n_k0 + n_res
        if   n_ar == 0:      verdict = "not_assignment_relevant"
        elif n_k0 == n_ar:   verdict = "k0"
        elif n_res == n_ar:  verdict = "resolvable"
        else:                verdict = "mixed"
        k0_pairs = ";".join(
            f"{p['chrom']}:{p['spanA'][0]}-{p['spanA'][1]}~{p['spanB'][0]}-{p['spanB'][1]}"
            for p, k in zip(classed, cls) if k in ("k0", "k0_strict"))
        rows.append(dict(family_id=fid, n_copies=int(ncopies), n_colocated_pairs=len(classed),
                         n_assignment_relevant=n_ar, n_resolvable=n_res, n_k0=n_k0, n_k0_strict=n_k0s,
                         family_verdict=verdict, k0_pairs=k0_pairs))
    cols = ["family_id","n_copies","n_colocated_pairs","n_assignment_relevant","n_resolvable",
            "n_k0","n_k0_strict","family_verdict","k0_pairs"]
    with open(OUT_TSV, "w") as fh:
        fh.write("\t".join(cols) + "\n")
        for r in rows:
            fh.write("\t".join(str(r[c]) for c in cols) + "\n")
    return rows, all_pairs

def check_census():
    rows, pairs = per_family_census()
    ar   = [p for p in pairs if _pair_class(p) in ("k0", "k0_strict", "resolvable")]  # n_xmap>=3
    n    = len(ar)
    res  = sum(1 for p in ar if _pair_class(p) == "resolvable")
    k0   = sum(1 for p in ar if _pair_class(p) in ("k0", "k0_strict"))
    reads = sum(p["n_xmap"] for p in ar); diff = sum(p["nm_diff"] for p in ar)
    assert n >= 200, f"expected a broad assignment-relevant set, got {n}"
    assert 0.70 <= res / n <= 0.85, f"resolvable fraction off: {res/n:.2f}"
    assert 0.78 <= diff / reads <= 0.88, f"read-level resolvable off: {diff/reads:.2f}"
    dnfam1 = [r for r in rows if r["family_id"] == "DNFAM1"]
    assert dnfam1, "DNFAM1 missing from census"
    assert dnfam1[0]["family_verdict"] in ("k0", "mixed"), f"DNFAM1 verdict={dnfam1[0]['family_verdict']}"
    assert dnfam1[0]["n_k0"] >= 2, f"DNFAM1 must carry >=2 K0 pairs (pair2/pair3), got n_k0={dnfam1[0]['n_k0']}"
    for r in rows:
        assert r["n_resolvable"] + r["n_k0"] == r["n_assignment_relevant"], f"{r['family_id']} verdict-count mismatch"
        assert r["n_assignment_relevant"] <= r["n_colocated_pairs"]
        assert r["n_k0_strict"] <= r["n_k0"]
    print(f"OK  - census: {len(rows)} families; assignment-relevant {n} pairs -> {res} resolvable "
          f"({100*res/n:.0f}%), {k0} K0; reads {100*diff/reads:.0f}% resolvable; "
          f"DNFAM1 verdict={dnfam1[0]['family_verdict']} n_k0={dnfam1[0]['n_k0']}")
    return rows

if __name__ == "__main__":
    check_census()
```

- [ ] **Step 4: Run the check**

Run: `python3 bench/copy_resolution_census.py 2>&1 | tail -4`
Expected (≈, the exact pair count may exceed 258 now that pairing is exhaustive rather than sampled):
`OK  - census: <N> families; assignment-relevant ~258 pairs -> ~200 resolvable (~77%), ~39 K0; reads ~83% resolvable; DNFAM1 verdict=k0 n_k0=2` (verdict may be `mixed` if pair0 lands resolvable). Runtime ≈ 2–5 min (fetches bounded by distinct loci via the cache). The first lines are the kept MAGEA sanity prints — fine.

- [ ] **Step 5: Confirm the TSV**

Run: `head -1 bench/copy_resolution_census.tsv; grep -P "^DNFAM1\t" bench/copy_resolution_census.tsv`
Expected: the 9-column header; a `DNFAM1` row with `family_verdict` `k0`/`mixed`, `n_k0>=2`, and `k0_pairs` listing the 164.38M~164.44M / 164.39M~164.42M spans.

- [ ] **Step 6: Commit**

```bash
git add bench/copy_resolution_census.py bench/copy_resolution_census.tsv
git commit -m "$(printf '%s\n' "bench(census): exhaustive per-family K=0 labels + per-family TSV" "" "Co-Authored-By: Claude Opus 4.8 (1M context) <noreply@anthropic.com>" "Claude-Session: https://claude.ai/code/session_0193po8vtwu5dbDNComynKAa")"
```

---

## Task 2: Splice-divergence negative-result probe (Tier-2 collapses into Tier-3)

> **REFRAMED during execution (user decision).** The original plan tried a per-read splice *resolver*
> (assign reads via copy-specific splice sites, targeting the prior "MAGEA pair0 ~33%"). Hands-on testing
> proved that is not achievable per-read: **minimap2's spliced aligner independently snaps every junction to
> the nearest canonical `GT-AG` at EACH copy**, so a read is fully canonical at *both* copies (intron lengths
> differing only by the intronic indel, e.g. 3702 nt at A vs 3705 nt at B). The 3 bp indel sits inside the
> long intron — absent from the spliced mature read, accommodated without an exonic edit — so NM stays tied
> (these are K=0 by construction) and there is **no per-read direction signal**. The prior "33%" was a
> reference-level junction count, not a per-read assignment. **Decision: collapse Tier-2 into Tier-3** — the
> deliverable is an honest negative-result probe + the boundary doc update.

**Files:**
- Modify: `bench/splice_divergence_resolver.py` (a prior commit `5f2517e` created it as a resolver; rewrite to the negative-result probe below)
- Modify: `bench/copy_resolution_boundary.md` (collapse Tier-2 into Tier-3)

**Interfaces:**
- Produces: `probe_pair(bam_path, fasta_path, locusA, locusB) -> dict` where each locus is `(chrom, start, end)`; returns `{n_reads, n_junction_reads, n_resolved, resolved_fraction, n_chain_divergent, chain_divergent_fraction, distinguishing_junctions}`.

**What the probe shows (the negative result, made rigorous):** for each cross-mapping junction-read it measures TWO things — (1) `resolved_fraction`: per-read direction via canonicality (fully-canonical at one copy, degraded at the other) → **~0 for ALL pairs incl. pair0**, because the aligner snaps; (2) `chain_divergent_fraction`: whether the read's intron-length multiset differs between its A- and B-alignment (reflection preserves intron lengths, so a length difference is a genuine splice divergence) → **HIGH for pair0** (the divergence IS detectable) but **~0 for the splice-identical pair2/pair3**. Conclusion: pair0's divergence is *detectable but non-directional* — no per-read rescue — so the K=0 residual collapses into Tier-3 (co-quantify).

- [ ] **Step 1: Rewrite the script as the negative-result probe + check**

Replace the entire contents of `bench/splice_divergence_resolver.py` with:

```python
"""Splice-divergence NEGATIVE-RESULT probe (Tier-2 collapses into Tier-3).

For K=0 co-located inverted-dup copies, intronic divergence near a splice site (e.g. the MAGEA pair0 3 bp
intronic indel) is REAL in the reference but is NOT usable for per-read copy assignment: minimap2's spliced
aligner independently snaps every junction to the nearest canonical GT-AG site at EACH copy, so a read is
fully canonical at both copies (intron lengths differing only by the indel, e.g. 3702 nt at A vs 3705 nt at B).
The indel sits inside the long intron -- absent from the spliced mature read, accommodated without an exonic
edit -- so NM stays tied (these are K=0 by construction) and there is no per-read direction signal.

This probe quantifies that on the MAGEA panel:
- resolved_fraction        (per-read canonicality direction)  -> ~0 for ALL pairs incl. pair0 (the negative
                                                                 result: aligner snapping masks the divergence)
- chain_divergent_fraction (read's A-intron-length multiset != B's) -> HIGH at pair0 (divergence IS detectable)
                                                                 but ~0 at splice-identical pair2/pair3
=> detectable yet NON-DIRECTIONAL: the K=0 residual has no per-read splice rescue; it collapses into Tier-3."""
import pysam

CANONICAL = {("GT", "AG"), ("CT", "AC"), ("GC", "AG"), ("CT", "GC"), ("AT", "AC"), ("GT", "AT")}

def _introns(read):
    """Genomic introns as (donor, acceptor) 0-based half-open: donor=intron first base, acceptor=intron end."""
    out, pos = [], read.reference_start
    for op, ln in read.cigartuples:
        if op in (0, 2, 7, 8):       # M/D/=/X consume reference
            pos += ln
        elif op == 3:                 # N = intron
            out.append((pos, pos + ln)); pos += ln
    return out

def _splice(fasta, chrom, d, a):
    return (fasta.fetch(chrom, d, d + 2).upper(), fasta.fetch(chrom, a - 2, a).upper())

def _canon(introns, fasta, chrom):
    return sum(1 for (d, a) in introns if _splice(fasta, chrom, d, a) in CANONICAL)

def probe_pair(bam_path, fasta_path, locusA, locusB):
    bam = pysam.AlignmentFile(bam_path, "rb"); fasta = pysam.FastaFile(fasta_path)
    cA, sA, eA = locusA; cB, sB, eB = locusB

    def reads_at(c, s, e):
        d = {}
        for r in bam.fetch(c, s, e):
            if r.is_unmapped or r.is_supplementary:
                continue
            if min(r.reference_end, e) - max(r.reference_start, s) < 200:
                continue
            d.setdefault(r.query_name, r)   # one alignment per query at this locus
        return d

    A, B = reads_at(cA, sA, eA), reads_at(cB, sB, eB)
    common = set(A) & set(B)
    distinguishing, n_junction, n_resolved, n_chain_div = set(), 0, 0, 0
    for q in common:
        ia, ib = _introns(A[q]), _introns(B[q])
        if not ia and not ib:
            continue
        n_junction += 1
        # (1) DETECTABLE divergence: intron-length multiset differs between the two alignments.
        #     Reflection preserves intron lengths, so a length difference is a genuine splice divergence.
        if sorted(a - d for d, a in ia) != sorted(a - d for d, a in ib):
            n_chain_div += 1
        # (2) per-read DIRECTION via canonicality -> ~0: minimap2 snaps BOTH copies to canonical sites.
        ca, cb = _canon(ia, fasta, cA), _canon(ib, fasta, cB)
        if ia and ca == len(ia) and ib and cb < len(ib):        # A fully canonical, B degraded -> copy A
            n_resolved += 1
            for (d, a) in ib:
                if _splice(fasta, cB, d, a) not in CANONICAL:
                    distinguishing.add((cB, d, a))
        elif ib and cb == len(ib) and ia and ca < len(ia):      # B fully canonical, A degraded -> copy B
            n_resolved += 1
            for (d, a) in ia:
                if _splice(fasta, cA, d, a) not in CANONICAL:
                    distinguishing.add((cA, d, a))
    return dict(n_reads=len(common), n_junction_reads=n_junction, n_resolved=n_resolved,
                resolved_fraction=(n_resolved / n_junction if n_junction else 0.0),
                n_chain_divergent=n_chain_div,
                chain_divergent_fraction=(n_chain_div / n_junction if n_junction else 0.0),
                distinguishing_junctions=sorted(distinguishing))


PANEL = {
    "pair0": (("NC_073247.2", 161251228, 161257000), ("NC_073247.2", 161458538, 161464324)),
    "pair2": (("NC_073247.2", 164381222, 164384848), ("NC_073247.2", 164442447, 164446101)),
    "pair3": (("NC_073247.2", 164397061, 164401095), ("NC_073247.2", 164426194, 164430228)),
}

def check_probe():
    BAM   = "/home/juanfra/winloci_scratch/GGO.bam"
    FASTA = "/home/juanfra/winloci_scratch/GGO.fasta"
    res = {k: probe_pair(BAM, FASTA, a, b) for k, (a, b) in PANEL.items()}
    for k, r in res.items():
        print(f"    {k}: junction_reads={r['n_junction_reads']} "
              f"per-read-resolved={r['resolved_fraction']:.2f} "
              f"chain-divergent={r['chain_divergent_fraction']:.2f} "
              f"distinguishing_junc={len(r['distinguishing_junctions'])}")
    # NEGATIVE RESULT: per-read splice resolution is ~0 for ALL pairs incl. pair0 -- aligner junction-snapping
    # masks the intronic divergence (the prior "33%" was a reference-level junction count, not per-read).
    for k in PANEL:
        assert res[k]["resolved_fraction"] <= 0.05, f"{k} per-read resolved should be ~0 (aligner-masked), got {res[k]['resolved_fraction']:.2f}"
    # ...but the divergence IS real and DETECTABLE as intron-length divergence at pair0, and ABSENT at the
    # splice-identical pairs -- detectable yet non-directional (so it cannot rescue per-read assignment).
    assert res["pair0"]["chain_divergent_fraction"] >= 0.15, f"pair0 divergence should be detectable, got {res['pair0']['chain_divergent_fraction']:.2f}"
    assert res["pair2"]["chain_divergent_fraction"] <= 0.05, f"pair2 splice-identical, got {res['pair2']['chain_divergent_fraction']:.2f}"
    assert res["pair3"]["chain_divergent_fraction"] <= 0.05, f"pair3 splice-identical, got {res['pair3']['chain_divergent_fraction']:.2f}"
    print(f"OK  - splice probe (NEGATIVE per-read): pair0 per-read-resolved={res['pair0']['resolved_fraction']:.2f} "
          f"but chain-divergent={res['pair0']['chain_divergent_fraction']:.2f} (detectable, non-directional); "
          f"pair2/3 resolved~0 chain-divergent~0 (splice-identical). Tier-2 -> Tier-3.")
    return res

if __name__ == "__main__":
    check_probe()
```

- [ ] **Step 2: Run the check**

Run: `python3 bench/splice_divergence_resolver.py 2>&1 | tail -6`
Expected: three per-pair lines (all `per-read-resolved=0.00`; pair0 `chain-divergent` high ≥0.15, pair2/pair3 ~0), then `OK  - splice probe (NEGATIVE per-read): ... Tier-2 -> Tier-3.`. Load-bearing asserts: every pair `resolved_fraction <= 0.05`; pair0 `chain_divergent_fraction >= 0.15`; pair2/pair3 `chain_divergent_fraction <= 0.05`. The run is fast (small windows). If pair0's `chain_divergent_fraction` comes in below 0.15, that contradicts the observed 3702-vs-3705 intron divergence — recheck `_introns` before concluding; do NOT weaken the assert.

- [ ] **Step 3: Collapse Tier-2 into Tier-3 in `bench/copy_resolution_boundary.md`**

Make these edits (the section headings are stable match anchors; rewrite each named region, preserving everything else):

(a) **Bottom-line** — the sentence "they are recent inverted segmental duplications … so no spliced read carries copy information." is fine; but in the same "Bottom line" section, if any clause credits splice sites with rescuing the residual, drop it. Then in the **"### Tier 2 …"** section (heading currently `### Tier 2 — K=0 residual, partly rescued by copy-specific splice sites (~15% of pairs)`) REPLACE the entire section (heading + body) with:

```markdown
### Tier 2 — K=0 residual: splice divergence is real in the reference but per-read-masked (~15-18% of pairs)

The exons are identical, so per-base PSVs give 0. The *introns* do diverge (4–60 SNPs/indels across the
duplicated block), and where that lands at a splice site it shifts a junction (e.g. MAGEA pair0's 3 bp
intronic indel makes copy A's intron 3702 nt and copy B's 3705 nt). One might hope a read using that junction
betrays its copy — but it does **not, per read**: minimap2's spliced aligner independently snaps every junction
to the nearest canonical `GT-AG` at *each* copy, so the read is fully canonical at **both** copies and, because
the indel lives inside the long intron (never in the spliced mature read, and absorbed without an exonic edit),
NM stays tied. Measured on the panel (`bench/splice_divergence_resolver.py`): per-read splice **resolution = 0
for every pair, including pair0** — the divergence is *detectable* (pair0's reads carry a different intron-length
chain at A vs B) but **non-directional** (it cannot say which copy a read came from). The earlier "pair0 ~33%"
was a reference-level count of junctions whose homologous site differs, not a per-read assignment.

So there is **no Tier-2 per-read rescue**: the K=0 residual collapses into Tier-3. (The reference-level splice
distinguishability remains relevant to copy-*model* / family-graph work — interest #3 — but not to assigning an
individual mature read.)
```

(b) **The "## The three tiers" heading** → change to `## The two tiers`. **The Tier-3 heading** (currently `### Tier 3 — The irreducible core (~3% strict K=0)`) → broaden to `### Tier 3 — The K=0 core: co-quantify (no per-read assignment)`, and add one sentence at the end of that section: "With Tier-2 shown to be per-read-masked, the *entire* K=0 residual (not only the strict core) lands here: report a co-quantified ambiguity set, not forced per-read assignments."

(c) **The boundary diagram** (the fenced ```` ``` ```` block under `## The boundary as the result`) — REPLACE its body with:

```
multi-copy family with cross-mapping reads
        │
        ├─ exonic PSVs (NM_A≠NM_B)         → Tier 1  ~77%  → per-read PSV assignment (Thm 2/3)        [interest #2]
        │
        └─ K=0 (exons identical)            ~15-18%
              → splice divergence is reference-real but per-read-masked (minimap2 snaps junctions)
              → Tier 2 has NO per-read rescue → co-quantified ambiguity set (Tier 3)                  [interest #3]
```

(d) **Reproducibility** — under the `## Reproducibility` or `## Open / next` list, replace the "Productize the donor/acceptor … splice-divergence scan as the opt-in Tier-2 resolver." bullet with: "Splice-divergence is per-read-masked by aligner junction-snapping (`bench/splice_divergence_resolver.py` — negative-result probe); the K=0 residual co-quantifies (Tier-3)." Leave the other bullets.

After editing, sanity-check: `grep -n "Tier 2\|two tiers\|per-read-masked\|33%" bench/copy_resolution_boundary.md` — there should be no surviving claim that splice sites *rescue* a fraction of reads, and no stale "three tiers".

- [ ] **Step 4: Commit**

```bash
git add bench/splice_divergence_resolver.py bench/copy_resolution_boundary.md
git commit -m "$(printf '%s\n' "bench(splice): negative-result probe — per-read splice resolution is aligner-masked; collapse Tier-2 into Tier-3" "" "minimap2 snaps junctions to canonical at both copies; pair0 intronic indel is detectable (chain divergence) but non-directional (resolution 0). Boundary doc updated to two tiers." "" "Co-Authored-By: Claude Opus 4.8 (1M context) <noreply@anthropic.com>" "Claude-Session: https://claude.ai/code/session_0193po8vtwu5dbDNComynKAa")"
```

---

## Task 3: Tier-3 co-quantification proposition + check

**Files:**
- Modify: `bench/copy_assignment_theory.md` (add the Tier-3 Proposition + proof after the §6 Corollary)
- Modify: `bench/copy_assignment_theory_checks.py` (add `check_tier3_coquant_unidentifiable`)

**Interfaces:**
- Consumes: nothing new (self-contained likelihood witness).
- Produces (checks): `check_tier3_coquant_unidentifiable`.

- [ ] **Step 1: Find the check-registration idiom**

Run: `grep -n "^def check_\|CHECKS\|^if __name__\|for .* in .*check" bench/copy_assignment_theory_checks.py | head -40`
Note how checks are collected and run (e.g. a `CHECKS = [...]` list, or the runner enumerates `check_*` functions). You will register the new check the same way.

- [ ] **Step 2: Add the executable witness**

Add to `bench/copy_assignment_theory_checks.py` (place the function with the other `check_*` defs; register it the SAME way they are — append to the `CHECKS` list if one exists, otherwise add it wherever the runner enumerates checks):

```python
import math as _math

def check_tier3_coquant_unidentifiable():
    """Tier-3 witness: under K=0 (copies identical over the transcript) the mixture log-likelihood is CONSTANT
    over the per-copy abundance simplex (per-copy split unidentifiable), while the aggregate N is fixed."""
    N = 10
    p = 0.3  # P(read | copy) — identical across copies because the copies are identical
    def loglik(a0, a1):
        assert abs((a0 + a1) - N) < 1e-9
        per_read = _math.log((a0 / N) * p + (a1 / N) * p)  # = log(p), independent of (a0, a1)
        return N * per_read
    splits = [(N, 0), (N / 2, N / 2), (1, N - 1), (3, N - 3)]
    vals = [loglik(*s) for s in splits]
    assert max(vals) - min(vals) < 1e-9, f"likelihood must be flat over the simplex, spread={max(vals)-min(vals)}"
    assert abs(loglik(N, 0) - loglik(0, N)) < 1e-9
    return "Tier-3: K=0 per-copy split unidentifiable (flat likelihood over the simplex); aggregate identifiable"
```
Then register it (e.g. `CHECKS.append(check_tier3_coquant_unidentifiable)` if the file uses a `CHECKS` list, matching what Step 1 found).

- [ ] **Step 3: Run the suite**

Run: `python3 bench/copy_assignment_theory_checks.py 2>&1 | tail -5`
Expected: all checks `OK` including the new Tier-3 line, process exits 0. Confirm the new check is actually invoked (the count of OK lines increased by one).

- [ ] **Step 4: Write the Proposition + proof into the note**

In `bench/copy_assignment_theory.md`, after the §6 Corollary (paths/isoforms), add:

```markdown
## §6b Tier-3: co-quantification of the irreducible core

For copies identical over the entire transcribed region (the K=0 / Strong-Separation-fails regime), per-read
assignment is impossible (no distinguishing column). The natural fallback is per-copy *quantification* — but the
same identical-sequence condition makes that unidentifiable too, in a precise sense.

> **Proposition (Tier-3 unidentifiability).** Let a family have copies $c_1,\ldots,c_K$ identical over the
> transcribed region, and let reads $R$ ($|R| = N$) each be consistent with every copy. Under the mixture
> likelihood $L(a) = \prod_{r\in R}\sum_{k=1}^{K}\frac{a_k}{N}\,P(r\mid c_k)$ with per-copy abundances
> $a=(a_1,\ldots,a_K)$, $a_k\ge 0$, $\sum_k a_k = N$, the likelihood $L(a)$ is **constant** on the simplex
> $\{a:\sum_k a_k=N\}$: the per-copy split is statistically **unidentifiable** from RNA. The aggregate $N=|R|$
> is identifiable (a sufficient statistic). Under a copy-number / dosage prior $\pi(a)$, the MAP estimate is
> well-posed and equals the mode of $\pi$ scaled to $N$ — RNA contributes nothing to the per-copy direction.

*Proof.* The copies are identical, so $P(r\mid c_k)=p_r$ is the same for every $k$. Each read's mixture term is
$\sum_k\frac{a_k}{N}p_r = p_r\frac{\sum_k a_k}{N}=p_r$, independent of $a$. Hence $L(a)=\prod_r p_r$ is constant
in $a$; it factors through $N$ alone, so $N$ is identifiable and $a$ is not. With prior $\pi$, the posterior is
$\propto \pi(a)\,L(a)\propto\pi(a)$ on the simplex, so the MAP is the mode of $\pi$ (scaled to $N$). $\square$

This is the identifiability theorem's K=0 floor in the quantification frame: the same emptiness of the
distinguishing-column set that forbids per-read assignment (Theorem 2 corollary) flattens the per-copy
likelihood. The honest Tier-3 output is therefore the **family aggregate** (identifiable) plus a
prior-conditioned per-copy split whose uncertainty set is the entire simplex. (Read error perturbs each $p_r$
but not the $a$-independence of the mixture term, so the unidentifiability is robust to symmetric error.)
`check_tier3_coquant_unidentifiable` is the executable witness: every apportionment of a fixed $N$ yields
identical likelihood.
```

- [ ] **Step 5: Re-run the suite (safety after the markdown edit)**

Run: `python3 bench/copy_assignment_theory_checks.py 2>&1 | tail -3`
Expected: exit 0, all `OK`. Spot-read §6b: the proposition statement matches the witness (flat likelihood / identifiable aggregate / prior-conditioned MAP).

- [ ] **Step 6: Commit**

```bash
git add bench/copy_assignment_theory.md bench/copy_assignment_theory_checks.py
git commit -m "$(printf '%s\n' "theory(copy-assign): Tier-3 co-quantification proposition + witness check" "" "Per-copy split unidentifiable under K=0 (flat likelihood over the simplex); aggregate N identifiable; prior-conditioned MAP." "" "Co-Authored-By: Claude Opus 4.8 (1M context) <noreply@anthropic.com>" "Claude-Session: https://claude.ai/code/session_0193po8vtwu5dbDNComynKAa")"
```

---

## Self-Review

**Spec coverage:** Component 1 (exhaustive per-family census + per-family TSV + check) → Task 1; Component 2 (splice-divergence resolver `resolve_pair` + panel-pinned check) → Task 2; Component 3 (Tier-3 Proposition + proof in `copy_assignment_theory.md` + `check_tier3_coquant_unidentifiable`) → Task 3. All three spec components covered.

**Placeholder scan:** none — full code for every function and check, exact commands with expected output. The conditional instructions ("match the file's check-registration idiom", "if pair0 underperforms, diagnose") are real integration steps resolved by reading the file / output, not deliverable gaps.

**Type/name consistency:**
- `classify_pair` returns `frac_same, nm_diff, n_xmap, spanA, spanB` — used identically by `_pair_class`, `per_family_census`, `check_census`. K=0 thresholds (`frac_same >= 0.95` / `== 1.0`, `n_xmap < 3`) match the Global Constraints and the existing census.
- (Task 2 reframed mid-execution — see the box at its head.) `probe_pair` returns `resolved_fraction, n_chain_divergent, chain_divergent_fraction, …` — `check_probe` pins the NEGATIVE result (per-read `resolved_fraction ≤ 0.05` for all pairs incl. pair0, since the aligner snaps junctions) plus the detectable-but-non-directional divergence (pair0 `chain_divergent_fraction ≥ 0.15`, pair2/3 ≤ 0.05). Tier-2 collapses into Tier-3 in `copy_resolution_boundary.md`.
- The fetch-cache wrapper relies on Python late binding (`classify_pair` resolves the global `fetch_locus` at call time) — correct, since `classify_pair` references `fetch_locus` by name, not by a captured reference.
- `check_tier3_coquant_unidentifiable` registered via the existing idiom (Step 1 confirms); the proposition's symbols ($a$, $N$, $p_r$, $\pi$) match the witness's variables (`a0/a1`, `N`, `p`).

**Risk note (flagged, not a gap):** Task 1 lowers the threshold to `n_xmap >= 1`, so the exhaustive pair set is larger than the original 258; the `check_census` bands (`n >= 200`, `res/n ∈ [0.70,0.85]`, `reads ∈ [0.78,0.88]`) are deliberately tolerant of the count shift while still pinning the established ≈77%/≈83% headline. The `win=2_000_000` window + sorted-`break` keeps mega-families bounded; the fetch cache keeps runtime bounded by distinct loci.
```
