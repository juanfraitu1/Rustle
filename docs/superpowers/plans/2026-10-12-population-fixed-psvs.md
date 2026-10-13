# Population-fixed PSVs Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Mark every chrY PSV site as population-FIXED / POLYMORPHIC / UNCOVERED / NOT_ASSESSED from 43+ human Y
assemblies, and let O2 report a portability tier per contested molecule (`--psv-population`) or decide on fixed PSVs
only (`--psv-fixed-only`).

**Architecture:** A Python bench (`bench/popfix/`) downloads the Hallast Y assemblies and builds
`psv_population.tsv` with a k-mer layer (all columns, via the existing `kscan`) and an alignment layer (copies with
unique flanks, via `minimap2`). In Rust, the default read-star two-form path keeps each candidate's genome position per
column; a pure `fixed_only_certificate` re-runs O2's pairwise certificate on columns whose sites are FIXED for both
copies; a process-wide side table (the codebase's pattern for opt-in levers, cf. `CORE_HULLS`) carries the population
table and the tier results to the binary, which writes `<out>.population_tier.tsv`.

**Tech Stack:** Rust (rustle lib `src/family.rs`, binary `src/copy_assign_bin.rs`), Python 3 (pytest), `kscan`
(`bench/psv_sun/kscan.c`), minimap2 2.30, samtools.

**Spec:** `docs/superpowers/specs/2026-10-12-population-fixed-psvs-design.md`

## Global Constraints

- Data dir: `/mnt/c/Users/jfris/Desktop/data/hgsvc_sigY/` (C:); outputs `/mnt/linuxdisk/tmp/popfix/` (small files only; `/mnt/linuxdisk` has ~5.6 GB free).
- Source: `https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/HGSVC3/working/20230412_sigY_assembly_data/assemblies/verkko_1.0_release/chrY/*.chrY.fasta`, md5 from `MANIFEST_20230412_sigY_assembly_data.tsv`.
- k = 31; allele k-mers = the 31 31-mers containing the site; flank k-mers start 40 bp left / right of the site.
- "Half the men" = ceil(N / 2) with N = number of downloaded assemblies.
- Alignment layer: +- 10 kb flanks, `minimap2 -cx asm5 --cs`; ortholog = the single alignment with qs < 10,000 and qe > qlen - 10,000.
- Certificate constants unchanged: alpha 1e-3 (`AssignParams::alpha`), per-column error e/3 with e = `error_rate` (0.003).
- Without `--psv-population` / `--psv-fixed-only` every output must be byte-identical to today (`assignments`, `quant`, `families`, `star_reads` format excepted only when `--dump-star` is on: it gains a `positions` column).
- Build: `CARGO_TARGET_DIR=/mnt/linuxdisk/home/juanfraitu/rustle_target TMPDIR=/mnt/linuxdisk/tmp`, always `--release`, output to a file; every cargo / copy_assign / minimap2-on-genome run under `bash tools/rlock.sh heavy`, small Python under `bash tools/rlock.sh light`.
- Positions are 0-based everywhere (catalog exons, `psv-discover` `genome_pos`, the table, Rust `u64`).

## Review Focus

- Table chrom names that do not match the catalog (`NC_060948.1` vs `chrY`) -> 0 fixed sites, every molecule UNASSIGNED with no error; expected: the binary warns loudly when no table site lies on any catalog copy's chromosome (test in Task 7).
- An off-by-one between the table (0-based) and read-star genome positions -> plausible but wrong matches; expected: a site at a known 0-based PSV position is found for a read aligned over it (test in Task 6).
- A `-` strand copy in the unit (spliced) form -> genome position must come through `exon_map`, which reverses for `-`; expected: correct position for a `-` copy (test in Task 6).
- `--psv-fixed-only` without `--psv-population`, or with a read-star path other than the default two-form genomic one -> silently ignored; expected: exit with an error (test in Task 7).
- A column where a candidate has no genome position (an insertion column, an unaligned stretch) -> must never count as fixed (test in Task 5).

---

### Task 1: Fetch the population Y assemblies

**Files:**
- Create: `bench/popfix/fetch.sh`
- Create (data): `/mnt/c/Users/jfris/Desktop/data/hgsvc_sigY/SOURCE.md` (written by the script)

**Interfaces:**
- Produces: `/mnt/c/Users/jfris/Desktop/data/hgsvc_sigY/*.chrY.fasta` (+ `.fai`), `SOURCE.md` listing N and sample ids.

- [ ] **Step 1: Write the script**

```bash
#!/bin/bash
# Fetch the Hallast et al. 2023 (HGSVC) per-man chrY assemblies and verify md5 (spec 2026-10-12-population-fixed-psvs §1).
#   fetch.sh OUTDIR
set -euo pipefail
OUT=${1:?OUTDIR}
BASE=https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/HGSVC3/working/20230412_sigY_assembly_data
mkdir -p "$OUT"
curl -s --fail "$BASE/MANIFEST_20230412_sigY_assembly_data.tsv" -o "$OUT/MANIFEST.tsv"
awk -F'\t' '$1 ~ /verkko_1.0_release\/chrY\// && $1 !~ /\.fai$/ {print $1"\t"$3}' "$OUT/MANIFEST.tsv" > "$OUT/want.tsv"
while IFS=$'\t' read -r path md5; do
  f=$(basename "$path")
  if [ ! -s "$OUT/$f" ] || [ "$(md5sum < "$OUT/$f" | cut -d' ' -f1)" != "$md5" ]; then
    curl -s --fail --retry 3 "$BASE/${path#20230412_sigY_assembly_data/}" -o "$OUT/$f.part" && mv "$OUT/$f.part" "$OUT/$f"
  fi
  got=$(md5sum < "$OUT/$f" | cut -d' ' -f1)
  [ "$got" = "$md5" ] || { echo "md5 mismatch: $f" >&2; exit 1; }
  samtools faidx "$OUT/$f"
done < "$OUT/want.tsv"
n=$(wc -l < "$OUT/want.tsv")
{
  echo "# Hallast et al. 2023 Y assemblies (HGSVC, Verkko 1.0), fetched $(date -I)"
  echo "source: $BASE/assemblies/verkko_1.0_release/chrY/"
  echo "n_assemblies: $n (md5-verified against MANIFEST.tsv)"
  echo "samples:"; cut -f1 "$OUT/want.tsv" | xargs -n1 basename | cut -d. -f1 | sed 's/^/  /'
} > "$OUT/SOURCE.md"
echo "$n assemblies OK in $OUT"
```

- [ ] **Step 2: Run it (network; ~2.8 GB to C:)**

Run: `bash bench/popfix/fetch.sh /mnt/c/Users/jfris/Desktop/data/hgsvc_sigY`
Expected: `48 assemblies OK in ...` (the number in `want.tsv`); re-running downloads nothing.

- [ ] **Step 3: Commit**

```bash
git add bench/popfix/fetch.sh
git commit -m "popfix: fetch the Hallast 2023 Y assemblies (md5-verified)"
```

---

### Task 2: k-mer layer — column status from allele / flank k-mers

**Files:**
- Create: `bench/popfix/kmer_layer.py`
- Test: `bench/popfix/test_kmer_layer.py`

**Interfaces:**
- Consumes: `utilities psv-discover` outputs `PREFIX.psv_pos.tsv` (family_id, col_index, copy_index, genome_pos, base) and `PREFIX.psv_copies.tsv` (family_id, copy_index, copy_tid, alleles, strand); the catalog `copies.tsv` (tid -> chrom); the reference FASTA; `kscan` binary.
- Produces: `columns(pos_tsv, copies_tsv_dump, catalog_tsv) -> dict[(fam, col)] -> list[Site]` with `Site = (copy_index, chrom, pos)`; `kmer_status(columns, genome, men_counts, n_men) -> dict[(fam, col)] -> dict(status, n_cover, n_carry_all)`; CLI writing `OUT/kmer_status.tsv` (family_id, col_index, status, n_men, n_cover, n_carry_all).

- [ ] **Step 1: Write the failing tests**

```python
"""Unit tests for the k-mer layer (run: python3 -m pytest -q bench/popfix)."""
import os, sys
sys.path.insert(0, os.path.dirname(__file__))
import kmer_layer as K

G = {"chrY": "A" * 100 + "ACGTTGCATGCAAGTCCGATAGCTAGGCTTACGATCGTAGCTAGCTTAGCAT" + "C" * 100}


def test_allele_kmers_are_the_31_windows_containing_the_site():
    ks = K.allele_kmers(G, "chrY", 120)
    assert len(ks) == 31
    assert all(len(k) == 31 for k in ks)
    raw = G["chrY"][90:151]
    assert {K.canon(raw[i:i + 31]) for i in range(31)} == set(ks)


def test_flank_kmers_sit_40_bp_from_the_site():
    left, right = K.flank_kmers(G, "chrY", 120)
    assert left == K.canon(G["chrY"][50:81])      # ends at 120 - 40
    assert right == K.canon(G["chrY"][160:191])   # starts at 120 + 40


def test_fixed_needs_half_the_men_covering_and_every_allele_carried():
    sites = {"A": ["a1"], "C": ["c1"]}           # allele -> its sites' allele k-mer sets (ids)
    cover = [True, True, True, False]            # man covers the column
    carry = [{"A", "C"}, {"A", "C"}, {"A", "C"}, set()]
    assert K.decide(sites, cover, carry, n_men=4) == "FIXED"
    carry[1] = {"A"}                             # one covering man lost allele C (homogenised)
    assert K.decide(sites, cover, carry, n_men=4) == "POLYMORPHIC"
    assert K.decide(sites, [True, False, False, False], carry, n_men=4) == "UNCOVERED"
```

- [ ] **Step 2: Run tests to verify they fail**

Run: `cd /mnt/linuxdisk/home/juanfraitu/rustle_m2_soto && python3 -m pytest -q bench/popfix/test_kmer_layer.py`
Expected: FAIL with `ModuleNotFoundError: No module named 'kmer_layer'`.

- [ ] **Step 3: Implement**

```python
#!/usr/bin/env python3
"""k-mer layer of spec 2026-10-12-population-fixed-psvs §2: a PSV column is FIXED when every man whose assembly covers
it still carries every allele of it.

    kmer_layer.py PREFIX CATALOG_COPIES_TSV GENOME_FA ASSEMBLY_DIR KSCAN OUTDIR
PREFIX = `utilities psv-discover` output prefix. Columns with a site off chrY are NOT_ASSESSED.
"""
import collections, csv, glob, math, os, subprocess, sys

K, FLANK_GAP = 31, 40
COMP = str.maketrans("ACGTN", "TGCAN")


def canon(k):
    r = k.translate(COMP)[::-1]
    return k if k <= r else r


def fetch(genome, chrom, s, e):
    """genome: dict chrom -> str (tests) or a callable(chrom, s, e) -> str (CLI)."""
    if callable(genome):
        return genome(chrom, s, e)
    return genome[chrom][max(0, s):e]


def allele_kmers(genome, chrom, pos):
    w = fetch(genome, chrom, pos - (K - 1), pos + K).upper()
    return sorted({canon(w[i:i + K]) for i in range(len(w) - K + 1) if set(w[i:i + K]) <= set("ACGT")})


def flank_kmers(genome, chrom, pos):
    left = fetch(genome, chrom, pos - FLANK_GAP - (K - 1), pos - FLANK_GAP + 1).upper()
    right = fetch(genome, chrom, pos + FLANK_GAP, pos + FLANK_GAP + K).upper()
    return canon(left), canon(right)


def decide(allele_sites, cover, carry, n_men):
    """allele_sites: allele -> sites; cover[m]: man m covers the column; carry[m]: alleles man m carries."""
    covering = [m for m in range(len(cover)) if cover[m]]
    if len(covering) < math.ceil(n_men / 2):
        return "UNCOVERED"
    alleles = set(allele_sites)
    return "FIXED" if all(alleles <= carry[m] for m in covering) else "POLYMORPHIC"


def tsv(path):
    with open(path) as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def columns(prefix, catalog_tsv):
    chrom_of = {r["tid"]: r["chrom"] for r in tsv(catalog_tsv)}
    tid = {(r["family_id"], r["copy_index"]): r["copy_tid"] for r in tsv(prefix + ".psv_copies.tsv")}
    cols = collections.defaultdict(list)
    for r in tsv(prefix + ".psv_pos.tsv"):
        c = chrom_of[tid[(r["family_id"], r["copy_index"])]]
        cols[(r["family_id"], int(r["col_index"]))].append((int(r["copy_index"]), c, int(r["genome_pos"])))
    return cols


def main(prefix, catalog_tsv, genome_fa, asm_dir, kscan, out):
    sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "psv_sun"))
    import psv_sun
    F = psv_sun.Fasta(genome_fa)
    genome = lambda c, s, e: F.fetch(c, s, e)
    os.makedirs(out, exist_ok=True)
    cols = columns(prefix, catalog_tsv)
    plan, want = {}, set()
    for key, sites in cols.items():
        if any(c != "chrY" for _, c, _ in sites):
            plan[key] = None
            continue
        al = {i: allele_kmers(genome, c, p) for i, c, p in sites}
        fl = {i: flank_kmers(genome, c, p) for i, c, p in sites}
        base = {i: F.fetch(c, p, p + 1).upper() for i, c, p in sites}
        plan[key] = (al, fl, base)
        for ks in al.values():
            want.update(ks)
        for l, r in fl.values():
            want.update((l, r))
    kfile = os.path.join(out, "kmers.txt")
    with open(kfile, "w") as fh:
        fh.write("\n".join(sorted(want)) + "\n")
    men = sorted(glob.glob(os.path.join(asm_dir, "*.chrY.fasta")))
    present = []
    for fa in men:
        res = subprocess.run([kscan, str(K), kfile, fa], capture_output=True, text=True, check=True).stdout
        present.append({k for k, n in (l.split("\t") for l in res.splitlines()) if int(n) > 0})
    with open(os.path.join(out, "kmer_status.tsv"), "w") as fh:
        fh.write("family_id\tcol_index\tstatus\tn_men\tn_cover\tn_carry_all\n")
        for (fam, col), p in sorted(plan.items()):
            if p is None:
                fh.write(f"{fam}\t{col}\tNOT_ASSESSED\t{len(men)}\t0\t0\n")
                continue
            al, fl, base = p
            allele_sites = collections.defaultdict(list)
            for i, b in base.items():
                allele_sites[b].append(i)
            cover = [all(l in pm and r in pm for l, r in fl.values()) for pm in present]
            carry = [{b for b, ids in allele_sites.items() if any(any(k in pm for k in al[i]) for i in ids)}
                     for pm in present]
            st = decide(allele_sites, cover, carry, len(men))
            n_cover = sum(cover)
            n_all = sum(1 for m in range(len(men)) if cover[m] and set(allele_sites) <= carry[m])
            fh.write(f"{fam}\t{col}\t{st}\t{len(men)}\t{n_cover}\t{n_all}\n")
    print(f"{len(plan)} columns, {len(men)} assemblies -> {out}/kmer_status.tsv")


if __name__ == "__main__":
    main(*sys.argv[1:7])
```

- [ ] **Step 4: Run tests to verify they pass**

Run: `python3 -m pytest -q bench/popfix/test_kmer_layer.py`
Expected: `3 passed`.

- [ ] **Step 5: Commit**

```bash
git add bench/popfix/kmer_layer.py bench/popfix/test_kmer_layer.py
git commit -m "popfix: k-mer layer (allele/flank 31-mers, FIXED iff every covering man carries every allele)"
```

---

### Task 3: alignment layer — copy-level check where flanks are unique

**Files:**
- Create: `bench/popfix/align_layer.py`
- Test: `bench/popfix/test_align_layer.py`

**Interfaces:**
- Consumes: `kmer_layer.columns`, `kmer_layer.canon`, `psv_sun.Fasta`, `psv_sun.run_paf`, `psv_sun.anchor_pairs`, `kscan`.
- Produces: `ortholog_base(paf_recs, qlen, flank, site_offset) -> str | None` and CLI writing `OUT/align_status.tsv` (family_id, col_index, status ∈ {FIXED, POLYMORPHIC, UNCOVERED, NO_ALIGN}, n_men, n_ortholog, n_mismatch). NO_ALIGN = some site of the column lacks unique flanks (the k-mer verdict stands).

- [ ] **Step 1: Write the failing test**

```python
"""Unit test for the alignment layer's ortholog base read-through (needs minimap2)."""
import os, random, sys, tempfile
sys.path.insert(0, os.path.dirname(__file__))
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "psv_sun"))
import align_layer as A
import psv_sun as P


def test_ortholog_base_reads_the_mans_base_at_the_site():
    rng = random.Random(1)
    q = "".join(rng.choice("ACGT") for _ in range(30_000))      # locus +- 10 kb, site at 15,000
    man = list(rng.choice("ACGT") for _ in range(5_000)) + list(q) + list(rng.choice("ACGT") for _ in range(5_000))
    alt = "A" if q[15_000] != "A" else "C"
    man[5_000 + 15_000] = alt                                    # this man carries another base at the site
    with tempfile.TemporaryDirectory() as d:
        open(f"{d}/q.fa", "w").write(f">q\n{q}\n")
        open(f"{d}/m.fa", "w").write(">m\n" + "".join(man) + "\n")
        recs = P.run_paf(["minimap2", "-cx", "asm5", "--cs", f"{d}/m.fa", f"{d}/q.fa"])
        tseq = "".join(man)
    assert A.ortholog_base(recs, len(q), 10_000, 15_000, lambda name, s, e: tseq[s:e]) == alt


def test_no_single_flank_spanning_alignment_is_no_ortholog():
    assert A.ortholog_base([], 30_000, 10_000, 15_000, lambda n, s, e: "") is None
```

- [ ] **Step 2: Run tests to verify they fail**

Run: `python3 -m pytest -q bench/popfix/test_align_layer.py`
Expected: FAIL with `ModuleNotFoundError: No module named 'align_layer'`.

- [ ] **Step 3: Implement**

```python
#!/usr/bin/env python3
"""Alignment layer of spec 2026-10-12-population-fixed-psvs §2: for copies whose +- 10 kb flanks are unique in the
reference, each man's ortholog is the single alignment of locus + flanks spanning both flanks; his base at every PSV
position is read through it. A column gets an alignment verdict only if ALL its sites qualify.

    align_layer.py PREFIX CATALOG_COPIES_TSV GENOME_FA ASSEMBLY_DIR KSCAN OUTDIR
"""
import collections, glob, math, os, subprocess, sys
sys.path.insert(0, os.path.dirname(__file__))
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "psv_sun"))
import kmer_layer as KL
import psv_sun as P

FLANK, STEP = 10_000, 100


def ortholog_base(recs, qlen, flank, site_offset, fetch_target):
    """The man's base at query offset `site_offset`, through the single alignment spanning both flanks, else None."""
    span = [r for r in recs if r["tp"] == "P" and r["qs"] < flank and r["qe"] > qlen - flank]
    if len(span) != 1:
        return None
    r = span[0]
    t = P.paf_qmap(r).get(site_offset)
    if t is None:
        return None
    b = fetch_target(r["tname"], t, t + 1).upper()
    return b.translate(KL.COMP) if r["strand"] == "-" else b


def unique_flanks(F, kscan, genome_fa, sites, workdir):
    """site -> True if each +-10 kb flank has a 31-mer occurring once in the reference (sampled every STEP bp)."""
    probe = {}
    for (c, p) in sites:
        for side, (s, e) in (("L", (p - FLANK, p - 1000)), ("R", (p + 1000, p + FLANK))):
            seq = F.fetch(c, max(0, s), e).upper()
            probe[(c, p, side)] = [KL.canon(seq[i:i + KL.K]) for i in range(0, len(seq) - KL.K, STEP)
                                   if set(seq[i:i + KL.K]) <= set("ACGT")]
    kf = os.path.join(workdir, "flank_kmers.txt")
    with open(kf, "w") as fh:
        fh.write("\n".join(sorted({k for v in probe.values() for k in v})) + "\n")
    res = subprocess.run([kscan, str(KL.K), kf, genome_fa], capture_output=True, text=True, check=True).stdout
    once = {k for k, n in (l.split("\t") for l in res.splitlines()) if int(n) == 1}
    return {(c, p): any(k in once for k in probe[(c, p, "L")]) and any(k in once for k in probe[(c, p, "R")])
            for (c, p) in sites}


def main(prefix, catalog_tsv, genome_fa, asm_dir, kscan, out):
    os.makedirs(out, exist_ok=True)
    F = P.Fasta(genome_fa)
    cols = KL.columns(prefix, catalog_tsv)
    ycols = {k: v for k, v in cols.items() if all(c == "chrY" for _, c, _ in v)}
    sites = sorted({(c, p) for v in ycols.values() for _, c, p in v})
    ok = unique_flanks(F, kscan, genome_fa, sites, out)
    qual = [s for s in sites if ok[s]]
    qfa = os.path.join(out, "loci.fa")
    with open(qfa, "w") as fh:
        for c, p in qual:
            fh.write(f">{c}:{p}\n{F.fetch(c, p - FLANK, p + FLANK + 1)}\n")
    men = sorted(glob.glob(os.path.join(asm_dir, "*.chrY.fasta")))
    base = collections.defaultdict(dict)   # site -> man -> base or None
    for mi, fa in enumerate(men):
        T = P.Fasta(fa)
        recs = collections.defaultdict(list)
        for r in P.run_paf(["minimap2", "-cx", "asm5", "--cs", "-t", "2", fa, qfa]):
            recs[r["qname"]].append(r)
        for c, p in qual:
            base[(c, p)][mi] = ortholog_base(recs.get(f"{c}:{p}", []), 2 * FLANK + 1, FLANK, FLANK,
                                             lambda n, s, e: T.fetch(n, s, e))
    half = math.ceil(len(men) / 2)
    with open(os.path.join(out, "align_status.tsv"), "w") as fh:
        fh.write("family_id\tcol_index\tstatus\tn_men\tn_ortholog\tn_mismatch\n")
        for (fam, col), v in sorted(ycols.items()):
            if not all(ok[(c, p)] for _, c, p in v):
                fh.write(f"{fam}\t{col}\tNO_ALIGN\t{len(men)}\t0\t0\n")
                continue
            have = [m for m in range(len(men)) if all(base[(c, p)].get(m) for _, c, p in v)]
            mism = sum(1 for m in have if any(base[(c, p)][m] != F.fetch(c, p, p + 1).upper() for _, c, p in v))
            st = "UNCOVERED" if len(have) < half else ("FIXED" if mism == 0 else "POLYMORPHIC")
            fh.write(f"{fam}\t{col}\t{st}\t{len(men)}\t{len(have)}\t{mism}\n")
    print(f"{len(qual)}/{len(sites)} sites with unique flanks, {len(men)} assemblies -> {out}/align_status.tsv")


if __name__ == "__main__":
    main(*sys.argv[1:7])
```

- [ ] **Step 4: Run tests to verify they pass**

Run: `python3 -m pytest -q bench/popfix/test_align_layer.py`
Expected: `2 passed`.

- [ ] **Step 5: Commit**

```bash
git add bench/popfix/align_layer.py bench/popfix/test_align_layer.py
git commit -m "popfix: alignment layer (unique-flank copies, ortholog base read-through)"
```

---

### Task 4: build the population table for the A119b chrY families

**Files:**
- Create: `bench/popfix/table.py`
- Test: `bench/popfix/test_table.py`

**Interfaces:**
- Consumes: `OUT/kmer_status.tsv`, `OUT/align_status.tsv`, `PREFIX.psv_pos.tsv`, `PREFIX.psv_copies.tsv`, catalog copies table.
- Produces: `psv_population.tsv` — header `chrom\tpos\tfamily_id\tcol_index\tstatus\tmethod\tn_men\tn_cover\tn_mismatch`, one row per (chrom, pos) site, 0-based `pos`; `combine(kmer_row, align_row) -> (status, method)`.

- [ ] **Step 1: Write the failing test**

```python
import os, sys
sys.path.insert(0, os.path.dirname(__file__))
import table as T


def test_alignment_verdict_overrides_kmer_verdict():
    assert T.combine({"status": "FIXED"}, {"status": "POLYMORPHIC"}) == ("POLYMORPHIC", "A")
    assert T.combine({"status": "POLYMORPHIC"}, {"status": "FIXED"}) == ("FIXED", "A")


def test_kmer_verdict_stands_without_alignment():
    assert T.combine({"status": "FIXED"}, {"status": "NO_ALIGN"}) == ("FIXED", "K")
    assert T.combine({"status": "FIXED"}, None) == ("FIXED", "K")
    assert T.combine({"status": "NOT_ASSESSED"}, None) == ("NOT_ASSESSED", "-")


def test_alignment_uncovered_falls_back_to_kmer():
    assert T.combine({"status": "FIXED"}, {"status": "UNCOVERED"}) == ("FIXED", "K")
```

- [ ] **Step 2: Run tests to verify they fail**

Run: `python3 -m pytest -q bench/popfix/test_table.py`
Expected: FAIL with `ModuleNotFoundError: No module named 'table'`.

- [ ] **Step 3: Implement**

```python
#!/usr/bin/env python3
"""psv_population.tsv (spec §2): per (chrom, pos) site, the alignment verdict where it exists, else the k-mer verdict.

    table.py PREFIX CATALOG_COPIES_TSV KMER_STATUS ALIGN_STATUS OUT.tsv
"""
import os, sys
sys.path.insert(0, os.path.dirname(__file__))
import kmer_layer as KL


def combine(krow, arow):
    if krow["status"] == "NOT_ASSESSED":
        return "NOT_ASSESSED", "-"
    if arow is not None and arow["status"] in ("FIXED", "POLYMORPHIC"):
        return arow["status"], "A"
    return krow["status"], "K"


def main(prefix, catalog_tsv, kmer_tsv, align_tsv, out):
    km = {(r["family_id"], int(r["col_index"])): r for r in KL.tsv(kmer_tsv)}
    al = {(r["family_id"], int(r["col_index"])): r for r in KL.tsv(align_tsv)} if os.path.exists(align_tsv) else {}
    cols = KL.columns(prefix, catalog_tsv)
    rank = {"POLYMORPHIC": 0, "UNCOVERED": 1, "NOT_ASSESSED": 2, "FIXED": 3}
    site = {}
    for key, sites in cols.items():
        st, meth = combine(km[key], al.get(key))
        a = al.get(key, {})
        for _, c, p in sites:
            row = (st, meth, key, km[key]["n_men"], km[key]["n_cover"], a.get("n_mismatch", "0"))
            if (c, p) not in site or rank[st] < rank[site[(c, p)][0]]:   # a site in two columns: the stricter wins
                site[(c, p)] = row
    with open(out, "w") as fh:
        fh.write("chrom\tpos\tfamily_id\tcol_index\tstatus\tmethod\tn_men\tn_cover\tn_mismatch\n")
        for (c, p), (st, meth, (fam, col), n, nc, nm) in sorted(site.items()):
            fh.write(f"{c}\t{p}\t{fam}\t{col}\t{st}\t{meth}\t{n}\t{nc}\t{nm}\n")
    print(f"{len(site)} sites -> {out}")


if __name__ == "__main__":
    main(*sys.argv[1:6])
```

- [ ] **Step 4: Run tests to verify they pass**

Run: `python3 -m pytest -q bench/popfix`
Expected: `8 passed`.

- [ ] **Step 5: Build the table on real data**

```bash
cd /mnt/linuxdisk/home/juanfraitu/rustle_m2_soto
O=/mnt/linuxdisk/tmp/popfix; H=/mnt/linuxdisk/tmp/rustle_figures/runs/human_A119b/human_A119b
G=/mnt/linuxdisk/home/juanfraitu/winloci_data/chm13v2.0.fa; A=/mnt/c/Users/jfris/Desktop/data/hgsvc_sigY
U=/mnt/linuxdisk/home/juanfraitu/rustle_target/release/utilities; KS=/mnt/linuxdisk/tmp/psv_sun/bin/kscan
mkdir -p $O
awk -F'\t' 'NR>1 && $4=="chrY"{print $1}' $H.fam.copies.tsv | sort -u > $O/yfams.txt
bash tools/rlock.sh light $U psv-discover --families $H.fam.copies.tsv --copies-fa $H.fam.copies.fa --only $O/yfams.txt --out $O/y --threads 2
bash tools/rlock.sh heavy python3 bench/popfix/kmer_layer.py $O/y $H.fam.copies.tsv $G $A $KS $O
bash tools/rlock.sh heavy python3 bench/popfix/align_layer.py $O/y $H.fam.copies.tsv $G $A $KS $O
python3 bench/popfix/table.py $O/y $H.fam.copies.tsv $O/kmer_status.tsv $O/align_status.tsv $O/psv_population.tsv
cut -f5,6 $O/psv_population.tsv | sort | uniq -c
```
Expected: a table with counts per status/method; every row's `pos` is a `psv-discover` `genome_pos`.

- [ ] **Step 6: Commit**

```bash
git add bench/popfix/table.py bench/popfix/test_table.py
git commit -m "popfix: psv_population.tsv (alignment verdict overrides k-mer verdict; stricter status wins per site)"
```

---

### Task 5: Rust — `fixed_only_certificate` (pure) and the population side table

**Files:**
- Modify: `src/family.rs` (module `copy_assign::copy_assign_pipeline`, next to `copy_pair_significance` users; tests in its `tests` module)

**Interfaces:**
- Produces:
  - `pub struct PsvPopulation { pub fixed: DetHashSet<(String, u64)>, pub strict: bool }`
  - `pub fn parse_psv_population(text: &str, strict: bool) -> anyhow::Result<PsvPopulation>` (FIXED rows only; 0-based `pos`)
  - `pub fn register_psv_population(t: PsvPopulation)` / `pub fn psv_population() -> Option<&'static PsvPopulation>`
  - `pub struct FixedCert { pub status: AssignStatus, pub best: usize, pub p_read: f64, pub min_p: f64, pub margin: f64, pub n_dec: usize }`
  - `pub fn fixed_only_certificate(obs: &[Option<u8>], profiles: &[CopyProfile], fixed: &[Vec<bool>], partner: &[bool], p: &AssignParams) -> FixedCert` — `fixed[k][j]` = candidate k's site at column j is FIXED; `best` is a candidate index (into `profiles`).

- [ ] **Step 1: Write the failing tests** (in `copy_assign_pipeline::tests`)

```rust
            fn prof(alleles: &[u8]) -> CopyProfile {
                CopyProfile { copy_id: 0, alleles: alleles.iter().map(|&b| Some(b)).collect(), junctions: vec![] }
            }

            #[test]
            fn fixed_only_certificate_assigns_on_fixed_columns_and_ignores_the_rest() {
                // 3 columns; the read matches copy 0 at all; column 2 is NOT fixed and is where copy 1 would differ
                let obs = vec![Some(b'A'), Some(b'C'), Some(b'G')];
                let profiles = vec![prof(b"ACG"), prof(b"TGT")];
                let all_fixed = vec![vec![true; 3], vec![true; 3]];
                let c = fixed_only_certificate(&obs, &profiles, &all_fixed, &[false, false], &AssignParams::default());
                assert_eq!((c.status, c.best, c.n_dec), (AssignStatus::Assigned, 0, 3));
                let one_fixed = vec![vec![true, false, false], vec![true, false, false]];
                let c = fixed_only_certificate(&obs, &profiles, &one_fixed, &[false, false], &AssignParams::default());
                assert_eq!(c.n_dec, 1, "only the fixed column counts");
                assert_eq!(c.status, AssignStatus::Ambiguous, "one column never assigns (m = 1)");
            }

            #[test]
            fn a_column_without_a_fixed_site_on_both_copies_never_counts() {
                // Review Focus: a candidate with no genome position (insertion column) is never fixed
                let obs = vec![Some(b'A'), Some(b'C')];
                let profiles = vec![prof(b"AC"), prof(b"TG")];
                let fixed = vec![vec![true, true], vec![false, false]];
                let c = fixed_only_certificate(&obs, &profiles, &fixed, &[false, false], &AssignParams::default());
                assert_eq!((c.status, c.n_dec), (AssignStatus::Tied, 0));
            }

            #[test]
            fn psv_population_keeps_fixed_rows_only() {
                let t = "chrom\tpos\tfamily_id\tcol_index\tstatus\tmethod\tn_men\tn_cover\tn_mismatch\n\
                         chrY\t100\tF\t0\tFIXED\tK\t48\t40\t0\nchrY\t200\tF\t1\tPOLYMORPHIC\tA\t48\t44\t3\n";
                let pop = parse_psv_population(t, false).unwrap();
                assert!(pop.fixed.contains(&("chrY".to_string(), 100)));
                assert!(!pop.fixed.contains(&("chrY".to_string(), 200)));
                assert!(parse_psv_population("chrom\tpos\n", false).is_err(), "no status column");
            }
```

- [ ] **Step 2: Run tests to verify they fail**

Run: `CARGO_TARGET_DIR=/mnt/linuxdisk/home/juanfraitu/rustle_target TMPDIR=/mnt/linuxdisk/tmp bash tools/rlock.sh heavy bash -c "cargo test --release --lib fixed_only_certificate a_column_without psv_population_keeps > /mnt/linuxdisk/tmp/popfix_t5_red.log 2>&1"; grep -E '^error\[' /mnt/linuxdisk/tmp/popfix_t5_red.log | sort | uniq -c`
Expected: `cannot find function fixed_only_certificate / parse_psv_population` compile errors.

- [ ] **Step 3: Implement** (place after `copy_pair_significance`'s module neighbours, inside `copy_assign_pipeline`)

```rust
        /// spec 2026-10-12-population-fixed-psvs §3: the population table — FIXED sites only, 0-based (chrom, pos).
        pub struct PsvPopulation {
            pub fixed: DetHashSet<(String, u64)>,
            /// `--psv-fixed-only`: the fixed-only certificate IS the decision
            pub strict: bool,
        }

        pub fn parse_psv_population(text: &str, strict: bool) -> anyhow::Result<PsvPopulation> {
            let mut lines = text.lines();
            let head: Vec<&str> = lines.next().unwrap_or("").split('\t').collect();
            let col = |n: &str| head.iter().position(|h| *h == n).ok_or_else(|| anyhow::anyhow!("--psv-population: no `{n}` column"));
            let (ic, ip, is) = (col("chrom")?, col("pos")?, col("status")?);
            let mut fixed = DetHashSet::default();
            for l in lines.filter(|l| !l.is_empty()) {
                let f: Vec<&str> = l.split('\t').collect();
                if f.get(is) == Some(&"FIXED") {
                    fixed.insert((f[ic].to_string(), f[ip].parse::<u64>()?));
                }
            }
            Ok(PsvPopulation { fixed, strict })
        }

        static PSV_POPULATION: std::sync::OnceLock<PsvPopulation> = std::sync::OnceLock::new();
        pub fn register_psv_population(t: PsvPopulation) {
            let _ = PSV_POPULATION.set(t);
        }
        pub fn psv_population() -> Option<&'static PsvPopulation> {
            PSV_POPULATION.get()
        }

        /// The read-star pairwise certificate (§6fc) on the columns where BOTH copies' sites are population-FIXED
        /// (`fixed[k][j]`): best = maximin pairwise LLR over non-partner candidates (ties: fixed-column PSV score, then
        /// lower index); Tied if some rival shares no fixed distinguishing column; Assigned iff every rival's
        /// `copy_pair_significance` p < alpha / #rivals and its LLR > 0; else Ambiguous. Junctions are not used.
        pub struct FixedCert {
            pub status: AssignStatus,
            pub best: usize,
            pub p_read: f64,
            pub min_p: f64,
            pub margin: f64,
            pub n_dec: usize,
        }

        pub fn fixed_only_certificate(
            obs: &[Option<u8>],
            profiles: &[CopyProfile],
            fixed: &[Vec<bool>],
            partner: &[bool],
            p: &AssignParams,
        ) -> FixedCert {
            let n = profiles.len();
            let ks: Vec<usize> = (0..n).filter(|&k| !partner.get(k).copied().unwrap_or(false)).collect();
            let lr = ((1.0 - p.error_rate) / (p.error_rate / 3.0)).ln();
            let masked = |a: usize, b: usize| -> Vec<Option<u8>> {
                (0..obs.len())
                    .map(|j| if fixed[a].get(j).copied().unwrap_or(false) && fixed[b].get(j).copied().unwrap_or(false) { obs[j] } else { None })
                    .collect()
            };
            let pair = |a: usize, b: usize| -> (usize, f64) {
                let o = masked(a, b);
                let (mut dec, mut llr) = (0usize, 0.0f64);
                for j in 0..o.len() {
                    if let (Some(x), Some(ba), Some(ca)) = (o[j], profiles[a].alleles[j], profiles[b].alleles[j]) {
                        if ba != ca {
                            dec += 1;
                            if x == ba { llr += lr } else if x == ca { llr -= lr }
                        }
                    }
                }
                (dec, llr)
            };
            let score = |a: usize| -> i64 {
                let o: Vec<Option<u8>> = (0..obs.len()).map(|j| if fixed[a].get(j).copied().unwrap_or(false) { obs[j] } else { None }).collect();
                profiles[a].alleles.iter().zip(o.iter()).filter_map(|(&c, &x)| match (c, x) {
                    (Some(c), Some(x)) if c == x => Some(1i64),
                    (Some(_), Some(_)) => Some(-1),
                    _ => None,
                }).sum()
            };
            let maximin = |a: usize| ks.iter().filter(|&&b| b != a).map(|&b| pair(a, b).1).fold(f64::INFINITY, f64::min);
            let best = *ks
                .iter()
                .max_by(|&&a, &&b| {
                    maximin(a).partial_cmp(&maximin(b)).unwrap_or(std::cmp::Ordering::Equal)
                        .then(score(a).cmp(&score(b)))
                        .then(b.cmp(&a)) // lower index wins ties
                })
                .unwrap_or(&0);
            let rivals: Vec<usize> = ks.iter().copied().filter(|&b| b != best).collect();
            if rivals.is_empty() {
                return FixedCert { status: AssignStatus::Tied, best, p_read: 1.0, min_p: 1.0, margin: 0.0, n_dec: 0 };
            }
            let thr = p.alpha / rivals.len() as f64;
            let (mut p_read, mut min_p, mut margin, mut n_dec, mut k0) = (0.0f64, 1.0f64, f64::INFINITY, usize::MAX, false);
            for &b in &rivals {
                let feats = ReadFeatures { psv_obs: masked(best, b), psv_qual: Vec::new(), junctions: Vec::new() };
                let (pk, mp, _) = copy_pair_significance(&feats, &profiles[best], &profiles[b], p, &[]);
                let (dec, llr) = pair(best, b);
                k0 |= dec == 0;
                p_read = p_read.max(pk);
                min_p = min_p.min(mp);
                margin = margin.min(llr);
                n_dec = n_dec.min(dec);
            }
            let status = if k0 {
                AssignStatus::Tied
            } else if p_read < thr && margin > 0.0 {
                AssignStatus::Assigned
            } else {
                AssignStatus::Ambiguous
            };
            FixedCert { status, best, p_read, min_p, margin, n_dec }
        }
```

- [ ] **Step 4: Run tests to verify they pass**

Run: `CARGO_TARGET_DIR=/mnt/linuxdisk/home/juanfraitu/rustle_target TMPDIR=/mnt/linuxdisk/tmp bash tools/rlock.sh heavy bash -c "cargo test --release --lib > /mnt/linuxdisk/tmp/popfix_t5.log 2>&1"; grep -E 'test result|FAILED' /mnt/linuxdisk/tmp/popfix_t5.log`
Expected: `test result: ok.` with 3 more tests than before (1,036 passed), 0 failed.

- [ ] **Step 5: Commit**

```bash
git add src/family.rs
git commit -m "popfix: fixed_only_certificate + population side table (FIXED sites only)"
```

---

### Task 6: Rust — genome positions of read-star columns, tier side table, strict override

**Files:**
- Modify: `src/family.rs` — the two-form branch `let alns: Vec<ReadStarObs> = if genomic_spans.is_some() && p.read_star_two_form {` (≈ line 24216), the consumer `.zip(alns.par_iter())` (≈ line 24282), the status decision `let status = if single || k0 {` (≈ line 24580), `StarProof` (≈ line 21265).

**Interfaces:**
- Consumes: Task 5's `psv_population()`, `fixed_only_certificate`, `FixedCert`.
- Produces:
  - `pub fn star_column_genome_positions(cols: &[u32], cand: &[usize], per_pos: &[Option<Vec<Option<u32>>>], to_genome: &dyn Fn(usize, u32) -> Option<u64>) -> Vec<Vec<Option<u64>>>` — `[k][j]`, candidate order.
  - `pub struct PopulationTier { pub status_all: AssignStatus, pub copy_all: usize, pub status_fixed: AssignStatus, pub copy_fixed: usize, pub n_dec_fixed: usize, pub cand_tids: Vec<String> }` registered per molecule by `register_population_tier(read_name, tier)`, drained by `take_population_tiers() -> DetHashMap<String, Vec<PopulationTier>>`.
  - `StarProof` gains `pub gpos: Vec<Vec<Option<u64>>>` (candidate order, per column).

- [ ] **Step 1: Write the failing tests**

```rust
            #[test]
            fn star_column_positions_map_genomic_and_unit_forms_both_strands() {
                // copy 0: genomic form, locus starts at 1,000 -> target offset 5 is genome 1,005
                // copy 1: unit form of a '-' copy with exons [2,000,2,010) -> unit offset 0 is genome 2,009 (exon_map)
                let c1 = copy_tx("B", 2_000, 2_010, '-', &[], vec![b'A'; 10]);
                let emap1 = exon_map(&c1);
                let per_pos = vec![Some(vec![None, Some(5)]), Some(vec![None, Some(0)])];
                let to_g = |ci: usize, off: u32| -> Option<u64> {
                    if ci == 0 { Some(1_000 + off as u64) } else { emap1.get(off as usize).copied() }
                };
                let g = star_column_genome_positions(&[1], &[0, 1], &per_pos, &to_g);
                assert_eq!(g, vec![vec![Some(1_005)], vec![Some(2_009)]]);
                // a read position a candidate does not cover has no position (Review Focus)
                let g0 = star_column_genome_positions(&[0], &[0, 1], &per_pos, &to_g);
                assert_eq!(g0, vec![vec![None], vec![None]]);
            }

            #[test]
            fn a_site_at_a_known_zero_based_psv_position_is_found() {
                // Review Focus: table pos and read-star positions are both 0-based
                let pop = parse_psv_population(
                    "chrom\tpos\tfamily_id\tcol_index\tstatus\tmethod\tn_men\tn_cover\tn_mismatch\nchrY\t1005\tF\t0\tFIXED\tK\t48\t40\t0\n",
                    false,
                ).unwrap();
                let per_pos = vec![Some(vec![None, Some(5)])];
                let g = star_column_genome_positions(&[1], &[0], &per_pos, &|_, off| Some(1_000 + off as u64));
                assert!(pop.fixed.contains(&("chrY".to_string(), g[0][0].unwrap())));
            }
```

- [ ] **Step 2: Run tests to verify they fail**

Run: `CARGO_TARGET_DIR=/mnt/linuxdisk/home/juanfraitu/rustle_target TMPDIR=/mnt/linuxdisk/tmp bash tools/rlock.sh heavy bash -c "cargo test --release --lib star_column_positions a_site_at_a_known > /mnt/linuxdisk/tmp/popfix_t6_red.log 2>&1"; grep -E '^error\[' /mnt/linuxdisk/tmp/popfix_t6_red.log | sort | uniq -c`
Expected: `cannot find function star_column_genome_positions`.

- [ ] **Step 3: Implement the helper and the side table**

```rust
        /// spec §3: candidate `cand[k]`'s genome position at each read-star column (`cols` = read positions):
        /// `per_pos[copy][read pos]` = the target offset, `to_genome(copy, offset)` maps it (genomic form: locus start +
        /// offset; unit form: `exon_map`). `None` where the candidate has no base at that read position.
        pub fn star_column_genome_positions(
            cols: &[u32],
            cand: &[usize],
            per_pos: &[Option<Vec<Option<u32>>>],
            to_genome: &dyn Fn(usize, u32) -> Option<u64>,
        ) -> Vec<Vec<Option<u64>>> {
            cand.iter()
                .map(|&ci| {
                    cols.iter()
                        .map(|&r| {
                            per_pos.get(ci).and_then(|v| v.as_ref()).and_then(|v| v.get(r as usize).copied().flatten())
                                .and_then(|off| to_genome(ci, off))
                        })
                        .collect()
                })
                .collect()
        }

        /// spec §3: one contested molecule's portability tier (`--psv-population`), keyed like `StarProof`.
        pub struct PopulationTier {
            pub status_all: AssignStatus,
            pub copy_all: usize,
            pub status_fixed: AssignStatus,
            pub copy_fixed: usize,
            pub n_dec_fixed: usize,
            pub cand_tids: Vec<String>,
        }
        static POPULATION_TIERS: std::sync::OnceLock<std::sync::Mutex<DetHashMap<String, Vec<PopulationTier>>>> =
            std::sync::OnceLock::new();
        pub fn register_population_tier(read_name: &str, t: PopulationTier) {
            POPULATION_TIERS.get_or_init(Default::default).lock().unwrap().entry(read_name.to_string()).or_default().push(t);
        }
        pub fn take_population_tiers() -> DetHashMap<String, Vec<PopulationTier>> {
            POPULATION_TIERS.get().map(|m| std::mem::take(&mut *m.lock().unwrap())).unwrap_or_default()
        }
```

- [ ] **Step 4: Thread the positions through the two-form branch**

Replace the two-form branch's `.map(|(ri, (g, u))| { ... build_star_obs(...) })` and `.collect()` so it also returns positions (computed only when `psv_population()` is registered or `p.dump_star`), and collect into `alns` + `star_gpos`:

```rust
                let want_gpos = psv_population().is_some() || p.dump_star;
                let unit_maps: Vec<Vec<u64>> = if want_gpos { copies.iter().map(|c| exon_map(c)).collect() } else { Vec::new() };
                let (alns, star_gpos): (Vec<ReadStarObs>, Vec<Vec<Vec<Option<u64>>>>) = if genomic_spans.is_some() && p.read_star_two_form {
                    // ⭐ §6fp: both forms, the better one per (molecule, candidate)
                    let unit_seqs: Vec<&[u8]> = copies.iter().map(|c| c.seq.as_slice()).collect();
                    let raw_g = read_star_raw(&copy_seqs, &read_seqs, preset, &windows);
                    let raw_u = read_star_raw(&unit_seqs, &read_seqs, "map-hifi", &[]);
                    raw_g
                        .into_iter()
                        .zip(raw_u)
                        .enumerate()
                        .map(|(ri, (g, u))| {
                            let ((pc, pp, ed, ev), form) = merge_star_raw(g, u);
                            let o = build_star_obs(read_seqs[ri], &pc, &pp, &ed, &[], &ev, Some(&form), indel_min_len);
                            let gp = if want_gpos {
                                let to_g = |ci: usize, off: u32| -> Option<u64> {
                                    if form[ci] {
                                        unit_maps[ci].get(off as usize).copied() // unit form: spliced offset -> genome
                                    } else {
                                        Some(target_of(copies[ci]).0 + off as u64) // genomic form: locus start + offset
                                    }
                                };
                                star_column_genome_positions(&o.5, &o.0, &pp, &to_g)
                            } else {
                                Vec::new()
                            };
                            (o, gp)
                        })
                        .unzip()
                } else {
                    let o = read_star_observations(
                        &copy_seqs,
                        &read_seqs,
                        if genomic_spans.is_some() { &[] } else { &copy_bounds },
                        preset,
                        &windows,
                        indel_min_len,
                    );
                    let n = o.len();
                    (o, vec![Vec::new(); n]) // no positions off the default path (the binary refuses the flags there)
                };
```

In the consumer, change `.zip(alns.par_iter())` to `.zip(alns.par_iter().zip(star_gpos.par_iter()))` and the closure pattern to `|(&ri, ((cand, obs, profiles, unit_edits, covered, star_cols), gpos))|`.

- [ ] **Step 5: Tier + strict override at the decision**

Directly after `let status = if single || k0 { ... };` insert (and make `bk`, `status`, `p_read`, `min_p`, `margin_min`, `n_dec_min` mutable at their `let`s — they are only reassigned here):

```rust
                if let Some(pop) = psv_population() {
                    if !gpos.is_empty() {
                        let fixed: Vec<Vec<bool>> = gpos
                            .iter()
                            .enumerate()
                            .map(|(k, v)| {
                                let chrom = &copies[cand[k]].chrom;
                                v.iter().map(|g| g.map_or(false, |g| pop.fixed.contains(&(chrom.clone(), g)))).collect()
                            })
                            .collect();
                        let partner: Vec<bool> = (0..cand.len()).map(|k| partner_k.contains(&k)).collect();
                        let fc = fixed_only_certificate(&feats.psv_obs, &profiles, &fixed, &partner, p);
                        register_population_tier(
                            &names[ri],
                            PopulationTier {
                                status_all: status,
                                copy_all: cand[bk],
                                status_fixed: fc.status,
                                copy_fixed: cand[fc.best],
                                n_dec_fixed: if fc.n_dec == usize::MAX { 0 } else { fc.n_dec },
                                cand_tids: cand.iter().map(|&k| copies[k].tid.clone()).collect(),
                            },
                        );
                        if pop.strict {
                            // --psv-fixed-only: the fixed-only certificate IS the decision (later gates still apply)
                            status = fc.status;
                            bk = fc.best;
                            p_read = fc.p_read;
                            min_p = fc.min_p;
                            margin_min = fc.margin;
                            n_dec_min = fc.n_dec;
                        }
                    }
                }
```

In the `if p.dump_star { register_star_proof(..., StarProof { ... }) }` block add `gpos: gpos.clone(),` and add the field `pub gpos: Vec<Vec<Option<u64>>>,` to `StarProof` (update the `star_proofs_of_one_molecule_in_two_families_are_both_kept` test's `mk` to set `gpos: vec![]`).

- [ ] **Step 6: Run the whole library suite**

Run: `CARGO_TARGET_DIR=/mnt/linuxdisk/home/juanfraitu/rustle_target TMPDIR=/mnt/linuxdisk/tmp bash tools/rlock.sh heavy bash -c "cargo test --release --lib > /mnt/linuxdisk/tmp/popfix_t6.log 2>&1"; grep -E 'test result|FAILED' /mnt/linuxdisk/tmp/popfix_t6.log`
Expected: `test result: ok.` (1,038 passed), 0 failed.

- [ ] **Step 7: Commit**

```bash
git add src/family.rs
git commit -m "popfix: read-star column genome positions, population tier side table, --psv-fixed-only override"
```

---

### Task 7: binary flags, tier output, star-dump positions, guards

**Files:**
- Modify: `src/copy_assign_bin.rs` — the `Args` struct (next to `dump_star`, ≈ line 1039), argument validation (next to `if args.copies_fa.is_some()` checks, ≈ line 2147), the params writer (`row("psv_engine", ...)`), the `--dump-star` writer (≈ line 10470), and a new writer after it.

**Interfaces:**
- Consumes: Task 5 `parse_psv_population`, `register_psv_population`; Task 6 `take_population_tiers`, `PopulationTier`, `StarProof::gpos`.
- Produces: flags `--psv-population FILE`, `--psv-fixed-only`; output `<out>.population_tier.tsv` with header `read_name\tfamily_id\tstatus_all\tcopy_all\tstatus_fixed\tcopy_fixed\tn_dec_fixed\ttier`; `star_reads.tsv` gains a 10th column `positions` (per column `c0/c1/...` genome positions in candidate order, `.` = none, columns separated by `,`); `pub fn check_population_flags(population: bool, fixed_only: bool, two_form_genomic: bool) -> Result<(), String>` (lib-free helper in the binary, unit-tested there).

- [ ] **Step 1: Write the failing test** (binary unit test module at the end of `src/copy_assign_bin.rs`, `#[cfg(test)] mod popfix_flag_tests`)

```rust
#[cfg(test)]
mod popfix_flag_tests {
    use super::check_population_flags;

    #[test]
    fn fixed_only_requires_the_table_and_the_default_read_star_path() {
        assert!(check_population_flags(false, false, true).is_ok());
        assert!(check_population_flags(true, true, true).is_ok());
        assert!(check_population_flags(false, true, true).unwrap_err().contains("--psv-population"));
        assert!(check_population_flags(true, false, false).unwrap_err().contains("read-star"));
    }

    #[test]
    fn tier_names_follow_the_spec() {
        use rustle::copy_assign::copy_assign_pipeline::AssignStatus as S;
        assert_eq!(super::tier_name(S::Assigned, 3, S::Assigned, 3), "PORTABLE");
        assert_eq!(super::tier_name(S::Assigned, 3, S::Ambiguous, 3), "REFERENCE_ONLY");
        assert_eq!(super::tier_name(S::Assigned, 3, S::Assigned, 4), "REFERENCE_ONLY");
        assert_eq!(super::tier_name(S::Tied, 3, S::Assigned, 4), "UNASSIGNED");
    }
}
```

- [ ] **Step 2: Run to verify it fails**

Run: `CARGO_TARGET_DIR=/mnt/linuxdisk/home/juanfraitu/rustle_target TMPDIR=/mnt/linuxdisk/tmp bash tools/rlock.sh heavy bash -c "cargo test --release --bin copy_assign popfix_flag_tests > /mnt/linuxdisk/tmp/popfix_t7_red.log 2>&1"; grep -E '^error\[' /mnt/linuxdisk/tmp/popfix_t7_red.log | sort | uniq -c`
Expected: `cannot find function check_population_flags / tier_name`.

- [ ] **Step 3: Implement flags, validation, helpers**

In `Args` (next to `dump_star`):

```rust
    /// spec 2026-10-12-population-fixed-psvs: a `psv_population.tsv` (bench/popfix/table.py). Writes
    /// `<out>.population_tier.tsv`: per contested molecule the normal call, the call on population-FIXED PSVs only, and
    /// the tier PORTABLE / REFERENCE_ONLY / UNASSIGNED. Needs the default read-star path (genomic, two-form).
    #[arg(long)]
    psv_population: Option<String>,
    /// With `--psv-population`: decide on population-FIXED PSVs only (assignments.tsv = portable calls).
    #[arg(long, default_value_t = false)]
    psv_fixed_only: bool,
```

Free functions (top level of the binary):

```rust
fn check_population_flags(population: bool, fixed_only: bool, two_form_genomic: bool) -> Result<(), String> {
    if fixed_only && !population {
        return Err("--psv-fixed-only needs --psv-population FILE".into());
    }
    if population && !two_form_genomic {
        return Err("--psv-population needs the default read-star path (genomic, two-form); drop the read-star escape flags".into());
    }
    Ok(())
}

fn tier_name(
    s_all: rustle::copy_assign::copy_assign_pipeline::AssignStatus,
    c_all: usize,
    s_fixed: rustle::copy_assign::copy_assign_pipeline::AssignStatus,
    c_fixed: usize,
) -> &'static str {
    use rustle::copy_assign::copy_assign_pipeline::AssignStatus as S;
    match (s_all, s_fixed) {
        (S::Assigned, S::Assigned) if c_all == c_fixed => "PORTABLE",
        (S::Assigned, _) => "REFERENCE_ONLY",
        _ => "UNASSIGNED",
    }
}
```

Directly after `let params = AssignParams { ... };` (`src/copy_assign_bin.rs` ≈ line 6202; `region_families` from
`load_supplied_families` ≈ line 6120 is in scope):

```rust
    check_population_flags(args.psv_population.is_some(), args.psv_fixed_only, params.read_star_genomic && params.read_star_two_form)
        .map_err(|e| anyhow::anyhow!(e))?;
    if let Some(path) = args.psv_population.as_deref() {
        let text = std::fs::read_to_string(path).with_context(|| format!("reading --psv-population {path}"))?;
        let pop = rustle::copy_assign::copy_assign_pipeline::parse_psv_population(&text, args.psv_fixed_only)?;
        // Review Focus: chrom names must match the catalog
        let catalog_chroms: std::collections::HashSet<String> = region_families
            .iter()
            .flat_map(|m| m.values().flatten().flat_map(|f| f.copies.iter().map(|c| c.chrom.clone())))
            .collect();
        let on_catalog = pop.fixed.iter().filter(|(c, _)| catalog_chroms.contains(c.as_str())).count();
        if on_catalog == 0 {
            eprintln!("[copy_assign] ⚠ --psv-population: 0 of {} FIXED sites lie on a catalog copy's chromosome (chrom names?); every molecule will be UNASSIGNED", pop.fixed.len());
        }
        rustle::copy_assign::copy_assign_pipeline::register_psv_population(pop);
    }
```

Params rows (next to `psv_engine`):

```rust
        row("psv_population", args.psv_population.clone().unwrap_or_else(|| "NONE".into()))?;
        row("psv_fixed_only", format!("{}", args.psv_fixed_only))?;
```

Tier writer (after the `--dump-star` writer), reusing `tid_family`:

```rust
    if args.psv_population.is_some() {
        let tiers = rustle::copy_assign::copy_assign_pipeline::take_population_tiers();
        let mut th = std::fs::File::create(format!("{}.population_tier.tsv", args.out))?;
        writeln!(th, "read_name\tfamily_id\tstatus_all\tcopy_all\tstatus_fixed\tcopy_fixed\tn_dec_fixed\ttier")?;
        let mut seen: std::collections::HashSet<(String, String)> = std::collections::HashSet::new();
        for r in &assign_rows {
            let Some(t) = tiers.get(&r.read_name).and_then(|v| v.iter().find(|t| t.cand_tids.first().and_then(|x| tid_family.get(x)) == Some(&r.family_id))) else { continue };
            if !seen.insert((r.read_name.clone(), r.family_id.clone())) { continue; }
            let st = |s| format!("{s:?}").to_lowercase();
            writeln!(th, "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}", r.read_name, r.family_id, st(t.status_all), t.copy_all, st(t.status_fixed), t.copy_fixed, t.n_dec_fixed, tier_name(t.status_all, t.copy_all, t.status_fixed, t.copy_fixed))?;
        }
    }
```

and make `tid_family` filled whenever `args.dump_star || args.psv_population.is_some()` (change the `if args.dump_star {` guard where it is filled).

Star-dump positions column: in the `--dump-star` writer add `\tpositions` to the header and append per column the candidates' positions:

```rust
            let positions: Vec<String> = (0..pf.cols.len())
                .map(|j| pf.gpos.iter().map(|v| v.get(j).copied().flatten().map_or(".".to_string(), |g| g.to_string())).collect::<Vec<_>>().join("/"))
                .collect();
```

written as a 10th field `positions.join(",")`.

- [ ] **Step 4: Run tests**

Run: `CARGO_TARGET_DIR=/mnt/linuxdisk/home/juanfraitu/rustle_target TMPDIR=/mnt/linuxdisk/tmp bash tools/rlock.sh heavy bash -c "cargo test --release --bin copy_assign popfix_flag_tests > /mnt/linuxdisk/tmp/popfix_t7.log 2>&1 && cargo build --release >> /mnt/linuxdisk/tmp/popfix_t7.log 2>&1"; grep -E 'test result|^error' /mnt/linuxdisk/tmp/popfix_t7.log`
Expected: `test result: ok. 2 passed`, build OK.

- [ ] **Step 5: Byte-identity without the flags (spec §4 d)**

```bash
O=/mnt/linuxdisk/tmp/psv_sun; R=/mnt/linuxdisk/tmp/popfix/nofl; mkdir -p $R
bash tools/rlock.sh heavy /mnt/linuxdisk/home/juanfraitu/rustle_target/release/copy_assign --bam /mnt/linuxdisk/home/juanfraitu/fibroblasts/GCA_029281585.2_flnc_mm.bam --fasta /mnt/linuxdisk/home/juanfraitu/_from_wsl/winloci_scratch/GGO.fasta --regions $O/run/regions.txt --families $O/run/copies.tsv --copies-fa $O/run/copies.fa --dump-psv --threads 2 --out $R/s > $R/s.log 2>&1
for f in assignments quant families psv_cols psv_copies; do cmp -s $O/run_default/s.$f.tsv $R/s.$f.tsv && echo "$f identical" || echo "$f DIFFERS"; done
```
Expected: all five `identical`.

- [ ] **Step 6: Commit**

```bash
git add src/copy_assign_bin.rs src/family.rs
git commit -m "popfix: copy_assign --psv-population / --psv-fixed-only, population_tier.tsv, star-dump positions, flag guards"
```

---

### Task 8: pre-register and run the validation (spec §4)

**Files:**
- Create: `docs/PREREG_popfix_2026-10-12.md` (commit BEFORE any run below)
- Modify: `bench/yag_o2/split_half.py` (fixed-only mode)
- Create: `bench/popfix/sanity.py`
- Results appended to the prereg.

**Interfaces:**
- Consumes: `psv_population.tsv`; `copy_assign` with the new flags; `star_reads.tsv` `positions` column.
- Produces: `split_half.py` env `POPULATION=psv_population.tsv` -> columns used only where both candidates' positions are FIXED (same rule as `fixed_only_certificate`).

- [ ] **Step 1: Write the prereg** — bars copied from spec §4: (b) simulation `--psv-fixed-only`: 0 wrong among assigned, 0 of the undecidable assigned; (c) A119b testis, YAG families: split-half concordance on fixed PSVs >= 0.98 with >= 30 molecules whose halves both assign; PORTABLE count reported; (d) byte-identity (Task 7 Step 5) and full suite (3 known `pipeline_driver_dispatch` failures excepted); (a) descriptive. Commit it.

- [ ] **Step 2: Fixed-only split-half** — in `split_half.py`, when `POPULATION` is set, load FIXED `(chrom, pos)`; parse the `positions` field; in `certify` treat a column as present for pair (a, b) only if both candidates' positions at it are FIXED (chrom = the family's chrom from the catalog). Test: add to `bench/yag_o2/test_split_half.py`:

```python
import os, sys
sys.path.insert(0, os.path.dirname(__file__))
import split_half as S


def test_fixed_only_mask_drops_non_fixed_columns():
    cols = [(10, "A", "AC"), (20, "G", "GT")]
    pos = [["1000", "2000"], ["1010", "2010"]]          # per column: candidates' positions
    fixed = {("chrY", 1000), ("chrY", 2000)}             # column 1 not fixed
    st, cp = S.certify([0, 1], S.mask_fixed(cols, pos, "chrY", fixed))
    assert st == "ambiguous"                               # one fixed column (m = 1) never assigns
```

with `mask_fixed(cols, pos, chrom, fixed)` returning the columns whose every non-`.` candidate position is fixed.

- [ ] **Step 3: Runs** — (a) `python3 bench/popfix/sanity.py /mnt/linuxdisk/tmp/popfix/psv_population.tsv /mnt/linuxdisk/tmp/yag_o2/yag_human/yag_families.json` printing status x method counts per YAG family; (b) the 10-06 chrY simulation command (as in `docs/PREREG_yag_o2_2026-10-12.md` Arm A) plus `--psv-population /mnt/linuxdisk/tmp/popfix/psv_population.tsv --psv-fixed-only --dump-star`, scored with `bench/yag_o2/arm_a.py` (YAG families); (c) the A119b YAG run (`/mnt/linuxdisk/tmp/yag_o2/yag_human/` inputs) plus `--psv-population ... --dump-star`, then `POPULATION=... python3 bench/yag_o2/split_half.py ...` and PORTABLE counts from `population_tier.tsv`; (d) `cargo test --release --no-fail-fast`.

- [ ] **Step 4: Append results + verdicts to the prereg, commit, push.**

---

## Self-review notes

- Spec §1 -> Task 1; §2 k-mer -> Task 2, alignment -> Task 3, table/statuses -> Task 4; §3 positions -> Task 6, tier/strict -> Tasks 5-7, separate file -> Task 7, byte-identity -> Task 7 Step 5; §4 -> Task 8.
- Interpretation recorded: an alignment verdict needs ALL sites of a column to have unique flanks (else the k-mer verdict stands); a site in two columns takes the stricter status.
- Review Focus tests: chrom mismatch warning (Task 7 Step 3), 0-based positions (Task 6 test), '-' unit form (Task 6 test), flag guards (Task 7 test), no-position columns (Task 5 test).
