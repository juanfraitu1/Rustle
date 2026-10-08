# O3 Maternal- and Paternal-Reference Study Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Measure what happens to the reads of copies one parental haplotype lacks (unmapped / tied / absorbed), run the truth-free recovery chain (IsoCon and the in-house `o3_candidates`) from that haplotype's alignment alone, score both against the other haplotype, and show the two haplotypes side by side in one artifact: where the copy is missing, what the chain recovers; where it is present, it is simply in the reference.

**Architecture:** A package `bench/o3_maternal/` of pure, unit-tested functions (`common.py`) plus thin CLI steps that reuse the registered machinery unchanged (`bench/rna_allele/control_test.py net/outputs/contigs`, `merge_test.py`, `panel_to_copies.py`, `truth_lift.py`, the `o3_candidates` binary). One environment switch, `O3_REF` (`mat` = primary run, `pat` = the reverse), selects the reference haplotype; the other haplotype is the truth. Reads and their alignments to both haplotypes are shared; every other output lives under `W/<REF>/`. Heavy jobs run in the foreground under `tools/rlock.sh heavy`, one at a time. Work directory `/mnt/linuxdisk/tmp/o3_mat/` (`W`).

**Tech Stack:** Python 3 (`pysam` 0.23, `unittest`), minimap2 2.30, samtools, IsoCon 0.3.3 (conda env `isocon`), the Rust binary `o3_candidates`, a hand-written inline-SVG artifact page.

**Spec:** `docs/PREREG_o3_maternal_reference_2026-10-08.md` — the approved body plus **Amendment 1** (paternal-reference run and side-by-side; written 2026-10-08 at the user's request, pending the user's review together with this plan). Executors read both first; section numbers below (S3, S5, ...) refer to the prereg body.

## Global Constraints

- Substrate: KB3781 fibroblast Iso-Seq only. Reference of a run = `O3_REF` (`mat`: `GCA_028885495.2`; `pat`: `GCA_028885475.2`); truth = the other haplotype; indexes `winloci_data/mGorGor1.{mat,pat}.splice.mmi`, FASTAs `gorilla_haps/{mat,pat}.fa`. Testis excluded. No human numbers.
- Order (Amendment 1): the `mat` run is finished and its verdicts written to `W/mat/VERDICTS.txt` BEFORE the `pat` run starts. The `pat` run applies the same rules with no tuning.
- Registered constants, **no new ones**: allele cutoff `DELTA = 0.00958`; hit coverage `0.80`; hit identity `0.90`; tie rule `0.98` of the primary alignment score; expressed `>= 3` reads, LARGE `>= 20`; recovery hit `identity x coverage >= 0.999`; flag floor `>= 2` transcripts; family net `<= 1,000` reads, seed 1; per-locus read cap `2,000`, seed 1.
- Read alignment uses the fibroblast BAM's own command: `minimap2 -ax splice:hq -uf --eqx -Y -N 50 -p 0.1 --secondary=yes -K 100M -t 4` (minimap2 2.30 here vs 2.31 for the baseline BAM: disclose). Always record `-p` and `-N` next to any identity or copy count.
- chrY copy = SEX CONTROL (absent from `mat` by sex): excluded from every bar. p12 = reported on its own line, excluded from R1/R2.
- Crash rule (WSL2): foreground only, serial, one heavy job at a time (`tools/rlock.sh heavy`, default timeout 600 s, override with `RLOCK_TIMEOUT`); no `pkill -f`; no `nohup`; `ps` for orphans first. Light Python (< 2 GB, < 3 min) may use `tools/rlock.sh light`.
- Disk: `/mnt/linuxdisk` had ~25 GB free; never write a new minimap2 index; each renamed `W/<REF>.idx.fa` (3.6 GB) is deleted when its run's in-house batches are done.
- Build rule if the Rust binary must be rebuilt: `CARGO_TARGET_DIR=/mnt/linuxdisk/home/juanfraitu/rustle_target`, `--release`, cargo output to a FILE.
- Max 5 agents per workflow. Verifiers get frozen read-only copies, never the live work dir.
- Tests: `python3 -B -m unittest bench/o3_maternal/test_<name>.py` (repo convention, `bench/hierarchy`). Commit only `bench/o3_maternal/` and the `docs/` files named in a task; the worktree has unrelated uncommitted files: never `git add -A`. Do not push; do not share the artifact beyond its private link unless asked.
- **Run environment.** Every run step starts with this block (each Bash call is a fresh shell). `O3_REF` defaults to `mat`; the paternal pass (Task 10) sets `O3_REF=pat` before it.

```bash
cd /mnt/linuxdisk/home/juanfraitu/rustle_m2_soto
export O3_REF=${O3_REF:-mat}; REF=$O3_REF; OTHER=$([ "$REF" = mat ] && echo pat || echo mat)
export TMPDIR=/mnt/linuxdisk/tmp/o3_mat/tmp; W=/mnt/linuxdisk/tmp/o3_mat; WR=$W/$REF; mkdir -p $TMPDIR $WR
IDX=/mnt/linuxdisk/home/juanfraitu/winloci_data/mGorGor1; IDXREF=$IDX.$REF.splice.mmi; IDXOTH=$IDX.$OTHER.splice.mmi
commit() { git commit -F - <<MSG
$1

Co-Authored-By: Claude Sonnet 5.5 <noreply@anthropic.com>
Claude-Session: https://claude.ai/code/session_01W5su8r24NdRauw9JKTv34k
MSG
}
```

## Review Focus

Failure modes the spec implies but no happy-path test exercises, most likely first. Each has a pinning test in the owning task.

1. A read that is in both the 34-family set and an LRPAP1 locus (family collision): keeps its 34-family label and is still usable for the LRPAP1 loci (Task 2 `test_overlap_read_keeps_family_label`).
2. Haplotype contig names that are not `chrN_<hap>_hsaX` (unplaced scaffolds): pass through unchanged instead of `KeyError` (Task 1 `test_unknown_name_passes_through`).
3. A locus whose reads are all tied/unplaced on the truth haplotype has 0 labelled reads: the fate table and bars must not divide by zero (Task 5 `test_zero_reads_has_no_verdict`).
4. R34 has 0 unmapped reads by construction (every read mapped on `_pri` first): UNMAPPED must count only reads with no primary record, never be inferred from absence in a BAM (Task 1 `test_no_record_is_unmapped`, Task 5 selection note).
5. A flagged candidate set that is empty, or a chrY candidate: the scorer must print R1 FAIL / NOT TESTABLE and must never count the sex locus (Task 9 `test_empty_candidates_fail_r1`, `test_sex_locus_excluded`).
6. A chromosome present on one haplotype only (chrY in `pat`; chrX in `mat`) or a lift below 50%: the locus has NO counterpart, not a zero-length one (Task 4 `test_chry_is_pat_only_and_sex`, `test_poor_lift_means_no_counterpart`).
7. Pairing a read's two placements when one haplotype has no primary for it: the read is skipped, never paired with `None` (Task 11 `test_pairs_need_a_primary_on_both`).

---

## File Structure

```
bench/o3_maternal/
  common.py            constants (O3_REF switch), alias/accession, PAF + BAM readers, place(), classify_fate()
  testutil.py          write_bam() for synthetic BAMs in tests
  extract_reads.py     R_LRP (LRPAP1 net) and R_unm (unmapped primaries)          [shared by both runs]
  map_reads.sh         the fibroblast BAM's minimap2 command against a haplotype index [shared]
  truth.py             absent loci (S3) and read labels (S5), per reference
  fate.py              Q1 table: fates, verdicts, nearest reference paralog, per reference
  chain_inputs.py      panel, labels, scored.fa, renamed FASTA, per reference
  iso_batch.sh         resumable IsoCon per family (copy of refabsent/iso_batch.sh, H from env)
  adapt_isocon.py      IsoCon arm -> cands.tsv / cands.fa
  adapt_inhouse.py     in-house arm -> cands.tsv / cands.fa / cands_all.tsv
  score.py             R1-R4, p12 line, per reference
  side.py              paired reads + divergence pile, per reference          (Amendment 1)
  compare.py           the chain's view of the same copy in the two runs      (Amendment 1)
  artifact_data.py     data.json + index.html for both directions
  template.html        the artifact page
  test_<name>.py       one per module above (common, extract_reads, truth, fate, chain_inputs, score, side, compare)
docs/O3_MATERNAL_REFERENCE_2026-10-08.md     results write-up (Task 13)
docs/REGISTER_DRAFTS_o3_maternal.md          proposed register rows (Task 13)
```

---

### Task 1: Shared helpers, the reference switch and the fate rule

**Files:**
- Create: `bench/o3_maternal/common.py`, `bench/o3_maternal/testutil.py`, `bench/o3_maternal/test_common.py`

**Interfaces:**
- Produces (used by every later task):
  - constants `TRUTH, REF, OTHER, W, WR, HAP_FA, HAP_IDX, DELTA, COV_MIN, TIE` (`REF` from env `O3_REF`, default `mat`; `OTHER` the other haplotype; `WR = W/REF`)
  - `Rec = namedtuple("Rec", "primary mapq score qcov ref start end de")`
  - `alias() -> {(hap, num): accession}`; `accession(name, al) -> str`
  - `read_records(path, al=None) -> {read: [Rec]}` (primary first, supplementary dropped, unmapped -> `[]`)
  - `place(recs) -> Rec | None`; `classify_fate(recs, paralog) -> "UNMAPPED"|"PARTIAL"|"TIED"|"ABSORBED_NEAREST"|"ABSORBED_OTHER"`
  - `paf_hits(path, al, idmin=0.90, covmin=0.80) -> {query: [(acc, s, e, ident, cov)]}`; `best_hits(path, al) -> {query: (score, acc, s, e, ident)}`
  - `testutil.write_bam(path, refs, segs)`

- [ ] **Step 1: Write the failing tests**

`bench/o3_maternal/testutil.py`:

```python
"""Synthetic BAMs for the o3_maternal tests."""
import os

import pysam


def write_bam(path, refs, segs):
    """refs: {name: length}; segs: [dict(name, flag=0, ref, start, cigar='100M', mapq=60, AS=100, de=0.01, seq='A'*100)].
    Written coordinate-sorted and indexed."""
    tmp = path + ".unsorted"
    header = {"HD": {"VN": "1.6"}, "SQ": [{"SN": n, "LN": l} for n, l in refs.items()]}
    ids = {n: i for i, n in enumerate(refs)}
    with pysam.AlignmentFile(tmp, "wb", header=header) as f:
        for s in segs:
            a = pysam.AlignedSegment(f.header)
            a.query_name = s["name"]
            a.flag = s.get("flag", 0)
            unm = bool(a.flag & 4)
            a.reference_id = -1 if unm else ids[s["ref"]]
            a.reference_start = -1 if unm else s["start"]
            a.mapping_quality = s.get("mapq", 60)
            if not unm:
                a.cigarstring = s.get("cigar", "100M")
            a.query_sequence = s.get("seq", "A" * 100)
            if not unm:
                a.set_tag("AS", s.get("AS", 100))
                a.set_tag("de", s.get("de", 0.01))
            f.write(a)
    pysam.sort("-o", path, tmp)
    pysam.index(path)
    os.remove(tmp)
```

`bench/o3_maternal/test_common.py`:

```python
#!/usr/bin/env python3
"""Run: python3 -B -m unittest bench/o3_maternal/test_common.py"""
import os
import sys
import tempfile
import unittest

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import common as C  # noqa: E402
import testutil  # noqa: E402
from common import Rec  # noqa: E402


def rec(primary=True, mapq=60, score=1000, qcov=1.0, ref="chr1", start=100, end=200, de=0.01):
    return Rec(primary, mapq, score, qcov, ref, start, end, de)


class Fate(unittest.TestCase):
    def test_no_record_is_unmapped(self):
        self.assertEqual(C.classify_fate([], None), "UNMAPPED")

    def test_only_secondary_is_unmapped(self):
        self.assertEqual(C.classify_fate([rec(primary=False)], None), "UNMAPPED")

    def test_low_coverage_primary_is_partial(self):
        self.assertEqual(C.classify_fate([rec(qcov=0.79)], None), "PARTIAL")

    def test_coverage_boundary_is_mapped(self):
        self.assertEqual(C.classify_fate([rec(qcov=0.80)], None), "ABSORBED_OTHER")

    def test_mapq0_is_tied(self):
        self.assertEqual(C.classify_fate([rec(mapq=0)], None), "TIED")

    def test_close_secondary_is_tied(self):
        self.assertEqual(C.classify_fate([rec(score=1000), rec(primary=False, score=981)], None), "TIED")

    def test_far_secondary_is_not_tied(self):
        self.assertEqual(C.classify_fate([rec(score=1000), rec(primary=False, score=979)], None), "ABSORBED_OTHER")

    def test_absorbed_on_nearest_paralog(self):
        self.assertEqual(C.classify_fate([rec(ref="A", start=100, end=200)], ("A", 150, 400)), "ABSORBED_NEAREST")

    def test_other_chromosome_is_other(self):
        self.assertEqual(C.classify_fate([rec(ref="B", start=100, end=200)], ("A", 150, 400)), "ABSORBED_OTHER")


class Place(unittest.TestCase):
    def test_untied_best(self):
        p = C.place([rec(score=1000, start=1), rec(primary=False, score=900, start=2)])
        self.assertEqual(p.start, 1)

    def test_tied_is_none(self):
        self.assertIsNone(C.place([rec(score=1000), rec(primary=False, score=985)]))

    def test_empty_is_none(self):
        self.assertIsNone(C.place([]))


class Names(unittest.TestCase):
    def test_index_name_maps_to_accession(self):
        self.assertEqual(C.accession("chr3_mat_hsa4", {("mat", "3"): "CM1"}), "CM1")

    def test_unknown_name_passes_through(self):
        self.assertEqual(C.accession("scaffold_12", {}), "scaffold_12")
        self.assertEqual(C.accession("chr9_mat_hsa9", {}), "chr9_mat_hsa9")


class Bam(unittest.TestCase):
    def test_read_records(self):
        with tempfile.TemporaryDirectory() as d:
            p = os.path.join(d, "t.bam")
            testutil.write_bam(p, {"chr1_mat_hsa1": 5000}, [
                dict(name="r1", ref="chr1_mat_hsa1", start=100, AS=190, de=0.01),
                dict(name="r1", flag=256, ref="chr1_mat_hsa1", start=900, mapq=0, AS=185, de=0.02),
                dict(name="r1", flag=2048, ref="chr1_mat_hsa1", start=900, cigar="50S50M", mapq=0, AS=90),
                dict(name="r2", flag=4),
                dict(name="r3", ref="chr1_mat_hsa1", start=300, cigar="20S80M", AS=70),
            ])
            recs = C.read_records(p, {("mat", "1"): "CM1"})
            self.assertEqual(len(recs["r1"]), 2)
            self.assertTrue(recs["r1"][0].primary)
            self.assertEqual(recs["r1"][0].score, 190)
            self.assertEqual(recs["r1"][0].ref, "CM1")
            self.assertAlmostEqual(recs["r1"][0].qcov, 1.0)
            self.assertFalse(recs["r1"][1].primary)
            self.assertEqual(recs["r2"], [])
            self.assertAlmostEqual(recs["r3"][0].qcov, 0.8)


class Paf(unittest.TestCase):
    def test_hits_and_best(self):
        line = lambda q, t, ts, te, m, aln, qs, qe, ql: "\t".join(
            [q, str(ql), str(qs), str(qe), "+", t, "9999", str(ts), str(te), str(m), str(aln), "60"])
        with tempfile.TemporaryDirectory() as d:
            p = os.path.join(d, "t.paf")
            open(p, "w").write("\n".join([
                line("q1", "chr1_mat_hsa1", 10, 1010, 990, 1000, 0, 1000, 1000),      # ident .99, cov 1.0
                line("q1", "chr2_mat_hsa2", 50, 550, 480, 500, 0, 500, 1000),         # cov .5: not a hit
                line("q2", "chr1_mat_hsa1", 70, 1070, 800, 1000, 0, 1000, 1000),      # ident .80: not a hit
            ]) + "\n")
            al = {("mat", "1"): "CM1", ("mat", "2"): "CM2"}
            h = C.paf_hits(p, al)
            self.assertEqual(h["q1"], [("CM1", 10, 1010, 0.99, 1.0)])
            self.assertNotIn("q2", h)
            b = C.best_hits(p, al)
            self.assertEqual(b["q1"][1:4], ("CM1", 10, 1010))
            self.assertAlmostEqual(b["q1"][0], 0.99)


if __name__ == "__main__":
    unittest.main()
```

- [ ] **Step 2: Run the tests to verify they fail**

Run: `cd /mnt/linuxdisk/home/juanfraitu/rustle_m2_soto && python3 -B -m unittest bench/o3_maternal/test_common.py`
Expected: `ModuleNotFoundError: No module named 'common'`.

- [ ] **Step 3: Write `bench/o3_maternal/common.py`**

```python
#!/usr/bin/env python3
"""Shared helpers of the maternal-reference study (docs/PREREG_o3_maternal_reference_2026-10-08.md). Pure functions are unit-tested in test_common.py."""
import collections
import csv
import os
import re

import pysam

TRUTH = "/mnt/linuxdisk/tmp/rna_allele"
REF = os.environ.get("O3_REF", "mat")      # the reference haplotype of this run: mat (copies only the father has are missing) or pat (the reverse)
assert REF in ("mat", "pat"), REF
OTHER = "pat" if REF == "mat" else "mat"   # the truth haplotype: where the missing copies are present
W = "/mnt/linuxdisk/tmp/o3_mat"            # shared by both runs: reads/, map/, artifact/
WR = f"{W}/{REF}"                          # per-run outputs: truth/, fate/, isoc/, inhouse/, score/
HAP_FA = "/mnt/linuxdisk/home/juanfraitu/gorilla_haps/{}.fa"
HAP_IDX = "/mnt/linuxdisk/home/juanfraitu/winloci_data/mGorGor1.{}.splice.mmi"
DELTA = 0.00958        # merge_test.DELTA: the registered allele cutoff (Amendments 7-10)
COV_MIN = 0.80         # registered hit-coverage rule
TIE = 0.98             # registered tie rule (refabsent_truth.express)
IDX = re.compile(r"chr(\w+?)_(mat|pat)_hsa[^_]*")

Rec = collections.namedtuple("Rec", "primary mapq score qcov ref start end de")


def alias():
    """(hap, chromosome number) -> GenBank accession, from TRUTH/{mat,pat}.len.tsv"""
    out = {}
    for h in ("mat", "pat"):
        for acc, num, _ in csv.reader(open(f"{TRUTH}/{h}.len.tsv"), delimiter="\t"):
            out[(h, num)] = acc
    return out


def accession(name, al):
    """haplotype index name chrN_<hap>_hsaX -> accession; any other name (an unplaced scaffold, an unknown number) passes through"""
    m = IDX.fullmatch(name)
    if m and (m.group(2), m.group(1)) in al:
        return al[(m.group(2), m.group(1))]
    return name


def read_records(path, al=None):
    """{read: [Rec, ...]} from a BAM: supplementary records dropped, primary first; an unmapped read maps to []"""
    al = al or {}
    out = {}
    with pysam.AlignmentFile(path) as bam:
        for rd in bam.fetch(until_eof=True):
            if rd.is_supplementary:
                continue
            recs = out.setdefault(rd.query_name, [])
            if rd.is_unmapped:
                continue
            n = rd.infer_read_length() or 1
            recs.append(Rec(not rd.is_secondary, rd.mapping_quality, rd.get_tag("AS") if rd.has_tag("AS") else 0,
                            rd.query_alignment_length / n, accession(rd.reference_name, al), rd.reference_start,
                            rd.reference_end, rd.get_tag("de") if rd.has_tag("de") else None))
    for recs in out.values():
        recs.sort(key=lambda r: not r.primary)
    return out


def place(recs):
    """the read's single best placement over all its records, or None when there is none or the runner-up scores within TIE of it"""
    rs = sorted(recs, key=lambda r: -r.score)
    if not rs or (len(rs) > 1 and rs[1].score >= TIE * rs[0].score):
        return None
    return rs[0]


def classify_fate(recs, paralog):
    """Fate of one read on the reference (prereg S5). recs: Rec list (primary first); paralog: (acc, start, end) of the
    locus' nearest reference paralog, or None."""
    prim = [r for r in recs if r.primary]
    if not prim:
        return "UNMAPPED"
    p = prim[0]
    if p.qcov < COV_MIN:
        return "PARTIAL"
    if p.mapq == 0 or any((not r.primary) and r.score >= TIE * p.score for r in recs):
        return "TIED"
    if paralog and p.ref == paralog[0] and p.start < paralog[2] and paralog[1] < p.end:
        return "ABSORBED_NEAREST"
    return "ABSORBED_OTHER"


def paf_hits(path, al, idmin=0.90, covmin=0.80):
    """{query: [(acc, start, end, identity, query coverage)]} of the PAF records at identity >= idmin and coverage >= covmin"""
    out = collections.defaultdict(list)
    for ln in open(path):
        f = ln.rstrip("\n").split("\t")
        ident = int(f[9]) / max(1, int(f[10]))
        cov = (int(f[3]) - int(f[2])) / max(1, int(f[1]))
        if ident >= idmin and cov >= covmin:
            out[f[0]].append((accession(f[5], al), int(f[7]), int(f[8]), ident, cov))
    return dict(out)


def best_hits(path, al):
    """{query: (score, acc, start, end, identity)}: best PAF record per query by identity x query coverage (control_test.best_hits)"""
    b = {}
    for ln in open(path):
        f = ln.rstrip("\n").split("\t")
        ident = int(f[9]) / max(1, int(f[10]))
        s = ident * (int(f[3]) - int(f[2])) / max(1, int(f[1]))
        if f[0] not in b or s > b[f[0]][0]:
            b[f[0]] = (s, accession(f[5], al), int(f[7]), int(f[8]), ident)
    return b
```

- [ ] **Step 4: Run the tests to verify they pass**

Run: `cd /mnt/linuxdisk/home/juanfraitu/rustle_m2_soto && python3 -B -m unittest bench/o3_maternal/test_common.py`
Expected: `Ran 16 tests ... OK`.

- [ ] **Step 5: Commit**

```bash
# run environment block first
git add bench/o3_maternal/common.py bench/o3_maternal/testutil.py bench/o3_maternal/test_common.py docs/PREREG_o3_maternal_reference_2026-10-08.md docs/superpowers/plans/2026-10-08-o3-maternal-reference.md
commit "o3_maternal: prereg (+ Amendment 1), plan, shared helpers, reference switch and the fate rule"
```

---

### Task 2: Extract the reads (R_LRP, R_unm) — shared by both runs

**Files:**
- Create: `bench/o3_maternal/extract_reads.py`, `bench/o3_maternal/test_extract_reads.py`
- Outputs (not committed): `W/reads/R_LRP.fa`, `W/reads/R_LRP.names.tsv`, `W/reads/R_unm.fa`

**Interfaces:**
- Consumes: `common.W`, `common.TRUTH`, `testutil.write_bam`.
- Produces: `loci() -> [(cid, name, chrom, lo0, hi)]` (the 11 LRPAP1 loci, pri coordinates); `net_names(bam_path, loci_, cap=2000, seed=1) -> {read: [cid]}`; `merge_labels(r34, lrp) -> (rows, overlap)`; constants `BAM, LRP, R34`; file `W/reads/R_LRP.names.tsv` with columns `read, cids, in_r34`.

- [ ] **Step 1: Write the failing tests**

`bench/o3_maternal/test_extract_reads.py`:

```python
#!/usr/bin/env python3
"""Run: python3 -B -m unittest bench/o3_maternal/test_extract_reads.py"""
import os
import sys
import tempfile
import unittest

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import extract_reads as E  # noqa: E402
import testutil  # noqa: E402


class Net(unittest.TestCase):
    def bam(self, d):
        p = os.path.join(d, "t.bam")
        testutil.write_bam(p, {"chrA": 100000}, [
            dict(name="a", ref="chrA", start=1000),
            dict(name="b", ref="chrA", start=1500),
            dict(name="b", flag=256, ref="chrA", start=40000, mapq=0),       # secondary on locus 2 only
            dict(name="c", ref="chrA", start=90000),                          # outside every locus
            dict(name="s", flag=2048, ref="chrA", start=1200, cigar="50S50M"),  # supplementary: ignored
        ])
        return p

    def test_reads_on_a_locus_primary_or_secondary(self):
        with tempfile.TemporaryDirectory() as d:
            n = E.net_names(self.bam(d), [("c1", "g1", "chrA", 500, 3000), ("c2", "g2", "chrA", 39000, 42000)])
            self.assertEqual(n, {"a": ["c1"], "b": ["c1", "c2"]})

    def test_cap_is_per_locus_and_seeded(self):
        with tempfile.TemporaryDirectory() as d:
            n1 = E.net_names(self.bam(d), [("c1", "g1", "chrA", 500, 3000)], cap=1)
            n2 = E.net_names(self.bam(d), [("c1", "g1", "chrA", 500, 3000)], cap=1)
            self.assertEqual(len(n1), 1)
            self.assertEqual(n1, n2)


class Merge(unittest.TestCase):
    def test_overlap_read_keeps_family_label(self):
        r34 = {"x": "GWFAM9"}
        rows, overlap = E.merge_labels(r34, {"x": ["p12"], "y": ["p12", "p14"]})
        self.assertEqual(overlap, ["x"])
        self.assertEqual(rows, [("y", "LRPAP1", "p12,p14")])


if __name__ == "__main__":
    unittest.main()
```

- [ ] **Step 2: Run to verify failure**

Run: `cd /mnt/linuxdisk/home/juanfraitu/rustle_m2_soto && python3 -B -m unittest bench/o3_maternal/test_extract_reads.py`
Expected: `ModuleNotFoundError: No module named 'extract_reads'`.

- [ ] **Step 3: Write `bench/o3_maternal/extract_reads.py`**

```python
#!/usr/bin/env python3
"""R_LRP and R_unm of docs/PREREG_o3_maternal_reference_2026-10-08.md section 4.

    extract_reads.py lrp    # reads with a primary/secondary record on one of the 11 LRPAP1 loci -> W/reads/R_LRP.{fa,names.tsv}
    extract_reads.py unm    # unmapped primaries of the fibroblast BAM -> W/reads/R_unm.fa
"""
import collections
import csv
import os
import random
import subprocess
import sys

import pysam

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import common as C  # noqa: E402

BAM = "/mnt/linuxdisk/home/juanfraitu/fibroblasts/GCA_029281585.2_flnc_mm.bam"
LRP = "/mnt/linuxdisk/tmp/lrpap1"
R34 = f"{C.TRUTH}/refabsent/labels.tsv"
CAP = 2000


def loci():
    """the 11 LRPAP1 loci (8 full-length + 3 fragments): (cid, name, chrom, lo0, hi), pri coordinates"""
    out = []
    for f in ("lrpap1.copies.tsv", "partial.copies.tsv"):
        for r in csv.DictReader(open(f"{LRP}/{f}"), delimiter="\t"):
            out.append((r["cid"], r["name"], r["chrom"], int(r["terr_lo0"]), int(r["terr_hi"])))
    return out


def net_names(bam_path, loci_, cap=CAP, seed=1):
    """{read: [cid, ...]}: reads with a primary or secondary record on a locus, at most `cap` per locus (seeded shuffle of the sorted names)"""
    rng = random.Random(seed)
    out = collections.defaultdict(list)
    with pysam.AlignmentFile(bam_path) as bam:
        for cid, _name, chrom, lo, hi in loci_:
            names = sorted({rd.query_name for rd in bam.fetch(chrom, lo, hi) if not (rd.is_unmapped or rd.is_supplementary)})
            rng.shuffle(names)
            for n in names[:cap]:
                out[n].append(cid)
    return dict(out)


def merge_labels(r34, lrp):
    """r34: {read: family} of the 34-family set; lrp: {read: [cid]}. -> (rows, overlap): rows = [(read, 'LRPAP1', 'cid,cid')] for reads
    NOT in the 34-family set; overlap = reads in both (they keep their 34-family label and are still LRPAP1-locus reads via names.tsv)"""
    rows, overlap = [], []
    for n, cids in sorted(lrp.items()):
        if n in r34:
            overlap.append(n)
        else:
            rows.append((n, "LRPAP1", ",".join(sorted(cids))))
    return rows, overlap


def sequences(bam_path, names, out_fa):
    """primary-record sequences in the original read orientation (samtools fasta restores the strand); ONE scan of the BAM"""
    nf = out_fa + ".names"
    open(nf, "w").write("\n".join(sorted(names)) + "\n")
    subprocess.run(f"samtools view -b -F 2308 -N {nf} -@ 4 {bam_path} | samtools fasta -@ 2 - > {out_fa}", shell=True, check=True)


def lrp():
    os.makedirs(f"{C.W}/reads", exist_ok=True)
    r34 = {r["read"]: r["family"] for r in csv.DictReader(open(R34), delimiter="\t")}
    names = net_names(BAM, loci())
    rows, overlap = merge_labels(r34, names)
    with open(f"{C.W}/reads/R_LRP.names.tsv", "w") as o:
        o.write("read\tcids\tin_r34\n")
        for n, cids in sorted(names.items()):
            o.write(f"{n}\t{','.join(sorted(cids))}\t{int(n in r34)}\n")
    sequences(BAM, [r[0] for r in rows], f"{C.W}/reads/R_LRP.fa")
    got = sum(1 for ln in open(f"{C.W}/reads/R_LRP.fa") if ln[0] == ">")
    print(f"LRPAP1 net reads {len(names)}; new (not in the 34-family set) {len(rows)}; sequences written {got}; in both {len(overlap)}")


def unm():
    os.makedirs(f"{C.W}/reads", exist_ok=True)
    out = f"{C.W}/reads/R_unm.fa"
    subprocess.run(f"samtools view -b -f 4 {BAM} '*' | samtools fasta - > {out}", shell=True, check=True)
    print("unmapped primaries", sum(1 for ln in open(out) if ln[0] == ">"))


if __name__ == "__main__":
    {"lrp": lrp, "unm": unm}[sys.argv[1]]()
```

- [ ] **Step 4: Run the tests to verify they pass**

Run: `cd /mnt/linuxdisk/home/juanfraitu/rustle_m2_soto && python3 -B -m unittest bench/o3_maternal/test_extract_reads.py`
Expected: `Ran 3 tests ... OK`.

- [ ] **Step 5: Run the real extraction (one BAM scan, ~7 min)**

```bash
# run environment block first
ps -eo pid,etimes,rss,cmd --sort=-rss | head -5          # no orphans
RLOCK_TIMEOUT=900 tools/rlock.sh heavy python3 bench/o3_maternal/extract_reads.py lrp
tools/rlock.sh light python3 bench/o3_maternal/extract_reads.py unm
```
Expected: `unmapped primaries 959`; the LRPAP1 line reports roughly 2,500-3,500 net reads and `sequences written` equal to `new`. If `sequences written` < `new`, stop and report.

- [ ] **Step 6: Commit**

```bash
# run environment block first
git add bench/o3_maternal/extract_reads.py bench/o3_maternal/test_extract_reads.py
commit "o3_maternal: extract the LRPAP1 net reads and the unmapped primaries"
```

---

### Task 3: Map the new reads to `mat` and `pat` — shared by both runs

**Files:**
- Create: `bench/o3_maternal/map_reads.sh`
- Outputs (not committed): `W/map/R_LRP.{mat,pat}.bam`, `W/map/R_unm.{mat,pat}.bam`, `W/map/reads.{mat,pat}.all.bam` (R34 + R_LRP merged)

**Interfaces:**
- Consumes: `W/reads/R_LRP.fa`, `W/reads/R_unm.fa`; existing `TRUTH/refabsent/reads.{mat,pat}.bam` (32,219 net reads of the 34 families, same command).
- Produces: `W/map/reads.{mat,pat}.all.bam` (indexed; reference names `chrN_<hap>_hsaX`); `W/map/R_unm.{mat,pat}.bam`.

- [ ] **Step 1: Write `bench/o3_maternal/map_reads.sh`**

```bash
#!/bin/bash
# map_reads.sh <mat|pat> <reads.fa> <out.bam>: the fibroblast BAM's own minimap2 command (@PG of GCA_029281585.2_flnc_mm.bam)
# against a haplotype splice index. Run under tools/rlock.sh heavy (loads a 13 GB index, ~15 GB RSS).
set -euo pipefail
hap=${1:?mat|pat}; fa=${2:?reads.fa}; out=${3:?out.bam}
IDX=/mnt/linuxdisk/home/juanfraitu/winloci_data/mGorGor1.$hap.splice.mmi
[ -s "$IDX" ] || { echo "missing $IDX" >&2; exit 2; }
minimap2 -ax splice:hq -uf --eqx -Y -N 50 -p 0.1 --secondary=yes -K 100M -t 4 "$IDX" "$fa" 2> "$out.log" \
  | samtools sort -@ 1 -m 500M -o "$out" -
samtools index "$out"
echo "$out: $(samtools view -c "$out") records, $(samtools view -c -F 2308 "$out") primaries, $(samtools view -c -f 4 "$out") unmapped"
```

- [ ] **Step 2: Verify the minimap2 and index preconditions**

```bash
# run environment block first
minimap2 --version                                   # expect 2.30-r1287
ls -l $IDX.mat.splice.mmi $IDX.pat.splice.mmi
df -h /mnt/linuxdisk | tail -1                       # expect >= 15 GB free
```

- [ ] **Step 3: Map R_LRP and R_unm to both haplotypes, serial, foreground**

```bash
# run environment block first
mkdir -p $W/map
for hap in mat pat; do tools/rlock.sh heavy bash bench/o3_maternal/map_reads.sh $hap $W/reads/R_LRP.fa $W/map/R_LRP.$hap.bam; done
for hap in mat pat; do tools/rlock.sh heavy bash bench/o3_maternal/map_reads.sh $hap $W/reads/R_unm.fa $W/map/R_unm.$hap.bam; done
```
Expected: four lines `... primaries ...`; R_LRP primaries equal the number of reads in `R_LRP.fa`.

- [ ] **Step 4: Merge R34 and R_LRP per haplotype and verify coverage of every labelled read**

```bash
# run environment block first
T=/mnt/linuxdisk/tmp/rna_allele/refabsent
for hap in mat pat; do
  samtools merge -f $W/map/reads.$hap.all.bam $T/reads.$hap.bam $W/map/R_LRP.$hap.bam && samtools index $W/map/reads.$hap.all.bam
  echo "$hap primaries: $(samtools view -c -F 2308 $W/map/reads.$hap.all.bam)  expected $(( $(grep -c '>' $W/reads/R_LRP.fa) + 32219 ))"
done
```
Expected: both counts equal the `expected` value. Any mismatch = a read name present in both inputs or a dropped read: stop and report.

- [ ] **Step 5: Commit**

```bash
# run environment block first
git add bench/o3_maternal/map_reads.sh
commit "o3_maternal: map_reads.sh (the baseline BAM's minimap2 command against mat/pat)"
```

---

### Task 4: Absent-locus truth and read labels (per reference)

**Files:**
- Create: `bench/o3_maternal/truth.py`, `bench/o3_maternal/test_truth.py`
- Outputs: `WR/truth/loci.tsv`, `WR/truth/lrpap1_loci.tsv`, `WR/truth/labels.tsv`

**Interfaces:**
- Consumes: `common.*`, `extract_reads.loci/LRP/R34`, `rna_allele/truth_lift.py`, `TRUTH/refabsent/bonly.tsv`, `TRUTH/chrmap.tsv`, `TRUTH/out/*.paf`, `LRP/{copies8,partial3}.{mat,pat}.paf`, `W/map/reads.<OTHER>.all.bam`.
- Produces:
  - `catalog_loci(path, other=None) -> [dict(locus, kind='catalog', family, chrom, start, end, name)]` (rows of `bonly.tsv` with `hap == other`)
  - `hap_intervals(loci_, chrmap, lift, hits, qof, al) -> {cid: dict(name, mat, pat, sex)}` (`mat`/`pat` = `(acc, s, e)` or `None`)
  - `label_reads(placements, loci_, fam_of, lrp_of) -> {read: locus | 'shared' | 'ambiguous'}`
  - files `WR/truth/loci.tsv` (columns `locus kind family chrom start end name`; `chrom` = accession on the truth haplotype), `lrpap1_loci.tsv` (columns `cid name sex pat_acc pat_s pat_e mat_acc mat_s mat_e`), `labels.tsv` (columns `read label`)

- [ ] **Step 1: Write the failing tests**

`bench/o3_maternal/test_truth.py`:

```python
#!/usr/bin/env python3
"""Run: python3 -B -m unittest bench/o3_maternal/test_truth.py"""
import os
import sys
import tempfile
import unittest

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import truth as T  # noqa: E402
from common import Rec  # noqa: E402


class Catalog(unittest.TestCase):
    def test_only_pat_rows(self):
        with tempfile.TemporaryDirectory() as d:
            p = os.path.join(d, "bonly.tsv")
            open(p, "w").write("family\tlocus\thap\tchrom\tstart\tend\n"
                               "F1\tF1_B0\tpat\tCM1\t10\t50\nF1\tF1_B1\tmat\tCM2\t5\t9\n")
            rows = T.catalog_loci(p, other="pat")
            self.assertEqual([r["locus"] for r in rows], ["F1_B0"])
            self.assertEqual([r["locus"] for r in T.catalog_loci(p, other="mat")], ["F1_B1"])
            self.assertEqual((rows[0]["chrom"], rows[0]["start"], rows[0]["end"], rows[0]["kind"]), ("CM1", 10, 50, "catalog"))


class Lrpap1(unittest.TestCase):
    loci = [("c0", "L0", "NC1", 100, 200), ("c1", "L1", "NC1", 300, 400), ("c2", "L2", "NC2", 100, 200),
            ("c3", "L3", "NCY", 100, 200), ("c4", "L4", "NC3", 100, 200)]
    chrmap = {"NC1": dict(same_hap="pat", same_name="CMP1"), "NC2": dict(same_hap="mat", same_name="CMM2"),
              "NC3": dict(same_hap="pat", same_name="CMP3")}
    lift = {"c0": dict(lift_frac="1.0", B_chrom="CMM1", B_start="1000", B_end="1100"),
            "c1": dict(lift_frac="1.0", B_chrom="CMM1", B_start="2000", B_end="2100"),
            "c2": dict(lift_frac="1.0", B_chrom="CMP2", B_start="50", B_end="150"),
            "c4": dict(lift_frac="0.2", B_chrom="CMM3", B_start="1", B_end="2")}
    hits = {"mat": {"q0": [("CMM1", 1010, 1090, 0.99, 1.0)], "q1": [("CMM9", 5, 50, 0.99, 1.0)], "q4": [("CMM3", 1, 2, 0.99, 1.0)]},
            "pat": {"q2": [("CMP2", 60, 140, 0.99, 1.0)]}}
    qof = {"c0": "q0", "c1": "q1", "c2": "q2", "c3": "q3", "c4": "q4"}

    def iv(self):
        return T.hap_intervals(self.loci, self.chrmap, self.lift, self.hits, self.qof, {("pat", "Y"): "CMY"})

    def test_pat_chromosome_with_a_mat_ortholog(self):
        d = self.iv()["c0"]
        self.assertEqual((d["pat"], d["mat"]), (("CMP1", 100, 200), ("CMM1", 1000, 1100)))

    def test_pat_chromosome_without_a_mat_ortholog_is_mat_absent(self):
        d = self.iv()["c1"]
        self.assertEqual((d["pat"], d["mat"]), (("CMP1", 300, 400), None))

    def test_mat_chromosome_with_a_pat_ortholog(self):
        d = self.iv()["c2"]
        self.assertEqual((d["mat"], d["pat"]), (("CMM2", 100, 200), ("CMP2", 50, 150)))

    def test_chry_is_pat_only_and_sex(self):
        d = self.iv()["c3"]
        self.assertEqual((d["pat"], d["mat"], d["sex"]), (("CMY", 100, 200), None, True))

    def test_poor_lift_means_no_counterpart(self):
        self.assertIsNone(self.iv()["c4"]["mat"])

    def test_query_of_picks_the_overlapping_query(self):
        qs = ["NC1:101-200", "NC1:5000-6000", "NC2:101-200"]
        out = T.query_of([("c0", "L0", "NC1", 100, 200), ("c2", "L2", "NC2", 100, 200)], qs)
        self.assertEqual(out, {"c0": "NC1:101-200", "c2": "NC2:101-200"})


class Labels(unittest.TestCase):
    loci = [dict(locus="F1_B0", kind="catalog", family="F1", chrom="CM1", start=100, end=200),
            dict(locus="LRPAP1_c9", kind="lrpap1", family="LRPAP1", chrom="CM2", start=100, end=200)]

    def rec(self, ref, s, e):
        return Rec(True, 60, 100, 1.0, ref, s, e, 0.0)

    def test_catalog_needs_family_match(self):
        lab = T.label_reads({"a": self.rec("CM1", 120, 180), "b": self.rec("CM1", 120, 180)}, self.loci, {"a": "F1", "b": "F2"}, {})
        self.assertEqual(lab, {"a": "F1_B0", "b": "shared"})

    def test_lrpap1_needs_net_membership(self):
        lab = T.label_reads({"a": self.rec("CM2", 120, 180), "b": self.rec("CM2", 120, 180)}, self.loci, {}, {"a": ["c9"]})
        self.assertEqual(lab, {"a": "LRPAP1_c9", "b": "shared"})

    def test_tied_is_ambiguous_and_elsewhere_is_shared(self):
        lab = T.label_reads({"a": None, "b": self.rec("CM7", 1, 9)}, self.loci, {}, {})
        self.assertEqual(lab, {"a": "ambiguous", "b": "shared"})


if __name__ == "__main__":
    unittest.main()
```

- [ ] **Step 2: Run to verify failure**

Run: `cd /mnt/linuxdisk/home/juanfraitu/rustle_m2_soto && python3 -B -m unittest bench/o3_maternal/test_truth.py`
Expected: `ModuleNotFoundError: No module named 'truth'`.

- [ ] **Step 3: Write `bench/o3_maternal/truth.py`**

```python
#!/usr/bin/env python3
"""Truth of docs/PREREG_o3_maternal_reference_2026-10-08.md section 3 (+ Amendment 1) and the read labels of section 5, for the run's reference
haplotype C.REF (env O3_REF, default mat): a locus is ABSENT iff it exists on C.OTHER and has no counterpart on C.REF.

    O3_REF=mat truth.py loci     # WR/truth/loci.tsv and lrpap1_loci.tsv
    O3_REF=mat truth.py labels   # WR/truth/labels.tsv (needs W/map/reads.<OTHER>.all.bam and WR/truth/loci.tsv)
"""
import collections
import csv
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.join(HERE, "..", "rna_allele"))
import common as C  # noqa: E402
import extract_reads as E  # noqa: E402


def catalog_loci(path=f"{C.TRUTH}/refabsent/bonly.tsv", other=None):
    """haplotype-only loci of the 378-family catalog lying on the truth haplotype `other` (default C.OTHER): the rows of bonly.tsv with
    hap == other (chromosomes `_pri` took from the REFERENCE haplotype); 0-based half-open on `other`"""
    other = other or C.OTHER
    out = []
    for r in csv.DictReader(open(path), delimiter="\t"):
        if r["hap"] == other:
            out.append(dict(locus=r["locus"], kind="catalog", family=r["family"], chrom=r["chrom"],
                            start=int(r["start"]), end=int(r["end"]), name=r["locus"]))
    return out


def fasta_names(path):
    return [ln[1:].strip() for ln in open(path) if ln[0] == ">"]


def query_of(loci_, queries):
    """cid -> the PAF query ('NC:lo+1-hi') on the same chromosome overlapping the locus the most"""
    out = {}
    for cid, _name, chrom, lo, hi in loci_:
        best = None
        for q in queries:
            c, rng = q.rsplit(":", 1)
            s, e = (int(x) for x in rng.split("-"))
            ov = min(hi, e) - max(lo, s - 1)
            if c == chrom and ov > 0 and (best is None or ov > best[0]):
                best = (ov, q)
        if best:
            out[cid] = best[1]
    return out


def hap_intervals(loci_, chrmap, lift, hits, qof, al):
    """cid -> dict(name, mat=(acc, s, e)|None, pat=(acc, s, e)|None, sex=bool).
    On a chromosome `_pri` took from haplotype H, `_pri` coordinates ARE H coordinates (chrmap seq_check identical) and the other haplotype's
    interval is the lift (B) iff lift_frac >= .5 and a hit of the locus body (identity >= .90, coverage >= .80) overlaps it.
    hits = {'mat': {query: [(acc, s, e, ident, cov)]}, 'pat': {...}}. A chromosome missing from chrmap is chrY: pat only (pri chrY = pat chrY)."""
    out = {}
    for cid, name, chrom, lo, hi in loci_:
        d = dict(name=name, mat=None, pat=None, sex=False)
        row, L = chrmap.get(chrom), lift.get(cid)
        if row is None:
            d["pat"], d["sex"] = (al[("pat", "Y")], lo, hi), True
        else:
            h = row["same_hap"]
            o = "mat" if h == "pat" else "pat"
            d[h] = (row["same_name"], lo, hi)
            if L is not None and float(L["lift_frac"]) >= 0.5 and any(
                    x[0] == L["B_chrom"] and x[1] < int(L["B_end"]) and int(L["B_start"]) < x[2] for x in hits[o].get(qof.get(cid), [])):
                d[o] = (L["B_chrom"], int(L["B_start"]), int(L["B_end"]))
        out[cid] = d
    return out


def label_reads(placements, loci_, fam_of, lrp_of):
    """placements: {read: Rec | None} the read's untied best placement on the TRUTH haplotype (accessions); loci_: loci dicts (start/end ints);
    fam_of: {read: family} of the 34-family set; lrp_of: {read: [cid]} of the LRPAP1 net.
    -> {read: locus id | 'shared' (untied placement elsewhere) | 'ambiguous' (tied or unplaced)}"""
    out = {}
    for n, p in placements.items():
        if p is None:
            out[n] = "ambiguous"
            continue
        hit = None
        for L in loci_:
            if p.ref != L["chrom"] or not (p.start < L["end"] and L["start"] < p.end):
                continue
            if L["kind"] == "catalog" and fam_of.get(n) != L["family"]:
                continue
            if L["kind"] in ("lrpap1", "sex") and n not in lrp_of:
                continue
            hit = L["locus"]
            break
        out[n] = hit or "shared"
    return out


def cmd_loci():
    import truth_lift
    os.makedirs(f"{C.WR}/truth", exist_ok=True)
    al = C.alias()
    L11 = E.loci()
    chrmap = {r["pri"]: r for r in csv.DictReader(open(f"{C.TRUTH}/chrmap.tsv"), delimiter="\t")}
    genes = f"{C.WR}/truth/lrpap1.genes.tsv"
    with open(genes, "w") as o:
        o.write("gene_id\tchrom\tstrand\texons\n")
        for cid, _n, chrom, lo, hi in L11:
            if chrom in chrmap:
                o.write(f"{cid}\t{chrom}\t+\t{lo}-{hi}\n")
    truth_lift.main(["--chrmap", f"{C.TRUTH}/chrmap.tsv", "--paf-dir", f"{C.TRUTH}/out", "--genes", genes,
                     "--out", f"{C.WR}/truth/lrpap1.lift.tsv"])
    lift = {r["gene_id"]: r for r in csv.DictReader(open(f"{C.WR}/truth/lrpap1.lift.tsv"), delimiter="\t")}
    queries = fasta_names(f"{E.LRP}/copies8.fa") + fasta_names(f"{E.LRP}/partial3.fa")
    hits = {"mat": {}, "pat": {}}
    for h in hits:
        for f in (f"copies8.{h}.paf", f"partial3.{h}.paf"):
            for q, v in C.paf_hits(f"{E.LRP}/{f}", al).items():      # asm20 -N 50 -p 0.5 (copies8) as produced on 10-04
                hits[h].setdefault(q, []).extend(v)
    iv = hap_intervals(L11, chrmap, lift, hits, query_of(L11, queries), al)
    with open(f"{C.WR}/truth/lrpap1_loci.tsv", "w") as o:
        o.write("cid\tname\tsex\tpat_acc\tpat_s\tpat_e\tmat_acc\tmat_s\tmat_e\n")
        for cid, d in iv.items():
            p, m = d["pat"] or ("", "", ""), d["mat"] or ("", "", "")
            o.write("\t".join(str(x) for x in [cid, d["name"], int(d["sex"]), *p, *m]) + "\n")
    rows = catalog_loci()
    for cid, d in iv.items():
        if d[C.REF] is None and d[C.OTHER] is not None:
            rows.append(dict(locus=f"LRPAP1_{cid}", kind="sex" if d["sex"] else "lrpap1", family="LRPAP1",
                             chrom=d[C.OTHER][0], start=d[C.OTHER][1], end=d[C.OTHER][2], name=d["name"]))
    with open(f"{C.WR}/truth/loci.tsv", "w") as o:
        o.write("locus\tkind\tfamily\tchrom\tstart\tend\tname\n")
        for r in rows:
            o.write("\t".join(str(r[k]) for k in ("locus", "kind", "family", "chrom", "start", "end", "name")) + "\n")
    print(f"reference {C.REF}, truth {C.OTHER}: LRPAP1 loci present on mat {sum(1 for d in iv.values() if d['mat'])}, on pat "
          f"{sum(1 for d in iv.values() if d['pat'])} of {len(iv)}; absent from {C.REF}: {[c for c, d in iv.items() if d[C.REF] is None]}")
    print(f"{C.REF}-absent loci: {len(rows)} ({sum(1 for r in rows if r['kind'] == 'catalog')} catalog, "
          f"{sum(1 for r in rows if r['kind'] == 'lrpap1')} LRPAP1, {sum(1 for r in rows if r['kind'] == 'sex')} sex control)")


def cmd_labels():
    al = C.alias()
    loci_ = list(csv.DictReader(open(f"{C.WR}/truth/loci.tsv"), delimiter="\t"))
    for L in loci_:
        L["start"], L["end"] = int(L["start"]), int(L["end"])
    fam_of = {r["read"]: r["family"] for r in csv.DictReader(open(E.R34), delimiter="\t")}
    lrp_of = {r["read"]: r["cids"].split(",") for r in csv.DictReader(open(f"{C.W}/reads/R_LRP.names.tsv"), delimiter="\t")}
    recs = C.read_records(f"{C.W}/map/reads.{C.OTHER}.all.bam", al)
    labels = label_reads({n: C.place(rs) for n, rs in recs.items()}, loci_, fam_of, lrp_of)
    with open(f"{C.WR}/truth/labels.tsv", "w") as o:
        o.write("read\tlabel\n")
        for n, g in sorted(labels.items()):
            o.write(f"{n}\t{g}\n")
    cnt = collections.Counter(labels.values())
    old = {r["locus"]: int(r["n_reads"]) for r in csv.DictReader(open(f"{C.TRUTH}/refabsent/bonly_expressed.tsv"), delimiter="\t")}
    print(f"reads {len(labels)}: shared {cnt['shared']}, ambiguous {cnt['ambiguous']}")
    print(f"locus\tkind\treads(new, best untied placement on {C.OTHER})\treads(Amendment 10 express, both haplotypes)")
    for L in loci_:
        if cnt[L["locus"]] or old.get(L["locus"]):
            print(f"{L['locus']}\t{L['kind']}\t{cnt[L['locus']]}\t{old.get(L['locus'], '-')}")


if __name__ == "__main__":
    {"loci": cmd_loci, "labels": cmd_labels}[sys.argv[1]]()
```

- [ ] **Step 4: Run the tests to verify they pass**

Run: `cd /mnt/linuxdisk/home/juanfraitu/rustle_m2_soto && python3 -B -m unittest bench/o3_maternal/test_truth.py`
Expected: `Ran 10 tests ... OK`.

- [ ] **Step 5: Build the truth loci (this step is run once per reference; here `mat`) and check them against the prereg**

```bash
# run environment block first (O3_REF=mat)
tools/rlock.sh light python3 bench/o3_maternal/truth.py loci
cat $WR/truth/lrpap1_loci.tsv
```
Expected for `mat` (prereg S3): 35 catalog loci; LRPAP1 absent set = `p12` (`LRPAP1_p12`) plus `c07` as `sex`; every other LRPAP1 locus has a `mat` interval. If a different set comes out, **stop and report to the user** (the prereg was written from the artifact's flags; the PAFs are the authority).

- [ ] **Step 6: Label the reads and cross-check against Amendment 10**

```bash
# run environment block first (O3_REF=mat)
tools/rlock.sh light python3 bench/o3_maternal/truth.py labels
```
Expected: the table lists GWFAM175_B0 near 281, GWFAM205_B0 near 77, `LRPAP1_p12` near 83. Differences against Amendment 10's `express` counts (which used both haplotype BAMs jointly) are expected to be small; a locus differing by more than 25% is reported in the results write-up, not silently accepted.

- [ ] **Step 7: Commit**

```bash
# run environment block first
git add bench/o3_maternal/truth.py bench/o3_maternal/test_truth.py
commit "o3_maternal: absent-locus truth and read labels, symmetric in the reference haplotype"
```

---

### Task 5: Q1 — the fate table (per reference)

**Files:**
- Create: `bench/o3_maternal/fate.py`, `bench/o3_maternal/test_fate.py`
- Outputs: `WR/truth/loci.fa`, `WR/truth/loci.ref.paf`, `WR/fate/fate.tsv`, `WR/fate/fate.json`, `W/unm.txt`

**Interfaces:**
- Consumes: `common.*`, `WR/truth/{loci.tsv,labels.tsv}`, `W/map/reads.<REF>.all.bam`, `W/map/R_unm.{mat,pat}.bam`, `gorilla_haps/<OTHER>.fa`.
- Produces: `fractions(fates, n) -> {'absorbed','unmapped','tied'} | None`; `verdict(fr) -> str`; `fate_rows(recs, labels, paralog) -> {group: {'n','fates','de','reads'}}`; `WR/fate/fate.json` = `{"loci": {locus: {kind, n, fates, fractions, verdict, de_median, paralog: [acc,s,e], paralog_identity, reads: [[read, fate, ref, start, de]]}}, "shared": {...}}` (used by Tasks 9, 11, 12).

- [ ] **Step 1: Write the failing tests**

`bench/o3_maternal/test_fate.py`:

```python
#!/usr/bin/env python3
"""Run: python3 -B -m unittest bench/o3_maternal/test_fate.py"""
import collections
import os
import sys
import unittest

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import fate as F  # noqa: E402
from common import Rec  # noqa: E402


def rec(primary=True, mapq=60, score=1000, qcov=1.0, ref="A", start=0, end=100, de=0.01):
    return Rec(primary, mapq, score, qcov, ref, start, end, de)


class Verdict(unittest.TestCase):
    def fr(self, a, u, t):
        return {"absorbed": a, "unmapped": u, "tied": t}

    def test_refuted(self):
        self.assertEqual(F.verdict(self.fr(0.9, 0.05, 0.05)), "REFUTED")

    def test_partly(self):
        self.assertEqual(F.verdict(self.fr(0.7, 0.2, 0.1)), "PARTLY")

    def test_supported_when_lost_reads_reach_half(self):
        self.assertEqual(F.verdict(self.fr(0.4, 0.3, 0.3)), "SUPPORTED")

    def test_between_bars(self):
        self.assertEqual(F.verdict(self.fr(0.85, 0.0, 0.15)), "BETWEEN_BARS")

    def test_zero_reads_has_no_verdict(self):
        self.assertEqual(F.verdict(F.fractions(collections.Counter(), 0)), "NO_READS")


class Rows(unittest.TestCase):
    def test_counts_and_de(self):
        recs = {"r1": [rec(ref="P", start=10, end=90, de=0.002)],           # absorbed on the nearest paralog
                "r2": [rec(ref="Q", de=0.05)],                                # absorbed elsewhere
                "r3": [rec(mapq=0)],                                           # tied
                "r4": [],                                                      # unmapped
                "r5": [rec(qcov=0.5)],                                         # partial
                "r6": [rec()]}                                                 # shared, not counted for the locus
        labels = {"r1": "L", "r2": "L", "r3": "L", "r4": "L", "r5": "L", "r6": "shared", "r7": "ambiguous"}
        out = F.fate_rows(recs, labels, {"L": ("P", 0, 100)})
        self.assertEqual(out["L"]["n"], 5)
        self.assertEqual(dict(out["L"]["fates"]), {"ABSORBED_NEAREST": 1, "ABSORBED_OTHER": 1, "TIED": 1, "UNMAPPED": 1, "PARTIAL": 1})
        self.assertEqual(sorted(out["L"]["de"]), [0.002, 0.05])
        self.assertEqual(out["shared"]["n"], 1)
        self.assertNotIn("ambiguous", out)

    def test_read_missing_from_the_bam_is_skipped(self):
        self.assertEqual(F.fate_rows({}, {"x": "L"}, {}), {})

    def test_fractions_merge_partial_into_unmapped(self):
        fr = F.fractions(collections.Counter({"UNMAPPED": 1, "PARTIAL": 1, "TIED": 2, "ABSORBED_OTHER": 6}), 10)
        self.assertEqual(fr, {"absorbed": 0.6, "unmapped": 0.2, "tied": 0.2})


if __name__ == "__main__":
    unittest.main()
```

- [ ] **Step 2: Run to verify failure**

Run: `cd /mnt/linuxdisk/home/juanfraitu/rustle_m2_soto && python3 -B -m unittest bench/o3_maternal/test_fate.py`
Expected: `ModuleNotFoundError: No module named 'fate'`.

- [ ] **Step 3: Write `bench/o3_maternal/fate.py`**

```python
#!/usr/bin/env python3
"""Q1 of docs/PREREG_o3_maternal_reference_2026-10-08.md: the fate on the reference haplotype C.REF (env O3_REF, default mat) of the reads of
the loci absent from it.

    O3_REF=mat fate.py fasta    # WR/truth/loci.fa: the C.OTHER sequence of every absent locus
    O3_REF=mat fate.py run      # WR/fate/fate.{tsv,json} (needs WR/truth/loci.ref.paf, see Task 5 step 5)
    fate.py unm                 # descriptive: where the 959 unmapped primaries go on mat / pat (reference-independent) -> W/unm.txt
"""
import collections
import csv
import json
import os
import statistics
import sys

import pysam

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import common as C  # noqa: E402

FATES = ("UNMAPPED", "PARTIAL", "TIED", "ABSORBED_NEAREST", "ABSORBED_OTHER")


def fractions(fates, n):
    """PARTIAL is summed into unmapped (prereg S5). None when the locus has no reads."""
    if n == 0:
        return None
    return {"absorbed": (fates["ABSORBED_NEAREST"] + fates["ABSORBED_OTHER"]) / n,
            "unmapped": (fates["UNMAPPED"] + fates["PARTIAL"]) / n, "tied": fates["TIED"] / n}


def verdict(fr):
    """The bar of the 09-23 registration, per LARGE locus. 'unmapped or tied' is read as the two lost classes together:
    >= 50% lost -> SUPPORTED, 20-50% -> PARTLY, >= 80% absorbed with <= 10% unmapped and <= 10% tied -> REFUTED, else BETWEEN_BARS."""
    if fr is None:
        return "NO_READS"
    lost = fr["unmapped"] + fr["tied"]
    if lost >= 0.5:
        return "SUPPORTED"
    if fr["absorbed"] >= 0.8 and fr["unmapped"] <= 0.1 and fr["tied"] <= 0.1:
        return "REFUTED"
    if 0.2 <= lost < 0.5:
        return "PARTLY"
    return "BETWEEN_BARS"


def fate_rows(recs, labels, paralog):
    """recs: {read: [Rec]} on mat (accessions); labels: {read: locus|'shared'|'ambiguous'}; paralog: {locus: (acc, s, e)}.
    -> {group: {'n', 'fates': Counter, 'de': [floats of absorbed reads], 'reads': [[read, fate, ref, start, de]]}}"""
    out = {}
    for n, g in labels.items():
        if g == "ambiguous" or n not in recs:
            continue
        r = recs[n]
        f = C.classify_fate(r, paralog.get(g))
        d = out.setdefault(g, {"n": 0, "fates": collections.Counter(), "de": [], "reads": []})
        d["n"] += 1
        d["fates"][f] += 1
        p = r[0] if r else None
        if f.startswith("ABSORBED") and p.de is not None:
            d["de"].append(p.de)
        d["reads"].append([n, f, p.ref if p else None, p.start if p else None, p.de if p else None])
    return out


def cmd_fasta():
    fa = pysam.FastaFile(C.HAP_FA.format(C.OTHER))
    with open(f"{C.WR}/truth/loci.fa", "w") as o:
        for r in csv.DictReader(open(f"{C.WR}/truth/loci.tsv"), delimiter="\t"):
            o.write(f">{r['locus']}\n{fa.fetch(r['chrom'], int(r['start']), int(r['end'])).upper()}\n")
    print("loci.fa:", sum(1 for ln in open(f"{C.WR}/truth/loci.fa") if ln[0] == ">"), "sequences")


def cmd_run():
    al = C.alias()
    loci_ = {r["locus"]: r for r in csv.DictReader(open(f"{C.WR}/truth/loci.tsv"), delimiter="\t")}
    labels = {r["read"]: r["label"] for r in csv.DictReader(open(f"{C.WR}/truth/labels.tsv"), delimiter="\t")}
    best = C.best_hits(f"{C.WR}/truth/loci.ref.paf", al)
    paralog = {k: tuple(v[1:4]) for k, v in best.items()}
    ident = {k: v[4] for k, v in best.items()}
    recs = C.read_records(f"{C.W}/map/reads.{C.REF}.all.bam", al)
    res = fate_rows(recs, labels, paralog)
    out = {"loci": {}, "shared": None}
    os.makedirs(f"{C.WR}/fate", exist_ok=True)
    with open(f"{C.WR}/fate/fate.tsv", "w") as o:
        o.write("group\tkind\tn\t" + "\t".join(FATES) + "\tde_median\tparalog\tparalog_identity\tverdict\n")
        for g in list(loci_) + ["shared"]:
            d = res.get(g, {"n": 0, "fates": collections.Counter(), "de": [], "reads": []})
            fr = fractions(d["fates"], d["n"])
            kind = loci_[g]["kind"] if g in loci_ else "control"
            large = g in loci_ and kind != "sex" and d["n"] >= 20
            v = verdict(fr) if large else ""
            dm = statistics.median(d["de"]) if d["de"] else None
            row = dict(kind=kind, n=d["n"], fates={f: d["fates"][f] for f in FATES}, fractions=fr, verdict=v, de_median=dm,
                       paralog=list(paralog[g]) if g in paralog else None, paralog_identity=ident.get(g), reads=d["reads"])
            if g in loci_:
                out["loci"][g] = row
            else:
                out["shared"] = row
            o.write("\t".join(str(x) for x in [g, kind, d["n"], *[d["fates"][f] for f in FATES],
                                               "" if dm is None else f"{dm:.4f}", ":".join(str(x) for x in paralog.get(g, ())),
                                               "" if g not in ident else f"{ident[g]:.4f}", v]) + "\n")
    json.dump(out, open(f"{C.WR}/fate/fate.json", "w"))
    print(open(f"{C.WR}/fate/fate.tsv").read())
    print("selection: every read of the 34-family and LRPAP1 sets was mapped on `_pri` first, so UNMAPPED/PARTIAL there means "
          f"'mapped on _pri, lost on {C.REF}'; reads unmapped on _pri are only in R_unm (fate.py unm).")


def cmd_unm():
    al = C.alias()
    lines = []
    for hap in ("mat", "pat"):
        recs = C.read_records(f"{C.W}/map/R_unm.{hap}.bam", al)
        mapped = sum(1 for rs in recs.values() if any(r.primary and r.qcov >= C.COV_MIN for r in rs))
        lines.append(f"R_unm on {hap}: {len(recs)} reads, {mapped} with a primary at query coverage >= {C.COV_MIN}")
    open(f"{C.W}/unm.txt", "w").write("\n".join(lines) + "\n")
    print("\n".join(lines))


if __name__ == "__main__":
    {"fasta": cmd_fasta, "run": cmd_run, "unm": cmd_unm}[sys.argv[1]]()
```

- [ ] **Step 4: Run the tests to verify they pass**

Run: `cd /mnt/linuxdisk/home/juanfraitu/rustle_m2_soto && python3 -B -m unittest bench/o3_maternal/test_fate.py`
Expected: `Ran 8 tests ... OK`.

- [ ] **Step 5: Nearest reference paralog of every locus (heavy, ~1 min; once per reference, here `mat`)**

```bash
# run environment block first (O3_REF=mat)
tools/rlock.sh light python3 bench/o3_maternal/fate.py fasta
tools/rlock.sh heavy bash -c "minimap2 -c -x asm20 --secondary=yes -N 50 -p 0.1 -t 4 $IDXREF $WR/truth/loci.fa > $WR/truth/loci.ref.paf 2> $WR/truth/loci.ref.log"
wc -l $WR/truth/loci.ref.paf
```
Expected: a non-empty PAF (the command is the 10-04 LRPAP1 one with `-p 0.1`; record `-N 50 -p 0.1` next to any identity quoted).

- [ ] **Step 6: Run the fate table**

```bash
# run environment block first (O3_REF=mat)
tools/rlock.sh light python3 bench/o3_maternal/fate.py run
tools/rlock.sh light python3 bench/o3_maternal/fate.py unm        # reference-independent; run once
```
Expected: `fate.tsv` has one row per locus plus `shared`; LARGE loci carry a verdict. Compare with the prediction in prereg S5 (REFUTED at every LARGE locus; GWFAM175_B0 absorbed at `de` ~0.066; p12 absorbed on the p14 ortholog). A prediction that fails is **reported as failed**, never reworded.

- [ ] **Step 7: Commit**

```bash
# run environment block first
git add bench/o3_maternal/fate.py bench/o3_maternal/test_fate.py
commit "o3_maternal: Q1 fate table (unmapped / tied / absorbed) with the registered bar"
```

---

### Task 6: Chain inputs (per reference)

**Files:**
- Create: `bench/o3_maternal/chain_inputs.py`, `bench/o3_maternal/test_chain_inputs.py`, `bench/o3_maternal/iso_batch.sh`
- Outputs: `WR/isoc/{panel.json,labels.tsv,scored.fa,R0.bam}`, `WR/inhouse/{panel.json,labels.tsv,R0.bam,M.copies.tsv,M.copies.fa,M.regions}`, `W/<REF>.idx.fa` (+ `.fai`)

**Interfaces:**
- Consumes: `common.*`, `extract_reads.R34`, `TRUTH/refabsent/{copies.<REF>.paf,fams_bonly.txt,scored.fa,labels.tsv}`, `LRP/{copies8,partial3}.<REF>.paf`, `W/reads/R_LRP.{fa,names.tsv}`, `W/map/reads.<REF>.all.bam`, `rna_allele/panel_to_copies.py`.
- Produces: `merge_loci(hits, gap=5000) -> [(chrom, s, e)]`; `build_panel(copy_loci, fams) -> [dict(fam, mask, keep)]` (panel layout of `control_test.panel`: each copy `[chrom, start, end, gene]`, chromosome = **index name** `chrN_<REF>_hsaX`); `name_map(fai_rows, sq) -> {fasta name: index name}`; a work dir `WR/isoc` that `control_test.py net --w WR/isoc --l WR/isoc` accepts unchanged.

- [ ] **Step 1: Write the failing tests**

`bench/o3_maternal/test_chain_inputs.py`:

```python
#!/usr/bin/env python3
"""Run: python3 -B -m unittest bench/o3_maternal/test_chain_inputs.py"""
import os
import sys
import unittest

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import chain_inputs as I  # noqa: E402


class Loci(unittest.TestCase):
    def test_merge_within_gap(self):
        self.assertEqual(I.merge_loci([("c1", 100, 200), ("c1", 4000, 4500), ("c1", 90000, 91000), ("c2", 10, 20)]),
                         [("c1", 100, 4500), ("c1", 90000, 91000), ("c2", 10, 20)])

    def test_merge_is_order_independent(self):
        a = [("c1", 4000, 4500), ("c1", 100, 200)]
        self.assertEqual(I.merge_loci(a), I.merge_loci(list(reversed(a))))

    def test_panel_mask_first_then_keep(self):
        p = I.build_panel({"F1": [("c1", 100, 200), ("c2", 5, 9)], "F2": [("c1", 1, 2)]}, ["F1", "F2", "F3"])
        self.assertEqual([x["fam"] for x in p], ["F1", "F2"])          # F3 has no mat locus: dropped
        self.assertEqual(p[0]["mask"], ["c1", 100, 200, "F1:0"])
        self.assertEqual(p[0]["keep"], [["c2", 5, 9, "F1:1"]])
        self.assertEqual(p[1]["keep"], [])


class Names(unittest.TestCase):
    def test_name_map_pairs_by_order_and_length(self):
        self.assertEqual(I.name_map([("CM1", 10), ("CM2", 20)], [("chr1_mat_hsa1", 10), ("chr2_mat_hsa2", 20)]),
                         {"CM1": "chr1_mat_hsa1", "CM2": "chr2_mat_hsa2"})

    def test_name_map_refuses_a_length_mismatch(self):
        with self.assertRaises(AssertionError):
            I.name_map([("CM1", 10)], [("chr1_mat_hsa1", 11)])


if __name__ == "__main__":
    unittest.main()
```

- [ ] **Step 2: Run to verify failure**

Run: `cd /mnt/linuxdisk/home/juanfraitu/rustle_m2_soto && python3 -B -m unittest bench/o3_maternal/test_chain_inputs.py`
Expected: `ModuleNotFoundError: No module named 'chain_inputs'`.

- [ ] **Step 3: Write `bench/o3_maternal/chain_inputs.py` and `iso_batch.sh`**

```python
#!/usr/bin/env python3
"""Inputs of the recovery chain with C.REF (env O3_REF, default mat) as the reference (docs/PREREG_o3_maternal_reference_2026-10-08.md section 6).

    O3_REF=mat chain_inputs.py panel     # WR/isoc/{panel.json,labels.tsv,scored.fa,R0.bam}; WR/inhouse/{panel.json,labels.tsv,R0.bam}
    O3_REF=mat chain_inputs.py fasta     # W/<REF>.idx.fa: gorilla_haps/<REF>.fa with the haplotype-index sequence names (chrN_<REF>_hsaX), so that
                                         # BAM, FASTA and .mmi agree for panel_to_copies.py and o3_candidates
"""
import collections
import csv
import json
import os
import shutil
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import common as C  # noqa: E402
import extract_reads as E  # noqa: E402

REF = f"{C.TRUTH}/refabsent"


def merge_loci(hits, gap=5000):
    """hits: [(chrom, start, end)] -> merged [(chrom, start, end)] (hits within `gap` bp on one chromosome join)"""
    out = []
    for c, s, e in sorted(hits):
        if out and out[-1][0] == c and s <= out[-1][2] + gap:
            out[-1][2] = max(out[-1][2], e)
        else:
            out.append([c, s, e])
    return [tuple(x) for x in out]


def build_panel(copy_loci, fams):
    """copy_loci: {family: [(chrom, start, end)]} (merged or not); control_test.panel layout, mask = first copy, keep = the rest;
    families without a reference locus are dropped"""
    out = []
    for f in fams:
        cps = [(c, s, e, f"{f}:{k}") for k, (c, s, e) in enumerate(merge_loci(copy_loci.get(f, [])))]
        if cps:
            out.append(dict(fam=f, mask=list(cps[0]), keep=[list(x) for x in cps[1:]]))
    return out


def name_map(fai_rows, sq):
    """fai_rows: [(name, length)] of the FASTA; sq: [(name, length)] of the BAM header in index order -> {fasta name: index name}"""
    assert len(fai_rows) == len(sq), (len(fai_rows), len(sq))
    for (a, la), (b, lb) in zip(fai_rows, sq):
        assert la == lb, (a, la, b, lb)
    return {a: b for (a, _), (b, _) in zip(fai_rows, sq)}


def families():
    return [ln.strip() for ln in open(f"{REF}/fams_bonly.txt") if ln.strip()] + ["LRPAP1"]


def ref_copy_loci():
    """family -> raw hits (index names) on the reference haplotype of the family's copies at identity >= .90, coverage >= .80 (-p/-N of the
    source PAFs: copies.<hap>.paf = Amendment 10; LRPAP1 = asm20 -N 50 -p 0.5)"""
    out = collections.defaultdict(list)
    for q, hs in C.paf_hits(f"{REF}/copies.{C.REF}.paf", {}).items():
        out[q.split(":")[0]].extend((h[0], h[1], h[2]) for h in hs)
    for f in (f"copies8.{C.REF}.paf", f"partial3.{C.REF}.paf"):
        for _q, hs in C.paf_hits(f"{E.LRP}/{f}", {}).items():
            out["LRPAP1"].extend((h[0], h[1], h[2]) for h in hs)
    return out


def cmd_panel():
    fams = families()
    panel = build_panel(ref_copy_loci(), fams)
    missing = sorted(set(fams) - {p["fam"] for p in panel})
    for d in ("isoc", "inhouse"):
        os.makedirs(f"{C.WR}/{d}", exist_ok=True)
        json.dump(panel, open(f"{C.WR}/{d}/panel.json", "w"), indent=0)
        bam = f"{C.W}/map/reads.{C.REF}.all.bam"
        for ext in ("", ".bai"):
            link = f"{C.WR}/{d}/R0.bam{ext}"
            if os.path.lexists(link):
                os.remove(link)
            os.symlink(bam + ext, link)
    r34 = list(csv.DictReader(open(E.R34), delimiter="\t"))
    seen = {r["read"] for r in r34}
    with open(f"{C.WR}/isoc/labels.tsv", "w") as o:
        o.write("read\tfamily\trole\tcopy\n")
        for r in r34:
            o.write(f"{r['read']}\t{r['family']}\t{r['role']}\t{r['copy']}\n")
        new = 0
        for r in csv.DictReader(open(f"{C.W}/reads/R_LRP.names.tsv"), delimiter="\t"):
            if r["read"] not in seen:
                o.write(f"{r['read']}\tLRPAP1\tS\t{r['cids']}\n")
                new += 1
    with open(f"{C.WR}/isoc/scored.fa", "w") as o:
        for src in (f"{REF}/scored.fa", f"{C.W}/reads/R_LRP.fa"):
            with open(src) as f:
                shutil.copyfileobj(f, o)
    shutil.copy(f"{C.WR}/isoc/labels.tsv", f"{C.WR}/inhouse/labels.tsv")
    print(f"reference {C.REF}: families in the panel {len(panel)} of {len(fams)}; dropped (no {C.REF} locus): {missing}; copies per family: "
          f"{collections.Counter(1 + len(p['keep']) for p in panel)}; labels {len(r34)} + {new} LRPAP1-only reads")


def cmd_fasta():
    import pysam
    src = C.HAP_FA.format(C.REF)
    fai = [(r[0], int(r[1])) for r in csv.reader(open(src + ".fai"), delimiter="\t")]
    with pysam.AlignmentFile(f"{C.W}/map/reads.{C.REF}.all.bam") as b:
        sq = [(s["SN"], s["LN"]) for s in b.header.to_dict()["SQ"]]
    m = name_map(fai, sq)
    out = f"{C.W}/{C.REF}.idx.fa"
    with open(src) as f, open(out, "w") as o:
        for ln in f:
            o.write(">" + m[ln[1:].split()[0]] + "\n" if ln[0] == ">" else ln)
    subprocess.run(["samtools", "faidx", out], check=True)
    print(f"{out}: {len(m)} sequences renamed")


if __name__ == "__main__":
    {"panel": cmd_panel, "fasta": cmd_fasta}[sys.argv[1]]()
```

`bench/o3_maternal/iso_batch.sh`:

```bash
#!/bin/bash
# resumable: IsoCon per family until the deadline; skips families with final_candidates.fa. H = the work dir (env, required).
# Copy of TRUTH/refabsent/iso_batch.sh with H from the environment.   H=W/isoc BUDGET=540 tools/rlock.sh heavy bash iso_batch.sh
H=${H:?work dir}; ISO=/home/juanfra/miniforge3/envs/isocon/bin/IsoCon
DEADLINE=$(( $(date +%s) + ${BUDGET:-540} ))
for fa in $(ls -S -r $H/fam/*.fa); do
  f=$(basename $fa .fa); out=$H/iso/$f
  [ -s $out/final_candidates.fa ] && continue
  [ -e $out.empty ] && continue
  [ $(grep -c ">" $fa) -lt 2 ] && { mkdir -p $H/iso; touch $out.empty; continue; }
  now=$(date +%s); [ $now -ge $(( DEADLINE - 30 )) ] && { echo "deadline"; exit 75; }
  rm -rf $out; mkdir -p $H/iso
  timeout $(( DEADLINE - now )) $ISO pipeline -fl_reads $fa -outfolder $out --nr_cores 4 > $out.log 2>&1
  rc=$?
  if [ $rc -ne 0 ]; then echo "$f rc=$rc"; [ $rc -eq 124 ] && { rm -rf $out; exit 75; }; touch $out.empty; fi
done
echo ALL_DONE
```

- [ ] **Step 4: Run the tests to verify they pass**

Run: `cd /mnt/linuxdisk/home/juanfraitu/rustle_m2_soto && python3 -B -m unittest bench/o3_maternal/test_chain_inputs.py`
Expected: `Ran 5 tests ... OK`.

- [ ] **Step 5: Build the panel and inputs, then the renamed FASTA (once per reference, here `mat`)**

```bash
# run environment block first (O3_REF=mat)
tools/rlock.sh light python3 bench/o3_maternal/chain_inputs.py panel
df -h /mnt/linuxdisk | tail -1                                   # >= 10 GB free before the next line
tools/rlock.sh heavy python3 bench/o3_maternal/chain_inputs.py fasta
```
Expected: `families in the panel` 35 (or fewer with the dropped list printed; a dropped family has no reference copy at identity >= .90 and is excluded from the chain with a note in the write-up); `<REF>.idx.fa: 225 sequences renamed` (`mat`; the `pat` assembly has 24). The `name_map` assertion firing means the index order differs from the FASTA order: stop and report.

- [ ] **Step 6: Build the in-house copies table (uses the registered `panel_to_copies.py`)**

```bash
# run environment block first (O3_REF=mat)
tools/rlock.sh light python3 bench/rna_allele/panel_to_copies.py --all --panel $WR/inhouse/panel.json --bam $WR/inhouse/R0.bam --fasta $W/$REF.idx.fa --out $WR/inhouse/M
head -3 $WR/inhouse/M.copies.tsv | cut -c1-200; grep -c '>' $WR/inhouse/M.copies.fa
```
Expected: `M.copies.tsv` rows = total copies of the panel; `M.copies.fa` has the same number of records, none all-N.

- [ ] **Step 7: Commit**

```bash
# run environment block first
git add bench/o3_maternal/chain_inputs.py bench/o3_maternal/test_chain_inputs.py bench/o3_maternal/iso_batch.sh
commit "o3_maternal: chain inputs for a reference haplotype (panel, labels, renamed FASTA)"
```

---

### Task 7: Arm I — IsoCon, truth-free, to candidates (per reference)

**Files:**
- Create: `bench/o3_maternal/adapt_isocon.py`
- Outputs: `WR/isoc/{fam/*.fa,iso/*,outputs.fa,outputs.pri.paf,outputs.other.paf,contigs*.fa,contigs.tsv,merge/*,cands.tsv,cands.fa,cands.ref.paf,cands.other.paf}`

**Interfaces:**
- Consumes: Task 6 outputs; `rna_allele/control_test.py net|outputs|contigs`, `control_test.comps_at`; `merge_test`.
- Produces: `WR/isoc/cands.tsv` (columns `candidate family n_transcripts contigs`; `contigs` = comma-separated contig names), `WR/isoc/cands.fa` (the new-copy contigs), the hit tables `cands.{ref,other}.paf` (contigs vs the reference / the truth haplotype) used by Task 9, and `outputs.pri.paf` (ALL IsoCon outputs vs the reference — the file name `control_test.contigs` expects) and `outputs.other.paf` (all outputs vs the truth haplotype) used by Tasks 9 and 11.

- [ ] **Step 1: Write `bench/o3_maternal/adapt_isocon.py`**

```python
#!/usr/bin/env python3
"""IsoCon arm -> the common candidate tables.   adapt_isocon.py <workdir>   (workdir = W/isoc, after control_test.py contigs)

cands.tsv: candidate, family, n_transcripts (= contigs in the merge component, Amendment 8's rule), contigs; cands.fa: the new-copy contigs."""
import collections
import os
import shutil
import sys
import types

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.join(HERE, "..", "rna_allele"))
import common as C  # noqa: E402
import control_test  # noqa: E402


def main(w):
    comp, fam_of = control_test.comps_at(types.SimpleNamespace(w=w), C.DELTA)
    members = collections.defaultdict(list)
    for c, cid in comp.items():
        members[cid].append(c)
    with open(f"{w}/cands.tsv", "w") as o:
        o.write("candidate\tfamily\tn_transcripts\tcontigs\n")
        for cid, cs in sorted(members.items()):
            o.write(f"{cid}\t{fam_of[cs[0]]}\t{len(cs)}\t{','.join(sorted(cs))}\n")
    if os.path.exists(f"{w}/contigs_L.fa"):
        shutil.copy(f"{w}/contigs_L.fa", f"{w}/cands.fa")
    else:
        open(f"{w}/cands.fa", "w").close()
    flagged = sum(1 for cs in members.values() if len(cs) >= 2)
    print(f"candidates {len(members)} in {len({fam_of[cs[0]] for cs in members.values()})} families; flagged (>= 2 transcripts) {flagged}")


if __name__ == "__main__":
    main(sys.argv[1])
```

- [ ] **Step 2: Net and IsoCon input (light; once per reference, here `mat`)**

```bash
# run environment block first (O3_REF=mat)
tools/rlock.sh light python3 bench/rna_allele/control_test.py net --w $WR/isoc --l $WR/isoc
```
Expected: `IsoCon inputs: <= 35 families, total ~20-25k, median ..., max 1000`.

- [ ] **Step 3: IsoCon per family, resumable foreground calls**

```bash
# run environment block first (O3_REF=mat)
H=$WR/isoc BUDGET=540 tools/rlock.sh heavy bash bench/o3_maternal/iso_batch.sh
```
Repeat the identical command until it prints `ALL_DONE` (exit 75 + `deadline` = call again; any `<fam> rc=<n>` line is kept and reported). Check: `ls $WR/isoc/iso/*/final_candidates.fa | wc -l`.

- [ ] **Step 4: Outputs, flag step against the reference, contigs (link + merge)**

```bash
# run environment block first (O3_REF=mat)
python3 bench/rna_allele/control_test.py outputs --w $WR/isoc --l $WR/isoc
tools/rlock.sh heavy bash -c "minimap2 -c -x splice:hq -uf -N 20 -t 4 $IDXREF $WR/isoc/outputs.fa > $WR/isoc/outputs.pri.paf 2> $WR/isoc/outputs.pri.log"
tools/rlock.sh light python3 bench/rna_allele/control_test.py contigs --w $WR/isoc --l $WR/isoc
python3 bench/o3_maternal/adapt_isocon.py $WR/isoc
```
(`outputs.pri.paf` is the file name `control_test.contigs` reads; here it holds the alignment to the run's **reference**.) Expected: `outputs N; flagged ...; linked back ...; kept as new copies ...`, then the `candidates ... flagged ...` line.

- [ ] **Step 5: Hits of the candidates and of all outputs on both haplotypes (heavy, serial)**

```bash
# run environment block first (O3_REF=mat)
tools/rlock.sh heavy bash -c "minimap2 -c -x splice:hq -uf -N 20 -t 4 $IDXREF $WR/isoc/cands.fa > $WR/isoc/cands.ref.paf 2> $WR/isoc/cands.ref.log"
tools/rlock.sh heavy bash -c "minimap2 -c -x splice:hq -uf -N 20 -t 4 $IDXOTH $WR/isoc/cands.fa > $WR/isoc/cands.other.paf 2> $WR/isoc/cands.other.log"
tools/rlock.sh heavy bash -c "minimap2 -c -x splice:hq -uf -N 20 -t 4 $IDXOTH $WR/isoc/outputs.fa > $WR/isoc/outputs.other.paf 2> $WR/isoc/outputs.other.log"
```
(An empty `cands.fa` makes minimap2 print nothing: skip those two mappings and let Task 9 report R1 FAIL / NOT TESTABLE.)

- [ ] **Step 6: Commit**

```bash
# run environment block first
git add bench/o3_maternal/adapt_isocon.py
commit "o3_maternal: IsoCon arm adapter (candidates from the registered merge components)"
```

---

### Task 8: Arm H — the in-house `o3_candidates` stage (per reference)

**Files:**
- Create: `bench/o3_maternal/adapt_inhouse.py`
- Outputs: `WR/inhouse/cand_g{0..3}.*`, `WR/inhouse/cands.tsv`, `cands_all.tsv`, `cands.fa`, `cands.{ref,other}.paf`

**Interfaces:**
- Consumes: Task 6 outputs; binary `/mnt/linuxdisk/home/juanfraitu/rustle_target/release/o3_candidates`.
- Produces: `WR/inhouse/cands.tsv` (same columns as Task 7; `n_transcripts` = the stage's `n_clusters`; only rows with `flagged == 1`; `contigs` = the candidate name), `cands_all.tsv` (every stage candidate: `candidate family n_clusters flagged nearest_locus d`; used by Task 11), `cands.fa`, `cands.{ref,other}.paf`.

- [ ] **Step 1: Provenance of the binary (prereg S6: "as shipped at main b29afa55, Amendment 15 not applied"); once**

```bash
# run environment block first
ls -l --time-style=full-iso /mnt/linuxdisk/home/juanfraitu/rustle_target/release/o3_candidates
git diff --stat b29afa55 HEAD -- src | tail -15
git log --format='%h %cI %s' b29afa55..HEAD -- src | cut -c1-110
```
Decision rule (no new constants): if no commit after `b29afa55` changed `o3_candidates` behaviour (consolidation/move commits only, Amendment 15 not implemented), use the existing binary and record `sha256sum` of it in the write-up. If behaviour changed, **stop and ask the user** whether to build `b29afa55` in a separate worktree (cargo to a file, `CARGO_TARGET_DIR=/mnt/linuxdisk/home/juanfraitu/rustle_target_b29`, ~1 GB) or accept HEAD with disclosure.

- [ ] **Step 2: Write `bench/o3_maternal/adapt_inhouse.py`**

```python
#!/usr/bin/env python3
"""In-house arm -> the common candidate tables.   adapt_inhouse.py <workdir>   (workdir = W/inhouse, after the o3_candidates batches)

Reads cand_g*.candidates.tsv (columns family candidate n_clusters n_reads flagged union_len nearest_locus d n_net n_used) and
cand_g*.contigs.fa; also writes cands_all.tsv. A candidate is a new-copy candidate iff flagged == 1 (distance to the nearest locus > delta); n_transcripts = n_clusters."""
import csv
import glob
import os
import sys


def read_fa(path):
    seq, cur = {}, None
    for ln in open(path):
        if ln[0] == ">":
            cur = ln[1:].strip()
            seq[cur] = []
        elif cur:
            seq[cur].append(ln.strip())
    return {k: "".join(v) for k, v in seq.items()}


def main(w):
    rows, seqs = [], {}
    for f in sorted(glob.glob(f"{w}/cand_g*.candidates.tsv")):
        rows += list(csv.DictReader(open(f), delimiter="\t"))
        fa = f.replace(".candidates.tsv", ".contigs.fa")
        if os.path.exists(fa):
            seqs.update(read_fa(fa))
    new = [r for r in rows if r["flagged"] == "1"]
    with open(f"{w}/cands_all.tsv", "w") as o:        # every stage candidate, flagged or not (used by compare.py: "found in the reference")
        o.write("candidate\tfamily\tn_clusters\tflagged\tnearest_locus\td\n")
        for r in rows:
            o.write(f"{r['candidate']}\t{r['family']}\t{r['n_clusters']}\t{r['flagged']}\t{r['nearest_locus']}\t{r['d']}\n")
    with open(f"{w}/cands.tsv", "w") as o:
        o.write("candidate\tfamily\tn_transcripts\tcontigs\n")
        for r in new:
            o.write(f"{r['candidate']}\t{r['family']}\t{r['n_clusters']}\t{r['candidate']}\n")
    with open(f"{w}/cands.fa", "w") as o:
        for r in new:
            if r["candidate"] in seqs:
                o.write(f">{r['candidate']}\n{seqs[r['candidate']]}\n")
    print(f"stage candidates {len(rows)} in {len({r['family'] for r in rows})} families; new-copy (flagged==1) {len(new)}; "
          f">= 2 clusters {sum(1 for r in new if int(r['n_clusters']) >= 2)}")


if __name__ == "__main__":
    main(sys.argv[1])
```

- [ ] **Step 3: Run the stage in 4 foreground batches (each < 10 min under `rlock heavy`; once per reference, here `mat`)**

```bash
# run environment block first (O3_REF=mat)
O3=/mnt/linuxdisk/home/juanfraitu/rustle_target/release/o3_candidates
python3 - <<'PY'
import json, os
wr = f"/mnt/linuxdisk/tmp/o3_mat/{os.environ['O3_REF']}/inhouse"
fams = [p["fam"] for p in json.load(open(f"{wr}/panel.json"))]
n = (len(fams) + 3) // 4
for i in range(4):
    open(f"{wr}/batch{i}.txt", "w").write(",".join(fams[i * n:(i + 1) * n]))
PY
for i in 0 1 2 3; do
  tools/rlock.sh heavy $O3 --bam $WR/inhouse/R0.bam --fasta $W/$REF.idx.fa --copies $WR/inhouse/M.copies.tsv --copies-fa $WR/inhouse/M.copies.fa \
    --index $IDXREF --out $WR/inhouse/cand_g$i --families "$(cat $WR/inhouse/batch$i.txt)" > $WR/inhouse/cand_g$i.log 2>&1 || echo "batch $i rc=$?"
done
ls $WR/inhouse/cand_g*.candidates.tsv
```
A batch that hits the 600 s timeout (`rc=124`): split it into smaller `--families` lists and rerun only that batch; never raise the timeout past the foreground limit.

- [ ] **Step 4: Adapt, map to both haplotypes, free the disk**

```bash
# run environment block first (O3_REF=mat)
python3 bench/o3_maternal/adapt_inhouse.py $WR/inhouse
tools/rlock.sh heavy bash -c "minimap2 -c -x splice:hq -uf -N 20 -t 4 $IDXREF $WR/inhouse/cands.fa > $WR/inhouse/cands.ref.paf 2> $WR/inhouse/cands.ref.log"
tools/rlock.sh heavy bash -c "minimap2 -c -x splice:hq -uf -N 20 -t 4 $IDXOTH $WR/inhouse/cands.fa > $WR/inhouse/cands.other.paf 2> $WR/inhouse/cands.other.log"
rm -f $W/$REF.idx.fa $W/$REF.idx.fa.fai
```
Expected: the adapter line, and non-empty PAFs unless the stage flagged nothing (then Task 9 reports R1 FAIL).

- [ ] **Step 5: Commit**

```bash
# run environment block first
git add bench/o3_maternal/adapt_inhouse.py
commit "o3_maternal: in-house arm adapter (o3_candidates as shipped, Amendment 15 not applied)"
```

---

### Task 9: The scorer — R1-R4 and the p12 line (per reference)

**Files:**
- Create: `bench/o3_maternal/score.py`, `bench/o3_maternal/test_score.py`
- Outputs: `WR/score/{isocon,inhouse}.json` and printed tables

**Interfaces:**
- Consumes: `common.*`, `WR/truth/{loci.tsv,lrpap1_loci.tsv}`, `WR/fate/fate.json`, `WR/<arm>/{cands.tsv,cands.ref.paf,cands.other.paf,panel.json}`, `WR/isoc/{outputs.pri.paf,outputs.other.paf}`.
- Produces:
  - `recovered_loci(contigs, other_best, loci_, min_score=0.999) -> set(locus)`
  - `candidate_class(contigs, ref_best, other_best, recovers) -> 'a_recovered'|'ref'|'b_other'|'c_unmatched'`
  - `max_matching(edges, n_left) -> int` (Kuhn)
  - `evaluate(cands, loci_, fam_net, ref_best, other_best, floor=2) -> dict(R1, R2, R3, R4, rows)` (`rows[i]` has `candidate, family, n_transcripts, contigs, recovers, cls`)
  - `type_hits(mat_best, pat_best, targets, min_score=0.999) -> {target: [sequence]}`
  - JSON per arm with keys `R1, R2, R3, R4, rows, candidates` (and `p12` for IsoCon when `O3_REF=mat`) consumed by Tasks 11-12.

- [ ] **Step 1: Write the failing tests**

`bench/o3_maternal/test_score.py`:

```python
#!/usr/bin/env python3
"""Run: python3 -B -m unittest bench/o3_maternal/test_score.py"""
import os
import sys
import unittest

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import score as S  # noqa: E402

DELTA = 0.00958


def L(locus, kind="catalog", family="F1", chrom="CM1", start=100, end=200, n=50, ident=0.91):
    return dict(locus=locus, kind=kind, family=family, chrom=chrom, start=start, end=end, n=n, ident=ident)


def cand(name, fam, n, contigs=None):
    return dict(candidate=name, family=fam, n_transcripts=n, contigs=contigs or [name])


def hit(score, acc, s, e):
    return (score, acc, s, e, score)


class Recovery(unittest.TestCase):
    def test_recovers_by_overlap_at_the_floor(self):
        pat = {"c1": hit(0.9995, "CM1", 120, 180)}
        self.assertEqual(S.recovered_loci(["c1"], pat, [L("A")]), {"A"})

    def test_below_the_floor_or_elsewhere_does_not(self):
        self.assertEqual(S.recovered_loci(["c1"], {"c1": hit(0.998, "CM1", 120, 180)}, [L("A")]), set())
        self.assertEqual(S.recovered_loci(["c1"], {"c1": hit(0.9995, "CM2", 120, 180)}, [L("A")]), set())

    def test_class(self):
        self.assertEqual(S.candidate_class(["c"], {}, {}, {"A"}), "a_recovered")
        self.assertEqual(S.candidate_class(["c"], {"c": hit(1.0, "M", 1, 2)}, {"c": hit(1.0, "P", 1, 2)}, set()), "ref")
        self.assertEqual(S.candidate_class(["c"], {"c": hit(0.99, "M", 1, 2)}, {"c": hit(1.0, "P", 1, 2)}, set()), "b_other")
        self.assertEqual(S.candidate_class(["c"], {"c": hit(0.95, "M", 1, 2)}, {}, set()), "c_unmatched")
        self.assertEqual(S.candidate_class(["c"], {}, {}, set()), "c_unmatched")


class Matching(unittest.TestCase):
    def test_kuhn(self):
        self.assertEqual(S.max_matching([(0, 0), (1, 0), (1, 1)], 2), 2)
        self.assertEqual(S.max_matching([(0, 0), (1, 0)], 2), 1)
        self.assertEqual(S.max_matching([], 0), 0)


class Evaluate(unittest.TestCase):
    loci = [L("BIG", n=281, ident=0.91), L("SMALL", n=13, ident=0.92), L("NEAR", n=77, ident=0.995),
            L("SEX", kind="sex", family="LRPAP1", chrom="CMY", n=30, ident=0.98),
            L("LRPAP1_p12", kind="lrpap1", family="LRPAP1", chrom="CM12", n=83, ident=0.9984)]
    fams = ["F1", "F2", "F3", "LRPAP1"]

    def run_eval(self, cands, pat):
        return S.evaluate(cands, self.loci, self.fams, {}, pat)

    def test_r1_passes_when_the_big_beyond_delta_locus_is_recovered(self):
        out = self.run_eval([cand("k1", "F1", 3)], {"k1": hit(1.0, "CM1", 120, 180)})
        self.assertEqual(out["R1"]["verdict"], "PASS")
        self.assertEqual(out["R1"]["targets"], ["BIG"])

    def test_empty_candidates_fail_r1(self):
        out = self.run_eval([], {})
        self.assertEqual(out["R1"]["verdict"], "FAIL")
        self.assertEqual(out["R4"]["flagged"], 0)

    def test_one_transcript_is_not_a_flag(self):
        out = self.run_eval([cand("k1", "F1", 1)], {"k1": hit(1.0, "CM1", 120, 180)})
        self.assertEqual(out["R1"]["verdict"], "FAIL")

    def test_r2_flags_a_recovered_within_delta_locus(self):
        pat = {"k1": hit(1.0, "CM1", 120, 180)}
        loci = [L("NEAR", ident=0.995, n=77)]
        out = S.evaluate([cand("k1", "F1", 2)], loci, self.fams, {}, pat)
        self.assertEqual(out["R2"]["verdict"], "FAIL")

    def test_sex_locus_excluded(self):
        out = self.run_eval([cand("k1", "LRPAP1", 2)], {"k1": hit(1.0, "CMY", 120, 180)})
        self.assertNotIn("SEX", out["R1"]["targets"])
        self.assertEqual(out["R4"]["expressed"], 4)            # BIG, SMALL, NEAR, p12; SEX not counted

    def test_p12_not_in_r1_or_r2(self):
        out = self.run_eval([cand("k1", "LRPAP1", 2)], {"k1": hit(1.0, "CM12", 120, 180)})
        self.assertNotIn("LRPAP1_p12", out["R1"]["targets"])
        self.assertEqual(out["R2"]["verdict"], "PASS")

    def test_r3_counts_false_flags_in_families_without_a_truth_locus(self):
        out = self.run_eval([cand("k1", "F2", 2)], {"k1": hit(1.0, "CM9", 1, 9)})
        self.assertEqual(out["R3"]["denominator"], 2)          # F2, F3 (F1 and LRPAP1 hold expressed loci)
        self.assertEqual(out["R3"]["false_families"], ["F2"])
        self.assertEqual(out["R3"]["fraction"], 0.5)


class Types(unittest.TestCase):
    def test_distinct_targets(self):
        mat = {"s1": hit(1.0, "M14", 10, 90)}
        pat = {"s2": hit(1.0, "P12", 10, 90), "s3": hit(0.9, "P14", 10, 90)}
        t = {"p12@pat": ("pat", "P12", 0, 100), "p14@pat": ("pat", "P14", 0, 100), "p14@mat": ("mat", "M14", 0, 100)}
        self.assertEqual(S.type_hits(mat, pat, t), {"p12@pat": ["s2"], "p14@pat": [], "p14@mat": ["s1"]})


if __name__ == "__main__":
    unittest.main()
```

- [ ] **Step 2: Run to verify failure**

Run: `cd /mnt/linuxdisk/home/juanfraitu/rustle_m2_soto && python3 -B -m unittest bench/o3_maternal/test_score.py`
Expected: `ModuleNotFoundError: No module named 'score'`.

- [ ] **Step 3: Write `bench/o3_maternal/score.py`**

```python
#!/usr/bin/env python3
"""R1-R4 and the p12 line of docs/PREREG_o3_maternal_reference_2026-10-08.md section 6, per arm and per run (reference C.REF, env O3_REF).

    O3_REF=mat score.py isocon | inhouse   # reads WR/<dir>/{cands.tsv,cands.ref.paf,cands.other.paf}; writes WR/score/<arm>.json, prints verdicts

Loci passed to evaluate() are dicts: locus, kind (catalog|lrpap1|sex), family, chrom (accession on the TRUTH haplotype C.OTHER), start, end,
n (reads), ident (identity of the nearest reference paralog). Truth is used only here."""
import csv
import json
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import common as C  # noqa: E402

MIN_SCORE = 0.999       # registered recovery hit: identity x coverage


def recovered_loci(contigs, other_best, loci_, min_score=MIN_SCORE):
    """loci whose interval on the truth haplotype a contig's best hit there (score >= min_score) overlaps"""
    out = set()
    for c in contigs:
        h = other_best.get(c)
        if h and h[0] >= min_score:
            for L in loci_:
                if h[1] == L["chrom"] and h[2] < L["end"] and L["start"] < h[3]:
                    out.add(L["locus"])
    return out


def candidate_class(contigs, ref_best, other_best, recovers):
    """a_recovered (hits a truth locus) | ref (best hit >= .999 on the reference, the reference winning ties) | b_other (>= .999 on the truth
    haplotype, no truth locus) | c_unmatched"""
    if recovers:
        return "a_recovered"
    best = None
    for c in contigs:
        for hap, B in (("ref", ref_best), ("other", other_best)):
            h = B.get(c)
            if h and (best is None or h[0] > best[0]):
                best = (h[0], hap)
    if best is None or best[0] < MIN_SCORE:
        return "c_unmatched"
    return "ref" if best[1] == "ref" else "b_other"


def max_matching(edges, n_left):
    """Kuhn's maximum bipartite matching; edges = [(left, right)]"""
    adj = [[] for _ in range(n_left)]
    for a, b in edges:
        adj[a].append(b)
    match = {}

    def try_(u, seen):
        for v in adj[u]:
            if v in seen:
                continue
            seen.add(v)
            if v not in match or try_(match[v], seen):
                match[v] = u
                return True
        return False
    return sum(1 for u in range(n_left) if try_(u, set()))


def evaluate(cands, loci_, fam_net, ref_best, other_best, floor=2):
    """cands: dicts candidate, family, n_transcripts, contigs[list]; loci_: truth loci (see module doc); fam_net: every family that entered the chain.
    R1: each expressed (>= 3), LARGE (>= 20), beyond-delta catalog locus is recovered by a flagged candidate. R2: no expressed within-delta catalog
    locus is recovered by a flagged candidate. R3: families without an expressed truth locus having a flagged non-recovering candidate / those
    families (bar <= 20%). R4 (reported): sensitivity, precision, bipartite matching. kind 'sex' is never counted; kind 'lrpap1' (p12) is
    reported on the p12 line, not in R1/R2."""
    expressed = [L for L in loci_ if L["kind"] != "sex" and L["n"] >= 3]
    flagged = []
    for c in cands:
        if c["n_transcripts"] < floor:
            continue
        rec = recovered_loci(c["contigs"], other_best, loci_)
        flagged.append(dict(c, recovers=sorted(rec), cls=candidate_class(c["contigs"], ref_best, other_best, rec)))
    got = {l for c in flagged for l in c["recovers"]}
    beyond = 1 - C.DELTA
    t1 = [L["locus"] for L in expressed if L["kind"] == "catalog" and L["n"] >= 20 and L["ident"] < beyond]
    r1 = dict(targets=t1, recovered=[l for l in t1 if l in got], verdict=("PASS" if t1 and all(l in got for l in t1) else "NOT TESTABLE" if not t1 else "FAIL"))
    t2 = [L["locus"] for L in expressed if L["kind"] == "catalog" and L["ident"] >= beyond]
    flagged_within = [l for l in t2 if l in got]
    r2 = dict(targets=t2, flagged=flagged_within, verdict="FAIL" if flagged_within else "PASS")
    fam_with = {L["family"] for L in expressed}
    den = [f for f in fam_net if f not in fam_with]
    bad = sorted({c["family"] for c in flagged if c["cls"] != "a_recovered" and c["family"] in den})
    frac = len(bad) / len(den) if den else None
    r3 = dict(denominator=len(den), false_families=bad, fraction=frac, verdict="NOT TESTABLE" if frac is None else ("PASS" if frac <= 0.20 else "FAIL"))
    idx = {L["locus"]: i for i, L in enumerate(expressed)}
    edges = [(i, idx[l]) for i, c in enumerate(flagged) for l in c["recovers"] if l in idx]
    m = max_matching(edges, len(flagged))
    r4 = dict(expressed=len(expressed), flagged=len(flagged), matched=m,
              sensitivity=(m / len(expressed)) if expressed else None, precision=(m / len(flagged)) if flagged else None)
    return dict(R1=r1, R2=r2, R3=r3, R4=r4, rows=flagged)


def type_hits(mat_best, pat_best, targets, min_score=MIN_SCORE):
    """targets: {name: (hap, acc, s, e)} -> {name: [sequences whose best hit on that haplotype scores >= min_score and overlaps the target]}"""
    out = {}
    for name, (hap, acc, s, e) in targets.items():
        B = mat_best if hap == "mat" else pat_best
        out[name] = sorted(q for q, h in B.items() if h[0] >= min_score and h[1] == acc and h[2] < e and s < h[3])
    return out


def load_loci():
    fate = json.load(open(f"{C.WR}/fate/fate.json"))["loci"]
    out = []
    for r in csv.DictReader(open(f"{C.WR}/truth/loci.tsv"), delimiter="\t"):
        f = fate.get(r["locus"], {})
        out.append(dict(locus=r["locus"], kind=r["kind"], family=r["family"], chrom=r["chrom"], start=int(r["start"]), end=int(r["end"]),
                        n=f.get("n", 0), ident=f.get("paralog_identity") if f.get("paralog_identity") is not None else 0.0))
    return out


def p12_targets():
    rows = {r["name"]: r for r in csv.DictReader(open(f"{C.WR}/truth/lrpap1_loci.tsv"), delimiter="\t")}
    t = {}
    for key, name in (("p12", "LOC134756368"), ("p14", "LOC115932954")):
        r = rows[name]
        if r["pat_acc"]:
            t[f"{key}@pat"] = ("pat", r["pat_acc"], int(r["pat_s"]), int(r["pat_e"]))
        if r["mat_acc"]:
            t[f"{key}@mat"] = ("mat", r["mat_acc"], int(r["mat_s"]), int(r["mat_e"]))
    return t


def main(arm):
    al = C.alias()
    d = f"{C.WR}/{'isoc' if arm == 'isocon' else 'inhouse'}"
    cands = [dict(candidate=r["candidate"], family=r["family"], n_transcripts=int(r["n_transcripts"]), contigs=r["contigs"].split(","))
             for r in csv.DictReader(open(f"{d}/cands.tsv"), delimiter="\t")]
    best = lambda p: C.best_hits(p, al) if os.path.exists(p) and os.path.getsize(p) else {}
    ref_best, other_best = best(f"{d}/cands.ref.paf"), best(f"{d}/cands.other.paf")
    fam_net = [p["fam"] for p in json.load(open(f"{d}/panel.json"))]
    res = evaluate(cands, load_loci(), fam_net, ref_best, other_best)
    res["candidates"] = len(cands)
    if arm == "isocon" and C.REF == "mat":      # the p12 line exists only with the mother as the reference (p12 is mother-absent)
        om, op = best(f"{d}/outputs.pri.paf"), best(f"{d}/outputs.other.paf")
        fam = lambda B: {k: v for k, v in B.items() if k.startswith("LRPAP1|")}
        res["p12"] = type_hits(fam(om), fam(op), p12_targets())
    os.makedirs(f"{C.WR}/score", exist_ok=True)
    json.dump(res, open(f"{C.WR}/score/{arm}.json", "w"), indent=1)
    for k in ("R1", "R2", "R3", "R4"):
        print(k, json.dumps(res[k]))
    print("p12 line (sequences whose best hit is >= .999 on each target):", json.dumps(res.get("p12", "n/a for this arm / reference")))
    print("flagged candidates by class:", {c: sum(1 for r in res["rows"] if r["cls"] == c) for c in ("a_recovered", "ref", "b_other", "c_unmatched")})


if __name__ == "__main__":
    main(sys.argv[1])
```

- [ ] **Step 4: Run the tests to verify they pass**

Run: `cd /mnt/linuxdisk/home/juanfraitu/rustle_m2_soto && python3 -B -m unittest bench/o3_maternal/test_score.py`
Expected: `Ran 12 tests ... OK`.

- [ ] **Step 5: Score both arms (once per reference, here `mat`) and record the verdicts before anything else is run**

```bash
# run environment block first (O3_REF=mat)
python3 bench/o3_maternal/score.py isocon
python3 bench/o3_maternal/score.py inhouse
{ cat $WR/fate/fate.tsv; python3 bench/o3_maternal/score.py isocon; python3 bench/o3_maternal/score.py inhouse; } > $WR/VERDICTS.txt
```
Expected: R1-R4 verdict lines for each arm. Write nothing about winners; report each bar per arm. The p12 line prints for IsoCon only (it needs all outputs); for the in-house arm state "n/a" in the write-up. **`W/mat/VERDICTS.txt` must exist before Task 10 starts** (Amendment 1).

- [ ] **Step 6: Commit**

```bash
# run environment block first
git add bench/o3_maternal/score.py bench/o3_maternal/test_score.py
commit "o3_maternal: R1-R4 scorer and the p12 line"
```

---

### Task 10: The paternal-reference pass (runs only, no new code)

**Files:** none created. Outputs under `W/pat/`.

**Interfaces:**
- Consumes: Tasks 1-9 code; `W/mat/VERDICTS.txt` (must exist and be non-empty).
- Produces: `W/pat/{truth,fate,isoc,inhouse,score}/...` with the same names as `W/mat/`; `W/pat/VERDICTS.txt`.

- [ ] **Step 1: Gate — the mat verdicts are on record**

```bash
test -s /mnt/linuxdisk/tmp/o3_mat/mat/VERDICTS.txt && echo "mat verdicts recorded" || echo "STOP: finish and score the mat run first"
```

- [ ] **Step 2: Repeat, with `export O3_REF=pat` before the run environment block, exactly the run steps of**
  - Task 4 steps 5-6 (`truth.py loci`, `labels`). Expected: `reference pat, truth mat`; catalog loci = the `hap = mat` rows of `bonly.tsv` (92); LRPAP1 loci absent from `pat`: the PAFs decide (expected none); `LRPAP1_c07` (chrY) is **present** in `pat`.
  - Task 5 steps 5-6 (`fate.py fasta`, the nearest-paralog `minimap2`, `fate.py run`); skip `fate.py unm`.
  - Task 6 steps 5-6 (`chain_inputs.py panel`, `fasta` — 24 sequences for `pat`, `panel_to_copies.py`).
  - Task 7 steps 2-5 (net, IsoCon batches until `ALL_DONE`, outputs and flag step against `$IDXREF` = the `pat` index, contigs, adapter, the three hit-table mappings).
  - Task 8 steps 3-4 (four in-house batches on the `pat` panel; adapter; mappings; delete `W/pat.idx.fa`).
  - Task 9 step 5 (`score.py isocon`, `inhouse`; write `W/pat/VERDICTS.txt`).

- [ ] **Step 3: Check the pass is complete and symmetric**

```bash
export O3_REF=pat
# run environment block
ls $WR/fate/fate.json $WR/score/isocon.json $WR/score/inhouse.json $WR/isoc/outputs.other.paf $WR/inhouse/cands_all.tsv
grep -c . $WR/truth/loci.tsv
```
Expected: all files exist. A missing file means a step above was skipped: do it, do not continue.

---

### Task 11: Side-by-side — paired reads, divergence pile, the chain's view of the same copy

**Files:**
- Create: `bench/o3_maternal/side.py`, `bench/o3_maternal/test_side.py`, `bench/o3_maternal/compare.py`, `bench/o3_maternal/test_compare.py`
- Outputs: `W/mat/fate/side.json`, `W/pat/fate/side.json`, `W/compare_mat.json`, `W/compare_pat.json`

**Interfaces:**
- Consumes: `W/<REF>/fate/fate.json`, `W/<REF>/truth/loci.tsv`, both `W/map/reads.<hap>.all.bam`; for `compare.py`: `W/<A>/{truth/loci.tsv,isoc/outputs.other.paf,isoc/contigs.tsv,score/*.json}` and `W/<B>/{isoc/outputs.pri.paf,isoc/contigs.tsv,inhouse/cands_all.tsv}`.
- Produces:
  - `side.add_alignment(cov, mis, ref_start, cigar, lo, hi)`, `side.pair_rows(reads, recs_ref, recs_other, fate_of) -> [[read, fate, de_ref, mapq_ref, de_other, mapq_other]]`, `side.track(bam_path, names, target, al, nbins=60) -> {target, cov, mis}`
  - `compare.outputs_on(best, loci_, min_score)`, `compare.new_outputs(contig_rows)`, `compare.recovering(score_rows, locus)`, `compare.inhouse_near(rows, loci_, al)`
  - `W/<REF>/fate/side.json` = `{locus: {n, kind, pairs, de_ref_median, de_other_median, mapq_ref_median, mapq_other_median, [ref_track, other_track]}}`; `W/compare_<A>.json` = `{a, b, loci: {locus: {kind, isocon: {run_a: {outputs, flagged_new, recovering}, run_b: {outputs, flagged_new}}, inhouse: {run_a: {recovering}, run_b: {near}}}}}`.

- [ ] **Step 1: Write the failing tests**

`bench/o3_maternal/test_side.py`:

```python
#!/usr/bin/env python3
"""Run: python3 -B -m unittest bench/o3_maternal/test_side.py"""
import os
import sys
import tempfile
import unittest

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import side as S  # noqa: E402
import testutil  # noqa: E402
from common import Rec  # noqa: E402


class Track(unittest.TestCase):
    def run_one(self, cigar, start=0, lo=0, hi=100, n=2):
        cov, mis = [0] * n, [0] * n
        S.add_alignment(cov, mis, start, cigar, lo, hi)
        return cov, mis

    def test_mismatch_lands_in_its_bin(self):
        self.assertEqual(self.run_one([(7, 10), (8, 1), (7, 89)]), ([50, 50], [1, 0]))

    def test_insertion_adds_a_mismatch_without_coverage(self):
        self.assertEqual(self.run_one([(7, 60), (1, 3), (7, 40)]), ([50, 50], [0, 1]))

    def test_deletion_is_covered_and_mismatched(self):
        self.assertEqual(self.run_one([(7, 40), (2, 20), (7, 40)]), ([50, 50], [10, 10]))

    def test_intron_covers_nothing(self):
        self.assertEqual(self.run_one([(7, 20), (3, 500), (7, 20)], hi=1000), ([20, 20], [0, 0]))

    def test_outside_the_interval_is_ignored(self):
        self.assertEqual(self.run_one([(7, 100)], start=500), ([0, 0], [0, 0]))


class Pairs(unittest.TestCase):
    def rec(self, de, mapq, primary=True):
        return Rec(primary, mapq, 100, 1.0, "A", 0, 10, de)

    def test_pairs_need_a_primary_on_both(self):
        ref = {"a": [self.rec(0.05, 60)], "b": [self.rec(0.05, 60)], "c": []}
        oth = {"a": [self.rec(0.001, 39)], "b": [], "c": [self.rec(0.0, 60)]}
        rows = S.pair_rows(["a", "b", "c"], ref, oth, {"a": "ABSORBED_OTHER", "b": "TIED", "c": "UNMAPPED"})
        self.assertEqual(rows, [["a", "ABSORBED_OTHER", 0.05, 60, 0.001, 39]])


class BamTrack(unittest.TestCase):
    def test_track_counts_only_the_named_primaries(self):
        with tempfile.TemporaryDirectory() as d:
            p = os.path.join(d, "t.bam")
            testutil.write_bam(p, {"chr1_mat_hsa1": 5000}, [
                dict(name="a", ref="chr1_mat_hsa1", start=100, cigar="50=1X49="),
                dict(name="b", ref="chr1_mat_hsa1", start=100, cigar="100="),                    # not named
                dict(name="a", flag=256, ref="chr1_mat_hsa1", start=100, cigar="100=", mapq=0),   # secondary: ignored
            ])
            t = S.track(p, {"a"}, ("CM1", 100, 200), {("mat", "1"): "CM1"}, nbins=2)
            self.assertEqual(t["cov"], [50, 50])
            self.assertEqual(t["mis"], [0, 1])
            self.assertEqual(t["target"], ["CM1", 100, 200])


if __name__ == "__main__":
    unittest.main()
```

`bench/o3_maternal/test_compare.py`:

```python
#!/usr/bin/env python3
"""Run: python3 -B -m unittest bench/o3_maternal/test_compare.py"""
import os
import sys
import unittest

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import compare as K  # noqa: E402

LOCI = [dict(locus="A", chrom="CM1", start=100, end=200), dict(locus="B", chrom="CM1", start=500, end=600)]


class Outputs(unittest.TestCase):
    def test_outputs_on_the_locus_at_the_floor(self):
        best = {"o1": (1.0, "CM1", 120, 180, 1.0), "o2": (0.998, "CM1", 120, 180, 0.998), "o3": (1.0, "CM2", 120, 180, 1.0),
                "o4": (0.9995, "CM1", 550, 590, 0.9995)}
        self.assertEqual(K.outputs_on(best, LOCI), {"A": ["o1"], "B": ["o4"]})

    def test_new_outputs_are_the_unlinked_ones(self):
        rows = [dict(output="F|x one", linked="0"), dict(output="F|y", linked="1")]
        self.assertEqual(K.new_outputs(rows), {"F|x"})


class Inhouse(unittest.TestCase):
    def test_recovering(self):
        rows = [dict(candidate="k1", n_transcripts=3, recovers=["A"]), dict(candidate="k2", n_transcripts=2, recovers=["B"])]
        self.assertEqual(K.recovering(rows, "A"), [dict(candidate="k1", n_transcripts=3)])

    def test_near_uses_accessions_and_overlap(self):
        al = {("pat", "5"): "CM1"}
        rows = [dict(candidate="c1", n_clusters="4", d="0.0003", flagged="0", nearest_locus="chr5_pat_hsa5:150-400"),
                dict(candidate="c2", n_clusters="2", d="0.2", flagged="1", nearest_locus="chr5_pat_hsa5:1000-2000")]
        out = K.inhouse_near(rows, LOCI, al)
        self.assertEqual(out["A"], [dict(candidate="c1", n_clusters=4, d=0.0003, flagged=0)])
        self.assertEqual(out["B"], [])


if __name__ == "__main__":
    unittest.main()
```

- [ ] **Step 2: Run to verify failure**

Run: `cd /mnt/linuxdisk/home/juanfraitu/rustle_m2_soto && python3 -B -m unittest bench/o3_maternal/test_side.py bench/o3_maternal/test_compare.py`
Expected: `ModuleNotFoundError: No module named 'side'` (and `compare`).

- [ ] **Step 3: Write `bench/o3_maternal/side.py` and `bench/o3_maternal/compare.py`**

`bench/o3_maternal/side.py`:

```python
#!/usr/bin/env python3
"""Side-by-side read view of one run (docs/PREREG_o3_maternal_reference_2026-10-08.md, Amendment 1): the reads of each absent locus on the
reference haplotype (where the copy is missing) and on the other haplotype (where it is present).

    O3_REF=mat side.py     # WR/fate/side.json (needs WR/fate/fate.json, WR/truth/loci.tsv and both W/map/reads.<hap>.all.bam)

Per locus with >= 3 reads: `pairs` = [read, fate on the reference, de ref, MAPQ ref, de other, MAPQ other]; for LARGE loci (>= 20 reads) also
`ref_track` (the nearest reference paralog, where the reads land) and `other_track` (the locus itself): per-bin coverage and mismatch counts of
the primary alignments, the interval cut in NBINS equal bins (positions are fractions of the interval, not shared coordinates)."""
import csv
import json
import os
import statistics
import sys

import pysam

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import common as C  # noqa: E402

NBINS = 60


def add_alignment(cov, mis, ref_start, cigar, lo, hi):
    """Add one alignment to the per-bin coverage / mismatch counts of the interval [lo, hi) cut in len(cov) equal bins.
    cigar: pysam (op, length) tuples with --eqx operators: 7 '=', 8 'X', 0 'M' (covered, mismatch unknown), 1 'I', 2 'D', 3 'N' (intron: nothing).
    A deleted base counts as covered and mismatched; an insertion adds one mismatch at the current reference base and no coverage."""
    n = len(cov)

    def b(p):
        return min(n - 1, (p - lo) * n // (hi - lo))
    pos = ref_start
    for op, ln in cigar:
        if op in (0, 7, 8, 2):
            for p in range(max(pos, lo), min(pos + ln, hi)):
                k = b(p)
                cov[k] += 1
                if op in (8, 2):
                    mis[k] += 1
            pos += ln
        elif op == 3:
            pos += ln
        elif op == 1 and lo <= pos < hi:
            mis[b(pos)] += 1


def pair_rows(reads, recs_ref, recs_other, fate_of):
    """[[read, fate, de ref, MAPQ ref, de other, MAPQ other]] for the reads with a primary record on both haplotypes"""
    out = []
    for n in reads:
        r, o = recs_ref.get(n), recs_other.get(n)
        if r and o and r[0].primary and o[0].primary:
            out.append([n, fate_of[n], r[0].de, r[0].mapq, o[0].de, o[0].mapq])
    return out


def track(bam_path, names, target, al, nbins=NBINS):
    """per-bin coverage and mismatches of the primary alignments of `names` overlapping target = (accession, lo, hi)"""
    cov, mis = [0] * nbins, [0] * nbins
    with pysam.AlignmentFile(bam_path) as bam:
        idx = {C.accession(s["SN"], al): s["SN"] for s in bam.header.to_dict()["SQ"]}
        for rd in bam.fetch(idx[target[0]], target[1], target[2]):
            if rd.is_unmapped or rd.is_secondary or rd.is_supplementary or rd.query_name not in names:
                continue
            add_alignment(cov, mis, rd.reference_start, rd.cigartuples, target[1], target[2])
    return dict(target=list(target), cov=cov, mis=mis)


def median(xs):
    return statistics.median(xs) if xs else None


def main():
    al = C.alias()
    fate = json.load(open(f"{C.WR}/fate/fate.json"))["loci"]
    loci_ = {r["locus"]: r for r in csv.DictReader(open(f"{C.WR}/truth/loci.tsv"), delimiter="\t")}
    bam_ref, bam_oth = f"{C.W}/map/reads.{C.REF}.all.bam", f"{C.W}/map/reads.{C.OTHER}.all.bam"
    recs_ref, recs_oth = C.read_records(bam_ref, al), C.read_records(bam_oth, al)
    out = {}
    for k, v in fate.items():
        if v["n"] < 3:
            continue
        fate_of = {r[0]: r[1] for r in v["reads"]}
        pairs = pair_rows(sorted(fate_of), recs_ref, recs_oth, fate_of)
        row = dict(n=v["n"], kind=v["kind"], pairs=pairs, de_ref_median=median([p[2] for p in pairs if p[2] is not None]),
                   de_other_median=median([p[4] for p in pairs if p[4] is not None]),
                   mapq_ref_median=median([p[3] for p in pairs]), mapq_other_median=median([p[5] for p in pairs]))
        if v["n"] >= 20 and v["paralog"]:
            L = loci_[k]
            row["ref_track"] = track(bam_ref, set(fate_of), tuple(v["paralog"]), al)
            row["other_track"] = track(bam_oth, set(fate_of), (L["chrom"], int(L["start"]), int(L["end"])), al)
        out[k] = row
    json.dump(out, open(f"{C.WR}/fate/side.json", "w"))
    print(f"reference {C.REF}: {len(out)} loci with >= 3 reads; tracks for {sum(1 for r in out.values() if 'ref_track' in r)}")
    for k, r in out.items():
        print(f"{k}\tn={r['n']}\tde {C.REF} {r['de_ref_median']}\tde {C.OTHER} {r['de_other_median']}\tMAPQ {r['mapq_ref_median']} / {r['mapq_other_median']}")


if __name__ == "__main__":
    main()
```

`bench/o3_maternal/compare.py`:

```python
#!/usr/bin/env python3
"""Chain side-by-side of the two runs (docs/PREREG_o3_maternal_reference_2026-10-08.md, Amendment 1).

For the loci absent from the reference of run A (A = mat: copies only the father has; A = pat: the reverse) show the SAME copies as seen by
run A (the reference lacks them: are their transcripts flagged as new, are they recovered?) and by run B (the reference has them: are their
transcripts simply found in it?). Locus coordinates are on the truth haplotype of run A, which is the reference of run B.

    compare.py mat | pat      # writes W/compare_<A>.json; needs W/<A>/{truth,isoc,inhouse,score} and W/<B>/{isoc,inhouse}
"""
import csv
import json
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import common as C  # noqa: E402

MIN_SCORE = 0.999


def outputs_on(best, loci_, min_score=MIN_SCORE):
    """{locus: [output names]}: outputs whose best hit (identity x coverage >= min_score) overlaps the locus; best = common.best_hits() on the
    haplotype the locus lies on"""
    out = {L["locus"]: [] for L in loci_}
    for q, h in sorted(best.items()):
        if h[0] < min_score:
            continue
        for L in loci_:
            if h[1] == L["chrom"] and h[2] < L["end"] and L["start"] < h[3]:
                out[L["locus"]].append(q)
    return out


def new_outputs(contig_rows):
    """names of the IsoCon outputs the link step kept as new copies (contigs.tsv rows with linked == '0')"""
    return {r["output"].split()[0] for r in contig_rows if r["linked"] == "0"}


def recovering(score_rows, locus):
    """[(candidate, n_transcripts)] of the flagged candidates of a score JSON that recover the locus"""
    return [dict(candidate=r["candidate"], n_transcripts=r["n_transcripts"]) for r in score_rows if locus in r["recovers"]]


def inhouse_near(rows, loci_, al):
    """{locus: [dict(candidate, n_clusters, d, flagged)]}: stage candidates (cands_all.tsv rows) whose nearest locus overlaps the locus"""
    out = {L["locus"]: [] for L in loci_}
    for r in rows:
        chrom, rng = r["nearest_locus"].rsplit(":", 1)
        s, e = (int(x) for x in rng.split("-"))
        acc = C.accession(chrom, al)
        for L in loci_:
            if acc == L["chrom"] and s < L["end"] and L["start"] < e:
                out[L["locus"]].append(dict(candidate=r["candidate"], n_clusters=int(r["n_clusters"]), d=float(r["d"]), flagged=int(r["flagged"])))
    return out


def tsv(path):
    return list(csv.DictReader(open(path), delimiter="\t"))


def main(a):
    b = "pat" if a == "mat" else "mat"
    al = C.alias()
    WA, WB = f"{C.W}/{a}", f"{C.W}/{b}"
    loci_ = tsv(f"{WA}/truth/loci.tsv")
    for L in loci_:
        L["start"], L["end"] = int(L["start"]), int(L["end"])
    best = lambda p: C.best_hits(p, al) if os.path.exists(p) and os.path.getsize(p) else {}
    on_a = outputs_on(best(f"{WA}/isoc/outputs.other.paf"), loci_)      # run A: outputs against the truth haplotype
    on_b = outputs_on(best(f"{WB}/isoc/outputs.pri.paf"), loci_)        # run B: outputs against ITS reference = the same haplotype
    new_a, new_b = new_outputs(tsv(f"{WA}/isoc/contigs.tsv")), new_outputs(tsv(f"{WB}/isoc/contigs.tsv"))
    sc = {arm: json.load(open(f"{WA}/score/{arm}.json")) for arm in ("isocon", "inhouse") if os.path.exists(f"{WA}/score/{arm}.json")}
    near = inhouse_near(tsv(f"{WB}/inhouse/cands_all.tsv"), loci_, al) if os.path.exists(f"{WB}/inhouse/cands_all.tsv") else {}
    out = {}
    for L in loci_:
        k = L["locus"]
        out[k] = dict(
            kind=L["kind"],
            isocon=dict(run_a=dict(outputs=len(on_a[k]), flagged_new=sum(1 for o in on_a[k] if o in new_a),
                                   recovering=recovering(sc.get("isocon", {}).get("rows", []), k)),
                        run_b=dict(outputs=len(on_b[k]), flagged_new=sum(1 for o in on_b[k] if o in new_b))),
            inhouse=dict(run_a=dict(recovering=recovering(sc.get("inhouse", {}).get("rows", []), k)),
                         run_b=dict(near=near.get(k, []))))
    json.dump(dict(a=a, b=b, loci=out), open(f"{C.W}/compare_{a}.json", "w"), indent=1)
    print(f"run A = {a} (reference lacks the copy), run B = {b} (reference has it)")
    print("locus\tkind\tIsoCon outputs on the locus: A (flagged new) | B (flagged new)\trecovered by IsoCon / in-house candidate (A)\tin-house near the locus (B)")
    for k, v in out.items():
        i, h = v["isocon"], v["inhouse"]
        print(f"{k}\t{v['kind']}\t{i['run_a']['outputs']} ({i['run_a']['flagged_new']}) | {i['run_b']['outputs']} ({i['run_b']['flagged_new']})\t"
              f"{len(i['run_a']['recovering'])} / {len(h['run_a']['recovering'])}\t{len(h['run_b']['near'])}")


if __name__ == "__main__":
    main(sys.argv[1])
```

- [ ] **Step 4: Run the tests to verify they pass**

Run: `cd /mnt/linuxdisk/home/juanfraitu/rustle_m2_soto && python3 -B -m unittest bench/o3_maternal/test_side.py bench/o3_maternal/test_compare.py`
Expected: `Ran 11 tests ... OK`.

- [ ] **Step 5: Run both directions**

```bash
# run environment block first
for r in mat pat; do O3_REF=$r tools/rlock.sh light python3 bench/o3_maternal/side.py; done
python3 bench/o3_maternal/compare.py mat
python3 bench/o3_maternal/compare.py pat
```
Expected: a per-locus table for each reference (paired-read medians; tracks for LARGE loci) and the chain table. Compare with prereg Amendment 1 predictions S1-S2 and record each as held or failed in the write-up. If `side.py` takes more than 3 minutes (reading both merged BAMs), run it under `tools/rlock.sh heavy`.

- [ ] **Step 6: Commit**

```bash
# run environment block first
git add bench/o3_maternal/side.py bench/o3_maternal/test_side.py bench/o3_maternal/compare.py bench/o3_maternal/test_compare.py
commit "o3_maternal: side-by-side paired reads, divergence pile and the chain's view of the same copy"
```

---

### Task 12: Data file and artifact (both directions, side by side)

**Files:**
- Create: `bench/o3_maternal/artifact_data.py`, `bench/o3_maternal/template.html`
- Outputs: `W/artifact/{data.json,index.html}`; a published private artifact

**Interfaces:**
- Consumes: both runs' `fate.json`, `side.json`, `score/*.json`, `truth/loci.tsv`; `W/compare_{mat,pat}.json`; `W/unm.txt`; `/home/juanfra/winloci_scratch/o3_excise/per_family2.json` (162 families: `fam, n, unaln, conc, is_sister, ratio, dest, mig_de`).
- Produces: `data.json` = `{directions: {mat|pat: {ref, other, title, loci[], shared, arms, p12, compare}}, context[], unm[]}` with `loci[]` items `id, label, kind, n, fates, fractions, verdict, de_median, paralog_identity, reads, pairs, de_other_median, mapq_ref_median, mapq_other_median, ref_track, other_track`.

- [ ] **Step 1: Write `bench/o3_maternal/artifact_data.py`**

```python
#!/usr/bin/env python3
"""W/artifact/data.json and index.html for the artifact. Every number comes from the two runs' fate/side/score/compare files and the excision
table; nothing is typed in. Needs both W/mat and W/pat to be complete (Tasks 4-9 for each reference, then side.py and compare.py)."""
import csv
import json
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import common as C  # noqa: E402

EXC = "/home/juanfra/winloci_scratch/o3_excise/per_family2.json"
LABEL = {"LRPAP1_p12": "LRPAP1 5′ fragment, chr12 23.07 Mb (LOC134756368)", "LRPAP1_c07": "LRPAP1 chrY copy (sex control)"}
TITLE = {"mat": "The mother's genome is the reference: copies only the father has are missing",
         "pat": "The father's genome is the reference: copies only the mother has are missing"}


def load(path, default=None):
    return json.load(open(path)) if os.path.exists(path) else default


def direction(a):
    """one direction: reference haplotype `a`; `b` = the haplotype that has the copies"""
    b = "pat" if a == "mat" else "mat"
    WA = f"{C.W}/{a}"
    fate = load(f"{WA}/fate/fate.json")
    side = load(f"{WA}/fate/side.json", {})
    kinds = {r["locus"]: r for r in csv.DictReader(open(f"{WA}/truth/loci.tsv"), delimiter="\t")}
    loci = []
    for k, v in fate["loci"].items():
        if v["n"] == 0:
            continue
        s = side.get(k, {})
        loci.append(dict(id=k, label=LABEL.get(k, f"{kinds[k]['family']} · {k}"), kind=v["kind"], n=v["n"], fates=v["fates"], fractions=v["fractions"],
                         verdict=v["verdict"], de_median=v["de_median"], paralog_identity=v["paralog_identity"], reads=v["reads"],
                         pairs=s.get("pairs", []), de_other_median=s.get("de_other_median"), mapq_ref_median=s.get("mapq_ref_median"),
                         mapq_other_median=s.get("mapq_other_median"), ref_track=s.get("ref_track"), other_track=s.get("other_track")))
    loci.sort(key=lambda x: (x["kind"] == "sex", -x["n"]))
    arms = {arm: load(f"{WA}/score/{arm}.json") for arm in ("isocon", "inhouse")}
    cmp_ = load(f"{C.W}/compare_{a}.json", {}).get("loci", {})
    return dict(ref=a, other=b, title=TITLE[a], loci=loci, shared=fate["shared"],
                arms={k: {x: v[x] for x in ("R1", "R2", "R3", "R4", "candidates")} for k, v in arms.items() if v}, p12=(arms.get("isocon") or {}).get("p12"),
                compare=cmp_)


def main():
    out = dict(directions={a: direction(a) for a in ("mat", "pat") if os.path.exists(f"{C.W}/{a}/fate/fate.json")},
               context=[dict(fam=r["fam"], unaln=r["unaln"], conc=r["conc"], mig_de=r["mig_de"]) for r in json.load(open(EXC))],
               unm=open(f"{C.W}/unm.txt").read().strip().split("\n") if os.path.exists(f"{C.W}/unm.txt") else [])
    os.makedirs(f"{C.W}/artifact", exist_ok=True)
    json.dump(out, open(f"{C.W}/artifact/data.json", "w"))
    html = open(f"{HERE}/template.html").read().replace("/*DATA*/null", json.dumps(out))
    open(f"{C.W}/artifact/index.html", "w").write(html)
    print(f"data.json: directions {list(out['directions'])}, loci {[len(d['loci']) for d in out['directions'].values()]}, "
          f"{len(out['context'])} excision families; index.html {len(html) // 1024} KB")


if __name__ == "__main__":
    main()
```

- [ ] **Step 2: Load the design skills, then write `bench/o3_maternal/template.html`**

Before writing the file run the Skill tool for `artifact-design` and `dataviz` (both are required by the Artifact tool), and follow them: title of two to four words, color tokens on `:root` with dark-mode redefinitions, external scripts only from the allowed CDNs (none are needed), 16 px gutters, no horizontal page scroll at phone width. The page is self-contained: one `const DATA = /*DATA*/null;` line that `artifact_data.py` replaces; the page `<script>` starts with the rendering core below.

Content contract (every panel reads `DATA`, never hard-coded numbers):

1. **Header + direction toggle** (two buttons, `aria-pressed`, default the mother's genome): title `Missing copies, side by side`; the toggle label is `DATA.directions[k].title`. Tiles for the chosen direction: number of absent loci with reads, total reads, `% absorbed`, `% tied`, `% unmapped` over the LARGE loci (sex control excluded), the two arms' R1 verdicts.
2. **Where the reads go** (`stackedBars`, one bar per locus; the chrY row dashed and labelled "sex control"); a final row `shared copies (control)` from `direction.shared`.
3. **Side by side** (the heart of the page): a locus selector (one button per LARGE locus) and, for the selected locus, two columns **`<reference> reference: copy missing` | `<other> reference: copy present`**:
   - top: the fate bar on the left; on the right the text "the copy is on the reference: the reads sit on it with median de X, MAPQ Y" (from `de_other_median`, `mapq_other_median`);
   - middle: `pairedDe` (one line per read from its de on the reference where the copy is missing to its de where it is present) and, beside it, `trackPair` (coverage and mismatch rate along the paralog the reads land on vs along the locus itself);
   - bottom: `sideChain` (IsoCon transcripts of this copy kept as new and recovered on the left; found in the reference on the right; the in-house stage likewise).
4. **Spectrum** (`spectrum`) for the mother's-genome direction only: x = divergence to the nearest relative, y = fraction of reads unmapped; natural loci filled, excision families hollow ("synthetic deletion, same animal, 2026-08-14"); left half is context only.
5. **p12 read strip**: one small square per p12 read, colored by fate, tooltip with `de`.
6. **Can the chain recover them?** `chainCard` per arm and direction (R1/R2/R3 chips, R4 numbers). The in-house card carries the line "run without the Amendment 15 fix (known false-flag defect)". The p12 line (which of `p12@pat`, `p14@pat`, `p14@mat` have a sequence at >= .999) under the mother's-genome IsoCon card.
7. **Footer**: caveats from the prereg S4/S7/S8 in one paragraph (one animal, fibroblast only, selection through `_pri`, number of expressed loci, `-p 0.1 -N 50` for every identity quoted, DELTA, IsoCon version, the pat run as the held-back substrate), plus the prereg path.

Palette (Okabe-Ito, works in both themes), as CSS custom properties: `--c-near #0072B2`, `--c-oth #56B4E9`, `--c-tie #E69F00`, `--c-unm #D55E00`, neutral `#9AA0AA`; classes `.chip.ok/.bad/.na`, `.card`, `.rd` (a 9 px square), `.lab/.num/.tk/.axis/.cap` for SVG text and lines.

The rendering core (syntax-checked with `node --check`; the design skill's tokens, tooltips and responsive CSS are added around it):

```javascript
const DATA = /*DATA*/null;
const FATE = [["UNMAPPED","unmapped"],["PARTIAL","unmapped"],["TIED","tied"],["ABSORBED_NEAREST","nearest"],["ABSORBED_OTHER","other"]];
const SEG = {unmapped:"var(--c-unm)", tied:"var(--c-tie)", nearest:"var(--c-near)", other:"var(--c-oth)"};
const NS = "http://www.w3.org/2000/svg";
const el = (n, a = {}, t) => { const e = document.createElementNS(NS, n); for (const k in a) e.setAttribute(k, a[k]); if (t != null) e.textContent = t; return e; };
const html = (tag, cls, txt) => { const e = document.createElement(tag); if (cls) e.className = cls; if (txt != null) e.textContent = txt; return e; };
const fateColour = f => f === "TIED" ? SEG.tied : f === "ABSORBED_NEAREST" ? SEG.nearest : f === "ABSORBED_OTHER" ? SEG.other : SEG.unmapped;
const fmt = (x, d = 4) => x == null ? "–" : x.toFixed(d);

function stackedBars(host, rows) {
  const W = 860, rowH = 30, lab = 250, bar = 360, h = rows.length * rowH + 8;
  const svg = el("svg", {viewBox: `0 0 ${W} ${h}`, role: "img", "aria-label": "Fate of the reads of each absent locus"});
  rows.forEach((r, i) => {
    const y = 4 + i * rowH;
    svg.append(el("text", {x: 0, y: y + 17, class: "lab"}, r.label));
    let x = lab;
    const by = {unmapped: 0, tied: 0, nearest: 0, other: 0};
    for (const [f, k] of FATE) by[k] += (r.fates[f] || 0);
    for (const k of ["unmapped", "tied", "nearest", "other"]) {
      const w = r.n ? bar * by[k] / r.n : 0;
      if (w > 0) { const s = el("rect", {x, y, width: w, height: 20, fill: SEG[k]}); s.append(el("title", {}, `${k}: ${by[k]} of ${r.n}`)); svg.append(s); }
      x += w;
    }
    svg.append(el("text", {x: lab + bar + 10, y: y + 15, class: "num"}, `n=${r.n} · id ${fmt(r.paralog_identity)} · de ${fmt(r.de_median)} · ${r.verdict || ""}`));
  });
  host.replaceChildren(svg);
}

function spectrum(host, natural, context) {
  const W = 640, H = 300, m = {l: 48, r: 12, t: 12, b: 36};
  const X = d => m.l + (W - m.l - m.r) * Math.min(d, 0.30) / 0.30, Y = f => m.t + (H - m.t - m.b) * (1 - f);
  const svg = el("svg", {viewBox: `0 0 ${W} ${H}`, role: "img", "aria-label": "Fraction of reads unmapped against divergence to the nearest relative"});
  svg.append(el("line", {x1: m.l, y1: Y(0), x2: W - m.r, y2: Y(0), class: "axis"}), el("line", {x1: m.l, y1: Y(0), x2: m.l, y2: Y(1), class: "axis"}));
  for (const f of [0, .25, .5, .75, 1]) svg.append(el("text", {x: m.l - 6, y: Y(f) + 4, class: "tk", "text-anchor": "end"}, `${Math.round(f * 100)}%`));
  for (const d of [0, .05, .1, .15, .2, .25, .3]) svg.append(el("text", {x: X(d), y: H - 14, class: "tk", "text-anchor": "middle"}, d.toFixed(2)));
  svg.append(el("text", {x: W / 2, y: H - 1, class: "lab", "text-anchor": "middle"}, "divergence to the nearest relative"));
  for (const c of context) if (c.mig_de != null) svg.append(el("circle", {cx: X(c.mig_de), cy: Y(c.unaln), r: 3.5, fill: "none", stroke: "var(--c-oth)", "stroke-width": 1.2}));
  for (const p of natural) if (p.fractions && p.paralog_identity != null) svg.append(el("circle", {cx: X(1 - p.paralog_identity), cy: Y(p.fractions.unmapped), r: 6, fill: "var(--c-near)"}));
  host.replaceChildren(svg);
}

// one line per read: its divergence (de) on the haplotype that lacks the copy (left) and on the one that has it (right), coloured by the fate on the left
function pairedDe(host, pairs, refName, otherName) {
  const W = 420, H = 260, m = {l: 44, r: 44, t: 22, b: 28};
  const vals = pairs.flatMap(p => [p[2] ?? 0, p[4] ?? 0]);
  const maxde = Math.min(0.12, Math.max(0.01, ...vals));
  const Y = d => m.t + (H - m.t - m.b) * (1 - Math.min(d, maxde) / maxde), xl = m.l, xr = W - m.r;
  const svg = el("svg", {viewBox: `0 0 ${W} ${H}`, role: "img", "aria-label": `Divergence of each read on the ${refName} and ${otherName} haplotypes`});
  for (const t of [0, .25, .5, .75, 1]) svg.append(el("text", {x: m.l - 6, y: Y(maxde * t) + 4, class: "tk", "text-anchor": "end"}, (maxde * t).toFixed(3)));
  svg.append(el("text", {x: xl, y: 12, class: "lab", "text-anchor": "middle"}, `${refName}: copy missing`), el("text", {x: xr, y: 12, class: "lab", "text-anchor": "middle"}, `${otherName}: copy present`));
  for (const p of pairs) if (p[2] != null && p[4] != null)
    svg.append(el("line", {x1: xl, y1: Y(p[2]), x2: xr, y2: Y(p[4]), stroke: fateColour(p[1]), "stroke-opacity": .35, "stroke-width": 1}));
  host.replaceChildren(svg);
}

// coverage (up) and mismatch rate per base (down, 10% = full height) in 60 bins, for the locus the reads land on (reference) and the locus itself (other)
function trackPair(host, refTrack, otherTrack, refName, otherName) {
  const W = 420, H = 150, m = {l: 8, r: 8, t: 16, b: 6};
  const mk = (tr, title) => {
    const svg = el("svg", {viewBox: `0 0 ${W} ${H}`, role: "img", "aria-label": `${title}: coverage and mismatch rate along the locus`});
    const n = tr.cov.length, bw = (W - m.l - m.r) / n, cmax = Math.max(1, ...tr.cov), half = 0.5 * (H - m.t - m.b), mid = m.t + half;
    svg.append(el("text", {x: m.l, y: 11, class: "lab"}, title));
    tr.cov.forEach((c, i) => {
      const x = m.l + i * bw, h1 = half * c / cmax, h2 = half * Math.min(1, (c ? tr.mis[i] / c : 0) / 0.1);
      svg.append(el("rect", {x, y: mid - h1, width: Math.max(1, bw - .5), height: h1, fill: "var(--c-near)", "fill-opacity": .5}));
      svg.append(el("rect", {x, y: mid + 2, width: Math.max(1, bw - .5), height: h2, fill: "var(--c-unm)", "fill-opacity": .85}));
    });
    return svg;
  };
  host.replaceChildren(mk(refTrack, `${refName}: where the reads land`), mk(otherTrack, `${otherName}: the copy itself`));
}

function chip(txt, cls) { return html("span", "chip " + cls, txt); }
function chainCard(host, title, arm, note) {
  const d = html("div", "card"); d.append(html("h3", null, title));
  if (!arm) { d.append(html("p", "cap", "not run")); host.append(d); return; }
  const row = html("p");
  for (const k of ["R1", "R2", "R3"]) row.append(chip(`${k} ${arm[k].verdict}`, arm[k].verdict === "PASS" ? "ok" : arm[k].verdict === "FAIL" ? "bad" : "na"), document.createTextNode(" "));
  const r4 = arm.R4;
  d.append(row, html("p", null, `flagged ${r4.flagged} · matched ${r4.matched} of ${r4.expressed} expressed loci · sensitivity ${fmt(r4.sensitivity, 2)} · precision ${fmt(r4.precision, 2)}`));
  if (note) d.append(html("p", "cap", note));
  host.append(d);
}

// the same copy seen by the two runs: IsoCon outputs on it, flagged as new (reference lacks it) or found in the reference; in-house candidates
function sideChain(host, cmp, refName, otherName) {
  const a = html("div", "card"), b = html("div", "card");
  a.append(html("h3", null, `${refName} reference: copy missing`));
  b.append(html("h3", null, `${otherName} reference: copy present`));
  if (cmp) {
    const i = cmp.isocon, h = cmp.inhouse;
    a.append(html("p", null, `IsoCon: ${i.run_a.outputs} transcripts of this copy, ${i.run_a.flagged_new} kept as a new copy`),
             html("p", null, i.run_a.recovering.length ? `recovered by ${i.run_a.recovering.map(r => `${r.candidate} (${r.n_transcripts} transcripts)`).join(", ")}` : "not recovered by a flagged IsoCon candidate"),
             html("p", null, h.run_a.recovering.length ? `in-house stage: recovered by ${h.run_a.recovering.map(r => r.candidate).join(", ")}` : "in-house stage: not recovered"));
    b.append(html("p", null, `IsoCon: ${i.run_b.outputs} transcripts of this copy match the reference, ${i.run_b.flagged_new} kept as new`),
             html("p", null, h.run_b.near.length ? `in-house stage: ${h.run_b.near.length} candidates sit on it, closest d ${fmt(Math.min(...h.run_b.near.map(x => x.d)))}, ${h.run_b.near.filter(x => x.flagged).length} flagged` : "in-house stage: no candidate on it"));
  }
  host.replaceChildren(a, b);
}
```

- [ ] **Step 3: Build the data and render-test locally**

```bash
# run environment block first
python3 bench/o3_maternal/artifact_data.py
```
Render `W/artifact/index.html` with Windows Chrome headless from WSL (memory `reference_wsl_chrome_headless_render`) at 1200 px and 390 px wide, light and dark, both directions; look at the PNGs. Checks: no horizontal page scroll at 390 px, every bar's segments sum to the bar width, the chrY row is visibly a control, the left and right columns of the side-by-side card have the same axis range for `pairedDe`, no text clipped, tooltips present, verdict chips match `fate.tsv` / `score/*.json`.

- [ ] **Step 4: Publish (private) and commit**

Publish `W/artifact/index.html` with the Artifact tool (new artifact, `icon` = `chart`, one-sentence description). Give the user the link; do not share beyond the private link.

```bash
# run environment block first
git add bench/o3_maternal/artifact_data.py bench/o3_maternal/template.html
commit "o3_maternal: artifact data builder and page template (both directions, side by side)"
```

---

### Task 13: Write-up, register drafts, memory, cleanup, final review

**Files:**
- Create: `docs/O3_MATERNAL_REFERENCE_2026-10-08.md`, `docs/REGISTER_DRAFTS_o3_maternal.md`
- Modify: memory `project_o3_maternal_reference_2026-10-08.md` (+ its `MEMORY.md` line)

- [ ] **Step 1: Write `docs/O3_MATERNAL_REFERENCE_2026-10-08.md`**

Sections, in this order, each number copied from `fate.tsv` / `score/*.json` / `compare_*.json` (never retyped from memory): (1) result in one paragraph per question and per direction; (2) fate table per locus with the prediction column from prereg S5 and PASS/FAIL against it; (3) the unmapped-set line (`W/unm.txt`) and the selection caveat; (4) the side-by-side findings against Amendment 1 predictions S1-S3 (held / failed, per locus); (5) the arms table (R1-R4 per arm and direction) and the p12 line; (6) what failed or was not as predicted, stated plainly; (7) provenance: minimap2 2.30 vs 2.31, `-p 0.1 -N 50`, IsoCon 0.3.3, `sha256sum` of the `o3_candidates` binary, the Task 4 label cross-check against Amendment 10, the mat-before-pat order (`W/mat/VERDICTS.txt` timestamp); (8) not covered (prereg S8 and A1.1). No pooled percentage; no winner between arms.

- [ ] **Step 2: Draft the register rows** in `docs/REGISTER_DRAFTS_o3_maternal.md` (one row per decided clause: Q1 per LARGE locus and direction, R1-R3 per arm and direction, S1-S3, the p12 line), following the format of `docs/NEGATIVE_RESULTS_REGISTER.md` (read its last 10 rows first). Do not edit the register itself.

- [ ] **Step 3: Cleanup**

```bash
# run environment block first
rm -f $W/mat.idx.fa $W/mat.idx.fa.fai $W/pat.idx.fa $W/pat.idx.fa.fai
rm -rf $W/tmp
du -sh $W; df -h /mnt/linuxdisk | tail -1
```
Keep `W/{mat,pat}/{truth,fate,score}`, `W/artifact`, `W/compare_*.json` and the per-arm candidate tables; the merged BAMs may be deleted after the write-up if space is needed (ask first).

- [ ] **Step 4: Final whole-branch review (fresh reviewer)**

Dispatch ONE fresh reviewer agent on the most capable model (`feedback_final_review_catches_what_task_review_misses`), giving it a frozen read-only copy: `cp -r bench/o3_maternal /mnt/linuxdisk/tmp/o3_mat_review/code` plus copies of `W/{mat,pat}/{truth,fate,score}`, `W/compare_*.json`, the prereg and the write-up. Brief: "check every number in the write-up against the tables, every clause of the prereg (body and Amendment 1) against the code, and the seven Review Focus items; report findings, change nothing." Fix what it confirms, re-run the affected unit tests.

- [ ] **Step 5: Update memory and commit**

Update `project_o3_maternal_reference_2026-10-08.md` with the results (what held, what failed, artifact link, paths) and refresh its `MEMORY.md` line. Then:

```bash
# run environment block first
git add docs/O3_MATERNAL_REFERENCE_2026-10-08.md docs/REGISTER_DRAFTS_o3_maternal.md
commit "o3_maternal: results write-up and register drafts"
```
Do not push; tell the user the branch has local commits.

---

## Self-Review

1. **Spec coverage.** S1 Q1 -> Tasks 4-5, Q2 -> Tasks 6-9. S2 substrate -> Global Constraints. S3 truth (catalog, LRPAP1, chrY control) -> Task 4; the reverse direction (S3 item 4, superseded by Amendment 1) -> Task 10. S4 reads -> Tasks 2-3. S5 labels + fate + bar + context panel -> Tasks 4, 5, 12. S6 chain, both arms, scoring, R1-R4, p12 line -> Tasks 6-9. S7 outcomes -> Task 13. S8 not covered -> Task 13. S9 compute plan -> Global Constraints + per-task commands. Amendment 1: A1.1 symmetric runs + order -> Global Constraints, Task 9 step 5, Task 10; A1.2 (a) paired reads, (b) divergence pile -> Task 11 `side.py`; (c) the chain's view -> Task 11 `compare.py`; A1.3 predictions -> Task 11 step 5 and Task 13; A1.4 visual contract -> Task 12; A1.5 cost -> Task 10.
2. **Placeholders.** None: every code step has code; run steps have commands and expected output. The artifact template's visual styling is delegated to the `artifact-design`/`dataviz` skills by contract (Task 12 Step 2), with the data-driven rendering core written out and syntax-checked.
3. **Type consistency.** All Python in this plan was extracted into a scratch directory and run: 65 unit tests pass and every module compiles (`render_core.js` passes `node --check`). Names agree across tasks: `Rec`, `read_records`, `place`, `classify_fate`, `paf_hits`, `best_hits` (Task 1) are used with the same signatures in Tasks 4, 5, 9, 11; loci dicts written by Task 4 (`locus kind family chrom start end name`) are read by Tasks 5, 9, 11; `cands.tsv` columns are identical in Tasks 7 and 8 and parsed in Tasks 9 and 11; `fate.json` keys written in Task 5 are read in Tasks 9, 11, 12; `side.json`/`compare_*.json` keys written in Task 11 are read in Task 12.
4. **Review Focus.** The seven items map to tests: 1 -> `test_overlap_read_keeps_family_label`; 2 -> `test_unknown_name_passes_through`; 3 -> `test_zero_reads_has_no_verdict`; 4 -> `test_no_record_is_unmapped` plus the selection note printed by `fate.py run`; 5 -> `test_empty_candidates_fail_r1`, `test_sex_locus_excluded`; 6 -> `test_chry_is_pat_only_and_sex`, `test_poor_lift_means_no_counterpart`; 7 -> `test_pairs_need_a_primary_on_both`.
5. **Decisions this plan takes that the prereg left open (for the user's review):** (a) the in-house arm's flag floor is `n_clusters >= 2` (the analogue of "at least 2 transcripts"); (b) R3's "class b/c" is every flagged candidate that does not recover a truth locus, because with one haplotype as the only reference the registered `_pri`-based b/c split does not apply; (c) the bar's "unmapped or tied" is read as the two lost classes together; (d) the shared-copy control is every read with an untied placement on the truth haplotype that is on no absent locus (a superset of "copies with an ortholog"); (e) the divergence pile bins each interval into 60 equal fractions, so the two panels share a relative axis, not coordinates; (f) T1-candidate copies are not in either run's truth set.
