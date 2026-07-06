#!/usr/bin/env python
"""Self-check for bench/family_rna_refine.py (the DEFAULT-ON RNA-only refinement stage).

Asserts:
  (a) LEGACY opt-out (--legacy / RUSTLE_RNA_ORACLE=0) writes nothing; shipped path untouched.
  (b) DEFAULT (no flag) writes the refined catalog and is deterministic (== --rna-oracle == env=1).
  (c) The allele demotion removes DHRSX + LOC129530050 (DNA-confirmed).
  (d) The residual-removed count matches RNA_ONLY_EDGE_ORACLE.md's recall-preserving row:
      residual remaining {allele 0, oversize 3, multifam 4}; shipped total 12; 6/12 named FP removed.
  (e) REPEAT-HUB GATE negative controls: the fam17 16-family repeat-bridge hub is SHATTERED (no
      longer one block) by default, while GSTM2 + MAGE families are SPARED (largest block byte-
      identical membership with vs without the gate).
  (f) --no-repeat-gate / RUSTLE_NO_REPEAT_GATE=1 ABLATION recovers the pre-gate catalog (fam17 one
      block again) and is itself byte-identical across runs.
  (g) recall cost: EVERY gene-pair the repeat gate cuts is a VG-'genuine' over-merge -- 0 TP and 0
      truthbar (real / borderline-real) edges are cut.
  (h) RNA-only / LIBRARY-FREE guard: repeat-hub features = {min_shared_mult, cyclic}, disjoint from
      every soft-mask / RepeatMasker / DNA column.
  (i) DEFAULT (no flag) catalog is BYTE-IDENTICAL to the golden TSV (md5) and its JSON records
      gamma=0.20 with NO high_precision block.
  (j) --high-precision (== env RUSTLE_HIGH_PRECISION=1) swaps gamma 0.20 -> 0.40
      (PRECISION_RECALL_FRONTIER.md): JSON records the active gamma, the catalog is DIFFERENT +
      deterministic, n_families -> 623 (frontier gamma=0.40 row), the collapsed-array OVERSIZE blobs
      MPHOSPH8 + LOC134758618 are removed (oversize residual 3 -> 1), and the HONEST off-oracle
      KRAB-ZNF + MAGE-floor caveats are carried in the JSON (not dropped).
  (k) --legacy still opts out even with --high-precision; --high-precision composes with
      --no-repeat-gate (gamma=0.40 AND the gate ablated).
  (p) ANTISENSE / reciprocal-overlap gate (4th default-on FP gate): DEFAULT (gate ON) dissolves the
      named same-locus opposite-strand nested-gene over-merges (CCNH+RASA1, KAT5+RNASEH2C, EXOSC3+
      TRMT10B, ARHGEF39+CCDC107, HDGFL3+TM6SF1 no longer co-membered) while SPARING the real antisense
      pair MPDU1/MPDU1-AS1 (reciprocal-overlap 0.4855 < 0.50) and, via the MEGA-SPAN guard, the GSTM2
      array (largest block byte-identical with vs without the gate).  --no-antisense-gate /
      RUSTLE_NO_ANTISENSE_GATE=1 recovers the pre-antisense catalog (00c4ff2e) BYTE-IDENTICAL, and the
      gate composes with every other --no-* flag.
  (q) COLINEAR-MERGE gate (5th default-on stage): post-demote family merge by exon colinearity +
      adaptive adjacent-junction floor.  --no-colinear-merge / RUSTLE_NO_COLINEAR_MERGE=1 recovers the
      pre-merge catalog (548029ad) BYTE-IDENTICAL.  DEFAULT keeps MAGEA9 + GSTM1/2/4/5 merged while
      blocking the ANKRD18 + ANKRD36C domain-share over-merge.

Run: /home/juanfra/miniforge3/bin/python bench/test_family_rna_refine.py
"""
import hashlib
import json
import os
import subprocess
import sys
from collections import defaultdict

BENCH = os.path.dirname(os.path.abspath(__file__))
SCRIPT = os.path.join(BENCH, "family_rna_refine.py")
OUT_TSV = os.path.join(BENCH, "family_rna_refine.tsv")
OUT_JSON = os.path.join(BENCH, "family_rna_refine.json")
SHIPPED = os.path.join(BENCH, "denovo_families.py")
PY = sys.executable

# md5 of the DEFAULT catalog (gamma=0.20, repeat-hub + ANTISENSE + recombinant-split + multi-repeat-bridge
# + COLINEAR-MERGE ALL FIVE stages ON, 553 families).  The 5th stage (colinear merge) is DEFAULT-ON and
# merges split families by exon colinearity after demote.  The --high-precision flag must NOT perturb the
# default: `family_rna_refine.py` with no flag reproduces this.
GOLDEN_DEFAULT_TSV_MD5 = "991913da4345e88bf56fbd939d753292"
# md5 of the PRE-COLINEAR-MERGE catalog (--no-colinear-merge: 566 families, the 4-gate default before
# the colinear-merge stage was added).  Ablating ONLY the colinear-merge stage recovers this EXACTLY
# (byte-for-byte) -- the byte-identical-when-OFF discipline.
GOLDEN_NO_COLINEAR_MERGE_TSV_MD5 = "548029ade983645f71c772a9d9257334"
# md5 of the PRE-ANTISENSE catalog (--no-antisense-gate: gamma=0.20, repeat-hub + recombinant-split +
# multi-repeat-bridge + colinear-merge ON, antisense OFF, 560 families) = the shipped default WITHOUT the
# antisense gate.  The antisense gate is DEFAULT-ON; ablating ONLY it must recover this EXACTLY.
GOLDEN_NO_ANTISENSE_TSV_MD5 = "00c4ff2ea0d7afb8b8505d68487404a2"
# md5 of the pre-repeat-bridge catalog WITH the antisense gate ON (--no-repeat-bridge-gate: gamma=0.20,
# repeat-hub + antisense + recombinant-split + colinear-merge ON, multi-repeat-bridge OFF, 581 families).
GOLDEN_ANTISENSE_NO_REPEAT_BRIDGE_TSV_MD5 = "0da31d5925e1087d051f1133660ec93e"
# md5 of the pre-repeat-bridge catalog with BOTH the antisense AND repeat-bridge gates OFF
# (--no-antisense-gate --no-repeat-bridge-gate: 590 families).
GOLDEN_NO_REPEAT_BRIDGE_TSV_MD5 = "5f438104f2f9225b734c3fba1695f1c3"
# md5 of the pre-split catalog with antisense + repeat-bridge OFF (--no-antisense-gate
# --no-split-recombinants --no-repeat-bridge-gate, gamma=0.20, repeat-hub + colinear-merge ON, 591 families).
GOLDEN_NOSPLIT_TSV_MD5 = "8f79718a79768b45b011efd883a2b207"

# fam17 = the 27-gene / 16-protein-family EXTREME repeat-bridge hub (VG_REPEAT_CATALOG.md sec 3).
FAM17 = ["C9H11orf65", "COPS9", "DCAF16", "DCP1B", "GIMAP4", "KLHL33", "LIPA", "LOC109024560",
         "LOC109025730", "LOC109028669", "LOC129530096", "LOC129530131", "LOC134756589", "LRFN5",
         "NCALD", "NCR3LG1", "NFXL1", "OXCT1", "PDCD5", "SLC11A2", "SMYD3", "SPATA33", "SPC25",
         "THNSL1", "TSPAN8", "UCHL1", "TMEM38B"]
# negative controls: GSTM2 (protein DOMAIN, per-edge max mult 9) + MAGE (cardinality, max mult 8).
GSTM2 = ["GSTM2", "LOC115930164", "LOC115930576", "LOC101126097", "LOC134756922", "LOC101129940",
         "SEC22B", "LOC109023809", "TRIP13", "LOC129532045", "LOC134757399", "LOC101131274",
         "LOC101135585", "LOC115932984", "LOC115933235", "LOC115933241", "LOC129525330",
         "LOC129525331", "LOC129525599", "LOC134756861"]
MAGE = ["MAGEA1", "MAGEA4", "MAGEA12", "LOC129529976", "LOC129529978", "LOC129529983",
        "LOC129529986", "MAGEA9", "LOC129530018", "MAGEB6", "MAGEB6B"]
# ANTISENSE gate: the 5 NAMED genuine over-merges (same-locus opposite-strand nested genes) the gate must
# DISSOLVE (FP_EXCLUSION_DISCRIMINATORS.md) + the REAL antisense pair it must SPARE (reciprocal-overlap
# 0.4855 < 0.50 floor).  GSTM2 (1.18 Mb array span >= the 500 kb MEGA-SPAN guard) must be UNTOUCHED.
ANTISENSE_NAMED_FP = [("CCNH", "RASA1"), ("KAT5", "RNASEH2C"), ("EXOSC3", "TRMT10B"),
                      ("ARHGEF39", "CCDC107"), ("HDGFL3", "TM6SF1")]
ANTISENSE_SPARED_PAIR = ("MPDU1", "MPDU1-AS1")   # reciprocal-overlap 0.4855 -> below the 0.50 floor


# WSL /mnt/c (drvfs/9p) can lag file existence + content visibility after a CHILD process writes or
# removes a file -- subprocess.run returns before the parent's immediate os.path.exists / json.load sees
# the new state.  These are filesystem visibility races, NOT logic defects.  _settle() / _writes() below
# poll until the just-written (or just-removed) outputs are visible + self-consistent so the suite is
# deterministic on this filesystem.  (No-op on POSIX where close() already guarantees visibility.)
import time

_SETTLE_TRIES = 300
_SETTLE_DELAY = 0.05


def _writes(args, env):
    """True iff this invocation is EXPECTED to write the catalog (i.e. NOT the legacy opt-out / --no-write)."""
    if "--legacy" in args or "--no-write" in args:
        return False
    if env.get("RUSTLE_RNA_ORACLE") == "0":
        return False
    return True


def _settle(path, want_exists, tries=_SETTLE_TRIES, delay=_SETTLE_DELAY):
    """Poll until os.path.exists(path) == want_exists (bounded).  Returns the final observed state."""
    for _ in range(tries):
        if os.path.exists(path) == want_exists:
            return want_exists
        time.sleep(delay)
    return os.path.exists(path)


def _settle_json(path, tries=_SETTLE_TRIES, delay=_SETTLE_DELAY):
    """Poll until `path` exists AND parses as JSON (a child mid-write can leave a truncated file)."""
    for _ in range(tries):
        try:
            with open(path) as fh:
                json.load(fh)
            return True
        except (OSError, ValueError):
            time.sleep(delay)
    return False


def _run(args, env_extra=None):
    env = dict(os.environ)
    env.pop("RUSTLE_RNA_ORACLE", None)
    env.pop("RUSTLE_NO_REPEAT_GATE", None)
    env.pop("RUSTLE_NO_SPLIT_RECOMBINANTS", None)
    env.pop("RUSTLE_NO_REPEAT_BRIDGE_GATE", None)
    env.pop("RUSTLE_NO_ANTISENSE_GATE", None)
    if env_extra:
        env.update(env_extra)
    r = subprocess.run([PY, SCRIPT] + args, cwd=BENCH, env=env,
                       capture_output=True, text=True)
    # WSL visibility settle: wait for the freshly-written outputs to be visible + parseable before the
    # parent reads them (mirrors what the shipped repeat-hub / recombinant-split gate tests rely on).
    if r.returncode == 0 and _writes(args, env):
        _settle(OUT_TSV, True)
        _settle_json(OUT_JSON)
    return r


def _remove_settled(*paths):
    """Remove `paths` and poll until each is actually gone (drvfs unlink-visibility lag)."""
    for p in paths:
        if os.path.exists(p):
            os.remove(p)
    for p in paths:
        _settle(p, False)


def _md5(path, tries=_SETTLE_TRIES, delay=_SETTLE_DELAY):
    for _ in range(tries):
        try:
            with open(path, "rb") as fh:
                return hashlib.md5(fh.read()).hexdigest()
        except OSError:
            time.sleep(delay)
    with open(path, "rb") as fh:
        return hashlib.md5(fh.read()).hexdigest()


def _fam_to_genes(path):
    fam = defaultdict(set)
    with open(path) as fh:
        hdr = fh.readline().rstrip("\n").split("\t")
        ix = {h: i for i, h in enumerate(hdr)}
        for ln in fh:
            f = ln.rstrip("\n").split("\t")
            fam[f[ix["family_id"]]].add(f[ix["member_gene"]])
    return fam


def _blocks_with(fam, genes):
    gs = set(genes)
    return [members for members in fam.values() if members & gs]


def test_legacy_opt_out_writes_nothing():
    _remove_settled(OUT_TSV, OUT_JSON)              # settle the unlink (drvfs visibility) before asserting absence
    shipped_before = _md5(SHIPPED)
    r = _run(["--legacy"])                          # legacy opt-out (flag form)
    assert r.returncode == 0, f"legacy exit={r.returncode}"
    assert r.stdout.startswith("legacy:"), repr(r.stdout)
    assert not os.path.exists(OUT_TSV) and not os.path.exists(OUT_JSON), "legacy wrote outputs"
    assert _md5(SHIPPED) == shipped_before, "denovo_families.py (shipped path) was modified"
    _remove_settled(OUT_TSV, OUT_JSON)
    r2 = _run([], env_extra={"RUSTLE_RNA_ORACLE": "0"})   # legacy opt-out (env form)
    assert r2.returncode == 0 and r2.stdout.startswith("legacy:"), repr(r2.stdout)
    assert not os.path.exists(OUT_TSV), "env-legacy wrote outputs"
    print("(a) legacy opt-out (--legacy / env=0) writes nothing / shipped untouched : OK")


def test_default_deterministic():
    r1 = _run([])                                  # DEFAULT: no flag -> refined catalog
    assert r1.returncode == 0, r1.stderr
    assert os.path.exists(OUT_TSV) and os.path.exists(OUT_JSON), "default did not write the catalog"
    tsv1, json1 = _md5(OUT_TSV), _md5(OUT_JSON)
    r2 = _run(["--rna-oracle"])                     # deprecated no-op == default
    assert r2.returncode == 0, r2.stderr
    assert _md5(OUT_TSV) == tsv1 and _md5(OUT_JSON) == json1, "default != --rna-oracle"
    r3 = _run([], env_extra={"RUSTLE_RNA_ORACLE": "1"})   # explicit-enable == default
    assert r3.returncode == 0, r3.stderr
    assert _md5(OUT_TSV) == tsv1, "TSV not byte-identical across runs"
    assert _md5(OUT_JSON) == json1, "JSON not byte-identical across runs"
    print("(b) default deterministic (default == --rna-oracle == env=1; byte-identical) : OK")


def test_default_byte_identical_to_golden():
    r = _run([])            # DEFAULT: no flag (repeat-hub + antisense + recombinant-split + multi-repeat-bridge + colinear-merge ALL ON)
    assert r.returncode == 0, r.stderr
    got = _md5(OUT_TSV)
    assert got == GOLDEN_DEFAULT_TSV_MD5, \
        f"default catalog drifted from the golden {GOLDEN_DEFAULT_TSV_MD5}: {got}"
    d = json.load(open(OUT_JSON))
    assert d["rule"]["gamma"] == 0.2, f"default JSON gamma != 0.20: {d['rule']['gamma']}"
    assert d["guards"]["gamma"] == 0.2, f"default guards gamma != 0.20: {d['guards']['gamma']}"
    assert d["n_families"] == 553, f"default n_families != 553: {d['n_families']}"
    assert "high_precision" not in d, "default JSON must NOT carry the high_precision block"
    # recombinant-split gate is DEFAULT-ON and splits the 2 recall-safe HIGH-confidence mosaics
    rs = d["recombinant_split"]
    assert rs["enabled"] is True and rs["n_families_split"] == 2, \
        f"default recombinant-split not ON/2: enabled={rs['enabled']} n_split={rs['n_families_split']}"
    assert d["rule"]["split_recombinants_enabled"] is True
    # multi-repeat-bridge gate (3rd VG-native gate) is DEFAULT-ON at T=8/C=2 and cuts the DISCONNECTED
    # repeat-bridge FP class (all cut families are no-full-shared-exon by construction).  The antisense
    # gate runs UPSTREAM (edge-level) so the count is 59 (was 61 before the antisense gate).
    mrb = d["multi_repeat_bridge"]
    assert mrb["enabled"] is True, "default multi-repeat-bridge gate not enabled"
    assert d["rule"]["repeat_bridge_gate_enabled"] is True
    assert (mrb["mult_min"], mrb["count_min"], mrb["exon_id"]) == (8, 2, 0.70), \
        f"repeat-bridge T/C/exon_id drifted: {mrb['mult_min']}/{mrb['count_min']}/{mrb['exon_id']}"
    assert mrb["n_families_cut"] == 59, f"default repeat-bridge n_families_cut != 59: {mrb['n_families_cut']}"
    assert mrb["n_disconnected_removed"] == mrb["n_families_cut"], "disconnected_removed != n_families_cut"
    assert len(mrb["families_cut_named"]) == mrb["n_families_cut"]
    # antisense / reciprocal-overlap gate (4th default-on FP gate) is DEFAULT-ON and cuts the same-locus
    # opposite-strand nested-gene over-merges (16 cross-gene pairs incl. the 5 named FP)
    ag = d["antisense_overlap_gate"]
    assert ag["enabled"] is True, "default antisense gate not enabled"
    assert d["rule"]["antisense_gate_enabled"] is True
    assert (d["rule"]["antisense_recip_overlap_min"], d["rule"]["mega_span_max"]) == (0.5, 500000), \
        f"antisense thresholds drifted: {d['rule']['antisense_recip_overlap_min']}/{d['rule']['mega_span_max']}"
    assert ag["n_cross_gene_pairs_cut"] == 16, f"antisense n_cross_gene_pairs_cut != 16: {ag['n_cross_gene_pairs_cut']}"
    for a, b in ANTISENSE_NAMED_FP:
        assert "+".join(sorted((a, b))) in ag["pairs_cut"], f"named antisense FP {a}+{b} not in pairs_cut"
    # colinear-merge stage (5th default-on stage) is DEFAULT-ON and records merge edges
    cm = d["colinear_merge"]
    assert cm["enabled"] is True, "default colinear merge not enabled"
    assert cm["n_merge_edges"] >= 1, f"colinear merge produced no edges: {cm['n_merge_edges']}"
    assert cm["min_colinear"] == 2 and cm["min_adjacent_junctions"] == 2
    print(f"(i) default (no flag) byte-identical to golden TSV (md5 {GOLDEN_DEFAULT_TSV_MD5}); "
          f"gamma=0.20, {d['n_families']} families, recombinant-split ON (2 split), multi-repeat-bridge ON "
          f"({mrb['n_families_cut']} cut, T=8/C=2), antisense ON ({ag['n_cross_gene_pairs_cut']} pairs cut), "
          f"colinear-merge ON ({cm['n_merge_edges']} edges), no high_precision block : OK")


def test_recombinant_split_gate():
    """recombinant-split DEFAULT-ON breaks the fid-210 mosaic (GALNT17 + LOC101126070 no longer
    co-membered; GALNT17's single-locus remnant drops via the >=2-loci multi-copy filter, LOC101126070
    survives); --no-split-recombinants (+ --no-repeat-bridge-gate + --no-antisense-gate) recovers the
    591-family pre-split golden byte-identically; env==flag.  Byte-identity to the pre-split golden IS
    the recall-safety proof (the split never separates same-gene loci)."""
    # PRE-SPLIT reference (split + bridge + antisense all OFF; colinear-merge stays DEFAULT-ON):
    # fid co-members GALNT17 + LOC101126070.
    rp = _run(["--no-split-recombinants", "--no-repeat-bridge-gate", "--no-antisense-gate"])
    assert rp.returncode == 0, rp.stderr
    pre = _fam_to_genes(OUT_TSV)
    assert any({"GALNT17", "LOC101126070"} <= gs for gs in pre.values()), \
        "pre-split: GALNT17 + LOC101126070 not co-membered (expected the fid-210 mosaic)"
    off = _md5(OUT_TSV);  d_off = json.load(open(OUT_JSON))
    assert off == GOLDEN_NOSPLIT_TSV_MD5, f"--no-split-recombinants drifted from pre-split golden: {off}"
    assert d_off["n_families"] == 591 and d_off["recombinant_split"]["enabled"] is False
    assert d_off["multi_repeat_bridge"]["enabled"] is False and d_off["antisense_overlap_gate"]["enabled"] is False
    # env form == flag form (byte-identical) -- proves same-gene loci never separated (recall-safe)
    r2 = _run(["--no-repeat-bridge-gate", "--no-antisense-gate"],
              env_extra={"RUSTLE_NO_SPLIT_RECOMBINANTS": "1"})
    assert r2.returncode == 0 and _md5(OUT_TSV) == off, "RUSTLE_NO_SPLIT_RECOMBINANTS=1 != flag form"
    # DEFAULT (split ON): the fid-210 mosaic is BROKEN -> GALNT17 + LOC101126070 no longer co-membered
    r = _run([]);  assert r.returncode == 0, r.stderr
    fam = _fam_to_genes(OUT_TSV)
    assert not any({"GALNT17", "LOC101126070"} <= gs for gs in fam.values()), \
        "recombinant mosaic fid 210 NOT split (GALNT17 + LOC101126070 still co-membered)"
    d = json.load(open(OUT_JSON))
    assert d["recombinant_split"]["enabled"] is True and d["recombinant_split"]["n_families_split"] == 2
    print("(l) recombinant-split ON breaks the fid-210 mosaic (GALNT17|LOC101126070 no longer "
          "co-membered); --no-split-recombinants (+--no-repeat-bridge-gate +--no-antisense-gate) "
          "recovers 591 pre-split golden; env==flag : OK")
    _run([])   # restore DEFAULT


def test_antisense_gate():
    """(p) ANTISENSE / reciprocal-overlap gate (4th default-on FP gate, edge-level GENOME-ARCHITECTURE):
    (a) --no-antisense-gate (flag == env) recovers the PRE-ANTISENSE golden 00c4ff2e BYTE-IDENTICAL
        (560 families; the byte-identical-when-OFF discipline); (b) DEFAULT (gate ON) dissolves the 5
        named same-locus opposite-strand nested-gene over-merges while SPARING the real antisense pair
        MPDU1/MPDU1-AS1 (recip-overlap 0.4855 < 0.50) and the GSTM2 array (MEGA-SPAN guard -> largest
        block byte-identical); (c) deterministic (default byte-identical across runs)."""
    # (a) ablation (flag) recovers the pre-antisense golden 00c4ff2e, byte-identical
    r = _run(["--no-antisense-gate"]);  assert r.returncode == 0, r.stderr
    off_tsv = _md5(OUT_TSV);  off = _fam_to_genes(OUT_TSV);  off_json = json.load(open(OUT_JSON))
    assert off_tsv == GOLDEN_NO_ANTISENSE_TSV_MD5, \
        f"--no-antisense-gate drifted from the pre-antisense golden {GOLDEN_NO_ANTISENSE_TSV_MD5}: {off_tsv}"
    assert off_json["antisense_overlap_gate"]["enabled"] is False
    assert off_json["rule"]["antisense_gate_enabled"] is False
    assert off_json["antisense_overlap_gate"]["n_cross_gene_pairs_cut"] == 0, "ablation still cut edges"
    assert off_json["n_families"] == 560, f"--no-antisense-gate n_families != 560: {off_json['n_families']}"
    # env form == flag form (byte-identical) -- the byte-identical-when-OFF discipline
    r2 = _run([], env_extra={"RUSTLE_NO_ANTISENSE_GATE": "1"});  assert r2.returncode == 0, r2.stderr
    assert _md5(OUT_TSV) == off_tsv, "RUSTLE_NO_ANTISENSE_GATE=1 != --no-antisense-gate"
    # DEFAULT (gate ON) -- (c) deterministic across two runs
    r3 = _run([]);  assert r3.returncode == 0, r3.stderr
    on_tsv = _md5(OUT_TSV);  on = _fam_to_genes(OUT_TSV);  on_json = json.load(open(OUT_JSON))
    r4 = _run([]);  assert r4.returncode == 0, r4.stderr
    assert _md5(OUT_TSV) == on_tsv, "default not byte-identical across runs (non-deterministic)"
    assert on_tsv == GOLDEN_DEFAULT_TSV_MD5 and on_tsv != off_tsv, "antisense gate had no effect / drifted"
    # (b) the 5 named genuine over-merges are DISSOLVED (no longer co-membered) by DEFAULT
    for a, b in ANTISENSE_NAMED_FP:
        assert any({a, b} <= gs for gs in off.values()), \
            f"pre-antisense: {a}+{b} expected co-membered (the over-merge to dissolve)"
        assert not any({a, b} <= gs for gs in on.values()), \
            f"antisense gate did NOT dissolve the {a}+{b} over-merge (still co-membered)"
    # (b) the REAL antisense pair MPDU1/MPDU1-AS1 (recip 0.4855 < 0.50) is SPARED -- still co-membered
    sa, sb = ANTISENSE_SPARED_PAIR
    assert any({sa, sb} <= gs for gs in on.values()), \
        f"antisense gate wrongly cut the real pair {sa}+{sb} (recip 0.4855 < 0.50 floor)"
    # (b) GSTM2 MEGA-SPAN guard: its largest block is byte-identical membership with vs without the gate
    bo = max(_blocks_with(off, GSTM2), key=len)
    bn = max(_blocks_with(on, GSTM2), key=len)
    assert bo == bn, f"GSTM2 largest block CHANGED by the antisense gate (mega-span guard failed): {bo ^ bn}"
    ag = on_json["antisense_overlap_gate"]
    print(f"(p) antisense gate: --no-antisense-gate recovers pre-antisense golden 00c4ff2e byte-identical "
          f"(flag==env, 560 fam); DEFAULT dissolves the 5 named FP ({ag['n_cross_gene_pairs_cut']} pairs "
          f"cut), spares MPDU1/MPDU1-AS1 + GSTM2 (mega-span); deterministic : OK")
    _run([])   # restore DEFAULT


def test_colinear_merge_gate():
    """(q) COLINEAR-MERGE stage (5th default-on): post-demote family merge by exon colinearity with an
    adaptive adjacent-junction floor.  --no-colinear-merge recovers the pre-merge catalog 548029ad
    byte-identical (flag == env).  DEFAULT keeps known split-family siblings (MAGEA9 rejoins the MAGEA
    family; GSTM1/2/4/5 rejoin) while blocking the ANKRD18 + ANKRD36C domain-share over-merge."""
    # (a) ablation (flag) recovers the pre-colinear-merge golden 548029ad, byte-identical
    r = _run(["--no-colinear-merge"]);  assert r.returncode == 0, r.stderr
    off_tsv = _md5(OUT_TSV);  off = _fam_to_genes(OUT_TSV);  off_json = json.load(open(OUT_JSON))
    assert off_tsv == GOLDEN_NO_COLINEAR_MERGE_TSV_MD5, \
        f"--no-colinear-merge drifted from the pre-merge golden {GOLDEN_NO_COLINEAR_MERGE_TSV_MD5}: {off_tsv}"
    assert off_json["colinear_merge"]["enabled"] is False
    assert off_json["n_families"] == 566, f"--no-colinear-merge n_families != 566: {off_json['n_families']}"
    # env form == flag form (byte-identical)
    r2 = _run([], env_extra={"RUSTLE_NO_COLINEAR_MERGE": "1"});  assert r2.returncode == 0, r2.stderr
    assert _md5(OUT_TSV) == off_tsv, "RUSTLE_NO_COLINEAR_MERGE=1 != --no-colinear-merge"
    # DEFAULT (merge ON): deterministic across two runs
    r3 = _run([]);  assert r3.returncode == 0, r.stderr
    on_tsv = _md5(OUT_TSV);  on = _fam_to_genes(OUT_TSV);  on_json = json.load(open(OUT_JSON))
    r4 = _run([]);  assert r4.returncode == 0, r.stderr
    assert _md5(OUT_TSV) == on_tsv, "default not byte-identical across runs (non-deterministic)"
    assert on_tsv == GOLDEN_DEFAULT_TSV_MD5 and on_tsv != off_tsv, "colinear merge had no effect / drifted"
    cm = on_json["colinear_merge"]
    assert cm["enabled"] is True and cm["n_merge_edges"] >= 1
    # (b) MAGEA9 rejoins the MAGEA family with the merge ON
    assert any({"MAGEA9", "MAGEA1"} <= gs for gs in on.values()), \
        "DEFAULT: MAGEA9 did not rejoin the MAGEA family"
    # (c) GSTM1/2/4/5 are co-membered in the GSTM block with the merge ON
    assert any({"GSTM1", "GSTM2"} <= gs for gs in on.values()), \
        "DEFAULT: GSTM1 and GSTM2 are not co-membered after colinear merge"
    # (d) ANKRD18 + ANKRD36C domain-share is BLOCKED (not co-membered) with the merge ON
    assert not any({"ANKRD18A", "ANKRD18B"} & gs and "ANKRD36C" in gs for gs in on.values()), \
        "DEFAULT: ANKRD18 + ANKRD36C domain-share over-merge was not blocked"
    print(f"(q) colinear merge: --no-colinear-merge recovers pre-merge golden 548029ad byte-identical "
          f"(flag==env, 566 fam); DEFAULT ({cm['n_merge_edges']} edges) keeps MAGEA9 + GSTM1/2/4/5 "
          f"merged, blocks ANKRD18 + ANKRD36C; deterministic : OK")
    _run([])   # restore DEFAULT


def test_allele_demote_removes_known_fp():
    d = json.load(open(OUT_JSON))
    demoted = {x["gene"] for x in d["alleles_demoted"]}
    conf = {x["gene"] for x in d["alleles_demoted"] if x["dna_confirmed"]}
    assert "DHRSX" in demoted and "LOC129530050" in demoted, f"missing allele demotions: {demoted}"
    assert {"DHRSX", "LOC129530050"} <= conf, f"DHRSX/LOC129530050 not flagged DNA-confirmed: {conf}"
    print(f"(c) allele demotion removes DHRSX + LOC129530050 : OK  (demoted={sorted(demoted)})")


def test_residual_removed_matches_recall_preserving_row():
    # the recall-preserving row is the HISTORICAL PRE-repeat-bridge catalog (--no-repeat-bridge-gate
    # --no-colinear-merge), which the multi-repeat-bridge gate does not touch: residual
    # {allele 0, oversize 3, multifam 4}, 6/12 removed.
    r = _run(["--no-repeat-bridge-gate", "--no-colinear-merge"]);  assert r.returncode == 0, r.stderr
    d = json.load(open(OUT_JSON))
    rem = d["residual_fp"]["residual_remaining"]
    assert rem == {"allele": 0, "oversize": 3, "multifam": 4}, f"residual remaining drift: {rem}"
    assert d["residual_fp"]["shipped_residual_total"] == 12, d["residual_fp"]["shipped_residual_total"]
    assert d["residual_fp"]["named_removed"] == 6, d["residual_fp"]["named_removed"]
    bd = d["residual_fp"]["named_removed_breakdown"]
    assert bd == {"allele": 2, "oversize": 2, "multifam": 2}, f"removed breakdown drift: {bd}"
    print("(d) historical pre-repeat-bridge residual {allele 0, oversize 3, multifam 4}; 6/12 removed : OK")
    # the DEFAULT (multi-repeat-bridge ON) removes 2 MORE oversize blobs (MPHOSPH8-hub + LOC134758618)
    r2 = _run([]);  assert r2.returncode == 0, r2.stderr
    d2 = json.load(open(OUT_JSON))
    rem2 = d2["residual_fp"]["residual_remaining"]
    assert rem2 == {"allele": 0, "oversize": 1, "multifam": 4}, f"default residual drift: {rem2}"
    assert d2["residual_fp"]["named_removed"] == 7, d2["residual_fp"]["named_removed"]
    print("(d') DEFAULT (multi-repeat-bridge + colinear-merge ON) residual {allele 0, oversize 1, multifam 4}; 7/12 "
          "removed (MPHOSPH8-hub + LOC134758618 oversize blobs cut as DISCONNECTED; colinear-merge "
          "suppresses one multifam FP) : OK")


def test_repeat_gate_shatters_fam17_spares_controls():
    # isolate the REPEAT-HUB gate: hold the recombinant-split, multi-repeat-bridge AND antisense gates OFF
    # in both arms so the only difference is the repeat-hub gate (each other gate is exercised separately).
    base = ["--no-split-recombinants", "--no-repeat-bridge-gate", "--no-antisense-gate"]
    r = _run(base)   # repeat-hub ON, split+bridge+antisense OFF
    assert r.returncode == 0, r.stderr
    on = _fam_to_genes(OUT_TSV);  on_json = json.load(open(OUT_JSON))
    on_tsv = _md5(OUT_TSV)
    # ablation = repeat-hub gate OFF (flag form) -- run twice for byte-identical determinism
    r2 = _run(["--no-repeat-gate"] + base)
    assert r2.returncode == 0, r2.stderr
    off = _fam_to_genes(OUT_TSV);  off_json = json.load(open(OUT_JSON))
    off_tsv = _md5(OUT_TSV)
    r3 = _run(base, env_extra={"RUSTLE_NO_REPEAT_GATE": "1"})   # env form == flag form
    assert r3.returncode == 0 and _md5(OUT_TSV) == off_tsv, "env ablation != flag ablation / not deterministic"

    # the gate MUST change the catalog, and only via the repeat gate
    assert on_tsv != off_tsv, "repeat gate had no effect on the catalog"
    assert on_json["edges"]["n_dn_cross_gene_cut_by_repeat_gate"] > 0, "gate cut nothing (ON)"
    assert off_json["edges"]["n_dn_cross_gene_cut_by_repeat_gate"] == 0, "ablation still cut edges"

    # (a) fam17: ONE block pre-gate -> shattered (>=2 blocks / no longer one block) with the gate
    off_blocks, on_blocks = _blocks_with(off, FAM17), _blocks_with(on, FAM17)
    assert len(off_blocks) == 1, f"fam17 not a single block pre-gate: {len(off_blocks)} blocks"
    assert len(on_blocks) >= 2, f"repeat gate did NOT shatter fam17: still {len(on_blocks)} block(s)"

    # (a) GSTM2 + MAGE SPARED: their largest block is byte-identical membership with vs without the gate
    for name, gl in (("GSTM2", GSTM2), ("MAGE", MAGE)):
        bo = max(_blocks_with(off, gl), key=len)
        bn = max(_blocks_with(on, gl), key=len)
        assert bo == bn, f"{name} largest block CHANGED by the repeat gate (not spared): {bo ^ bn}"
    print(f"(e) repeat gate shatters fam17 (1 -> {len(on_blocks)} blocks); GSTM2 + MAGE spared "
          f"(largest block identical) : OK")
    print("(f) --no-repeat-gate / env ablation recovers the one-block fam17 + byte-identical : OK")
    # restore DEFAULT (gate ON) state for any downstream test that reads the catalog
    _run([])


def test_repeat_gate_recall_cost_is_genuine_only():
    """Every gene-pair the repeat gate cuts is a VG-'genuine' over-merge; 0 TP / 0 truthbar cut."""
    sys.path.insert(0, BENCH)
    import family_rna_refine as R
    gene_of_dn = R.RO.load_gene_of_dn()
    pair_core = R.RO.load_pair_core(gene_of_dn)
    univ_aln = R.RO.load_universal_aln()
    rep = R.load_repeat_mult()
    # VG per-edge class
    cls = {}
    with open(os.path.join(BENCH, "vg_repeat_catalog.tsv")) as fh:
        in_e, ix = False, None
        for ln in fh:
            if ln.startswith("# SECTION edges"):
                in_e, ix = True, None; continue
            if not in_e:
                continue
            if ln.startswith("gene_a\t"):
                hdr = ln.rstrip("\n").split("\t"); ix = {h: i for i, h in enumerate(hdr)}; continue
            if ix is None:
                continue
            f = ln.rstrip("\n").split("\t")
            cls[frozenset((f[ix["gene_a"]], f[ix["gene_b"]]))] = f[ix["cls"]]

    def core_aln_keep(k):
        return (pair_core.get(k) or 0.0) >= R.CORE_MIN and (univ_aln.get(k) or 0.0) >= R.ALN_MIN

    cut = [k for k, m in rep.items() if m >= R.REPEAT_MULT_MIN and core_aln_keep(k)]
    labels = [cls.get(k, "NA") for k in cut]
    assert cut, "repeat gate cut no core+aln-surviving edges"
    assert all(l == "genuine" for l in labels), \
        f"repeat gate cut a real/borderline edge: {[sorted(k) for k, l in zip(cut, labels) if l != 'genuine']}"
    print(f"(g) recall cost: all {len(cut)} repeat-gate-cut pairs are VG-'genuine' over-merge "
          f"(0 TP / 0 truthbar cut) : OK")


def test_rna_only_guard():
    d = json.load(open(OUT_JSON))
    g = d["guards"]
    assert g["edge_decision_features"] == ["core_recip", "aln_frac"], g["edge_decision_features"]
    assert g["repeat_gate_features"] == ["min_shared_mult", "cyclic"], g["repeat_gate_features"]
    assert g["antisense_gate_features"] == ["gene_strand", "gene_span"], g["antisense_gate_features"]
    assert g["demote_features"] == ["balanced_frac", "copy_like"], g["demote_features"]
    assert g["no_dna_in_inference"] is True and g["gamma"] == 0.2 and g["seed"] == 0
    assert g["repeat_gate_library_free"] is True and g["no_softmask_in_repeat_gate"] is True
    assert g["antisense_gate_genome_architecture"] is True and g["no_dna_in_antisense_gate"] is True
    # source-level: the repeat-hub feature set is library-free (VG multiplicity), NOT soft-mask/DNA;
    # the antisense feature set is genome-architecture (coordinates + strand), NOT DNA/soft-mask.
    sys.path.insert(0, BENCH)
    import family_rna_refine as R
    R.rna_only_guard()   # must not raise
    rep = set(R.REPEAT_GATE_FEATURES)
    assert rep == {"min_shared_mult", "cyclic"}, rep
    assert not (rep & (R.DNA_FORBIDDEN | R.LIBRARY_FORBIDDEN)), "soft-mask/DNA leaked into repeat gate"
    anti = set(R.ANTISENSE_GATE_FEATURES)
    assert anti == {"gene_strand", "gene_span"}, anti
    assert not (anti & (R.DNA_FORBIDDEN | R.LIBRARY_FORBIDDEN)), "DNA/soft-mask leaked into antisense gate"
    assert "softmask" in R.LIBRARY_FORBIDDEN and "mask" in R.LIBRARY_FORBIDDEN
    print("(h) RNA-only guard (repeat gate = min_shared_mult/cyclic; antisense gate = gene_strand/gene_span "
          "genome-architecture; no soft-mask/DNA in either) : OK")


def test_high_precision_gamma_040():
    """--high-precision swaps gamma 0.20 -> 0.40 (PRECISION_RECALL_FRONTIER.md): DIFFERENT +
    deterministic catalog, JSON records the active gamma, n_families increases (frontier direction),
    the two collapsed-array OVERSIZE blobs (MPHOSPH8 + LOC134758618) are removed, and the HONEST
    off-oracle KRAB-ZNF + MAGE-floor caveats are carried in the JSON (not dropped).  env == flag."""
    # ISOLATE gamma: hold the multi-repeat-bridge AND antisense gates OFF in all arms so ONLY gamma
    # (0.20->0.40) moves.  With colinear-merge DEFAULT-ON, the reference is 590 families and high-precision
    # is 606 families.  The default (all gates) md5/count is checked in test_default_byte_identical_to_golden;
    # compose with the bridge gate in test_repeat_bridge_composes and with the antisense gate in
    # test_antisense_gate.
    iso = ["--no-repeat-bridge-gate", "--no-antisense-gate"]
    # default reference (repeat-bridge + antisense OFF)
    r = _run(iso);  assert r.returncode == 0, r.stderr
    default_tsv = _md5(OUT_TSV);  d_def = json.load(open(OUT_JSON))
    # --high-precision (flag), repeat-bridge + antisense OFF
    r1 = _run(["--high-precision"] + iso);  assert r1.returncode == 0, r1.stderr
    hp_tsv, hp_json = _md5(OUT_TSV), _md5(OUT_JSON);  d = json.load(open(OUT_JSON))
    # (j) records gamma=0.40 everywhere it reports gamma
    assert d["rule"]["gamma"] == 0.4, f"rule gamma != 0.40: {d['rule']['gamma']}"
    assert d["guards"]["gamma"] == 0.4, f"guards gamma != 0.40: {d['guards']['gamma']}"
    assert d["high_precision"]["active_gamma"] == 0.4, d["high_precision"]["active_gamma"]
    assert d["high_precision"]["default_gamma"] == 0.2, d["high_precision"]["default_gamma"]
    # (j) DIFFERENT catalog vs (repeat-bridge+antisense-OFF) reference (gamma actually applied)
    assert hp_tsv != default_tsv, "high-precision catalog identical to default (gamma not applied)"
    # (j) deterministic across two runs
    r2 = _run(["--high-precision"] + iso);  assert r2.returncode == 0, r2.stderr
    assert _md5(OUT_TSV) == hp_tsv and _md5(OUT_JSON) == hp_json, "high-precision not byte-identical across runs"
    # (j) env form == flag form
    r3 = _run(iso, env_extra={"RUSTLE_HIGH_PRECISION": "1"})
    assert r3.returncode == 0, r3.stderr
    assert _md5(OUT_TSV) == hp_tsv and _md5(OUT_JSON) == hp_json, "RUSTLE_HIGH_PRECISION=1 != --high-precision"
    # (j) precision impact: n_families increases toward the frontier; oversize residual drops (blobs removed)
    assert d["n_families"] == 606, f"high-precision n_families {d['n_families']} != 606"
    assert d["n_families"] > d_def["n_families"], "high-precision did not split toward the frontier"
    rem_hp, rem_def = d["residual_fp"]["residual_remaining"], d_def["residual_fp"]["residual_remaining"]
    assert rem_hp["oversize"] < rem_def["oversize"], \
        f"high-precision did not reduce oversize residual: {rem_def['oversize']} -> {rem_hp['oversize']}"
    assert rem_hp["oversize"] == 1, f"expected oversize residual 1 (MAGE X-array only), got {rem_hp['oversize']}"
    ex_hp = " ".join(d["residual_fp"]["residual_examples"]["oversize"])
    assert "MPHOSPH8" not in ex_hp and "LOC134758618" not in ex_hp, f"OVERSIZE blobs not removed: {ex_hp}"
    assert "LOC129529978" in ex_hp, f"MAGE X-array DNA-only floor missing (should survive): {ex_hp}"
    # (j) HONEST caveats carried, not dropped
    hp = d["high_precision"]
    for key in ("precision_impact", "offoracle_krabznf_cost", "mage_floor", "over_split_cost",
                "frontier_row_gamma040", "live_precision_signal", "source"):
        assert hp.get(key), f"high_precision field missing: {key}"
    assert "KRAB-ZNF" in hp["offoracle_krabznf_cost"] and "MAGE" in hp["mage_floor"], "caveats hollowed out"
    print(f"(j) --high-precision gamma=0.40 (JSON records it); n_families {d_def['n_families']} -> "
          f"{d['n_families']} (frontier direction); oversize residual {rem_def['oversize']} -> {rem_hp['oversize']} "
          f"(MPHOSPH8 + LOC134758618 removed, MAGE floor survives); deterministic; env==flag; caveats carried : OK")
    _run([])   # restore DEFAULT catalog for any downstream test


def test_flags_compose():
    """--legacy still opts out even with --high-precision; --high-precision composes with --no-repeat-gate."""
    # --legacy wins (opt-out) even alongside --high-precision -> nothing written
    _remove_settled(OUT_TSV, OUT_JSON)
    r = _run(["--legacy", "--high-precision"])
    assert r.returncode == 0 and r.stdout.startswith("legacy:"), repr(r.stdout)
    assert not os.path.exists(OUT_TSV) and not os.path.exists(OUT_JSON), "legacy+high-precision wrote outputs"
    # --high-precision composes with --no-repeat-gate: gamma=0.40 AND the repeat-hub gate ablated
    r = _run(["--high-precision", "--no-repeat-gate"]);  assert r.returncode == 0, r.stderr
    d = json.load(open(OUT_JSON))
    assert d["rule"]["gamma"] == 0.4, f"gamma not 0.40 under compose: {d['rule']['gamma']}"
    assert d["rule"]["repeat_gate_enabled"] is False, "--no-repeat-gate not applied under --high-precision"
    assert d["edges"]["n_dn_cross_gene_cut_by_repeat_gate"] == 0, "repeat-hub gate still cut under --no-repeat-gate"
    assert "high_precision" in d, "high_precision block missing when composing with --no-repeat-gate"
    assert d["multi_repeat_bridge"]["enabled"] is True, "multi-repeat-bridge gate lost under compose"
    # SANITY (isolate the repeat-HUB gate with the repeat-bridge gate held OFF): at gamma=0.40 the
    # DEFAULT-ON repeat-bridge gate SUBSUMES the single-edge repeat-hub gate (it removes the same
    # DISCONNECTED families), so --no-repeat-gate is a NET no-op on the FINAL catalog with the bridge
    # gate on -- the repeat-hub gate's independent effect is only visible with --no-repeat-bridge-gate.
    r = _run(["--high-precision", "--no-repeat-gate", "--no-repeat-bridge-gate"])
    assert r.returncode == 0, r.stderr
    hp_norg_nb = _md5(OUT_TSV)
    r = _run(["--high-precision", "--no-repeat-bridge-gate"]);  assert r.returncode == 0, r.stderr
    hp_nb = _md5(OUT_TSV)
    assert hp_norg_nb != hp_nb, "--no-repeat-gate had no effect under --high-precision (bridge held off)"
    print("(k) --legacy opts out even with --high-precision; --high-precision + --no-repeat-gate + "
          "--no-repeat-bridge-gate compose (repeat-bridge subsumes repeat-hub at gamma=0.40) : OK")
    _run([])   # restore DEFAULT catalog


def test_repeat_bridge_gate():
    """(a) --no-repeat-bridge-gate recovers the pre-repeat-bridge golden (0da31d59 / 5f438104) byte-identically
    (flag == env, deterministic); (b) DEFAULT (gate ON) cuts >=1 DISCONNECTED repeat-bridge over-merge
    (named -- fam17 hub + PDIA3+PRKAB2 etc.) and RECORDS the count; (c) GSTM2 (fid 9/13) + MAGE
    (546/548/553) STAY MERGED (largest block byte-identical membership with vs without the gate);
    (g) determinism (default byte-identical across runs)."""
    # (a) ablation (flag; antisense + colinear-merge still DEFAULT-ON) recovers the antisense-ON
    # pre-repeat-bridge catalog byte-identical.
    r = _run(["--no-repeat-bridge-gate"]);  assert r.returncode == 0, r.stderr
    off_tsv = _md5(OUT_TSV);  off = _fam_to_genes(OUT_TSV);  off_json = json.load(open(OUT_JSON))
    assert off_tsv == GOLDEN_ANTISENSE_NO_REPEAT_BRIDGE_TSV_MD5, \
        f"--no-repeat-bridge-gate drifted from the antisense-ON pre-bridge golden " \
        f"{GOLDEN_ANTISENSE_NO_REPEAT_BRIDGE_TSV_MD5}: {off_tsv}"
    assert off_json["multi_repeat_bridge"]["enabled"] is False
    assert off_json["multi_repeat_bridge"]["n_families_cut"] == 0, "ablation still cut families"
    assert off_json["n_families"] == 581
    assert off_json["antisense_overlap_gate"]["enabled"] is True, "antisense gate wrongly disabled"
    # additionally ablating the antisense gate recovers the pre-repeat-bridge golden 5f438104.
    rh = _run(["--no-repeat-bridge-gate", "--no-antisense-gate"]);  assert rh.returncode == 0, rh.stderr
    assert _md5(OUT_TSV) == GOLDEN_NO_REPEAT_BRIDGE_TSV_MD5, \
        f"--no-repeat-bridge-gate --no-antisense-gate != pre-bridge golden {GOLDEN_NO_REPEAT_BRIDGE_TSV_MD5}"
    assert json.load(open(OUT_JSON))["n_families"] == 590
    # (e) env form == flag form (byte-identical)
    r2 = _run([], env_extra={"RUSTLE_NO_REPEAT_BRIDGE_GATE": "1"});  assert r2.returncode == 0, r2.stderr
    assert _md5(OUT_TSV) == off_tsv, "RUSTLE_NO_REPEAT_BRIDGE_GATE=1 != --no-repeat-bridge-gate"
    # DEFAULT (gate ON) -- (g) deterministic across two runs
    r3 = _run([]);  assert r3.returncode == 0, r3.stderr
    on_tsv = _md5(OUT_TSV);  on = _fam_to_genes(OUT_TSV);  on_json = json.load(open(OUT_JSON))
    r4 = _run([]);  assert r4.returncode == 0, r4.stderr
    assert _md5(OUT_TSV) == on_tsv, "default not byte-identical across runs (non-deterministic)"
    assert on_tsv == GOLDEN_DEFAULT_TSV_MD5 and on_tsv != off_tsv, "gate had no effect / drifted from golden"
    # (b) gate cuts >=1 DISCONNECTED and records the count + named list; every cut is DISCONNECTED
    mrb = on_json["multi_repeat_bridge"]
    assert mrb["enabled"] is True
    assert mrb["n_families_cut"] >= 1, "gate cut nothing"
    assert mrb["n_families_cut"] == 59, f"default n_families_cut != 59: {mrb['n_families_cut']}"
    assert mrb["n_disconnected_removed"] == mrb["n_families_cut"], "not every cut is DISCONNECTED"
    named = " ".join(mrb["families_cut_named"])
    # a named repeat-bridge over-merge is split: the fam17 Alu hub (C9H11orf65...) + PDIA3+PRKAB2
    assert "C9H11orf65" in named, "fam17 repeat-bridge hub not among the cut families"
    assert "PDIA3+PRKAB2" in mrb["families_cut_named"], "PDIA3+PRKAB2 over-merge not cut"
    # (c) GSTM2 + MAGE control families STAY MERGED: largest block identical membership ON vs OFF
    for name, gl in (("GSTM2", GSTM2), ("MAGE", MAGE)):
        bo = max(_blocks_with(off, gl), key=len)
        bn = max(_blocks_with(on, gl), key=len)
        assert bo == bn, f"{name} largest block CHANGED by the repeat-bridge gate (not spared): {bo ^ bn}"
    print(f"(m) --no-repeat-bridge-gate recovers 0da31d59 byte-identical (+--no-antisense-gate -> 5f438104); "
          f"DEFAULT gate ON cuts {mrb['n_families_cut']} DISCONNECTED (fam17 + PDIA3+PRKAB2 named); "
          f"GSTM2/MAGE spared (largest block identical); deterministic : OK")


def test_repeat_bridge_r_oracle_held():
    """(d) R_oracle 51/57 = 0.8947 HELD on the NEW default catalog (antisense + colinear-merge ON),
    measured VIA bench/family_level_pr_current.py (the shipped diploid-DNA-oracle recall +
    candidate_generation_recall relabel).  E_p purity is 0.8879; diploid distinct-FP stays 4, oversize 1."""
    r = _run([]);  assert r.returncode == 0, r.stderr          # ensure DEFAULT catalog on disk
    assert _md5(OUT_TSV) == GOLDEN_DEFAULT_TSV_MD5
    pr_json = os.path.join(BENCH, "family_level_pr_current.json")
    pr = subprocess.run([PY, os.path.join(BENCH, "family_level_pr_current.py")], cwd=BENCH,
                        capture_output=True, text=True)
    assert pr.returncode == 0, pr.stderr
    assert _settle_json(pr_json), "family_level_pr_current.json not visible/parseable after run"
    d = json.load(open(pr_json))
    o = d["truth3_diploid_oracle"]["current"]
    # RECALL NEUTRAL: the antisense gate cuts only non-oracle same-locus opposite-strand edges, so the
    # 51/57 (with the candidate_generation_recall relabel) oracle recovery is UNCHANGED vs pre-antisense.
    assert o["oracle_genes_recovered_multicopy"] == 51, \
        f"R_oracle numerator moved off 51 (recall NOT neutral): {o['oracle_genes_recovered_multicopy']}"
    assert o["oracle_multicopy_genes"] == 57, o["oracle_multicopy_genes"]
    assert o["recall_oracle"] == 0.8947, f"R_oracle != 0.8947: {o['recall_oracle']}"
    # E_p purity; diploid distinct-FP stays 4, oversize 1
    assert d["truth1_protein_Ep"]["current"]["precision_Ep"] == 0.8879, \
        d["truth1_protein_Ep"]["current"]["precision_Ep"]
    assert o["distinct_fp_blocks"] == 4, f"distinct-FP not <= 4: {o['distinct_fp_blocks']}"
    assert o["n_oversize"] == 1, f"oversize FP not 1: {o['n_oversize']}"
    print(f"(n) R_oracle {o['oracle_genes_recovered_multicopy']}/{o['oracle_multicopy_genes']} = "
          f"{o['recall_oracle']} HELD (recall NEUTRAL, via family_level_pr_current); P_Ep -> "
          f"{d['truth1_protein_Ep']['current']['precision_Ep']} (distinct-FP "
          f"{o['distinct_fp_blocks']}, oversize {o['n_oversize']}) : OK")


def test_repeat_bridge_composes():
    """(f) the multi-repeat-bridge gate COMPOSES with the other ablation flags: with --no-repeat-gate,
    --no-split-recombinants and --high-precision the gate is still DEFAULT-ON (cuts families); adding
    --no-repeat-bridge-gate to each disables ONLY it (recovers the matching pre-repeat-bridge catalog)."""
    # compose with --no-repeat-gate (repeat-hub OFF, bridge ON)
    r = _run(["--no-repeat-gate"]);  assert r.returncode == 0, r.stderr
    d = json.load(open(OUT_JSON))
    assert d["rule"]["repeat_gate_enabled"] is False and d["multi_repeat_bridge"]["enabled"] is True
    assert d["multi_repeat_bridge"]["n_families_cut"] >= 1, "bridge gate inert under --no-repeat-gate"
    # compose with --no-split-recombinants (split OFF, bridge + antisense ON); +--no-repeat-bridge-gate
    # --no-antisense-gate recovers the 591 pre-split golden (8f79718a).
    r = _run(["--no-split-recombinants"]);  assert r.returncode == 0, r.stderr
    d = json.load(open(OUT_JSON))
    assert d["recombinant_split"]["enabled"] is False and d["multi_repeat_bridge"]["enabled"] is True
    assert d["antisense_overlap_gate"]["enabled"] is True
    assert d["multi_repeat_bridge"]["n_families_cut"] >= 1, "bridge gate inert under --no-split-recombinants"
    on_nosplit = _md5(OUT_TSV)
    r = _run(["--no-split-recombinants", "--no-repeat-bridge-gate", "--no-antisense-gate"])
    assert r.returncode == 0, r.stderr
    assert _md5(OUT_TSV) == GOLDEN_NOSPLIT_TSV_MD5 and _md5(OUT_TSV) != on_nosplit, \
        "--no-split-recombinants --no-repeat-bridge-gate --no-antisense-gate did not recover the 591 pre-split golden"
    # compose with --high-precision (gamma 0.40, bridge ON)
    r = _run(["--high-precision"]);  assert r.returncode == 0, r.stderr
    d = json.load(open(OUT_JSON))
    assert d["rule"]["gamma"] == 0.4 and d["multi_repeat_bridge"]["enabled"] is True
    assert d["multi_repeat_bridge"]["n_families_cut"] >= 1, "bridge gate inert under --high-precision"
    assert "high_precision" in d, "high_precision block missing when composing with the bridge gate"
    hp_bridge = _md5(OUT_TSV)
    r = _run(["--high-precision", "--no-repeat-bridge-gate"]);  assert r.returncode == 0, r.stderr
    assert _md5(OUT_TSV) != hp_bridge, "--no-repeat-bridge-gate had no effect under --high-precision"
    print("(o) multi-repeat-bridge composes with --no-repeat-gate / --no-split-recombinants / "
          "--high-precision (default-ON in each; --no-repeat-bridge-gate disables only it) : OK")
    _run([])   # restore DEFAULT catalog


if __name__ == "__main__":
    test_legacy_opt_out_writes_nothing()
    test_default_deterministic()
    test_default_byte_identical_to_golden()
    test_antisense_gate()
    test_colinear_merge_gate()
    test_allele_demote_removes_known_fp()
    test_residual_removed_matches_recall_preserving_row()
    test_recombinant_split_gate()
    test_repeat_gate_shatters_fam17_spares_controls()
    test_repeat_gate_recall_cost_is_genuine_only()
    test_rna_only_guard()
    test_repeat_bridge_gate()
    test_repeat_bridge_r_oracle_held()
    test_repeat_bridge_composes()
    test_high_precision_gamma_040()
    test_flags_compose()
    _run([])   # leave the DEFAULT (all-gates) catalog on disk
    print("\nALL TESTS PASSED")
