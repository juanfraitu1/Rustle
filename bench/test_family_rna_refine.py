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
  (i) DEFAULT (no flag) catalog is BYTE-IDENTICAL to the pre-change golden TSV (md5) and its JSON
      records gamma=0.20 with NO high_precision block.
  (j) --high-precision (== env RUSTLE_HIGH_PRECISION=1) swaps gamma 0.20 -> 0.40
      (PRECISION_RECALL_FRONTIER.md): JSON records the active gamma, the catalog is DIFFERENT +
      deterministic, n_families -> 623 (frontier gamma=0.40 row), the collapsed-array OVERSIZE blobs
      MPHOSPH8 + LOC134758618 are removed (oversize residual 3 -> 1), and the HONEST off-oracle
      KRAB-ZNF + MAGE-floor caveats are carried in the JSON (not dropped).
  (k) --legacy still opts out even with --high-precision; --high-precision composes with
      --no-repeat-gate (gamma=0.40 AND the gate ablated).

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

# md5 of the DEFAULT catalog (gamma=0.20, repeat-hub + recombinant-split + multi-repeat-bridge ALL ON,
# 573 families).  The --high-precision flag must NOT perturb the default: `family_rna_refine.py` with no
# flag reproduces this.  The multi-repeat-bridge gate (3rd VG-native gate, T=8/C=2) is DEFAULT-ON.
GOLDEN_DEFAULT_TSV_MD5 = "dca64cbdc1de0617c130253221a0b6bd"
# md5 of the pre-repeat-bridge catalog (--no-repeat-bridge-gate: gamma=0.20, repeat-hub + recombinant-
# split ON, multi-repeat-bridge OFF, 605 families) = the CURRENT shipped default BEFORE this gate.
# The multi-repeat-bridge gate is DEFAULT-ON; ablating it must recover this exactly (byte-for-byte).
GOLDEN_NO_REPEAT_BRIDGE_TSV_MD5 = "5e58378ad3920b94c0048d383dfada29"
# md5 of the pre-split catalog (--no-split-recombinants --no-repeat-bridge-gate, gamma=0.20, repeat-hub
# ON, 606 families).  The recombinant-split + multi-repeat-bridge gates are both DEFAULT-ON; ablating
# BOTH recovers this exactly (byte-for-byte).
GOLDEN_NOSPLIT_TSV_MD5 = "f94f387e4f3a53a69e9d8a2d0b1f497a"

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
    r = _run([])            # DEFAULT: no flag (repeat-hub + recombinant-split + multi-repeat-bridge ALL ON)
    assert r.returncode == 0, r.stderr
    got = _md5(OUT_TSV)
    assert got == GOLDEN_DEFAULT_TSV_MD5, \
        f"default catalog drifted from the golden {GOLDEN_DEFAULT_TSV_MD5}: {got}"
    d = json.load(open(OUT_JSON))
    assert d["rule"]["gamma"] == 0.2, f"default JSON gamma != 0.20: {d['rule']['gamma']}"
    assert d["guards"]["gamma"] == 0.2, f"default guards gamma != 0.20: {d['guards']['gamma']}"
    assert d["n_families"] == 573, f"default n_families != 573: {d['n_families']}"
    assert "high_precision" not in d, "default JSON must NOT carry the high_precision block"
    # recombinant-split gate is DEFAULT-ON and splits the 2 recall-safe HIGH-confidence mosaics
    rs = d["recombinant_split"]
    assert rs["enabled"] is True and rs["n_families_split"] == 2, \
        f"default recombinant-split not ON/2: enabled={rs['enabled']} n_split={rs['n_families_split']}"
    assert d["rule"]["split_recombinants_enabled"] is True
    # multi-repeat-bridge gate (3rd VG-native gate) is DEFAULT-ON at T=8/C=2 and cuts the DISCONNECTED
    # repeat-bridge FP class (all cut families are no-full-shared-exon by construction)
    mrb = d["multi_repeat_bridge"]
    assert mrb["enabled"] is True, "default multi-repeat-bridge gate not enabled"
    assert d["rule"]["repeat_bridge_gate_enabled"] is True
    assert (mrb["mult_min"], mrb["count_min"], mrb["exon_id"]) == (8, 2, 0.70), \
        f"repeat-bridge T/C/exon_id drifted: {mrb['mult_min']}/{mrb['count_min']}/{mrb['exon_id']}"
    assert mrb["n_families_cut"] == 61, f"default repeat-bridge n_families_cut != 61: {mrb['n_families_cut']}"
    assert mrb["n_disconnected_removed"] == mrb["n_families_cut"], "disconnected_removed != n_families_cut"
    assert len(mrb["families_cut_named"]) == mrb["n_families_cut"]
    print(f"(i) default (no flag) byte-identical to golden TSV (md5 {GOLDEN_DEFAULT_TSV_MD5}); "
          f"gamma=0.20, 573 families, recombinant-split ON (2 split), multi-repeat-bridge ON "
          f"({mrb['n_families_cut']} cut, T=8/C=2), no high_precision block : OK")


def test_recombinant_split_gate():
    """recombinant-split DEFAULT-ON breaks the fid-210 mosaic (GALNT17 + LOC101126070 no longer
    co-membered; GALNT17's single-locus remnant drops via the >=2-loci multi-copy filter, LOC101126070
    survives); --no-split-recombinants (+ --no-repeat-bridge-gate to isolate the 3rd gate) recovers the
    606-family pre-split golden byte-identically; env==flag.  Byte-identity to the pre-split golden IS
    the recall-safety proof (the split never separates same-gene loci)."""
    # PRE-SPLIT reference (split + bridge both OFF): fid 210 co-members GALNT17 + LOC101126070
    rp = _run(["--no-split-recombinants", "--no-repeat-bridge-gate"]);  assert rp.returncode == 0, rp.stderr
    pre = _fam_to_genes(OUT_TSV)
    assert any({"GALNT17", "LOC101126070"} <= gs for gs in pre.values()), \
        "pre-split: GALNT17 + LOC101126070 not co-membered (expected the fid-210 mosaic)"
    off = _md5(OUT_TSV);  d_off = json.load(open(OUT_JSON))
    assert off == GOLDEN_NOSPLIT_TSV_MD5, f"--no-split-recombinants drifted from pre-split golden: {off}"
    assert d_off["n_families"] == 606 and d_off["recombinant_split"]["enabled"] is False
    assert d_off["multi_repeat_bridge"]["enabled"] is False
    # env form == flag form (byte-identical) -- proves same-gene loci never separated (recall-safe)
    r2 = _run(["--no-repeat-bridge-gate"], env_extra={"RUSTLE_NO_SPLIT_RECOMBINANTS": "1"})
    assert r2.returncode == 0 and _md5(OUT_TSV) == off, "RUSTLE_NO_SPLIT_RECOMBINANTS=1 != flag form"
    # DEFAULT (split ON): the fid-210 mosaic is BROKEN -> GALNT17 + LOC101126070 no longer co-membered
    r = _run([]);  assert r.returncode == 0, r.stderr
    fam = _fam_to_genes(OUT_TSV)
    assert not any({"GALNT17", "LOC101126070"} <= gs for gs in fam.values()), \
        "recombinant mosaic fid 210 NOT split (GALNT17 + LOC101126070 still co-membered)"
    d = json.load(open(OUT_JSON))
    assert d["recombinant_split"]["enabled"] is True and d["recombinant_split"]["n_families_split"] == 2
    print("(l) recombinant-split ON breaks the fid-210 mosaic (GALNT17|LOC101126070 no longer "
          "co-membered); --no-split-recombinants (+--no-repeat-bridge-gate) recovers 606 pre-split "
          "golden; env==flag : OK")
    _run([])   # restore DEFAULT


def test_allele_demote_removes_known_fp():
    d = json.load(open(OUT_JSON))
    demoted = {x["gene"] for x in d["alleles_demoted"]}
    conf = {x["gene"] for x in d["alleles_demoted"] if x["dna_confirmed"]}
    assert "DHRSX" in demoted and "LOC129530050" in demoted, f"missing allele demotions: {demoted}"
    assert {"DHRSX", "LOC129530050"} <= conf, f"DHRSX/LOC129530050 not flagged DNA-confirmed: {conf}"
    print(f"(c) allele demotion removes DHRSX + LOC129530050 : OK  (demoted={sorted(demoted)})")


def test_residual_removed_matches_recall_preserving_row():
    # the recall-preserving row is the PRE-repeat-bridge catalog (--no-repeat-bridge-gate), which the
    # multi-repeat-bridge gate does not touch: residual {allele 0, oversize 3, multifam 4}, 6/12 removed.
    r = _run(["--no-repeat-bridge-gate"]);  assert r.returncode == 0, r.stderr
    d = json.load(open(OUT_JSON))
    rem = d["residual_fp"]["residual_remaining"]
    assert rem == {"allele": 0, "oversize": 3, "multifam": 4}, f"residual remaining drift: {rem}"
    assert d["residual_fp"]["shipped_residual_total"] == 12, d["residual_fp"]["shipped_residual_total"]
    assert d["residual_fp"]["named_removed"] == 6, d["residual_fp"]["named_removed"]
    bd = d["residual_fp"]["named_removed_breakdown"]
    assert bd == {"allele": 2, "oversize": 2, "multifam": 2}, f"removed breakdown drift: {bd}"
    print("(d) pre-repeat-bridge residual {allele 0, oversize 3, multifam 4}; 6/12 removed : OK")
    # the DEFAULT (multi-repeat-bridge ON) removes 2 MORE oversize blobs (MPHOSPH8-hub + LOC134758618)
    r2 = _run([]);  assert r2.returncode == 0, r2.stderr
    d2 = json.load(open(OUT_JSON))
    rem2 = d2["residual_fp"]["residual_remaining"]
    assert rem2 == {"allele": 0, "oversize": 1, "multifam": 4}, f"default residual drift: {rem2}"
    assert d2["residual_fp"]["named_removed"] == 8, d2["residual_fp"]["named_removed"]
    print("(d') DEFAULT (multi-repeat-bridge ON) residual {allele 0, oversize 1, multifam 4}; 8/12 "
          "removed (MPHOSPH8-hub + LOC134758618 oversize blobs cut as DISCONNECTED) : OK")


def test_repeat_gate_shatters_fam17_spares_controls():
    # isolate the REPEAT-HUB gate: hold the recombinant-split AND multi-repeat-bridge gates OFF in both
    # arms so the only difference is the repeat-hub gate (each other gate is exercised separately).
    r = _run(["--no-split-recombinants", "--no-repeat-bridge-gate"])   # repeat-hub ON, split+bridge OFF
    assert r.returncode == 0, r.stderr
    on = _fam_to_genes(OUT_TSV);  on_json = json.load(open(OUT_JSON))
    on_tsv = _md5(OUT_TSV)
    # ablation = repeat-hub gate OFF (flag form) -- run twice for byte-identical determinism
    r2 = _run(["--no-repeat-gate", "--no-split-recombinants", "--no-repeat-bridge-gate"])
    assert r2.returncode == 0, r2.stderr
    off = _fam_to_genes(OUT_TSV);  off_json = json.load(open(OUT_JSON))
    off_tsv = _md5(OUT_TSV)
    r3 = _run(["--no-split-recombinants", "--no-repeat-bridge-gate"],
              env_extra={"RUSTLE_NO_REPEAT_GATE": "1"})   # env form == flag form
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
    assert g["demote_features"] == ["balanced_frac", "copy_like"], g["demote_features"]
    assert g["no_dna_in_inference"] is True and g["gamma"] == 0.2 and g["seed"] == 0
    assert g["repeat_gate_library_free"] is True and g["no_softmask_in_repeat_gate"] is True
    # source-level: the repeat-hub feature set is library-free (VG multiplicity), NOT soft-mask/DNA
    sys.path.insert(0, BENCH)
    import family_rna_refine as R
    R.rna_only_guard()   # must not raise
    rep = set(R.REPEAT_GATE_FEATURES)
    assert rep == {"min_shared_mult", "cyclic"}, rep
    assert not (rep & (R.DNA_FORBIDDEN | R.LIBRARY_FORBIDDEN)), "soft-mask/DNA leaked into repeat gate"
    assert "softmask" in R.LIBRARY_FORBIDDEN and "mask" in R.LIBRARY_FORBIDDEN
    print("(h) RNA-only / library-free guard (repeat gate = min_shared_mult/cyclic; no soft-mask/DNA) : OK")


def test_high_precision_gamma_040():
    """--high-precision swaps gamma 0.20 -> 0.40 (PRECISION_RECALL_FRONTIER.md): DIFFERENT +
    deterministic catalog, JSON records the active gamma, n_families -> 623 (frontier row), the two
    collapsed-array OVERSIZE blobs (MPHOSPH8 + LOC134758618) are removed, and the HONEST off-oracle
    KRAB-ZNF + MAGE-floor caveats are carried in the JSON (not dropped).  env == flag."""
    # ISOLATE gamma: hold the multi-repeat-bridge gate OFF in all arms so ONLY gamma (0.20->0.40) moves,
    # preserving the frontier semantics (605->622, oversize 3->1).  The default (all gates) md5/count is
    # checked in test_default_byte_identical_to_golden; compose with the bridge gate in test_repeat_bridge_composes.
    # default reference (repeat-bridge OFF)
    r = _run(["--no-repeat-bridge-gate"]);  assert r.returncode == 0, r.stderr
    default_tsv = _md5(OUT_TSV);  d_def = json.load(open(OUT_JSON))
    # --high-precision (flag), repeat-bridge OFF
    r1 = _run(["--high-precision", "--no-repeat-bridge-gate"]);  assert r1.returncode == 0, r1.stderr
    hp_tsv, hp_json = _md5(OUT_TSV), _md5(OUT_JSON);  d = json.load(open(OUT_JSON))
    # (j) records gamma=0.40 everywhere it reports gamma
    assert d["rule"]["gamma"] == 0.4, f"rule gamma != 0.40: {d['rule']['gamma']}"
    assert d["guards"]["gamma"] == 0.4, f"guards gamma != 0.40: {d['guards']['gamma']}"
    assert d["high_precision"]["active_gamma"] == 0.4, d["high_precision"]["active_gamma"]
    assert d["high_precision"]["default_gamma"] == 0.2, d["high_precision"]["default_gamma"]
    # (j) DIFFERENT catalog vs (repeat-bridge-OFF) reference (gamma actually applied)
    assert hp_tsv != default_tsv, "high-precision catalog identical to default (gamma not applied)"
    # (j) deterministic across two runs
    r2 = _run(["--high-precision", "--no-repeat-bridge-gate"]);  assert r2.returncode == 0, r2.stderr
    assert _md5(OUT_TSV) == hp_tsv and _md5(OUT_JSON) == hp_json, "high-precision not byte-identical across runs"
    # (j) env form == flag form
    r3 = _run(["--no-repeat-bridge-gate"], env_extra={"RUSTLE_HIGH_PRECISION": "1"})
    assert r3.returncode == 0, r3.stderr
    assert _md5(OUT_TSV) == hp_tsv and _md5(OUT_JSON) == hp_json, "RUSTLE_HIGH_PRECISION=1 != --high-precision"
    # (j) precision impact: n_families toward the frontier (623 pre-split; 622 with the default-on
    # recombinant-split gate, which removes one net family), oversize residual drops (blobs removed)
    assert d["n_families"] == 622, f"high-precision n_families {d['n_families']} != 622 (frontier 623 - 1 recombinant split)"
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
          f"{d['n_families']} (frontier 623 - 1 recombinant split); oversize residual {rem_def['oversize']} -> {rem_hp['oversize']} "
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
    """(a) --no-repeat-bridge-gate recovers the pre-repeat-bridge golden (5e58378a) byte-identically
    (flag == env, deterministic); (b) DEFAULT (gate ON) cuts >=1 DISCONNECTED repeat-bridge over-merge
    (named -- fam17 hub + PDIA3+PRKAB2 etc.) and RECORDS the count; (c) GSTM2 (fid 9/13) + MAGE
    (546/548/553) STAY MERGED (largest block byte-identical membership with vs without the gate);
    (g) determinism (default byte-identical across runs)."""
    # (a) ablation (flag) recovers the pre-repeat-bridge golden, byte-identical
    r = _run(["--no-repeat-bridge-gate"]);  assert r.returncode == 0, r.stderr
    off_tsv = _md5(OUT_TSV);  off = _fam_to_genes(OUT_TSV);  off_json = json.load(open(OUT_JSON))
    assert off_tsv == GOLDEN_NO_REPEAT_BRIDGE_TSV_MD5, \
        f"--no-repeat-bridge-gate drifted from the pre-repeat-bridge golden {GOLDEN_NO_REPEAT_BRIDGE_TSV_MD5}: {off_tsv}"
    assert off_json["multi_repeat_bridge"]["enabled"] is False
    assert off_json["multi_repeat_bridge"]["n_families_cut"] == 0, "ablation still cut families"
    assert off_json["n_families"] == 605
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
    assert mrb["n_families_cut"] == 61, f"default n_families_cut != 61: {mrb['n_families_cut']}"
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
    print(f"(m) --no-repeat-bridge-gate recovers 5e58378a byte-identical (flag==env); DEFAULT gate ON "
          f"cuts {mrb['n_families_cut']} DISCONNECTED (fam17 + PDIA3+PRKAB2 named); GSTM2/MAGE spared "
          f"(largest block identical); deterministic : OK")


def test_repeat_bridge_r_oracle_held():
    """(d) R_oracle 50/57 HELD on the NEW default catalog, measured VIA bench/family_level_pr_current.py
    (the shipped diploid-DNA-oracle recall) -- and E_p purity IMPROVES 0.8694 -> 0.8918."""
    r = _run([]);  assert r.returncode == 0, r.stderr          # ensure DEFAULT catalog on disk
    assert _md5(OUT_TSV) == GOLDEN_DEFAULT_TSV_MD5
    pr_json = os.path.join(BENCH, "family_level_pr_current.json")
    pr = subprocess.run([PY, os.path.join(BENCH, "family_level_pr_current.py")], cwd=BENCH,
                        capture_output=True, text=True)
    assert pr.returncode == 0, pr.stderr
    assert _settle_json(pr_json), "family_level_pr_current.json not visible/parseable after run"
    d = json.load(open(pr_json))
    o = d["truth3_diploid_oracle"]["current"]
    assert o["oracle_genes_recovered_multicopy"] == 50, \
        f"R_oracle numerator dropped below 50: {o['oracle_genes_recovered_multicopy']}"
    assert o["oracle_multicopy_genes"] == 57, o["oracle_multicopy_genes"]
    assert o["recall_oracle"] == 0.8772, f"R_oracle != 0.8772: {o['recall_oracle']}"
    # E_p purity improves (diploid distinct-FP shrinks 6 -> 4; oversize 3 -> 1)
    assert d["truth1_protein_Ep"]["current"]["precision_Ep"] == 0.8918, \
        d["truth1_protein_Ep"]["current"]["precision_Ep"]
    assert o["distinct_fp_blocks"] == 4, f"distinct-FP not reduced to 4: {o['distinct_fp_blocks']}"
    assert o["n_oversize"] == 1, f"oversize FP not reduced to 1: {o['n_oversize']}"
    print(f"(n) R_oracle {o['oracle_genes_recovered_multicopy']}/{o['oracle_multicopy_genes']} = "
          f"{o['recall_oracle']} HELD (via family_level_pr_current); P_Ep -> "
          f"{d['truth1_protein_Ep']['current']['precision_Ep']} (distinct-FP {o['distinct_fp_blocks']}, "
          f"oversize {o['n_oversize']}) : OK")


def test_repeat_bridge_composes():
    """(f) the multi-repeat-bridge gate COMPOSES with the other ablation flags: with --no-repeat-gate,
    --no-split-recombinants and --high-precision the gate is still DEFAULT-ON (cuts families); adding
    --no-repeat-bridge-gate to each disables ONLY it (recovers the matching pre-repeat-bridge catalog)."""
    # compose with --no-repeat-gate (repeat-hub OFF, bridge ON)
    r = _run(["--no-repeat-gate"]);  assert r.returncode == 0, r.stderr
    d = json.load(open(OUT_JSON))
    assert d["rule"]["repeat_gate_enabled"] is False and d["multi_repeat_bridge"]["enabled"] is True
    assert d["multi_repeat_bridge"]["n_families_cut"] >= 1, "bridge gate inert under --no-repeat-gate"
    # compose with --no-split-recombinants (split OFF, bridge ON) -> 574 families; +--no-repeat-bridge-gate -> 606
    r = _run(["--no-split-recombinants"]);  assert r.returncode == 0, r.stderr
    d = json.load(open(OUT_JSON))
    assert d["recombinant_split"]["enabled"] is False and d["multi_repeat_bridge"]["enabled"] is True
    assert d["multi_repeat_bridge"]["n_families_cut"] >= 1, "bridge gate inert under --no-split-recombinants"
    on574 = _md5(OUT_TSV)
    r = _run(["--no-split-recombinants", "--no-repeat-bridge-gate"]);  assert r.returncode == 0, r.stderr
    assert _md5(OUT_TSV) == GOLDEN_NOSPLIT_TSV_MD5 and _md5(OUT_TSV) != on574, \
        "--no-split-recombinants --no-repeat-bridge-gate did not recover the 606 pre-split golden"
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
