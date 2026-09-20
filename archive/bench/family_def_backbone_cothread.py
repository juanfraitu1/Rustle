#!/usr/bin/env python3
"""family_def_backbone_cothread.py — A3: is the VARIATION GRAPH (its junction
backbone) the right STRUCTURE for family DEFINITION, not just copy resolution?

HYPOTHESIS (structural, orthogonal to cDNA identity):
  In a REAL whole-gene-duplication family, every copy is a parallel PATH through
  ONE shared bubble-chain VG: they CO-THREAD a homologous splice-junction backbone
  (same number of introns, same exon-length signature, same phase). A
  domain-sharer shares at most a SUB-PATH (one bubble / a partial backbone). A
  RETROCOPY is INTRONLESS (0 junctions) so it cannot thread the backbone at all,
  even though its cDNA id ~ 1.0.

So the family test is: "do the candidate members co-thread a shared junction
backbone (parallel paths in one bubble chain)?" -- a structural criterion the
de-tie read-conflict graph and ~B coverage do not directly encode.

WHAT THIS SCRIPT MEASURES (annotation-derived, no `rustle --vg`, no BAM):
  For each panel case, take each member's REPRESENTATIVE transcript (the mRNA with
  the most exons; ties -> longest span). Its junction backbone = ordered tuple of
  (intron, then the exon-length signature). We then build the family VG SPINE by
  aligning member backbones to the member with the most introns (the reference
  spine) and asking what fraction of the reference's introns is CO-THREADED by all
  members within tolerance (matching exon-length signature => homologous junction).

  shared_backbone_fraction(case)
     = (# reference introns co-threaded by ALL members) / (# reference introns)
  Also: per-member co-thread fraction, intron counts, and the bridge test (does any
  member thread <50% of the backbone => it is a sub-path sharer, not a family member).

This is a MODEL-level (not read-level) test, complementary to the existing
family_def_junction_concordance.py which threads READS.

Run: /home/juanfra/miniforge3/bin/python bench/family_def_backbone_cothread.py
Output: bench/family_def_backbone_cothread.json
"""
import collections
import json
import os
import subprocess

HERE = os.path.dirname(os.path.abspath(__file__))
GFF = "/mnt/c/Users/jfris/Desktop/GGO_genomic.gff"
# Pre-extracted caches (built once with one GFF pass; see docstring at bottom).
REGION_CACHE = "/tmp/panel_gene_regions.txt"     # name<TAB>chrom<TAB>start<TAB>end
FEATURE_CACHE = "/tmp/panel_features.gff"        # mRNA/exon lines in panel windows

# Panel: case -> (kind, [member gene Names])   (LOC names = retrocopies)
PANEL = [
    ("RABL2",            "real",       ["RABL2A", "RABL2B"]),
    ("APOBEC3",          "real",       ["APOBEC3D", "APOBEC3F"]),
    ("RFPL",             "real",       ["RFPL1", "RFPL2", "RFPL3"]),
    ("CASP8~FLACC1",     "domain",     ["CASP8", "FLACC1"]),
    ("ASDURF~ASNSD1",    "domain",     ["ASDURF", "ASNSD1"]),
    ("CREB1~METTL21A",   "domain",     ["CREB1", "METTL21A"]),
    ("GPR39~LYPD1",      "domain",     ["GPR39", "LYPD1"]),
    ("EEF1A1_retro",     "retro",      ["EEF1A1", "LOC109023808"]),
    ("CNN2_retro",       "retro",      ["CNN2", "LOC129524764"]),
]

ABS_TOL = 25      # bp tolerance on exon-block length match
REL_TOL = 0.15


_REGIONS = None
_FEATURES = None  # parsed once: mrna_parent + exons-by-mrna


def _load_caches():
    """One-time load of the small pre-extracted panel caches. If absent, rebuild
    them with a single GFF pass (slow ~1 min, but only once)."""
    global _REGIONS, _FEATURES
    if _REGIONS is not None:
        return
    if not (os.path.exists(REGION_CACHE) and os.path.exists(FEATURE_CACHE)):
        _build_caches()
    _REGIONS = {}
    for line in open(REGION_CACHE):
        name, chrom, s, e = line.rstrip("\n").split("\t")
        _REGIONS[name] = (chrom, int(s), int(e))
    mrna_parent = {}
    exons = collections.defaultdict(list)
    for line in open(FEATURE_CACHE):
        f = line.rstrip("\n").split("\t")
        typ, s, e, attr = f[2], int(f[3]), int(f[4]), f[8]
        if typ == "mRNA":
            mrna_parent[_attr(attr, "ID")] = _attr(attr, "Parent")
        elif typ == "exon":
            exons[_attr(attr, "Parent")].append((s, e))
    _FEATURES = (mrna_parent, exons)


def _build_caches():
    names = sorted({g for _, _, ms in PANEL for g in ms})
    pat = "Name=(" + "|".join(names) + ");"
    gl = subprocess.run(["grep", "-P", r"\tgene\t", GFF], capture_output=True, text=True).stdout
    with open(REGION_CACHE, "w") as o:
        for line in gl.splitlines():
            f = line.split("\t")
            if f[2] != "gene":
                continue
            nm = _attr(f[8], "Name")
            if nm in names:
                o.write(f"{nm}\t{f[0]}\t{f[3]}\t{f[4]}\n")
    # one awk pass for mRNA/exon features in panel windows
    awk = (
        'NR==FNR { name[FNR]=$1; chr[FNR]=$2; s[FNR]=$3; e[FNR]=$4; n=FNR; next }'
        '($3=="mRNA"||$3=="exon"){for(i=1;i<=n;i++){'
        'if($1==chr[i]&&$4>=s[i]-5&&$5<=e[i]+5){print;break}}}')
    with open(FEATURE_CACHE, "w") as o:
        subprocess.run(["awk", "-F", "\t", awk, REGION_CACHE, GFF], stdout=o)


def gene_region(name):
    _load_caches()
    return _REGIONS.get(name)


def transcripts_for_gene(name):
    """Return ({mRNA_id: sorted [(exon_start, exon_end), ...]}, region) for gene."""
    _load_caches()
    reg = _REGIONS.get(name)
    if reg is None:
        return {}, None
    mrna_parent, exons = _FEATURES
    want_gene = f"gene-{name}"
    tx = {}
    for mid, exs in exons.items():
        if mrna_parent.get(mid) == want_gene:
            tx[mid] = sorted(exs)
    return tx, reg


def _attr(attr, key):
    for kv in attr.split(";"):
        if kv.startswith(key + "="):
            return kv[len(key) + 1:]
    return None


def representative(tx):
    """Pick the mRNA with the most exons; tie -> longest span."""
    best = None
    for mid, exs in tx.items():
        span = exs[-1][1] - exs[0][0]
        key = (len(exs), span)
        if best is None or key > best[0]:
            best = (key, mid, exs)
    return best[2] if best else None


def exon_lengths(exs):
    return [b - a + 1 for a, b in exs]


def intron_lengths(exs):
    return [exs[i + 1][0] - exs[i][1] - 1 for i in range(len(exs) - 1)]


def close(x, y):
    return abs(x - y) <= max(ABS_TOL, REL_TOL * max(x, y))


def best_offset_match(ref_exlen, mem_exlen):
    """Align member exon-length chain to reference chain allowing a start OFFSET
    (a copy may be a contiguous sub-chain of the spine). Return (best_offset,
    matched_intron_count) where matched_intron_count counts consecutive exon-pair
    matches => co-threaded introns. We slide the member over the reference and take
    the offset with the most matched introns."""
    nref, nmem = len(ref_exlen), len(mem_exlen)
    best = (0, 0)
    # member must overlap; slide so member[0] aligns to ref[off]
    for off in range(-(nmem - 1), nref):
        matched_exons = 0
        # count matched consecutive exon-block lengths
        run = 0
        best_run = 0
        for j in range(nmem):
            i = off + j
            if 0 <= i < nref and close(ref_exlen[i], mem_exlen[j]):
                run += 1
                best_run = max(best_run, run)
            else:
                run = 0
        # introns co-threaded = (consecutive matched exon blocks - 1), at least 0
        introns = max(0, best_run - 1)
        if introns > best[1]:
            best = (off, introns)
    return best


def run_case(case, kind, members):
    recs = []
    for m in members:
        tx, reg = transcripts_for_gene(m)
        rep = representative(tx) if tx else None
        if rep is None:
            recs.append(dict(gene=m, exons=0, introns=0, exlen=[], inlen=[], reg=reg))
            continue
        exl = exon_lengths(rep)
        inl = intron_lengths(rep)
        recs.append(dict(gene=m, exons=len(rep), introns=len(inl),
                         exlen=exl, inlen=inl, reg=reg,
                         span=rep[-1][1] - rep[0][0]))
    # reference spine = member with most introns
    ref = max(recs, key=lambda r: r["introns"])
    ref_intr = ref["introns"]
    per_member = []
    co_all = ref_intr  # introns co-threaded by ALL members (start = ref's own)
    for r in recs:
        if r is ref:
            per_member.append(dict(gene=r["gene"], introns=r["introns"],
                                   co_threaded=ref_intr,
                                   frac=1.0 if ref_intr else 0.0, is_ref=True))
            continue
        if r["introns"] == 0 or ref_intr == 0:
            ct = 0
        else:
            _, ct = best_offset_match(ref["exlen"], r["exlen"])
        per_member.append(dict(gene=r["gene"], introns=r["introns"], co_threaded=ct,
                               frac=ct / ref_intr if ref_intr else 0.0, is_ref=False))
        co_all = min(co_all, ct)
    shared_frac = co_all / ref_intr if ref_intr else 0.0
    # bridge test: any member threading < 50% of backbone => sub-path sharer
    min_member_frac = min(p["frac"] for p in per_member)
    any_intronless = any(r["introns"] == 0 for r in recs)
    return dict(case=case, kind=kind, ref_gene=ref["gene"], ref_introns=ref_intr,
                shared_backbone_introns=co_all, shared_backbone_fraction=round(shared_frac, 3),
                min_member_fraction=round(min_member_frac, 3),
                any_intronless_member=any_intronless,
                members=per_member,
                member_intron_counts={r["gene"]: r["introns"] for r in recs})


def main():
    results = [run_case(*c) for c in PANEL]
    print("=== A3: family-VG junction-BACKBONE co-threading ===\n")
    print(f"{'case':18} {'kind':7} {'refIntr':>7} {'shared':>7} {'minMemb':>7} "
          f"{'intronless?':>11}  members(introns,co-thread,frac)")
    for r in results:
        ms = "  ".join(f"{p['gene']}({p['introns']},{p['co_threaded']},{p['frac']:.2f})"
                       for p in r["members"])
        print(f"{r['case']:18} {r['kind']:7} {r['ref_introns']:>7} "
              f"{r['shared_backbone_fraction']:>7.2f} {r['min_member_fraction']:>7.2f} "
              f"{str(r['any_intronless_member']):>11}  {ms}")

    # separation summary
    real = [r["shared_backbone_fraction"] for r in results if r["kind"] == "real"]
    domain = [r["shared_backbone_fraction"] for r in results if r["kind"] == "domain"]
    retro = [r["shared_backbone_fraction"] for r in results if r["kind"] == "retro"]
    print("\n--- shared_backbone_fraction by class ---")
    print(f"  REAL    : {sorted(real)}   range [{min(real):.2f},{max(real):.2f}]")
    print(f"  DOMAIN  : {sorted(domain)}  range [{min(domain):.2f},{max(domain):.2f}]")
    print(f"  RETRO   : {sorted(retro)}   range [{min(retro):.2f},{max(retro):.2f}]")
    sep_real_domain = min(real) - max(domain)
    sep_real_retro = min(real) - max(retro)
    print(f"\n  margin REAL_min - DOMAIN_max = {sep_real_domain:+.3f}  (>0 => clean separation)")
    print(f"  margin REAL_min - RETRO_max  = {sep_real_retro:+.3f}")

    json.dump(dict(panel=results,
                   summary=dict(real=sorted(real), domain=sorted(domain), retro=sorted(retro),
                                margin_real_minus_domain=round(sep_real_domain, 3),
                                margin_real_minus_retro=round(sep_real_retro, 3))),
              open(os.path.join(HERE, "family_def_backbone_cothread.json"), "w"), indent=2)
    print(f"\n[+] wrote {os.path.join(HERE, 'family_def_backbone_cothread.json')}")


if __name__ == "__main__":
    main()
