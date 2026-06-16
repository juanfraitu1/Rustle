#!/usr/bin/env python
"""
compara_fetch.py — Fetch human Ensembl Compara paralogy for the NAMED genes of
the Rustle paralog universe and cache it for a NON-CIRCULAR validation phase.

This phase ONLY fetches/caches + computes mapping coverage. The actual comparison
(does Ensembl-Compara paralogy agree with our family grouping?) is a later phase.

Run with: /home/juanfra/miniforge3/bin/python

Deterministic given the cache (compara_cache.json). Reruns load the cache and
fetch ONLY missing entries (resumable). To force a full refetch, delete the cache.
"""
import csv
import json
import os
import re
import sys
import time
import urllib.error
import urllib.request

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------
BASE = os.path.dirname(os.path.abspath(__file__))
UNIVERSE_TSV = os.path.join(BASE, "copy_recovery_eval", "results", "universe.tsv")
CACHE_JSON = os.path.join(BASE, "compara_cache.json")
RELATION_JSON = os.path.join(BASE, "compara_paralog_relation.json")
REPORT_MD = os.path.join(BASE, "compara_fetch.md")

REST = "https://rest.ensembl.org"

# Network robustness knobs
TIMEOUT = 30.0          # 30s per request
MAX_RETRIES = 5         # up to 5 attempts
BACKOFF_BASE = 1.5      # seconds, exponential
SPACING = 0.4           # ~0.4s between successful requests (Ensembl rate limit)

LOC_RE = re.compile(r"^LOC[0-9]+")

# ---------------------------------------------------------------------------
# Cache (keyed by "endpoint|symbol")
# ---------------------------------------------------------------------------
def load_cache():
    if os.path.exists(CACHE_JSON):
        with open(CACHE_JSON) as f:
            return json.load(f)
    return {}


def save_cache(cache):
    tmp = CACHE_JSON + ".tmp"
    with open(tmp, "w") as f:
        json.dump(cache, f, indent=1, sort_keys=True)
    os.replace(tmp, CACHE_JSON)


def cache_key(endpoint, symbol):
    return "%s|%s" % (endpoint, symbol)


# ---------------------------------------------------------------------------
# Robust GET with retries/backoff. Returns dict:
#   {"status": "ok",     "data": <parsed json>}
#   {"status": "notfound","http": 400}          (symbol not in Ensembl)
#   {"status": "error",  "note": "..."}         (persistent failure)
# These are cached verbatim so reruns do not re-hit the API for resolved keys.
# A persistent transient "error" is NOT cached, so a later rerun retries it.
# ---------------------------------------------------------------------------
def robust_get(url):
    last_note = None
    for attempt in range(1, MAX_RETRIES + 1):
        try:
            req = urllib.request.Request(
                url, headers={"Content-Type": "application/json",
                              "User-Agent": "rustle-compara-fetch/1.0"}
            )
            with urllib.request.urlopen(req, timeout=TIMEOUT) as resp:
                body = resp.read().decode("utf-8")
                data = json.loads(body) if body.strip() else None
                return {"status": "ok", "data": data}
        except urllib.error.HTTPError as e:
            code = e.code
            if code in (400, 404):
                # Symbol genuinely not resolvable in Ensembl -> deterministic.
                return {"status": "notfound", "http": code}
            if code == 429:
                # Rate limited: honor Retry-After if present.
                ra = e.headers.get("Retry-After")
                wait = float(ra) if ra else BACKOFF_BASE ** attempt
                last_note = "http 429 (rate limit)"
                time.sleep(wait)
                continue
            # 500/502/503/504 etc -> transient, retry.
            last_note = "http %d" % code
            time.sleep(BACKOFF_BASE ** attempt)
        except (urllib.error.URLError, TimeoutError, OSError) as e:
            last_note = "neterr %s" % (e,)
            time.sleep(BACKOFF_BASE ** attempt)
        except json.JSONDecodeError as e:
            last_note = "jsonerr %s" % (e,)
            time.sleep(BACKOFF_BASE ** attempt)
    return {"status": "error", "note": last_note or "unknown"}


def fetch_cached(cache, endpoint, symbol, url):
    """Return cached resolved result, else fetch. Caches ok/notfound (terminal).
    Does NOT cache persistent transient 'error' so reruns retry."""
    key = cache_key(endpoint, symbol)
    cached = cache.get(key)
    if cached is not None and cached.get("status") in ("ok", "notfound"):
        return cached, False  # (result, did_fetch)
    res = robust_get(url)
    if res.get("status") in ("ok", "notfound"):
        cache[key] = res
        save_cache(cache)
    return res, True


# ---------------------------------------------------------------------------
# Universe parsing
# ---------------------------------------------------------------------------
def parse_universe():
    genes = set()
    families = set()
    fam_genes = {}  # family_id -> set(gene_id)
    with open(UNIVERSE_TSV) as f:
        r = csv.DictReader(f, delimiter="\t")
        for row in r:
            g = row["gene_id"].strip()
            fam = row["family_id"].strip()
            genes.add(g)
            families.add(fam)
            fam_genes.setdefault(fam, set()).add(g)
    return genes, families, fam_genes


def is_named(g):
    return not LOC_RE.match(g)


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
def main():
    genes, families, fam_genes = parse_universe()
    named = sorted(g for g in genes if is_named(g))
    loc_only = sorted(g for g in genes if not is_named(g))

    print("universe.tsv: %s" % UNIVERSE_TSV)
    print("total genes      : %d" % len(genes))
    print("named genes      : %d" % len(named))
    print("LOC-only genes   : %d" % len(loc_only))
    print("families         : %d" % len(families))
    print()

    cache = load_cache()
    fetched = 0

    # gene_symbol -> {"ensg": str|None, "paralogues": [ENSG...], "status": ...,
    #                 "para_status": ...}
    gene_info = {}

    for i, sym in enumerate(named, 1):
        # 1) lookup ENSG id
        lk_url = ("%s/lookup/symbol/homo_sapiens/%s?content-type=application/json"
                  % (REST, sym))
        lk, did = fetch_cached(cache, "lookup", sym, lk_url)
        if did:
            fetched += 1
            time.sleep(SPACING)

        ensg = None
        lk_status = lk.get("status")
        if lk_status == "ok" and isinstance(lk.get("data"), dict):
            ensg = lk["data"].get("id")
        # notfound -> unmapped (symbol not in Ensembl)

        # 2) paralogues (condensed) — fetch regardless; keyed by symbol
        pa_url = ("%s/homology/symbol/homo_sapiens/%s?type=paralogues;"
                  "compara=vertebrates;content-type=application/json;format=condensed"
                  % (REST, sym))
        para_ids = []
        para_status = "skip"
        if lk_status != "notfound":
            pa, did = fetch_cached(cache, "homology_paralogues", sym, pa_url)
            if did:
                fetched += 1
                time.sleep(SPACING)
            para_status = pa.get("status")
            if para_status == "ok":
                data = pa.get("data") or {}
                homs = []
                d = data.get("data") if isinstance(data, dict) else None
                if isinstance(d, list) and d:
                    homs = d[0].get("homologies", []) or []
                para_ids = sorted({h.get("id") for h in homs if h.get("id")})
            elif para_status == "notfound":
                para_ids = []
        else:
            para_status = "notfound"

        gene_info[sym] = {
            "ensg": ensg,
            "lookup_status": lk_status,
            "paralogues": para_ids,
            "para_status": para_status,
        }
        flag = "OK " if ensg else "UNMAPPED"
        print("[%2d/%2d] %-14s ensg=%-18s paralogues=%-4d (%s)"
              % (i, len(named), sym, ensg or "-", len(para_ids), flag))

    save_cache(cache)

    # -----------------------------------------------------------------------
    # Build symmetric within-universe paralog relation (restricted to OUR genes)
    # X ~ Y iff ENSID(Y) in paralogues(X) OR ENSID(X) in paralogues(Y)
    # -----------------------------------------------------------------------
    mapped_syms = [s for s in named if gene_info[s]["ensg"]]
    ensg_to_sym = {}
    for s in mapped_syms:
        ensg_to_sym.setdefault(gene_info[s]["ensg"], []).append(s)

    relation_pairs = set()
    for x in mapped_syms:
        ex = gene_info[x]["ensg"]
        px = set(gene_info[x]["paralogues"])
        for y in mapped_syms:
            if x >= y:
                continue
            ey = gene_info[y]["ensg"]
            py = set(gene_info[y]["paralogues"])
            if (ey in px) or (ex in py):
                relation_pairs.add((x, y))

    # -----------------------------------------------------------------------
    # Coverage on families: families with >=2 NAMED genes that BOTH mapped
    # (i.e. have a checkable within-family pair via Ensembl ids).
    # -----------------------------------------------------------------------
    fam_named_genes = {}      # family -> sorted named genes in family
    fam_mapped_genes = {}     # family -> named genes that got an ENSG
    for fam, gs in fam_genes.items():
        nm = sorted(g for g in gs if is_named(g))
        fam_named_genes[fam] = nm
        fam_mapped_genes[fam] = [g for g in nm if gene_info.get(g, {}).get("ensg")]

    fam_ge2_named = [f for f, nm in fam_named_genes.items() if len(nm) >= 2]
    fam_ge2_mapped = [f for f, mg in fam_mapped_genes.items() if len(mg) >= 2]

    # Coverage numbers
    n_ensg = sum(1 for s in named if gene_info[s]["ensg"])
    n_para = sum(1 for s in named if gene_info[s]["para_status"] == "ok")
    n_para_nonempty = sum(1 for s in named if gene_info[s]["paralogues"])
    n_unmapped = sum(1 for s in named if gene_info[s]["lookup_status"] == "notfound")
    n_lookup_err = sum(1 for s in named if gene_info[s]["lookup_status"] == "error")
    n_para_err = sum(1 for s in named
                     if gene_info[s]["para_status"] == "error")

    # -----------------------------------------------------------------------
    # Write relation summary JSON
    # -----------------------------------------------------------------------
    relation = {
        "meta": {
            "total_genes": len(genes),
            "named_genes": len(named),
            "loc_only_genes": len(loc_only),
            "families": len(families),
            "named_genes_with_ensg": n_ensg,
            "named_genes_with_paralogue_data": n_para,
            "named_genes_with_nonempty_paralogues": n_para_nonempty,
            "named_genes_unmapped": n_unmapped,
            "lookup_errors": n_lookup_err,
            "paralogue_errors": n_para_err,
            "families_with_ge2_named": len(fam_ge2_named),
            "families_with_ge2_mapped": len(fam_ge2_mapped),
            "within_universe_paralog_pairs": len(relation_pairs),
            "rest_endpoint": REST,
            "compara": "vertebrates",
        },
        "genes": {
            s: {
                "ensg": gene_info[s]["ensg"],
                "lookup_status": gene_info[s]["lookup_status"],
                "para_status": gene_info[s]["para_status"],
                "paralogues": gene_info[s]["paralogues"],
            }
            for s in named
        },
        "families": {
            f: {
                "named_genes": fam_named_genes[f],
                "mapped_genes": fam_mapped_genes[f],
                "checkable_pair": len(fam_mapped_genes[f]) >= 2,
            }
            for f in sorted(families)
        },
        "within_universe_paralog_relation": sorted(
            [list(p) for p in relation_pairs]
        ),
    }
    with open(RELATION_JSON, "w") as f:
        json.dump(relation, f, indent=2, sort_keys=True)

    # -----------------------------------------------------------------------
    # Write report markdown
    # -----------------------------------------------------------------------
    unmapped_syms = sorted(s for s in named
                           if gene_info[s]["lookup_status"] == "notfound")
    err_syms = sorted(s for s in named
                      if gene_info[s]["lookup_status"] == "error"
                      or gene_info[s]["para_status"] == "error")
    md = []
    md.append("# Ensembl Compara paralogy fetch — coverage report\n")
    md.append("Non-circular validation prep. This phase ONLY fetches/caches "
              "Ensembl Compara paralogy for the NAMED universe genes; the "
              "comparison vs our family grouping is the next phase.\n")
    md.append("- Cache: `bench/compara_cache.json` (keyed by `endpoint|symbol`, "
              "resumable; reruns fetch only missing).\n")
    md.append("- Relation summary: `bench/compara_paralog_relation.json`\n")
    md.append("- Source universe: `bench/copy_recovery_eval/results/universe.tsv`\n")
    md.append("\n## Universe gene inventory\n")
    md.append("| metric | count |")
    md.append("| --- | --- |")
    md.append("| total distinct gene_ids | %d |" % len(genes))
    md.append("| named genes (not `^LOC[0-9]+`) | %d |" % len(named))
    md.append("| LOC-only genes | %d |" % len(loc_only))
    md.append("| families | %d |" % len(families))
    md.append("\n## Mapping coverage (named genes)\n")
    md.append("| metric | count | of named |")
    md.append("| --- | --- | --- |")
    md.append("| got ENSG id | %d | %d |" % (n_ensg, len(named)))
    md.append("| got paralogue data (HTTP ok) | %d | %d |" % (n_para, len(named)))
    md.append("| non-empty paralogue set | %d | %d |"
              % (n_para_nonempty, len(named)))
    md.append("| unmapped (symbol not in Ensembl) | %d | %d |"
              % (n_unmapped, len(named)))
    md.append("| persistent fetch errors | %d | %d |"
              % (len(err_syms), len(named)))
    md.append("\n## Family-level checkability\n")
    md.append("| metric | count | of 62 families |")
    md.append("| --- | --- | --- |")
    md.append("| families with >=2 NAMED genes | %d | %d |"
              % (len(fam_ge2_named), len(families)))
    md.append("| families with >=2 MAPPED genes (checkable within-family pair) "
              "| %d | %d |" % (len(fam_ge2_mapped), len(families)))
    md.append("\nWithin-universe symmetric paralog pairs (both genes mapped, "
              "Compara-linked): **%d**\n" % len(relation_pairs))
    if unmapped_syms:
        md.append("\n### Unmapped symbols (not in Ensembl)\n")
        md.append(", ".join("`%s`" % s for s in unmapped_syms) + "\n")
    if err_syms:
        md.append("\n### Symbols with persistent fetch errors (rerun to retry)\n")
        md.append(", ".join("`%s`" % s for s in err_syms) + "\n")
    md.append("\nDeterministic given the cache. Fetched %d new API responses "
              "this run.\n" % fetched)
    with open(REPORT_MD, "w") as f:
        f.write("\n".join(md) + "\n")

    # -----------------------------------------------------------------------
    # Console summary
    # -----------------------------------------------------------------------
    print()
    print("=== COVERAGE ===")
    print("named genes                         : %d" % len(named))
    print("  got ENSG id                       : %d" % n_ensg)
    print("  got paralogue data (HTTP ok)      : %d" % n_para)
    print("  non-empty paralogue set           : %d" % n_para_nonempty)
    print("  unmapped (not in Ensembl)         : %d" % n_unmapped)
    print("  persistent fetch errors           : %d" % len(err_syms))
    print("families                            : %d" % len(families))
    print("  with >=2 NAMED genes              : %d" % len(fam_ge2_named))
    print("  with >=2 MAPPED genes (checkable) : %d" % len(fam_ge2_mapped))
    print("within-universe paralog pairs       : %d" % len(relation_pairs))
    print("new API responses fetched this run  : %d" % fetched)
    print()
    print("wrote: %s" % CACHE_JSON)
    print("wrote: %s" % RELATION_JSON)
    print("wrote: %s" % REPORT_MD)


if __name__ == "__main__":
    sys.exit(main())
