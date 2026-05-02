#!/usr/bin/env python3
"""Classify ST + rustle GTFs against a panel of multi-copy ref families.

For each family in the panel, runs gffcompare with that family's ref
GFF, then counts how many distinct ref gene_ids are caught (matched at
locus level) by each input GTF. Produces a single wide table summary.

Usage: classify_panel.py <st.gtf> <em.gtf>

Reads ref_<FAMILY>.gff in the working directory.
"""
import os, re, sys, subprocess, tempfile

if len(sys.argv) != 3:
    sys.exit("usage: classify_panel.py <st.gtf> <em.gtf>")

GFFCOMPARE = "/storage/work/jxi21/.pixi_cache/pkgs/gffcompare-0.12.10-h9948957_0/bin/gffcompare"
queries = {"ST": sys.argv[1], "rustle-EM": sys.argv[2]}
families = ["GOLGA6", "GOLGA8", "AMY", "TBC1D3", "NBPF", "OR", "MUC", "HBB_HBA", "KRABZNF"]

def n_ref_genes(reffile):
    """Count ref gene records (denominator)."""
    n = 0
    with open(reffile) as fh:
        for line in fh:
            if line.startswith("#"): continue
            p = line.rstrip("\n").split("\t")
            if len(p) >= 9 and p[2] == "gene":
                n += 1
    return n

def caught_genes(loci_path):
    """Return set of ref gene_ids that are caught (XLOC has both ref and query)."""
    caught = set()
    if not os.path.exists(loci_path): return caught
    with open(loci_path) as fh:
        for line in fh:
            p = line.rstrip("\n").split("\t")
            if len(p) < 4: continue
            ref_col, q_col = p[2], p[3]
            if ref_col == "-" or q_col == "-": continue
            for entry in ref_col.split(","):
                entry = entry.strip()
                if not entry: continue
                gene = entry.split("|", 1)[0]
                caught.add(gene)
    return caught

def parse_stats(stats_path):
    """Pull the key gffcompare summary numbers."""
    out = {}
    if not os.path.exists(stats_path): return out
    with open(stats_path) as fh:
        for line in fh:
            line = line.rstrip()
            for k, lbl in [("Locus level", "locus"), ("Intron chain level", "ichain"),
                           ("Transcript level", "tx")]:
                m = re.match(rf"\s*{re.escape(k)}\s*:\s*([\d.]+)\s*\|", line)
                if m: out[f"sn_{lbl}"] = float(m.group(1))
            m = re.match(r"\s*Matching transcripts:\s*(\d+)", line)
            if m: out["match_tx"] = int(m.group(1))
            m = re.match(r"\s*Matching intron chains:\s*(\d+)", line)
            if m: out["match_ichain"] = int(m.group(1))
            m = re.match(r"\s*Matching loci:\s*(\d+)", line)
            if m: out["match_loci"] = int(m.group(1))
    return out

# Compute everything
results = {}
for fam in families:
    ref = f"ref_{fam}.gff"
    if not os.path.exists(ref):
        continue
    n_ref = n_ref_genes(ref)
    fam_data = {"n_ref": n_ref, "by_query": {}}
    with tempfile.TemporaryDirectory() as tmp:
        for label, q in queries.items():
            prefix = os.path.join(tmp, f"cmp_{fam}_{label.replace('+','_').replace('-','_')}")
            try:
                subprocess.run([GFFCOMPARE, "-r", ref, "-o", prefix, "--no-merge", q],
                               stdout=subprocess.PIPE, stderr=subprocess.PIPE, check=True)
            except subprocess.CalledProcessError:
                continue
            n_caught = len(caught_genes(prefix + ".loci"))
            stats = parse_stats(prefix + ".stats")
            fam_data["by_query"][label] = {
                "n_caught": n_caught,
                "match_tx": stats.get("match_tx", 0),
                "match_ichain": stats.get("match_ichain", 0),
                "match_loci": stats.get("match_loci", 0),
            }
    results[fam] = fam_data

# Print panel summary
print("=" * 110)
print(f"{'family':10}  {'ref_genes':>9}  | {'ST_caught':>10} {'EM_caught':>10}  Δ  | {'ST_tx':>6} {'EM_tx':>6}  Δ_tx")
print("=" * 110)
totals = {"n_ref":0, "ST_caught":0, "EM_caught":0, "ST_tx":0, "EM_tx":0}
for fam in families:
    if fam not in results:
        print(f"{fam:10}  (skipped — no ref subset)")
        continue
    r = results[fam]
    nref = r["n_ref"]
    st = r["by_query"].get("ST", {"n_caught":0, "match_tx":0})
    em = r["by_query"].get("rustle-EM", {"n_caught":0, "match_tx":0})
    delta = em["n_caught"] - st["n_caught"]
    delta_tx = em["match_tx"] - st["match_tx"]
    print(f"{fam:10}  {nref:>9}  | {st['n_caught']:>10} {em['n_caught']:>10}  {delta:+}  | "
          f"{st['match_tx']:>6} {em['match_tx']:>6}  {delta_tx:+}")
    totals["n_ref"] += nref
    totals["ST_caught"] += st["n_caught"]
    totals["EM_caught"] += em["n_caught"]
    totals["ST_tx"] += st["match_tx"]
    totals["EM_tx"] += em["match_tx"]

print("-" * 110)
delta_t = totals["EM_caught"] - totals["ST_caught"]
delta_tx_t = totals["EM_tx"] - totals["ST_tx"]
print(f"{'TOTAL':10}  {totals['n_ref']:>9}  | "
      f"{totals['ST_caught']:>10} {totals['EM_caught']:>10}  {delta_t:+}  | "
      f"{totals['ST_tx']:>6} {totals['EM_tx']:>6}  {delta_tx_t:+}")
