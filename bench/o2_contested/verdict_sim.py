#!/usr/bin/env python3
"""r1705 (docs/PREREG_placement_verdicts_2026-11-13.md, section Row 1705): planted cases with known read origin for the
placement verdicts. Each case = a famsim scenario (make, verify, align with the shipped minimap2 command), as_table,
placement_verdicts.py verdicts at 0.8, then verdicts against each read's home copy.
Usage: verdict_sim.py OUTDIR AS_TABLE_BINARY [CASE ...]     (needs minimap2, samtools, pysam on PATH)"""
import collections, csv, json, os, subprocess, sys

HERE = os.path.dirname(os.path.abspath(__file__)); BENCH = os.path.dirname(HERE)
snp = lambda r: {"op": "snp", "rate": r}
conv = lambda *ex: [{"op": "convert", "from": "A2", "exon": e} for e in ex]
absent = lambda cid, ops: {"id": cid, "in_reference": False, "ops": ops}
# case -> (copies, {absent copy: home copy or None})
CASES = {
    "identical": ([{"id": "A"}, {"id": "A2", "ops": []}], {}),
    "psv_0.001": ([{"id": "A"}, {"id": "A2", "ops": [snp(0.001)]}], {}),
    "psv_0.005": ([{"id": "A"}, {"id": "A2", "ops": [snp(0.005)]}], {}),
    "psv_0.02": ([{"id": "A"}, {"id": "A2", "ops": [snp(0.02)]}], {}),
    "retro_0.002": ([{"id": "A"}, {"id": "AR", "ops": [{"op": "retro"}, snp(0.002)]}], {}),
    "retro_0.01": ([{"id": "A"}, {"id": "AR", "ops": [{"op": "retro"}, snp(0.01)]}], {}),
    "unexpressed": ([{"id": "A"}, {"id": "A2", "expression": 0, "ops": [snp(0.005)]}], {}),
    "absent_copy": ([{"id": "A"}, {"id": "A2", "ops": [snp(0.005)]}, absent("A3", [snp(0.01)])], {"A3": None}),
    "allele_private": ([{"id": "A"}, {"id": "A2", "ops": [snp(0.005)]}, absent("Astar", [snp(0.003)])], {"Astar": "A"}),
    "allele_minor": ([{"id": "A"}, {"id": "A2", "ops": [snp(0.01)]}, absent("Astar", conv(3))], {"Astar": "A"}),
    "allele_minor6": ([{"id": "A"}, {"id": "A2", "ops": [snp(0.01)]}, absent("Astar", conv(6))], {"Astar": "A"}),  # amendment: exon 3 carried no SNP
    "allele_major": ([{"id": "A"}, {"id": "A2", "ops": [snp(0.01)]}, absent("Astar", conv(2, 3, 4, 5, 6))], {"Astar": "A"}),
}
READS = {"perfect": {"err": 0.0, "indel": 0.0}, "hifi": {"err": 0.001, "indel": 0.0003}}
VERD = ("excluded", "identical", "undecided", "over_primary")

def sh(cmd, log):
    with open(log, "a") as fh:
        if subprocess.run(cmd, shell=True, stdout=fh, stderr=fh, cwd=os.path.dirname(BENCH)).returncode: sys.exit(f"FAILED: {cmd} (see {log})")

def run(out, as_table, case, mode):
    d = os.path.join(out, f"{case}.{mode}"); os.makedirs(d, exist_ok=True); copies, home = CASES[case]
    spec = {"name": case, "seed": 5, "background": {"source": "random", "length": 200000},
            "template": {"source": "synthetic", "exons": [220, 150, 300, 130, 90, 600], "introns": [1200, 2500, 800, 1500, 900]},
            "copies": copies, "reads": dict(READS[mode], per_copy=200, jitter=30)}
    json.dump(spec, open(d + ".json", "w")); log = d + ".log"; open(log, "w").close()
    sh(f"python3 bench/famsim make {d}.json --out {d}", log)
    subprocess.run(f"python3 bench/famsim verify {d}", shell=True, stdout=open(log, "a"), stderr=subprocess.STDOUT, cwd=os.path.dirname(BENCH))
    fails = [l.split("\t")[0] for l in open(os.path.join(d, "verify.tsv")) if "\tFAIL\t" in l]
    sh(f"python3 bench/famsim align {d}", log)
    sh(f"{as_table} --bam {d}/reads.bam --out {d}/reads.as_table.tsv && python3 {HERE}/placement_verdicts.py verdicts {d}/reads.bam {d}/genome.ref.fa {d}/reads.as_table.tsv 0.8 > {d}/verdicts.tsv", log)
    span = {r["copy_id"]: (r["contig"], int(r["start"]), int(r["end"])) for r in csv.DictReader(open(os.path.join(d, "copies.tsv")), delimiter="\t") if r["in_reference"] == "1"}
    c = collections.Counter(); seen = set()
    for r in csv.DictReader(open(os.path.join(d, "verdicts.tsv")), delimiter="\t"):
        cp = r["read"].split("|")[0]; h = home.get(cp, cp)
        on = h is not None and h in span and span[h][0] == r["chrom"] and span[h][1] < int(r["end"]) and int(r["pos"]) - 1 < span[h][2]
        where = "nohome" if h is None else "home" if on else "off"
        if r["primary"] == "1": c[f"primary {where}"] += 1; seen.add(r["read"])
        else: c[f"{where} {r['verdict']}"] += 1
    nreads = sum(1 for l in open(os.path.join(d, "reads.truth.tsv"))) - 1
    return dict(case=case, mode=mode, verify="FAIL:" + ",".join(sorted(set(fails))) if fails else "PASS", reads=nreads, judged=len(seen), c=c)

def main(out, as_table, *only):
    os.makedirs(out, exist_ok=True); rows = [run(out, as_table, k, m) for k in (only or CASES) for m in READS]
    with open(os.path.join(out, "summary.tsv"), "w") as fh:
        hdr = ["case", "reads", "verify", "simulated", "with_verdicts", "primary_home", "primary_off", "primary_nohome"] + [f"{w}_{v}" for w in ("home", "off", "nohome") for v in VERD]
        for o in (fh, sys.stdout): print(*hdr, sep="\t", file=o)
        for r in rows:
            c = r["c"]; line = [r["case"], r["mode"], r["verify"], r["reads"], r["judged"], c["primary home"], c["primary off"], c["primary nohome"]] + [c[f"{w} {v}"] for w in ("home", "off", "nohome") for v in VERD]
            for o in (fh, sys.stdout): print(*line, sep="\t", file=o)
    ok = lambda r: r["verify"] == "PASS"; hs = lambda r: sum(r["c"][f"home {v}"] for v in VERD); offs = lambda r: sum(r["c"][f"off {v}"] for v in VERD)
    b1 = [(r["case"], r["c"]["home excluded"]) for r in rows if r["mode"] == "perfect" and r["case"] not in ("allele_major", "allele_minor6") and ok(r)]
    b2 = [(r["case"], r["c"]["home excluded"], hs(r)) for r in rows if r["mode"] == "hifi" and r["case"] not in ("allele_major", "allele_minor6") and ok(r)]
    b3 = [(r["case"], r["c"]["off excluded"], offs(r)) for r in rows if r["mode"] == "perfect" and r["case"] in ("psv_0.005", "psv_0.02", "unexpressed", "retro_0.01") and ok(r)]
    print("B1 perfect reads, excluded on home (bar 0 each):", b1, "PASS" if all(n == 0 for _, n in b1) else "FAIL")
    print("B2 errors, excluded on home / home secondaries (bar <= 1% each):", b2, "PASS" if all(n <= 0.01 * t for _, n, t in b2) else "FAIL")
    print("B3 perfect reads, off-home excluded / off-home secondaries (bar >= 90% each):", b3, "PASS" if all(t and n >= 0.9 * t for _, n, t in b3) else "FAIL")
    print("voided (famsim verify failed):", [(r["case"], r["mode"], r["verify"]) for r in rows if not ok(r)])

if __name__ == "__main__": main(*sys.argv[1:])
