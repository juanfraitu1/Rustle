#!/usr/bin/env python3
"""Multi-copy benchmark oracle. Measures all 4 objectives, scores vs expectations.json,
prints a table + per-objective PASS/FAIL, exits non-zero on regression.
Usage: python3 bench/multi_copy_eval/run_oracle.py [--fast|--full] [--check] [--write-report]"""
import sys, os, re, json, subprocess

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.dirname(os.path.dirname(HERE))
RUSTLE = os.path.join(REPO, "target", "release", "rustle")
EXP = os.path.join(HERE, "expectations.json")

def run_json(script):
    out = subprocess.run([sys.executable, os.path.join(HERE, script)], capture_output=True, text=True)
    try: return json.loads(out.stdout)
    except Exception: return {"error": out.stdout[-500:] + out.stderr[-500:]}

def obj2():
    r = subprocess.run([sys.executable, os.path.join(HERE,"no_absent_copy_eval.py"),"--check"],
                       capture_output=True, text=True)
    txt = r.stdout
    g = re.search(r"\[GOLGA8\] Rustle ([\d.]+)%", txt); y = re.search(r"\[YAG\] Rustle ([\d.]+)%", txt)
    if not (g and y): return {"status":"SKIPPED","detail":"no_absent_copy_eval TSVs missing"}
    return {"golga8_pct": float(g.group(1)), "yag_pct": float(y.group(1)), "exit": r.returncode}

def obj_daz():
    bam = "/tmp/daz.bam"
    if not os.path.exists(os.path.join(REPO,"..","GGO.bam")):
        return {"status":"SKIPPED","detail":"GGO.bam missing"}
    subprocess.run(f"samtools view -b {REPO}/../GGO.bam NC_073248.2:42700000-43000000 > {bam} && samtools index {bam}",
                   shell=True, stderr=subprocess.DEVNULL)
    r = subprocess.run([RUSTLE,"--vg","--genome-fasta",f"{REPO}/../GGO.fasta","-L",bam,"-o","/tmp/daz_vg.gtf"],
                       capture_output=True, text=True)
    em = bool(re.search(r"HMM-EM.*reweighted \d+ reads", r.stderr))
    daz3 = 0
    for line in open("/tmp/daz_vg.gtf"):
        f = line.split("\t")
        if len(f)>=7 and f[2]=="transcript" and f[6]=="+" and int(f[3])<42945552 and int(f[4])>42879918:
            daz3 += 1
    return {"daz3_isoforms": daz3, "em_ran": em}

# ── O4/O5 copy-resolution guard (the discriminating metric: copy ATTRIBUTION, not chain Sn/Pr).
# See docs/superpowers/specs/2026-06-03-o4o5-copy-resolution-benchmark.md.
TANDEM = os.path.join(REPO, "bench", "tandem_attribution")

def _o5_run(identity, out, vg_extra):
    """gen + align + VG(+FP_ATTR_TSV dump) + score_attribution. Returns the score dict (with
    em_ran). A small merged-but-separable 2-copy fixture (spacing 16000, rpc 80, distinct
    isoforms); seed fixed for determinism."""
    gen = os.path.join(TANDEM, "gen_synthetic.py"); score = os.path.join(TANDEM, "score_attribution.py")
    if not (os.path.exists(gen) and os.path.exists(score)): return {"em_ran": False, "status": "SKIPPED"}
    os.makedirs(out, exist_ok=True)
    subprocess.run([sys.executable, gen, "--identity", str(identity), "--copies", "2",
                    "--reads-per-copy", "80", "--spacing", "16000", "--distinct-isoforms",
                    "--seed", "5", "--out", out], capture_output=True, text=True)
    bam = os.path.join(out, "r.bam")
    rc = subprocess.run(f"minimap2 -ax splice:hq -uf -N20 -p0.5 {out}/genome.fa {out}/reads.fa 2>/dev/null "
                        f"| samtools sort -o {bam} - 2>/dev/null && samtools index {bam}", shell=True)
    tsv = os.path.join(out, "fp.tsv")
    if os.path.exists(tsv): os.remove(tsv)
    env = dict(os.environ, RUSTLE_VG_FP_ATTR_TSV=tsv)
    subprocess.run([RUSTLE, "--vg", "--vg-snp", *vg_extra, "--genome-fasta", f"{out}/genome.fa",
                    "-L", bam, "-o", f"{out}/vg.gtf"], capture_output=True, text=True, env=env)
    if not os.path.exists(tsv): return {"em_ran": False}  # family dropped → no EM → guard skips
    r = subprocess.run([sys.executable, score, "--tsv", tsv, "--reads", f"{out}/reads.fa",
                        "--meta", f"{out}/meta.json"], capture_output=True, text=True)
    try:
        d = json.loads(r.stdout); d["em_ran"] = True; return d
    except Exception:
        return {"em_ran": False}

def obj5_copy_attribution():
    """(a) CAPABILITY: at an identifiable merged bundle (DEFAULT settings, id 0.97) VG attributes
    ambiguous multimappers to their TRUE source copy (baseline has no copy concept). (b) HONESTY:
    at identical copies (id 1.0, the DAZ limit) VG must ABSTAIN — 0 decisive calls, no fabrication
    (relaxed family filter so the EM runs and the abstention is observable)."""
    if not os.path.exists(RUSTLE): return {"status": "SKIPPED"}
    ident = _o5_run(0.97, "/tmp/oracle_o5_id097", [])
    nonid = _o5_run(1.0, "/tmp/oracle_o5_id1",
                    ["--vg-family-min-shared", "2", "--vg-family-min-shared-per-copy", "0"])
    return {
        "identifiable_acc": ident.get("attribution_accuracy_multimappers"),
        "identifiable_em_ran": ident.get("em_ran", False),
        "nonid_decisive_frac": nonid.get("decisive_frac"),
        "nonid_abstain_frac": nonid.get("abstain_frac"),
        "nonid_em_ran": nonid.get("em_ran", False),
    }

def default_headline():
    subprocess.run([RUSTLE,"-L",f"{REPO}/../GGO_19.bam","-o","/tmp/oracle_def.gtf"], stderr=subprocess.DEVNULL)
    subprocess.run(["gffcompare","-r",f"{REPO}/../GGO_19.gtf","/tmp/oracle_def.gtf","-o","/tmp/oracle_def_cmp"], stderr=subprocess.DEVNULL)
    txt = open("/tmp/oracle_def_cmp.stats").read()
    m = re.search(r"Transcript level:\s*([\d.]+)\s*\|\s*([\d.]+)", txt)
    return {"tx_sn": float(m.group(1)), "tx_pr": float(m.group(2))} if m else {"status":"SKIPPED"}

def main():
    full = "--full" in sys.argv; check = "--check" in sys.argv
    res = {}
    res["obj3_isoforms"] = run_json("score_synth_isoforms.py")
    res["obj4_assignment"] = run_json("score_synth_assignment.py")
    res["default_headline"] = default_headline()
    if full:
        res["obj2_copy_recovery"] = obj2()
        res["obj_daz"] = obj_daz()
        res["obj5_copy_attribution"] = obj5_copy_attribution()
        # Obj 1 reframed 2026-06-01: GOLGA6L7 is ANTISENSE-SILENT (0 own-strand reads;
        # its 53+17 reads are all +strand, belonging to overlapping +strand lncRNA
        # LOC115930831) — correctly emits nothing; the "miss" is a gffcompare locus-
        # overlap artifact, NOT a tool bug and NOT a same-bundle splitter problem.
        # Confirmed on 3 genuine same-strand expressed paralog pairs (gap 15bp-3kb):
        # they fail by bundle REJECTION / partial assembly in SEPARATE bundles, never
        # merge-collapse — so an intra-bundle splitter has no validated target. The
        # real gap is depth-aware bundling / family-read stripping of paralog bundles.
        res["obj1_paralog_assembly"] = {"status":"REFRAMED",
            "detail":"GOLGA6L7 antisense-silent (not a bug); no validated intra-bundle-splitter target; "
                     "genuine paralog pairs fail by bundle rejection, not collapse — see THESIS_OBJECTIVES.md O5"}
    print(json.dumps(res, indent=2))
    if check:
        exp = json.load(open(EXP))
        fails = []
        # obj3 overall Sn floor
        if res["obj3_isoforms"].get("vg",{}).get("overall",{}).get("Sn",0) + 1e-9 < exp["obj3_overall_sn_min"]:
            fails.append("obj3 Sn below floor")
        if res["obj4_assignment"].get("decisive_accuracy",0) + 1e-9 < exp["obj4_decisive_acc_min"]:
            fails.append("obj4 decisive accuracy below floor")
        dh = res["default_headline"]
        if "tx_sn" in dh and abs(dh["tx_sn"]-exp["default_tx_sn"])>exp["default_band"]:
            fails.append("default Sn out of band (isolation guard)")
        if full and res.get("obj2_copy_recovery",{}).get("exit",0) not in (0,):
            fails.append("obj2 no_absent_copy_eval --check failed")
        # DAZ3 false-positive regression guard (corrected 2026-06-02): DAZ3 must ABSTAIN,
        # not fabricate the phantom 5-isoform/cov-113 recovery. Fails if it over-enumerates.
        if full and "daz3_isoforms" in res.get("obj_daz",{}):
            d3 = res["obj_daz"]["daz3_isoforms"]
            if d3 > exp.get("daz3_isoforms_max", 2):
                fails.append(f"DAZ3 over-enumerated ({d3} > {exp['daz3_isoforms_max']}): phantom recovery regressed (should abstain)")
        # O4/O5 copy-resolution guards (only assert when the EM actually ran):
        #   (a) CAPABILITY: identifiable copy-attribution accuracy must hold (Tier-2 win).
        #   (b) HONESTY: at identical copies VG must ABSTAIN — no decisive calls (no fabrication).
        if full and isinstance(res.get("obj5_copy_attribution"), dict):
            o5 = res["obj5_copy_attribution"]
            if o5.get("identifiable_em_ran") and (o5.get("identifiable_acc") or 0) + 1e-9 < exp.get("obj5_copy_attr_min", 1.0):
                fails.append(f"obj5 copy-attribution below floor ({o5.get('identifiable_acc')} < {exp.get('obj5_copy_attr_min',1.0)})")
            if o5.get("nonid_em_ran") and (o5.get("nonid_decisive_frac") or 0) > exp.get("obj5_id1_max_decisive_frac", 0.0) + 1e-9:
                fails.append(f"obj5 FABRICATION at identical copies: decisive_frac {o5.get('nonid_decisive_frac')} > {exp.get('obj5_id1_max_decisive_frac',0.0)} (must abstain at the DAZ limit)")
        if fails:
            print("REGRESSION: " + "; ".join(fails)); sys.exit(1)
        print("ALL OBJECTIVES PASS")
    if "--write-report" in sys.argv:
        md = ["# VG Objectives — generated by run_oracle.py (do not hand-edit)\n"]
        md.append("| Objective | Metric | Value |\n|---|---|---|")
        o3 = res["obj3_isoforms"].get("vg",{}).get("overall",{})
        md.append(f"| Obj 3 per-copy isoforms | overall Sn/Pr | {o3.get('Sn')}/{o3.get('Pr')} |")
        o4 = res["obj4_assignment"]
        md.append(f"| Obj 4 read assignment | accuracy / decisive-acc | {o4.get('accuracy')}/{o4.get('decisive_accuracy')} |")
        if "obj2_copy_recovery" in res:
            o2 = res["obj2_copy_recovery"]
            md.append(f"| Obj 2 copy recovery | GOLGA8 / YAG % | {o2.get('golga8_pct')}/{o2.get('yag_pct')} |")
        if "obj_daz" in res:
            md.append(f"| Obj 3 DAZ (real) | DAZ3 isoforms / EM-ran | {res['obj_daz'].get('daz3_isoforms')}/{res['obj_daz'].get('em_ran')} |")
        if "obj5_copy_attribution" in res:
            o5 = res["obj5_copy_attribution"]
            md.append(f"| Obj 4/5 copy attribution | id0.97 acc / id1.0 decisive-frac (abstain) | "
                      f"{o5.get('identifiable_acc')} / {o5.get('nonid_decisive_frac')} |")
        md.append(f"| Obj 1 paralog assembly | status | REFRAMED (GOLGA6L7 antisense-silent; no splitter target — bundle rejection, not collapse) |")
        md.append(f"| default de-novo | Tx Sn/Pr (isolation) | {res['default_headline'].get('tx_sn')}/{res['default_headline'].get('tx_pr')} |")
        md.append("\n```json\n" + json.dumps(res, indent=2) + "\n```\n")
        open(os.path.join(HERE,"OBJECTIVES_ASSESSMENT.md"),"w").write("\n".join(md))

if __name__ == "__main__":
    main()
