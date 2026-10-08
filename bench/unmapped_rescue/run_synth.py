#!/usr/bin/env python3
"""Amendment 10 pipeline on the synthetic world (docs/PREREG_unmapped_rescue_2026-10-08.md). Run with /home/juanfra/miniforge3/bin/python (edlib, pyabpoa).

    run_synth.py <E0|E1|E2> <map|cluster|truth|attribute|augment|report>
World from synth_world.py at W0; per-variant work in W0/<variant>/. Stages are resumable; 'cluster' exits 75 between minimap2 rounds."""
import collections
import csv
import json
import os
import statistics
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import attribute as A  # noqa: E402
import augment as G  # noqa: E402
import ends as E  # noqa: E402
import graph as GR  # noqa: E402
import run_augment as RA  # noqa: E402
import score as S  # noqa: E402
import seeds as SD  # noqa: E402

W0 = os.environ.get("SYNTH_DIR", "/mnt/linuxdisk/tmp/o3_rescue/synth")
DELTA = 0.00958
BLAST = "/home/juanfra/miniforge3/envs/blast/bin"
RC = str.maketrans("ACGTacgt", "TGCAtgca")


def rc(s):
    return s.translate(RC)[::-1]


def copies():
    return {r["copy"]: r for r in csv.DictReader(open(f"{W0}/copies.tsv"), delimiter="\t")}


def stage_map(V):
    d = f"{W0}/{V}"
    os.makedirs(d, exist_ok=True)
    cp = copies()
    lines = open(f"{W0}/reads.{V}.fq").read().splitlines()
    seqs = {}
    for i in range(0, len(lines), 4):
        seqs[lines[i][1:].split()[0].replace("|", ".")] = lines[i + 1]
    SD.write_fa(f"{d}/reads.fa", seqs, sorted(seqs))
    subprocess.run(f"minimap2 {RA.MM2} -t 4 {W0}/genome.ref.fa {d}/reads.fa > {d}/ref.sam", shell=True, check=True)
    prim = RA.primaries(open(f"{d}/ref.sam"))
    rows, net = [], []
    for n in sorted(seqs):
        c = n.split(".")[0]
        role = "E" if cp[c]["role"] == "E" else "S"
        rec = prim.get(n)
        innet = rec is None or rec["de"] > DELTA
        rows.append((n, c, cp[c]["family"], role, cp[c]["D"], int(innet), "" if rec is None else rec["de"]))
        if innet:
            net.append(n)
    with open(f"{d}/labels.tsv", "w") as o:
        o.write("read\tcopy\tfamily\trole\tD\tin_net\tde\n")
        for r in rows:
            o.write("\t".join(str(x) for x in r) + "\n")
    SD.write_fa(f"{d}/pool.fa", seqs, net)
    json.dump(prim, open(f"{d}/genome.json", "w"))
    by = collections.Counter((r[3], r[5]) for r in rows)
    print(f"{V}: {len(seqs)} reads, net {len(net)} (erased-copy reads in net {by[('E', 1)]} of {by[('E', 1)] + by[('E', 0)]}; survivor reads in net {by[('S', 1)]})")


def labels(V):
    return {r["read"]: r for r in csv.DictReader(open(f"{W0}/{V}/labels.tsv"), delimiter="\t")}


def stage_cluster(V):
    d = f"{W0}/{V}"
    lab = labels(V)
    seqs = SD.read_fa(f"{d}/pool.fa")
    truth = {n: lab[n]["copy"] for n in seqs}
    lens = {n: len(s) for n, s in seqs.items()}
    os.makedirs(f"{d}/rounds", exist_ok=True)
    try:
        comp = SD.run_rounds(lens, SD.minimap_map_fn(seqs, f"{d}/rounds", DELTA, 0.5), n_seeds=2000, min_size=3, max_rounds=6)
    except SD.Pause as e:
        print(f"round {e} finished; run again")
        sys.exit(75)
    cl = GR.clusters(comp, 3)
    maj = {c: S.majority(rs, truth) for c, rs in cl.items()}
    with open(f"{d}/clusters.tsv", "w") as o:
        o.write("read\tcluster\tsize\tmajority\n")
        for c, rs in cl.items():
            for r in rs:
                o.write(f"{r}\t{c}\t{len(rs)}\t{maj[c] or ''}\n")
    subprocess.run([sys.executable, f"{HERE}/consensus.py", f"{d}/pool.fa", f"{d}/clusters.tsv", f"{d}/cons.fa"], check=True)
    print(f"{V}: {len(cl)} clusters, {sum(len(r) for r in cl.values())} of {len(seqs)} net reads")


def cons_seqs(V):
    return RA_read_cons(f"{W0}/{V}/cons.fa")


def RA_read_cons(path):
    import run_polish as RP
    return RP.read_cons(path)


def cluster_info(V):
    """{cluster key 'clX': dict(majority copy, size, purity)}"""
    d = f"{W0}/{V}"
    lab = labels(V)
    cl = collections.defaultdict(list)
    for r in csv.DictReader(open(f"{d}/clusters.tsv"), delimiter="\t"):
        cl["cl" + r["cluster"]].append(r["read"])
    out = {}
    for k, rs in cl.items():
        cnt = collections.Counter(lab[r]["copy"] for r in rs)
        m, n = cnt.most_common(1)[0]
        out[k] = dict(majority=m, size=len(rs), purity=n / len(rs), family=copies()[m]["family"], role=copies()[m]["role"])
    return out


def stage_truth(V, cons_file="cons.fa", out="truth"):
    d = f"{W0}/{V}"
    info = cluster_info(V)
    cons = RA_read_cons(f"{d}/{cons_file}")
    tx = SD.read_fa(f"{W0}/transcripts.fa")
    cp = copies()
    res, oriented = {}, {}
    for name, seq in cons.items():
        k = name.split("|")[0]
        t = tx[info[k]["majority"]]
        a, b = E.compare(seq, t), E.compare(rc(seq), t)
        best, s2 = (a, seq) if a["identity"] >= b["identity"] else (b, rc(seq))   # orientation chosen by the global identity
        oriented[name] = s2
        res[k] = dict(best, orientation="+" if s2 is seq else "-")
    SD.write_fa(f"{d}/cons.oriented.{out}.fa", oriented, list(oriented))
    paf = f"{d}/cons.{out}.paf"
    subprocess.run(f"minimap2 -c --cs -x splice:hq -uf -N 5 -t 4 {W0}/genome.truth.fa {d}/cons.oriented.{out}.fa > {paf}", shell=True, check=True)
    best = {}
    for ln in open(paf):
        f = ln.rstrip("\n").split("\t")
        q = f[0].split("|")[0]
        if q not in best or int(f[9]) > best[q]["m"]:
            best[q] = dict(m=int(f[9]), ident=int(f[9]) / max(1, int(f[10])), cov=(int(f[3]) - int(f[2])) / int(f[1]), ref=f[5], start=int(f[7]), end=int(f[8]),
                           clip5=int(f[2]), clip3=int(f[1]) - int(f[3]), nm=int([t for t in f[12:] if t.startswith("NM:i:")][0][5:]))
    for k, r in res.items():
        c = cp[info[k]["majority"]]
        h = best.get(k)
        on = bool(h) and h["ref"] == c["contig"] and h["start"] < int(c["end"]) and int(c["pos0"]) < h["end"]
        r.update(on_true_copy=on, idcov=(h["ident"] * h["cov"]) if on else None, clip5=h["clip5"] if on else None, clip3=h["clip3"] if on else None, nm=h["nm"] if on else None)
    json.dump(res, open(f"{d}/{out}.json", "w"), indent=1)
    print(f"{V}: compared {len(res)} consensus sequences ({cons_file}) with their true transcripts and genome copies -> {out}.json")


def stage_attribute(V):
    d = f"{W0}/{V}"
    if not os.path.exists(f"{W0}/targets_db.nhr"):
        subprocess.run(f"{BLAST}/makeblastdb -in {W0}/targets.fa -dbtype nucl -out {W0}/targets_db", shell=True, check=True, capture_output=True)
    out = f"{d}/cons.blastn.tsv"
    subprocess.run(f"{BLAST}/blastn -task dc-megablast -query {d}/cons.fa -db {W0}/targets_db -evalue 1e-5 -num_threads 4 -max_target_seqs 5000 "
                   f"-outfmt '6 qseqid sseqid bitscore evalue length pident qstart qend' -out {out}", shell=True, check=True)
    hs = list(A.read_blastn(out))
    assert not A.capped_queries(hs, 5000)
    att = A.attribute_cover(A.cover_scores(hs), 1.10)
    info = cluster_info(V)
    res = {k: (att[k][0] if k in att else None) for k in info}
    json.dump(res, open(f"{d}/attribution.json", "w"), indent=1)
    print(f"{V}: attributed {sum(v is not None for v in res.values())} of {len(res)} clusters")


def stage_augment(V):
    d = f"{W0}/{V}"
    lab = labels(V)
    genome = json.load(open(f"{d}/genome.json"))
    sam = f"{d}/aug.sam"
    if not os.path.exists(sam + ".done"):
        subprocess.run(f"minimap2 {RA.MM2} -t 4 {d}/cons.fa {d}/reads.fa > {sam}", shell=True, check=True)
        open(sam + ".done", "w").write("ok")
    cons = RA.primaries(open(sam))
    info = cluster_info(V)
    att = json.load(open(f"{d}/attribution.json"))
    cf_att = {n: att.get(n.split("|")[0]) for n in cons_seqs(V)}
    cf_true = {n: info[n.split("|")[0]]["family"] for n in cons_seqs(V)}
    rows = []
    for n, r in lab.items():
        cls = "S" if r["role"] == "S" else ("E_net" if r["in_net"] == "1" else "E_abs")
        c = cons.get(n)
        rows.append((cls, r["family"], genome.get(n), c, c["ref"] if c else None, r["D"]))
    res = {}
    for rule in ("score", "divergence"):
        res[rule] = {}
        for scope, cf in (("attributed", cf_att), ("true_majority", cf_true)):
            res[rule][scope] = {}
            for D in sorted({r[5] for r in rows}, key=float):
                sub = [(a, b, c, e, f) for a, b, c, e, f, g in rows if g == D]
                res[rule][scope][D] = G.move_metrics(sub, cf, rule=rule)
    json.dump(res, open(f"{d}/augment.json", "w"), indent=1)
    print(f"{V}: augmented-reference moves written")


def stage_trim(V):
    """rule T (Amendment 13): trim a 1-3 bp pure-G 5' prefix not covered by the attributed family's best surviving copy"""
    import trim5g as T
    d = f"{W0}/{V}"
    att = json.load(open(f"{d}/attribution.json"))
    hs = collections.defaultdict(list)
    for q, t, a, b in A.read_blastn(f"{d}/cons.blastn.tsv"):
        hs[q].append((q, t, a, b))
    cons = cons_seqs(V)
    out, log = {}, {}
    for name, seq in cons.items():
        k = name.split("|")[0]
        fam = att.get(k)
        if fam is None:
            out[name], log[k] = seq, dict(trimmed=0, reason="abstained", prefix=None)
            continue
        prefix = T.best_copy_prefix(hs.get(k, []), fam)
        new, n, why = T.decide(seq, prefix)
        out[name] = new
        log[k] = dict(trimmed=n, reason=why, prefix=prefix)
    SD.write_fa(f"{d}/cons.trimmed.fa", out, list(out))
    json.dump(log, open(f"{d}/trim.json", "w"), indent=1)
    print(f"{V}: " + ", ".join(f"{r} {c}" for r, c in collections.Counter(v["reason"] for v in log.values()).most_common()))


def stage_trim2(V):
    """rule T2 (Amendment 14): the 5' prefix is the consensus start of its best alignment to the attributed family's surviving copies"""
    import trim5g as T
    d = f"{W0}/{V}"
    att = json.load(open(f"{d}/attribution.json"))
    cons = cons_seqs(V)
    paf = f"{d}/cons.targets.paf"
    subprocess.run(f"minimap2 -c -x splice:hq -uf -N 5 -t 4 {W0}/targets.fa {d}/cons.fa > {paf}", shell=True, check=True, stderr=subprocess.DEVNULL)
    fams = {n.split("|")[0]: att.get(n.split("|")[0]) for n in cons}
    pref = T.prefixes_from_paf(open(paf).read().splitlines(), fams)
    out, log = {}, {}
    for name, seq in cons.items():
        k = name.split("|")[0]
        if fams[k] is None:
            out[name], log[k] = seq, dict(trimmed=0, reason="abstained", prefix=None)
            continue
        new, n, why = T.decide(seq, pref.get(k))
        out[name], log[k] = new, dict(trimmed=n, reason=why, prefix=pref.get(k))
    SD.write_fa(f"{d}/cons.trimmed2.fa", out, list(out))
    json.dump(log, open(f"{d}/trim2.json", "w"), indent=1)
    print(f"{V}: " + ", ".join(f"{r} {c}" for r, c in collections.Counter(v["reason"] for v in log.values()).most_common()))


def erased_clusters(V):
    """per erased copy the largest pure cluster of >= 3 reads: {copy: cluster key}"""
    info = cluster_info(V)
    best = {}
    for k, i in info.items():
        if i["role"] == "E" and i["purity"] == 1.0 and i["size"] >= 3 and (i["majority"] not in best or i["size"] > info[best[i["majority"]]]["size"]):
            best[i["majority"]] = k
    return best


def stage_trimreport(V, tag="trimmed", trimjson="trim.json"):
    d = f"{W0}/{V}"
    cp = copies()
    tr0, tr1 = json.load(open(f"{d}/truth.json")), json.load(open(f"{d}/truth.{tag}.json"))
    trim = json.load(open(f"{d}/{trimjson}"))
    best = erased_clusters(V)
    out = {}
    for D in sorted({r["D"] for r in cp.values()}, key=float):
        keys = [k for c, k in best.items() if cp[c]["D"] == D]
        row = dict(clusters=len(keys))
        for tag, tr in (("untrimmed", tr0), ("trimmed", tr1)):
            ic = [tr[k]["idcov"] for k in keys if tr[k]["idcov"] is not None]
            row[tag] = dict(idcov_ge_0_999=sum(x >= 0.999 for x in ic), idcov_n=len(ic), left_in_window=sum(-4 <= tr[k]["left"] <= 1 for k in keys if tr[k]["left"] is not None),
                            left_below_minus4=sum(tr[k]["left"] < -4 for k in keys if tr[k]["left"] is not None), lefts=sorted(tr[k]["left"] for k in keys if tr[k]["left"] is not None))
        row["decisions"] = dict(collections.Counter(trim[k]["reason"] for k in keys))
        out[D] = row
    json.dump(out, open(f"{d}/trimreport.{tag}.json", "w"), indent=1)
    print(f"{'D':>6} {'n':>3} | untrimmed idcov>=.999 | trimmed idcov>=.999 | left in [-4,+1] (before/after) | over-trim (<-4, after) | decisions")
    tot = collections.Counter()
    for D, r in out.items():
        u, t = r["untrimmed"], r["trimmed"]
        print(f"{D:>6} {r['clusters']:>3} | {u['idcov_ge_0_999']:>3}/{u['idcov_n']:<3}            | {t['idcov_ge_0_999']:>3}/{t['idcov_n']:<3}          | {u['left_in_window']:>3} / {t['left_in_window']:<3}                  | {t['left_below_minus4']:>3}                   | {r['decisions']}")
        if float(D) >= 0.01:
            for key, v in (("n", r["clusters"]), ("u_ok", u["idcov_ge_0_999"]), ("t_ok", t["idcov_ge_0_999"]), ("u_win", u["left_in_window"]), ("t_win", t["left_in_window"]), ("t_over", t["left_below_minus4"]),
                           ("trimmed", r["decisions"].get("trimmed", 0))):
                tot[key] += v
    print("D >= 1% totals:", dict(tot))


def stage_report(V):
    d = f"{W0}/{V}"
    lab = labels(V)
    info, att = cluster_info(V), json.load(open(f"{d}/attribution.json"))
    tr = json.load(open(f"{d}/truth.json"))
    aug = json.load(open(f"{d}/augment.json"))["score"]["attributed"]
    augd = json.load(open(f"{d}/augment.json"))["divergence"]["attributed"]
    cp = copies()
    classes = sorted({r["D"] for r in lab.values()}, key=float)
    out = {}
    for D in classes:
        erased = [c for c, r in cp.items() if r["role"] == "E" and r["D"] == D]
        ereads = [r for r in lab.values() if r["role"] == "E" and r["D"] == D]
        # clusters whose majority copy is an erased copy of this class, one per copy (the largest)
        best = {}
        for k, i in info.items():
            if i["majority"] in erased and (i["majority"] not in best or i["size"] > info[best[i["majority"]]]["size"]):
                best[i["majority"]] = k
        s2 = [k for k in best.values() if info[k]["purity"] == 1.0 and info[k]["size"] >= 3]
        ident = [tr[k]["core_identity"] for k in s2 if tr[k]["core_identity"] is not None]
        idcov = [tr[k]["idcov"] for k in s2 if tr[k]["idcov"] is not None]
        e_net = [r for r in ereads if r["in_net"] == "1"]
        in_cl = {r for k in s2 for r in []}
        # S5: erased-copy net reads in a pure cluster that is attributed to the right family
        cl_reads = collections.defaultdict(list)
        for r in csv.DictReader(open(f"{d}/clusters.tsv"), delimiter="\t"):
            cl_reads["cl" + r["cluster"]].append(r["read"])
        right = wrong = 0
        for k, rs in cl_reads.items():
            a = att.get(k)
            for r in rs:
                if lab[r]["role"] != "E" or lab[r]["D"] != D:
                    continue
                if a is None:
                    continue
                right += a == lab[r]["family"]
                wrong += a != lab[r]["family"]
        m = aug[D]
        out[D] = dict(erased_copies=len(erased), erased_reads=len(ereads), erased_reads_in_net=len(e_net), frac_in_net=len(e_net) / max(1, len(ereads)),
                      copies_with_pure_cluster=len(s2), consensus_identity_ge_0_999=sum(x >= 0.999 for x in ident), idcov_ge_0_999=sum(x >= 0.999 for x in idcov),
                      idcov_n=len(idcov), s5_right=right, s5_wrong=wrong, s5_right_frac_of_net=right / max(1, len(e_net)),
                      s6_E_net_moved_own=m.get("E_net", {}).get("moved_to_own_family"), s6_E_net_n=m.get("E_net", {}).get("n"),
                      s6_S_moved=m.get("S", {}).get("moved"), s6_S_n=m.get("S", {}).get("n"),
                      s6d_E_net_moved_own=augd[D].get("E_net", {}).get("moved_to_own_family"), s6d_E_abs_moved_own=augd[D].get("E_abs", {}).get("moved_to_own_family"),
                      s6d_E_abs_n=augd[D].get("E_abs", {}).get("n"), s6d_S_moved=augd[D].get("S", {}).get("moved"),
                      ends_left=[tr[k]["left"] for k in s2], ends_right=[tr[k]["right"] for k in s2],
                      clip5=[tr[k]["clip5"] for k in s2 if tr[k]["clip5"] is not None], clip3=[tr[k]["clip3"] for k in s2 if tr[k]["clip3"] is not None])
    json.dump(out, open(f"{d}/report.json", "w"), indent=1)
    print(f"{'D':>6} {'copies':>6} {'inNet':>7} {'pure':>5} {'id>=.999':>8} {'idcov':>7} {'S5 right/net':>13} {'S5 wrong':>8} {'S6 E moved':>11} {'S6 S moved':>11}")
    for D, o in out.items():
        print(f"{D:>6} {o['erased_copies']:>6} {o['frac_in_net']:>7.3f} {o['copies_with_pure_cluster']:>5} {o['consensus_identity_ge_0_999']:>8} {o['idcov_ge_0_999']:>4}/{o['idcov_n']:<2} "
              f"{o['s5_right']:>5}/{o['erased_reads_in_net']:<7} {o['s5_wrong']:>8} {o['s6_E_net_moved_own']!s}/{o['s6_E_net_n']!s:>5} {o['s6_S_moved']}/{o['s6_S_n']}   [divergence rule: E_net own {o['s6d_E_net_moved_own']}, E_abs own {o['s6d_E_abs_moved_own']}/{o['s6d_E_abs_n']}, S moved {o['s6d_S_moved']}]")
    for D, o in out.items():
        print(D, "ends 5' offsets", sorted(o["ends_left"]), "3' offsets", sorted(o["ends_right"]), "genome clips 5'/3' medians",
              (statistics.median(o["clip5"]) if o["clip5"] else None), (statistics.median(o["clip3"]) if o["clip3"] else None))


if __name__ == "__main__":
    {"map": stage_map, "cluster": stage_cluster, "truth": stage_truth, "attribute": stage_attribute, "augment": stage_augment, "report": stage_report,
     "trim": stage_trim, "trim2": stage_trim2, "trimreport": stage_trimreport}[sys.argv[2]](sys.argv[1], *sys.argv[3:])
