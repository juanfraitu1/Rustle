"""score: every stage's output against the planted truth -> DIR/score.tsv (long), DIR/summary.md, DIR/reads.placement.tsv.

Matchings are one-to-one (loci <-> copies by shared exonic bp, scipy linear_sum_assignment; metric trap: two truth
copies must never share one node). Read attribution is by the planted interval of the copy the primary alignment lies
in (copies never overlap here, so interval = exon attribution).

Objectives
  alignment  what minimap2 did with the reads (reported, not judged)
  o1_loci    assembled loci vs copies: found / split / merged / missed; junction P/R of the best transcript; exact chain
  o1_family  de novo (run.fam.clusters.tsv) and guided (guided.clusters.tsv) families: pairwise S/P/F, bipartite F,
             family_intact, guided_minus_denovo
  o2         run.assign.assignments.tsv: correct / wrong / abstain per read, tied vs unique, per copy
  o3         run.flag.missing_copy.tsv: verdict at the locus that absorbed each reference-absent copy's reads
"""
import collections
import os
import re

from lib import bipartite_families, merge, pairwise  # noqa: E402

from .chromosome import load_planted
from .reads import load_truth


class Score:
    def __init__(self, scenario):
        self.scenario = scenario
        self.rows = []

    def add(self, objective, metric, copy, value, note=""):
        self.rows.append((self.scenario, objective, metric, copy, value if isinstance(value, str) else f"{value:.4f}" if isinstance(value, float) else str(value), note))

    def write(self, path):
        with open(path, "w") as fh:
            fh.write("scenario\tobjective\tmetric\tcopy\tvalue\tnote\n")
            for r in self.rows:
                fh.write("\t".join(r) + "\n")

    def get(self, objective, metric, copy="ALL"):
        for r in self.rows:
            if r[1] == objective and r[2] == metric and r[3] == copy:
                return r[4]
        return ""


def _ov(a0, a1, b0, b1):
    return max(0, min(a1, b1) - max(a0, b0))


# ---------------------------------------------------------------- alignment
def alignment_report(out_dir, planted, truth, S):
    import pysam
    byid = {p.id: p for p in planted}
    iv = [(p.contig, p.pos, p.end, p.id) for p in planted]

    def copy_at(contig, s, e):
        best = None
        for c, a, b, cid in iv:
            if c == contig:
                o = _ov(s, e, a, b)
                if o > 0 and (best is None or o > best[1]):
                    best = (cid, o)
        return best[0] if best else None

    recs = collections.defaultdict(list)
    bam = os.path.join(out_dir, "reads.bam")
    for r in pysam.AlignmentFile(bam):
        recs[r.query_name].append(r)
    rows = []
    per = collections.defaultdict(collections.Counter)
    jrec = collections.defaultdict(lambda: [0, 0])
    clip = collections.defaultdict(list)
    for name, t in truth.items():
        rs = recs.get(name, [])
        prim = [r for r in rs if not r.is_secondary and not r.is_supplementary and not r.is_unmapped]
        secs = [r for r in rs if r.is_secondary]
        sups = [r for r in rs if r.is_supplementary]
        cid = t["copy"]
        if not prim:
            rows.append((name, cid, "unmapped", "", -1, 0, 0, 0, 0.0, 0, 0)); per[cid]["unmapped"] += 1
            continue
        p = prim[0]
        at = copy_at(p.reference_name, p.reference_start, p.reference_end)
        if at is None:
            cls = "outside"
        elif at == cid:
            cls = "own"
        elif not byid[cid].in_reference:
            cls = f"absorbed_by:{at}"
        else:
            cls = f"other:{at}"
        pas = p.get_tag("AS") if p.has_tag("AS") else 0
        tied = int(p.mapping_quality == 0 or any((r.get_tag("AS") if r.has_tag("AS") else 0) >= 0.98 * pas for r in secs))
        introns = set()
        rp = p.reference_start
        for op, L in p.cigartuples:
            if op in (0, 7, 8, 2):
                rp += L
            elif op == 3:
                introns.add((rp, rp + L)); rp += L
        sc = sum(L for op, L in p.cigartuples if op == 4)
        big_ins = max([L for op, L in p.cigartuples if op == 1] or [0])
        jr = len(t["junctions"] & introns) if cls == "own" else -1
        if cls == "own":
            jrec[cid][0] += jr; jrec[cid][1] += len(t["junctions"])
        clip[cid].append(sc / max(1, p.query_length + sc))
        per[cid][cls.split(":")[0]] += 1
        per[cid]["mapq0"] += p.mapping_quality == 0
        per[cid]["tied"] += tied
        per[cid]["supplementary"] += bool(sups)
        per[cid]["softclip50"] += sc >= 50
        per[cid]["insertion50"] += big_ins >= 50
        rows.append((name, cid, cls, f"{p.reference_name}:{p.reference_start}-{p.reference_end}", p.mapping_quality, tied, len(secs),
                     len(sups), sc, big_ins, jr, len(t["junctions"])))
    with open(os.path.join(out_dir, "reads.placement.tsv"), "w") as fh:
        fh.write("read\tcopy\tclass\tprimary\tmapq\ttied\tn_secondary\tn_supplementary\tsoftclip\tmax_insertion\tjunctions_recovered\tjunctions_true\n")
        for r in rows:
            fh.write("\t".join(str(x) for x in r) + "\n")
    for p in planted:
        n = sum(1 for t in truth.values() if t["copy"] == p.id)
        if n == 0:
            continue
        c = per[p.id]
        S.add("alignment", "reads", p.id, n)
        for k in ("own", "other", "absorbed_by", "outside", "unmapped", "mapq0", "tied", "supplementary", "softclip50", "insertion50"):
            S.add("alignment", k, p.id, c[k] / n)
        if jrec[p.id][1]:
            S.add("alignment", "junction_recall", p.id, jrec[p.id][0] / jrec[p.id][1], "own-copy primaries only")
        S.add("alignment", "mean_clip_frac", p.id, sum(clip[p.id]) / len(clip[p.id]) if clip[p.id] else 0.0)
    return per


# ---------------------------------------------------------------- O1 loci
def gtf_loci(path):
    """{gene_id: {"contig", "strand", "exons": {(s0,e1)}, "transcripts": {tid: [(s0,e1)]}, "reads": {tid: n}}}"""
    loci = {}
    for line in open(path):
        f = line.rstrip("\n").split("\t")
        if len(f) < 9 or f[2] not in ("exon", "transcript"):
            continue
        gid = re.search(r'gene_id "([^"]+)"', f[8]).group(1)
        tid = re.search(r'transcript_id "([^"]+)"', f[8]).group(1)
        L = loci.setdefault(gid, {"contig": f[0], "strand": f[6], "exons": set(), "transcripts": {}, "reads": {}})
        if f[2] == "transcript":
            m = re.search(r'reads "(\d+)"', f[8])
            L["reads"][tid] = int(m.group(1)) if m else 0
        else:
            L["exons"].add((int(f[3]) - 1, int(f[4])))
            L["transcripts"].setdefault(tid, []).append((int(f[3]) - 1, int(f[4])))
    for L in loci.values():
        for t in L["transcripts"].values():
            t.sort()
    return loci


def match_loci(loci, planted):
    """One-to-one loci <-> reference-present copies by shared exonic bp. Returns (copy->locus, locus->copy, overlap)."""
    import numpy as np
    from scipy.optimize import linear_sum_assignment
    cands = [p for p in planted if p.in_reference]
    lids = sorted(loci)
    ovl = {}
    for i, p in enumerate(cands):
        ex = [(s, e) for s, e, *_ in p.dna_exon_intervals()]
        for j, l in enumerate(lids):
            L = loci[l]
            if L["contig"] != p.contig:
                continue
            o = sum(_ov(a, b, s, e) for a, b in ex for s, e in merge(sorted(L["exons"])))
            if o > 0:
                ovl[(p.id, l)] = o
    M = np.zeros((len(cands), len(lids)), dtype=int)
    for (cid, l), o in ovl.items():
        M[[k for k, p in enumerate(cands) if p.id == cid][0], lids.index(l)] = o
    c2l, l2c = {}, {}
    if M.size:
        r, c = linear_sum_assignment(-M)
        for i, j in zip(r, c):
            if M[i, j] > 0:
                c2l[cands[i].id] = lids[j]; l2c[lids[j]] = cands[i].id
    return c2l, l2c, ovl


def o1_loci(out_dir, planted, S):
    gtf = os.path.join(out_dir, "run.gtf")
    if not os.path.exists(gtf):
        return None
    loci = gtf_loci(gtf)
    c2l, l2c, ovl = match_loci(loci, planted)
    status = {}
    for p in planted:
        if not p.in_reference:
            continue
        exbp = sum(e - s for s, e, *_ in p.dna_exon_intervals())
        touching = [l for (cid, l), o in ovl.items() if cid == p.id and o >= 100]
        l = c2l.get(p.id)
        if p.expression == 0:
            st = "unexpressed"
        elif l is None or ovl.get((p.id, l), 0) < 100:
            st = "missed"
        else:
            others = [cid for (cid, l2), o in ovl.items() if l2 == l and cid != p.id and o >= 100]
            if others:
                st = "merged:" + ",".join(sorted(others))
            elif len(touching) >= 2:
                st = f"split:{len(touching)}"
            else:
                st = "found"
        status[p.id] = st
        S.add("o1_loci", "status", p.id, st)
        S.add("o1_loci", "n_loci_touching", p.id, len(touching))
        if l is not None:
            L = loci[l]
            S.add("o1_loci", "locus", p.id, l)
            S.add("o1_loci", "transcripts", p.id, len(L["transcripts"]))
            lex = sum(e - s for s, e in merge(sorted(L["exons"])))
            S.add("o1_loci", "exonic_recall", p.id, ovl[(p.id, l)] / exbp if exbp else 0.0, "shared exonic bp / copy's DNA exonic bp")
            S.add("o1_loci", "exonic_precision", p.id, ovl[(p.id, l)] / lex if lex else 0.0, "shared exonic bp / locus exonic bp")
            truth_j = p.junctions()
            best = None
            for tid, ex in L["transcripts"].items():
                j = {(a, b) for (_, a), (b, _) in zip(ex, ex[1:])}
                tp = len(j & truth_j)
                pr = tp / len(j) if j else (1.0 if not truth_j else 0.0)
                rc_ = tp / len(truth_j) if truth_j else (1.0 if not j else 0.0)
                f = 2 * pr * rc_ / (pr + rc_) if pr + rc_ else 0.0
                if best is None or (f, tp) > (best[0], best[1]):
                    best = (f, tp, pr, rc_, j == truth_j, tid)
            S.add("o1_loci", "junction_precision", p.id, best[2]); S.add("o1_loci", "junction_recall", p.id, best[3])
            S.add("o1_loci", "exact_chain", p.id, int(best[4]), best[5])
    n_exp = sum(1 for p in planted if p.in_reference and p.expression != 0)
    S.add("o1_loci", "found_frac", "ALL", sum(1 for s in status.values() if s == "found") / n_exp if n_exp else 0.0)
    S.add("o1_loci", "n_loci", "ALL", len(loci))
    return c2l, l2c, loci


# ---------------------------------------------------------------- O1 families
def read_clusters(path):
    """[(cluster_id, chrom, start1, end1)] from mcl_families clusters.tsv."""
    out = []
    for line in open(path):
        f = line.rstrip("\n").split("\t")
        if f[0] == "cluster_id" or len(f) < 8:
            continue
        out.append((f[0], f[5], int(f[6]), int(f[7])))
    return out


def o1_family(out_dir, planted, S, mode, c2l=None, loci=None):
    path = os.path.join(out_dir, "run.fam.clusters.tsv" if mode == "denovo" else "guided.clusters.tsv")
    if not os.path.exists(path):
        return None
    members = read_clusters(path)
    # universe: reference-present copies and decoys; de novo additionally needs expression
    uni = [p for p in planted if p.in_reference and (mode == "guided" or p.expression != 0)]
    label = {}
    for cid, chrom, s, e in members:
        best = None
        for p in uni:
            if p.contig != chrom:
                continue
            o = _ov(s - 1, e, p.pos, p.end)
            if o > 0 and (best is None or o > best[1]):
                best = (p, o)
        if best is None:
            S.add(f"o1_family_{mode}", "member_unmatched", f"{chrom}:{s}-{e}", cid, "a cluster member overlapping no planted copy")
            continue
        p = best[0]
        if p.id in label and label[p.id] != cid:
            S.add(f"o1_family_{mode}", "copy_in_two_clusters", p.id, f"{label[p.id]},{cid}")
        label[p.id] = cid
    pred = [label.get(p.id, f"single:{p.id}") for p in uni]
    true = [p.family for p in uni]
    for p in uni:
        S.add(f"o1_family_{mode}", "cluster", p.id, label.get(p.id, "singleton"))
    sens, prec = pairwise(pred, true)
    f = 2 * sens * prec / (sens + prec) if sens == sens and prec == prec and sens + prec else 0.0
    S.add(f"o1_family_{mode}", "pair_sensitivity", "ALL", sens); S.add(f"o1_family_{mode}", "pair_precision", "ALL", prec)
    S.add(f"o1_family_{mode}", "pair_F", "ALL", f)
    truth_fams = collections.defaultdict(list); pred_fams = collections.defaultdict(list)
    for p, t, pr in zip(uni, true, pred):
        truth_fams[t].append(p.id)
        if not pr.startswith("single:"):
            pred_fams[pr].append(p.id)
    tf = {k: v for k, v in truth_fams.items() if len(v) >= 2}
    bp = bipartite_families(tf, dict(pred_fams)) if tf else None
    if bp:
        for k, v in bp.items():
            S.add(f"o1_family_{mode}", "bipartite_F", k, v["f"], f"sens {v['sens']} prec {v['prec']} cluster {v['cluster']}")
    a_copies = [p for p in uni if p.family == "A"]
    cl = {label.get(p.id) for p in a_copies}
    intact = len(a_copies) >= 2 and len(cl) == 1 and None not in cl and all(label.get(d.id) not in cl for d in uni if d.kind == "decoy")
    S.add(f"o1_family_{mode}", "family_intact", "ALL", "n/a" if len(a_copies) < 2 else int(intact),
          f"{len(a_copies)} copies of A in {len(cl)} cluster(s)" + ("" if len(a_copies) >= 2 else " (fewer than 2 scoreable copies in this mode)"))
    S.add(f"o1_family_{mode}", "n_clusters", "ALL", len({m[0] for m in members}))
    return f


# ---------------------------------------------------------------- O2
def o2(out_dir, planted, truth, placement, S):
    path = os.path.join(out_dir, "run.assign.assignments.tsv")
    if not os.path.exists(path):
        return None
    byid = {p.id: p for p in planted}
    # catalog copy -> planted copy (the families copy table rows are locus representatives)
    cat = {}
    tab = None
    for cand in ("run.fam.copies.tsv", "run.cat.copies.tsv"):
        if os.path.exists(os.path.join(out_dir, cand)):
            tab = os.path.join(out_dir, cand); break
    if tab:
        with open(tab) as fh:
            hdr = fh.readline().rstrip("\n").split("\t")
            for line in fh:
                f = dict(zip(hdr, line.rstrip("\n").split("\t")))
                s, e = int(f["start"]), int(f["end"])
                best = None
                for p in planted:
                    if p.contig == f["chrom"]:
                        o = _ov(s, e, p.pos, p.end)
                        if o > 0 and (best is None or o > best[1]):
                            best = (p.id, o)
                cat[(f["family_id"], f["copy_idx"])] = best[0] if best else None
    rows = collections.defaultdict(list)
    with open(path) as fh:
        hdr = fh.readline().rstrip("\n").split("\t")
        for line in fh:
            f = dict(zip(hdr, line.rstrip("\n").split("\t")))
            rows[f["read_name"]].append(f)
    tied_of, class_of = {}, {}
    for line in open(os.path.join(out_dir, "reads.placement.tsv")):
        f = line.rstrip("\n").split("\t")
        if f[0] != "read":
            tied_of[f[0]] = f[5] == "1"; class_of[f[0]] = f[2]
    per = collections.defaultdict(collections.Counter)
    for name, t in truth.items():
        cid = t["copy"]
        rs = rows.get(name)
        tied = tied_of.get(name, False)
        band = "tied" if tied else "unique"
        if not rs or (not tied and not any(r["status"] == "assigned" for r in rs)):
            # unique molecules are skipped by the AS-tied gate: the aligner's placement is the answer (sim.py tandem's rule)
            if tied:
                per[cid]["tied_no_row"] += 1
            elif class_of.get(name) == "own":
                per[cid]["unique_placed_correct"] += 1
            elif class_of.get(name, "").startswith("absorbed_by"):
                per[cid]["unique_placed_absent"] += 1
            elif class_of.get(name) == "unmapped":
                per[cid]["unique_unmapped"] += 1
            else:
                per[cid]["unique_placed_wrong"] += 1
            continue
        assigned = [r for r in rs if r["status"] == "assigned" and r.get("origin_rejected", "0") == "0"]
        if not assigned:
            per[cid][f"{band}_abstain"] += 1; continue
        targets = {cat.get((r["family_id"], r["catalog_copy_idx"])) for r in assigned}
        if not byid[cid].in_reference:
            per[cid][f"{band}_absent_assigned"] += 1          # the true copy is not in the roster: any assignment is wrong
        elif targets == {cid}:
            per[cid][f"{band}_correct"] += 1
        else:
            per[cid][f"{band}_wrong"] += 1
    tot = collections.Counter()
    for p in planted:
        n = sum(1 for t in truth.values() if t["copy"] == p.id)
        if n == 0:
            continue
        c = per[p.id]
        for k in sorted(c):
            S.add("o2", k, p.id, c[k] / n); tot[k] += c[k]
        tot["n"] += n
    n = tot["n"] or 1
    for k in ("tied_correct", "tied_wrong", "tied_abstain", "tied_absent_assigned", "unique_placed_correct", "unique_placed_wrong",
              "unique_placed_absent", "unique_correct", "unique_wrong"):
        S.add("o2", k, "ALL", tot[k] / n)
    S.add("o2", "tied_frac", "ALL", sum(v for k, v in tot.items() if k.startswith("tied_")) / n)
    ass = sum(v for k, v in tot.items() if k.endswith("_correct") or k.endswith("_wrong") or k.endswith("_absent_assigned"))
    wrong = sum(v for k, v in tot.items() if k.endswith("_wrong") or k.endswith("_absent_assigned"))
    S.add("o2", "wrong_among_assigned", "ALL", wrong / ass if ass else 0.0, "assigned = certificate assignments + unique placements")
    return tot


# ---------------------------------------------------------------- O3
def o3(out_dir, planted, placement, S):
    path = os.path.join(out_dir, "run.flag.missing_copy.tsv")
    if not os.path.exists(path):
        return None
    byid = {p.id: p for p in planted}
    flags = []
    with open(path) as fh:
        hdr = fh.readline().rstrip("\n").split("\t")
        for line in fh:
            f = dict(zip(hdr, line.rstrip("\n").split("\t")))
            flags.append(f)
    hosts = {}
    for p in planted:
        if p.in_reference or p.expression == 0:
            continue
        absorb = collections.Counter()
        for line in open(os.path.join(out_dir, "reads.placement.tsv")):
            f = line.rstrip("\n").split("\t")
            if f[0] != "read" and f[1] == p.id and f[2].startswith("absorbed_by:"):
                absorb[f[2].split(":")[1]] += 1
        hosts[p.id] = absorb.most_common(1)[0][0] if absorb else None
        S.add("o3", "absorbing_copy", p.id, hosts[p.id] or "none")
    for p in planted:
        if p.expression == 0:
            continue
        target = byid[hosts[p.id]] if (not p.in_reference and hosts.get(p.id)) else (p if p.in_reference else None)
        if target is None:
            S.add("o3", "verdict", p.id, "no_host"); continue
        expect_flag = (not p.in_reference) or (p.id in hosts.values())
        hit = [f for f in flags if f["chrom"] == target.contig and _ov(int(f["start"]), int(f["end"]), target.pos, target.end) > 0]
        v = hit[0]["verdict"] if hit else "not_flagged"
        S.add("o3", "verdict", p.id, v, f"at {target.id}'s locus; expected {'reference_absent_candidate' if expect_flag else 'not_flagged'}"
              + (" (hosts an absent copy's reads)" if p.in_reference and expect_flag else ""))
        S.add("o3", "flag_correct", p.id, int((v == "reference_absent_candidate") == expect_flag))
        if hit:
            for k in ("class", "m", "delta", "n_psv", "shared_frac", "conf_truth_identity", "conf_truth"):
                if k in hit[0]:
                    S.add("o3", k, p.id, hit[0][k])
    return flags


# ---------------------------------------------------------------- summary
def run(out_dir, log=print):
    planted, man = load_planted(out_dir)
    truth = load_truth(out_dir)
    S = Score(man["spec"].get("name", os.path.basename(out_dir.rstrip("/"))))
    placement = alignment_report(out_dir, planted, truth, S) if os.path.exists(os.path.join(out_dir, "reads.bam")) else None
    res = o1_loci(out_dir, planted, S)
    c2l, loci = (res[0], res[2]) if res else (None, None)
    fd = o1_family(out_dir, planted, S, "denovo", c2l, loci)
    fg = o1_family(out_dir, planted, S, "guided")
    if fd is not None and fg is not None:
        S.add("o1_family", "guided_minus_denovo_F", "ALL", fg - fd)
    if placement is not None:
        o2(out_dir, planted, truth, placement, S)
        o3(out_dir, planted, placement, S)
    S.write(os.path.join(out_dir, "score.tsv"))
    write_summary(out_dir, planted, S)
    log(f"score: {len(S.rows)} rows -> score.tsv, summary.md")
    return S


def headline(S):
    """The ladder's one-row view of a scenario."""
    return {
        "denovo_intact": S.get("o1_family_denovo", "family_intact"), "guided_intact": S.get("o1_family_guided", "family_intact"),
        "denovo_F": S.get("o1_family_denovo", "pair_F"), "guided_F": S.get("o1_family_guided", "pair_F"),
        "loci_found": S.get("o1_loci", "found_frac"),
        "o2_wrong": S.get("o2", "wrong_among_assigned"), "o2_tied_correct": S.get("o2", "tied_correct"), "o2_tied_abstain": S.get("o2", "tied_abstain"),
        "o3": ";".join(r[4] for r in S.rows if r[1] == "o3" and r[2] == "verdict") or "",
    }


def write_summary(out_dir, planted, S):
    lines = [f"# {S.scenario}\n"]
    objs = []
    for r in S.rows:
        if r[1] not in objs:
            objs.append(r[1])
    for o in objs:
        rows = [r for r in S.rows if r[1] == o]
        metrics = []
        for r in rows:
            if r[2] not in metrics:
                metrics.append(r[2])
        copies = []
        for r in rows:
            if r[3] not in copies:
                copies.append(r[3])
        lines.append(f"\n## {o}\n")
        lines.append("| metric | " + " | ".join(copies) + " |")
        lines.append("|---|" + "---|" * len(copies))
        for m in metrics:
            vals = {r[3]: r[4] for r in rows if r[2] == m}
            lines.append(f"| {m} | " + " | ".join(vals.get(c, "") for c in copies) + " |")
        notes = {r[5] for r in rows if r[5]}
        if notes and len(notes) <= 6:
            lines.append("\n" + "; ".join(sorted(notes)))
    with open(os.path.join(out_dir, "summary.md"), "w") as fh:
        fh.write("\n".join(lines) + "\n")
