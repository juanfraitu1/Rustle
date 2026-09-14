#!/usr/bin/env python3
"""Prereg Addendum Y: a per-copy mechanism layer. Each family member is classified against the family PARENT (the
protein-coding member with the most introns) as GENOMIC (shares >= 1 kb non-exonic, non-repeat-masked sequence with the
parent), RETROCOPY (parent introns lost at exon-exon junctions plus a 3' poly(A) and/or target-site duplication), or
UNRESOLVED; independently, a protein-coding gene is TE-DERIVED when >= 50% of its CDS lies in RepeatMasker
LINE/SINE/LTR/Retroposon records.

usage: copy_mechanism.py --gff full.gff.gz --genome softmasked.fa --rmsk hs1.repeatMasker.out.gz --outdir DIR
                         [--amy-truth truth.tsv] [--families A,B] [--threads 2]
The parsed GFF, RepeatMasker intervals and per-family alignments are cached in --outdir.
"""
import argparse
import collections
import csv
import gzip
import itertools
import os
import pickle
import random
import re
import subprocess
import sys

import pysam

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import duplicon_origin as do  # noqa: E402

MIN_ID, NONEXONIC, CAP = 0.80, 1000, 40
MIN_INTRON, JUNC_IN, LOST_MAX, RETAINED_MIN = 70, 20, 30, 70
POLYA_W, POLYA_MIN, POLYA_REACH = 15, 12, 100
TSD_MIN, TSD_MAXFRAC, CHANCE_MAX = 10, 0.60, 0.10
RETRO_CLASSES = {"LINE", "SINE", "LTR", "Retroposon"}

FAMILIES = {
    # RETRO positives (held out)
    "RPL7": (r"^RPL7(P\d+)?$", "retro_pos"), "RPL23A": (r"^RPL23A(P\d+)?$", "retro_pos"),
    "HNRNPA1": (r"^HNRNPA1(P\d+)?$", "retro_pos"), "PTMA": (r"^PTMA(P\d+)?$", "retro_pos"),
    "TPT1": (r"^TPT1(P\d+)?$", "retro_pos"), "NACA": (r"^NACA(P\d+)?$", "retro_pos"),
    "KRT8": (r"^KRT8(P\d+)?$", "retro_pos"), "FTH1": (r"^FTH1(P\d+)?$", "retro_pos"),
    "PGK": (r"^PGK[12]$", "retro_pos"), "GLUD": (r"^GLUD[12]$", "retro_pos"), "UTP14": (r"^UTP14[AC]$", "retro_pos"),
    "TAF1L": (r"^TAF1L?$", "retro_pos"), "RPL10L": (r"^RPL10L?$", "retro_pos"), "PABPC": (r"^PABPC[13]$", "retro_pos"),
    "POU5F1": (r"^POU5F1B?$", "retro_pos"), "NANOG": (r"^NANOG(P8)?$", "retro_pos"),
    # SD negatives (held out for this layer)
    "NOTCH2NL": (r"^NOTCH2NL[A-Z]?$|^NOTCH2$", "sd_neg"), "SRGAP2": (r"^SRGAP2[A-D]?$", "sd_neg"),
    "ARHGAP11": (r"^ARHGAP11[AB]$", "sd_neg"), "HYDIN": (r"^HYDIN\d?$", "sd_neg"), "FAM72": (r"^FAM72[A-D]$", "sd_neg"),
    "SMN": (r"^SMN[12]$", "sd_neg"), "SERF1": (r"^SERF1[AB]$", "sd_neg"), "ZNG1": (r"^ZNG1[A-F]$", "sd_neg"),
    "NCF1": (r"^NCF1[BC]?$", "sd_neg"), "DEFB4": (r"^DEFB4[AB]$", "sd_neg"), "GTF2H2": (r"^GTF2H2C?$", "sd_neg"),
    "SPDYE": (r"^SPDYE\d+[A-Z]?$", "sd_neg"), "USP17L": (r"^USP17L\d+$", "sd_neg"),
    # development
    **{f: (p, "dev_" + k) for f, (p, k) in do.FAMILIES.items() if k in ("core", "retro", "retro_ho")},
}
TE_GROUPS = {
    "ERV-env": (r"^ERVW-1$|^ERVFRD-1$|^ERVV-[12]$|^ERVH48-1$|^ERVMER34-1$|^ERVK3-1$|^ERV3-1$", "deciding"),
    "Ty3/gypsy": (r"^PEG10$|^RTL\d+[A-Z]?$|^ARC$|^ASPRV1$|^NYNRIN$", "reported"),
    "PNMA": (r"^PNMA\d+[A-Z]?$", "reported"), "L1TD1": (r"^L1TD1$", "reported"),
}


def lower_runs(seq, offset=0):
    return [(m.start() + offset, m.end() + offset) for m in re.finditer(r"[a-z]+", seq)]


def polya_windows(seq, start, stop, base):
    """Start positions p in [start, stop) (seq coordinates) whose 15-bp window holds >= 12 copies of `base`."""
    return [p for p in range(max(0, start), min(len(seq) - POLYA_W, stop - 1) + 1)
            if seq[p:p + POLYA_W].count(base) >= POLYA_MIN]


def tsd(w1, w2):
    """Longest exact direct repeat between w1 and w2 whose base composition passes; returns its length or 0."""
    w1, w2 = w1.upper(), w2.upper()
    best = 0
    for i in range(len(w1)):
        for j in range(len(w2)):
            if w1[i] != w2[j] or (i and j and w1[i - 1] == w2[j - 1]):
                continue
            k = 0
            while i + k < len(w1) and j + k < len(w2) and w1[i + k] == w2[j + k]:
                k += 1
            rep = w1[i:i + k]
            if k >= TSD_MIN and max(rep.count(b) for b in "ACGT") <= TSD_MAXFRAC * k:
                best = max(best, k)
    return best


def hallmarks(genome, chrom, clen, ts, te, strand, apply_polya=True):
    """POLYA and TSD for an insertion aligned at genomic [ts, te) on `strand` (3' end = te on +, ts on -)."""
    pad = 120
    lo, hi = max(0, ts - pad), min(clen, te + pad)
    seq = genome.fetch(chrom, lo, hi)
    up = seq.upper()
    L, R, pa = ts, te, False
    if apply_polya:
        if strand == "+":
            ws = polya_windows(up, te - lo - 10, te - lo + 50, "A")
            if ws:
                pa, R = True, lo + max(ws) + POLYA_W
        else:
            # windows ending inside (ts - 50, ts + 10]: starts in [ts - 50 - 15 + 1 ... ] mirrored
            ws = [p for p in polya_windows(up, ts - lo - 50 - POLYA_W + 1, ts - lo + 10 - POLYA_W + 1, "T")]
            if ws:
                pa, L = True, lo + min(ws)
    w1 = up[max(0, L - lo - 40):L - lo + 10]
    w2 = up[R - lo - 10:min(len(up), R - lo + 40)]
    return pa, tsd(w1, w2)


def q2t_map(rec):
    """Aligned (query forward coordinate, target coordinate) M blocks of one PAF record."""
    plus = rec["strand"] == "+"
    t, q = rec["ts"], (rec["qs"] if plus else rec["qe"])
    blocks = []
    for n, op in re.findall(r"(\d+)([MIDNSHP=X])", rec["cg"]):
        n = int(n)
        if op in "M=X":
            if plus:
                blocks.append((q, t, n, 1))
                q += n
            else:
                blocks.append((q - 1, t, n, -1))  # query base q-1 aligns to t, q-2 to t+1, ...
                q -= n
            t += n
        elif op in "DN":
            t += n
        elif op == "I":
            q += n if plus else -n
    return blocks


def t_of(blocks, x):
    """Target coordinate of query base x (nearest aligned base when x is unaligned)."""
    best = None
    for q0, t0, n, d in blocks:
        lo, hi = (q0, q0 + n - 1) if d == 1 else (q0 - n + 1, q0)
        if lo <= x <= hi:
            return t0 + (x - q0) * d
        dist, edge = (lo - x, lo) if x < lo else (x - hi, hi)
        if best is None or dist < best[0]:
            best = (dist, t0 + (edge - q0) * d)
    return best[1]


def junction_calls(rec, junctions):
    """(lost, retained, ambiguous) over parent junctions (query position J, intron length) inside the record."""
    blocks = q2t_map(rec)
    lost = retained = amb = 0
    for J, ilen in junctions:
        if ilen < MIN_INTRON or not (rec["qs"] + JUNC_IN <= J <= rec["qe"] - JUNC_IN):
            continue
        gap = abs(t_of(blocks, J + 9) - t_of(blocks, J - 10)) - 19
        if gap <= LOST_MAX:
            lost += 1
        elif gap >= RETAINED_MIN:
            retained += 1
        else:
            amb += 1
    return lost, retained, amb


def parse_paf(text):
    out = []
    for line in text.splitlines():
        f = line.split("\t")
        out.append({"q": f[0], "qlen": int(f[1]), "qs": int(f[2]), "qe": int(f[3]), "strand": f[4], "t": f[5],
                    "tlen": int(f[6]), "ts": int(f[7]), "te": int(f[8]), "nm": int(f[9]), "bl": int(f[10]),
                    "cg": next((x[5:] for x in f[12:] if x.startswith("cg:Z:")), "")})
    return out


def cached_run(path, cmd):
    if not os.path.exists(path):
        out = subprocess.run(cmd, capture_output=True, text=True, check=True).stdout
        open(path + ".tmp", "w").write(out)
        os.replace(path + ".tmp", path)
    return open(path).read()


def load_rmsk(path):
    iv = collections.defaultdict(list)
    with gzip.open(path, "rt") as fh:
        for i, line in enumerate(fh):
            if i < 3:
                continue
            f = line.split()
            if len(f) < 11:
                continue
            cls = f[10].split("/")[0].rstrip("?")
            if cls in RETRO_CLASSES:
                iv[f[4]].append((int(f[5]) - 1, int(f[6])))
    return {c: do.merge(v) for c, v in iv.items()}


def overlap_bp(ivs, merged):
    ends = [b for _, b in merged]
    return sum(min(b, e) - max(a, s) for a, b in ivs for s, e in do._hits(merged, ends, a, b))


def parent_mrna(genome, g):
    tx = max(g["transcripts"] or [[(g["start0"], g["end"])]], key=lambda t: (len(t), sum(e - s for s, e in t)))
    order = tx if g["strand"] == "+" else tx[::-1]
    seq = "".join(genome.fetch(g["chrom"], s, e).upper() for s, e in tx)  # genomic order; reverse-complemented below
    if g["strand"] == "-":
        seq = seq.translate(str.maketrans("ACGTN", "TGCAN"))[::-1]
    junctions, pos = [], 0
    for k in range(len(order) - 1):
        pos += order[k][1] - order[k][0]
        ilen = (order[k + 1][0] - order[k][1]) if g["strand"] == "+" else (order[k][0] - order[k + 1][1])
        junctions.append((pos, ilen))
    return seq, junctions


def chance_rates(genome, n=1000, seed=13):
    rng = random.Random(seed)
    chroms = [(c, l) for c, l in zip(genome.references, genome.lengths) if re.fullmatch(r"chr(\d+|X)", c)]
    total = sum(l for _, l in chroms)
    pa = ts_ = done = 0
    while done < n:
        x = rng.randrange(total)
        for c, l in chroms:
            if x < l:
                break
            x -= l
        strand = rng.choice("+-")
        if x < 200 or x + 1200 > l:
            continue
        s = genome.fetch(c, x - 150, x + 1150)
        if "N" in s.upper() or not s[150].isupper():
            continue
        p, t = hallmarks(genome, c, l, x, x + 1000, strand)
        pa += p
        ts_ += t >= TSD_MIN
        done += 1
    return pa / n, ts_ / n


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--gff", required=True)
    ap.add_argument("--genome", required=True)
    ap.add_argument("--rmsk", required=True)
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--amy-truth")
    ap.add_argument("--families")
    ap.add_argument("--threads", type=int, default=2)
    a = ap.parse_args()
    os.makedirs(a.outdir, exist_ok=True)
    genome = pysam.FastaFile(a.genome)
    clen = dict(zip(genome.references, genome.lengths))
    gc = f"{a.outdir}/genes.pkl"
    if not os.path.exists(gc):
        pickle.dump(do.load_genes(a.gff), open(gc, "wb"))
    genes = pickle.load(open(gc, "rb"))
    by_name = collections.defaultdict(list)
    for g in genes.values():
        by_name[g["name"]].append(g)

    # ---------------- chance-rate gate ----------------
    cr = f"{a.outdir}/chance.pkl"
    if not os.path.exists(cr):
        pickle.dump(chance_rates(genome), open(cr, "wb"))
    pa_rate, tsd_rate = pickle.load(open(cr, "rb"))
    use_pa, use_tsd = pa_rate <= CHANCE_MAX, tsd_rate <= CHANCE_MAX
    print(f"CHANCE RATES (1,000 random pseudo-insertions): POLYA {pa_rate:.3f} -> {'used' if use_pa else 'DROPPED'}; "
          f"TSD {tsd_rate:.3f} -> {'used' if use_tsd else 'DROPPED'}")

    wanted = set(a.families.split(",")) if a.families else set(FAMILIES) | {"AMY"}
    fams = {f: (k, [g for n, gs in by_name.items() if re.match(p, n) for g in gs if g["chrom"] in clen])
            for f, (p, k) in FAMILIES.items() if f in wanted}
    if a.amy_truth and "AMY" in wanted:
        names = {(r["name"], r["chrom"]) for r in csv.DictReader(open(a.amy_truth), delimiter="\t")}
        fams["AMY"] = ("dev_amy", [g for (n, c) in sorted(names) for g in by_name.get(n, []) if g["chrom"] == c])

    member_rows, fam_rows = [], []
    for fam, (kind, members) in fams.items():
        pcs = [g for g in members if g["biotype"] == "protein_coding"] or members
        parent = min(pcs, key=lambda g: (-g["introns"], g["name"], g["chrom"], g["start0"]))
        others = sorted([g for g in members if g is not parent], key=lambda g: (g["name"], g["chrom"], g["start0"]))
        dropped = max(0, len(others) - CAP)
        others = others[:CAP]

        def region(g):
            L = g["end"] - g["start0"]
            return max(0, g["start0"] - L), min(clen[g["chrom"]], g["end"] + L)

        ps, pe = region(parent)
        pseq = genome.fetch(parent["chrom"], ps, pe)
        tfa, qfa = f"{a.outdir}/{fam}.parent.fa", f"{a.outdir}/{fam}.members.fa"
        open(tfa, "w").write(f">parent\n{pseq.upper()}\n")
        keys = {}
        with open(qfa, "w") as fh:
            for i, g in enumerate(others):
                s, e = region(g)
                keys[f"m{i}"] = (g, s, e)
                fh.write(f">m{i}\n{genome.fetch(g['chrom'], s, e).upper()}\n")
        paf = cached_run(f"{a.outdir}/{fam}.region.paf",
                         ["minimap2", "-c", "-x", "asm20", "-N", "50", "-p", "0.1", "-t", str(a.threads), tfa, qfa])
        chains = do.chains_from_paf(paf.splitlines())
        pt0, pt1 = parent["start0"] - ps, parent["end"] - ps
        t_block = do.merge([(x - ps, y - ps) for x, y in parent["exons"]] + lower_runs(pseq))
        mrna, junctions = parent_mrna(genome, parent)
        mfa = f"{a.outdir}/{fam}.mrna.fa"
        open(mfa, "w").write(f">mrna\n{mrna}\n")

        counts = collections.Counter()
        for key, (g, s, e) in keys.items():
            gq0, gq1 = g["start0"] - s, g["end"] - s
            row = {"family": fam, "kind": kind, "parent": parent["name"], "parent_introns": parent["introns"],
                   "member": g["name"], "chrom": g["chrom"], "start0": g["start0"], "end": g["end"], "biotype": g["biotype"],
                   "member_introns": g["introns"], "region_identity": float("nan"), "shared_nonexonic_unique": 0,
                   "mrna_identity": float("nan"), "mrna_qcov": 0.0, "lost": 0, "retained": 0, "ambiguous": 0,
                   "polya": False, "tsd_len": 0, "class": "UNASSESSED"}
            region_ok = False
            cand = [c for c in chains if c["q"] == key and c["qs"] < gq1 and gq0 < c["qe"] and c["ts"] < pt1 and pt0 < c["te"]]
            if cand:
                c = max(cand, key=lambda c: c["aligned"])
                row["region_identity"] = c["nm"] / c["bl"]
                if c["nm"] / c["bl"] >= MIN_ID:
                    region_ok = True
                    mseq = genome.fetch(g["chrom"], s, e)
                    q_block = do.merge([(x - s, y - s) for x, y in g["exons"]] + lower_runs(mseq))
                    row["shared_nonexonic_unique"] = do.nonexonic_bp(c["recs"], t_block, q_block)
            mtfa = f"{a.outdir}/{fam}.{key}.target.fa"
            if not os.path.exists(f"{a.outdir}/{fam}.{key}.splice.paf"):
                open(mtfa, "w").write(f">{key}\n{genome.fetch(g['chrom'], s, e).upper()}\n")
            sp = parse_paf(cached_run(f"{a.outdir}/{fam}.{key}.splice.paf",
                                      ["minimap2", "-c", "-x", "splice", "-uf", "-N", "50", "-p", "0.1", "-t", "1", mtfa, mfa]))
            sp = [r for r in sp if r["ts"] < gq1 and gq0 < r["te"]]
            mrna_ok = False
            if sp:
                r = max(sp, key=lambda r: (r["qe"] - r["qs"], r["nm"]))
                row["mrna_identity"], row["mrna_qcov"] = r["nm"] / r["bl"], (r["qe"] - r["qs"]) / len(mrna)
                if r["nm"] / r["bl"] >= MIN_ID and r["qe"] - r["qs"] >= 100:
                    mrna_ok = True
                    lost, ret, amb = junction_calls(r, junctions)
                    row.update(lost=lost, retained=ret, ambiguous=amb)
                    reach = (r["qe"] >= len(mrna) - POLYA_REACH)
                    gts, gte = s + r["ts"], s + r["te"]
                    pa, tl = hallmarks(genome, g["chrom"], clen[g["chrom"]], gts, gte, r["strand"], apply_polya=reach)
                    row["polya"], row["tsd_len"] = pa, tl
            if row["shared_nonexonic_unique"] >= NONEXONIC:
                row["class"] = "GENOMIC"
            elif mrna_ok or region_ok:
                pa = row["polya"] and use_pa
                ts_ = row["tsd_len"] >= TSD_MIN and use_tsd
                scorable = row["lost"] + row["retained"] + row["ambiguous"] > 0
                if mrna_ok and row["lost"] >= 1 and row["retained"] == 0 and (pa or ts_ or not (use_pa or use_tsd)):
                    row["class"] = "RETROCOPY"
                elif mrna_ok and not scorable and (use_pa or use_tsd) and (pa or not use_pa) and (ts_ or not use_tsd):
                    row["class"] = "RETROCOPY"
                else:
                    row["class"] = "UNRESOLVED"
            counts[row["class"]] += 1
            member_rows.append(row)
        assessed = sum(v for k, v in counts.items() if k != "UNASSESSED")
        fam_rows.append({"family": fam, "kind": kind, "members": len(members), "dropped": dropped, "parent": parent["name"],
                         "parent_introns": parent["introns"], "assessed": assessed, **{k: counts[k] for k in
                         ("GENOMIC", "RETROCOPY", "UNRESOLVED", "UNASSESSED")},
                         "retro_derived": assessed > 0 and counts["RETROCOPY"] / assessed >= 0.5,
                         "genomic_derived": assessed > 0 and counts["GENOMIC"] / assessed >= 0.5})

        # report-only: retrocopy then genomic duplication
        rc = [(k, v) for k, v in keys.items() if any(r["member"] == v[0]["name"] and r["chrom"] == v[0]["chrom"] and
              r["start0"] == v[0]["start0"] and r["class"] == "RETROCOPY" for r in member_rows)]
        fam_rows[-1]["retro_then_sd"] = 0
        if len(rc) >= 2:
            rfa = f"{a.outdir}/{fam}.retro.fa"
            if not os.path.exists(f"{a.outdir}/{fam}.retro.paf"):
                with open(rfa, "w") as fh:
                    for k, (g, s, e) in rc:
                        fh.write(f">{k}\n{genome.fetch(g['chrom'], s, e).upper()}\n")
            rp = parse_paf(cached_run(f"{a.outdir}/{fam}.retro.paf",
                                      ["minimap2", "-c", "-x", "asm20", "-N", "50", "-p", "0.1", "-t", str(a.threads), rfa, rfa]))
            info = dict(rc)
            hit = set()
            for (qk, tk), recs in itertools.groupby(sorted((r for r in rp if r["q"] != r["t"]), key=lambda r: (r["q"], r["t"])),
                                                      key=lambda r: (r["q"], r["t"])):
                recs = list(recs)
                if sum(r["nm"] for r in recs) / max(1, sum(r["bl"] for r in recs)) < MIN_ID:
                    continue
                (gq, sq, _), (gt, st, _) = info[qk], info[tk]
                qseq, tseq = genome.fetch(gq["chrom"], *info[qk][1:]), genome.fetch(gt["chrom"], *info[tk][1:])
                qb = do.merge([(x - sq, y - sq) for x, y in gq["exons"]] + lower_runs(qseq))
                tb = do.merge([(x - st, y - st) for x, y in gt["exons"]] + lower_runs(tseq))
                if do.nonexonic_bp(recs, tb, qb) >= NONEXONIC:
                    hit.add(qk)
            fam_rows[-1]["retro_then_sd"] = len(hit)

    # ---------------- TE-derived ----------------
    rc_path = f"{a.outdir}/rmsk_retro.pkl"
    if not os.path.exists(rc_path):
        pickle.dump(load_rmsk(a.rmsk), open(rc_path, "wb"))
    rmsk = pickle.load(open(rc_path, "rb"))
    te_rows = []

    def te_row(group, kind, g):
        cds = g["cds"]
        L = sum(e - s for s, e in cds)
        ov = overlap_bp(cds, rmsk.get(g["chrom"], [])) if L else 0
        return {"group": group, "kind": kind, "gene": g["name"], "chrom": g["chrom"], "cds_bp": L, "te_bp": ov,
                "te_frac": ov / L if L else float("nan"), "te_derived": bool(L) and ov / L >= 0.5}
    for grp, (pat, kind) in TE_GROUPS.items():
        for n, gs in sorted(by_name.items()):
            if re.match(pat, n):
                for g in gs:
                    if g["chrom"] in clen and g["biotype"] == "protein_coding":
                        te_rows.append(te_row(grp, kind, g))
    if not a.families:
        for fam, (kind, members) in fams.items():
            for g in members:
                if g["biotype"] == "protein_coding":
                    te_rows.append(te_row(fam, "negative:" + kind, g))

    # ---------------- outputs ----------------
    tag = "" if not a.families else "." + "_".join(sorted(wanted))
    for name, rows in (("members", member_rows), ("te", te_rows)):
        if rows:
            with open(f"{a.outdir}/{name}{tag}.tsv", "w") as fh:
                cols = list(rows[0].keys())
                fh.write("\t".join(cols) + "\n")
                for r in rows:
                    fh.write("\t".join(f"{r[k]:.4f}" if isinstance(r[k], float) else str(r[k]) for k in cols) + "\n")
    print("family\tkind\tmembers\tdropped\tparent\tparent_introns\tassessed\tGENOMIC\tRETROCOPY\tUNRESOLVED\tUNASSESSED\t"
          "retro_derived\tgenomic_derived\tretro_then_sd")
    for r in fam_rows:
        print("\t".join(str(r[k]) for k in ("family", "kind", "members", "dropped", "parent", "parent_introns", "assessed",
                                            "GENOMIC", "RETROCOPY", "UNRESOLVED", "UNASSESSED", "retro_derived",
                                            "genomic_derived", "retro_then_sd")))
    k = lambda kind, field: (sum(1 for r in fam_rows if r["kind"] == kind and r[field]), sum(1 for r in fam_rows if r["kind"] == kind))
    rp_r, sd_r = k("retro_pos", "retro_derived"), k("sd_neg", "retro_derived")
    rp_g, sd_g = k("retro_pos", "genomic_derived"), k("sd_neg", "genomic_derived")
    if not a.families:
        print(f"\nR1 RETROCOPY: retro positives RETRO-DERIVED {rp_r[0]}/{rp_r[1]}; SD negatives {sd_r[0]}/{sd_r[1]} -> "
              f"{'SUPPORTED' if rp_r[0] >= 12 and sd_r[0] <= 1 else 'NOT SUPPORTED'}")
        print(f"R2 GENOMIC (repeats excluded): SD negatives GENOMIC-DERIVED {sd_g[0]}/{sd_g[1]}; retro positives "
              f"{rp_g[0]}/{rp_g[1]} -> {'SUPPORTED' if sd_g[0] >= 10 and rp_g[0] <= 1 else 'NOT SUPPORTED'}")
        for dev in ("dev_core", "dev_retro", "dev_retro_ho", "dev_amy"):
            print(f"DEVELOPMENT {dev}: RETRO-DERIVED {k(dev, 'retro_derived')}, GENOMIC-DERIVED {k(dev, 'genomic_derived')}")
        erv = [r for r in te_rows if r["group"] == "ERV-env"]
        neg = [r for r in te_rows if r["kind"].startswith("negative:") and r["cds_bp"]]
        ne, nn = sum(r["te_derived"] for r in erv), sum(r["te_derived"] for r in neg)
        print(f"R3 TE-DERIVED: ERV-env {ne}/{len(erv)}; negatives {nn}/{len(neg)} ({nn / max(1, len(neg)):.3f}) -> "
              f"{'SUPPORTED' if ne >= 6 and nn <= 0.05 * len(neg) else 'NOT SUPPORTED'}")
        for grp in TE_GROUPS:
            rows = [r for r in te_rows if r["group"] == grp]
            print(f"  {grp}: " + " ".join(f"{r['gene']}:{r['te_frac']:.2f}" for r in rows))
        print("  negatives TE-DERIVED: " + " ".join(f"{r['gene']}({r['group']}):{r['te_frac']:.2f}" for r in neg if r["te_derived"]))


if __name__ == "__main__":
    main()
