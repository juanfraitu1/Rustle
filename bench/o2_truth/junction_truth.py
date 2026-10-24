#!/usr/bin/env python3
"""Arm 2 of docs/archive/2026-10/merged/PREREG_o2_anchor_2026-10-13.md (Amendment 2): junction truth from rival copies that cannot splice.

    junction_truth.py RUN_PREFIX CATALOG_COPIES.tsv BAM GENOME_FA OUTDIR

For every contested molecule with a read-star row whose genome primary alignment has >= 1 intron (N >= 30 bp): the read
(original orientation) is aligned to each read-star candidate's genomic span +- 1 kb (minimap2 -c -x splice:hq -uf); a
candidate BREAKS a read junction q if its best hit has an intron within 5 bp of q with a non-canonical motif (transcript
orientation; canonical = GT-AG, GC-AG, AT-AC), MAKES it if canonical. EXCLUDED = breaks >= 1 junction; PRODUCES = makes
all. Scores O2 (J1 falsification, J2 determined truth) and the aligner's primary. Molecules named `*|trim` are skipped.
"""
import collections, csv, json, os, re, subprocess, sys

MIN_INTRON = 30
TOL = 5
PAD = 1000
CANON = {("GT", "AG"), ("GC", "AG"), ("AT", "AC")}
COMP = str.maketrans("ACGTNacgtn", "TGCANtgcan")
CIG = re.compile(r"(\d+)([MIDNSHP=X])")


def revcomp(s):
    return s.translate(COMP)[::-1]


def read_junctions(cigar, flag, qlen):
    """Read positions (original read orientation) of the introns (N >= MIN_INTRON) of a BAM record. A junction at q lies
    between read bases q-1 and q. `qlen` = the full read length (SEQ + hard clips)."""
    q, out = 0, []
    for n, op in CIG.findall(cigar):
        n = int(n)
        if op in "MIS=XH":
            q += n
        elif op == "N" and n >= MIN_INTRON:
            out.append(q)
    if int(flag) & 16:
        out = sorted(qlen - x for x in out)
    return out


def paf_introns(cg, strand, qs, qe, ts):
    """[(read junction position, target intron start, end)] from a PAF `cg` cigar (forward query coordinates)."""
    q = qs if strand == "+" else qe
    t, out = ts, []
    for n, op in CIG.findall(cg):
        n = int(n)
        if op in "M=X":
            q += n if strand == "+" else -n
            t += n
        elif op == "I":
            q += n if strand == "+" else -n
        elif op == "D":
            t += n
        elif op == "N":
            out.append((q, t, t + n))
            t += n
    return out


def motif(target_seq, t0, t1, strand):
    """Intron dinucleotides (donor, acceptor) in transcript orientation."""
    if strand == "+":
        return target_seq[t0:t0 + 2].upper(), target_seq[t1 - 2:t1].upper()
    intr = revcomp(target_seq[t0:t1].upper())
    return intr[:2], intr[-2:]


def canonical(m):
    return m in CANON


def classify(juncs, introns):
    """introns = [(read pos, t0, t1, canonical)] of a candidate's best hit -> EXCLUDED / PRODUCES / UNDETERMINED."""
    made = 0
    for q in juncs:
        near = [x for x in introns if abs(x[0] - q) <= TOL]
        if any(not x[3] for x in near):
            return "EXCLUDED"
        made += bool(near)
    return "PRODUCES" if made == len(juncs) else "UNDETERMINED"


def overlaps(a, b):
    return a[0] == b[0] and a[1] < b[2] and b[1] < a[2]


def one_locus(spans):
    return bool(spans) and all(overlaps(a, b) for i, a in enumerate(spans) for b in spans[i + 1:])


def tsv(path):
    csv.field_size_limit(sys.maxsize)
    with open(path) as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def ref_end(pos0, cigar):
    return pos0 + sum(int(n) for n, op in CIG.findall(cigar) if op in "MDN=X")


def fasta_fetch(genome, chrom, s, e):
    out = subprocess.run(["samtools", "faidx", genome, f"{chrom}:{s + 1}-{e}"], capture_output=True, text=True,
                         check=True).stdout.splitlines()
    return "".join(out[1:])


def main(prefix, catalog_tsv, bam, genome, out):
    os.makedirs(out, exist_ok=True)
    span = {(r["family_id"], r["copy_idx"]): (r["chrom"], int(r["start"]), int(r["end"])) for r in tsv(catalog_tsv)}
    join = {(r["family_id"], r["copy_index"]): r["catalog_copy_idx"] for r in tsv(prefix + ".family_join.tsv")}
    asg = {(r["read_name"], r["family_id"]): r for r in tsv(prefix + ".assignments.tsv")}
    mols = {}
    for s in tsv(prefix + ".star_reads.tsv"):
        k = (s["read_name"], s["family_id"])
        a = asg.get(k)
        if a is None or a["contested"] != "1" or s["read_name"].endswith("|trim") or not s["candidates"]:
            continue
        cands = [join[(s["family_id"], c)] for c in s["candidates"].split(",")]
        mols[k] = {"cands": cands, "status": a["status"], "copy": a["catalog_copy_idx"], "cols": s["columns"]}
    c = collections.Counter(population_rows=len(mols))
    names = sorted({k[0] for k in mols})
    fams = sorted({k[1] for k in mols})
    # primary records of the molecules, from their families' copy windows (+- 50 kb)
    with open(f"{out}/names.txt", "w") as fh:
        fh.write("\n".join(names) + "\n")
    win = sorted({(span[(f, i)][0], max(0, span[(f, i)][1] - 50000), span[(f, i)][2] + 50000)
                  for (f, i) in span if f in set(fams)})
    with open(f"{out}/windows.bed", "w") as fh:
        fh.writelines(f"{ch}\t{s}\t{e}\n" for ch, s, e in win)
    prim = {}
    res = subprocess.run(["samtools", "view", "-F", "2308", "-L", f"{out}/windows.bed", "-N", f"{out}/names.txt", bam],
                         capture_output=True, text=True, check=True).stdout
    for line in res.splitlines():
        f = line.split("\t")
        if f[0] in prim or f[9] == "*":
            continue
        hard = sum(int(n) for n, op in CIG.findall(f[5]) if op == "H")
        seq = f[9] if not int(f[1]) & 16 else revcomp(f[9])
        prim[f[0]] = {"juncs": read_junctions(f[5], f[1], len(f[9]) + hard), "seq": seq,
                      "span": (f[2], int(f[3]) - 1, ref_end(int(f[3]) - 1, f[5]))}
    # per family: align the spliced molecules to their candidates' genomic spans
    per = {}
    for fam in fams:
        ks = [k for k in mols if k[1] == fam and k[0] in prim and prim[k[0]]["juncs"]]
        c["no_primary"] += sum(1 for k in mols if k[1] == fam and k[0] not in prim)
        c["unspliced"] += sum(1 for k in mols if k[1] == fam and k[0] in prim and not prim[k[0]]["juncs"])
        if not ks:
            continue
        tids = sorted({x for k in ks for x in mols[k]["cands"]}, key=int)
        tseq = {}
        with open(f"{out}/_t.fa", "w") as fh:
            for t in tids:
                ch, s, e = span[(fam, t)]
                s0 = max(0, s - PAD)
                tseq[t] = fasta_fetch(genome, ch, s0, e + PAD)
                fh.write(f">{t}\n{tseq[t]}\n")
        with open(f"{out}/_r.fa", "w") as fh:
            fh.writelines(f">{k[0]}\n{prim[k[0]]['seq']}\n" for k in ks)
        paf = subprocess.run(["minimap2", "-c", "-x", "splice:hq", "-uf", "--secondary=yes", "-N", "50", "-p", "0.1",
                              "-t", "2", f"{out}/_t.fa", f"{out}/_r.fa"], capture_output=True, text=True, check=True).stdout
        best = {}
        for l in paf.splitlines():
            f = l.split("\t")
            tags = dict(t.split(":", 2)[0::2] for t in f[12:] if t.count(":") >= 2)
            aswt = int(tags.get("AS", "0"))
            key = (f[0], f[5])
            if key not in best or aswt > best[key][0]:
                best[key] = (aswt, f[4], int(f[2]), int(f[3]), int(f[7]), tags.get("cg", ""))
        for k in ks:
            cls, detail = {}, {}
            for t in set(mols[k]["cands"]):
                h = best.get((k[0], t))
                if h is None:
                    cls[t] = "UNDETERMINED"
                    continue
                _, strand, qs, qe, ts, cg = h
                intr = [(q, t0, t1, canonical(motif(tseq[t], t0, t1, strand)))
                        for q, t0, t1 in paf_introns(cg, strand, qs, qe, ts)]
                cls[t] = classify(prim[k[0]]["juncs"], intr)
                detail[t] = [(q, "-".join(motif(tseq[t], t0, t1, strand))) for q, t0, t1, _ in
                             [(q, t0, t1, None) for q, t0, t1 in paf_introns(cg, strand, qs, qe, ts)]]
            per[k] = (cls, detail)
    for f in ("_t.fa", "_r.fa"):
        if os.path.exists(f"{out}/{f}"):
            os.remove(f"{out}/{f}")
    # scoring
    ex = {"wrong_J1": [], "wrong_J2": [], "right_J2": []}
    fam_c = collections.defaultdict(collections.Counter)
    for k, (cls, detail) in per.items():
        m = mols[k]
        c["spliced_scored"] += 1
        excl = {t for t, v in cls.items() if v == "EXCLUDED"}
        keep = [t for t in cls if t not in excl]
        p = prim[k[0]]["span"]
        al = max(((t, max(0, min(p[2], span[(k[1], t)][2]) - max(p[1], span[(k[1], t)][1])))
                  for t in cls if span[(k[1], t)][0] == p[0]), key=lambda x: x[1], default=(None, 0))
        aligner_copy = al[0] if al[1] > 0 else None
        if not excl or not keep:
            continue
        c["informative"] += 1
        fam_c[k[1]]["informative"] += 1
        c["aligner_on_excluded"] += aligner_copy in excl
        assigned = m["status"] == "assigned"
        if assigned:
            c["J1_assigned"] += 1
            fam_c[k[1]]["J1_assigned"] += 1
            if m["copy"] in excl:
                c["J1_wrong"] += 1
                ex["wrong_J1"].append({"read": k[0], "family": k[1], "o2_copy": m["copy"], "classes": cls,
                                       "motifs": detail, "columns": m["cols"][:300]})
        spans = [span[(k[1], t)] for t in keep]
        if one_locus(spans) and all(cls[t] == "PRODUCES" for t in keep):
            c["determined"] += 1
            if assigned:
                c["J2_assigned"] += 1
                ok = any(overlaps(span[(k[1], m["copy"])], s) for s in spans)
                c["J2_right" if ok else "J2_wrong"] += 1
                rec = {"read": k[0], "family": k[1], "o2_copy": m["copy"], "truth": keep, "classes": cls,
                       "motifs": detail, "aligner_copy": aligner_copy, "columns": m["cols"][:300]}
                if ok and len(ex["right_J2"]) < 20:
                    ex["right_J2"].append(rec)
                if not ok:
                    ex["wrong_J2"].append(rec)
            c["determined_aligner_right"] += aligner_copy is not None and any(
                overlaps(span[(k[1], aligner_copy)], s) for s in spans)

    def verdict(w, n):
        return "UNDERPOWERED" if n < 30 else ("PASS" if w / n <= 0.01 else "FAIL")
    resj = {"counts": dict(c),
            "J1_falsification": {"assigned": c["J1_assigned"], "wrong": c["J1_wrong"],
                                 "verdict": verdict(c["J1_wrong"], c["J1_assigned"])},
            "J2_determined": {"assigned": c["J2_assigned"], "right": c["J2_right"], "wrong": c["J2_wrong"],
                              "recall": round(c["J2_right"] / c["determined"], 4) if c["determined"] else None,
                              "verdict": verdict(c["J2_wrong"], c["J2_assigned"])},
            "aligner": {"on_excluded_of_informative": c["aligner_on_excluded"],
                        "right_of_determined": c["determined_aligner_right"]},
            "per_family": {f: dict(v) for f, v in sorted(fam_c.items(), key=lambda x: -x[1]["informative"])},
            "examples": ex}
    json.dump(resj, open(f"{out}/junction_truth.json", "w"), indent=1)
    print(json.dumps({k: v for k, v in resj.items() if k not in ("per_family", "examples")}, indent=1))


if __name__ == "__main__":
    main(*sys.argv[1:6])
