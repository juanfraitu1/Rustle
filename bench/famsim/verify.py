"""verify: re-derive the simulated condition from the PRODUCTS (FASTA, truth GTF, copies.fa, reads) and write
DIR/verify.tsv — `claim  copy  expected  observed  status  note`, status PASS | FAIL | INFO.

The generator's bookkeeping (manifest ops) says WHAT was intended; every claim here is measured on the files an
outside reader has: sequences read back from genome.truth.fa, motifs at the truth junctions, minimap2 identities between
the planted sequences, k-mer presence/orientation of exons, read counts and read junctions. INFO rows carry measured
values that have no pass bar (e.g. identity between structurally different copies).
"""
import collections
import json
import os
import subprocess
import tempfile

from .chromosome import load_planted
from .model import rc

K = 11


class Report:
    def __init__(self):
        self.rows = []

    def add(self, claim, copy, expected, observed, ok, note=""):
        status = "INFO" if ok is None else ("PASS" if ok else "FAIL")
        self.rows.append((claim, copy, str(expected), str(observed), status, note))

    def write(self, path):
        with open(path, "w") as fh:
            fh.write("claim\tcopy\texpected\tobserved\tstatus\tnote\n")
            for r in self.rows:
                fh.write("\t".join(r) + "\n")

    @property
    def failed(self):
        return [r for r in self.rows if r[4] == "FAIL"]


# ---------------------------------------------------------------- helpers
def kmers(s, k=K):
    return {s[i:i + k] for i in range(len(s) - k + 1)}


def kmer_fraction(query, target_kmers, k=K):
    q = kmers(query, k)
    return (sum(1 for x in q if x in target_kmers) / len(q)) if q else 0.0


def read_fasta(path):
    out, name = {}, None
    for line in open(path):
        if line.startswith(">"):
            name = line[1:].split()[0]; out[name] = []
        else:
            out[name].append(line.strip())
    return {k: "".join(v) for k, v in out.items()}


def minimap2_paf(query_fa, target_fa, preset="-cx asm20 --eqx", extra="", threads=1):
    cmd = f"minimap2 {preset} {extra} -t {threads} {target_fa} {query_fa} 2>/dev/null"
    out = subprocess.run(cmd, shell=True, stdout=subprocess.PIPE, text=True).stdout
    rows = []
    for line in out.splitlines():
        f = line.split("\t")
        if len(f) < 12:
            continue
        de = next((float(x[5:]) for x in f[12:] if x.startswith("de:f:")), None)
        rows.append({"q": f[0], "qlen": int(f[1]), "qs": int(f[2]), "qe": int(f[3]), "strand": f[4], "t": f[5], "tlen": int(f[6]),
                     "ts": int(f[7]), "te": int(f[8]), "nmatch": int(f[9]), "alen": int(f[10]), "mapq": int(f[11]), "de": de})
    return rows


def nw_identity(a, b):
    """Global alignment identity (matches / columns), linear gaps; the fallback when minimap2 finds no hit."""
    n, m = len(a), len(b)
    if n == 0 or m == 0:
        return 0.0
    prev = list(range(0, -m - 1, -1)); tb = [[0] * (m + 1) for _ in range(n + 1)]
    for j in range(1, m + 1):
        tb[0][j] = 2
    for i in range(1, n + 1):
        cur = [-i] + [0] * m; tb[i][0] = 1; ai = a[i - 1]
        for j in range(1, m + 1):
            d = prev[j - 1] + (1 if ai == b[j - 1] else -1); u = prev[j] - 1; l = cur[j - 1] - 1
            if d >= u and d >= l:
                cur[j] = d; tb[i][j] = 0
            elif u >= l:
                cur[j] = u; tb[i][j] = 1
            else:
                cur[j] = l; tb[i][j] = 2
        prev = cur
    i, j, match, cols = n, m, 0, 0
    while i > 0 or j > 0:
        t = tb[i][j]
        if t == 0:
            match += a[i - 1] == b[j - 1]; i -= 1; j -= 1
        elif t == 1:
            i -= 1
        else:
            j -= 1
        cols += 1
    return match / cols


def gtf_transcripts(path):
    """{tid: {"gene", "contig", "strand", "exons": [(s0, e1, label)] in file order, "n": int}}"""
    tx = {}
    for line in open(path):
        f = line.rstrip("\n").split("\t")
        if len(f) < 9 or f[2] != "exon":
            continue
        import re
        tid = re.search(r'transcript_id "([^"]+)"', f[8]).group(1)
        gid = re.search(r'gene_id "([^"]+)"', f[8]).group(1)
        lab = re.search(r'exon_label "([^"]+)"', f[8])
        t = tx.setdefault(tid, {"gene": gid, "contig": f[0], "strand": f[6], "exons": []})
        t["exons"].append((int(f[3]) - 1, int(f[4]), lab.group(1) if lab else None))
    return tx


# ---------------------------------------------------------------- the claims
def run(out_dir, threads=1, allow_background_homology=False, log=print):
    import pysam
    planted, man = load_planted(out_dir)
    byid = {p.id: p for p in planted}
    template_seq = json.load(open(os.path.join(out_dir, "template.json")))
    template_exons = {e["label"]: template_seq["seq"][e["start"]:e["end"]] for e in template_seq["exons"]}
    template_motifs = set(man.get("template_motifs") or [])
    n_template_exons = len(template_seq["exons"])
    R = Report()
    truth_fa = pysam.FastaFile(os.path.join(out_dir, "genome.truth.fa"))
    copies = read_fasta(os.path.join(out_dir, "copies.fa"))
    tx = gtf_transcripts(os.path.join(out_dir, "truth.gtf"))
    spec_reads = man["spec"].get("reads", {})

    # 1. planted sequences are in the genome, exon by exon, at the truth coordinates
    for p in planted:
        g = truth_fa.fetch(p.contig, p.pos, p.end).upper()
        if p.strand == "-":
            g = rc(g)
        same_span = g == copies[p.id]
        bad = 0
        for tid, t in tx.items():
            if t["gene"] != p.id:
                continue
            for s, e, lab in t["exons"]:
                es = truth_fa.fetch(p.contig, s, e).upper()
                if p.strand == "-":
                    es = rc(es)
                ex = next((x for x in p.model.exons if x.label == lab), None)
                if ex is None or p.model.seq[ex.start:ex.end] != es:
                    bad += 1
        R.add("planted_in_genome", p.id, "span and every truth exon read back identical", f"span={'ok' if same_span else 'DIFF'} exon_mismatches={bad}",
              same_span and bad == 0)

    # 2. every truth junction canonical (or the template's own non-canonical motif)
    for tid, t in tx.items():
        p = byid[t["gene"]]
        ex = sorted(t["exons"])
        bad, motifs = [], []
        for (s1, e1, _), (s2, e2, _) in zip(ex, ex[1:]):
            d = truth_fa.fetch(t["contig"], e1, e1 + 2).upper(); a = truth_fa.fetch(t["contig"], s2 - 2, s2).upper()
            motif = (d + a) if t["strand"] == "+" else (rc(a) + rc(d))
            motifs.append(motif)
            if motif != "GTAG" and motif not in template_motifs:
                bad.append(f"{e1}-{s2}:{motif}")
        R.add("junctions_canonical", tid, "GTAG at every intron", f"{len(motifs) - len(bad)}/{len(motifs)} GTAG" + (f" bad={bad}" if bad else ""),
              not bad, "non-canonical motifs inherited from the template are allowed" if template_motifs - {"GTAG"} else "")

    # 3. exon counts (RNA, iso0) follow the ops; length bookkeeping
    for p in planted:
        if p.kind != "copy":
            continue
        n_rna, n_dna, dlen = n_template_exons, n_template_exons, 0
        for o in p.model.ops:
            op = o["op"]
            if op == "exon_delete":
                n_rna -= 1; n_dna -= 1; dlen -= o["removed_bp"]
            elif op == "exon_insert":
                n_rna += 1; n_dna += 1; dlen += o["length"] + 4
            elif op == "splice_kill":
                n_rna -= 1
            elif op == "invert" and not o.get("whole"):
                if not o.get("keep_in_rna"):
                    n_rna -= len(o.get("exons_inverted", []))
            elif op == "truncate":
                n_rna -= (n_dna - o["exons_left"]); n_dna = o["exons_left"]; dlen -= o["removed_bp"]
            elif op == "indel":
                dlen += sum(len(ev.get("ins", "")) - ev.get("del", 0) for ev in o["events"])
            elif op == "intron_resize":
                dlen += o["new"] - o["old"]
            elif op == "convert" and "bp" in o and "exon" in o:
                src = byid.get(o["source"])
                mine = next((x for x in p.model.exons if x.label == o["label"]), None)
                # exon conversion can change the length when the donor's exon differs in length
                dlen += 0 if mine is None else 0
        got_rna = len(tx[f"{p.id}.iso0"]["exons"]) if f"{p.id}.iso0" in tx else 0
        R.add("rna_exon_count", p.id, n_rna, got_rna, got_rna == n_rna, "iso0 exons in truth.gtf vs template + op deltas")
        exp_len = len(template_seq["seq"]) + dlen
        has_conv = any(o["op"] == "convert" for o in p.model.ops)
        R.add("length_bookkeeping", p.id, exp_len, len(copies[p.id]), None if has_conv else len(copies[p.id]) == exp_len,
              "template length + op deltas" + (" (INFO: conversion may change the length)" if has_conv else ""))

    # 4. pairwise identity between planted sequences (genomic) and between RNA chains
    cfa = os.path.join(out_dir, "copies.fa")
    ava = minimap2_paf(cfa, cfa, extra="-X -N 50 -p 0.1 --secondary=yes", threads=threads)
    best = {}
    for r in ava:
        if r["q"] == r["t"]:
            continue
        key = (r["q"], r["t"])
        if key not in best or r["nmatch"] > best[key]["nmatch"]:
            best[key] = r
    with tempfile.TemporaryDirectory() as td:
        chain_fa = os.path.join(td, "chains.fa")
        with open(chain_fa, "w") as fh:
            for p in planted:
                fh.write(f">{p.id}\n{p.model.chain_seq()}\n")
        cava = minimap2_paf(chain_fa, chain_fa, extra="-k15 -w5 -X -N 50 -p 0.1 --secondary=yes", threads=threads)
    cbest = {}
    for r in cava:
        if r["q"] != r["t"]:
            key = (r["q"], r["t"])
            if key not in cbest or r["nmatch"] > cbest[key]["nmatch"]:
                cbest[key] = r
    A = [p for p in planted if p.kind == "copy"]
    pair_ident = {}
    for i, p in enumerate(A):
        for q in A[i + 1:]:
            r = best.get((p.id, q.id)) or best.get((q.id, p.id))
            cr = cbest.get((p.id, q.id)) or cbest.get((q.id, p.id))
            gid = r["nmatch"] / r["alen"] if r else None
            gcov = r["alen"] / min(len(copies[p.id]), len(copies[q.id])) if r else 0.0
            if cr:
                cid_ = cr["nmatch"] / cr["alen"]; src = "minimap2"
            else:
                a, b = p.model.chain_seq(), q.model.chain_seq()
                cid_ = nw_identity(a, b) if len(a) * len(b) <= 25_000_000 else None; src = "nw" if cid_ is not None else "none"
            pair_ident[(p.id, q.id)] = (gid, cid_)
            ops_p = [o["op"] for o in p.model.ops]; ops_q = [o["op"] for o in q.model.ops]
            snp_only = all(o in ("snp",) for o in ops_p + ops_q)
            rate = sum(o["rate"] for o in p.model.ops + q.model.ops if o["op"] == "snp")
            if not ops_p and not ops_q:
                ok = gid is not None and gid >= 0.9999 and cid_ is not None and cid_ >= 0.9999
                exp = "1.0000 (identical copies)"
            elif snp_only:
                tol = 0.35 * rate + 0.003
                ok = gid is not None and abs((1 - gid) - rate) <= tol
                exp = f"~{1 - rate:.4f} +- {tol:.4f}"
            else:
                ok = None; exp = "measured"
            R.add("pair_identity_genomic", f"{p.id}~{q.id}", exp,
                  "no hit" if gid is None else f"{gid:.4f} (gap-compressed {1 - r['de']:.4f}, cov {gcov:.2f})", ok,
                  "minimap2 -cx asm20 --eqx: nmatch/alen of the best hit (gaps count as columns) and 1 - de")
            R.add("pair_identity_chain", f"{p.id}~{q.id}", exp if ok is not None else "measured",
                  "none" if cid_ is None else f"{cid_:.4f} ({src})", None if ok is None else (cid_ is not None and (abs((1 - cid_) - rate) <= 0.35 * rate + 0.006 if snp_only else cid_ >= 0.9999)),
                  "spliced chains")
    # decoys unrelated to A
    for d in [p for p in planted if p.kind == "decoy"]:
        hits = [r for r in ava if (r["q"] == d.id and byid[r["t"]].kind == "copy") or (r["t"] == d.id and byid[r["q"]].kind == "copy")]
        strong = [r for r in hits if r["alen"] >= 300 and r["nmatch"] / r["alen"] >= 0.8]
        R.add("decoy_unrelated", d.id, "no hit >= 300 bp at >= 0.80 identity to any copy of A",
              "none" if not strong else "; ".join(f"{r['q']}~{r['t']} {r['nmatch'] / r['alen']:.3f}x{r['alen']}" for r in strong[:3]), not strong)

    # 5. exon-level structural claims, from k-mers of the template exon vs the copy's DNA, and labels in the RNA chain
    for p in A:
        ck = kmers(copies[p.id]); ckr = kmers(rc(copies[p.id]))
        iso0 = tx.get(f"{p.id}.iso0", {"exons": []})
        rna_labels = [lab for _, _, lab in (iso0["exons"] if p.strand == "+" else iso0["exons"])]
        for o in p.model.ops:
            op = o["op"]
            if op == "exon_delete":
                f = kmer_fraction(template_exons[o["label"]], ck | ckr)
                R.add("deleted_exon_absent", f"{p.id}:{o['label']}", "< 0.20 of its 11-mers in the copy (a present exon at <= 10% divergence keeps >= 0.3)",
                      f"{f:.3f}", f < 0.20, "short exons share a few 11-mers with the rest of the gene by chance")
                R.add("deleted_exon_not_in_rna", f"{p.id}:{o['label']}", "label absent from iso0", "absent" if o["label"] not in rna_labels else "PRESENT",
                      o["label"] not in rna_labels)
            elif op == "splice_kill":
                dn = copies[p.id][o["donor_at"]:o["donor_at"] + 2]
                R.add("splice_killed", f"{p.id}:{o['label']}", "donor CT in DNA, exon absent from iso0",
                      f"donor={dn} rna={'absent' if o['label'] not in rna_labels else 'PRESENT'}", dn == "CT" and o["label"] not in rna_labels)
                ex = next(x for x in p.model.exons if x.label == o["label"])
                f = kmer_fraction(template_exons[o["label"]], ck) if o["label"] in template_exons else None
                R.add("killed_exon_in_dna", f"{p.id}:{o['label']}", ">= 0.05 of the template exon's 11-mers still in the copy", f"{f:.3f}" if f is not None else "n/a",
                      f is None or f >= 0.05)
            elif op == "exon_insert":
                ex = next(x for x in p.model.exons if x.label == o["label"])
                present = o["label"] in rna_labels
                flank = copies[p.id][ex.start - 2:ex.start] + ".." + copies[p.id][ex.end:ex.end + 2]
                R.add("inserted_exon", f"{p.id}:{o['label']}", "in iso0, AG..GT flanks", f"rna={'present' if present else 'ABSENT'} flanks={flank}",
                      present and flank == "AG..GT")
            elif op == "exon_shuffle":
                la, lb = o["labels"]
                ia, ib = (rna_labels.index(la) if la in rna_labels else -1), (rna_labels.index(lb) if lb in rna_labels else -1)
                R.add("exons_shuffled", f"{p.id}:{la}<->{lb}", f"{lb} before {la} in iso0", f"order={rna_labels}", 0 <= ib < ia)
            elif op == "invert" and not o.get("whole"):
                seg = copies[p.id][o["start"]:o["end"]]
                tk = kmers(template_seq["seq"])
                fwd, rcf = kmer_fraction(seg, tk), kmer_fraction(rc(seg), tk)
                R.add("segment_inverted", f"{p.id}:{o['kind']}", "rc 11-mers >= 0.05 and >= 5x forward", f"fwd={fwd:.3f} rc={rcf:.3f}",
                      rcf >= 0.05 and rcf >= 5 * fwd)
                for lab in o.get("exons_inverted", []):
                    R.add("inverted_exon_not_in_rna", f"{p.id}:{lab}", "absent from iso0" if not o.get("keep_in_rna") else "present",
                          "absent" if lab not in rna_labels else "present", (lab not in rna_labels) != bool(o.get("keep_in_rna")))
            elif op == "truncate":
                R.add("truncated", p.id, f"{o['exons_left']} DNA exons, {o['removed_bp']} bp removed from {o['side']}'",
                      f"{len(p.model.exons)} exons, len {len(copies[p.id])}", len(p.model.exons) == o["exons_left"])
            elif op == "convert":
                src = byid.get(o["source"])
                if src is not None and "label" in o:
                    me = next(x for x in p.model.exons if x.label == o["label"]); se = next(x for x in src.model.exons if x.label == o["label"])
                    same = copies[p.id][me.start:me.end] == copies[src.id][se.start:se.end]
                    R.add("exon_converted", f"{p.id}:{o['label']}<-{src.id}", "identical to the donor's exon", "identical" if same else "DIFFERENT", same)
        if p.strand == "-" or any(o.get("whole") for o in p.model.ops if o["op"] == "invert"):
            g = truth_fa.fetch(p.contig, p.pos, p.end).upper()
            R.add("whole_inverted", p.id, "genome holds the reverse complement of the model", "yes" if rc(g) == copies[p.id] else "NO", rc(g) == copies[p.id])

    # 6. reference-absent copies and background homology (copies vs the two genomes)
    ref_fa, tru_fa = os.path.join(out_dir, "genome.ref.fa"), os.path.join(out_dir, "genome.truth.fa")
    vs_ref = minimap2_paf(cfa, ref_fa, extra="-N 50 -p 0.1 --secondary=yes", threads=threads) if os.path.getsize(ref_fa) > 0 else []
    vs_tru = minimap2_paf(cfa, tru_fa, extra="-N 50 -p 0.1 --secondary=yes", threads=threads)
    planted_iv = [(p.contig, p.pos, p.end, p.id) for p in planted]

    def best_of(rows, q):
        rs = [r for r in rows if r["q"] == q]
        return max(rs, key=lambda r: r["nmatch"]) if rs else None

    for p in planted:
        if p.in_reference:
            continue
        br, bt = best_of(vs_ref, p.id), best_of(vs_tru, p.id)
        ir = br["nmatch"] / br["alen"] if br else 0.0; it = bt["nmatch"] / bt["alen"] if bt else 0.0
        nr, nt = (br["nmatch"] if br else 0), (bt["nmatch"] if bt else 0)
        R.add("absent_from_reference", p.id, "best hit in genome.ref.fa matches fewer bases than in genome.truth.fa (= itself, every base)",
              f"ref {nr} matched bases, identity {ir:.4f}" + (f" at {br['t']}:{br['ts']}-{br['te']}" if br else "") + f" | truth {nt}, identity {it:.4f}",
              nr < nt, "an identical absent copy is undetectable by sequence; its FAIL here is the honest verdict" if nr >= nt else "")
    for p in planted:
        stray = []
        for r in vs_tru:
            if r["q"] != p.id or r["alen"] < 300 or r["nmatch"] / r["alen"] < 0.8:
                continue
            inside = any(r["t"] == c and r["ts"] < e and r["te"] > s for c, s, e, _ in planted_iv)
            if not inside:
                stray.append(f"{r['t']}:{r['ts']}-{r['te']} {r['nmatch'] / r['alen']:.3f}x{r['alen']}")
        R.add("background_clean", p.id, "no hit >= 300 bp at >= 0.80 identity outside the planted intervals",
              "none" if not stray else "; ".join(stray[:3]), None if allow_background_homology else not stray,
              "unplanned homology in the background would change the condition" if stray else "")

    # 7. reads: counts and junctions
    counts = collections.Counter()
    for line in open(os.path.join(out_dir, "reads.fq")):
        if line.startswith("@") and "|" in line:
            counts[line[1:].split("|")[0]] += 1
    per_copy = int(spec_reads.get("per_copy", 50))
    for p in planted:
        exp = per_copy if p.expression is None else int(p.expression)
        got = counts.get(p.id, 0)
        R.add("read_count", p.id, exp, got, abs(got - exp) <= max(1, len(p.isoform_list()) - 1), "isoform rounding allows +-1 per extra isoform")
    truth_j = {}
    for tid, t in tx.items():
        ex = sorted(t["exons"])
        truth_j[tid] = {(e1, s2) for (s1, e1, _), (s2, e2, _) in zip(ex, ex[1:])}
    bad = collections.Counter(); tot = collections.Counter()
    with open(os.path.join(out_dir, "reads.truth.tsv")) as fh:
        hdr = fh.readline().rstrip("\n").split("\t")
        for line in fh:
            f = dict(zip(hdr, line.rstrip("\n").split("\t")))
            js = {tuple(int(x) for x in j.split("-")) for j in f["junctions"].split(",") if j}
            tot[f["copy"]] += 1
            if not js <= truth_j.get(f"{f['copy']}.{f['isoform']}", set()):
                bad[f["copy"]] += 1
    for p in planted:
        if tot[p.id]:
            R.add("read_junctions_in_truth", p.id, "every read's junctions within its isoform's", f"{tot[p.id] - bad[p.id]}/{tot[p.id]}", bad[p.id] == 0)

    R.write(os.path.join(out_dir, "verify.tsv"))
    n_fail = len(R.failed)
    log(f"verify: {len(R.rows)} claims, {sum(1 for r in R.rows if r[4] == 'PASS')} PASS, {n_fail} FAIL, "
        f"{sum(1 for r in R.rows if r[4] == 'INFO')} INFO -> verify.tsv")
    for r in R.failed:
        log(f"  FAIL {r[0]} {r[1]}: expected {r[2]}, observed {r[3]}")
    return R
