#!/usr/bin/env python3
"""protein_attach.py — the MANUAL extra-sensitive step of the de novo families: attach loci that the RNA family rule
left out to an existing family by protein homology. Never part of the default pipeline (driver stage
`families-protein`, not in `all`); the RNA families are read, never modified; no two families are ever merged
(proposed merges are only reported). Pre-registered in docs/archive/2026-09/PREREG_protein_attach_2026-09-25.md (every constant
below is fixed there; changing one is off the pre-registration).

    python3 tools/protein_attach.py --fam PREFIX.fam --fasta GENOME.fa --out PREFIX.fam_protein [--threads 4]
                                    [--budget-s S] [--shard-size N] [--no-null] [--no-paf]

Input: the driver `families` stage outputs PREFIX.fam.{clusters.tsv, loci.gff3, loci.tsv[, loci.paf]}.
Steps (PREREG section 2):
  1. protein of each locus: the representative's exons (the positional exon sum) spliced from the genome,
     oriented on the transcribed strand; longest stop-to-stop frame (no ATG needed; an N codon ends a frame) among
     the frames less than half soft-masked, >= 100 aa (both strands only for a locus without one);
  2. BLASTP (bench/truth.py's binary and command, -evalue 1e-5) of every locus protein against the proteins of the
     family MEMBERS, in resumable query shards (exit 75 = shards remain, re-run the same command);
  3. base hit (truth.edges_from style): greedy non-overlapping HSPs by bit score on the longer protein, coverage
     >= 0.30 of it; identity = sum nident / sum length of the chosen HSPs;
  4. family calibration: I_F = min over F's members of their best base-hit identity to another member of F (F needs
     >= 2 such members, else it is uncalibrated and receives nothing);
  5. an unattached locus u joins F* (the family of its best BLASTP hit) iff F* is calibrated, u has a base hit to
     F* with identity >= max(I_F*, 0.60), and no other family passes the same test (else `ambiguous`);
  6. null arm: the same with every candidate protein shuffled (fixed per-locus seed).
Outputs (OUT = --out): .orfs.tsv .proteins.faa .members.faa .calibration.tsv .candidates.tsv .candidate_hits.tsv
(every candidate x family it hits) .attached.tsv
.clusters.tsv (the RNA rows unchanged + one row per attached locus, column added_by) .merges.tsv (proposed, never
applied) .null.tsv .params.tsv, and the BLASTP shard directories OUT.blastp.shards/ OUT.null.blastp.shards/.
"""
from __future__ import annotations

import argparse
import collections
import hashlib
import json
import os
import random
import re
import shutil
import subprocess
import sys
import time
import zlib

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(REPO, "bench"))
import lib    # noqa: E402  (CODON, rc, merge)
import truth  # noqa: E402  (BLAST, blastp_version, the command and the edge rule this step mirrors)

# ---------------------------------------------------------------- pre-registered constants (PREREG section 2)
MIN_ORF_AA = 100          # r1027's floor (its codon-shuffled null gave 0 edges at this floor)
MAX_MASKED_FRAC = 0.50    # a frame at least half soft-masked is transposable-element sequence
MIN_COV = truth.MIN_COV   # 0.30 of the longer protein (§6ko, unchanged; r1030)
SCOPE_FLOOR = 0.60        # §6t1 divergence cliff of RNA edges; r1096 protein-tier precision below it
EVALUE = "1e-5"           # every protein arm of the register
NULL_SEED = 20260925
ORF_VERSION = "1"
BLAST_FMT = "6 qseqid sseqid nident length qstart qend sstart send bitscore"
VERDICTS = ["attached", "ambiguous", "below_scope_floor", "below_family_identity", "below_coverage",
            "family_uncalibrated", "no_hit", "overlaps_member"]


def log(msg):
    print(f"[protein_attach] {msg}", file=sys.stderr, flush=True)


def fingerprint(path):
    try:
        st = os.stat(path)
        return f"{os.path.abspath(path)}\t{st.st_size}\t{int(st.st_mtime)}"
    except OSError:
        return f"{path}\tabsent"


def write_atomic(path, text):
    with open(path + ".tmp", "w") as fh:
        fh.write(text)
    os.replace(path + ".tmp", path)


def write_if_changed(path, text):
    """Keep the file (and its mtime / md5-keyed caches) when the content is unchanged."""
    if os.path.exists(path) and open(path).read() == text:
        return
    write_atomic(path, text)


def span_key(s):
    c, se = s.rsplit(":", 1)
    a, b = se.split("-")
    return c, int(a), int(b)


# ---------------------------------------------------------------- inputs: loci, families
def read_loci(gff3):
    """loci.gff3 -> list of dicts (order of the file): id, contig, start, end (1-based closed), strand, exons."""
    loci, by_id = [], {}
    for line in open(gff3):
        if line.startswith("#"):
            continue
        f = line.rstrip("\n").split("\t")
        if len(f) < 9:
            continue
        if f[2] == "gene":
            m = re.search(r"ID=gene-([^;]+)", f[8])
            d = {"id": m.group(1), "contig": f[0], "start": int(f[3]), "end": int(f[4]), "strand": f[6],
                 "exons": []}
            by_id[d["id"]] = d
            loci.append(d)
        elif f[2] == "exon":
            m = re.search(r"Parent=gene-([^;]+)", f[8])
            by_id[m.group(1)]["exons"].append((int(f[3]), int(f[4])))
    for i, d in enumerate(loci):
        d["exons"].sort()
        d["idx"] = i
    return loci


def read_families(clusters, loci_tsv):
    """{(contig, start, end): family} for member spans, the fold map, the RNA header and rows (verbatim)."""
    rows = open(clusters).read().splitlines()
    head = rows[0].split("\t")
    ci = {h: i for i, h in enumerate(head)}
    fam_of_span = {}
    for r in rows[1:]:
        f = r.split("\t")
        if len(f) < len(head):
            continue
        k = (f[ci["chrom"]], int(f[ci["start"]]), int(f[ci["end"]]))
        if k in fam_of_span and fam_of_span[k] != f[ci["cluster_id"]]:
            raise RuntimeError(f"{clusters}: span {k} in two families")
        fam_of_span[k] = f[ci["cluster_id"]]
    fold = {}
    if loci_tsv and os.path.exists(loci_tsv):
        for line in open(loci_tsv).read().splitlines()[1:]:
            a, r = line.split("\t")
            fold[span_key(a)] = span_key(r)
    return fam_of_span, fold, rows


class ExonIndex:
    """Representative exons of a set of loci, binned for overlap queries (1-based closed intervals)."""
    B = 100_000

    def __init__(self, loci):
        self.bins = collections.defaultdict(list)
        for d in loci:
            for a, b in d["exons"]:
                for k in range(a // self.B, b // self.B + 1):
                    self.bins[(d["contig"], k)].append((a, b, d["idx"]))

    def hits(self, d):
        out = set()
        for a, b in d["exons"]:
            for k in range(a // self.B, b // self.B + 1):
                for x, y, j in self.bins.get((d["contig"], k), ()):
                    if j != d["idx"] and x <= b and a <= y:
                        out.add(j)
        return out


def exons_overlap(p, q):
    if p["contig"] != q["contig"] or p["end"] < q["start"] or q["end"] < p["start"]:
        return False
    i = j = 0
    while i < len(p["exons"]) and j < len(q["exons"]):
        (a, b), (x, y) = p["exons"][i], q["exons"][j]
        if a <= y and x <= b:
            return True
        if b < y:
            i += 1
        else:
            j += 1
    return False


# ---------------------------------------------------------------- step 1: the protein of a locus
def orf_runs(seq):
    """Every maximal stop-free run of the three forward frames of an oriented, soft-masked sequence:
    (aa_len, frame, aa_start, aa_end, masked_frac, protein)."""
    up = seq.upper()
    cum = [0]
    for ch in seq:
        cum.append(cum[-1] + (ch.islower()))
    out = []
    for fr in range(3):
        aa = [lib.CODON.get(up[i:i + 3], "*") for i in range(fr, len(up) - 2, 3)]
        start = 0
        for j in range(len(aa) + 1):
            if j == len(aa) or aa[j] == "*":
                if j > start:
                    b0, b1 = fr + 3 * start, fr + 3 * j
                    out.append((j - start, fr, start, j, (cum[b1] - cum[b0]) / (b1 - b0), "".join(aa[start:j])))
                start = j + 1
    return out


def locus_protein(fa, d):
    """(status, protein, info) for one locus (PREREG 2.2)."""
    seq = "".join(fa.fetch(d["contig"], a - 1, b) for a, b in d["exons"])
    strands = ["+", "-"] if d["strand"] not in "+-" else [d["strand"]]
    runs = []
    for st in strands:
        s = lib.rc(seq) if st == "-" else seq
        runs += [(r, st) for r in orf_runs(s)]
    ok = [(r, st) for r, st in runs if r[4] < MAX_MASKED_FRAC]
    best = max(ok, key=lambda x: x[0][0], default=None)      # first maximum: frame order, then position
    longest_any = max((r[0] for r, _ in runs), default=0)
    info = {"exonic_nt": len(seq), "longest_any_aa": longest_any, "orf_aa": best[0][0] if best else 0,
            "orf_masked": best[0][4] if best else float("nan"), "orf_strand": best[1] if best else ".",
            "orf_frame": best[0][1] if best else -1}
    if best and best[0][0] >= MIN_ORF_AA:
        return "ok", best[0][5], info
    if any(r[0] >= MIN_ORF_AA for r, _ in runs):
        return "te_majority", None, info
    if best and best[0][0] >= 50:
        return "short_50_99", None, info
    return "short", None, info


# ---------------------------------------------------------------- step 2: sharded BLASTP (query set vs member db)
def blastp_vs(query_faa, db_prefix, db_md5, prefix, threads, shard_size, budget_s, t0, tag):
    """truth.blastp_sharded with a separate database: contiguous query shards PREFIX.blastp.shards/<s>-<e>.tsv
    written through .tmp + rename, keyed by a manifest (md5 of the queries, of the database FASTA, BLASTP
    version, arguments); a changed key discards the shards. --budget-s: no shard starts that the slowest shard
    so far would not finish in time; a running shard is killed at the budget. Returns (paths, done, n)."""
    sd = prefix + ".blastp.shards"
    recs, name = [], None
    for line in open(query_faa):
        if line.startswith(">"):
            name = line[1:].strip()
        else:
            recs.append((name, line.strip()))
    n = len(recs)
    cmd = ["-evalue", EVALUE, "-max_target_seqs", "100000", "-num_threads", str(threads), "-outfmt", BLAST_FMT]
    key = {"query_md5": hashlib.md5(open(query_faa, "rb").read()).hexdigest(), "db_md5": db_md5, "n_queries": n,
           "blastp": truth.blastp_version(), "args": cmd[:4] + cmd[6:]}
    man = os.path.join(sd, "manifest.json")
    old = json.load(open(man)) if os.path.exists(man) else None
    if old != key:
        if old is not None:
            log(f"{sd}: cache key changed; discarding its shards")
        shutil.rmtree(sd, ignore_errors=True)
        os.makedirs(sd)
        write_atomic(man, json.dumps(key, indent=1, sort_keys=True))
    have = {}
    for fn in os.listdir(sd):
        m = re.fullmatch(r"(\d+)-(\d+)\.tsv", fn)
        if m:
            have[int(m.group(1))] = (int(m.group(2)), os.path.join(sd, fn))
    times = os.path.join(sd, "times.tsv")
    slowest = max((float(l.split("\t")[2]) for l in open(times) if not l.startswith("start")), default=0.0) \
        if os.path.exists(times) else 0.0
    done, start = [], 0
    while start < n:
        if start in have:
            end, path = have[start]
            done.append(path)
            start = end
            continue
        end = min(n, start + shard_size)
        timeout = None
        if budget_s:
            elapsed = time.time() - t0
            if elapsed + slowest > budget_s:
                log(f"budget: {elapsed:.0f} s used, slowest shard {slowest:.0f} s -> stopping before {tag} "
                    f"queries {start}-{end}")
                break
            timeout = budget_s - elapsed
        out = os.path.join(sd, f"{start:07d}-{end:07d}.tsv")
        qf = out[:-4] + ".query.faa"
        with open(qf, "w") as fh:
            for nm, sq in recs[start:end]:
                fh.write(f">{nm}\n{sq}\n")
        ts = time.time()
        try:
            with open(out + ".tmp", "w") as fh:
                subprocess.run([truth.BLAST + "/blastp", "-query", qf, "-db", db_prefix] + cmd, stdout=fh,
                               check=True, timeout=timeout)
        except subprocess.TimeoutExpired:
            os.remove(out + ".tmp")
            log(f"budget: {tag} shard {start}-{end} killed after {time.time() - ts:.0f} s; it restarts next call")
            break
        wall = time.time() - ts
        os.replace(out + ".tmp", out)
        os.remove(qf)
        new = not os.path.exists(times)
        with open(times, "a") as fh:
            if new:
                fh.write("start\tend\twall_s\tbytes\n")
            fh.write(f"{start}\t{end}\t{wall:.1f}\t{os.path.getsize(out)}\n")
        slowest = max(slowest, wall)
        done.append(out)
        log(f"{tag} queries {start}-{end} of {n}: {wall:.0f} s")
        start = end
    return done, start, n


def make_db(members_faa, prefix):
    """makeblastdb of the member proteins, rebuilt only when the FASTA changes (md5 stamp)."""
    md5 = hashlib.md5(open(members_faa, "rb").read()).hexdigest()
    stamp = prefix + "_db.md5"
    if not (os.path.exists(stamp) and open(stamp).read() == md5):
        for p in os.listdir(os.path.dirname(os.path.abspath(prefix))):
            if p.startswith(os.path.basename(prefix) + "_db."):
                os.remove(os.path.join(os.path.dirname(os.path.abspath(prefix)), p))
        subprocess.run([truth.BLAST + "/makeblastdb", "-dbtype", "prot", "-in", members_faa, "-out", prefix + "_db"],
                       stdout=subprocess.DEVNULL, check=True)
        write_atomic(stamp, md5)
    return prefix + "_db", md5


# ---------------------------------------------------------------- step 3: base hits, streamed per query block
def pair_stats(rows, lq, ls):
    """truth.edges_from's rule for one ordered pair: greedy non-overlapping HSPs by bit score on the LONGER protein
    (half-open), coverage = their union / longer length, identity = sum nident / sum length; plus summed bits."""
    longer_is_q = lq >= ls
    taken, N, L, B = [], 0, 0, 0.0
    for bits, nid, ln, q0, q1, s0, s1 in sorted(rows, reverse=True):
        iv = (q0, q1) if longer_is_q else (s0, s1)
        if any(iv[0] < y and x < iv[1] for x, y in taken):
            continue
        taken.append(iv)
        N += nid
        L += ln
        B += bits
    cov = sum(y - x for x, y in lib.merge(taken)) / max(lq, ls)
    return (N / L if L else 0.0), cov, B


def stream_hits(paths, plen_q, plen_s):
    """Yield (query id, {subject id: (identity, coverage, bits)}) per query block of the concatenated shards."""
    cur, block = None, None
    seen = set()
    for path in paths:
        for line in open(path):
            q, s, nid, ln, q0, q1, s0, s1, bits = line.rstrip("\n").split("\t")
            if q != cur:
                if cur is not None:
                    yield cur, {k: pair_stats(v, plen_q[cur], plen_s[k]) for k, v in block.items()}
                    seen.add(cur)
                if q in seen:
                    raise ValueError(f"{path}: query {q} has two output blocks")
                cur, block = q, collections.defaultdict(list)
            if q == s:
                continue
            block[s].append((float(bits), int(nid), int(ln), int(q0) - 1, int(q1), int(s0) - 1, int(s1)))
    if cur is not None:
        yield cur, {k: pair_stats(v, plen_q[cur], plen_s[k]) for k, v in block.items()}


# ---------------------------------------------------------------- steps 4-6
def per_family(hits, fam_of, member_locus, loci, u):
    """{family: [best bits (any hit), best base identity, its coverage, its bits, its member]} of one query."""
    agg = {}
    for s, (ident, cov, bits) in hits.items():
        m = member_locus[s]
        if exons_overlap(loci[u], loci[m]):
            continue
        F = fam_of[m]
        a = agg.setdefault(F, [0.0, None, None, None, None, None])
        if bits > a[0]:
            a[0], a[5] = bits, m
        if cov >= MIN_COV and (a[1] is None or ident > a[1] or (ident == a[1] and bits > a[3])):
            a[1], a[2], a[3], a[4] = ident, cov, bits, m
    return agg


def verdict(agg, calib):
    """PREREG 2.6 steps 1-6 on one query's per-family aggregate -> (verdict, F*, n families passing, other passes)."""
    if not agg:
        return "no_hit", None, 0, False

    def passes(F):
        a = agg[F]
        return F in calib and a[1] is not None and a[1] >= max(calib[F], SCOPE_FLOOR)
    fstar = max(agg, key=lambda F: (agg[F][0], F))
    passing = [F for F in agg if passes(F)]
    other = any(F != fstar for F in passing)
    a = agg[fstar]
    if fstar not in calib:
        return "family_uncalibrated", fstar, len(passing), other
    if a[1] is None:
        return "below_coverage", fstar, len(passing), other
    if a[1] < calib[fstar]:
        return "below_family_identity", fstar, len(passing), other
    if a[1] < SCOPE_FLOOR:
        return "below_scope_floor", fstar, len(passing), other
    if other:
        return "ambiguous", fstar, len(passing), other
    return "attached", fstar, len(passing), other


def shuffled(pid, prot):
    r = random.Random(zlib.crc32(pid.encode()) ^ NULL_SEED)
    p = list(prot)
    r.shuffle(p)
    return "".join(p)


def nt_records(paf, pairs_wanted):
    """{(u, family)} for which PREFIX.fam.loci.paf holds any record between locus u and a member of the family.
    pairs_wanted: {span string: [(role, locus idx)]} built by the caller."""
    found = set()
    with open(paf) as fh:
        for line in fh:
            q, _, _, _, _, t = line.split("\t", 6)[:6]
            if q == t:
                continue
            for a, b in ((q, t), (t, q)):
                for ua in pairs_wanted["u"].get(a, ()):
                    for F in pairs_wanted["m"].get(b, ()):
                        found.add((ua, F))
    return found


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0],
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--fam", required=True, help="the families stage prefix (PREFIX.fam: .clusters.tsv, .loci.gff3, "
                    ".loci.tsv, .loci.paf)")
    ap.add_argument("--fasta", required=True, help="the genome FASTA (soft-masked; indexed)")
    ap.add_argument("--out", required=True, help="output prefix (e.g. PREFIX.fam_protein)")
    ap.add_argument("--threads", type=int, default=4)
    ap.add_argument("--shard-size", type=int, default=1000, help="queries per BLASTP shard (changing it keeps "
                    "finished shards)")
    ap.add_argument("--budget-s", type=float, default=0.0, help="wall-clock budget of this call; exit 75 = shards "
                    "remain, run the same command again")
    ap.add_argument("--no-null", action="store_true", help="skip the shuffled-protein null arm (PREREG 2.7 runs it)")
    ap.add_argument("--no-paf", action="store_true", help="skip the nucleotide-record label (large genome-wide PAF)")
    a = ap.parse_args(argv)
    t0 = time.time()
    import pysam
    P = a.fam
    clusters, gff3, loci_tsv, paf = P + ".clusters.tsv", P + ".loci.gff3", P + ".loci.tsv", P + ".loci.paf"
    for p in (clusters, gff3):
        if not os.path.exists(p):
            sys.exit(f"protein_attach: {p} missing (run the driver's families stage first)")
    os.makedirs(os.path.dirname(os.path.abspath(a.out)), exist_ok=True)

    # ---- loci and families
    loci = read_loci(gff3)
    fam_of_span, fold, rna_rows = read_families(clusters, loci_tsv)
    fam_of = {}
    for d in loci:
        k = (d["contig"], d["start"], d["end"])
        k = k if k in fam_of_span else fold.get(k)
        if k in fam_of_span:
            fam_of[d["idx"]] = fam_of_span[k]
    families = sorted(set(fam_of_span.values()), key=lambda x: (len(x), x))
    members_by_fam = collections.defaultdict(list)
    for i, F in fam_of.items():
        members_by_fam[F].append(i)
    log(f"{len(loci):,} loci, {len(families)} families, {len(fam_of):,} member loci")

    # ---- proteins (cached on the inputs' fingerprints)
    orfs_tsv, faa, mfaa = a.out + ".orfs.tsv", a.out + ".proteins.faa", a.out + ".members.faa"
    okey = f"v{ORF_VERSION}\n{fingerprint(gff3)}\n{fingerprint(clusters)}\n{fingerprint(loci_tsv)}\n" \
           f"{fingerprint(a.fasta)}\n{MIN_ORF_AA}\t{MAX_MASKED_FRAC}\n"
    kpath = a.out + ".orfs.key"
    prot, status = {}, {}
    if os.path.exists(kpath) and open(kpath).read() == okey and os.path.exists(orfs_tsv) and os.path.exists(faa):
        name = None
        for line in open(faa):
            if line.startswith(">"):
                name = int(line[2:].strip())
            else:
                prot[name] = line.strip()
        for r in open(orfs_tsv).read().splitlines()[1:]:
            f = r.split("\t")
            status[int(f[0])] = f[-1]
        log(f"proteins reused from {faa}")
    else:
        fa = pysam.FastaFile(a.fasta)
        rows = ["idx\tlocus\tcontig\tstart\tend\tstrand\tfamily\tn_exons\texonic_nt\torf_strand\torf_frame\torf_aa"
                "\torf_masked_frac\tlongest_any_aa\tstatus"]
        for d in loci:
            st, p, info = locus_protein(fa, d)
            status[d["idx"]] = st
            if p:
                prot[d["idx"]] = p
            rows.append("\t".join(map(str, [d["idx"], d["id"], d["contig"], d["start"], d["end"], d["strand"],
                                            fam_of.get(d["idx"], "-"), len(d["exons"]), info["exonic_nt"],
                                            info["orf_strand"], info["orf_frame"], info["orf_aa"],
                                            f"{info['orf_masked']:.3f}", info["longest_any_aa"], st])))
        write_if_changed(faa, "".join(f">L{i}\n{p}\n" for i, p in sorted(prot.items())))
        write_atomic(orfs_tsv, "\n".join(rows) + "\n")
        write_atomic(kpath, okey)
        log(f"proteins: {len(prot):,} of {len(loci):,} loci have a >= {MIN_ORF_AA} aa frame "
            f"({time.time() - t0:.0f} s)")
    write_if_changed(mfaa, "".join(f">L{i}\n{prot[i]}\n" for i in sorted(fam_of) if i in prot))
    member_idx = sorted(i for i in fam_of if i in prot)
    if not member_idx:
        sys.exit("protein_attach: no family member has a protein; nothing to search against")

    # candidate queries: unattached loci with a protein whose exons touch no member's exons
    mindex = ExonIndex([loci[i] for i in fam_of])
    overlapping = {d["idx"] for d in loci if d["idx"] not in fam_of and mindex.hits(d)}
    candidates = [i for i in sorted(prot) if i not in fam_of and i not in overlapping]
    null_faa = a.out + ".null.proteins.faa"
    write_if_changed(null_faa, "".join(f">L{i}\n{shuffled(f'L{i}', prot[i])}\n" for i in candidates))

    # ---- BLASTP (real: every protein; null: shuffled candidates), resumable
    dbp, dbmd5 = make_db(mfaa, a.out + ".members")
    paths, done, n = blastp_vs(faa, dbp, dbmd5, a.out, a.threads, a.shard_size, a.budget_s, t0, "real")
    npaths, ndone, nn = [], 0, 0
    if done == n and not a.no_null and candidates:
        npaths, ndone, nn = blastp_vs(null_faa, dbp, dbmd5, a.out + ".null", a.threads, a.shard_size, a.budget_s,
                                      t0, "null")
    if done < n or ndone < nn:
        print(f"protein_attach: INCOMPLETE — real {done}/{n}, null {ndone}/{nn} query proteins searched; run the "
              f"same command again")
        sys.exit(75)

    # ---- hits
    plen = {f"L{i}": len(p) for i, p in prot.items()}
    member_locus = {f"L{i}": i for i in member_idx}
    sib = collections.defaultdict(float)              # member -> best within-family base identity (-1 = none)
    sib_has = set()
    cross = collections.defaultdict(list)             # (A, B) -> [(identity, coverage, bits, a, b)]
    cand_agg = {}
    cand_set = set(candidates)
    for q, hits in stream_hits(paths, plen, plen):
        qi = int(q[1:])
        if qi in fam_of:
            A = fam_of[qi]
            for s, (ident, cov, bits) in hits.items():
                si = member_locus[s]
                if cov < MIN_COV or exons_overlap(loci[qi], loci[si]):
                    continue
                if fam_of[si] == A:
                    for x in (qi, si):
                        if x not in sib_has or ident > sib[x]:
                            sib[x] = ident
                            sib_has.add(x)
                elif ident >= SCOPE_FLOOR:
                    B = fam_of[si]
                    key = (A, B) if (len(A), A) < (len(B), B) else (B, A)
                    cross[key].append((ident, cov, bits, qi, si))
        elif qi in cand_set:
            cand_agg[qi] = per_family(hits, fam_of, member_locus, loci, qi)
    calib, calib_rows = {}, ["family\tmembers\tmembers_with_protein\tmembers_calibrated\tI_F\tcalibrated\tloosest_member"]
    for F in families:
        mem = members_by_fam[F]
        cal = [x for x in mem if x in sib_has]
        if len(cal) >= 2:
            loosest = min(cal, key=lambda x: (sib[x], x))
            calib[F] = sib[loosest]
        calib_rows.append("\t".join(map(str, [F, len(mem), sum(1 for x in mem if x in prot), len(cal),
                                              f"{calib[F]:.4f}" if F in calib else "NA",
                                              "yes" if F in calib else "no",
                                              loci[loosest]["id"] if F in calib else "-"])))
    write_atomic(a.out + ".calibration.tsv", "\n".join(calib_rows) + "\n")

    # ---- verdicts
    res = {}
    for i in candidates:
        agg = cand_agg.get(i, {})
        res[i] = (verdict(agg, calib), agg)
    null_att = []
    if npaths:
        for q, hits in stream_hits(npaths, plen, plen):
            qi = int(q[1:])
            agg = per_family(hits, fam_of, member_locus, loci, qi)
            v, F, _, _ = verdict(agg, calib)
            if v == "attached":
                null_att.append((qi, F, agg[F]))

    # ---- nucleotide records (sub-threshold alignments) between candidates and their best family's members
    ntfound = None
    if not a.no_paf and os.path.exists(paf):
        span = lambda d: f"{d['contig']}:{d['start']}-{d['end']}"   # noqa: E731  (loci.fa headers)
        want = {"u": collections.defaultdict(list), "m": collections.defaultdict(list)}
        for i, ((v, F, _, _), _) in res.items():
            if F is not None:
                want["u"][span(loci[i])].append(i)
        for i, F in fam_of.items():
            want["m"][span(loci[i])].append(F)
        found = nt_records(paf, want)
        ntfound = {i for (i, F) in found if res[i][0][1] == F}
        log(f"nucleotide records read from {paf}")

    # ---- outputs
    def nt(i):
        return "NA" if ntfound is None else ("yes" if i in ntfound else "no")
    crow = ["locus\tcontig\tstart\tend\tstrand\taa\tbest_family\tbest_hit_member\tbest_bits\tfamily_calibrated\tI_F"
            "\tbase_identity\tbase_coverage\tbase_bits\tbase_member\tfamilies_passing\tother_family_passes"
            "\tnt_record\tverdict"]
    arow = ["locus\tcontig\tstart\tend\tstrand\tfamily\tmember\tidentity\tcoverage\tbits\tI_F\tevidence"]
    counts = collections.Counter()
    ext = [rna_rows[0] + "\tadded_by"] + [r + "\trna" for r in rna_rows[1:]]
    for i in candidates:
        (v, F, npass, other), agg = res[i]
        d = loci[i]
        counts[v] += 1
        ag = agg.get(F) if F else None
        crow.append("\t".join(map(str, [
            d["id"], d["contig"], d["start"], d["end"], d["strand"], len(prot[i]), F or "-",
            loci[ag[5]]["id"] if ag else "-", f"{ag[0]:.1f}" if ag else "NA",
            ("yes" if F in calib else "no") if F else "NA", f"{calib[F]:.4f}" if F in calib else "NA",
            f"{ag[1]:.4f}" if ag and ag[1] is not None else "NA", f"{ag[2]:.4f}" if ag and ag[2] is not None else "NA",
            f"{ag[3]:.1f}" if ag and ag[3] is not None else "NA", loci[ag[4]]["id"] if ag and ag[4] is not None else "-",
            npass, "yes" if other else "no", nt(i) if F else "NA", v])))
        if v == "attached":
            ev = "NA" if ntfound is None else ("nucleotide-corroborated" if i in ntfound else "protein-only")
            arow.append("\t".join(map(str, [d["id"], d["contig"], d["start"], d["end"], d["strand"], F,
                                            loci[ag[4]]["id"], f"{ag[1]:.4f}", f"{ag[2]:.4f}", f"{ag[3]:.1f}",
                                            f"{calib[F]:.4f}", ev])))
            ext.append("\t".join([F, "NA", "NA", "NA", "NA", d["contig"], str(d["start"]), str(d["end"]), "protein"]))
    for i in sorted(overlapping):
        counts["overlaps_member"] += 1
    write_atomic(a.out + ".candidates.tsv", "\n".join(crow) + "\n")
    hrow = ["locus\tfamily\tbest_bits\tbase_identity\tbase_coverage\tbase_member\tpasses"]
    for i in candidates:
        for F, ag in sorted(res[i][1].items(), key=lambda x: -x[1][0]):
            ok = F in calib and ag[1] is not None and ag[1] >= max(calib[F], SCOPE_FLOOR)
            hrow.append("\t".join(map(str, [loci[i]["id"], F, f"{ag[0]:.1f}",
                                            f"{ag[1]:.4f}" if ag[1] is not None else "NA",
                                            f"{ag[2]:.4f}" if ag[2] is not None else "NA",
                                            loci[ag[4]]["id"] if ag[4] is not None else "-", "yes" if ok else "no"])))
    write_atomic(a.out + ".candidate_hits.tsv", "\n".join(hrow) + "\n")
    write_atomic(a.out + ".attached.tsv", "\n".join(arow) + "\n")
    write_atomic(a.out + ".clusters.tsv", "\n".join(ext) + "\n")

    mrow = ["family_a\tfamily_b\tvia\tn_pairs\tbest_identity\tbest_coverage\tlocus_a\tlocus_b\tI_a\tI_b\tstatus"]
    n_merge = 0
    for (A, B), prs in sorted(cross.items()):
        need = max([SCOPE_FLOOR] + [calib[x] for x in (A, B) if x in calib])
        ok = sorted((p for p in prs if p[0] >= need), reverse=True)
        pairs = {(min(p[3], p[4]), max(p[3], p[4])) for p in ok}
        if ok:
            b = ok[0]
            la, lb = (b[3], b[4]) if fam_of[b[3]] == A else (b[4], b[3])
            n_merge += 1
            mrow.append("\t".join(map(str, [A, B, "member_pair", len(pairs), f"{b[0]:.4f}", f"{b[1]:.4f}",
                                            loci[la]["id"], loci[lb]["id"],
                                            f"{calib[A]:.4f}" if A in calib else "NA",
                                            f"{calib[B]:.4f}" if B in calib else "NA", "proposed, never applied"])))
    for i in candidates:
        (v, F, _, _), agg = res[i]
        if v != "ambiguous":
            continue
        fams = sorted((G for G in agg if G in calib and agg[G][1] is not None
                       and agg[G][1] >= max(calib[G], SCOPE_FLOOR)), key=lambda G: -agg[G][1])
        mrow.append("\t".join(map(str, [fams[0], ",".join(fams[1:]), "ambiguous_locus", 1, f"{agg[fams[0]][1]:.4f}",
                                        f"{agg[fams[0]][2]:.4f}", loci[i]["id"], "-",
                                        f"{calib[fams[0]]:.4f}", "NA", "proposed, never applied"])))
    write_atomic(a.out + ".merges.tsv", "\n".join(mrow) + "\n")
    nrow = ["locus\tfamily\tmember\tidentity\tcoverage\tbits"]
    for qi, F, ag in null_att:
        nrow.append(f"{loci[qi]['id']}\t{F}\t{loci[ag[4]]['id']}\t{ag[1]:.4f}\t{ag[2]:.4f}\t{ag[3]:.1f}")
    write_atomic(a.out + ".null.tsv", "\n".join(nrow) + "\n")

    st = collections.Counter(status.values())
    params = [("prereg", "docs/archive/2026-09/PREREG_protein_attach_2026-09-25.md"), ("fam", P), ("fasta", a.fasta),
              ("min_orf_aa", MIN_ORF_AA), ("max_masked_frac", MAX_MASKED_FRAC), ("min_cov_longer", MIN_COV),
              ("scope_floor", SCOPE_FLOOR), ("evalue", EVALUE), ("blastp", truth.blastp_version()),
              ("loci", len(loci)), ("families", len(families)), ("member_loci", len(fam_of)),
              ("loci_protein_ok", st["ok"]), ("loci_te_majority", st["te_majority"]),
              ("loci_orf_50_99", st["short_50_99"]), ("loci_orf_lt50", st["short"]),
              ("member_loci_with_protein", len(member_idx)), ("families_calibrated", len(calib)),
              ("unattached_loci", len(loci) - len(fam_of)), ("candidates", len(candidates))]
    params += [(f"verdict_{v}", counts[v]) for v in VERDICTS]
    n_att_fam = len({res[i][0][1] for i in candidates if res[i][0][0] == "attached"})
    params += [("attached_families", n_att_fam),
               ("attached_nt_corroborated", "NA" if ntfound is None else
                sum(1 for i in candidates if res[i][0][0] == "attached" and i in ntfound)),
               ("null_arm", "skipped" if a.no_null else "run"), ("null_attached", len(null_att)),
               ("proposed_merges_member_pairs", n_merge), ("proposed_merges_ambiguous_loci", counts["ambiguous"]),
               ("wall_s", f"{time.time() - t0:.0f}")]
    write_atomic(a.out + ".params.tsv", "".join(f"{k}\t{v}\n" for k, v in params))
    print(f"protein_attach: {counts['attached']} of {len(candidates):,} candidate loci attached to "
          f"{n_att_fam} families ({len(calib)} of {len(families)} families calibrated); null {len(null_att)}; "
          f"{n_merge} proposed merges (never applied) -> {a.out}.attached.tsv")


if __name__ == "__main__":
    main()
