#!/usr/bin/env python3
"""Assembly parCN / famCN for Soto 2025's SD98 genes: QuicK-mer2's own k-mer rule, counted EXACTLY in complete
assemblies (2026-09-29, KEY=soto_parcn_asm; docs/archive/2026-09/PREREG_soto_parcn_assembly_2026-09-29.md, frozen sha1 2a76a930;
register rows 1171-1174).

WHY. Soto's paralog-specific copy number (Table S1E "Median parCN") is QuicK-mer2 on 2,504 1KGP short-read genomes
against CHM13 v1.0. Its whole-genome index needs a 2^32-slot hash: `search` peaks at ~52 GB and `count` at ~43 GB,
against 25 GiB + 16 GB swap here (register 1167; figs/soto_quickmer2.md), so the read-depth pipeline cannot be
re-run on this machine. What CAN be computed exactly is the same k-mer set counted in complete assemblies (HG002 v1.1
diploid, CHM13, chimpanzee / gorilla / orangutan), which is what this module does. Human and ape numbers are reported
per genome and never pooled.

THE RULE (QuicKmer.c, `-e 1 -w 500 -k 30` as Soto ran it):
  * k = 30, CANONICAL (min of the k-mer and its reverse complement in the A0 C1 G2 T3 code), case-insensitive
    (soft-masked bases count), non-ACGT resets the window, the all-A k-mer (code 0) is skipped.
  * a k-mer is PARALOG-SPECIFIC (SPEC) iff it occurs exactly once in the reference (CHM13 v2.0 WITHOUT chrY = v1.1
    content; Soto's v1.0 has no Y) AND its "edit depth" -- the sum over its 90 one-substitution neighbours of
    min(reference occurrences, 255) -- is < 100. The edit-depth filter drops 0.44% of the once-only k-mers.
  * parCN_G(gene) = s_G x median over SPEC(gene) of the exact count in genome G (s = 1 for a diploid sum such as
    HG002 MAT+PAT or gorilla mat+pat, 2 for a haploid primary: chimp, orangutan, CHM13). RESOLVED iff |SPEC| >= 100
    (one fifth of QuicK-mer2's 500-k-mer window), else reported as unresolved, never imputed.
  * famCN analogue (for their WSSD human-vs-ape calls only): FAM(gene) = the region's 30-mers with no soft-masked
    base in CHM13 v2.0 (RepeatMasker / TRF ~ WSSD's masked reference), any reference count; p_G = share present
    (count >= 1) in G; famCN_G = s_G x median count over the present positions if p_G >= 0.10 (an exact 30-mer
    survives divergence d with probability ~(1-d)^30), 0 if p_G < 0.10 (absent), unresolved if |FAM| < 100.
  * ASM-H (assembly human gain) iff famCN_HG002 > max(famCN over the apes); "duplicated" if that ape max < 2.5,
    else "expanded". Regions = S1E `SD98_v2.0` (gene ∩ SD98, BED half-open; 1,831 rows) + 300 RefSeq protein-coding
    controls (seed 20260929, first 20 kb, >= 1 Mb from any S1E region, names not in S1E/S1C).

REPRODUCED NUMBERS (analyze on the frozen count tables, 2026-09-29; the full Outcome is in the prereg):
  C0 controls: HG002 parCN = 2 for 297/299 (0.993); famCN = 2 for 0.983 / 0.983 / 0.976 of 297 (chimp / gorilla /
  orangutan). C5 resolved 1,163 / 1,831 (0.635). C1 Fixed 321/322 = 0.997 (mean statistic 0.994). C2 Nearly-Fixed
  408/629 = 0.649 (bar 0.70: fails; mean 0.692). C3 Polymorphic 83/212 = 0.392. C4 Spearman 0.414. H1 Duplicated in
  humans 105/109 = 0.963; H2 Expanded in humans 105/132 = 0.795; H3 non-calls ASM-H 131/631 = 0.208; H4 families
  105/118 = 0.890. F-cal Spearman 0.617 (human) / 0.886 (ape max).

SUBCOMMANDS (in order; every product is a plain file under --out-prefix / --work)
    regions     S1E SD98_v2.0 rows + 300 controls -> regions.tsv (rid, kind, s1e_row, gene, gene_id, chrom, start,
                end; `S<i>` = S1E data row i, `C<j>` = control j)
    kmers       every 30-mer of every region from the CHM13 v2.0 FASTA -> Q.u64 (sorted unique canonical codes),
                pos.npz (per position: index into Q, soft-mask flag; region offsets), regions.tsv with n_pos
    count       exact counts of the Q k-mers in one assembly -> <name>.i32 (int32 per Q entry). Via meryl
                (/home/juanfra/miniforge3/bin/meryl, `count k=30 threads=4 memory=16`, then
                `print [intersect genome.meryl Q.meryl]`); an existing meryl DB is reused with --db. The whole-
                genome DB is the cost (~3 G distinct 30-mers per genome, 15-30 GB on disk, ~10-30 min): the six
                tables of the 09-29 run were made with a C hash counter (`kc30`, identical to a Python brute force
                on 46,932 test k-mers) and live under
                /mnt/linuxdisk/tmp/rustle_figures_dev/soto_parcn_asm/counts/ -- `analyze` reads them as they are.
    edit-depth  QuicK-mer2's `-e 1` edit depth of every once-only Q k-mer (candidates) -> ed.u32, by streaming the
                reference once per batch of candidates (numpy; ~10-20 min per 2.4 M candidates, 3 batches for
                CHM13). The frozen ed0-2.u32 (kn30, identical to brute force on 2,998 test candidates) are reused.
    analyze     per-region SPEC / FAM statistics, parCN / famCN per genome, calls and the pre-registered clauses
                -> results_s1e.tsv, results_ctrl.tsv, results_named.tsv, summary.json
    The module-level functions (canonical_kmers, stream_kmer_counts, edit_depth_batch, region_stats) are what
    bench/soto/test_parcn_assembly.py exercises on toy data (`python3 bench/soto/test_parcn_assembly.py`).

INPUTS. In repo: bench/soto/soto_parCN_S1E.tsv (sha1 0c4cb730), soto_famCN_S1C.tsv (d008a179). Not in repo
(docs/DATA.md): chm13v2.0.fa (soft-masked, .fai), Reference/chm13v2.0_RefSeq_full.gff.gz (controls), the
assemblies (HG002 v1.1 `hg002v1.1.fasta.gz`; gorilla mGorGor1 mat.fa / pat.fa; mPanTro3 v2.0 pri; mPonPyg2 v2.0
pri), and the frozen k-mer products (Q.u64, pos.npz, regions.tsv, counts/*.i32, counts/ed*.u32, work/cand*.u64).

ENVIRONMENT. numpy and scipy (spearmanr); samtools on PATH (faidx); meryl for `count`. No pandas.
"""
import argparse, csv, gzip, json, os, random, re, statistics as st, subprocess, sys
from collections import defaultdict

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
S1E = os.path.join(HERE, "soto_parCN_S1E.tsv")
S1C = os.path.join(HERE, "soto_famCN_S1C.tsv")
K = 30
KMASK = (1 << (2 * K)) - 1
MINK, PFLOOR = 100, 0.10
MERYL = "/home/juanfra/miniforge3/bin/meryl"
SCALE = {"CHM13": 2, "HG002": 1, "GGO": 1, "GGOmat": 2, "GGOpat": 2, "PTR": 2, "PPY": 2}
NAMED = {"SRGAP2": ["ID_462"], "ARHGAP11": ["ID_145"], "FAM72": ["ID_354"], "NOTCH2NL": ["ID_400"],
         "NPIP": ["ID_149", "ID_151", "ID_152", "ID_153", "ID_154", "ID_155"], "TBC1D3": ["ID_468", "ID_469"]}
ACRO = {"chr13", "chr14", "chr15", "chr21", "chr22"}

_LUT = np.full(256, 4, dtype=np.uint64)
for _ch, _v in zip("ACGTacgt", [0, 1, 2, 3, 0, 1, 2, 3]):
    _LUT[ord(_ch)] = _v
_DEC = np.frombuffer(b"ACGT", dtype=np.uint8)


# ---------------------------------------------------------------------------------------------------------
# k-mer codes
# ---------------------------------------------------------------------------------------------------------

def revcomp_codes(x):
    """Reverse complement of an array (or int) of 2-bit K-mer codes (A0 C1 G2 T3: complement = 3 - base)."""
    x = np.asarray(x, dtype=np.uint64)
    r = np.zeros_like(x)
    for _ in range(K):
        r = (r << np.uint64(2)) | (np.uint64(3) - (x & np.uint64(3)))
        x = x >> np.uint64(2)
    return r


def canonical(x):
    return np.minimum(np.asarray(x, dtype=np.uint64), revcomp_codes(x))


def canonical_kmers(seq, with_mask=False):
    """All K-mers of `seq` (bytes / str) that contain no non-ACGT base, as canonical codes in sequence order (the
    all-A code 0 is KEPT here; QuicK-mer2's skip is applied by the callers that need it). With with_mask, also a
    boolean per kept K-mer: touches a soft-masked (lower-case) base."""
    b = np.frombuffer(seq.encode() if isinstance(seq, str) else seq, dtype=np.uint8)
    L = len(b) - K + 1
    if L <= 0:
        e = np.zeros(0, dtype=np.uint64)
        return (e, np.zeros(0, dtype=bool)) if with_mask else e
    code = _LUT[b]
    bad = np.concatenate([[0], np.cumsum(code == 4)])
    ok = (bad[K:] - bad[:-K]) == 0
    cc = np.where(code == 4, 0, code).astype(np.uint64)
    comp = (np.uint64(3) - cc).astype(np.uint64)
    fw, rv = np.zeros(L, dtype=np.uint64), np.zeros(L, dtype=np.uint64)
    for j in range(K):
        fw = (fw << np.uint64(2)) | cc[j:j + L]
    for j in range(K - 1, -1, -1):
        rv = (rv << np.uint64(2)) | comp[j:j + L]
    can = np.minimum(fw, rv)[ok]
    if not with_mask:
        return can
    low = np.concatenate([[0], np.cumsum(b >= 97)])
    return can, ((low[K:] - low[:-K]) > 0)[ok]


def decode_codes(codes):
    """Array of codes -> (n, K) uint8 array of ACGT characters."""
    codes = np.asarray(codes, dtype=np.uint64)
    out = np.empty((len(codes), K), dtype=np.uint8)
    for j in range(K):
        out[:, j] = _DEC[((codes >> np.uint64(2 * (K - 1 - j))) & np.uint64(3)).astype(np.int64)]
    return out


def encode_strings(chars):
    """(n, K) uint8 ACGT array -> forward codes."""
    code = _LUT[chars].astype(np.uint64)
    fw = np.zeros(len(chars), dtype=np.uint64)
    for j in range(K):
        fw = (fw << np.uint64(2)) | code[:, j]
    return fw


def fasta_records(path, exclude=()):
    """Yield (name, sequence bytes) from a FASTA (.gz allowed), skipping names in `exclude`."""
    op = gzip.open if path.endswith(".gz") else open
    name, buf = None, []
    with op(path, "rb") as fh:
        for ln in fh:
            if ln.startswith(b">"):
                if name is not None and name not in exclude:
                    yield name, b"".join(buf)
                name, buf = ln[1:].split()[0].decode(), []
            else:
                buf.append(ln.rstrip(b"\r\n"))
    if name is not None and name not in exclude:
        yield name, b"".join(buf)


def stream_kmer_counts(fasta, targets, exclude=(), chunk=1 << 24, skip_zero=True):
    """Exact occurrence counts of the sorted unique canonical codes `targets` over every record of `fasta`
    (case-insensitive, non-ACGT resets, the all-A k-mer skipped as QuicK-mer2 does). int64 per target."""
    counts = np.zeros(len(targets), dtype=np.int64)
    if len(targets) == 0:
        return counts
    for _name, seq in fasta_records(fasta, exclude):
        for s in range(0, max(1, len(seq) - K + 1), chunk):
            codes = canonical_kmers(seq[s:s + chunk + K - 1])
            if skip_zero:
                codes = codes[codes != 0]
            pos = np.searchsorted(targets, codes)
            pos[pos >= len(targets)] = 0
            hit = targets[pos] == codes
            if hit.any():
                counts += np.bincount(pos[hit], minlength=len(targets))
    return counts


def variants_1sub(cands):
    """(n, 90) canonical codes of every one-substitution neighbour of each candidate (QuicKmer.c Recurse_edit
    order: position p = 0..29 from the low bits, base + 1, + 2, + 3 mod 4)."""
    cands = np.asarray(cands, dtype=np.uint64)
    out = np.empty((len(cands), 3 * K), dtype=np.uint64)
    for p in range(K):
        b = (cands >> np.uint64(2 * p)) & np.uint64(3)
        cleared = cands & ~(np.uint64(3) << np.uint64(2 * p))
        for e in (1, 2, 3):
            y = cleared | (((b + np.uint64(e)) & np.uint64(3)) << np.uint64(2 * p))
            out[:, 3 * p + e - 1] = canonical(y)
    return out


def edit_depth_batch(cands, fasta, exclude=()):
    """QuicK-mer2 `-e 1` edit depth of each candidate: sum over its 90 neighbours of min(occurrences in the
    reference, 255); a neighbour equal to another neighbour's canonical form is summed once per neighbour."""
    var = variants_1sub(cands)
    uniq, inv = np.unique(var.ravel(), return_inverse=True)
    occ = stream_kmer_counts(fasta, uniq, exclude, skip_zero=False)
    return np.minimum(occ, 255)[inv].reshape(len(cands), 3 * K).sum(1).astype(np.uint32)


# ---------------------------------------------------------------------------------------------------------
# regions, kmers
# ---------------------------------------------------------------------------------------------------------

def read_tsv(path):
    with open(path, newline="") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def fnum(s):
    try:
        return float(s)
    except (TypeError, ValueError):
        return float("nan")


def regions_table(s1e_rows, s1c_rows, gff, n_controls=300, seed=20260929, control_bp=20000, exclude_bp=1_000_000):
    """S1E rows with a `SD98_v2.0` coordinate (BED half-open) + `n_controls` autosomal RefSeq protein-coding genes
    (first `control_bp`, >= `exclude_bp` from every S1E region, names absent from S1E / S1C), drawn with
    random.Random(seed).sample from the (chrom, start, end, name)-sorted eligible list."""
    rows = []
    for i, r in enumerate(s1e_rows):
        m = re.match(r"^(chr\w+):(\d+)-(\d+)$", str(r.get("SD98_v2.0", "")))
        if not m:
            continue
        rows.append(dict(rid=f"S{i}", kind="S1E", s1e_row=i, gene=r["Gene Name"], gene_id=r["Gene ID"],
                         chrom=m.group(1), start=int(m.group(2)), end=int(m.group(3))))
    names = {str(r["Gene Name"]) for r in s1e_rows} | {str(r["Gene Name"]) for r in s1c_rows}
    auto = {f"chr{i}" for i in range(1, 23)}
    bych = defaultdict(list)
    for r in rows:
        bych[r["chrom"]].append((r["start"], r["end"]))
    genes = []
    op = gzip.open if gff.endswith(".gz") else open
    with op(gff, "rt") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            t = line.rstrip("\n").split("\t")
            if len(t) < 9 or t[2] != "gene" or t[0] not in auto or "gene_biotype=protein_coding" not in t[8]:
                continue
            nm = re.search(r"(?:^|;)Name=([^;]+)", t[8]).group(1)
            genes.append((t[0], int(t[3]) - 1, int(t[4]), nm))
    def near(ch, s, e):
        return any(s < b + exclude_bp and e > a - exclude_bp for a, b in bych.get(ch, []))
    elig = sorted(g for g in genes if g[3] not in names and not near(g[0], g[1], g[2]))
    ctl = sorted(random.Random(seed).sample(elig, n_controls))
    for j, g in enumerate(ctl):
        rows.append(dict(rid=f"C{j}", kind="CTRL", s1e_row=-1, gene=g[3], gene_id="", chrom=g[0], start=g[1],
                         end=min(g[2], g[1] + control_bp)))
    return rows, len(elig)


REGION_COLS = ["rid", "kind", "s1e_row", "gene", "gene_id", "chrom", "start", "end"]


def write_regions(path, rows, extra=()):
    with open(path, "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(REGION_COLS + list(extra))
        for r in rows:
            w.writerow([r[c] for c in REGION_COLS] + [r[c] for c in extra])


def fetch(genome, chrom, start, end):
    out = subprocess.run(["samtools", "faidx", genome, f"{chrom}:{start + 1}-{end}"], capture_output=True,
                         text=True, check=True).stdout
    return "".join(out.split("\n")[1:])


def region_kmers(regions, genome):
    """Per region every K-mer fully inside it (ACGT only) -> Q (sorted unique canonical), idx (per position ->
    Q index, int32), mask (per position: touches a soft-masked base), offs (region offsets into idx / mask)."""
    allc, allm, offs = [], [], [0]
    for r in regions:
        seq = fetch(genome, r["chrom"], r["start"], r["end"])
        if len(seq) != r["end"] - r["start"]:
            sys.exit(f"{r['rid']}: fetched {len(seq)} bp for {r['end'] - r['start']}")
        can, msk = canonical_kmers(seq, with_mask=True)
        allc.append(can)
        allm.append(msk)
        offs.append(offs[-1] + len(can))
    can, msk = np.concatenate(allc), np.concatenate(allm)
    Q = np.unique(can)
    return Q, np.searchsorted(Q, can).astype(np.int32), msk, np.array(offs, dtype=np.int64)


# ---------------------------------------------------------------------------------------------------------
# meryl counting
# ---------------------------------------------------------------------------------------------------------

def write_kmer_fasta(codes, path):
    """One record per code (`>q` + the 30 bases), the query set for `meryl count`."""
    n = len(codes)
    rec = np.empty((n, K + 4), dtype=np.uint8)
    rec[:, 0], rec[:, 1], rec[:, 2], rec[:, K + 3] = ord(">"), ord("q"), ord("\n"), ord("\n")
    rec[:, 3:K + 3] = decode_codes(codes)
    rec.tofile(path)


def meryl_count_db(fasta, db, meryl=MERYL, threads=4, memory=16, log=None):
    cmd = [meryl, "count", f"k={K}", f"threads={threads}", f"memory={memory}", fasta, "output", db]
    with open(log or os.devnull, "w") as err:
        subprocess.run(cmd, stdout=err, stderr=err, check=True)


def meryl_lookup(genome_db, q_db, out_txt, meryl=MERYL, log=None):
    """`meryl print [intersect genome_db q_db]`: kmer<TAB>count for every query k-mer present in the genome."""
    with open(out_txt, "w") as out, open(log or os.devnull, "w") as err:
        subprocess.run([meryl, "print", "[", "intersect", genome_db, q_db, "]"], stdout=out, stderr=err, check=True)


def parse_meryl_print(path, Q):
    """kmer<TAB>count lines (meryl's own canonical form, which may be the reverse complement of ours) -> int32
    counts aligned to Q (0 where absent)."""
    counts = np.zeros(len(Q), dtype=np.int32)
    data = open(path, "rb").read()
    if not data:
        return counts
    lines = data.split(b"\n")
    if lines and lines[-1] == b"":
        lines.pop()
    chars = np.frombuffer(b"".join(l[:K] for l in lines), dtype=np.uint8).reshape(len(lines), K)
    codes = canonical(encode_strings(chars))
    vals = np.array([int(l[K + 1:]) for l in lines], dtype=np.int64)
    pos = np.searchsorted(Q, codes)
    pos[pos >= len(Q)] = 0
    hit = Q[pos] == codes
    if not hit.all():
        sys.exit(f"{int((~hit).sum())} printed k-mers are not in Q")
    counts[pos] = np.minimum(vals, np.iinfo(np.int32).max)
    return counts


def count_with_meryl(Q, genome, out, work, exclude=(), db=None, meryl=MERYL, threads=4, memory=16, keep_db=False):
    os.makedirs(work, exist_ok=True)
    fasta = genome
    if exclude:
        fasta = os.path.join(work, "genome.filtered.fa")
        names = [l.split("\t")[0] for l in open(genome + ".fai")]
        keep = [n for n in names if n not in exclude]
        with open(fasta, "w") as fh:
            subprocess.run(["samtools", "faidx", genome] + keep, stdout=fh, check=True)
    if db is None:
        db = os.path.join(work, "genome.meryl")
        meryl_count_db(fasta, db, meryl, threads, memory, os.path.join(work, "meryl_count.log"))
    qfa, qdb, hits = (os.path.join(work, x) for x in ("Q.fa", "Q.meryl", "hits.txt"))
    write_kmer_fasta(Q, qfa)
    meryl_count_db(qfa, qdb, meryl, threads, min(memory, 8), os.path.join(work, "meryl_q.log"))
    meryl_lookup(db, qdb, hits, meryl, os.path.join(work, "meryl_print.log"))
    counts = parse_meryl_print(hits, Q)
    counts.tofile(out)
    if not keep_db and db == os.path.join(work, "genome.meryl"):
        subprocess.run(["rm", "-rf", db])
    return counts


# ---------------------------------------------------------------------------------------------------------
# analysis
# ---------------------------------------------------------------------------------------------------------

def med(a):
    return float(np.median(a)) if len(a) else float("nan")


def region_stats(ii, mm, spec_k, noed_k, C, S, apes, min_k=MINK, p_floor=PFLOOR):
    """One region: ii = Q indices of its positions, mm = soft-mask flags, spec_k / noed_k = per-Q SPEC masks (with
    and without the edit-depth filter), C = {genome: counts}, S = {genome: scale}. Returns the per-genome parCN
    (median and mean statistic, NOED arm), SPEC presence, famCN analogue with its presence share p, the
    CHM13-private share of SPEC and the ape ratios r."""
    sp, no, fa = ii[spec_k[ii]], ii[noed_k[ii]], ii[~mm]
    row = dict(n_pos=len(ii), n_spec=len(sp), n_noed=len(no), n_fam=len(fa))
    for g, c in C.items():
        s = S[g]
        row[f"par_{g}"] = s * med(c[sp]) if len(sp) >= min_k else float("nan")
        row[f"parmean_{g}"] = s * float(c[sp].mean()) if len(sp) >= min_k else float("nan")
        row[f"parnoed_{g}"] = s * med(c[no]) if len(no) >= min_k else float("nan")
        row[f"specpres_{g}"] = float((c[sp] >= 1).mean()) if len(sp) else float("nan")
        if len(fa) >= min_k:
            cf = c[fa]
            pres = cf >= 1
            p = float(pres.mean())
            row[f"p_{g}"] = p
            row[f"fam_{g}"] = s * med(cf[pres]) if p >= p_floor else 0.0
        else:
            row[f"p_{g}"] = row[f"fam_{g}"] = float("nan")
    if len(sp):
        oth = np.zeros(len(sp), dtype=bool)
        for g in C:
            if g != "CHM13":
                oth |= C[g][sp] >= 1
        row["spec_chm13_private"] = float((~oth).mean())
    else:
        row["spec_chm13_private"] = float("nan")
    for a in apes:
        row[f"r_{a}"] = row[f"specpres_{a}"] / row[f"p_{a}"] if row.get(f"p_{a}") else float("nan")
    return row


def isnan(x):
    return x != x


def mean_of(flags):
    return float(np.mean(flags)) if len(flags) else float("nan")


def spearman(x, y):
    from scipy.stats import spearmanr
    return float(spearmanr(np.asarray(x, dtype=float), np.asarray(y, dtype=float)).statistic)


def margin_bin(m):
    if isnan(m):
        return "nan"
    for lo, hi in ((-1e9, 0.0), (0.0, 1.0), (1.0, 2.0), (2.0, 1e9)):
        if lo < m <= hi:
            return f"({lo}, {hi}]"
    return "nan"


def analyze(regions, idx, mask, offs, spec_k, noed_k, C, S, apes, s1e_rows, s1c_rows, min_k=MINK, p_floor=PFLOOR):
    """PREREG §2-§3 as written: every S1E row and control gets region_stats; S1E rows are joined to S1E (by data
    row) and S1C (first row of the gene with `In Table S1 = Yes`); clauses C0-C5, H1-H5, F-cal. Returns
    (s1e_result_rows, ctrl_rows, summary dict, named_rows)."""
    res = []
    for i, r in enumerate(regions):
        ii, mm = idx[offs[i]:offs[i + 1]], mask[offs[i]:offs[i + 1]]
        row = dict(r)
        row.update(region_stats(ii, mm, spec_k, noed_k, C, S, apes, min_k, p_floor))
        res.append(row)
    s1c_first = {}
    for r in s1c_rows:
        if r["In Table S1 (SD98 gene set)"] == "Yes" and r["Gene ID"] not in s1c_first:
            s1c_first[r["Gene ID"]] = r
    s, ctl = [], []
    for row in res:
        if row["kind"] != "S1E":
            ctl.append(row)
            continue
        e = s1e_rows[row["s1e_row"]]
        row["Family ID"], row["Family status"], row["parCN class"], row["Biotype"] = \
            e["Family ID"], e["Family status"], e["parCN class"], e["Biotype"]
        row["Median parCN"] = fnum(e["Median parCN"])
        c = s1c_first.get(row["gene_id"], {})
        row["Median famCN"] = fnum(c.get("Median famCN"))
        row["Max famCN Great Apes"] = fnum(c.get("Max famCN Great Apes"))
        row["famCN Status per Paralog"] = c.get("famCN Status per Paralog", "")
        row["Non-Syntenic with Chimpanzee"] = c.get("Non-Syntenic with Chimpanzee", "")
        s.append(row)
    R = {}
    # ---------- C0
    rs = [x for x in s if x["n_spec"] >= min_k]
    R["C0b_chm13_all2"] = bool(all(x["par_CHM13"] == 2 for x in rs))
    rc = [x for x in ctl if x["n_spec"] >= min_k]
    R["ctrl_resolved"] = len(rc)
    R["C0c_hg002_ctrl_eq2"] = mean_of([x["par_HG002"] == 2 for x in rc])
    fc = [x for x in ctl if x["n_fam"] >= min_k]
    R["ctrl_fam_resolved"] = len(fc)
    for a in apes + ["HG002"]:
        R[f"C0c_fam_{a}_ctrl_eq2"] = mean_of([x[f"fam_{a}"] == 2 for x in fc])
        R[f"ctrl_p_{a}_median"] = med([x[f"p_{a}"] for x in fc])
        R[f"ctrl_par_{a}_eq2"] = mean_of([x[f"par_{a}"] == 2 for x in rc])
    R["C0c_pass"] = R["C0c_hg002_ctrl_eq2"] >= 0.97 and all(R[f"C0c_fam_{a}_ctrl_eq2"] >= 0.90 for a in apes)
    # ---------- C1-C5
    for x in s:
        x["resolved"] = x["n_spec"] >= min_k
        x["agree"] = abs(x["par_HG002"] - x["Median parCN"]) <= 0.5
        x["agree_mean"] = abs(x["parmean_HG002"] - x["Median parCN"]) <= 0.5
        x["agree_noed"] = abs(x["parnoed_HG002"] - x["Median parCN"]) <= 0.5
        x["acro_parm"] = x["chrom"] in ACRO and x["start"] < 18_000_000
    R["n_s1e"] = len(s)
    R["resolved"] = sum(x["resolved"] for x in s)
    R["resolved_share"] = R["resolved"] / len(s)
    R["resolved_noed"] = sum(x["n_noed"] >= min_k for x in s)
    for cl in ("Fixed", "Nearly-Fixed", "Polymorphic"):
        sub = [x for x in s if x["resolved"] and x["parCN class"] == cl]
        R[f"{cl}_n"] = len(sub)
        R[f"{cl}_agree"] = mean_of([x["agree"] for x in sub])
        R[f"{cl}_agree_mean"] = mean_of([x["agree_mean"] for x in sub])
        R[f"{cl}_n_all"] = sum(x["parCN class"] == cl for x in s)
        sub2 = [x for x in sub if not x["acro_parm"]]
        R[f"{cl}_agree_noacro"] = mean_of([x["agree"] for x in sub2])
        R[f"{cl}_n_noacro"] = len(sub2)
        subn = [x for x in s if x["n_noed"] >= min_k and x["parCN class"] == cl]
        R[f"{cl}_agree_NOED"] = mean_of([x["agree_noed"] for x in subn])
        R[f"{cl}_n_NOED"] = len(subn)
        dist = defaultdict(int)
        for x in sub:
            dist[x["par_HG002"]] += 1
        R[f"{cl}_hg002_dist"] = {str(float(k)): v for k, v in sorted(dist.items())}
    rr = [x for x in s if x["resolved"]]
    R["C4_spearman"] = spearman([x["par_HG002"] for x in rr], [x["Median parCN"] for x in rr])
    R["C4_spearman_mean"] = spearman([x["parmean_HG002"] for x in rr], [x["Median parCN"] for x in rr])
    R["C1_pass"] = R["Fixed_agree"] >= 0.90
    R["C2_pass"] = R["Nearly-Fixed_agree"] >= 0.70
    R["C4_pass"] = R["C4_spearman"] >= 0.40

    def dis_reason(x):
        if x["spec_chm13_private"] > 0.5:
            return "CHM13-private"
        if x["acro_parm"]:
            return "acrocentric-p-arm"
        if x["n_spec"] < 300:
            return "SPEC<300"
        if not isnan(x["fam_HG002"]) and not isnan(x["fam_CHM13"]) and abs(x["fam_HG002"] - x["fam_CHM13"]) >= 1:
            return "HG002-CNV(famCN differs)"
        return "other"
    for x in s:
        x["c1_reason"] = dis_reason(x) if (x["resolved"] and not x["agree"]) else ""
    for cl in ("Fixed", "Nearly-Fixed", "Polymorphic"):
        cnt = defaultdict(int)
        for x in s:
            if x["resolved"] and not x["agree"] and x["parCN class"] == cl:
                cnt[x["c1_reason"]] += 1
        R[f"{cl}_disagree_reasons"] = dict(sorted(cnt.items(), key=lambda kv: -kv[1]))
    # ---------- H
    for x in s:
        x["ape_max"] = max(x[f"fam_{a}"] for a in apes) if not any(isnan(x[f"fam_{a}"]) for a in apes) else float("nan")
        x["evaluable"] = x["n_fam"] >= min_k
        x["ASM_H"] = bool(x["evaluable"] and x["fam_HG002"] > x["ape_max"])
        x["ASM_H_chm13"] = bool(x["evaluable"] and x["fam_CHM13"] > x["ape_max"])
        x["ASM_class"] = "NA" if not x["evaluable"] else ("dup" if x["ape_max"] < 2.5 else "exp") if x["ASM_H"] else "no"
        x["soto_margin"] = x["Median famCN"] - x["Max famCN Great Apes"]
    stl = lambda x: x["famCN Status per Paralog"]
    for lab, key in (("Duplicated in humans", "H1"), ("Expanded in humans", "H2"), ("Expanded in great apes", "H3"),
                     ("Undetermined", "U")):
        sub = [x for x in s if stl(x) == lab and x["evaluable"]]
        R[f"{key}_n"] = len(sub)
        R[f"{key}_n_all"] = sum(stl(x) == lab for x in s)
        R[f"{key}_asmH"] = mean_of([x["ASM_H"] for x in sub])
        R[f"{key}_asmH_chm13"] = mean_of([x["ASM_H_chm13"] for x in sub])
        cnt = defaultdict(int)
        for x in sub:
            cnt[x["ASM_class"]] += 1
        R[f"{key}_asm_class"] = dict(sorted(cnt.items(), key=lambda kv: -kv[1]))
    R["H1_pass"], R["H2_pass"], R["H3_pass"] = R["H1_asmH"] >= 0.80, R["H2_asmH"] >= 0.60, R["H3_asmH"] <= 0.25
    sub = [x for x in s if stl(x) == "Duplicated in humans" and x["evaluable"]]
    R["H1_apemax_lt2.5"] = mean_of([x["ape_max"] < 2.5 for x in sub])

    def h_reason(x):
        rs_ = []
        if not isnan(x["soto_margin"]) and x["soto_margin"] < 1:
            rs_.append("Soto-margin<1")
        if not isnan(x["Max famCN Great Apes"]) and abs(x["ape_max"] - x["Max famCN Great Apes"]) >= 1:
            rs_.append("ape-max-differs(ours %s Soto %.1f)" % (x["ape_max"], x["Max famCN Great Apes"]))
        if x["n_fam"] < 300:
            rs_.append("FAM<300")
        if abs(x["fam_HG002"] - x["fam_CHM13"]) >= 1:
            rs_.append("HG002!=CHM13x2")
        if any(x[f"p_{a}"] < p_floor for a in apes):
            rs_.append("ape-p-floor")
        return ";".join(rs_) if rs_ else "other"
    for x in s:
        flag = x["evaluable"] and ((stl(x) in ("Duplicated in humans", "Expanded in humans") and not x["ASM_H"])
                                   or (stl(x) == "Expanded in great apes" and x["ASM_H"]))
        x["h_reason"] = h_reason(x) if flag else ""
        x["h_reason_primary"] = re.sub(r"\(.*", "", x["h_reason"].split(";")[0]) if flag else ""
        x["_flag"] = flag
    for lab, key in (("Duplicated in humans", "H1"), ("Expanded in humans", "H2"), ("Expanded in great apes", "H3")):
        cnt = defaultdict(int)
        for x in s:
            if x["_flag"] and stl(x) == lab:
                cnt[x["h_reason_primary"]] += 1
        R[f"{key}_reasons_primary"] = dict(sorted(cnt.items(), key=lambda kv: -kv[1]))
        sub = [x for x in s if stl(x) == lab and x["evaluable"]]
        by = defaultdict(list)
        for x in sub:
            by[margin_bin(x["soto_margin"])].append(x["ASM_H"])
        R[f"{key}_asmH_by_margin"] = {b: [len(v), round(float(np.mean(v)), 3) if v else None] for b, v in sorted(by.items())}
    R["H3_soto_h_gt_a"] = sum(stl(x) == "Expanded in great apes" and x["soto_margin"] > 0 for x in s)
    fam = {}
    for x in s:
        f = x["Family ID"]
        if not f:
            continue
        d = fam.setdefault(f, dict(status=x["Family status"], sotoH=False, allU=True, asmH=False, asmH_c=False,
                                   anyeval=False))
        d["sotoH"] |= stl(x) in ("Duplicated in humans", "Expanded in humans")
        d["allU"] &= stl(x) == "Undetermined"
        d["asmH"] |= x["ASM_H"]
        d["asmH_c"] |= x["ASM_H_chm13"]
        d["anyeval"] |= x["evaluable"]
    f1 = [d for d in fam.values() if d["sotoH"] and d["anyeval"]]
    R["H4_n"], R["H4_n_all"] = len(f1), sum(d["sotoH"] for d in fam.values())
    R["H4_share"], R["H4_share_chm13"] = mean_of([d["asmH"] for d in f1]), mean_of([d["asmH_c"] for d in f1])
    R["H4_pass"] = R["H4_share"] >= 0.80
    f5 = [d for d in fam.values() if d["status"] == "Human duplicated gene family" and not d["sotoH"] and d["anyeval"]]
    R["H5_syntenyonly_n"], R["H5_syntenyonly_asmH"] = len(f5), mean_of([d["asmH"] for d in f5])
    f5a, f5b = [d for d in f5 if not d["allU"]], [d for d in f5 if d["allU"]]
    R["H5_syntenyonly_notallU_n"] = len(f5a)
    R["H5_syntenyonly_notallU_asmH"] = mean_of([d["asmH"] for d in f5a]) if f5a else None
    R["H5_syntenyonly_allU_n"] = len(f5b)
    R["H5_syntenyonly_allU_asmH"] = mean_of([d["asmH"] for d in f5b]) if f5b else None
    f6 = [d for d in fam.values() if d["status"] == "Undetermined" and d["anyeval"]]
    R["H5_undetermined_fam_n"], R["H5_undetermined_fam_asmH"] = len(f6), mean_of([d["asmH"] for d in f6])
    f7 = [d for d in f6 if not d["allU"]]
    R["H5_undetermined_fam_notallU_n"], R["H5_undetermined_fam_notallU_asmH"] = len(f7), mean_of([d["asmH"] for d in f7])
    fe = [x for x in s if stl(x) in ("Duplicated in humans", "Expanded in humans", "Expanded in great apes") and x["evaluable"]]
    R["Fcal_human_spearman"] = spearman([x["fam_HG002"] for x in fe], [x["Median famCN"] for x in fe])
    R["Fcal_human_spearman_chm13"] = spearman([x["fam_CHM13"] for x in fe], [x["Median famCN"] for x in fe])
    R["Fcal_ape_spearman"] = spearman([x["ape_max"] for x in fe], [x["Max famCN Great Apes"] for x in fe])
    R["Fcal_human_pass"], R["Fcal_ape_pass"] = R["Fcal_human_spearman"] >= 0.60, R["Fcal_ape_spearman"] >= 0.50
    ev = [x for x in s if x["evaluable"]]
    for a in apes:
        R[f"fam_p_median_{a}"] = med([x[f"p_{a}"] for x in ev])
        R[f"fam_pfloor_{a}"] = sum(x[f"p_{a}"] < p_floor for x in ev)
    named = []
    for k, ids in NAMED.items():
        for x in s:
            if x["Family ID"] in ids:
                named.append(dict(named=k, **x))
    for x in s:
        x.pop("_flag", None)
    return s, ctl, R, named


def write_rows(path, rows, cols=None):
    if not rows:
        open(path, "w").close()
        return
    cols = cols or [c for c in rows[0] if not c.startswith("_")]
    with open(path, "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(cols)
        for r in rows:
            w.writerow(["" if (v is None or (isinstance(v, float) and isnan(v))) else v for v in (r.get(c, "") for c in cols)])


# ---------------------------------------------------------------------------------------------------------
# subcommands
# ---------------------------------------------------------------------------------------------------------

def cmd_regions(a):
    rows, n_elig = regions_table(read_tsv(a.s1e), read_tsv(a.s1c), a.gff, a.n_controls, a.seed, a.control_bp,
                                 a.exclude_bp)
    write_regions(a.out, rows)
    print(f"[regions] {sum(r['kind'] == 'S1E' for r in rows)} S1E regions + {a.n_controls} controls (of {n_elig} "
          f"eligible) -> {a.out}", file=sys.stderr)


def cmd_kmers(a):
    regions = read_tsv(a.regions)
    for r in regions:
        r["start"], r["end"], r["s1e_row"] = int(r["start"]), int(r["end"]), int(r["s1e_row"])
    Q, idx, mask, offs = region_kmers(regions, a.genome)
    os.makedirs(a.out_dir, exist_ok=True)
    Q.tofile(os.path.join(a.out_dir, "Q.u64"))
    np.savez(os.path.join(a.out_dir, "pos.npz"), idx=idx, mask=mask, offs=offs)
    for i, r in enumerate(regions):
        r["n_pos"] = int(offs[i + 1] - offs[i])
    write_regions(os.path.join(a.out_dir, "regions.tsv"), regions, ["n_pos"])
    print(f"[kmers] {len(idx)} positions, {len(Q)} distinct canonical {K}-mers, masked share "
          f"{float(mask.mean()):.4f} -> {a.out_dir}/Q.u64, pos.npz, regions.tsv", file=sys.stderr)


def cmd_count(a):
    Q = np.fromfile(a.q, dtype=np.uint64)
    exclude = set(a.exclude.split(",")) if a.exclude else set()
    if a.meryl:
        counts = count_with_meryl(Q, a.genome, a.out, a.work or (a.out + "_work"), exclude, a.db, a.meryl,
                                  a.threads, a.memory, a.keep_db)
    else:
        counts = stream_kmer_counts(a.genome, Q, exclude).astype(np.int32)
        counts.tofile(a.out)
    print(f"[count] {a.genome}: {int((counts > 0).sum())} of {len(Q)} query k-mers present, "
          f"{int(counts.sum())} hits -> {a.out}", file=sys.stderr)


def cmd_edit_depth(a):
    Q = np.fromfile(a.q, dtype=np.uint64)
    c = np.fromfile(a.counts, dtype=np.int32)
    cands = Q[(c == 1) & (Q != 0)]
    exclude = set(a.exclude.split(",")) if a.exclude else set()
    if a.cands_out:
        cands.tofile(a.cands_out)
    out = []
    for i in range(0, len(cands), a.batch):
        out.append(edit_depth_batch(cands[i:i + a.batch], a.genome, exclude))
        print(f"  {min(i + a.batch, len(cands))}/{len(cands)} candidates", file=sys.stderr, flush=True)
    ed = np.concatenate(out) if out else np.zeros(0, dtype=np.uint32)
    ed.tofile(a.out)
    print(f"[edit-depth] {len(cands)} candidates; {int((ed >= 100).sum())} dropped at edit depth >= 100 -> {a.out}",
          file=sys.stderr)


def cmd_analyze(a):
    Q = np.fromfile(a.q, dtype=np.uint64)
    pz = np.load(a.pos)
    idx, mask, offs = pz["idx"], pz["mask"], pz["offs"]
    regions = read_tsv(a.regions)
    for r in regions:
        r["start"], r["end"], r["s1e_row"] = int(r["start"]), int(r["end"]), int(r["s1e_row"])
    C, S = {}, dict(SCALE)
    for spec in a.counts:
        name, _, files = spec.partition("=")
        C[name] = sum(np.fromfile(f, dtype=np.int32).astype(np.int64) for f in files.split("+"))
        if len(C[name]) != len(Q):
            sys.exit(f"{name}: {len(C[name])} counts for {len(Q)} k-mers")
    for spec in a.scale or ():
        name, _, v = spec.partition("=")
        S[name] = int(v)
    for g in C:
        if g not in S:
            sys.exit(f"no scale for {g} (--scale {g}=1|2)")
    if "CHM13" not in C or "HG002" not in C:
        sys.exit("--counts must include CHM13=... and HG002=...")
    apes = a.apes.split(",")
    candmask = (C["CHM13"] == 1) & (Q != 0)
    ed = np.full(len(Q), 1 << 30, dtype=np.int64)
    edv = np.concatenate([np.fromfile(f, dtype=np.uint32) for f in a.edit_depth.split(",")]).astype(np.int64)
    if len(edv) != int(candmask.sum()):
        sys.exit(f"edit depth has {len(edv)} values for {int(candmask.sum())} candidates")
    ed[candmask] = edv
    spec_k, noed_k = candmask & (ed <= 99), candmask
    s, ctl, R, named = analyze(regions, idx, mask, offs, spec_k, noed_k, C, S, apes, read_tsv(a.s1e), read_tsv(a.s1c),
                               a.min_k, a.p_floor)
    write_rows(a.out_prefix + "results_s1e.tsv", s)
    write_rows(a.out_prefix + "results_ctrl.tsv", ctl)
    write_rows(a.out_prefix + "results_named.tsv", named)
    json.dump(R, open(a.out_prefix + "summary.json", "w"), indent=1, default=str)
    print(json.dumps(R, indent=1, default=str))
    print(f"[analyze] C1 Fixed {R['Fixed_agree']:.3f} ({R['Fixed_n']}) C2 Nearly-Fixed {R['Nearly-Fixed_agree']:.3f} "
          f"({R['Nearly-Fixed_n']}) C3 {R['Polymorphic_agree']:.3f} C4 {R['C4_spearman']:.3f} | H1 {R['H1_asmH']:.3f} "
          f"H2 {R['H2_asmH']:.3f} H3 {R['H3_asmH']:.3f} H4 {R['H4_share']:.3f} -> {a.out_prefix}summary.json",
          file=sys.stderr)


def main(argv=None):
    ap = argparse.ArgumentParser(prog="parcn_assembly.py", description=__doc__.split("\n\n")[0],
                                 epilog="The module docstring has the rule, the numbers and the file layout.")
    sub = ap.add_subparsers(dest="cmd", required=True, metavar="SUBCOMMAND")

    p = sub.add_parser("regions", help="S1E SD98_v2.0 regions + 300 RefSeq controls -> regions.tsv")
    p.add_argument("--s1e", default=S1E)
    p.add_argument("--s1c", default=S1C)
    p.add_argument("--gff", required=True, help="Reference/chm13v2.0_RefSeq_full.gff.gz")
    p.add_argument("--n-controls", type=int, default=300)
    p.add_argument("--seed", type=int, default=20260929)
    p.add_argument("--control-bp", type=int, default=20000)
    p.add_argument("--exclude-bp", type=int, default=1_000_000)
    p.add_argument("--out", required=True)
    p.set_defaults(func=cmd_regions)

    p = sub.add_parser("kmers", help="every 30-mer of every region -> Q.u64, pos.npz, regions.tsv (+ n_pos)")
    p.add_argument("--regions", required=True)
    p.add_argument("--genome", required=True, help="chm13v2.0.fa (soft-masked, samtools-indexed)")
    p.add_argument("--out-dir", required=True)
    p.set_defaults(func=cmd_kmers)

    p = sub.add_parser("count", help="exact counts of the Q k-mers in one assembly -> .i32 (meryl or numpy)")
    p.add_argument("--q", required=True, help="Q.u64")
    p.add_argument("--genome", required=True, help="assembly FASTA (samtools-indexed when --exclude is used)")
    p.add_argument("--exclude", help="comma-separated sequence names to skip (e.g. chrY for the reference)")
    p.add_argument("--meryl", nargs="?", const=MERYL, help=f"count with meryl (default binary {MERYL}); omit for "
                                                            "the numpy streaming counter (slow, no DB)")
    p.add_argument("--db", help="an existing meryl DB of --genome (skips `meryl count`)")
    p.add_argument("--threads", type=int, default=4)
    p.add_argument("--memory", type=int, default=16, help="meryl memory= in GB")
    p.add_argument("--work", help="meryl work directory (default <out>_work)")
    p.add_argument("--keep-db", action="store_true", help="keep the genome meryl DB under --work")
    p.add_argument("--out", required=True)
    p.set_defaults(func=cmd_count)

    p = sub.add_parser("edit-depth", help="QuicK-mer2 -e 1 edit depth of the once-only Q k-mers -> .u32")
    p.add_argument("--q", required=True)
    p.add_argument("--counts", required=True, help="the reference's .i32 (candidates = count 1, code != 0)")
    p.add_argument("--genome", required=True, help="the reference FASTA")
    p.add_argument("--exclude", help="sequence names to skip (chrY)")
    p.add_argument("--batch", type=int, default=2_400_000)
    p.add_argument("--cands-out", help="also write the candidate codes (.u64)")
    p.add_argument("--out", required=True)
    p.set_defaults(func=cmd_edit_depth)

    p = sub.add_parser("analyze", help="parCN / famCN per genome, calls, clauses -> results_*.tsv, summary.json")
    p.add_argument("--regions", required=True, help="regions.tsv (from `kmers`)")
    p.add_argument("--q", required=True)
    p.add_argument("--pos", required=True)
    p.add_argument("--counts", required=True, nargs="+", metavar="NAME=FILE[+FILE]",
                   help="per-genome .i32 tables; CHM13 and HG002 required; a `+` sums haplotypes (GGO=mat+pat)")
    p.add_argument("--scale", nargs="*", metavar="NAME=1|2", help=f"override the default scales {SCALE}")
    p.add_argument("--apes", default="PTR,GGO,PPY", help="ape genome names among --counts (default %(default)s)")
    p.add_argument("--edit-depth", required=True, help="comma-separated .u32 files, in candidate order")
    p.add_argument("--s1e", default=S1E)
    p.add_argument("--s1c", default=S1C)
    p.add_argument("--min-k", type=int, default=MINK)
    p.add_argument("--p-floor", type=float, default=PFLOOR)
    p.add_argument("--out-prefix", required=True)
    p.set_defaults(func=cmd_analyze)

    a = ap.parse_args(argv)
    a.func(a)


if __name__ == "__main__":
    main()
