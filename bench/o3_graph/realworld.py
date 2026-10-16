#!/usr/bin/env python3
"""Simulation v2 (prereg docs/PREREG_o3_graph_real_2026-10-15.md): worlds from HANDPICKED REAL families (families.py) on the
real CHM13 v2.0 genome, every read simulated from a curated CAT isoform of a real copy.

    realworld.py mask OUT_DIR FAMILIES.json --set S --seed N     OUT_DIR/mask.bed + masked.fa: the REMOVED arm's genome
    realworld.py build OUT FAMILIES.json --set S --scenario SC --seed N [--genome FASTA] [--per-copy 60]

Per family one SPECIAL copy (seeded; its main isoform has >= 3 exons when any copy's does) plays the scenario's part:
  INTACT    every copy in the genome and the catalog
  REMOVED   the special copy hard-masked in the genome (`mask`) and dropped from the catalog; its reads still exist
  OUTSIDE   the special copy stays in the genome but leaves the catalog (an uncatalogued paralog)
  BRANCH    an extra ABSENT copy `<family>.B`: the special copy's isoforms with a novel 150-bp exon in the middle of its main
            isoform's middle intron (spliced into every isoform with an intron across that point) and 0.2% exonic SNPs
  SILENT    as INTACT, but a seeded third of the copies (at least one) express nothing; every copy stays in the catalog
            (round 7: what O1's RNA-supported seeding drops, and whether it drops an expressed copy)
  CROSS     as INTACT, and EVERY copy also expresses up to MAX_PROJ isoform structures annotated only on a paralog (review C1
            of register 1298: the catalog is not copy-specific): each paralog isoform's sense sequence spliced-aligned
            (minimap2 splice:hq -uf) to the copy's locus, kept when its primary alignment covers >= PROJ_COVER of it with the
            same exon count, canonical junctions and an intron chain none of the copy's own isoforms has; the copy's own
            genomic sequence over those exons is the transcript (truth isoform `proj:<paralog copy>|<isoform>`, weight
            OTHER; the copy's own-isoform reads are INTACT's). The catalog lists each copy's own isoforms only.
            `--projection liftover`: the structures are lifted through the two loci's genomic alignment instead (asm20),
            independent of the spliced projection union O2 itself uses.
Reads: per expressed copy PER_COPY reads over its isoforms (main 2.0, others 0.5), each the isoform's sense sequence with 0-30
bp end jitter and HiFi errors (bench/sim.py, 0.1% substitutions, 0.03% indels), from the UNMASKED genome. `--het RATE` (round 6):
a diploid donor, each copy with heterozygous sites at RATE of its exonic bases; a read comes from haplotype 1 or 2 (truth `hap`). Names are a
seeded shuffle `s0000001..` (they carry no truth); reads.truth.E0.tsv holds it, with `novel` = 1 for a read that
carries the novel exon (the error-free sequence before errors contains it whole).
Products (the layout detect.py / score.py read): genome.ref.fa (+ .fai; a link to the world's genome), reads.E0.fq,
reads.truth.E0.tsv, P.fam.copies.tsv / .fa (catalog; representative = the main isoform), isoforms.json, kinds.json,
copy_family.json, world.json, copies.tsv, regions.txt (each family's span per chromosome +- 20 kb, merged).
"""
import argparse
import collections
import importlib.util
import json
import math
import os
import random
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
BENCH = os.path.dirname(HERE)
sys.path.insert(0, BENCH)
sys.path.insert(0, os.path.join(BENCH, "o3_fewread"))

GENOME = "/mnt/linuxdisk/home/juanfraitu/npip_ladder/idx/target.fa"      # = chm13v2.0.fa (same contigs and sequence)
SCENARIOS = ("INTACT", "REMOVED", "OUTSIDE", "BRANCH", "CROSS", "SILENT")
MAIN_WEIGHT, OTHER_WEIGHT, NOVEL_LEN, BRANCH_SNP, REGION_PAD = 2.0, 0.5, 150, 0.002, 20000
MIN_READ = 200
MAX_PROJ, PROJ_COVER, PROJ_PAD, JUNCTION_TOL = 2, 0.95, 1000, 10
COMP = str.maketrans("ACGTNacgtn", "TGCANtgcan")


def splice(genome, exons, strand, snps, insert_at=None, novel=""):
    """An isoform's sequence in its sense orientation: genomic exons (0-based half-open; `genome` a string or fetch(a, b)),
    SNPs at genomic positions {pos: base}, and `novel` spliced in when one of the isoform's introns spans genomic position
    `insert_at` (a real new exon there joins every isoform that splices across it)."""
    get = (lambda a, b: genome[a:b]) if isinstance(genome, str) else genome
    ex = sorted(exons)
    parts = []
    for i, (a, b) in enumerate(ex):
        s = list(get(a, b).upper())
        for p, base in snps.items():
            if a <= p < b:
                s[p - a] = base
        parts.append("".join(s))
        if insert_at is not None and i + 1 < len(ex) and b <= insert_at < ex[i + 1][0]:
            parts.append(novel)
    seq = "".join(parts)
    return seq[::-1].translate(COMP) if strand == "-" else seq


def exons_from_cigar(tstart, cigar):
    """Reference blocks of a PAF alignment (`cg:Z:` CIGAR on the target's forward strand): a new exon at every N (intron);
    deletions stay inside their exon, insertions take no reference."""
    import re
    ex, pos, start = [], tstart, tstart
    for n, op in re.findall(r"(\d+)([MIDNSHP=X])", cigar):
        n = int(n)
        if op == "N":
            ex.append((start, pos))
            pos += n
            start = pos
        elif op in "MD=X":
            pos += n
    ex.append((start, pos))
    return ex


def canonical_junctions(genome, exons, strand):
    """Every intron GT-AG / GC-AG / AT-AC on the transcript strand (`genome`: a string or fetch(a, b), forward strand)."""
    get = (lambda a, b: genome[a:b]) if isinstance(genome, str) else genome
    ok = {("GT", "AG"), ("GC", "AG"), ("AT", "AC")} if strand == "+" else {("CT", "AC"), ("CT", "GC"), ("GT", "AT")}
    ex = sorted(exons)
    return all((get(a, a + 2).upper(), get(b - 2, b).upper()) in ok for (_, a), (b, _) in zip(ex, ex[1:]))


def project(iso_seq, locus, lo, strand, n_exons, cover=PROJ_COVER):
    """A paralog isoform's sense sequence onto a copy's locus (`locus` = genome[lo:...], forward strand): its exons in genome
    coordinates when the primary spliced alignment covers >= `cover` of it with `n_exons` exons and canonical junctions on
    `strand`, else None."""
    return project_many({"q": (iso_seq, n_exons)}, locus, lo, strand, cover)["q"]


def project_many(queries, locus, lo, strand, cover=PROJ_COVER):
    """project() for {name: (sense sequence, exon count)} in one minimap2 run -> {name: exons or None}."""
    import subprocess
    import tempfile
    with tempfile.TemporaryDirectory() as d:
        open(f"{d}/t.fa", "w").write(f">locus\n{locus}\n")
        with open(f"{d}/q.fa", "w") as fh:
            for name, (seq, _) in queries.items():
                fh.write(f">{name}\n{seq}\n")
        out = subprocess.run(["minimap2", "-c", "-x", "splice:hq", "-uf", "--secondary=no", f"{d}/t.fa", f"{d}/q.fa"],
                             capture_output=True, text=True, check=True).stdout
    res = {name: None for name in queries}
    done = set()
    for ln in out.splitlines():
        f = ln.split("\t")
        tags = {x[:5]: x[5:] for x in f[12:]}
        q = f[0]
        if q in done or tags.get("tp:A:") != "P" or (int(f[3]) - int(f[2])) < cover * int(f[1]):
            continue
        done.add(q)                                          # the first full-length primary decides
        if (f[4] == "+") != (strand == "+"):                 # -uf: the sense read lands on the copy's own strand
            continue
        ex = exons_from_cigar(int(f[7]), tags["cg:Z:"])
        if len(ex) == queries[q][1] and canonical_junctions(locus, ex, strand):
            res[q] = [(a + lo, b + lo) for a, b in ex]
    return res


def lift_positions(cigar, qs, ts, positions):
    """{query position: target position, or None outside the alignment} through a '+' PAF alignment's CIGAR (query from qs,
    target from ts). A query base inside an insertion lifts to the next target base."""
    import re
    want = sorted(set(positions))
    out, q, t, k = {}, qs, ts, 0
    for n, op in re.findall(r"(\d+)([MIDNSHP=X])", cigar):
        n = int(n)
        if op in "M=X":
            while k < len(want) and want[k] < q + n:
                out[want[k]] = t + (want[k] - q) if want[k] >= q else None
                k += 1
            q += n
            t += n
        elif op == "I":
            while k < len(want) and want[k] < q + n:
                out[want[k]] = t if want[k] >= q else None
                k += 1
            q += n
        elif op in "DN":
            t += n
    for x in want[k:]:
        out[x] = None
    return out


def liftover(p_locus, structures, c_locus):
    """Paralog isoform structures (exons on p_locus, both loci in their transcripts' orientation) lifted onto c_locus through
    the loci's genomic alignment (minimap2 -x asm20, the longest '+' primary): a projection independent of the spliced
    transcript alignment `project` uses. -> per structure its exons on c_locus, or None when an exon end falls outside the
    alignment or the exons do not stay ordered."""
    import subprocess
    import tempfile
    with tempfile.TemporaryDirectory() as d:
        open(f"{d}/t.fa", "w").write(f">c\n{c_locus}\n")
        open(f"{d}/q.fa", "w").write(f">p\n{p_locus}\n")
        paf = subprocess.run(["minimap2", "-c", "-x", "asm20", "--secondary=no", f"{d}/t.fa", f"{d}/q.fa"],
                             capture_output=True, text=True, check=True).stdout
    best = None
    for ln in paf.splitlines():
        f = ln.split("\t")
        tags = {x[:5]: x[5:] for x in f[12:]}
        if f[4] == "+" and tags.get("tp:A:") == "P" and (best is None or int(f[3]) - int(f[2]) > best[1] - best[0]):
            best = (int(f[2]), int(f[3]), int(f[7]), tags["cg:Z:"])
    out = []
    for ex in structures:
        if best is None:
            out.append(None)
            continue
        pos = [x for a, b in ex for x in (a, b - 1)]
        lift = lift_positions(best[3], best[0], best[2], pos)
        new = [(lift[a], lift[b - 1] + 1) if lift[a] is not None and lift[b - 1] is not None else None for a, b in ex]
        ok = all(e is not None and e[0] < e[1] for e in new) and all(x[1] <= y[0] for x, y in zip(new, new[1:]))
        out.append(new if ok else None)
    return out


def oriented(fetch, lo, hi, strand):
    """A locus genome[lo:hi] in its transcript's orientation, and a function taking its exons back to genome coordinates."""
    seq = fetch(lo, hi).upper()
    if strand == "-":
        return seq[::-1].translate(COMP), lambda ex: sorted((hi - b, hi - a) for a, b in ex)
    return seq, lambda ex: sorted((lo + a, lo + b) for a, b in ex)


def to_oriented(exons, lo, hi, strand):
    """Genomic exons -> coordinates on the oriented locus genome[lo:hi]."""
    return sorted((hi - b, hi - a) for a, b in exons) if strand == "-" else sorted((a - lo, b - lo) for a, b in exons)


def cross_structures(own, projected, max_proj=MAX_PROJ, rng=None):
    """{name: exons} of the projected paralog structures a copy expresses: intron chains new to its own isoforms (`own`; a
    chain whose every junction is within JUNCTION_TOL bp of a known chain's is that chain, a projection artifact), one name
    per new chain (the first in name order), at most max_proj of them (a seeded sample)."""
    def chain(ex):
        ex = sorted(ex)
        return tuple((a[1], b[0]) for a, b in zip(ex, ex[1:]))

    def known(c, seen):
        return any(len(k) == len(c) and all(abs(x - y) <= JUNCTION_TOL for i, j in zip(c, k) for x, y in zip(i, j))
                   for k in seen)
    seen, pick = [chain(v) for v in own.values()], {}
    for name in sorted(projected):
        ex = projected[name]
        if ex is None or not chain(ex) or known(chain(ex), seen):
            continue
        seen.append(chain(ex))
        pick[name] = ex
    names = sorted(pick)
    if max_proj is not None and len(names) > max_proj:
        names = sorted((rng or random.Random(0)).sample(names, max_proj))
    return {n: pick[n] for n in names}


def hp_runs(seq, min_run=3):
    """[(start, end)] of the runs of one base at least min_run long."""
    out, i = [], 0
    while i < len(seq):
        j = i
        while j < len(seq) and seq[j] == seq[i]:
            j += 1
        if j - i >= min_run:
            out.append((i, j))
        i = j
    return out


def _poisson(lam, rng):
    lim, k, p = math.exp(-lam), 0, 1.0
    while True:
        p *= rng.random()
        if p <= lim:
            return k
        k += 1


def _length(dist, rng):
    r, acc = rng.random(), 0.0
    for k, v in sorted(dist.items(), key=lambda kv: int(kv[0])):
        acc += v
        if r < acc:
            return int(k)
    return max(int(k) for k in dist)


def draw_edits(seq, rng, prof):
    """One read's sequencing errors at a calibrated profile's per-base rates (calibrate.py profile): substitutions;
    insertions and deletions, each inside a homopolymer run (>= 3) with the profile's homopolymer fraction, lengths drawn from
    its length distribution. [(position, 'S' | 'I' | 'D', base | bases | length)]."""
    n, runs, edits = len(seq), hp_runs(seq), []
    for _ in range(_poisson(prof["sub_error"] * n, rng)):
        p = rng.randrange(n)
        edits.append((p, "S", rng.choice([b for b in "ACGT" if b != seq[p]])))
    for kind, key in (("I", "ins"), ("D", "del")):
        for _ in range(_poisson(prof[f"{key}_event"] * n, rng)):
            k = _length(prof[f"{key}_len"], rng)
            if runs and rng.random() < prof[f"{key}_hp_frac"]:
                s0, e0 = rng.choice(runs)
                if kind == "I":
                    edits.append((rng.randrange(s0, e0), "I", seq[s0] * k))
                else:
                    k = min(k, e0 - s0 - 1)
                    edits.append((rng.randrange(s0, e0 - k + 1), "D", k))
            else:
                p = rng.randrange(n)
                edits.append((p, "I", "".join(rng.choice("ACGT") for _ in range(k))) if kind == "I" else (p, "D", k))
    return edits


def apply_edits(seq, edits):
    """The read after its edits, each at its position in the ORIGINAL sequence (applied from the 3' end)."""
    s = list(seq)
    for p, kind, x in sorted(edits, key=lambda e: (e[0], "IDS".index(e[1])), reverse=True):
        if kind == "S":
            if p < len(s):
                s[p] = x
        elif kind == "D":
            del s[p:p + x]
        else:
            s[p:p] = list(x)
    return "".join(s)


def draw_truncation(seq, pairs, rng, min_len=MIN_READ, tries=10):
    """A read's transcript after a calibrated (5', 3') truncation drawn from `pairs`, redrawn (up to `tries`) while it would leave
    < min_len; the whole transcript only when no draw fits (review of round 9: a single draw fell back to full length for 43-61%
    of FAM25's reads)."""
    for _ in range(tries):
        m5, m3 = rng.choice(pairs)
        out = truncate(seq, m5, m3, min_len)
        if len(out) < len(seq) or (m5 == 0 and m3 == 0) or int(m5 * len(seq)) + int(m3 * len(seq)) == 0:
            return out
    return seq


def truncate(seq, m5, m3, min_len=MIN_READ):
    """A transcript that lost fractions m5 of its 5' end and m3 of its 3' end; the whole transcript when < min_len would be left."""
    n5, n3 = int(m5 * len(seq)), int(m3 * len(seq))
    return seq[n5:len(seq) - n3] if len(seq) - n5 - n3 >= min_len else seq


def retain(exons, keep):
    """Exons (genome order) with the introns indexed in `keep` retained: each merges its two flanking exons."""
    out = [list(exons[0])]
    for i, (s0, e0) in enumerate(exons[1:]):
        if i in keep:
            out[-1][1] = e0
        else:
            out.append([s0, e0])
    return [tuple(x) for x in out]


def expr_mult(seed, cid, sigma):
    """A copy's expression multiplier: lognormal with mean 1 and log-sd sigma, deterministic per (seed, copy)."""
    if sigma <= 0:
        return 1.0
    return math.exp(sigma * random.Random(f"{seed}:{cid}:expr").gauss(0, 1) - sigma * sigma / 2)


def het_sites(positions, rate, rng, base_at):
    """{genomic position: alt base}: a copy's heterozygous sites (round 6: a diploid donor), `rate` of its exonic positions,
    each alt base different from the reference base (base_at(p))."""
    pos = sorted(set(positions))
    k = int(round(rate * len(pos)))
    if not k:
        return {}
    out = {}
    for p in sorted(rng.sample(pos, k)):
        ref = base_at(p).upper()
        out[p] = rng.choice([b for b in "ACGT" if b != ref])
    return out


def carries(novel, seq):
    """Does a read's error-free sequence hold the whole novel exon (it is reverse-complemented on a minus-strand copy)?"""
    return bool(novel) and (novel in seq or novel[::-1].translate(COMP) in seq)


def anon_names(rows, seed):
    """Names for `rows` (one per read): a seeded permutation of s0000001.. so a name carries nothing of its source."""
    order = list(range(1, len(rows) + 1))
    random.Random(f"anon:{seed}").shuffle(order)
    return [f"s{k:07d}" for k in order]


def special_copy(fam, seed):
    """The family's scenario copy: seeded among copies whose main isoform has >= 3 exons (else >= 2, else any)."""
    for need in (3, 2, 1):
        ok = [c for c in fam["copies"] if len(c["transcripts"][c["main"]]) >= need]
        if ok:
            return random.Random(f"{seed}:{fam['family']}").choice(ok)


def mid_junction(exons):
    """The genomic position of the novel exon: the middle of the main isoform's middle intron."""
    ex = sorted(exons)
    if len(ex) < 2:
        return None
    k = (len(ex) - 1) // 2
    return (ex[k][1] + ex[k + 1][0]) // 2


def cmd_mask(a):
    fj = json.load(open(a.families))
    os.makedirs(a.out_dir, exist_ok=True)
    with open(f"{a.out_dir}/mask.bed", "w") as fh:
        for fam in fj["sets"][a.set]:
            c = special_copy(fam, a.seed)
            fh.write(f"{c['chrom']}\t{c['start']}\t{c['end']}\t{fam['family']}:{c['gene_name']}\n")
    subprocess.run(f"bedtools maskfasta -fi {GENOME} -bed {a.out_dir}/mask.bed -fo {a.out_dir}/masked.fa && "
                   f"samtools faidx {a.out_dir}/masked.fa", shell=True, check=True)
    print(f"[realworld] {a.out_dir}/masked.fa", flush=True)
    return 0


def cmd_build(a):
    import pysam
    from sim import simulate_reads, stable_seed, write_fastq
    spec = importlib.util.spec_from_file_location("o3_fewread_world", os.path.join(BENCH, "o3_fewread", "world.py"))
    FW = importlib.util.module_from_spec(spec)            # by path: bench/o3_graph/world.py would shadow `import world`
    spec.loader.exec_module(FW)                           # (COPIES_HEADER lives in bench/o3_fewread/world.py)
    fj = json.load(open(a.families))
    fams = fj["sets"][a.set]
    prof = json.load(open(a.sim_profile)) if a.sim_profile else None   # a calibrated read profile (calibrate.py profile)
    out = a.out
    os.makedirs(out, exist_ok=True)
    gpath = a.genome or GENOME
    for ext in ("", ".fai"):
        dst = f"{out}/genome.ref.fa{ext}"
        if os.path.lexists(dst):
            os.unlink(dst)
        os.symlink(gpath + ext, dst)
    true_fa = pysam.FastaFile(GENOME)
    world_fa = pysam.FastaFile(gpath)
    world, kinds, cfam, rows_c, isoj, regions = {}, {}, {}, [], collections.defaultdict(list), []
    reads, cat_rows, cat_fa = [], [], []
    cross, n_proj = {}, collections.Counter()             # CROSS: {copy: {proj name: exons}}; projectable paralog isoforms
    n_het = collections.Counter()                         # --het: heterozygous sites per family
    for fam in fams:
        f = fam["family"]
        world[f] = a.scenario
        sp = special_copy(fam, a.seed)
        rng = random.Random(f"{a.seed}:{f}:snp")
        expressed = []
        idx = 0
        silent = set()
        if a.scenario == "SILENT":                          # a seeded third of the copies express nothing (all catalogued)
            ids = [f"{f}.{k}" for k in range(len(fam["copies"]))]
            silent = set(random.Random(f"{a.seed}:{f}:silent").sample(ids, max(1, len(ids) // 3)))
        for k, c in enumerate(fam["copies"]):
            cid = f"{f}.{k}"
            special = c is sp
            kind = "present"
            if special and a.scenario == "REMOVED":
                kind = "absent"
            elif special and a.scenario == "OUTSIDE":
                kind = "outside"
            elif cid in silent:
                kind = "silent"
            kinds[cid], cfam[cid] = kind, f
            in_ref = not (special and a.scenario == "REMOVED")
            catalogued = kind in ("present", "silent")
            rows_c.append([cid, f, c["gene_id"], c["gene_name"], c["chrom"], c["start"], c["end"], c["strand"], kind, int(in_ref),
                           int(catalogued)])
            if kind != "silent":
                expressed.append((cid, c, {}, None, ""))
            if catalogued:
                ex = [tuple(x) for x in c["transcripts"][c["main"]]]
                seq = "".join(world_fa.fetch(c["chrom"], s, e).upper() for s, e in ex)
                cat_rows.append([f, str(idx), cid, c["chrom"], str(ex[0][0]), str(ex[-1][1]), str(len(ex)), c["strand"], "0",
                                 ",".join(f"{s}-{e}" for s, e in ex), "1.000", "realworld", c["gene_name"], "NA", "1",
                                 str(len(seq)), "0.000", "kept", str(c["start"]), str(c["end"])])
                cat_fa.append(f">{f}|{idx}|{c['chrom']}:{ex[0][0]}-{ex[-1][1]}|{c['strand']}|nexon={len(ex)}\n{seq}\n")
                isoj[f].append(dict(copy=cid, idx=idx, contig=c["chrom"], start=c["start"], end=c["end"], strand=c["strand"],
                                    rep=c["main"], isoforms=c["transcripts"]))
                idx += 1
        hull = collections.defaultdict(list)                 # copy_assign wants a family's span on a chromosome in ONE region
        for c in fam["copies"]:
            hull[c["chrom"]] += [c["start"], c["end"]]
        regions += [(ch, max(0, min(v) - REGION_PAD), max(v) + REGION_PAD) for ch, v in hull.items()]
        if a.scenario == "CROSS":
            for k, c in enumerate(fam["copies"]):
                lo = max(0, c["start"] - PROJ_PAD)
                locus = true_fa.fetch(c["chrom"], lo, c["end"] + PROJ_PAD).upper()
                queries = {}
                for k2, p in enumerate(fam["copies"]):
                    if p is c:
                        continue
                    fetch_p = lambda s_, e_, ch=p["chrom"]: true_fa.fetch(ch, s_, e_)  # noqa: E731
                    for t, ex in p["transcripts"].items():
                        queries[f"proj:{f}.{k2}|{t}"] = (splice(fetch_p, [tuple(x) for x in ex], p["strand"], {}), len(ex))
                if a.projection == "liftover":
                    proj = {}
                    c_hi = c["end"] + PROJ_PAD
                    c_seq, back = oriented(lambda s_, e_, ch=c["chrom"]: true_fa.fetch(ch, s_, e_), lo, c_hi, c["strand"])
                    c_fetch = lambda s_, e_, ch=c["chrom"]: true_fa.fetch(ch, s_, e_)  # noqa: E731
                    for k2, p in enumerate(fam["copies"]):
                        if p is c:
                            continue
                        p_lo, p_hi = max(0, p["start"] - PROJ_PAD), p["end"] + PROJ_PAD
                        p_seq, _ = oriented(lambda s_, e_, ch=p["chrom"]: true_fa.fetch(ch, s_, e_), p_lo, p_hi, p["strand"])
                        ts = sorted(p["transcripts"])
                        lifted = liftover(p_seq, [to_oriented([tuple(x) for x in p["transcripts"][t]], p_lo, p_hi, p["strand"])
                                                  for t in ts], c_seq)
                        for t, ex in zip(ts, lifted):
                            g = back(ex) if ex else None
                            ok = (g is not None and len(g) == len(p["transcripts"][t]) and canonical_junctions(c_fetch, g, c["strand"]))
                            proj[f"proj:{f}.{k2}|{t}"] = g if ok else None
                else:
                    proj = project_many(queries, locus, lo, c["strand"]) if queries else {}
                own = {t: [tuple(x) for x in v] for t, v in c["transcripts"].items()}
                cross[f"{f}.{k}"] = cross_structures(own, proj, MAX_PROJ, random.Random(f"{a.seed}:{f}.{k}:cross"))
                n_proj[f] += sum(v is not None for v in proj.values())
        if a.scenario == "BRANCH":
            j = mid_junction(sp["transcripts"][sp["main"]])
            if j:
                exonic = sorted({p for v in sp["transcripts"].values() for s, e in v for p in range(s, e)})
                snps = {p: rng.choice([b for b in "ACGT" if b != true_fa.fetch(sp["chrom"], p, p + 1).upper()])
                        for p in rng.sample(exonic, int(round(BRANCH_SNP * len(exonic))))}
                novel = "".join(rng.choice("ACGT") for _ in range(NOVEL_LEN))
                cid = f"{f}.B"
                kinds[cid], cfam[cid] = "absent", f
                rows_c.append([cid, f, sp["gene_id"] + "+novel", sp["gene_name"], sp["chrom"], sp["start"], sp["end"], sp["strand"],
                               "absent", 0, 0])
                expressed.append((cid, sp, snps, j, novel))
        for cid, c, snps, j, novel in expressed:
            fetch = lambda s, e, ch=c["chrom"]: true_fa.fetch(ch, s, e)  # noqa: E731
            tx = c["transcripts"]
            w = {t: (MAIN_WEIGHT if t == c["main"] else OTHER_WEIGHT) for t in tx}
            tot = sum(w.values())
            het, hom = {}, {}
            het_rate = a.het if a.het > 0 else (2 * prof["het_per_base"] if prof else 0.0)
            exonic = {p for v in tx.values() for s_, e_ in v for p in range(s_, e_)}
            base_at = lambda p, ch=c["chrom"]: true_fa.fetch(ch, p, p + 1)  # noqa: E731
            if het_rate > 0:                                    # a diploid donor: haplotype 2 carries the copy's het alleles
                het = het_sites(exonic, het_rate, random.Random(f"{a.seed}:{cid}:het"), base_at)
            if prof:                                            # profile: the donor's homozygous alleles, on both haplotypes
                hom = het_sites(exonic, prof["hom_per_base"], random.Random(f"{a.seed}:{cid}:hom"), base_at)
            n_het[f] += len(het)
            mult = expr_mult(a.seed, cid, prof.get("expr_log_sd", 0.0)) if prof else 1.0
            cache = {}

            def body_of(ex, h, r, own_snps, nov, jj):
                """profile mode: retained introns, a calibrated (5', 3') truncation; the transcript sequence, then jitter."""
                ex = [tuple(x) for x in ex]
                if not nov and len(ex) > 1:
                    keep = {k for k in range(len(ex) - 1) if r.random() < prof["kept_intron_per_intron"]}
                    ex = retain(ex, keep) if keep else ex
                key = (tuple(ex), h)
                if key not in cache:
                    al = {**hom, **own_snps, **het} if h == 2 else {**hom, **own_snps}
                    cache[key] = splice(fetch, ex, c["strand"], al, jj, nov)
                if prof.get("truncation_redraw"):              # profile v2 (round 9 review I5): redraw a truncation too long
                    return draw_truncation(cache[key], prof["end_pairs"], r)
                m5, m3 = r.choice(prof["end_pairs"])            # profile v1: one draw (whole transcript when < MIN_READ)
                return truncate(cache[key], m5, m3)

            def hap_of(t, i):                                   # its own RNG: hap-1 reads equal the het-free world's
                return 2 if het and random.Random(stable_seed(str(a.seed), cid, t, str(i), "hap")).random() < 0.5 else 1
            for t, ex in tx.items():
                seqs = {h: splice(fetch, [tuple(x) for x in ex], c["strand"], {**snps, **het} if h == 2 else snps, j, novel)
                        for h in ((1, 2) if het else (1,))}
                n = int(round(a.per_copy * w[t] / tot * mult))
                for i in range(n):
                    h = hap_of(t, i)
                    r = random.Random(stable_seed(str(a.seed), cid, t, str(i)))
                    seq = body_of(ex, h, r, snps, novel, j) if prof else seqs[h]
                    lo, hi = r.randint(0, 30), len(seq) - r.randint(0, 30)
                    if hi - lo < 100:
                        lo, hi = 0, len(seq)
                    if prof:
                        b = seq[lo:hi]
                        b = apply_edits(b, draw_edits(b, r, prof))
                        rd = (b, "~" * len(b))
                    else:
                        rd = simulate_reads(seq[lo:hi], 1, err=0.001, indel=0.0003, seed=r.getrandbits(30))[0]
                    reads.append((cid, t, f, kinds[cid], rd, int(carries(novel, seq[lo:hi])), h))
            for t, ex in cross.get(cid, {}).items():            # CROSS: the paralog structures this copy also expresses
                seqs = {h: splice(fetch, ex, c["strand"], het if h == 2 else {}) for h in ((1, 2) if het else (1,))}
                for i in range(int(round(a.per_copy * OTHER_WEIGHT / tot * mult))):
                    h = hap_of(t, i)
                    r = random.Random(stable_seed(str(a.seed), cid, t, str(i)))
                    seq = body_of(ex, h, r, {}, "", None) if prof else seqs[h]
                    lo, hi = r.randint(0, 30), len(seq) - r.randint(0, 30)
                    if hi - lo < 100:
                        lo, hi = 0, len(seq)
                    if prof:
                        b = seq[lo:hi]
                        b = apply_edits(b, draw_edits(b, r, prof))
                        rd = (b, "~" * len(b))
                    else:
                        rd = simulate_reads(seq[lo:hi], 1, err=0.001, indel=0.0003, seed=r.getrandbits(30))[0]
                    reads.append((cid, t, f, kinds[cid], rd, 0, h))
    names = anon_names([(r[0], r[1]) for r in reads], a.seed)
    order = sorted(range(len(reads)), key=lambda i: names[i])
    with open(f"{out}/reads.E0.fq", "w") as fq, open(f"{out}/reads.truth.E0.tsv", "w") as tr:
        tr.write("read\tcopy\tkind\tfamily\tisoform\tnovel\thap\n")
        for i in order:
            cid, t, f, kd, rd, nv, h = reads[i]
            write_fastq(fq, names[i], rd)
            tr.write(f"{names[i]}\t{cid}\t{kd}\t{f}\t{t}\t{nv}\t{h}\n")
    with open(f"{out}/P.fam.copies.tsv", "w") as fh:
        fh.write(FW.COPIES_HEADER + "\n" + "".join("\t".join(r) + "\n" for r in cat_rows))
    open(f"{out}/P.fam.copies.fa", "w").write("".join(cat_fa))
    json.dump(isoj, open(f"{out}/isoforms.json", "w"), indent=1)
    json.dump(kinds, open(f"{out}/kinds.json", "w"), indent=1)
    json.dump(cfam, open(f"{out}/copy_family.json", "w"), indent=1)
    json.dump(world, open(f"{out}/world.json", "w"), indent=1)
    if a.scenario == "CROSS":
        json.dump(dict(projection=a.projection, projectable=n_proj, expressed=cross), open(f"{out}/cross.json", "w"), indent=1)
        print(f"[realworld] CROSS: {sum(len(v) for v in cross.values())} paralog structures expressed over {len(cross)} copies "
              f"({sum(n_proj.values())} paralog isoforms projectable)", flush=True)
    with open(f"{out}/copies.tsv", "w") as fh:
        fh.write("copy\tfamily\tgene_id\tgene_name\tchrom\tstart\tend\tstrand\tkind\tin_reference\tcatalogued\n")
        fh.writelines("\t".join(str(x) for x in r) + "\n" for r in rows_c)
    merged = []
    for ch, s, e in sorted(regions):
        if merged and merged[-1][0] == ch and s <= merged[-1][2]:
            merged[-1] = (ch, merged[-1][1], max(merged[-1][2], e))
        else:
            merged.append((ch, s, e))
    with open(f"{out}/regions.txt", "w") as fh:
        fh.writelines(f"{ch}:{s + 1}-{e}\n" for ch, s, e in merged)
    if a.het > 0:
        json.dump(dict(het=a.het, sites=n_het), open(f"{out}/het.json", "w"), indent=1)
    print(f"[realworld] {out}: {len(fams)} families, {len(cat_rows)} catalogued copies, {len(reads)} reads"
          + (f", {sum(n_het.values())} heterozygous sites (--het {a.het})" if a.het > 0 else ""), flush=True)
    return 0


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd", required=True)
    p = sub.add_parser("mask")
    p.add_argument("out_dir")
    p.add_argument("families")
    p.add_argument("--set", required=True)
    p.add_argument("--seed", type=int, required=True)
    p.set_defaults(func=cmd_mask)
    p = sub.add_parser("build")
    p.add_argument("out")
    p.add_argument("families")
    p.add_argument("--set", required=True)
    p.add_argument("--scenario", required=True, choices=SCENARIOS)
    p.add_argument("--seed", type=int, required=True)
    p.add_argument("--genome")
    p.add_argument("--per-copy", type=int, default=60)
    p.add_argument("--het", type=float, default=0.0,
                   help="a diploid donor: each copy's exonic heterozygous-site rate; reads come 50/50 from its two haplotypes")
    p.add_argument("--sim-profile", help="a calibrated read profile (calibrate.py profile + expr_log_sd): truncation, "
                   "homopolymer-biased errors, the donor's het and hom alleles, retained introns, per-copy expression spread")
    p.add_argument("--projection", choices=("spliced", "liftover"), default="spliced",
                   help="CROSS: paralog isoforms onto a copy by spliced transcript alignment (Amendment 1) or by liftover "
                        "through the loci's genomic alignment (independent of union O2's projection)")
    p.set_defaults(func=cmd_build)
    a = ap.parse_args(argv)
    return a.func(a)


if __name__ == "__main__":
    sys.exit(main())
