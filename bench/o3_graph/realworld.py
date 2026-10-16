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
  CROSS     as INTACT, and EVERY copy also expresses up to MAX_PROJ isoform structures annotated only on a paralog (review C1
            of register 1298: the catalog is not copy-specific): each paralog isoform's sense sequence spliced-aligned
            (minimap2 splice:hq -uf) to the copy's locus, kept when its primary alignment covers >= PROJ_COVER of it with the
            same exon count, canonical junctions and an intron chain none of the copy's own isoforms has; the copy's own
            genomic sequence over those exons is the transcript (truth isoform `proj:<paralog copy>|<isoform>`, weight
            OTHER; the copy's own-isoform reads are INTACT's). The catalog lists each copy's own isoforms only.
Reads: per expressed copy PER_COPY reads over its isoforms (main 2.0, others 0.5), each the isoform's sense sequence with 0-30
bp end jitter and HiFi errors (bench/sim.py, 0.1% substitutions, 0.03% indels), from the UNMASKED genome. Names are a
seeded shuffle `s0000001..` (they carry no truth); reads.truth.E0.tsv holds it, with `novel` = 1 for a read that
carries the novel exon (the error-free sequence before errors contains it whole).
Products (the layout detect.py / score.py read): genome.ref.fa (+ .fai; a link to the world's genome), reads.E0.fq,
reads.truth.E0.tsv, P.fam.copies.tsv / .fa (catalog; representative = the main isoform), isoforms.json, kinds.json,
copy_family.json, world.json, copies.tsv, regions.txt (each family's span per chromosome +- 20 kb, merged).
"""
import argparse
import collections
import json
import os
import random
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
BENCH = os.path.dirname(HERE)
sys.path.insert(0, BENCH)
sys.path.insert(0, os.path.join(BENCH, "o3_fewread"))

GENOME = "/mnt/linuxdisk/home/juanfraitu/npip_ladder/idx/target.fa"      # = chm13v2.0.fa (same contigs and sequence)
SCENARIOS = ("INTACT", "REMOVED", "OUTSIDE", "BRANCH", "CROSS")
MAIN_WEIGHT, OTHER_WEIGHT, NOVEL_LEN, BRANCH_SNP, REGION_PAD = 2.0, 0.5, 150, 0.002, 20000
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
    if len(names) > max_proj:
        names = sorted((rng or random.Random(0)).sample(names, max_proj))
    return {n: pick[n] for n in names}


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
    import world as FW                                   # bench/o3_fewread/world.py: COPIES_HEADER
    fj = json.load(open(a.families))
    fams = fj["sets"][a.set]
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
    for fam in fams:
        f = fam["family"]
        world[f] = a.scenario
        sp = special_copy(fam, a.seed)
        rng = random.Random(f"{a.seed}:{f}:snp")
        expressed = []
        idx = 0
        for k, c in enumerate(fam["copies"]):
            cid = f"{f}.{k}"
            special = c is sp
            kind = "present"
            if special and a.scenario == "REMOVED":
                kind = "absent"
            elif special and a.scenario == "OUTSIDE":
                kind = "outside"
            kinds[cid], cfam[cid] = kind, f
            in_ref = not (special and a.scenario == "REMOVED")
            catalogued = kind == "present"
            rows_c.append([cid, f, c["gene_id"], c["gene_name"], c["chrom"], c["start"], c["end"], c["strand"], kind, int(in_ref),
                           int(catalogued)])
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
            for t, ex in tx.items():
                seq = splice(fetch, [tuple(x) for x in ex], c["strand"], snps, j, novel)
                n = int(round(a.per_copy * w[t] / tot))
                for i in range(n):
                    r = random.Random(stable_seed(str(a.seed), cid, t, str(i)))
                    lo, hi = r.randint(0, 30), len(seq) - r.randint(0, 30)
                    if hi - lo < 100:
                        lo, hi = 0, len(seq)
                    rd = simulate_reads(seq[lo:hi], 1, err=0.001, indel=0.0003, seed=r.getrandbits(30))[0]
                    reads.append((cid, t, f, kinds[cid], rd, int(carries(novel, seq[lo:hi]))))
            for t, ex in cross.get(cid, {}).items():            # CROSS: the paralog structures this copy also expresses
                seq = splice(fetch, ex, c["strand"], {})
                for i in range(int(round(a.per_copy * OTHER_WEIGHT / tot))):
                    r = random.Random(stable_seed(str(a.seed), cid, t, str(i)))
                    lo, hi = r.randint(0, 30), len(seq) - r.randint(0, 30)
                    if hi - lo < 100:
                        lo, hi = 0, len(seq)
                    rd = simulate_reads(seq[lo:hi], 1, err=0.001, indel=0.0003, seed=r.getrandbits(30))[0]
                    reads.append((cid, t, f, kinds[cid], rd, 0))
    names = anon_names([(r[0], r[1]) for r in reads], a.seed)
    order = sorted(range(len(reads)), key=lambda i: names[i])
    with open(f"{out}/reads.E0.fq", "w") as fq, open(f"{out}/reads.truth.E0.tsv", "w") as tr:
        tr.write("read\tcopy\tkind\tfamily\tisoform\tnovel\n")
        for i in order:
            cid, t, f, kd, rd, nv = reads[i]
            write_fastq(fq, names[i], rd)
            tr.write(f"{names[i]}\t{cid}\t{kd}\t{f}\t{t}\t{nv}\n")
    with open(f"{out}/P.fam.copies.tsv", "w") as fh:
        fh.write(FW.COPIES_HEADER + "\n" + "".join("\t".join(r) + "\n" for r in cat_rows))
    open(f"{out}/P.fam.copies.fa", "w").write("".join(cat_fa))
    json.dump(isoj, open(f"{out}/isoforms.json", "w"), indent=1)
    json.dump(kinds, open(f"{out}/kinds.json", "w"), indent=1)
    json.dump(cfam, open(f"{out}/copy_family.json", "w"), indent=1)
    json.dump(world, open(f"{out}/world.json", "w"), indent=1)
    if a.scenario == "CROSS":
        json.dump(dict(projectable=n_proj, expressed=cross), open(f"{out}/cross.json", "w"), indent=1)
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
    print(f"[realworld] {out}: {len(fams)} families, {len(cat_rows)} catalogued copies, {len(reads)} reads", flush=True)
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
    p.set_defaults(func=cmd_build)
    a = ap.parse_args(argv)
    return a.func(a)


if __name__ == "__main__":
    sys.exit(main())
