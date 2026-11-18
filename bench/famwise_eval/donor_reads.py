#!/usr/bin/env python3
"""PREREG r1797 / r1798 (docs/PREREG_famwise_sim_2026-11-18.md): reads drawn from an individual's own copies of the
truth-family genes, each with its member, haplotype and isoform of origin.

  donor_reads.py families CAT.gff3.gz OUT.json                          (env as famtruth: PANEL_SPECIES, DENOVO_*)
  donor_reads.py queries FAM.json GENOME.fa OUT.fa [--flank 200000]
  donor_reads.py counterparts FAM.json QUERIES.fa OUT.tsv ACCOUNT.tsv HAP=PAF [HAP=PAF ...]
  donor_reads.py counts FAM.json SEL.bam OUT.tsv
  donor_reads.py build FAM.json COUNTERPARTS.tsv COUNTS.tsv SCENARIO PROFILE.json OUT_DIR HAP=FASTA [...] [--seed 2026]
  donor_reads.py build FAM.json - COUNTS.tsv SCENARIO PROFILE.json OUT_DIR ref=GENOME.fa [--seed 2026]    (reference-drawn control)

  donor_reads.py templates FAM.json SEL.bam OUT.tsv                                                       (r1809)
  donor_reads.py build-templates FAM.json QUERIES.fa ACCOUNT.tsv TEMPLATES.tsv SCENARIO PROFILE.json OUT_DIR HAP=PAF,FASTA [...] [--seed 2026]
  donor_reads.py build-templates FAM.json - - TEMPLATES.tsv SCENARIO PROFILE.json OUT_DIR ref=GENOME.fa   (reference-drawn control)

r1809: a read copies the aligned blocks of one real read's record at its member (templates), lifted onto a haplotype.
PAF: minimap2 -x asm20 -c --secondary=no of QUERIES.fa on a haplotype. SEL.bam: every record of the real reads the
expression is taken from. SCENARIO: observed | to-major | to-minor. OUT_DIR gets reads.fa, origin.tsv and build.json.
Borrowed, not edited: bench/o3_graph/donorsim.py (read_paf, lift_exons) and realworld.py (splice, retain, draw_edits)."""
import collections, csv, importlib.util, json, os, random, sys
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import records as R

def _load(name):
    spec = importlib.util.spec_from_file_location(f"fwsim_{name}", os.path.join(HERE, "..", "o3_graph", f"{name}.py"))
    mod = importlib.util.module_from_spec(spec); spec.loader.exec_module(mod)
    return mod

DS = _load("donorsim"); RW = DS.RW
MIN_MAPQ, NEAR, ANCHOR, MIN_READ, TRIES = 20, 10, 20, 200, 10
SCENARIOS = ("observed", "to-major", "to-minor")

# ---- families
def merge_abutting(exons):
    out = []
    for a, b in sorted(exons):
        if out and a <= out[-1][1]:
            out[-1][1] = max(out[-1][1], b)
        else:
            out.append([a, b])
    return [tuple(x) for x in out]

def canonical(tx):
    """most exons, then the longest span, then the first id (bench/o3_graph/families.py canonical)"""
    return min(tx, key=lambda t: (-len(tx[t]), -(tx[t][-1][1] - tx[t][0][0]), t))

def family_json(members, transcripts):
    """members = [(family, gene_id, name, chrom, strand, biotype)], transcripts = {gene_id: {tid: exons}} -> donorsim FAMILIES.json"""
    fams = collections.defaultdict(list)
    for fam, g, name, chrom, strand, bio in members:
        tx = {t: merge_abutting(e) for t, e in transcripts[g].items()}
        fams[fam].append(dict(chrom=chrom, start=min(e[0][0] for e in tx.values()), end=max(e[-1][1] for e in tx.values()), strand=strand,
                              main=canonical(tx), transcripts={t: [list(x) for x in e] for t, e in tx.items()}, biotype=bio, gene_id=g, gene_name=name))
    return {"sets": {"truth": [dict(family=f, copies=c) for f, c in sorted(fams.items())]}}

def copies(fam):
    """[(gene_id, family, copy)] with exons as tuples"""
    out = []
    for f in next(iter(fam["sets"].values())):
        for c in f["copies"]:
            out.append((c["gene_id"], f["family"], dict(c, transcripts={t: [tuple(e) for e in ex] for t, ex in c["transcripts"].items()})))
    return out

def cmd_families(gff_p, out_p):
    import pysam
    sys.path.insert(0, os.path.join(HERE, "..", "denovo_eval"))
    import famtruth as F
    G = F.load_genes(os.environ.get("PANEL_SPECIES", "human"))
    pairs, _ = F.truth_pairs(G, F.PARALOGS, set(F.LEVELS["simiiformes"]))
    gff, members, tx = pysam.TabixFile(gff_p), [], {}
    for comp in F.components(pairs):
        for g in sorted(comp):
            v = G[g]
            if not v["exons"]:
                continue
            tids, ex = set(), collections.defaultdict(list)
            for line in gff.fetch(v["chrom"], v["exons"][0][0], v["exons"][-1][1]):
                f = line.split("\t"); a = dict(x.split("=", 1) for x in f[8].split(";") if "=" in x)
                if f[2] == "transcript" and a.get("Parent") == g:
                    tids.add(a["ID"])
                elif f[2] == "exon":
                    ex[a.get("Parent")].append((int(f[3]) - 1, int(f[4])))
            tx[g] = {t: ex[t] for t in tids if ex.get(t)}
            assert tx[g], f"{g}: no transcript in {gff_p}"
            members.append((min(comp), g, v["name"], v["chrom"], v["strand"], v["biotype"]))
    out = family_json(members, tx)
    json.dump(out, open(out_p, "w"))
    print(f"{len(out['sets']['truth'])} families, {len(members)} members, {sum(map(len, tx.values()))} isoforms -> {out_p}")

def cmd_queries(fam_p, genome_p, out_p, flank):
    import pysam
    fa, n = pysam.FastaFile(genome_p), 0
    with open(out_p, "w") as fh:
        for g, _, c in copies(json.load(open(fam_p))):
            lo, hi = max(0, c["start"] - flank), min(fa.get_reference_length(c["chrom"]), c["end"] + flank)
            fh.write(f">{g}|{c['chrom']}|{lo}\n{fa.fetch(c['chrom'], lo, hi).upper()}\n"); n += 1
    print(f"{n} flanked member queries (flank {flank}) -> {out_p}")

# ---- counterparts
def pick_row(rows, chrom, gs, ge):
    """the primary row on the homologous chromosome overlapping most of the member's own span [gs, ge) of the query; tie: longer target"""
    ok = [r for r in rows if r["tp"] == "P" and r["tname"].split("_")[0] == chrom and min(r["qe"], ge) > max(r["qs"], gs)]
    return max(ok, key=lambda r: (min(r["qe"], ge) - max(r["qs"], gs), r["te"] - r["ts"])) if ok else None

def counterpart(rows, c, lo):
    """-> (counterpart or None, state, row). counterpart = {tname, strand (of the transcript on the haplotype), start, end, isoforms}"""
    best = pick_row(rows, c["chrom"], c["start"] - lo, c["end"] - lo)
    if best is None:
        return None, "no row on the homologous chromosome", None
    if best["mapq"] < MIN_MAPQ:
        return None, "mapq < 20", best
    iso = {}
    for t, ex in c["transcripts"].items():
        l = DS.lift_exons(best, [(a - lo, b - lo) for a, b in ex])
        if l:
            iso[t] = [tuple(x) for x in l]
    if not iso:
        return None, "no isoform lifts", best
    strand = c["strand"] if best["strand"] == "+" else "-+"[c["strand"] == "-"]
    return dict(tname=best["tname"], strand=strand, start=min(e[0][0] for e in iso.values()), end=max(e[-1][1] for e in iso.values()), isoforms=iso), "counterpart", best

def collapsed(cps, fam_of, ref):
    """members of one haplotype to drop: same family, lifts on one contig and strand overlapping >= 50% of the shorter, and
    not overlapping on the reference (ref = {member: (chrom, start, end)})"""
    bad, items = set(), sorted(cps.items())
    for i, (a, x) in enumerate(items):
        for b, y in items[i + 1:]:
            if fam_of[a] != fam_of[b] or x["tname"] != y["tname"] or x["strand"] != y["strand"]:
                continue
            ov = min(x["end"], y["end"]) - max(x["start"], y["start"])
            if 2 * ov < min(x["end"] - x["start"], y["end"] - y["start"]):
                continue
            (ca, sa, ea), (cb, sb, eb) = ref[a], ref[b]
            if not (ca == cb and sa < eb and sb < ea):
                bad |= {a, b}
    return bad

def cmd_counterparts(fam_p, queries_p, out_p, acct_p, haps):
    cs = copies(json.load(open(fam_p)))
    lo_of = {l[1:].split("|")[0]: (l[1:].strip(), int(l.strip().split("|")[2])) for l in open(queries_p) if l.startswith(">")}
    fam_of, ref = {g: f for g, f, _ in cs}, {g: (c["chrom"], c["start"], c["end"]) for g, _, c in cs}
    with open(out_p, "w") as out, open(acct_p, "w") as acct:
        out.write("member\thap\tcontig\tstrand\tisoforms\n")
        acct.write("member\tname\tfamily\thap\tstate\trow_target_span\tmember_coverage\tde\tmapq\tisoforms_lifted\tisoforms\n")
        for h, paf in haps.items():
            rows, cps, info = DS.read_paf(paf), {}, {}
            for g, f, c in cs:
                q, lo = lo_of[g]
                cp, state, best = counterpart(rows.get(q, []), c, lo)
                cov = (min(best["qe"], c["end"] - lo) - max(best["qs"], c["start"] - lo)) / (c["end"] - c["start"]) if best else 0.0
                info[g] = [state, best["te"] - best["ts"] if best else "", f"{cov:.3f}", best["de"] if best else "", best["mapq"] if best else "", len(cp["isoforms"]) if cp else 0]
                if cp:
                    cps[g] = cp
            for g in collapsed(cps, fam_of, ref):
                info[g][0] = "collapsed on the haplotype"; del cps[g]
            n = collections.Counter(v[0] for v in info.values())
            print(h, dict(sorted(n.items())))
            for g, f, c in cs:
                acct.write("\t".join(map(str, [g, c["gene_name"], f, h] + info[g] + [len(c["transcripts"])])) + "\n")
                if g in cps:
                    out.write(f"{g}\t{h}\t{cps[g]['tname']}\t{cps[g]['strand']}\t{json.dumps(cps[g]['isoforms'], separators=(',', ':'))}\n")

def reference_counterparts(fam):
    """the reference-drawn control: every member is its own counterpart on 'ref'"""
    return {g: {"ref": dict(tname=c["chrom"], strand=c["strand"], isoforms=c["transcripts"])} for g, _, c in copies(fam)}

def load_counterparts(path):
    out = collections.defaultdict(dict)
    for r in csv.DictReader(open(path), delimiter="\t"):
        out[r["member"]][r["hap"]] = dict(tname=r["contig"], strand=r["strand"], isoforms={t: [tuple(e) for e in ex] for t, ex in json.loads(r["isoforms"]).items()})
    return dict(out)

# ---- expression
def member_index(fam):
    """(chrom, strand) -> [(start, end, exon union, gene_id)]"""
    idx = collections.defaultdict(list)
    for g, _, c in copies(fam):
        idx[(c["chrom"], c["strand"])].append((c["start"], c["end"], merge_abutting([e for ex in c["transcripts"].values() for e in ex]), g))
    return idx

def blocks(rec):
    """a record's aligned blocks between its introns (records.Rec.chain is in gtf_names' convention: next exon start + 1)"""
    out, cur = [], rec.start
    for a, b in rec.chain:
        out.append((cur, a)); cur = b - 1
    return out + [(cur, rec.end)]

def member_of(rec, idx):
    """the member on the record's strand with the most overlapping exon bases (tie: lower gene id), or None"""
    bl, best = blocks(rec), (0, "")
    for s, e, ex, g in idx.get((rec.chrom, rec.strand), ()):
        if s < rec.end and rec.start < e:
            n = sum(max(0, min(b, y) - max(a, x)) for a, b in bl for x, y in ex)
            if n > best[0] or (n and n == best[0] and g < best[1]):
                best = (n, g)
    return best[1] or None

def isoform_weights(chains, tx):
    """chains = the junction chains ((exon end, next exon start), ...) of a member's primaries -> {tid: weight}: 1/k to each of
    the k isoforms a spliced read fits (its chain is the isoform's or a contiguous part); nothing fits: the canonical isoform"""
    ch = {t: tuple((ex[i][1], ex[i + 1][0]) for i in range(len(ex) - 1)) for t, ex in tx.items()}
    w = collections.Counter()
    for c in chains:
        fits = [t for t, full in ch.items() if c and any(full[i:i + len(c)] == tuple(c) for i in range(len(full) - len(c) + 1))]
        for t in fits:
            w[t] += 1.0 / len(fits)
    return dict(w) if w else {canonical(tx): 1.0}

def near_tied_members(recs, idx):
    """the members holding the read's places within NEAR of its best ms, when there are two or more and no place in that
    band lies off the members; else None (the read is decided, or left alone)"""
    places = []
    for r in sorted(recs, key=lambda r: -r.ms):
        if not any(p.chrom == r.chrom and p.strand == r.strand and p.start < r.end and r.start < p.end for p in places):
            places.append(r)                                    # the best record of each place stands for it
    band = [member_of(p, idx) for p in places if places[0].ms - p.ms <= NEAR]
    return set(band) if len(band) > 1 and None not in band and len(set(band)) > 1 else None

def scenario_members(reads, idx, weights_out=None):
    """reads = {name: [Rec]} -> {scenario: {read: member}} (reads whose primary is on no member are left out);
    weights_out, if given, collects {member: [junction chains of its primaries]}"""
    own, grouped, tied, parent = {}, set(), collections.defaultdict(collections.Counter), {}
    def find(x):
        while parent.setdefault(x, x) != x:
            parent[x] = parent[parent[x]]; x = parent[x]
        return x
    for q, recs in reads.items():
        p = next((r for r in recs if r.primary), None)
        m = member_of(p, idx) if p else None
        if m is None:
            continue
        own[q] = m
        if weights_out is not None:
            weights_out[m].append(tuple((a, b - 1) for a, b in p.chain))
        near = near_tied_members(recs, idx)
        if near is None or m not in near:
            continue
        for x in near:
            parent[find(x)] = find(m)
        grouped.add(q)
    groups = collections.defaultdict(collections.Counter)
    for q in grouped:
        groups[find(own[q])][own[q]] += 1
    win = {"observed": {}, "to-major": {g: min(c, key=lambda m: (-c[m], m)) for g, c in groups.items()},
           "to-minor": {g: min(c, key=lambda m: (c[m], m)) for g, c in groups.items()}}
    scenario_members.info = dict(off_members=len(reads) - len(own), near_tied=len(grouped), groups=len(groups))
    return {s: {q: win[s].get(find(m), m) if q in grouped else m for q, m in own.items()} for s in SCENARIOS}

def scenario_counts(reads, idx, weights_out=None):
    """-> {scenario: {member: reads}}"""
    mem = scenario_members(reads, idx, weights_out)
    scenario_counts.info = scenario_members.info
    return {s: dict(sorted(collections.Counter(m.values()).items())) for s, m in mem.items()}

def templates(reads, idx):
    """r1809 -> (scenario, read, member, blocks) per read: the record a simulated read copies. A read that moves copies its
    best-ms record on the member it moves to; with no record there it stays on its primary's member."""
    mem = scenario_members(reads, idx)
    templates.stayed = collections.Counter()
    for s in SCENARIOS:
        for q, m in mem[s].items():
            rec = next(r for r in reads[q] if r.primary)
            if m != mem["observed"][q]:
                there = [r for r in reads[q] if member_of(r, idx) == m]
                if there:
                    rec = max(there, key=lambda r: r.ms)
                else:
                    m = mem["observed"][q]; templates.stayed[s] += 1
            yield s, q, m, blocks(rec)

def cmd_templates(fam_p, bam_p, out_p):
    idx, n = member_index(json.load(open(fam_p))), collections.Counter()
    with open(out_p, "w") as fh:
        fh.write("scenario\tread\tmember\tblocks\n")
        for s, q, m, bl in templates(R.load_records(bam_p), idx):
            fh.write(f"{s}\t{q}\t{m}\t{','.join(f'{a}-{b}' for a, b in bl)}\n"); n[s] += 1
    print(f"templates per scenario {dict(n)}; moved reads with no record at their destination, left in place {dict(templates.stayed)}; {scenario_members.info} -> {out_p}")

def cmd_counts(fam_p, bam_p, out_p):
    fam = json.load(open(fam_p)); idx, chains = member_index(fam), collections.defaultdict(list)
    counts = scenario_counts(R.load_records(bam_p), idx, chains)
    tx = {g: c["transcripts"] for g, _, c in copies(fam)}
    with open(out_p, "w") as fh:
        fh.write("scenario\tmember\treads\tweights\n")
        for s in SCENARIOS:
            for m, n in counts[s].items():
                fh.write(f"{s}\t{m}\t{n}\t{json.dumps(isoform_weights(chains[m], tx[m]), sort_keys=True)}\n")
    spliced = [c for v in chains.values() for c in v if c]
    nofit = sum(1 for m, v in chains.items() for c in v if c and not any(c == f[i:i + len(c)] for ex in tx[m].values() for f in [tuple((ex[j][1], ex[j + 1][0]) for j in range(len(ex) - 1))] for i in range(len(f) - len(c) + 1)))
    print(f"reads per scenario {({s: sum(c.values()) for s, c in counts.items()})}; members with reads {({s: len(c) for s, c in counts.items()})}; "
          f"{scenario_counts.info}; spliced primaries {len(spliced)}, fitting no annotated isoform {nofit}")

def load_counts(path, scenario):
    return {r["member"]: (int(r["reads"]), json.loads(r["weights"])) for r in csv.DictReader(open(path), delimiter="\t") if r["scenario"] == scenario}

# ---- reads
def cells(n, weights, lifted, main):
    """reads per (haplotype, isoform): n split evenly over the haplotypes, then by weight over the isoforms lifted there
    (none of the weighted ones lifted: the member's canonical isoform if lifted, else the canonical of those that are)"""
    out = {}
    for h in sorted(lifted):
        w = {t: v for t, v in weights.items() if t in lifted[h] and v > 0}
        if not w:
            w = {main if main in lifted[h] else sorted(lifted[h])[0]: 1.0}
        for t, v in sorted(w.items()):
            k = int(round(n / len(lifted) * v / sum(w.values())))
            if k:
                out[(h, t)] = k
    return out

def cut(length, pairs, rng):
    """bases lost at the 5' and 3' end: a (5' fraction, 3' fraction) pair, redrawn while under MIN_READ bases would be left"""
    for _ in range(TRIES):
        m5, m3 = rng.choice(pairs)
        n5, n3 = int(m5 * length), int(m3 * length)
        if length - n5 - n3 >= MIN_READ:
            return n5, n3
    return 0, 0

def is_full(junctions, length, n5, n3, retained):
    """the read carries every junction of its isoform with ANCHOR bases beyond the outer ones, and no retained intron"""
    return bool(junctions) and not retained and n5 + ANCHOR <= junctions[0] and junctions[-1] <= length - n3 - ANCHOR

def pairs_for(prof, length):
    k = sum(length >= b for b in prof["classes"])
    for d in range(len(prof["classes"]) + 1):                       # an empty class takes the nearest
        for j in (k - d, k + d):
            if prof["end_pairs"].get(str(j)):
                return prof["end_pairs"][str(j)]
    raise SystemExit("the profile has no truncation pair")

def _stochastic(lam, rng):
    return int(lam) + (rng.random() < lam - int(lam))      # the measured per-read rates already hold the Poisson noise

def with_errors(body, rng, prof):
    """the read with one real read's event rates, rounded stochastically -> (sequence, events)"""
    poisson, RW._poisson = RW._poisson, _stochastic
    try:
        sub, ins, dele = rng.choice(prof["rates"])
        edits = RW.draw_edits(body, rng, dict(prof, sub_error=sub, ins_event=ins, del_event=dele))
    finally:
        RW._poisson = poisson
    return RW.apply_edits(body, edits), len(edits)

def named(out, seed):
    order = list(range(len(out))); random.Random(f"{seed}|names").shuffle(order)
    for name, i in enumerate(order):
        out[i]["read"] = f"s{name + 1:07d}"
    return sorted(out, key=lambda r: r["read"])

def exon_union(c):
    return ",".join(f"{a}-{b}" for a, b in merge_abutting([e for ex in c["transcripts"].values() for e in ex]))

def draw_reads(fam, cps, counts, prof, fetch, seed=2026):
    """-> [dict(read, seq, chrom, strand, blocks, member, family, hap, isoform, full, n5, n3, retained, edits)], names shuffled.
    fetch = {hap: f(contig, start, end)}. Every read takes one real read's event rates, rounded stochastically."""
    out = []
    for g, f, c in copies(fam):
        if g not in counts or g not in cps:
            continue
        n, w = counts[g]
        for (h, t), k in cells(n, w, {h: set(cp["isoforms"]) for h, cp in cps[g].items()}, c["main"]).items():
            cp = cps[g][h]; ex = sorted(cp["isoforms"][t]); cache = {}
            get = lambda a, b, h=h, contig=cp["tname"]: fetch[h](contig, a, b)
            for i in range(k):
                rng = random.Random(f"{seed}|{g}|{h}|{t}|{i}")
                keep = {j for j in range(len(ex) - 1) if rng.random() < prof["kept_intron_per_intron"]}
                e2 = tuple(RW.retain(ex, keep)) if keep else tuple(ex)
                if e2 not in cache:
                    cache[e2] = RW.splice(get, list(e2), cp["strand"], {})
                seq = cache[e2]
                n5, n3 = cut(len(seq), pairs_for(prof, len(seq)), rng)
                lens = [b - a for a, b in (e2 if cp["strand"] == "+" else e2[::-1])]
                junctions = [sum(lens[:j + 1]) for j in range(len(lens) - 1)]
                body, edits = with_errors(seq[n5:len(seq) - n3], rng, prof)
                out.append(dict(seq=body, chrom=c["chrom"], strand=c["strand"], blocks=exon_union(c), member=g, family=f, hap=h, isoform=t,
                                full=int(is_full(junctions, len(seq), n5, n3, bool(keep))), n5=n5, n3=n3, retained=len(keep), edits=edits))
    return named(out, seed)

# ---- r1809: reads that copy a real read's record
def lift_many(row, block_lists, max_growth=DS.MAX_EXON_GROWTH):
    """donorsim.lift_exons for many block lists in one walk of the row's alignment -> [lifted blocks or None]"""
    pos = sorted({x for bl in block_lists for a, b in bl for x in (a, b - 1)})
    if row["strand"] == "+":
        m = RW.lift_positions(row["cg"], row["qs"], row["ts"], pos)
        lift = lambda a, b: (m[a], m[b - 1] + 1) if m[a] is not None and m[b - 1] is not None else None
    else:
        m = RW.lift_positions(row["cg"], row["qlen"] - row["qe"], row["ts"], [row["qlen"] - 1 - p for p in pos])
        rc = lambda p: m[row["qlen"] - 1 - p]
        lift = lambda a, b: (rc(b - 1), rc(a) + 1) if rc(a) is not None and rc(b - 1) is not None else None
    out = []
    for bl in block_lists:
        new = [lift(a, b) for a, b in bl]
        ok = all(e is not None and e[0] < e[1] and abs((e[1] - e[0]) - (b - a)) <= max_growth for e, (a, b) in zip(new, bl))
        new = sorted(new) if ok else None
        out.append(new if new and all(x[1] <= y[0] for x, y in zip(new, new[1:])) else None)
    return out

def draw_template_reads(fam, tpl, lifts, prof, fetch, seed=2026):
    """tpl = [(real read, member, blocks on the reference)]; lifts = {member: {hap: (contig, strand of the transcript there,
    [lifted blocks or None, one per template of that member in tpl order])}} -> the reads, as draw_reads, with `template`.
    The haplotype is drawn evenly among those where the blocks lift; a read that lifts on none is not drawn."""
    by, out = {g: (f, c) for g, f, c in copies(fam)}, []
    chains = {g: {tuple((ex[i][1], ex[i + 1][0]) for i in range(len(ex) - 1)): t for t, ex in sorted(c["transcripts"].items(), reverse=True)} for g, (f, c) in by.items()}
    seen = collections.Counter()
    for q, g, bl in tpl:
        i = seen[g]; seen[g] += 1
        haps = sorted(h for h, (_, _, lifted) in lifts.get(g, {}).items() if lifted[i])
        if not haps:
            continue
        rng = random.Random(f"{seed}|{q}"); h = rng.choice(haps); contig, strand, lifted = lifts[g][h]
        f, c = by[g]
        chain = tuple((bl[j][1], bl[j + 1][0]) for j in range(len(bl) - 1))
        t = chains[g].get(chain) if chain else None
        full = int(t is not None and min(bl[0][1] - bl[0][0], bl[-1][1] - bl[-1][0]) >= ANCHOR)
        seq, edits = with_errors(RW.splice(lambda a, b: fetch[h](contig, a, b), list(lifted[i]), strand, {}), rng, prof)
        out.append(dict(seq=seq, chrom=c["chrom"], strand=c["strand"], blocks=exon_union(c), member=g, family=f, hap=h,
                        isoform=t if full else c["main"], full=full, template=q, edits=edits))
    return named(out, seed)

def cmd_build_templates(fam_p, queries_p, acct_p, tpl_p, scenario, prof_p, out_dir, haps, seed):
    import pysam
    assert scenario in SCENARIOS, scenario
    fam, prof = json.load(open(fam_p)), json.load(open(prof_p))
    by = {g: c for g, _, c in copies(fam)}
    tpl = [(r["read"], r["member"], [tuple(map(int, b.split("-"))) for b in r["blocks"].split(",")])
           for r in csv.DictReader(open(tpl_p), delimiter="\t") if r["scenario"] == scenario]
    per = collections.defaultdict(list)
    for _, g, bl in tpl:
        per[g].append(bl)
    lifts, fetch = collections.defaultdict(dict), {}
    if "ref" in haps:                                           # the reference-drawn control: the blocks as they are
        fa = pysam.FastaFile(haps["ref"]); fetch["ref"] = lambda contig, a, b, fa=fa: fa.fetch(contig, a, b)
        for g, bls in per.items():
            lifts[g]["ref"] = (by[g]["chrom"], by[g]["strand"], bls)
    else:
        lo_of = {l[1:].split("|")[0]: (l[1:].strip(), int(l.strip().split("|")[2])) for l in open(queries_p) if l.startswith(">")}
        has = {(r["member"], r["hap"]) for r in csv.DictReader(open(acct_p), delimiter="\t") if r["state"] == "counterpart"}
        for h, v in haps.items():
            paf, fasta = v.split(",")
            rows, fa = DS.read_paf(paf), pysam.FastaFile(fasta)
            fetch[h] = lambda contig, a, b, fa=fa: fa.fetch(contig, a, b)
            for g, bls in per.items():
                if (g, h) not in has:
                    continue
                q, lo = lo_of[g]; c = by[g]
                best = pick_row(rows.get(q, []), c["chrom"], c["start"] - lo, c["end"] - lo)
                strand = c["strand"] if best["strand"] == "+" else "-+"[c["strand"] == "-"]
                lifts[g][h] = (best["tname"], strand, lift_many(best, [[(a - lo, b - lo) for a, b in bl] for bl in bls]))
    reads = draw_template_reads(fam, tpl, lifts, prof, fetch, seed)
    os.makedirs(out_dir, exist_ok=True)
    cols = ("read", "chrom", "strand", "blocks", "member", "family", "hap", "isoform", "full", "template", "edits")
    with open(f"{out_dir}/reads.fa", "w") as fa, open(f"{out_dir}/origin.tsv", "w") as og:
        og.write("\t".join(cols) + "\tlength\n")
        for r in reads:
            fa.write(f">{r['read']}\n{r['seq']}\n")
            og.write("\t".join(str(r[k]) for k in cols) + f"\t{len(r['seq'])}\n")
    drawn = collections.Counter(r["member"] for r in reads); want = collections.Counter(g for _, g, _ in tpl)
    info = dict(scenario=scenario, seed=seed, model="templates", reads=len(reads), templates=len(tpl), not_drawn=len(tpl) - len(reads),
                members_emitting=len(drawn), members_without_a_lift=sum(1 for g in want if g not in drawn),
                largest_losses=sorted(((want[g] - drawn[g], by[g]["gene_name"]) for g in want if want[g] > drawn[g]), reverse=True)[:8],
                per_hap=dict(collections.Counter(r["hap"] for r in reads)), full=sum(r["full"] for r in reads))
    json.dump(info, open(f"{out_dir}/build.json", "w"), indent=1)
    print(info)

COLS = ("read", "chrom", "strand", "blocks", "member", "family", "hap", "isoform", "full", "n5", "n3", "retained", "edits")

def cmd_build(fam_p, cps_p, counts_p, scenario, prof_p, out_dir, haps, seed):
    import pysam
    assert scenario in SCENARIOS, scenario
    fam, prof, counts = json.load(open(fam_p)), json.load(open(prof_p)), load_counts(counts_p, scenario)
    cps = reference_counterparts(fam) if cps_p == "-" else load_counterparts(cps_p)
    fas = {h: pysam.FastaFile(p) for h, p in haps.items()}
    reads = draw_reads(fam, cps, counts, prof, {h: (lambda contig, a, b, fa=fa: fa.fetch(contig, a, b)) for h, fa in fas.items()}, seed)
    os.makedirs(out_dir, exist_ok=True)
    with open(f"{out_dir}/reads.fa", "w") as fa, open(f"{out_dir}/origin.tsv", "w") as og:
        og.write("\t".join(COLS) + "\tlength\n")
        for r in reads:
            fa.write(f">{r['read']}\n{r['seq']}\n")
            og.write("\t".join(str(r[k]) for k in COLS) + f"\t{len(r['seq'])}\n")
    lost = {g: n for g, (n, _) in counts.items() if g not in cps}
    info = dict(scenario=scenario, seed=seed, reads=len(reads), wanted=sum(n for n, _ in counts.values()), members_emitting=len({r["member"] for r in reads}),
                members_without_counterpart=len(lost), reads_of_members_without_counterpart=sum(lost.values()),
                per_hap=dict(collections.Counter(r["hap"] for r in reads)), full=sum(r["full"] for r in reads), with_retained_intron=sum(r["retained"] > 0 for r in reads))
    json.dump(info, open(f"{out_dir}/build.json", "w"), indent=1)
    print(info)

def main(argv):
    opt = {}
    for k, cast, default in (("--flank", int, 200000), ("--seed", int, 2026)):
        if k in argv:
            i = argv.index(k); opt[k] = cast(argv[i + 1]); argv = argv[:i] + argv[i + 2:]
        else:
            opt[k] = default
    haps = dict(a.split("=", 1) for a in argv[1:] if "=" in a)
    a = [x for x in argv if "=" not in x]
    if a[:1] == ["families"] and len(a) == 3:
        return cmd_families(*a[1:])
    if a[:1] == ["queries"] and len(a) == 4:
        return cmd_queries(*a[1:], opt["--flank"])
    if a[:1] == ["counterparts"] and len(a) == 5 and haps:
        return cmd_counterparts(*a[1:], haps)
    if a[:1] == ["counts"] and len(a) == 4:
        return cmd_counts(*a[1:])
    if a[:1] == ["templates"] and len(a) == 4:
        return cmd_templates(*a[1:])
    if a[:1] == ["build-templates"] and len(a) == 8 and haps:
        return cmd_build_templates(*a[1:], haps, opt["--seed"])
    if a[:1] == ["build"] and len(a) == 7 and haps:
        return cmd_build(*a[1:], haps, opt["--seed"])
    raise SystemExit(__doc__)

if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
