#!/usr/bin/env python3
"""SUNs vs our PSVs (M1) and a per-copy pileup check of PSVs (M2) — docs/archive/2026-10/merged/PREREG_psv_sun_pileup_2026-10-12.md.

    psv_sun.py prep            # subset S, its copy table / copies FASTA / regions -> OUT/run/
    (copy_assign --dump-psv on OUT/run/, see the prereg)
    psv_sun.py sanity          # S1: subset run vs the 2026-10-10 run
    psv_sun.py positions       # gate V: every copy's genome position per decisive column -> OUT/pos.tsv
    psv_sun.py kmers           # candidate 30-mers of private alleles -> OUT/kmers.txt (then kscan -> OUT/kmers.counts)
    psv_sun.py m1              # -> OUT/m1.json, OUT/m1_families.tsv
    psv_sun.py pileup CHROM    # per-copy pileup at the verified positions on CHROM -> OUT/pileup/CHROM.tsv
    psv_sun.py flags           # M2 flagged / control pairs -> OUT/m2_pairs.tsv
    psv_sun.py truth CHROM     # the other haplotype's base for the pairs on CHROM -> OUT/truth/CHROM.tsv
    psv_sun.py m2              # -> OUT/m2.json
Coordinates: copy exons 0-based half-open (fam.copies.tsv); copy sequences / dumped alleles in transcript orientation.
"""
import collections, csv, glob, json, os, re, subprocess, sys

R10 = "/mnt/linuxdisk/tmp/rustle_figures/runs/gorilla_KB3781/gorilla_KB3781"
BAM = "/mnt/linuxdisk/home/juanfraitu/fibroblasts/GCA_029281585.2_flnc_mm.bam"
GENOME = "/mnt/linuxdisk/home/juanfraitu/_from_wsl/winloci_scratch/GGO.fasta"
HAPS = "/mnt/linuxdisk/home/juanfraitu/gorilla_haps"
CHRMAP = "/mnt/linuxdisk/tmp/rna_allele/chrmap.tsv"
OUT = os.environ.get("PSV_SUN_OUT", "/mnt/linuxdisk/tmp/psv_sun")
SUBSET = f"{OUT}/run"  # the subset inputs written by `prep`
RUN = os.environ.get("PSV_SUN_RUN", f"{SUBSET}/s")  # a copy_assign --dump-psv prefix on that subset
RES = os.environ.get("PSV_SUN_RES", OUT)  # where `positions` writes pos.tsv / gate_v.json
K, FLANK = 30, 2000
MIN_READS, MIN_ALT, MIN_ALT_FRAC, CONTROL_REF = 10, 3, 0.20, 0.95
COMP = str.maketrans("ACGTNacgtn", "TGCANtgcan")
CS = re.compile(r"(:\d+|\*[a-z]{2}|\+[a-z]+|-[a-z]+|~[a-z]{2}\d+[a-z]{2})")


def rc(s):
    return s.translate(COMP)[::-1]


def comp(b):
    return b.translate(COMP)


def tsv(path):
    with open(path) as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def write_tsv(path, rows, cols):
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, "w") as fh:
        fh.write("\t".join(cols) + "\n")
        for r in rows:
            fh.write("\t".join(str(r[c]) for c in cols) + "\n")


class Fasta:
    """Random access through the .fai (0-based half-open fetch, upper case)."""

    def __init__(self, path):
        self.fh = open(path, "rb")
        self.idx = {}
        for line in open(path + ".fai"):
            f = line.split("\t")
            self.idx[f[0]] = (int(f[1]), int(f[2]), int(f[3]), int(f[4]))

    def fetch(self, chrom, s, e):
        n, off, lb, lw = self.idx[chrom]
        s, e = max(0, s), min(n, e)
        if s >= e:
            return ""
        a = off + (s // lb) * lw + s % lb
        b = off + ((e - 1) // lb) * lw + (e - 1) % lb + 1
        self.fh.seek(a)
        return self.fh.read(b - a).replace(b"\n", b"").replace(b"\r", b"").decode().upper()

    def length(self, chrom):
        return self.idx[chrom][0]


def parse_exons(s):
    return sorted(tuple(map(int, x.split("-"))) for x in s.split(","))


def tx_offset(exons, strand, g):
    """Genome position -> offset in the spliced transcript (transcript orientation); None if not exonic."""
    k = 0
    for s, e in exons:
        if s <= g < e:
            off = k + g - s
            L = sum(e2 - s2 for s2, e2 in exons)
            return L - 1 - off if strand == "-" else off
        k += e - s
    return None


def genome_pos(exons, strand, off):
    L = sum(e - s for s, e in exons)
    if not 0 <= off < L:
        return None
    if strand == "-":
        off = L - 1 - off
    for s, e in exons:
        if off < e - s:
            return s + off
        off -= e - s
    return None


def cs_pairs(cs, ts, q0):
    """Walk cs from target ts and query q0 (both increasing): yield (q, t) for aligned bases."""
    t, q = ts, q0
    for tok in CS.findall(cs):
        c = tok[0]
        if c == ":":
            n = int(tok[1:])
            for i in range(n):
                yield q + i, t + i
            q += n
            t += n
        elif c == "*":
            yield q, t
            q += 1
            t += 1
        elif c == "+":
            q += len(tok) - 1
        elif c == "-":
            t += len(tok) - 1
        elif c == "~":
            t += int(tok[3:-2])


def anchor_pairs(rec):
    """(query forward offset, target offset) over the aligned bases of one PAF record, either strand."""
    if rec["strand"] == "+":
        return cs_pairs(rec["cs"], rec["ts"], rec["qs"])
    # '-': cs walks the target forward against the reverse-complemented query
    qlen = rec["qlen"]
    return ((qlen - 1 - qr, t) for qr, t in cs_pairs(rec["cs"], rec["ts"], qlen - rec["qe"]))


def paf_qmap(rec):
    """{query forward offset: target offset} for one PAF record dict (strand-aware)."""
    return dict(anchor_pairs(rec))


def parse_paf(text):
    """PAF text (minimap2 stdout, or a saved copy of it) -> list of record dicts."""
    recs = []
    for line in text.splitlines():
        f = line.split("\t")
        tags = {t[:2]: t[5:] for t in f[12:]}
        recs.append({"qname": f[0], "qlen": int(f[1]), "qs": int(f[2]), "qe": int(f[3]), "strand": f[4],
                     "tname": f[5], "tlen": int(f[6]), "ts": int(f[7]), "te": int(f[8]), "tp": tags.get("tp", "P"),
                     "AS": int(tags.get("AS", 0)), "cs": tags.get("cs", "")})
    return recs


def run_paf(args):
    return parse_paf(subprocess.run(args, capture_output=True, text=True, check=True).stdout)


def canon(k):
    r = rc(k)
    return k if k < r else r


# ---------------------------------------------------------------- inputs

def subset():
    return sorted(r["family_id"] for r in tsv(R10 + ".assign.families.tsv")
                  if int(r["psv_cols"]) >= 1 and int(r["n_reads"]) >= 20 and int(r["n_copies"]) <= 20)


def copy_table():
    return {r["tid"]: r for r in tsv(R10 + ".fam.copies.tsv")}


def copies_fa(path=R10 + ".fam.copies.fa"):
    seqs, name = {}, None
    for line in open(path):
        line = line.strip()
        if line.startswith(">"):
            name = line[1:]
            seqs[name] = []
        elif name:
            seqs[name].append(line)
    out = {}
    for h, v in seqs.items():
        fam, idx = h.split("|")[:2]
        out[(fam, idx)] = "".join(v).upper()
    return out


def dump():
    cols = collections.defaultdict(dict)
    for r in tsv(RUN + ".psv_cols.tsv"):
        cols[r["family_id"]][int(r["col_index"])] = int(r["genome_pos"])
    copies = collections.defaultdict(list)
    for r in tsv(RUN + ".psv_copies.tsv"):
        copies[r["family_id"]].append(r)
    return cols, copies


def decisive(copies_rows, j):
    al = [c["alleles"][j] for c in copies_rows if j < len(c["alleles"]) and c["alleles"][j] != "."]
    return len(set(al)) >= 2


# ---------------------------------------------------------------- stages

def cmd_prep():
    S = set(subset())
    os.makedirs(f"{OUT}/run", exist_ok=True)
    with open(f"{OUT}/run/S.txt", "w") as fh:
        fh.write("\n".join(sorted(S)) + "\n")
    rows = tsv(R10 + ".fam.copies.tsv")
    cols = list(rows[0].keys())
    keep = [r for r in rows if r["family_id"] in S]
    write_tsv(f"{OUT}/run/copies.tsv", keep, cols)
    with open(R10 + ".fam.copies.fa") as fi, open(f"{OUT}/run/copies.fa", "w") as fo:
        on = False
        for line in fi:
            if line.startswith(">"):
                on = line[1:].split("|")[0] in S
            if on:
                fo.write(line)
    spans = collections.defaultdict(list)
    for r in keep:
        spans[r["chrom"]].append((int(r["start"]), int(r["end"])))
    n = 0
    with open(R10 + ".regions.txt") as fi, open(f"{OUT}/run/regions.txt", "w") as fo:
        for line in fi:
            chrom, rng = line.strip().rsplit(":", 1)
            a, b = map(int, rng.split("-"))
            if any(s < b and a - 1 < e for s, e in spans.get(chrom, ())):
                fo.write(line)
                n += 1
    print(f"S = {len(S)} families, {len(keep)} copies, {n} regions -> {OUT}/run/")


def cmd_sanity():
    S = set(open(f"{OUT}/run/S.txt").read().split())
    key = lambda r: (r["read_name"], r["family_id"], r["status"], r["assigned_copy"] if r["status"] == "assigned" else "-")
    old = collections.Counter(key(r) for r in tsv(R10 + ".assign.assignments.tsv") if r["family_id"] in S)
    new = collections.Counter(key(r) for r in tsv(RUN + ".assignments.tsv") if r["family_id"] in S)
    both = sum((old & new).values())
    res = {"rows_10_10": sum(old.values()), "rows_subset": sum(new.values()), "identical_rows": both,
           "agreement": round(both / max(sum(old.values()), 1), 4)}
    res["S1_pass"] = res["agreement"] >= 0.95
    json.dump(res, open(f"{OUT}/sanity.json", "w"), indent=1)
    print(json.dumps(res))


def cmd_positions():
    G, tab, seqs = Fasta(GENOME), copy_table(), copies_fa(f"{SUBSET}/copies.fa")
    cols, copies = dump()
    rows, stat, by_rel = [], collections.Counter(), collections.Counter()
    tmp = f"{OUT}/tmp"
    os.makedirs(tmp, exist_ok=True)
    for fam in sorted(copies):
        crs = copies[fam]
        info = []
        for c in crs:
            t = tab.get(c["copy_tid"])
            info.append(None if t is None else {"chrom": t["chrom"], "strand": t["strand"], "exons": parse_exons(t["exons"]),
                                                "seq": seqs.get((fam, t["copy_idx"]))})
        ncol = max((len(c["alleles"]) for c in crs), default=0)
        aln = {}  # anchor copy index -> {query copy index: {target offset: query offset}}

        def align_to(f):
            if f in aln:
                return aln[f]
            tq, qq = f"{tmp}/t.fa", f"{tmp}/q.fa"
            with open(tq, "w") as fh:
                fh.write(f">{f}\n{info[f]['seq']}\n")
            with open(qq, "w") as fh:
                for i, x in enumerate(info):
                    if i != f and x and x["seq"]:
                        fh.write(f">{i}\n{x['seq']}\n")
            # every primary/supplementary piece, either strand: copies differ in exon structure and orientation
            maps = collections.defaultdict(dict)  # query copy -> {anchor offset: (query forward offset, rel strand)}
            for r in sorted(run_paf(["minimap2", "-cx", "asm20", "--cs", "-t", "1", tq, qq]), key=lambda r: -r["AS"]):
                if r["tp"] == "P":
                    for qo, t in anchor_pairs(r):
                        maps[int(r["qname"])].setdefault(t, (qo, r["strand"]))
            aln[f] = maps
            return aln[f]

        for j in range(ncol):
            if not decisive(crs, j):
                continue
            pos = cols[fam].get(j, -1)
            have = [i for i, c in enumerate(crs) if j < len(c["alleles"]) and c["alleles"][j] != "."]
            # anchor: the first copy carrying the column whose exons contain the dumped position and whose base matches
            anchor = None
            for i in have:
                x = info[i]
                if x is None or pos < 0:
                    continue
                if tx_offset(x["exons"], x["strand"], pos) is None:
                    continue
                b = G.fetch(x["chrom"], pos, pos + 1)
                # the dump's strand column is the allele frame (the frame strand since the orientation fix)
                if (comp(b) if crs[i]["strand"] == "-" else b) == crs[i]["alleles"][j]:
                    anchor = i
                    break
            for i in have:
                c, x = crs[i], info[i]
                row = {"family_id": fam, "col": j, "copy": i, "tid": c["copy_tid"], "allele": c["alleles"][j],
                       "chrom": "", "gpos": -1, "strand": "", "status": "", "rel": ""}
                if x is None or not x["seq"]:
                    row["status"] = "no_copy_row"
                elif anchor is None:
                    row["status"] = "no_anchor"
                else:
                    if i == anchor:
                        g, row["rel"] = pos, "anchor"
                    else:
                        fo = tx_offset(info[anchor]["exons"], info[anchor]["strand"], pos)
                        hit = align_to(anchor).get(i, {}).get(fo)
                        row["rel"] = hit[1] if hit else ""
                        g = None if hit is None else genome_pos(x["exons"], x["strand"], hit[0])
                    if g is None:
                        row["status"] = "unmapped"
                    else:
                        b = G.fetch(x["chrom"], g, g + 1)
                        ok = (comp(b) if c["strand"] == "-" else b) == c["alleles"][j]
                        row.update(chrom=x["chrom"], gpos=g, strand=x["strand"], status="verified" if ok else "mismatch")
                stat[row["status"]] += 1
                by_rel[f"{row['status']}|{row['rel'] or 'none'}"] += 1
                rows.append(row)
    write_tsv(f"{RES}/pos.tsv", rows, ["family_id", "col", "copy", "tid", "allele", "chrom", "gpos", "strand", "status", "rel"])
    tot = sum(stat.values())
    res = {"pairs": tot, **stat, "verified_frac": round(stat["verified"] / max(tot, 1), 4), "by_rel": dict(by_rel)}
    res["V_pass"] = res["verified_frac"] >= 0.95
    json.dump(res, open(f"{RES}/gate_v.json", "w"), indent=1)
    print(json.dumps(res))


def homologous():
    """Amendment 1: decisive columns where every copy carrying a non-gap allele has a verified position."""
    ok = collections.defaultdict(lambda: True)
    for r in tsv(f"{RES}/pos.tsv"):
        key = (r["family_id"], int(r["col"]))
        ok[key] = ok[key] and r["status"] == "verified"
    return {k for k, v in ok.items() if v}


def verified():
    hom = homologous()
    return [r for r in tsv(f"{RES}/pos.tsv") if r["status"] == "verified" and (r["family_id"], int(r["col"])) in hom]


def private_pairs():
    """(family, col, copy) whose allele no other copy of the family carries at that decisive column."""
    _, copies = dump()
    out = set()
    for fam, crs in copies.items():
        for j in range(max(len(c["alleles"]) for c in crs)):
            al = collections.Counter(c["alleles"][j] for c in crs if j < len(c["alleles"]) and c["alleles"][j] != ".")
            if len(al) < 2:
                continue
            for i, c in enumerate(crs):
                if j < len(c["alleles"]) and c["alleles"][j] != "." and al[c["alleles"][j]] == 1:
                    out.add((fam, j, i))
    return out


def cmd_kmers():
    G, priv = Fasta(GENOME), private_pairs()
    rows, allk = [], set()
    for r in verified():
        key = (r["family_id"], int(r["col"]), int(r["copy"]))
        if key not in priv:
            continue
        g = int(r["gpos"])
        w = G.fetch(r["chrom"], g - (K - 1), g + K)
        ks = [canon(w[i:i + K]) for i in range(len(w) - K + 1) if set(w[i:i + K]) <= set("ACGT")]
        allk.update(ks)
        rows.append({"family_id": key[0], "col": key[1], "copy": key[2], "chrom": r["chrom"], "gpos": g, "kmers": ",".join(ks) or "-"})
    write_tsv(f"{RES}/private.tsv", rows, ["family_id", "col", "copy", "chrom", "gpos", "kmers"])
    with open(f"{RES}/kmers.txt", "w") as fh:
        fh.write("\n".join(sorted(allk)) + "\n")
    print(f"{len(rows)} private verified alleles, {len(allk)} distinct candidate {K}-mers -> {RES}/kmers.txt")


def cmd_m1():
    counts = {}
    for line in open(f"{RES}/kmers.counts"):
        k, n = line.split()
        counts[k] = int(n)
    _, copies = dump()
    sun, fam_stat = set(), collections.defaultdict(collections.Counter)
    for r in tsv(f"{RES}/private.tsv"):
        ks = [] if r["kmers"] == "-" else r["kmers"].split(",")
        f = fam_stat[r["family_id"]]
        if not ks:
            f["private_no_kmer"] += 1
            continue
        f["private"] += 1
        if any(counts.get(k, 0) == 1 for k in ks):
            f["sun"] += 1
            sun.add((r["family_id"], int(r["col"]), int(r["copy"])))
        else:
            f["private_not_genome_unique"] += 1
    hom = homologous()
    for fam, crs in copies.items():
        for j in range(max(len(c["alleles"]) for c in crs)):
            if decisive(crs, j):
                fam_stat[fam]["decisive_cols_all"] += 1
                if (fam, j) not in hom:
                    continue
                fam_stat[fam]["decisive_cols"] += 1
                al = collections.Counter(c["alleles"][j] for c in crs if j < len(c["alleles"]) and c["alleles"][j] != ".")
                if any(v == 1 for v in al.values()):
                    fam_stat[fam]["cols_with_private"] += 1
                if any((fam, j, i) in sun for i in range(len(crs))):
                    fam_stat[fam]["cols_with_sun"] += 1
    # R: O2-assigned rows carrying >= 1 SUN allele of their assigned copy
    alle = {(fam, i): c["alleles"] for fam, crs in copies.items() for i, c in enumerate(crs)}
    for r in tsv(RUN + ".psv_reads.tsv"):
        if r["status"] != "assigned":
            continue
        fam, x = r["family_id"], int(r["assigned_copy"])
        ca = alle.get((fam, x), "")
        carries = any(b != "." and j < len(ca) and b == ca[j] and (fam, j, x) in sun for j, b in enumerate(r["alleles"]))
        fam_stat[fam]["assigned_rows"] += 1
        fam_stat[fam]["assigned_rows_with_sun"] += int(carries)
    keys = ["decisive_cols_all", "decisive_cols", "cols_with_private", "cols_with_sun", "private", "sun",
            "private_not_genome_unique", "private_no_kmer", "assigned_rows", "assigned_rows_with_sun"]
    rows = [{"family_id": f, **{k: s[k] for k in keys}} for f, s in sorted(fam_stat.items())]
    write_tsv(f"{RES}/m1_families.tsv", rows, ["family_id"] + keys)
    tot = {k: sum(r[k] for r in rows) for k in keys}
    tot["homologous_of_decisive"] = round(tot["decisive_cols"] / max(tot["decisive_cols_all"], 1), 4)
    rel = collections.defaultdict(set)
    for r in tsv(f"{RES}/pos.tsv"):
        rel[r["family_id"]].add(r["rel"])
    for name, sel in (("strand_consistent", lambda f: "-" not in rel[f]), ("has_reverse_pairs", lambda f: "-" in rel[f])):
        a = sum(r["decisive_cols_all"] for r in rows if sel(r["family_id"]))
        h = sum(r["decisive_cols"] for r in rows if sel(r["family_id"]))
        tot[f"homologous_of_decisive_{name}"] = f"{h}/{a} = {round(h / max(a, 1), 4)}"
    tot["G"] = round(tot["private_not_genome_unique"] / max(tot["private"], 1), 4)
    tot["sun_of_private"] = round(tot["sun"] / max(tot["private"], 1), 4)
    tot["cols_with_sun_of_decisive"] = round(tot["cols_with_sun"] / max(tot["decisive_cols"], 1), 4)
    tot["R"] = round(tot["assigned_rows_with_sun"] / max(tot["assigned_rows"], 1), 4)
    G_fam = sorted(round(r["private_not_genome_unique"] / r["private"], 3) for r in rows if r["private"] >= 5)
    tot["G_per_family_median"] = G_fam[len(G_fam) // 2] if G_fam else None
    tot["families_with_ge5_private"] = len(G_fam)
    json.dump(tot, open(f"{RES}/m1.json", "w"), indent=1)
    print(json.dumps(tot, indent=1))


def chrom_pairs(chrom):
    return [r for r in verified() if r["chrom"] == chrom]


def cmd_pileup(chrom):
    """samtools mpileup (MAPQ>=20 primaries) at this chromosome's verified positions, one call per copy window."""
    pairs = chrom_pairs(chrom)
    bytid = collections.defaultdict(set)
    for r in pairs:
        bytid[r["tid"]].add(int(r["gpos"]))
    out = []
    tmp = f"{OUT}/tmp/pos.{chrom}.tsv"
    os.makedirs(os.path.dirname(tmp), exist_ok=True)
    seen = {}
    for tid, ps in bytid.items():
        lo, hi = min(ps), max(ps)
        with open(tmp, "w") as fh:
            fh.write("".join(f"{chrom}\t{p + 1}\n" for p in sorted(ps)))
        res = subprocess.run(["samtools", "mpileup", "-r", f"{chrom}:{lo + 1}-{hi + 1}", "-l", tmp, "-q", "20", "-Q", "0",
                              "--ff", "2308", "-B", "-d", "0", "-f", GENOME, BAM], capture_output=True, text=True, check=True).stdout
        for line in res.splitlines():
            f = line.split("\t")
            p, ref, bases = int(f[1]) - 1, f[2].upper(), f[4]
            if p in seen:
                continue
            seen[p] = 1
            c = collections.Counter()
            i = 0
            while i < len(bases):
                ch = bases[i]
                if ch == "^":
                    i += 2
                    continue
                if ch in "+-":
                    j = i + 1
                    while j < len(bases) and bases[j].isdigit():
                        j += 1
                    i = j + int(bases[i + 1:j])
                    continue
                if ch in ".,":
                    c[ref] += 1
                elif ch.upper() in "ACGT":
                    c[ch.upper()] += 1
                i += 1
            out.append({"chrom": chrom, "gpos": p, "ref": ref, **{b: c[b] for b in "ACGT"}})
    write_tsv(f"{RES}/pileup/{chrom}.tsv", out, ["chrom", "gpos", "ref", "A", "C", "G", "T"])
    print(f"{chrom}: {len(out)} positions piled up ({len(bytid)} copy windows)")


def cmd_flags():
    pile = {}
    for p in glob.glob(f"{RES}/pileup/*.tsv"):
        for r in tsv(p):
            pile[(r["chrom"], int(r["gpos"]))] = {b: int(r[b]) for b in "ACGT"}
    _, copies = dump()
    o2 = collections.defaultdict(collections.Counter)  # (fam, copy, col) -> allele counts (genome orientation later)
    for r in tsv(RUN + ".psv_reads.tsv"):
        if r["status"] == "assigned":
            for j, b in enumerate(r["alleles"]):
                if b != ".":
                    o2[(r["family_id"], int(r["assigned_copy"]), j)][b] += 1
    rows = []
    for r in verified():
        fam, j, x, g = r["family_id"], int(r["col"]), int(r["copy"]), int(r["gpos"])
        minus = r["strand"] == "-"
        ref_g = comp(r["allele"]) if minus else r["allele"]
        c = collections.Counter(pile.get((r["chrom"], g), {}))
        n_unique = sum(c.values())
        for b, n in o2.get((fam, x, j), {}).items():
            c[comp(b) if minus else b] += n
        n = sum(c.values())
        alts = [(v, b) for b, v in c.items() if b != ref_g]
        n_alt, alt = max(alts) if alts else (0, "")
        cls = ""
        if n >= MIN_READS and n_alt >= MIN_ALT and n_alt / n >= MIN_ALT_FRAC:
            cls = "flag"
        elif n >= MIN_READS and c[ref_g] / n >= CONTROL_REF:
            cls = "control"
        rows.append({"family_id": fam, "col": j, "copy": x, "tid": r["tid"], "chrom": r["chrom"], "gpos": g,
                     "strand": r["strand"], "ref_g": ref_g, "alt_g": alt or "-", "n": n, "n_unique": n_unique,
                     "n_o2": n - n_unique, "n_ref": c[ref_g], "n_alt": n_alt, "cls": cls or "-"})
    write_tsv(f"{RES}/m2_pairs.tsv", rows, ["family_id", "col", "copy", "tid", "chrom", "gpos", "strand", "ref_g", "alt_g",
                                            "n", "n_unique", "n_o2", "n_ref", "n_alt", "cls"])
    cnt = collections.Counter(r["cls"] for r in rows)
    print(json.dumps({"pairs": len(rows), "covered_ge10": sum(r["n"] >= MIN_READS for r in rows), **cnt}))


def chrmap():
    return {r["pri"]: r for r in tsv(CHRMAP)}


def cmd_truth(chrom):
    cm = chrmap().get(chrom)
    pairs = [r for r in tsv(f"{RES}/m2_pairs.tsv") if r["chrom"] == chrom and r["cls"] in ("flag", "control")]
    if not pairs or not cm or not cm["B_name"]:
        print(f"{chrom}: nothing to do (pairs {len(pairs)}, B {cm and cm['B_name']!r})")
        return
    G, tab = Fasta(GENOME), copy_table()
    tmp = f"{OUT}/tmp"
    bfa, qfa = f"{tmp}/B.{chrom}.fa", f"{tmp}/q.{chrom}.fa"
    with open(bfa, "w") as fh:
        subprocess.run(["samtools", "faidx", f"{HAPS}/{cm['B_hap']}.fa", cm["B_name"]], stdout=fh, check=True)
    subprocess.run(["samtools", "faidx", bfa], check=True)
    q0 = {}
    with open(qfa, "w") as fh:
        for tid in sorted({r["tid"] for r in pairs}):
            t = tab[tid]
            s = max(0, int(t["start"]) - FLANK)
            e = min(G.length(chrom), int(t["end"]) + FLANK)
            q0[tid] = s
            fh.write(f">{tid}\n{G.fetch(chrom, s, e)}\n")
    recs = collections.defaultdict(list)
    for r in run_paf(["minimap2", "-cx", "asm5", "--cs", "-t", "2", bfa, qfa]):
        if r["qe"] - r["qs"] >= 0.95 * r["qlen"]:
            recs[r["qname"]].append(r)
    B = Fasta(bfa)
    out, qmaps = [], {}
    for r in pairs:
        hits = recs.get(r["tid"], [])
        prim = [h for h in hits if h["tp"] == "P"]
        row = dict(r, b_base="-", truth="NO_COUNTERPART", n_hits=len(hits))
        if len(prim) == 1:
            h = prim[0]
            if r["tid"] not in qmaps:
                qmaps[r["tid"]] = paf_qmap(h)
            t = qmaps[r["tid"]].get(int(r["gpos"]) - q0[r["tid"]])
            if t is not None:
                b = B.fetch(h["tname"], t, t + 1)
                b = comp(b) if h["strand"] == "-" else b
                row["b_base"] = b
                row["truth"] = ("B_ALT" if b == r["alt_g"] else "B_SAME" if b == r["ref_g"] else "B_OTHER")
            else:
                row["truth"] = "UNALIGNED_BASE"
        out.append(row)
    cols = list(pairs[0].keys()) + ["b_base", "truth", "n_hits"]
    write_tsv(f"{RES}/truth/{chrom}.tsv", out, cols)
    for p in (bfa, bfa + ".fai", qfa):
        os.remove(p)
    print(f"{chrom}: {collections.Counter(r['truth'] + '/' + r['cls'] for r in out)}")


def cmd_m2():
    rows = tsv(f"{RES}/m2_pairs.tsv")
    truth = [r for p in glob.glob(f"{RES}/truth/*.tsv") for r in tsv(p)]
    cm = chrmap()
    res = {}
    for cls in ("flag", "control"):
        t = collections.Counter(r["truth"] for r in truth if r["cls"] == cls)
        withc = t["B_ALT"] + t["B_OTHER"] + t["B_SAME"]
        res[cls] = {**t, "with_counterpart": withc,
                    "B_ALT_frac": round(t["B_ALT"] / max(withc, 1), 4),
                    "B_differs_frac": round((t["B_ALT"] + t["B_OTHER"]) / max(withc, 1), 4)}
    f, c = res["flag"], res["control"]
    res["bar_pass"] = f["B_ALT_frac"] >= 0.5 and f["B_differs_frac"] >= 5 * max(c["B_differs_frac"], 1e-9)
    res["fold"] = round(f["B_differs_frac"] / c["B_differs_frac"], 2) if c["B_differs_frac"] else None
    hap = lambda ch: bool(cm.get(ch, {}).get("B_name"))
    cov = [r for r in rows if int(r["n"]) >= MIN_READS]
    for name, sel in (("autosome_like", lambda r: hap(r["chrom"])), ("no_second_haplotype", lambda r: not hap(r["chrom"]))):
        s = [r for r in cov if sel(r)]
        res[name] = {"covered_pairs": len(s), "flags": sum(r["cls"] == "flag" for r in s),
                     "flag_rate": round(sum(r["cls"] == "flag" for r in s) / max(len(s), 1), 4),
                     "chroms": sorted({r["chrom"] for r in s})}
    ts = {("A", "G"), ("G", "A"), ("C", "T"), ("T", "C")}
    fl = [r for r in rows if r["cls"] == "flag"]
    res["flag_transitions_frac"] = round(sum((r["ref_g"], r["alt_g"]) in ts for r in fl) / max(len(fl), 1), 4)
    res["flags_with_o2_reads"] = sum(int(r["n_o2"]) > 0 for r in fl)
    json.dump(res, open(f"{RES}/m2.json", "w"), indent=1)
    print(json.dumps(res, indent=1))


if __name__ == "__main__":
    cmd, args = sys.argv[1], sys.argv[2:]
    {"prep": cmd_prep, "sanity": cmd_sanity, "positions": cmd_positions, "kmers": cmd_kmers, "m1": cmd_m1,
     "pileup": cmd_pileup, "flags": cmd_flags, "truth": cmd_truth, "m2": cmd_m2}[cmd](*args)
