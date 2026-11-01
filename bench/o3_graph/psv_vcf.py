#!/usr/bin/env python3
"""PREREG PSV VCF REPORT (register 1484): O2's pileup and variant calls written out for checking by hand, and compared with
bcftools on the same reads.

  psv_vcf.py run --bam BAM --copies COPIES.tsv --net NET.tsv --genome G.fa --o2 ARM.o2.tsv --out PREFIX
                 [--families A,B] [--no-bcftools] [--threads N] [--mapq-reset Q]

O2 (copy_phasing / o3port.py, `--semi --guard`) is a pileup + variant-calling step per copy of a family: the records in the tie
window on a copy are its pile (`pile_records`), their aligned bases per reference column the pileup (`observations`), and a column
is VARIABLE when >= 2 alleles have >= MIN_READS (2) carriers at base quality >= MIN_Q (10) and every minor allele beats a Poisson
error test (lambda = sum of the carriers' per-read error / 3; Bonferroni ALPHA (0.01) / m over the columns with two such alleles;
`informative`). A tied read is placed (`semi_tied`) by the distinguishing sites it covers on both of its two best copies: reference
positions where the copies' bases differ, counted only when the read's base there is carried by >= MIN_SUPPORT (6) reads of the
leading copy's pile. This script rebuilds those piles with o3port.py's own functions and writes:
  PREFIX.vcf        O2's variable columns, one record per copy column: ID family:copy:POS (POS 1-based), REF the genome base,
                    ALT the other supported alleles, QUAL -10 log10 of the smallest minor-allele p (999 at p = 0), INFO
                    FAM, COPY, DP (bases at Q >= 10), AD (per listed allele, REF first when supported), M (Bonferroni m), P.
  PREFIX.reads.tsv  per assigned read of the arm: family, read, copy, path (`sites` when semi_tied reproduces the arm's copy,
                    `haplotype` otherwise = the forced / certified path, `certified` for new:k copies), the rival copy, the
                    deciding-site IDs (sites: family:copy:POS=base of each site that voted; haplotype: the variable columns the
                    read carries a supported allele at), votes for / against.
  PREFIX.sites.tsv  every distinguishing site a `sites` read voted with: ID on the leading copy, its base, the rival copy's ID
                    and base, the leading pile's carriers of each base, how many reads it decided.
  PREFIX.bcftools/  per copy with >= 1 variable column: exactly the pile's records (one per read, written out) and `bcftools mpileup -f G -r SPAN
                    -Q 10 -q 0 --ff UNMAP -a AD,DP -d 100000 | bcftools call -mv` (secondaries kept: the piles hold them).
  PREFIX.compare.tsv per copy: O2 variable columns, those bcftools also calls as an SNV sharing an ALT, bcftools SNV calls O2 does
                    not hold variable with O2's reason (one_allele = fewer than two alleles with >= 2 carriers at Q >= 10; bonferroni =
                    p under ALPHA but not under ALPHA / m; p=... = the minor allele is explained by error; no_pile), bcftools indels.
  PREFIX.report.txt the totals.
Nothing here changes a placement: the arm's o2.tsv is the input, the piles are rebuilt the way the arm built them."""
import argparse
import collections
import math
import os
import subprocess
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import o3port as P  # noqa: E402


def good_bases(xs):
    return [x for x in xs if x[2] >= P.MIN_Q and x[1] in "ACGT"]


def column_test(xs):
    """(good bases, counts, supported alleles, major, {minor allele: Poisson p}) for one column, informative()'s arithmetic."""
    good = good_bases(xs)
    cnt = collections.Counter(x[1] for x in good)
    alleles = [b for b, n in cnt.items() if n >= P.MIN_READS]
    if len(alleles) < 2:
        return good, cnt, alleles, None, {}
    lam = sum(P.err(x) / 3 for x in good)
    major = max(alleles, key=lambda b: (cnt[b], b))
    return good, cnt, alleles, major, {b: P.poisson_tail(cnt[b], lam) for b in alleles if b != major}


def vcf_rows(chrom, fam, copy, obs, ref_of):
    """O2's variable columns of one pile as VCF-shaped dicts (pos 1-based), exactly the columns informative(obs) returns."""
    real = P.informative(obs)
    m = max(1, sum(1 for xs in obs.values() if len(column_test(xs)[2]) >= 2))
    rows = []
    for col in sorted(real):
        good, cnt, alleles, major, ps = column_test(obs[col])
        ref = ref_of(col)
        alleles = sorted(real[col], key=lambda b: (-cnt[b], b))
        alt = [b for b in alleles if b != ref]
        p = min(ps[b] for b in alleles if b != major)
        qual = 999.0 if p <= 0 else min(999.0, -10 * math.log10(p))
        rows.append({"chrom": chrom, "pos": col + 1, "id": f"{fam}:{copy}:{col + 1}", "ref": ref, "alt": alt, "qual": qual,
                     "dp": len(good), "ad": {b: cnt[b] for b in alleles}, "m": m, "p": p, "fam": fam, "copy": copy, "alleles": set(alleles)})
    return rows


def reason(xs, m=1):
    """Why O2 holds a column not variable, from its own pile bases: no_pile, one_allele (fewer than two alleles with >= MIN_READS
    carriers at Q >= MIN_Q), bonferroni (the minor allele's p is under ALPHA but not under ALPHA / m), p=... (explained by error)."""
    if not xs:
        return "no_pile"
    good, cnt, alleles, major, ps = column_test(xs)
    if len(alleles) < 2:
        return "one_allele"
    p = min(ps.values())
    if p < P.ALPHA / m:
        return "variable"
    return "bonferroni" if p < P.ALPHA else f"p={p:.2g}"


def compare(o2_cols, calls, obs_at, m=1):
    """o2_cols {(chrom, pos1): alt set}; calls {(chrom, pos1): (ref, [alts], 'snv'|'indel')}; obs_at(chrom, col0) -> pile bases;
    m the pile's Bonferroni m."""
    called = [k for k, alts in o2_cols.items() if k in calls and calls[k][2] == "snv" and set(calls[k][1]) & alts]
    missed = sorted(k for k in o2_cols if k not in called)
    only = sorted((k, reason(obs_at(k[0], k[1] - 1), m)) for k, (_, _, kind) in calls.items() if kind == "snv" and k not in o2_cols)
    return {"o2_cols": len(o2_cols), "o2_called": len(called), "o2_missed": missed, "caller_only": only,
            "caller_indels": sum(1 for v in calls.values() if v[2] == "indel")}


def write_vcf(fh, rows, contigs):
    fh.write("##fileformat=VCFv4.2\n##source=psv_vcf.py (O2 variable columns, register 1484)\n")
    for c, ln in contigs:
        fh.write(f"##contig=<ID={c},length={ln}>\n")
    fh.write('##INFO=<ID=FAM,Number=1,Type=String,Description="family">\n##INFO=<ID=COPY,Number=1,Type=String,Description="copy">\n'
             '##INFO=<ID=DP,Number=1,Type=Integer,Description="pile bases at Q >= 10">\n'
             '##INFO=<ID=AD,Number=.,Type=Integer,Description="carriers per listed allele (REF first when supported, then ALT)">\n'
             '##INFO=<ID=M,Number=1,Type=Integer,Description="Bonferroni m: columns of the pile with two alleles of >= 2 carriers">\n'
             '##INFO=<ID=P,Number=1,Type=Float,Description="smallest minor-allele Poisson tail">\n'
             "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
    for r in rows:
        order = ([r["ref"]] if r["ref"] in r["ad"] else []) + r["alt"]
        ad = ",".join(str(r["ad"][b]) for b in order)
        fh.write(f"{r['chrom']}\t{r['pos']}\t{r['id']}\t{r['ref']}\t{','.join(r['alt'])}\t{r['qual']:.1f}\t.\t"
                 f"FAM={r['fam']};COPY={r['copy']};DP={r['dp']};AD={ad};M={r['m']};P={r['p']:.3g}\n")


def deciding_sites(r, ca, cb, piles, prim, genome, cache):
    """site_votes with the sites: [(ta, tb, ra, rb, read base in ca's strand, its carriers in ca's pile, vote 'a' | 'b', carriers of
    ra, carriers of rb)]; the same filters as o3port.site_votes (homologous flank, >= MIN_SUPPORT carriers of the read's base)."""
    a, b = (P.site_bases(piles[c][r], prim[r], genome) for c in (ca, cb))
    col = P.column_bases(piles, ca, prim, genome, cache)
    out = []
    for i in sorted(a.keys() & b.keys()):
        (ra, qa, ta), (rb, _, tb) = a[i], b[i]
        if ra != rb and ra in "ACGT" and rb in "ACGT" and P.homologous(i, a, b):
            fwd = P.rec_strand_base(qa, piles[ca][r], prim[r])
            if col[ta][fwd] < P.MIN_SUPPORT:
                continue
            if qa == ra or qa == rb:
                out.append((ta, tb, ra, rb, fwd, col[ta][fwd], "a" if qa == ra else "b", col[ta][ra], col[ta][rb]))
    return out


def ranked_copies(r, piles):
    key = (lambda rec: rec.get_tag("ms")) if P.TIE_MS is not None else (lambda rec: rec.get_tag("AS"))
    return sorted(((key(rec), -rec.get_tag("NM"), c) for c, p in piles.items() if r in p for rec in [p[r]]), reverse=True)


def bcftools_calls(bam, recs, chrom, s, e, genome, outdir, tag, threads, mapq=None):
    """bcftools mpileup | call on exactly O2's pile records (one record per read, the one in the tie window on this copy; the
    piles hold secondaries, so --ff keeps them), over the copy span; {(chrom, pos1): (ref, alts, kind)}. With `mapq`, every
    record's MAPQ is reset to it first: bcftools caps a base's quality at its read's MAPQ, and the piles' secondaries carry
    MAPQ 0 (as run on OR6737: 24% of pile records), so a tied read's bases count for nothing in its genotype likelihoods."""
    import pysam
    os.makedirs(outdir, exist_ok=True)
    bf, vf = (os.path.join(outdir, f"{tag}.{x}") for x in ("bam", "vcf"))
    with pysam.AlignmentFile(bf, "wb", template=bam) as out:
        for rec in sorted(recs, key=lambda r: (r.reference_start, r.query_name)):
            if mapq is not None:
                rec.mapping_quality = mapq
            out.write(rec)
    pysam.index(bf)
    region = f"{chrom}:{s + 1}-{e}"
    with open(vf, "w") as fh:
        mp = subprocess.Popen(["bcftools", "mpileup", "-f", genome, "-r", region, "-Q", "10", "-q", "0", "--ff", "UNMAP", "-a", "AD,DP",
                               "-d", "100000", "--threads", str(threads), bf], stdout=subprocess.PIPE, stderr=subprocess.DEVNULL)
        subprocess.run(["bcftools", "call", "-mv"], stdin=mp.stdout, stdout=fh, check=True, stderr=subprocess.DEVNULL)
        mp.wait()
    return parse_vcf(vf)


def parse_vcf(path):
    calls = {}
    for line in open(path):
        if line.startswith("#"):
            continue
        f = line.rstrip("\n").split("\t")
        alts = f[4].split(",")
        kind = "snv" if len(f[3]) == 1 and all(len(x) == 1 for x in alts) else "indel"
        calls[(f[0], int(f[1]))] = (f[3], alts, kind)
    return calls


def cmd_run(a):
    import pysam
    assigned = collections.defaultdict(dict)
    for line in open(a.o2):
        fam, read, copy, n = line.rstrip("\n").split("\t")
        if fam == "family" or copy in ("ABSTAIN", "HOPELESS"):
            continue
        assigned[fam][read] = copy
    fams = set(a.families.split(",")) if a.families else set(assigned)
    fams &= set(assigned)
    by_fam = P.load_copies(a.copies)
    net, _ = P.load_net(a.net, fams)
    bam = pysam.AlignmentFile(a.bam)
    genome = pysam.FastaFile(a.genome)
    wanted = set().union(*net.values()) if net else set()
    print(f"[psv_vcf] {len(fams)} families, {sum(len(v) for f, v in assigned.items() if f in fams)} assigned reads, "
          f"{len(wanted)} net reads; reading primaries", file=sys.stderr, flush=True)
    prim = P.read_primaries(bam, wanted)
    rows_all, sites_used, paths = [], collections.Counter(), collections.Counter()
    cmp_rows, tot = [], collections.Counter()
    only_reasons = collections.Counter()
    rf = open(f"{a.out}.reads.tsv", "w")
    rf.write("family\tread\tcopy\tpath\trival\tvotes_for\tvotes_against\tsites\n")
    for fam in sorted(fams):
        copies = by_fam[fam]
        piles = P.pile_records(bam, copies, net[fam], prim)
        span = {cid: (ch, s, e) for cid, ch, s, e in copies}
        cache = {}
        rows_by_copy, obs_by_copy, read_cols = {}, {}, collections.defaultdict(dict)
        for cid, pile in piles.items():
            if len(pile) < 2 * P.MIN_READS:
                continue
            ch, s, e = span[cid]
            obs, _ = P.observations(pile, prim, (s, e))
            obs_by_copy[cid] = obs
            rows = vcf_rows(ch, fam, cid, obs, lambda col, ch=ch: genome.fetch(ch, col, col + 1).upper())
            rows_by_copy[cid] = rows
            rows_all += rows
            for r in rows:
                for x in obs[r["pos"] - 1]:
                    if x[2] >= P.MIN_Q and x[1] in r["alleles"]:
                        read_cols[x[0]][cid] = read_cols[x[0]].get(cid, []) + [f"{r['id']}={x[1]}"]
        for read, copy in sorted(assigned[fam].items()):
            if copy.startswith("new:"):
                paths["certified"] += 1
                rf.write(f"{fam}\t{read}\t{copy}\tcertified\t.\t.\t.\t.\n")
                continue
            sc = ranked_copies(read, piles)
            semi = P.semi_tied(read, piles, prim, genome, cache)
            rival = sc[1][2] if len(sc) >= 2 and sc[0][2] == copy else next((c for _, _, c in sc if c != copy), ".")
            if semi == copy and len(sc) >= 2:
                ds = deciding_sites(read, copy, rival, piles, prim, genome, cache)
                va, vb = sum(1 for d in ds if d[6] == "a"), sum(1 for d in ds if d[6] == "b")
                ids = ";".join(f"{fam}:{copy}:{d[0] + 1}={d[4]}" for d in ds)
                for ta, tb, ra, rb, fwd, n, v, na, nb in ds:
                    sites_used[(fam, copy, ta, ra, rival, tb, rb, na, nb)] += 1
                paths["sites"] += 1
                rf.write(f"{fam}\t{read}\t{copy}\tsites\t{rival}\t{va}\t{vb}\t{ids}\n")
            else:
                paths["haplotype"] += 1
                ids = ";".join(read_cols.get(read, {}).get(copy, []))
                rf.write(f"{fam}\t{read}\t{copy}\thaplotype\t{rival}\t.\t.\t{ids}\n")
        if not a.no_bcftools:
            for cid, rows in rows_by_copy.items():
                if not rows:
                    continue
                ch, s, e = span[cid]
                calls = bcftools_calls(bam, list(piles[cid].values()), ch, s, e, a.genome, f"{a.out}.bcftools", f"{fam}.{cid}", a.threads, a.mapq_reset)
                o2c = {(ch, r["pos"]): set(r["alt"]) for r in rows}
                obs = obs_by_copy[cid]
                out = compare(o2c, calls, lambda chrom, col, obs=obs: obs.get(col, []), rows[0]["m"])
                tot["o2_cols"] += out["o2_cols"]
                tot["o2_called"] += out["o2_called"]
                tot["only"] += len(out["caller_only"])
                tot["indels"] += out["caller_indels"]
                tot["snv_calls"] += sum(1 for v in calls.values() if v[2] == "snv")
                for _, why in out["caller_only"]:
                    only_reasons[why.split("=")[0]] += 1
                cmp_rows.append((fam, cid, len(piles[cid]), out["o2_cols"], out["o2_called"], len(out["caller_only"]), out["caller_indels"],
                                 ";".join(f"{k[1]}:{w}" for k, w in out["caller_only"]), ";".join(str(k[1]) for k in out["o2_missed"])))
        print(f"[psv_vcf] {fam}: piles {len(piles)}, variable columns {sum(len(r) for r in rows_by_copy.values())}, "
              f"assigned {len(assigned[fam])}", file=sys.stderr, flush=True)
    rf.close()
    contigs = sorted({r["chrom"] for r in rows_all})
    with open(f"{a.out}.vcf", "w") as fh:
        write_vcf(fh, sorted(rows_all, key=lambda r: (r["chrom"], r["pos"], r["id"])), [(c, genome.get_reference_length(c)) for c in contigs])
    with open(f"{a.out}.sites.tsv", "w") as fh:
        fh.write("family\tcopy\tid\tbase\trival\trival_id\trival_base\tpile_carriers_base\tpile_carriers_rival_base\treads_decided\n")
        for (fam, c, ta, ra, cb, tb, rb, na, nb), n in sorted(sites_used.items()):
            fh.write(f"{fam}\t{c}\t{fam}:{c}:{ta + 1}\t{ra}\t{cb}\t{fam}:{cb}:{tb + 1}\t{rb}\t{na}\t{nb}\t{n}\n")
    with open(f"{a.out}.compare.tsv", "w") as fh:
        fh.write("family\tcopy\tpile_reads\to2_cols\tbcftools_called\tbcftools_only_snv\tbcftools_indels\tonly_reasons\to2_missed_pos\n")
        for r in cmp_rows:
            fh.write("\t".join(map(str, r)) + "\n")
    with open(f"{a.out}.report.txt", "w") as fh:
        fh.write(f"families {len(fams)}\nassigned reads {sum(paths.values())}: " + ", ".join(f"{k} {v}" for k, v in sorted(paths.items())) + "\n")
        fh.write(f"distinct deciding sites {len(sites_used)}\nO2 variable columns {len(rows_all)} over {len(cmp_rows)} copies with bcftools\n")
        if tot["o2_cols"]:
            fh.write(f"bcftools calls an SNV sharing an ALT at {tot['o2_called']} of {tot['o2_cols']} O2 columns ({100 * tot['o2_called'] / tot['o2_cols']:.1f}%)\n")
            fh.write(f"bcftools SNV calls {tot['snv_calls']}, not variable for O2 {tot['only']}: " + ", ".join(f"{k} {v}" for k, v in sorted(only_reasons.items())) + f"; indel calls {tot['indels']}\n")
    print(open(f"{a.out}.report.txt").read(), file=sys.stderr)


def main(argv):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd", required=True)
    r = sub.add_parser("run")
    for k in ("bam", "copies", "net", "genome", "o2", "out"):
        r.add_argument(f"--{k}", required=True)
    r.add_argument("--families", default="")
    r.add_argument("--no-bcftools", action="store_true")
    r.add_argument("--threads", type=int, default=4)
    r.add_argument("--mapq-reset", type=int, default=None, metavar="Q", help="reset every pile record's MAPQ to Q before bcftools "
                   "(bcftools caps base quality at MAPQ; the piles' secondaries carry MAPQ 0)")
    a = ap.parse_args(argv)
    cmd_run(a)


if __name__ == "__main__":
    main(sys.argv[1:])
