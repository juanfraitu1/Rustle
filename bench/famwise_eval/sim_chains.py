"""PREREG r1797 / r1798 (docs/PREREG_famwise_sim_2026-11-18.md), row 1798: the chain table. Which simulated isoform
chains each assembly input recovers, the ceiling C, R, and false copies.

  sim_chains.py FAM.json ORIGIN.tsv EXACT_READS.tsv OUT_PREFIX TOOL:INPUT=GTF [TOOL:INPUT=GTF ...] [--drop-families ID,ID,...]
  e.g. rustle:P=rP/x.gtf rustle:D=... rustle:O=... rustle:N=... stringtie:P=...   (P is required for every tool)

Chain = (chrom, strand, introns), the key of bench/famcmp/chain_prec.chains; single-exon isoforms are in no count.
Reference chain: >= 2 drawn reads with full == 1, summed over every member / isoform / haplotype with that chain.
Affected: >= 1 of those reads is off its origin (status not decidable, primary not in origin_loci, or not in the table).
OUT_PREFIX.chains.tsv   one row per reference chain; members / families = every member of FAM.json with that chain.
OUT_PREFIX.summary.tsv  per tool and subset (all, depth classes, without the --drop-families chains): affected chains
                        recovered per input, C = affected / P, X_over_P (R is O_over_P), ref_X = all reference chains
                        recovered. NA = the tool has no such input, or P recovered none.
OUT_PREFIX.false.tsv    per silent member (no drawn read, on no emitting member's span and strand): assembled chains on
                        its strand whose introns lie inside its span; TOTAL = distinct chains."""
import argparse, collections, csv, json, os, re, sys
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "famcmp"))
from chain_prec import chains

SUBSETS = ("all", "2-4", "5-19", "20+")
DROPPED = "without dropped families"

def chain_key(copy, exons):
    """the key chains() gives this isoform written as a GTF; None for a single exon"""
    introns = tuple((a, b) for (_, a), (b, _) in zip(exons, exons[1:]) if a < b)    # abutting exons are one exon
    return (copy["chrom"], copy["strand"], introns) if introns else None

def depth(full_reads):
    return "2-4" if full_reads < 5 else "5-19" if full_reads < 20 else "20+"

def load_family(path):
    """-> [(family, copy)] of the one set"""
    with open(path) as fh:
        sets = json.load(fh)["sets"]
    if len(sets) != 1:
        sys.exit(f"{path}: one set expected, found {sorted(sets)}")
    return [(str(f["family"]), c) for f in next(iter(sets.values())) for c in f["copies"]]

def at_origin(exact_tsv):
    """-> the reads whose primary is at their exact origin"""
    out = set()
    with open(exact_tsv) as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            m = re.search(r"(?:^|;)primary=(\d+)", r["picks"])
            if r["status"] == "decidable" and m and m.group(1) in r["origin_loci"].split(","):
                out.add(r["read"])
    return out

def reference(fam, origin_tsv, home):
    """-> ({chain: [full reads, of them off their origin]} with >= 2 full reads, the members that emitted a read)"""
    key = {(c["gene_id"], t): chain_key(c, ex) for _, c in fam for t, ex in c["transcripts"].items()}
    n, emitting = collections.defaultdict(lambda: [0, 0]), set()
    with open(origin_tsv) as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            emitting.add(r["member"])
            if (r["member"], r["isoform"]) not in key:
                sys.exit(f"{origin_tsv}: read {r['read']} of {r['member']} / {r['isoform']}, which the family file lacks")
            k = key[r["member"], r["isoform"]]
            if k and r["full"] == "1":
                n[k][0] += 1; n[k][1] += r["read"] not in home
    return {k: v for k, v in n.items() if v[0] >= 2}, emitting

def parse_arms(specs):
    """-> [(tool, label, gtf)] in the order given"""
    arms = []
    for s in specs:
        name, _, gtf = s.partition("="); tool, _, label = name.partition(":")
        if not (tool and label and gtf) or (tool, label) in [a[:2] for a in arms]:
            sys.exit(f"bad or repeated TOOL:INPUT=GTF: {s}")
        arms.append((tool, label, gtf))
    for tool in dict.fromkeys(a[0] for a in arms):
        if (tool, "P") not in [a[:2] for a in arms]:
            sys.exit(f"{tool}: input P is required (C and every ratio are against it)")
    return arms

def summary(table, arms, drop):
    """table rows: (families, depth class, affected, {(tool, label) recovered}) -> header and rows of the summary"""
    labels = list(dict.fromkeys(a[1] for a in arms))
    rest = [x for x in labels if x != "P"]
    out = [["tool", "subset", "reference", "affected", *labels, "C", *(f"{x}_over_P" for x in rest), *(f"ref_{x}" for x in labels)]]
    for tool in dict.fromkeys(a[0] for a in arms):
        has = {a[1] for a in arms if a[0] == tool}
        for sub in SUBSETS + ((DROPPED,) if drop else ()):
            rows = [r for r in table if sub == "all" or sub == r[1] or (sub == DROPPED and not r[0] & drop)]
            ref = {x: sum((tool, x) in r[3] for r in rows) for x in has}
            aff = {x: sum((tool, x) in r[3] for r in rows if r[2]) for x in has}
            n, p = sum(r[2] for r in rows), aff["P"]
            out.append([tool, sub, len(rows), n, *(aff.get(x, "NA") for x in labels), f"{n / p:.2f}" if p else "inf",
                        *(f"{aff[x] / p:.2f}" if p and x in has else "NA" for x in rest), *(ref.get(x, "NA") for x in labels)])
    return out

def false_copies(fam, emitting, got):
    """got: {(tool, label): chains} -> rows (member, gene_name, family, a count per arm), then TOTAL of distinct chains"""
    on = lambda a, b: a["chrom"] == b["chrom"] and a["strand"] == b["strand"] and a["start"] < b["end"] and b["start"] < a["end"]
    loud = [c for _, c in fam if c["gene_id"] in emitting]
    silent = [(f, c) for f, c in fam if c["gene_id"] not in emitting and not any(on(c, e) for e in loud)]
    by_strand = {arm: collections.defaultdict(list) for arm in got}
    for arm, cs in got.items():
        for k in cs:
            by_strand[arm][k[:2]].append(k)
    hit = {arm: [{k for k in d[c["chrom"], c["strand"]] if c["start"] <= k[2][0][0] and k[2][-1][1] <= c["end"]} for _, c in silent]
           for arm, d in by_strand.items()}
    rows = [[c["gene_id"], c["gene_name"], f, *(len(hit[arm][i]) for arm in got)] for i, (f, c) in enumerate(silent)]
    return rows + [["TOTAL", "", "", *(len(set().union(*hit[arm])) for arm in got)]]

def write(path, rows):
    with open(path, "w") as fh:
        fh.writelines("\t".join(map(str, r)) + "\n" for r in rows)

def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("fam"); ap.add_argument("origin"); ap.add_argument("exact"); ap.add_argument("out")
    ap.add_argument("arms", nargs="+", metavar="TOOL:INPUT=GTF"); ap.add_argument("--drop-families", default="")
    a = ap.parse_args(argv)
    arms, drop = parse_arms(a.arms), set(filter(None, a.drop_families.split(",")))
    fam = load_family(a.fam)
    ref, emitting = reference(fam, a.origin, at_origin(a.exact))
    got = {(tool, label): chains(gtf) for tool, label, gtf in arms}
    owners = collections.defaultdict(set)
    for f, c in fam:
        for ex in c["transcripts"].values():
            owners[chain_key(c, ex)].add((c["gene_id"], f))
    names = [f"{tool}:{label}" for tool, label in got]
    rows, table = [["chrom", "strand", "first_intron_start", "last_intron_end", "introns", "members", "families", "full_reads", "off_reads",
                    "depth", "affected", *names]], []
    for k in sorted(ref):
        full, off = ref[k]
        fams, rec = {f for _, f in owners[k]}, {arm for arm in got if k in got[arm]}
        table.append((fams, depth(full), off > 0, rec))
        rows.append([k[0], k[1], k[2][0][0], k[2][-1][1], ",".join(f"{s}-{e}" for s, e in k[2]), ",".join(sorted(g for g, _ in owners[k])),
                     ",".join(sorted(fams)), full, off, depth(full), int(off > 0), *(int(arm in rec) for arm in got)])
    write(a.out + ".chains.tsv", rows)
    sm = summary(table, arms, drop)
    write(a.out + ".summary.tsv", sm)
    fc = [["member", "gene_name", "family", *names]] + false_copies(fam, emitting, got)
    write(a.out + ".false.tsv", fc)
    for r in sm + [fc[0], fc[-1], ["silent members", len(fc) - 2]]:
        print("\t".join(map(str, r)))

if __name__ == "__main__":
    main()
