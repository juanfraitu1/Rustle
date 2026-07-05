#!/usr/bin/env python3
"""vg_family_prototype_domain_share.py — split VG families by whole-protein homology.

For each family that contains multiple annotated protein families (E_p-impure),
extract the protein sequences of its member genes, run an all-vs-all mmseqs2 search,
and split the family into sub-families where every cross-gene pair passes a
whole-protein reciprocal-coverage bar (similar to the shipped O1 catalog's anti-
sub-domain safeguard).

Run: PYTHONHASHSEED=0 /home/juanfra/miniforge3/bin/python bench/vg_family_prototype_domain_share.py
"""
import os
import sys

if os.environ.get("PYTHONHASHSEED") != "0":
    os.environ["PYTHONHASHSEED"] = "0"
    os.execv(sys.executable, [sys.executable] + sys.argv)

import csv
import json
import subprocess
import tempfile
from collections import defaultdict

sys.path.insert(0, os.path.dirname(__file__))
from family_level_pr_current import build_ctx

BENCH = os.path.dirname(os.path.abspath(__file__))
PROTOTSV = os.path.join(BENCH, "vg_family_prototype.tsv")
PROTEINS = "/home/juanfra/winloci_scratch/proteins.fa"
MMSEQS = "/home/juanfra/miniforge3/bin/mmseqs"
OUT_TSV = os.path.join(BENCH, "vg_family_prototype_protcov.tsv")
OUT_JSON = os.path.join(BENCH, "vg_family_prototype_protcov.json")

# whole-protein reciprocal-coverage bars (shipped O1 anti-sub-domain safeguard)
MIN_COV = 0.50
MAX_COV = 0.70
FIDENT = 0.30
EVAL = 1e-5


def load_catalog(path):
    fams = defaultdict(list)
    with open(path) as fh:
        r = csv.DictReader(fh, delimiter="\t")
        for row in r:
            fams[int(row["fam_id"])].append(row["member"])
    return fams


def load_protein_lengths(fa_path):
    lens = {}
    name = None
    seq = []
    with open(fa_path) as fh:
        for line in fh:
            if line.startswith(">"):
                if name is not None:
                    lens[name] = len("".join(seq))
                name = line[1:].split()[0]
                seq = []
            else:
                seq.append(line.strip())
        if name is not None:
            lens[name] = len("".join(seq))
    return lens


def split_by_protein_cov(members, gene_of_dn, g2f, mega, prot_lens, min_cov=MIN_COV, max_cov=MAX_COV):
    """Return list of sub-families (each >=2 members) after splitting by prot coverage.

    If there are no proteins or no pairs pass the bar, return [] (family dissolved).
    """
    genes = [gene_of_dn.get(m) for m in members]
    genes = [g for g in genes if g and g in prot_lens]
    if len(genes) < 2:
        return []

    with tempfile.TemporaryDirectory(prefix="vgfp_protcov_") as td:
        in_fa = os.path.join(td, "prots.fa")
        out_pref = os.path.join(td, "out")
        tmp = os.path.join(td, "tmp")
        # write subset fasta
        seqs = {}
        name = None
        seq = []
        with open(PROTEINS) as fh:
            for line in fh:
                if line.startswith(">"):
                    if name in genes:
                        seqs[name] = "".join(seq)
                    name = line[1:].split()[0]
                    seq = []
                else:
                    seq.append(line.strip())
            if name in genes:
                seqs[name] = "".join(seq)
        with open(in_fa, "w") as fh:
            for g in genes:
                if g in seqs:
                    fh.write(f">{g}\n{seqs[g]}\n")
        # mmseqs2 search
        cmd = [
            MMSEQS, "easy-search", in_fa, in_fa, out_pref + ".m8", tmp,
            "--min-seq-id", str(FIDENT),
            "-e", str(EVAL),
            "--cov-mode", "0",
            "-c", str(min_cov),
            "--threads", "4",
            "-v", "0",
        ]
        subprocess.run(cmd, check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)

        # parse self-vs-self, keep reciprocal good pairs
        good_pairs = set()
        with open(out_pref + ".m8") as fh:
            for line in fh:
                f = line.rstrip("\n").split("\t")
                if len(f) < 12:
                    continue
                q, t, ident, alen, mm, gap, qstart, qend, tstart, tend, evalue, bits = f
                if q == t:
                    continue
                qstart, qend, tstart, tend = int(qstart), int(qend), int(tstart), int(tend)
                qlen = prot_lens[q]
                tlen = prot_lens[t]
                qcov = (qend - qstart + 1) / qlen
                tcov = (tend - tstart + 1) / tlen
                if float(ident) >= FIDENT * 100 and min(qcov, tcov) >= min_cov and max(qcov, tcov) >= max_cov:
                    good_pairs.add(tuple(sorted([q, t])))

    # connected components on good pairs
    adj = {g: set() for g in genes}
    for a, b in good_pairs:
        adj[a].add(b)
        adj[b].add(a)
    seen = set()
    comps = []
    for g in genes:
        if g in seen:
            continue
        stack = [g]
        comp = []
        while stack:
            cur = stack.pop()
            if cur in seen:
                continue
            seen.add(cur)
            comp.append(cur)
            for nb in adj[cur]:
                if nb not in seen:
                    stack.append(nb)
        if len(comp) >= 2:
            # map genes back to dn members (a gene may map to multiple dn loci, take all)
            gene_to_members = defaultdict(list)
            for m in members:
                gg = gene_of_dn.get(m)
                if gg:
                    gene_to_members[gg].append(m)
            sub = []
            for gg in comp:
                sub.extend(gene_to_members.get(gg, []))
            if len(sub) >= 2:
                comps.append(sub)
    return comps


def main():
    print("[*] building ctx ...", flush=True)
    ctx, *_ = build_ctx()
    gene_of_dn = ctx["gene_of_dn"]
    g2f = ctx["g2f"]
    mega = ctx["mega"]

    print("[*] loading catalog + protein lengths ...", flush=True)
    catalog = load_catalog(PROTOTSV)
    prot_lens = load_protein_lengths(PROTEINS)

    # identify families with multiple protein families
    def n_prot_fams(members):
        fams = set()
        for m in members:
            g = gene_of_dn.get(m)
            if g and g in g2f and g2f[g] not in mega:
                fams.add(g2f[g])
        return len(fams)

    split_candidates = {fid: members for fid, members in catalog.items() if n_prot_fams(members) >= 2}
    print(f"    families with >=2 protein families: {len(split_candidates)}", flush=True)

    # build new catalog
    new_families = []
    n_split = 0
    n_dissolved = 0
    for fid, members in sorted(catalog.items()):
        if fid in split_candidates:
            subs = split_by_protein_cov(members, gene_of_dn, g2f, mega, prot_lens)
            if not subs:
                n_dissolved += 1
            else:
                n_split += 1
                new_families.extend(subs)
        else:
            new_families.append(members)

    print(f"    split {n_split} families, dissolved {n_dissolved}", flush=True)
    print(f"    new total families: {len(new_families)}", flush=True)

    # write TSV
    with open(OUT_TSV, "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["fam_id", "member"])
        for fid, members in enumerate(new_families):
            for m in members:
                w.writerow([fid, m])
    print(f"[write] {OUT_TSV}", flush=True)

    # summary
    summary = dict(
        original_families=len(catalog),
        split_candidates=len(split_candidates),
        n_split=n_split,
        n_dissolved=n_dissolved,
        new_families=len(new_families),
    )
    with open(OUT_JSON, "w") as fh:
        json.dump(summary, fh, indent=2, sort_keys=True)
    print(f"[write] {OUT_JSON}", flush=True)


if __name__ == "__main__":
    main()
