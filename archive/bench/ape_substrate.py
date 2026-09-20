#!/usr/bin/env python3
"""Evidence-agnostic families (spec 2026-09-14 §2): ape substrates = target chromosomes orthologous to human chr7/15/16/17,
found by Liftoff placement of human RefSeq genes (batched foreground calls); substrate FASTA; native GFF subset and the
50% native subsample (development/report species only)."""
import argparse
import collections
import glob
import os
import random
import subprocess
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

HUMAN_FA = "/mnt/linuxdisk/home/juanfraitu/winloci_data/chm13v2.0.fa"
HUMAN_GFF = "/mnt/linuxdisk/home/juanfraitu/winloci_data/Reference/chm13v2.0_RefSeq_full.gff.gz"
LIFTOFF = "/home/juanfra/miniforge3/envs/liftoff/bin/liftoff"
HUMAN_CHROMS = ("chr7", "chr15", "chr16", "chr17")
TOP = ("gene", "pseudogene")


def attrs(col):
    return dict(kv.split("=", 1) for kv in col.strip().split(";") if "=" in kv)


def gff_batches(lines, max_genes):
    top_of, order, members = {}, [], collections.defaultdict(list)
    for line in lines:
        if line.startswith("#"):
            continue
        f = line.split("\t")
        if len(f) < 9:
            continue
        a = attrs(f[8])
        if f[2] in TOP and "ID" in a:
            top_of[a["ID"]] = a["ID"]
            order.append(a["ID"])
            members[a["ID"]].append(line)
            continue
        p = a.get("Parent", "").split(",")[0]
        top = top_of.get(p)
        if top is None:
            continue
        if "ID" in a:
            top_of[a["ID"]] = top
        members[top].append(line)
    batches, cur, n = [], [], 0
    for gid in order:
        cur.extend(members[gid])
        n += 1
        if n == max_genes:
            batches.append(cur)
            cur, n = [], 0
    if cur:
        batches.append(cur)
    return batches


def select_substrate(counts, min_frac=0.10):
    tot = collections.Counter()
    for (h, t), n in counts.items():
        tot[h] += n
    keep = {t for (h, t), n in counts.items() if tot[h] and n / tot[h] >= min_frac}
    return sorted(keep)


def cmd_lift(a):
    os.makedirs(a.out, exist_ok=True)
    bdir = f"{a.out}/batches"
    os.makedirs(bdir, exist_ok=True)
    mgfile = f"{bdir}/max_genes.txt"
    if glob.glob(f"{bdir}/*.gff"):
        if os.path.exists(mgfile):
            prev = open(mgfile).read().strip()
            if prev != str(a.max_genes):
                print(f"[lift] batches were generated with --max-genes {prev}, not {a.max_genes}; "
                      f"delete {bdir} and rerun to regenerate")
                sys.exit(1)
        else:
            open(mgfile, "w").write(str(a.max_genes))
    else:
        # (re)generating batches: reused batch indices could otherwise be mistaken for already-lifted
        for pat in ("*.lifted.gff3", "*.unmapped.txt", "*.liftoff.log"):
            for p in glob.glob(f"{bdir}/{pat}"):
                os.remove(p)
        for h in HUMAN_CHROMS:
            lines = subprocess.run(["tabix", HUMAN_GFF, h], capture_output=True, text=True, check=True).stdout.splitlines(True)
            for i, b in enumerate(gff_batches(lines, a.max_genes)):
                open(f"{bdir}/{h}.{i:03d}.gff", "w").writelines(b)
        open(mgfile, "w").write(str(a.max_genes))
    from dup_evidence import run_budget
    mmi = a.target + ".mmi"
    if not os.path.exists(mmi):
        mmi_log = f"{mmi}.log"
        try:
            ok = run_budget(["minimap2", "-d", mmi + ".tmp", a.target, "-t", str(a.threads)], a.budget, log=mmi_log)
        except subprocess.CalledProcessError as e:
            if os.path.exists(mmi + ".tmp"):
                os.remove(mmi + ".tmp")
            print(f"[lift] index build failed (exit {e.returncode}); see {mmi_log}")
            sys.exit(1)
        if not ok:
            if os.path.exists(mmi + ".tmp"):
                os.remove(mmi + ".tmp")
            print("[lift] target index build did not finish; rerun")
            return
        os.replace(mmi + ".tmp", mmi)
    todo = [p for p in sorted(glob.glob(f"{bdir}/*.gff")) if not os.path.exists(p.replace(".gff", ".lifted.gff3"))]
    if not todo:
        print("[lift] all batches lifted")
        return
    import shutil
    b = todo[0]
    out = b.replace(".gff", ".lifted.gff3")
    blog = b.replace(".gff", ".liftoff.log")
    try:
        ok = run_budget([LIFTOFF, "-g", b, "-o", out + ".tmp", "-u", b.replace(".gff", ".unmapped.txt"), "-dir",
                         f"{a.out}/liftoff_tmp", "-p", str(a.threads), "-m", shutil.which("minimap2"), a.target, HUMAN_FA],
                        a.budget, log=blog)
    except subprocess.CalledProcessError as e:
        if os.path.exists(out + ".tmp"):
            os.remove(out + ".tmp")
        print(f"[lift] {os.path.basename(b)} failed (exit {e.returncode}); see {blog}")
        sys.exit(1)
    if not ok:
        if os.path.exists(out + ".tmp"):
            os.remove(out + ".tmp")
        print(f"[lift] {os.path.basename(b)} did not finish in {a.budget} s; regenerate batches with a smaller --max-genes")
        return
    os.replace(out + ".tmp", out)
    print(f"[lift] {os.path.basename(b)} done; {len(todo) - 1} batches left")


def cmd_select(a):
    counts = collections.Counter()
    with open(f"{a.out}/lifted.gff", "w") as fo:
        for p in sorted(glob.glob(f"{a.out}/batches/*.lifted.gff3")):
            h = os.path.basename(p).split(".")[0]
            for line in open(p):
                if line.startswith("#"):
                    continue
                f = line.rstrip("\n").split("\t")
                if len(f) < 9:
                    continue
                if f[2] in TOP:
                    counts[(h, f[0])] += 1
                fo.write("\t".join(f[:8] + [f[8] + f";human_source={h}"]) + "\n")
    sub = select_substrate(counts)
    with open(f"{a.out}/lift_matrix.tsv", "w") as fh:
        fh.write("human_chrom\ttarget_chrom\tgenes\n")
        for (h, t), n in sorted(counts.items()):
            fh.write(f"{h}\t{t}\t{n}\n")
    open(f"{a.out}/substrate.txt", "w").write("\n".join(sub) + "\n")
    print(f"[select] substrate: {','.join(sub)}")


def cmd_subset(a):
    sub = open(f"{a.out}/substrate.txt").read().split()
    fa = f"{a.out}/substrate.fa"
    if not os.path.exists(fa):
        with open(fa, "w") as fh:
            subprocess.run(["samtools", "faidx", a.target] + sub, stdout=fh, check=True)
        subprocess.run(["samtools", "faidx", fa], check=True)
    if a.native_gff:
        import guided_min
        s = set(sub)
        with open(f"{a.out}/native.gff", "w") as fh:
            for line in open(a.native_gff):
                if not line.startswith("#") and line.split("\t", 1)[0] in s:
                    fh.write(line)
        genes, _ = guided_min.load_genes(f"{a.out}/native.gff", s)
        names = sorted(genes)
        keep = sorted(random.Random(1).sample(names, round(0.5 * len(names))))
        open(f"{a.out}/native_subsample50.txt", "w").write("\n".join(keep) + "\n")
        print(f"[subset] native genes {len(names)}; subsample {len(keep)}")
    print(f"[subset] {fa}")


def main():
    ap = argparse.ArgumentParser()
    sub = ap.add_subparsers(dest="cmd", required=True)
    p = sub.add_parser("lift")
    p.add_argument("--species", required=True)
    p.add_argument("--target", required=True)
    p.add_argument("--out", required=True)
    p.add_argument("--max-genes", type=int, default=1500)
    p.add_argument("--threads", type=int, default=4)
    p.add_argument("--budget", type=int, default=560)
    p = sub.add_parser("select")
    p.add_argument("--out", required=True)
    p = sub.add_parser("subset")
    p.add_argument("--target", required=True)
    p.add_argument("--out", required=True)
    p.add_argument("--native-gff")
    a = ap.parse_args()
    {"lift": cmd_lift, "select": cmd_select, "subset": cmd_subset}[a.cmd](a)


if __name__ == "__main__":
    main()
