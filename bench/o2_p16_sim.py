# P16 (LEDGER_O2): draw 40 fresh GGO_gwcat families, simulate 30 reads per copy, write the copies table and the truth
import collections, csv, random
C = "/mnt/linuxdisk/home/juanfraitu/winloci_data/GGO_gwcat.copies.tsv"
F = "/mnt/linuxdisk/home/juanfraitu/winloci_data/GGO_gwcat.copies.fa"
used = set(open("/mnt/linuxdisk/tmp/p15/used_fams.txt").read().split())
rows = collections.defaultdict(list)
for r in csv.DictReader(open(C), delimiter="\t"):
    rows[r["family_id"]].append(r)
seqs, cur = {}, None
for l in open(F):
    if l[0] == ">": cur = l[1:].strip(); seqs[cur] = []
    else: seqs[cur].append(l.strip())
seqs = {tuple(k.split("|")[:2]): "".join(v).upper() for k, v in seqs.items()}
elig = sorted(f for f, rs in rows.items() if 2 <= len(rs) <= 8 and f not in used and any(int(r["n_exon"]) >= 2 for r in rs)
              and all((f, r["copy_idx"]) in seqs for r in rs))
key = lambda f: (f.rstrip("0123456789"), int(f[len(f.rstrip("0123456789")):]))
elig.sort(key=key)
rng = random.Random(20261025)
panel = sorted(rng.sample(elig, 40), key=key)
print(f"eligible {len(elig)}; panel {len(panel)}: {','.join(panel)}")
B = "ACGT"
with open("copies.tsv", "w") as ct, open("reads.fa", "w") as fa, open("truth.tsv", "w") as tr:
    ct.write("family_id\tcopy_idx\ttid\tchrom\tstart\tend\tn_exon\tstrand\texons\n")
    tr.write("read\tfamily\tcopy\tchrom\tstart\tend\n")
    k = 0
    for f in panel:
        for r in sorted(rows[f], key=lambda r: int(r["copy_idx"])):
            tid = f"{f}:{r['copy_idx']}"
            ct.write(f"{f}\t{r['copy_idx']}\t{tid}\t{r['chrom']}\t{r['start']}\t{r['end']}\t{r['n_exon']}\t{r['strand']}\t{r['start']}-{r['end']}\n")
            t = seqs[(f, r["copy_idx"])]
            for _ in range(30):
                o = []
                for ch in t:
                    x = rng.random()
                    if x < 0.0005: o.append(rng.choice([b for b in B if b != ch]))
                    elif x < 0.00075: o.append(ch); o.append(rng.choice(B))
                    elif x < 0.001: pass
                    else: o.append(ch)
                n = f"{tid}|r{k:05d}"; k += 1
                fa.write(f">{n}\n{''.join(o)}\n")
                tr.write(f"{n}\t{f}\t{tid}\t{r['chrom']}\t{r['start']}\t{r['end']}\n")
print(f"copies {sum(len(rows[f]) for f in panel)}; reads {k}")
