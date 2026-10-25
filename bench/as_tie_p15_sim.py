# P15 (LEDGER_O2, register 1404): simulate reads from 300 gorilla genes drawn from the real window-A reads; run in /mnt/linuxdisk/tmp/p15 (needs a_primaries.tsv)
import collections, random, bisect, subprocess
GFF = "/mnt/linuxdisk/home/juanfraitu/winloci_data/GGO_genomic.gff"
FA = "/mnt/linuxdisk/home/juanfraitu/_from_wsl/winloci_scratch/GGO.fasta"
TX_TYPES = {"mRNA", "lnc_RNA", "transcript"}
attr = lambda s, k: next((x[len(k) + 1:] for x in s.split(";") if x.startswith(k + "=")), None)
genes, tx_of, tx, exons = {}, collections.defaultdict(list), {}, collections.defaultdict(list)
for ln in open(GFF):
    if ln[0] == "#": continue
    f = ln.rstrip("\n").split("\t")
    if len(f) < 9: continue
    if f[2] == "gene":
        genes[attr(f[8], "ID")] = (f[0], int(f[3]) - 1, int(f[4]))
    elif f[2] in TX_TYPES:
        i, p = attr(f[8], "ID"), attr(f[8], "Parent")
        tx[i] = (f[0], int(f[3]) - 1, int(f[4]), f[6], p); tx_of[p].append(i)
    elif f[2] == "exon":
        p = attr(f[8], "Parent")
        exons[p].append((int(f[3]) - 1, int(f[4])))
by_chrom = collections.defaultdict(list)
for g, (c, s, e) in genes.items():
    if tx_of.get(g): by_chrom[c].append((s, e, g))
for c in by_chrom: by_chrom[c].sort()
cnt = collections.Counter()
starts = {c: [x[0] for x in v] for c, v in by_chrom.items()}
maxlen = max(e - s for v in by_chrom.values() for s, e, _ in v)
for ln in open("a_primaries.tsv"):
    c, p = ln.split("\t"); p = int(p) - 1
    v = by_chrom.get(c)
    if not v: continue
    k = bisect.bisect_right(starts[c], p)
    for s, e, g in v[max(0, bisect.bisect_left(starts[c], p - maxlen)):k]:
        if s <= p < e: cnt[g] += 1
elig = sorted(g for g, n in cnt.items() if n >= 5)
rng = random.Random(20261024)
panel = sorted(rng.sample(elig, 300))
print(f"eligible genes {len(elig)}; panel 300")
txs = []
for g in panel:
    t = sorted(i for i in tx_of[g] if exons.get(i))
    txs += [(g, i) for i in (rng.sample(t, 5) if len(t) > 5 else t)]
# sequences: one samtools faidx call over every exon
regs = sorted({f"{tx[i][0]}:{s + 1}-{e}" for _, i in txs for s, e in exons[i]})
open("regions.txt", "w").write("\n".join(regs) + "\n")
out = subprocess.run(["samtools", "faidx", FA, "-r", "regions.txt"], capture_output=True, text=True, check=True).stdout
seq, cur = {}, None
for ln in out.splitlines():
    if ln.startswith(">"): cur = ln[1:]; seq[cur] = []
    else: seq[cur].append(ln)
seq = {k: "".join(v).upper() for k, v in seq.items()}
comp = str.maketrans("ACGTN", "TGCAN")
B = "ACGT"
with open("reads.fa", "w") as fa, open("truth.tsv", "w") as tr:
    tr.write("read\tgene\ttranscript\tchrom\tstart\tend\tstrand\tn_exons\tlen\n")
    k = 0
    for g, i in txs:
        c, s, e, strand, _ = tx[i]
        ex = sorted(exons[i])
        t = "".join(seq[f"{c}:{a + 1}-{b}"] for a, b in ex)
        if strand == "-": t = t.translate(comp)[::-1]
        for r in range(20):
            o = []
            for ch in t:
                x = rng.random()
                if x < 0.0005: o.append(rng.choice([b for b in B if b != ch]))
                elif x < 0.00075: o.append(ch); o.append(rng.choice(B))
                elif x < 0.001: pass
                else: o.append(ch)
            name = f"sim{k:06d}"; k += 1
            fa.write(f">{name}\n{''.join(o)}\n")
            tr.write(f"{name}\t{g}\t{i}\t{c}\t{s}\t{e}\t{strand}\t{len(ex)}\t{len(t)}\n")
print(f"transcripts {len(txs)}; reads {k}")
