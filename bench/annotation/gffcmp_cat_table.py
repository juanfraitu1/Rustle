# The ASSEMBLY_POLISH 30-cell scorecard (5 gffcompare metrics x 6 chromosomes, ours >= StringTie at printed precision), RefSeq vs CAT
import re
B = "/mnt/linuxdisk/home/juanfraitu/bakeoff"
def stats(p):
    t = open(p).read(); g = lambda rx: re.search(rx, t)
    ic = g(r"Intron chain level:\s+([\d.]+)\s+\|\s+([\d.]+)"); tr = g(r"Transcript level:\s+([\d.]+)\s+\|\s+([\d.]+)")
    return dict(chains=int(g(r"Matching intron chains:\s+(\d+)").group(1)), tx=int(g(r"Matching transcripts:\s+(\d+)").group(1)),
                icSn=float(ic.group(1)), icPr=float(ic.group(2)), trSn=float(tr.group(1)), trPr=float(tr.group(2)),
                ref=int(g(r"Reference mRNAs :\s+(\d+)").group(1)), q=int(g(r"Query mRNAs :\s+(\d+)").group(1)))
M = ("chains", "icSn", "icPr", "trSn", "trPr")
for r, lab in (("ref", "RefSeq"), ("ref_cat", "CAT")):
    tot = 0; pool = {a: dict.fromkeys(("chains", "tx", "q", "ref"), 0) for a in ("SHIP", "st", "flair")}
    print(f"== {lab}")
    for c in (20, 11, 7, 14, 5, 9):
        s = {a: stats(f"{B}/human_chr{c}/gffship_cat/{a}_{r}.stats") for a in ("SHIP", "st", "flair")}
        for a in s:
            for k in pool[a]: pool[a][k] += s[a][k]
        win = sum(s["SHIP"][m] >= s["st"][m] for m in M); tot += win
        print(f"chr{c}\tref {s['SHIP']['ref']}\t" + "\t".join(f"{m} {s['SHIP'][m]}/{s['st'][m]}" for m in M) + f"\tq {s['SHIP']['q']}/{s['st']['q']}\tflair chains {s['flair']['chains']}\t{win}/5")
    print(f"cells {tot}/30; pooled chains " + " ".join(f"{a} {pool[a]['chains']}" for a in pool) + "; matching tx " + " ".join(f"{a} {pool[a]['tx']}" for a in pool)
          + "; tx Pr " + " ".join(f"{a} {100*pool[a]['tx']/pool[a]['q']:.1f}" for a in pool) + "; tx Sn " + " ".join(f"{a} {100*pool[a]['tx']/pool[a]['ref']:.2f}" for a in pool) + f"; ref mRNAs {pool['SHIP']['ref']}")
