#!/usr/bin/env python3
"""r1701 assembly bars (docs/PREREG_placement_verdicts_2026-11-13.md): spliced chains (chrom, strand, introns) lost and
new in ARM against BASE, exact matches to the annotation, overlap of the lost set with a control arm, single-exon
models, and secondary-only loci (r1671 table) left without any model.
Usage: pool_arm_score.py BASE.gtf ARM.gtf CTRL.gtf ANNOT.gff LOCI.tsv CONTIG[,CONTIG...]"""
import collections, csv, gzip, re, sys

def models(path, keep, parent=None):
    ex = collections.defaultdict(list); where = {}
    for l in (gzip.open if path.endswith(".gz") else open)(path, "rt"):
        f = l.rstrip("\n").split("\t")
        if l[0] == "#" or len(f) < 9 or f[2] != "exon" or f[0] not in keep: continue
        m = re.search(parent or r'transcript_id "([^"]+)"', f[8])
        if m: ex[m.group(1)].append((int(f[3]) - 1, int(f[4]))); where[m.group(1)] = (f[0], f[6])
    ch = set(); mono = []
    for t, e in ex.items():
        e.sort(); c = tuple((a[1], b[0]) for a, b in zip(e, e[1:]))
        if c: ch.add(where[t] + (c,))
        else: mono.append((where[t][0], e[0][0], e[0][1]))
    spans = collections.defaultdict(list)
    for t, e in ex.items(): spans[where[t][0]].append((min(a for a, _ in e), max(b for _, b in e)))
    return ch, mono, spans

def main(base, arm, ctrl, annot, loci, contigs):
    keep = set(contigs.split(",")); B, bm, _ = models(base, keep); A, am, asp = models(arm, keep); C, _, _ = models(ctrl, keep)
    R, _, _ = models(annot, keep, parent=r"Parent=([^;]+)")
    lost, new = B - A, A - B; bx = len(B & R); lx = len(lost & R); nx = len(new & R)
    print(f"BASE spliced chains {len(B)}, exact annotation match {bx} ({bx / max(len(B), 1):.4f}); single-exon models BASE {len(bm)} ARM {len(am)}")
    print(f"lost {len(lost)}, exact {lx} ({lx / max(len(lost), 1):.4f}); new {len(new)}, exact {nx} ({nx / max(len(new), 1):.4f})")
    print(f"A1 lost exact / BASE exact = {lx / max(bx, 1):.4f} (bar <= 0.01): {'PASS' if lx <= 0.01 * bx else 'FAIL'}")
    print(f"A2 exact share of lost {lx / max(len(lost), 1):.4f} < BASE share {bx / max(len(B), 1):.4f}: {'PASS' if lost and lx / len(lost) < bx / max(len(B), 1) else 'FAIL'}")
    print(f"A3 lost + new = {len(lost) + len(new)} (bar >= 50): {'PASS' if len(lost) + len(new) >= 50 else 'FAIL (inert)'}")
    both = len(lost & (B - C)); print(f"A4 lost chains also lost in the control {both} of {len(lost)} = {both / max(len(lost), 1):.3f} (bar < 0.80): {'PASS' if lost and both < 0.8 * len(lost) else 'FAIL'}")
    c = collections.Counter()
    for r in csv.DictReader(open(loci), delimiter="\t"):
        if r["chrom"] not in keep: continue
        hit = any(s < int(r["end"]) and int(r["start"]) < e for s, e in asp[r["chrom"]])
        c[(r["truth"], "kept" if hit else "gone")] += 1
    print("secondary-only loci on the panel (truth label, kept / gone in ARM):", dict(sorted(c.items())))

if __name__ == "__main__": main(*sys.argv[1:7])
