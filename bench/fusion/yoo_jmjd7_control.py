"""Positive control for the fusion classifier (Rust bridge_regroup::classify_bridge, mirrored): Yoo 2025's gorilla
JMJD7-PLA2G4B copies. For every assembled transcript that runs from the annotated JMJD7(-like) gene into the
downstream PLA2G4B part, the JMJD7->PLA2G4B junction = the first intron whose acceptor exon lies past the JMJD7 gene's
3' end; primary/secondary alignments carrying it are counted in the BAM; ORF/motif rules are the Rust ones."""
import re, sys, collections, pysam
gtf, gff, bam, fasta = sys.argv[1:5]
CH = "NC_073224.2"
jm = []
for ln in open(gff):
    f = ln.split("\t")
    if len(f) > 8 and f[0] == CH and f[2] == "gene" and "JMJD7" in f[8] and 15_600_000 < int(f[3]) < 29_200_000:
        jm.append((int(f[3]), int(f[4]), f[6], re.search(r"Name=([^;]+)", f[8]).group(1)))
tx = collections.defaultdict(lambda: {"ex": []})
for ln in open(gtf):
    f = ln.rstrip("\n").split("\t")
    if len(f) < 9 or f[0] != CH or f[2] != "exon":
        continue
    t = re.search(r'transcript_id "([^"]+)"', f[8]).group(1)
    tx[t]["ex"].append((int(f[3]), int(f[4])))
    tx[t]["strand"] = f[6]
fa = pysam.FastaFile(fasta)
bf = pysam.AlignmentFile(bam)
RC = str.maketrans("ACGT", "TGCA")
def orf(s):
    best = None
    for fr in range(3):
        op = None
        for i in range(fr, len(s) - 2, 3):
            c = s[i:i+3]
            if op is None and c == "ATG": op = i
            if c in ("TAA", "TAG", "TGA") and op is not None:
                L = i + 3 - op
                if best is None or L > best[1] - best[0] or (L == best[1] - best[0] and op < best[0]): best = (op, i + 3)
                op = None
    return best
STOPS = ("TAA", "TAG", "TGA")
def longest_open_run(x):
    best = 0
    for f in range(3):
        run = 0
        for i in range(f, len(x) - 2, 3):
            run = 0 if x[i:i+3] in STOPS else run + 1
            best = max(best, run)
    return best
def fusion_frame(t, o, p):
    """Rust bridge_regroup::fusion_frame: each partner read in its own dominant open frame (review I1)."""
    a, b = o
    if a + 3 > p or p + 6 > b:
        return False
    r0 = a % 3
    for i in range(a % 3, a - 2, 3):
        if t[i:i+3] in STOPS:
            r0 = i + 3
    before = (p - r0) // 3
    fb = (3 - (p - a) % 3) % 3
    after = (b - 3 - (p + fb)) // 3
    return before >= longest_open_run(t[:p]) and after >= longest_open_run(t[p:])
def support(s, e):
    p = q = 0
    for r in bf.fetch(CH, s - 2, s + 2):
        if r.is_unmapped or r.is_supplementary: continue
        pos = r.reference_start
        for op, n in r.cigartuples:
            if op == 3 and pos + 1 == s and pos + n == e:
                if r.is_secondary: q += 1
                else: p += 1
            if op in (0, 2, 3, 7, 8): pos += n
    return p, q
cache = {}
for (gs, ge, gstr, name) in sorted(jm):
    rows = []
    for t, v in tx.items():
        ex = sorted(v["ex"]); st = v["strand"]
        if st != gstr or ex[0][0] > ge or ex[-1][1] < gs: continue
        js = None
        if st == "+":
            for a, b in zip(ex, ex[1:]):
                if a[1] <= ge and b[0] > ge: js = (a[1] + 1, b[0] - 1); break
        else:
            for a, b in zip(ex, ex[1:]):
                if b[0] >= gs and a[1] < gs: js = (a[1] + 1, b[0] - 1)
        if js is None or not (gs <= (js[0] if st == "+" else js[1]) <= ge + 1 or True): continue
        seq = "".join(fa.fetch(CH, a - 1, b).upper() for a, b in ex)
        if st == "-": seq = seq.translate(RC)[::-1]
        tot = sum(b - a + 1 for a, b in ex); acc = 0; p = None
        for a, b in zip(ex, ex[1:]):
            acc += a[1] - a[0] + 1
            if a[1] + 1 == js[0] and b[0] - 1 == js[1]: p = tot - acc if st == "-" else acc
        o = orf(seq)
        aa = (o[1] - o[0]) // 3 - 1 if o else 0
        spans = bool(o) and p is not None and fusion_frame(seq, o, p)
        d = fa.fetch(CH, js[0] - 1, js[0] + 1).upper(); ac = fa.fetch(CH, js[1] - 2, js[1]).upper()
        mot = f"{d}-{ac}" if st == "+" else f"{ac.translate(RC)[::-1]}-{d.translate(RC)[::-1]}"
        if js not in cache: cache[js] = support(*js)
        pr, se = cache[js]
        canon = mot in ("GT-AG", "GC-AG", "AT-AC")
        cls = ("fusion_gene" if aa >= 100 and spans else "readthrough") if (pr >= 4 and canon) else "unsupported"
        rows.append((t, f"{js[0]}-{js[1]}", len(ex), pr, se, mot, aa, spans, cls, "contested" if se > pr else "."))
    print(f"== {name} {gs}-{ge} {gstr}: {len(rows)} transcripts run JMJD7 -> downstream")
    for r in sorted(rows, key=lambda r: -r[6])[:4]:
        print("   ", *r)
