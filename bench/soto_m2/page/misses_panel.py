"""Three read-level pictures of why our RNA pipeline missed Soto members (July attribution labels):
collapse-K0, mis-chain, seeding gap. misses_html(prefix) returns one self-contained block (unique SVG ids per prefix)."""

W, H = 360, 188


def t(x, y, s, size=10, fill="var(--muted)", anchor="start", weight=400):
    return f'<text x="{x}" y="{y}" font-size="{size}" fill="{fill}" text-anchor="{anchor}" font-weight="{weight}">{s}</text>'


def gene(x0, y, exons, col="var(--igv)", h=9):
    lo, hi = x0 + exons[0][0], x0 + exons[-1][1]
    s = f'<line x1="{lo}" x2="{hi}" y1="{y + h / 2}" y2="{y + h / 2}" stroke="{col}"/>'
    for a, b in exons:
        s += f'<rect x="{x0 + a}" y="{y}" width="{b - a}" height="{h}" rx="1.5" fill="{col}"/>'
    return s


def read(x0, y, blocks, col, op=1, dash_last=False):
    """A spliced read: thick blocks joined by thin intron lines."""
    s = ""
    for k, (a, b) in enumerate(blocks):
        s += f'<rect x="{x0 + a}" y="{y}" width="{b - a}" height="5" rx="2" fill="{col}" opacity="{op}"/>'
        if k:
            pa = blocks[k - 1][1]
            s += f'<line x1="{x0 + pa}" x2="{x0 + a}" y1="{y + 2.5}" y2="{y + 2.5}" stroke="{col}" stroke-width=".8" opacity="{op}"/>'
    return s


def locus_box(x, y, w, label):
    return (f'<rect x="{x}" y="{y}" width="{w}" height="16" rx="3" fill="var(--track)" stroke="var(--line)"/>'
            + t(x + 5, y + 11.5, label, 9, "var(--fg)"))


EX = [(10, 30), (55, 75), (100, 120)]  # a 3-exon gene, 130 px


def k0():
    s = locus_box(8, 4, 344, "chr7 · copy A and copy B, 60 kb apart (toy)")
    for x0, lab in ((20, "copy A"), (200, "copy B")):
        s += gene(x0, 34, EX) + t(x0 + 65, 30, lab, 10, "var(--fg)", "middle", 600)
        # PSV track: differences only in introns / flanks
        for px in (2, 38, 44, 84, 91, 126):
            s += f'<line x1="{x0 + px}" x2="{x0 + px}" y1="52" y2="60" stroke="var(--cross)" stroke-width="1.6"/>'
        for k in range(4):
            s += read(x0, 70 + k * 10, [(12 + k, 30), (55, 75), (100, 118 - k)], "var(--faint)")
    s += t(175, 59, "PSVs", 9, "var(--cross)", "middle", 600)
    s += t(180, 126, "the same 4 reads fit copy A and copy B equally (MAPQ 0)", 10, "var(--fg)", "middle")
    s += t(180, 140, "red ticks = differences between the copies: all in introns or flanks,", 10, anchor="middle")
    s += t(180, 152, "which spliced reads do not contain", 10, anchor="middle")
    s += (f'<rect x="100" y="162" width="160" height="18" rx="9" fill="var(--track)" stroke="var(--line)"/>'
          + t(180, 175, "RNA sees one copy, not two", 10, "var(--fg)", "middle", 600))
    return s


def mischain(prefix):
    s = locus_box(8, 4, 344, "chr7 · NCF1-like copies, ~300 kb apart")
    s += gene(20, 34, EX) + t(85, 30, "copy A", 10, "var(--fg)", "middle", 600)
    s += gene(210, 34, EX) + t(275, 30, "copy B", 10, "var(--fg)", "middle", 600)
    s += t(180, 43, "// 300 kb //", 10, "var(--muted)", "middle")
    # a correct read on copy A
    s += read(20, 58, [(10, 30), (55, 75), (100, 120)], "var(--s1)")
    s += t(152, 63, "correct read", 9, "var(--s1)")
    # a mis-chained read: exons 1-2 at copy A, exon 3 at copy B
    s += read(20, 76, [(10, 30), (55, 75)], "var(--s1)") + read(210, 76, [(100, 120)], "var(--s1)")
    s += (f'<path d="M{20 + 75} 81 C 150 104, 260 104, {210 + 100} 81" fill="none" stroke="var(--cross)" '
          f'stroke-width="1.6" stroke-dasharray="4 3"/>')
    s += t(200, 116, "fake 300-kb intron", 10, "var(--cross)", "middle", 600)
    s += t(180, 138, "minimap2 chains one read across two near-identical copies;", 10, "var(--fg)", "middle")
    s += t(180, 150, "our filter drops reads with giant, unsupported introns", 10, "var(--fg)", "middle")
    s += t(180, 168, "capping introns at 50 kb removes the fake intron, but also", 10, anchor="middle")
    s += t(180, 180, "cuts real long introns (CDH12 has 510-kb introns)", 10, anchor="middle")
    return s


def seeding():
    s = locus_box(8, 4, 168, "chr2 · ANKRD36B")
    s += locus_box(184, 4, 168, "chr13 · its sibling")
    s += gene(20, 30, [(5, 18), (30, 42), (58, 70), (86, 98), (112, 124), (138, 150)])
    s += gene(196, 30, [(5, 18), (30, 42), (58, 70), (86, 98), (112, 124), (138, 150)])
    import random
    r = random.Random(5)
    ex = [(5, 18), (30, 42), (58, 70), (86, 98), (112, 124), (138, 150)]
    for k in range(7):
        pick = sorted(r.sample(range(6), r.randint(2, 4)))
        s += read(20, 48 + k * 8, [ex[i] for i in pick], "var(--s1)")
    for k in range(4):
        s += read(196, 48 + k * 8, [ex[i] for i in (0, 1, 2, 3)], "var(--s1)")
    s += t(92, 116, "1,056 reads, MAPQ 60", 10, "var(--s1)", "middle", 600)
    s += t(92, 128, "scattered over 1,045 intron patterns", 10, "var(--muted)", "middle")
    s += (f'<rect x="30" y="136" width="124" height="16" rx="3" fill="none" stroke="var(--cross)" stroke-dasharray="4 3"/>'
          + t(92, 148, "our locus: none built", 10, "var(--cross)", "middle", 600))
    s += (f'<rect x="206" y="92" width="124" height="16" rx="3" fill="var(--accent)" opacity=".85"/>'
          + t(268, 104, "our locus", 10, "var(--surface)", "middle", 600))
    s += t(180, 174, "different chromosomes, unique reads: alignment is fine;", 10, "var(--fg)", "middle")
    s += t(180, 186, "the loss is in our own locus building", 10, "var(--fg)", "middle")
    return s


CARDS = [
    ("Collapse-K0", "7 members", "Two nearby copies whose exons are identical (K = 0 differing positions).",
     "No aligner can prefer one copy: the read carries no difference to use. Needs DNA evidence (read depth, long genomic reads).",
     "k0"),
    ("Mis-chain", "8 members (12 in the detection tab)", "One read is aligned half to one copy and half to a far-away near-identical copy.",
     "Tried in July: a 50-kb intron cap removed all fake introns but recovered 1 of 12 members (a repeat artifact) and broke 4 real copies. The copies still merged: a symptom of K = 0.",
     "mis"),
    ("Seeding gap", "5 members (11 in the detection tab, which also counts 6 \"genuine miss\")",
     "The copy is on another chromosome or more than 3 Mb away and has unique reads, but our pipeline never built a locus for it.",
     "The reads already align correctly, so re-aligning changes nothing. This one is ours to fix.",
     "seed"),
]


def misses_html(prefix):
    out = [f'<section class="miss" aria-label="Why RNA missed some Soto members">',
           '<div class="eyebrow">What the RNA miss reasons mean · read-level pictures (toy loci except where named)</div>',
           '<div class="miss-grid">']
    for title, n, what, fix, key in CARDS:
        body = {"k0": k0(), "mis": mischain(prefix), "seed": seeding()}[key]
        h = 188 if key != "seed" else 192
        out.append(f'<div class="mcard"><div class="mh"><b>{title}</b><span>{n}</span></div>'
                   f'<svg viewBox="0 0 {W} {h}" width="100%" role="img" aria-label="{title}: read-level picture">{body}</svg>'
                   f'<p><b>What:</b> {what}</p><p><b>Why re-aligning can\'t fix it:</b> {fix}</p></div>')
    out.append('</div></section>')
    return "\n".join(out)


if __name__ == "__main__":
    print(misses_html("x"))
