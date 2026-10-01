"""Builds the "How Soto builds families" tab: eight steps of Soto's released code (A_SD98_regions.md, B_SD98_families.ipynb)
drawn on one toy example. Prints the pane HTML to stdout."""

W = 760
PX = [10, 200, 390, 580]   # four locus panels (IGV multi-locus), 170 px each
PW = 170
LOCI = ["chrA:10,200-10,220 kb", "chrB:44,100-44,120 kb", "chrC:7,300-7,320 kb", "chrD:88,000-88,020 kb"]
SD = [(20, 150), (20, 150), (20, 88), (20, 110)]           # SD98 block per panel
GENES = {  # name, panel, biotype, exons (panel x), row
    "GENE1": (0, "pc", [(30, 50), (75, 95), (120, 140)], 0),
    "LNC1": (0, "lnc", [(28, 52)], 1),
    "GENE1B": (1, "pc", [(30, 50), (75, 95), (120, 140)], 0),
    "LNC2": (1, "lnc", [(118, 142)], 1),
    "GENE1P": (2, "up", [(30, 50), (75, 95)], 0),
    "GENE1P2": (3, "up", [(30, 50), (75, 95)], 0),
}
COL = {"pc": "var(--s1)", "up": "var(--s2)", "lnc": "var(--s3)"}
BT = {"pc": "protein-coding", "up": "unprocessed pseudogene", "lnc": "lncRNA"}
CN = {"GENE1": 6.1, "GENE1B": 6.4, "GENE1P": 13.8, "GENE1P2": 13.1}


def svg(h, body, label):
    return (f'<svg viewBox="0 0 {W} {h}" width="100%" role="img" aria-label="{label}" '
            f'xmlns="http://www.w3.org/2000/svg">{body}</svg>')


def txt(x, y, s, size=10.5, fill="var(--muted)", anchor="start", weight=400, mono=False):
    fam = ' font-family="var(--f-mono)"' if mono else ""
    return (f'<text x="{x}" y="{y}" font-size="{size}" fill="{fill}" text-anchor="{anchor}" '
            f'font-weight="{weight}"{fam}>{s}</text>')


def panels(y0, sd=True, genes=True, exon_mode="plain", hide=()):
    """IGV-style locus panels: locus box, SD98 track, gene tracks. exon_mode: plain | inside (colour exons fully
    inside SD98, grey the rest)."""
    s = ""
    for p, x in enumerate(PX):
        s += f'<rect x="{x}" y="{y0}" width="{PW}" height="18" rx="3" fill="var(--track)" stroke="var(--line)"/>'
        s += txt(x + 6, y0 + 13, LOCI[p], 9.5, "var(--fg)", mono=True)
        if sd:
            a, b = SD[p]
            s += f'<rect x="{x + a}" y="{y0 + 26}" width="{b - a}" height="10" rx="2" fill="var(--s4)" opacity=".75"/>'
        if genes:
            for g, (pp, bt, ex, row) in GENES.items():
                if pp != p or g in hide:
                    continue
                gy = y0 + 46 + row * 26
                lo, hi = ex[0][0], ex[-1][1]
                s += f'<line x1="{x + lo}" x2="{x + hi}" y1="{gy + 5}" y2="{gy + 5}" stroke="var(--igv)"/>'
                for a, b in ex:
                    inside = SD[p][0] <= a and b <= SD[p][1]
                    fill = "var(--igv)"
                    extra = ""
                    if exon_mode == "inside":
                        fill = "var(--accent)" if inside else "var(--faint)"
                        if not inside:
                            extra = (f'<path d="M{x + a + 4} {gy - 2} l12 14 M{x + a + 16} {gy - 2} l-12 14" '
                                     f'stroke="var(--cross)" stroke-width="1.6"/>')
                    s += f'<rect x="{x + a}" y="{gy}" width="{b - a}" height="10" rx="1.5" fill="{fill}"/>{extra}'
                s += txt(x + hi + 4, gy + 9, g, 9.5, COL[bt], weight=600)
    return s


def legend_row(y, items):
    s, x = "", 10
    for col, lab in items:
        s += f'<rect x="{x}" y="{y - 8}" width="10" height="10" rx="2" fill="{col}"/>' + txt(x + 15, y, lab, 10)
        x += 26 + len(lab) * 6.2
    return s


def step1():
    s = panels(4, genes=False)
    # arcs between SD98 blocks (links track)
    for a, b in [(0, 1), (0, 2), (0, 3), (1, 2), (2, 3)]:
        xa = PX[a] + (SD[a][0] + SD[a][1]) / 2
        xb = PX[b] + (SD[b][0] + SD[b][1]) / 2
        h = 10 + abs(xb - xa) * 0.05
        s += f'<path d="M{xa} 44 Q{(xa + xb) / 2} {44 + 2 * h} {xb} 44" fill="none" stroke="var(--s4)" stroke-width="1.6"/>'
    s += legend_row(112, [("var(--s4)", "SD98 block: a stretch whose copy elsewhere is ≥ 98% identical (SEDEF)")])
    return svg(120, s, "Four loci with their SD98 blocks linked by arcs")


def step2():
    s = panels(4, exon_mode="inside")
    s += txt(PX[2] + 100, 76, "exon crosses the block edge", 9.5, "var(--cross)")
    s += legend_row(118, [("var(--accent)", "exon kept (fully inside SD98)"), ("var(--faint)", "exon dropped"),
                          ("var(--s1)", "protein-coding"), ("var(--s2)", "unprocessed pseudogene"), ("var(--s3)", "lncRNA")])
    return svg(126, s, "Exons fully inside SD98 kept, of every biotype")


def step3():
    s = panels(4, exon_mode="inside", hide=("LNC1", "LNC2"))
    # GENE1 exon 1 aligns at every copy; LNC1's exon also aligns there
    sx = PX[0] + 40
    for p in (1, 2, 3):
        tx = PX[p] + 40
        s += (f'<path d="M{sx} 60 C{sx + 40} 118, {tx - 40} 118, {tx} 60" fill="none" stroke="var(--accent)" '
              f'stroke-width="1.4" marker-end="url(#ah)"/>')
    s = ('<defs><marker id="ah" viewBox="0 0 8 8" refX="7" refY="4" markerWidth="7" markerHeight="7" orient="auto">'
         '<path d="M0 0 L8 4 L0 8 z" fill="var(--accent)"/></marker></defs>') + s
    s += txt(PX[0] + 30, 124, "GENE1 exon 1 aligns at every copy; every other kept exon, lncRNA exons included, is mapped the same way", 10.5, "var(--fg)")
    return svg(132, s, "Each kept exon mapped back to the genome")


NODES = {"GENE1": (130, 32), "GENE1B": (130, 112), "LNC1": (330, 158), "LNC2": (36, 72),
         "GENE1P": (530, 32), "GENE1P2": (530, 112)}
EDGES = [("GENE1", "GENE1B"), ("GENE1", "GENE1P"), ("GENE1", "GENE1P2"), ("GENE1B", "GENE1P"), ("GENE1B", "GENE1P2"),
         ("GENE1P", "GENE1P2"), ("LNC1", "GENE1"), ("LNC1", "GENE1B"), ("LNC1", "GENE1P"), ("LNC1", "GENE1P2"),
         ("LNC2", "GENE1"), ("LNC2", "GENE1B")]


def bt_of(g):
    return GENES[g][1]


def graph(gate=False):
    s = ""
    for a, b in EDGES:
        (xa, ya), (xb, yb) = NODES[a], NODES[b]
        f = 0.5 if (xa == xb or ya == yb) else 0.27
        mx, my = xa + (xb - xa) * f, ya + (yb - ya) * f
        if not gate:
            s += f'<line x1="{xa}" y1="{ya}" x2="{xb}" y2="{yb}" stroke="var(--faint)" stroke-width="1.4"/>'
            continue
        if a in CN and b in CN:
            d = abs(CN[a] - CN[b])
            keep = d < 2
            s += (f'<line x1="{xa}" y1="{ya}" x2="{xb}" y2="{yb}" stroke="{"var(--fg)" if keep else "var(--cross)"}" '
                  f'stroke-width="{2.2 if keep else 1.4}"{"" if keep else " stroke-dasharray=\"5 4\""}/>')
            s += (f'<rect x="{mx - 22}" y="{my - 9}" width="44" height="15" rx="7" fill="var(--surface)"/>'
                  + txt(mx, my + 3, f"Δ {d:.1f}", 10, "var(--fg)" if keep else "var(--cross)", "middle", 600))
        else:
            s += f'<line x1="{xa}" y1="{ya}" x2="{xb}" y2="{yb}" stroke="var(--s3)" stroke-width="1.4"/>'
    for g, (x, y) in NODES.items():
        s += f'<circle cx="{x}" cy="{y}" r="9" fill="{COL[bt_of(g)]}" stroke="var(--surface)" stroke-width="2"/>'
        s += txt(x, y - 15, g, 10.5, "var(--fg)", "middle", 600)
    return s


def step4():
    s = graph()
    s += legend_row(196, [("var(--s1)", "protein-coding"), ("var(--s2)", "unprocessed pseudogene"), ("var(--s3)", "lncRNA")])
    s += txt(640, 70, "pair kept only if", 10.5, "var(--fg)") + txt(640, 84, "one side is coding", 10.5, "var(--fg)") + \
        txt(640, 98, "or an unprocessed", 10.5, "var(--fg)") + txt(640, 112, "pseudogene", 10.5, "var(--fg)")
    return svg(204, s, "Genes paired when one's exon covers the other's")


def strip_values(center, sd, n=268, seed=7):
    import random
    r = random.Random(seed)
    v = sorted(r.gauss(center, sd) for _ in range(n))
    med = (v[n // 2 - 1] + v[n // 2]) / 2
    return [x - med + center for x in v]


def step5():
    s, y = "", 16
    X0, SC = 150, 30
    # where one number comes from: 268 genomes -> one median per gene
    import random
    r = random.Random(3)
    for g, c, sd, yy in (("GENE1", 6.1, 0.35, 24), ("GENE1P", 13.8, 0.7, 56)):
        bt = bt_of(g)
        s += txt(X0 - 10, yy + 4, g, 10.5, COL[bt], "end", 600)
        for v in strip_values(c, sd, seed=11 if g == "GENE1" else 12):
            s += f'<circle cx="{X0 + v * SC:.1f}" cy="{yy + r.uniform(-7, 7):.1f}" r="1.6" fill="{COL[bt]}" opacity=".45"/>'
        s += (f'<line x1="{X0 + c * SC}" x2="{X0 + c * SC}" y1="{yy - 11}" y2="{yy + 11}" stroke="var(--fg)" stroke-width="2"/>'
              + txt(X0 + c * SC, yy - 13, f"median {c}", 10, "var(--fg)", "middle", 600))
    s += txt(X0, 84, "each dot = the gene's copy number in one of 268 genomes (illustrative spread); the median becomes the gene's famCN", 10)
    y = 100
    for t in range(0, 17, 4):
        x = X0 + t * SC
        s += f'<line x1="{x}" x2="{x}" y1="94" y2="258" stroke="var(--line)"/>' + txt(x, 272, t, 10, anchor="middle")
    for g in ["GENE1", "GENE1B", "GENE1P", "GENE1P2", "LNC1", "LNC2"]:
        bt = bt_of(g)
        s += txt(X0 - 10, y + 10, g, 10.5, COL[bt], "end", 600)
        if g in CN:
            s += (f'<path d="M{X0} {y} h{CN[g] * SC - 3} q3 0 3 3 v8 q0 3 -3 3 h-{CN[g] * SC - 3} z" fill="{COL[bt]}"/>'
                  + txt(X0 + CN[g] * SC + 6, y + 10, f"{CN[g]}", 10.5, "var(--fg)", weight=600))
        else:
            s += txt(X0 + 4, y + 10, "not measured (only coding genes and unprocessed pseudogenes get a copy number)", 10.5)
        y += 26
    s += txt(X0 + 8 * SC, 288, "famCN (diploid copies of that sequence in the genome)", 10.5, "var(--muted)", "middle")
    return svg(294, s, "Copy number per gene: the median over 268 genomes")


def mad_figure():
    X0, SC = 170, 15   # copy number 0..32
    s = ""
    for t in range(0, 33, 4):
        x = X0 + t * SC
        s += f'<line x1="{x}" x2="{x}" y1="18" y2="150" stroke="var(--line)"/>' + txt(x, 164, t, 10, anchor="middle")
    s += txt(X0 + 16 * SC, 180, "copy number (famCN)", 10.5, "var(--muted)", "middle")
    rows = [("GENE1 vs GENE1B", 6.1, 6.4, "toy"), ("GENE1 vs GENE1P", 6.1, 13.8, "toy"),
            ("NPIPB3 vs NPIPB4", 26.5, 29.1, "real, Soto's table")]
    for k, (lab, a, b, src) in enumerate(rows):
        y = 36 + k * 44
        med = (a + b) / 2
        mad = abs(a - b) / 2
        keep = mad < 1
        col = "var(--fg)" if keep else "var(--cross)"
        s += txt(X0 - 12, y + 4, lab, 10.5, "var(--fg)", "end", 600) + txt(X0 - 12, y + 17, src, 9.5, "var(--muted)", "end")
        s += f'<line x1="{X0 + a * SC}" x2="{X0 + b * SC}" y1="{y}" y2="{y}" stroke="{col}" stroke-width="2"/>'
        s += f'<line x1="{X0 + med * SC}" x2="{X0 + med * SC}" y1="{y - 9}" y2="{y + 9}" stroke="var(--muted)" stroke-dasharray="2 2"/>'
        for v in (a, b):
            s += f'<circle cx="{X0 + v * SC}" cy="{y}" r="5" fill="var(--surface)" stroke="{col}" stroke-width="2"/>'
        s += txt(X0 + max(a, b) * SC + 12, y + 4,
                 f"MAD = {mad:.2f} → {'kept' if keep else 'cut'}", 10.5, col, weight=650)
    s += txt(X0, 10, "dashed tick = the pair's median; each value sits half the gap away from it, so MAD = half the gap", 10)
    return svg(188, s, "MAD of a pair of copy numbers is half their gap")


def step6():
    s = graph(gate=True)
    s += (legend_row(196, [("var(--fg)", "kept: copy numbers differ by &lt; 2"), ("var(--cross)", "cut: differ by ≥ 2"),
                           ("var(--s3)", "kept: a gene without copy number always passes")]))
    return mad_figure() + svg(204, s, "Pairs cut when copy numbers differ by 2 or more")


def step7():
    s = ""
    fams = [("Family A", ["GENE1", "GENE1B", "LNC1", "LNC2"], 20), ("Family B", ["GENE1P", "GENE1P2", "LNC1"], 400)]
    for name, mem, x in fams:
        s += f'<rect x="{x}" y="8" width="340" height="104" rx="10" fill="var(--track)" stroke="var(--line)"/>'
        s += txt(x + 14, 30, name, 12, "var(--fg)", weight=650)
        for k, g in enumerate(mem):
            cx = x + 14 + (k % 2) * 160
            cy = 50 + (k // 2) * 30
            bt = bt_of(g)
            ring = ' stroke="var(--cross)" stroke-width="2"' if g == "LNC1" else ""
            s += f'<rect x="{cx}" y="{cy}" width="150" height="22" rx="11" fill="var(--surface)"{ring}/>'
            s += f'<circle cx="{cx + 13}" cy="{cy + 11}" r="5" fill="{COL[bt]}"/>' + txt(cx + 24, cy + 15, g, 10.5, "var(--fg)", weight=600)
    s += txt(380, 132, "LNC1 is in both families: a gene without copy number joins every family it pairs with", 10.5,
             "var(--cross)", "middle")
    return svg(140, s, "Two families; LNC1 belongs to both")


def step8():
    rows = [("A", "GENE1", "protein_coding", "Yes", "6.1"), ("A", "GENE1B", "protein_coding", "Yes", "6.4"),
            ("A", "LNC1", "lncRNA", "No", "N/A"), ("A", "LNC2", "lncRNA", "No", "N/A"),
            ("B", "GENE1P", "unprocessed_pseudogene", "Yes", "13.8"), ("B", "GENE1P2", "unprocessed_pseudogene", "Yes", "13.1"),
            ("B", "LNC1", "lncRNA", "No", "N/A")]
    tr = "".join(f'<tr><td class="mono">ID_{f}</td><td class="mono">{g}</td><td>{b}</td><td>{i}</td><td class="n">{c}</td></tr>'
                 for f, g, b, i, c in rows)
    return ('<div class="tbl"><table><thead><tr><th>Family ID</th><th>Gene Name</th><th>Biotype</th>'
            '<th>In Table S1 (SD98 gene set)</th><th>Median famCN</th></tr></thead><tbody>' + tr + '</tbody></table></div>')


STEPS = [
    ("Find the ≥ 98% duplications", "SEDEF segmental duplications on CHM13 v1.0 with ≥ 98% identity, merged; chrX, chrY and chrM left out.",
     "chm13.draft_v1.0_plus38Y.SDs-98.merged.bed", step1),
    ("Keep the exons that lie fully inside them", "Exons of every annotated gene (CAT v4), of any biotype, that sit entirely inside an SD98 block. 71 genes on a hand-made list are removed.",
     "bedtools intersect -wa -f 1 -a CHM13.combined.v4.exons.bed … | grep -Fvf genes_to_remove.txt", step2),
    ("Map every kept exon back to the genome", "Each kept exon's sequence is aligned to the whole genome (CHM13 v1.0), not only to the SD98 blocks, keeping up to 50 places. The SD98 blocks only decide which exons are queried (step 2) and which hits count: a hit links two genes only where it lands on another kept exon (step 4).",
     "minimap2 -c --end-bonus 5 --eqx -N 50 -p 0.5", step3),
    ("Pair genes that share an exon", "Two genes pair when an alignment of one gene's exon covers ≥ 99% of the other gene's exon, on the same strand. A pair needs at least one protein-coding gene or unprocessed pseudogene.",
     'bedtools intersect -f 0.99 -s … | grep "protein_coding\\|unprocessed_pseudogene"', step4),
    ("Measure each gene's copy number", "Read depth (WSSD) over the gene's SD98 part in each of 268 SGDP genomes (the paper says 269; their code drops one outlier); the median across genomes is the gene's famCN. Reads from every copy pile onto each copy, so famCN counts how many copies of that sequence the genome carries.",
     "SD98_WSSD.tsv (median famCN per gene)", step5),
    ("Cut pairs whose copy numbers disagree",
     "<b>What it compares:</b> the copy numbers of the two genes in a pair, with each other. Nothing else is involved: not the reference genome, "
     "not other people, not apes. Each gene's copy number is already a single value, the median over 268 genomes (step 5).<br>"
     "<b>How:</b> MAD is the median distance of the values from their own median (scipy <code>median_abs_deviation</code>, no rescaling). "
     "For two values it is exactly half their gap, so MAD &lt; 1 means the two copy numbers differ by less than 2.<br>"
     "<b>Details:</b> in their code a gene can give several values, one per piece of the gene inside SD98, and the MAD is taken over all of them. "
     "The cut-off is absolute: 2 copies apart whether a gene has 4 copies or 40.",
     "stats.median_abs_deviation(famCN) < 1", step6),
    ("Grow families through the kept pairs", "Starting from each protein-coding gene or unprocessed pseudogene, every kept pair that shares a protein-coding gene or unprocessed pseudogene with the family is merged in. A gene without copy number never pulls in anything else.",
     "B_SD98_families.ipynb, cells 10-11", step7),
    ("Publish the table", "One row per gene and family. Genes in no family are listed as Unassigned. Two families are merged by hand (\"Manual merge\").",
     "Table S1C", step8),
]


def pane():
    out = ['<div role="tabpanel" id="pane-how" aria-labelledby="tab-how" class="pane" hidden>',
           '<section class="how-intro"><div class="eyebrow">Soto et al. 2025 · family construction as in their released code '
           '(github.com/mydennislab/HSD_brain_evolution)</div><h2>How Soto builds a gene family</h2>'
           '<p class="note">One toy example runs through every step: two protein-coding genes, two unprocessed pseudogenes '
           'and two lncRNAs on four chromosomes. The names and numbers are illustrative; the rules are theirs.</p></section>',
           '<ol class="steps">']
    for k, (title, what, code, fn) in enumerate(STEPS, 1):
        out.append(f'<li class="step"><div class="sh"><span class="sn">{k}</span><div><h3>{title}</h3><p>{what}</p>'
                   f'<code class="sc">{code}</code></div></div><div class="sv">{fn()}</div></li>')
    out.append('</ol>')
    out.append('<section class="how-foot"><div class="facts3">'
               '<div><b>491</b> families</div><div><b>114</b> genes in no family (Unassigned)</div>'
               '<div><b>149</b> genes in more than one family, none of them coding</div>'
               '<div><b>2</b> families merged by hand</div><div><b>71</b> genes removed by hand before step 2</div></div>'
               '<p class="note">The paper\'s Methods describe mapping whole SD98 regions and one copy-number test per group. '
               'The released code maps exons and tests each pair. Its steps, with the family-growing loop run to completion, '
               'reproduce their table at ARI 0.97; the Methods wording gives 0.71. As released, the loop stops early and depends on '
               'hash order (0.82 to 0.87). This tab follows the code.</p>'
               '<div class="bar"><span class="note">See the same steps on real genes:</span>'
               '<button class="chip" data-real="ID_151">NPIPB3: split off by copy number</button>'
               '<button class="chip" data-real="ID_349">EEF1A1 and 2 retrocopies</button>'
               '<button class="chip" data-real="ID_468">TBC1D3</button>'
               '<button class="chip" data-real="ID_482">UBTFL: merged by hand</button></div></section>')
    out.append('</div>')
    return "\n".join(out)


if __name__ == "__main__":
    print(pane())
