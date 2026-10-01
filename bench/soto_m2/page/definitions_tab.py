"""The "Two definitions" tab: what "gene family" means in Soto 2025 and in the evolutionary sense used by the thesis.
Neutral, advisor-facing: Soto's definition is operational and stated explicitly; the two differ in scope. Prints the
pane HTML."""


def nested():
    s = []
    s.append('<rect x="8" y="8" width="604" height="300" rx="16" fill="var(--track)" stroke="var(--line)"/>')
    s.append('<text x="26" y="34" font-size="14" font-weight="650" fill="var(--fg)">Homology family</text>')
    s.append('<text x="26" y="51" font-size="11.5" fill="var(--muted)">evolutionary definition: genes descended from one ancestral gene by duplication</text>')
    s.append('<rect x="26" y="66" width="400" height="226" rx="12" fill="var(--s1)" fill-opacity=".10" stroke="var(--s1)"/>')
    s.append('<text x="42" y="90" font-size="13.5" font-weight="650" fill="var(--fg)">Sequence family</text>')
    s.append('<text x="42" y="107" font-size="11.5" fill="var(--muted)">genes sharing ≥ 98%-identical exons (SD98)</text>')
    s.append('<text x="42" y="122" font-size="11.5" fill="var(--muted)">Soto\'s own edge, before any copy-number cut</text>')
    for x, w, lab in ((40, 120, "Soto family"), (166, 120, "Soto family"), (292, 120, "Soto family")):
        s.append(f'<rect x="{x}" y="136" width="{w}" height="72" rx="10" fill="var(--s2)" fill-opacity=".18" stroke="var(--s2)"/>')
        s.append(f'<text x="{x + w / 2}" y="166" font-size="12.5" font-weight="600" fill="var(--fg)" text-anchor="middle">{lab}</text>')
        s.append(f'<text x="{x + w / 2}" y="184" font-size="10.5" fill="var(--muted)" text-anchor="middle">similar copy number</text>')
    s.append('<text x="226" y="236" font-size="11.5" fill="var(--fg)" text-anchor="middle">the sequence family is cut where copy numbers differ</text>')
    s.append('<text x="226" y="253" font-size="11.5" fill="var(--fg)" text-anchor="middle">(MAD &lt; 1, i.e. less than 2 copies apart)</text>')
    for cx, cy in ((472, 186), (502, 214), (548, 192), (484, 250), (530, 262), (566, 234)):
        s.append(f'<circle cx="{cx}" cy="{cy}" r="7" fill="var(--s3)" fill-opacity=".55" stroke="var(--s3)"/>')
    for k, (txt, size, weight, col) in enumerate((("older paralogs", 12, 600, "var(--fg)"), ("below 98% identity:", 11, 400, "var(--muted)"),
                                                   ("same homology family,", 11, 400, "var(--muted)"), ("outside SD98 by design", 11, 400, "var(--muted)"))):
        s.append(f'<text x="519" y="{96 + 16 * k}" font-size="{size}" font-weight="{weight}" fill="{col}" text-anchor="middle">{txt}</text>')
    return ('<svg viewBox="0 0 620 316" width="100%" role="img" aria-label="Nested sets: Soto families inside sequence '
            'families inside homology families, with older paralogs inside the homology family only">' + "".join(s) + '</svg>')


def axis():
    s = []
    W = 760
    x0, x1, x2, x3 = 40, 200, 470, 730
    s.append(f'<rect x="{x0}" y="70" width="{x1 - x0}" height="34" fill="var(--s2)" fill-opacity=".22"/>')
    s.append(f'<rect x="{x1}" y="70" width="{x2 - x1}" height="34" fill="var(--s1)" fill-opacity=".14"/>')
    s.append(f'<rect x="{x2}" y="70" width="{x3 - x2}" height="34" fill="var(--s3)" fill-opacity=".18"/>')
    s.append(f'<text x="{(x0 + x1) / 2}" y="92" font-size="11.5" font-weight="600" fill="var(--fg)" text-anchor="middle">&gt; 98% identical DNA</text>')
    s.append(f'<text x="{(x1 + x2) / 2}" y="92" font-size="11.5" font-weight="600" fill="var(--fg)" text-anchor="middle">DNA and RNA still alike</text>')
    s.append(f'<text x="{(x2 + x3) / 2}" y="92" font-size="11.5" font-weight="600" fill="var(--fg)" text-anchor="middle">DNA and RNA drift; protein persists</text>')
    # scope brackets
    s.append(f'<path d="M{x0} 52 V44 H{x1} V52" fill="none" stroke="var(--s2)" stroke-width="2"/>')
    s.append(f'<text x="{(x0 + x1) / 2}" y="36" font-size="12" font-weight="650" fill="var(--fg)" text-anchor="middle">Soto 2025</text>')
    s.append(f'<path d="M{x0} 24 V16 H{x3} V24" fill="none" stroke="var(--fg)" stroke-width="1.5"/>')
    s.append(f'<text x="{(x1 + x3) / 2 + 60}" y="11" font-size="12" font-weight="650" fill="var(--fg)" text-anchor="middle">evolutionary definition (this thesis)</text>')
    # time axis
    s.append(f'<line x1="{x0}" y1="118" x2="{x3}" y2="118" stroke="var(--muted)" stroke-width="1.5"/>')
    s.append(f'<path d="M{x3} 118 l-8 -4 v8 z" fill="var(--muted)"/>')
    s.append(f'<text x="{x0}" y="136" font-size="11.5" fill="var(--muted)">recent duplication</text>')
    s.append(f'<text x="{x3}" y="136" font-size="11.5" fill="var(--muted)" text-anchor="end">older duplication (identity falls)</text>')
    s.append(f'<text x="{(x0 + x1) / 2}" y="156" font-size="11" fill="var(--muted)" text-anchor="middle">evidence: shared exons</text>')
    s.append(f'<text x="{(x1 + x2) / 2}" y="156" font-size="11" fill="var(--muted)" text-anchor="middle">evidence: DNA / RNA alignment</text>')
    s.append(f'<text x="{(x2 + x3) / 2}" y="156" font-size="11" fill="var(--muted)" text-anchor="middle">evidence: protein similarity</text>')
    return (f'<svg viewBox="0 0 {W} 164" width="100%" role="img" aria-label="Which duplications each definition sees, '
            f'from recent to old">' + "".join(s) + '</svg>')


ROWS = [
    ("Question it serves", "which genes expanded specifically in humans", "which genes descend from a common ancestral gene"),
    ("Duplications it sees", "recent: &gt; 98% identical (SD98)", "all ages"),
    ("Evidence", "shared exons between annotated genes", "DNA / RNA similarity for recent copies; protein similarity once DNA has diverged"),
    ("Copy number", "part of the definition: homologs with different copy numbers are different families", "not part of the definition: copy number is a property of a family"),
    ("Non-coding features", "lncRNAs and processed pseudogenes are attached to families", "follows homology of the sequence"),
    ("A family is", "a sequence family cut by copy number", "a homology group"),
]

QUOTES = [
    ("SD98 genes were grouped into gene families based on shared exons.", "STAR Methods, Gene family clustering"),
    ("Based on sequence and famCN similarity, we clustered 1,679 of the paralogs into 491 multigene families.", "Results"),
    ("After comparing the median famCN values of SD98 genes with shared exons, groupings where the mean absolute deviation "
     "of the CN was less than one were selected.", "STAR Methods, Gene family clustering"),
    ("SD98 genes associated with other gene features, including lncRNAs and processed pseudogenes, were also assigned a "
     "gene family ID.", "STAR Methods, Gene family clustering"),
    ("A notable limitation of our study is its reliance on existing gene annotations, which we used to group "
     "human-duplicated paralogs into larger multigene families based on shared annotated sequences in SD98 regions.",
     "Limitations of the study"),
]

CLUSTERS = [('dev', 0.4863, 'ID_100,ID_192,ID_62,ID_99', 'AC016629.3, AC138393.2, AC139099.2, AL669831.4'), ('dev', 0.1699, 'ID_214,ID_396,ID_397,ID_400', 'AC239809.3, HYDIN, HYDIN2, NBPF1'), ('heldout', 0.3232, 'ID_23,ID_33,ID_34', 'AC006453.2, AC027612.1, AL356585.2, AP000550.3'), ('heldout', 0.3419, 'ID_126,ID_147', 'AC098826.2, AC125634.1, AL353626.3, ANKRD30BP2'), ('dev', 0.516, 'ID_104,ID_141,ID_142,ID_182,ID_279,ID_280', 'AC073464.1, AC119751.3, AC119751.5, AC137800.1'), ('heldout', 0.1372, 'ID_269,ID_275', 'AL591479.1, C2orf27AP3, CR382287.1'), ('heldout', 0.1541, 'ID_215,ID_250', 'AC239860.1, AC239860.2, AC241952.1, AC244394.2'), ('dev', 0.5953, 'ID_241,ID_363,ID_364', 'AL161457.2, FRG1BP, FRG1FP, FRG1GP'), ('dev', 0.3437, 'ID_403,ID_404', 'NF1P1, NF1P10, NF1P11, NF1P2'), ('heldout', 0.4664, 'ID_300,ID_301,ID_84', 'AC026273.1, BMS1P11, BMS1P12, BMS1P13'), ('dev', 0.3976, 'ID_156,ID_72', 'AC023310.3, AC091304.4, AC100756.1, AC100757.3'), ('dev', 0.2016, 'ID_113,ID_21,ID_368,ID_78,ID_88', 'AC006328.1, AC044860.1, AC091057.4, AC243562.1'), ('dev', 0.4397, 'ID_178,ID_179', 'CHRFAM7A, CHRNA7, ULK4P1, ULK4P2'), ('dev', 0.3647, 'ID_324,ID_481', 'CSPG4P10, CSPG4P11, CSPG4P12, CSPG4P4Y'), ('heldout', 0.375, 'ID_149,ID_154,ID_155,ID_41', 'AC009086.2, AC126755.1, AC138932.1, AC138969.1'), ('dev', -0.0288, 'ID_69,ID_76', 'AC022145.2, AC024257.2, AC113404.2, AL034380.2'), ('dev', 0.0163, 'ID_172,ID_184', 'AC133919.3, CICP19, CICP24, CICP26'), ('dev', 0.0489, 'ID_468,ID_469', 'TBC1D3, TBC1D3B, TBC1D3D, TBC1D3E'), ('heldout', 0.1722, 'ID_114,ID_14', 'AC005562.2, AC090616.5, AC090616.6, LRRC37A'), ('dev', 0.6214, 'ID_303,ID_359', 'BX088651.1, BX664615.1, BX664615.2, FGF7P1'), ('heldout', 0.2878, 'ID_395,ID_63', 'AC017002.3, AC083899.1, AC097527.1, ANAPC1'), ('heldout', 0.2377, 'ID_106,ID_107', 'AC087203.1, USP17L1, USP17L10, USP17L11'), ('dev', -0.1991, 'ID_163,ID_191', 'AC131392.1, AC138866.2, AC146949.1, AL021368.2'), ('heldout', 0.5626, 'ID_10,ID_11,ID_8', 'AC211476.3, AC211486.3, PMS2P1, PMS2P10'), ('heldout', 0.2108, 'ID_181,ID_187,ID_188,ID_93', 'AC055876.1, AC138649.5, FP700111.1, HERC2'), ('heldout', 0.8175, 'ID_355,ID_356', 'FAM86B1, FAM86B2, FAM90A10P, FAM90A11P'), ('heldout', -0.001, 'ID_96,ID_97', 'AC068587.6, AC068587.8, AC105233.3, AC134684.3'), ('heldout', 0.0, 'ID_271,ID_272', 'AL627230.1'), ('dev', 0.7372, 'ID_161,ID_330,ID_347,ID_62', 'AC016629.3, AC131281.1, CENPBD1P1, DUX4'), ('dev', 0.5065, 'ID_233,ID_270', 'AGAP10P, AGAP12P, AGAP13P, AGAP14P'), ('dev', 0.3518, 'ID_352,ID_485', 'FAM21EP, FAM21FP, WASHC2A, WASHC2C'), ('heldout', 0.1688, 'ID_289,ID_476', 'AP003122.1, AP004607.5, AP004607.8, AP005435.2'), ('heldout', 0.3101, 'ID_477,ID_478', 'TRIM49, TRIM49C, TRIM49D1, TRIM49D2')]  # docs/SOTO_CN_DUPLICON_CLUSTERS_2026-09-30.tsv (KEY=cnduplicon, frozen result)



def three_objects():
    s = []
    segs = [(120, 220, "var(--s1)", "D1"), (220, 340, "var(--s2)", "D2"), (340, 440, "var(--s3)", "D3"), (440, 550, "var(--s4)", "D4"),
            (550, 650, "var(--s2)", "D2"), (650, 740, "var(--s1)", "D1")]
    s.append('<path d="M120 40 V32 H738 V40" fill="none" stroke="var(--fg)" stroke-width="1.5"/>')
    s.append('<text x="429" y="24" font-size="12.5" font-weight="650" fill="var(--fg)" text-anchor="middle">SD block: one duplicated stretch (SEDEF; Soto\'s "SD Unit")</text>')
    for a, b, c, lab in segs:
        s.append(f'<rect x="{a}" y="48" width="{b - a - 2}" height="26" rx="4" fill="{c}" fill-opacity=".75"/>')
        s.append(f'<text x="{(a + b) / 2}" y="66" font-size="12" font-weight="650" fill="#fff" text-anchor="middle">{lab}</text>')
    s.append('<text x="110" y="66" font-size="11.5" fill="var(--muted)" text-anchor="end">duplicons</text>')
    genes = [(356, 532, "fusion gene: crosses the D3 | D4 boundary", 96), (130, 210, "gene of family A (on D1)", 120), (230, 330, "gene of family B (on D2)", 144)]
    for a, b, lab, y in genes:
        s.append(f'<line x1="{a}" x2="{b}" y1="{y}" y2="{y}" stroke="var(--igv)" stroke-width="2"/>')
        for x in range(a, b - 10, 32):
            s.append(f'<rect x="{x}" y="{y - 5}" width="16" height="10" rx="2" fill="var(--igv)"/>')
        s.append(f'<path d="M{b} {y} l-7 -5 v10 z" fill="var(--igv)"/>')
        if lab:
            s.append(f'<text x="{b + 10}" y="{y + 4}" font-size="11.5" fill="var(--fg)">{lab}</text>')
    s.append('<line x1="440" x2="440" y1="76" y2="106" stroke="var(--cross)" stroke-width="1.5" stroke-dasharray="3 3"/>')
    return ('<svg viewBox="0 0 840 160" width="100%" role="img" aria-label="An SD block made of duplicons, with genes riding on '
            'them and a fusion gene crossing a duplicon boundary">' + "".join(s) + '</svg>')


def famcn_diagram():
    """Why famCN tracks duplicons: read depth counts every copy of a gene's sequence, so genes on the same duplicons share
    a famCN and a gene on a more-copied duplicon gets a higher one. Toy numbers, diploid copies per genome."""
    s = []
    dup = {"D1": ("var(--s1)", 4), "D2": ("var(--s2)", 10), "D3": ("var(--s3)", 8)}
    s.append('<text x="0" y="14" font-size="12.5" font-weight="650" fill="var(--fg)">Copies of each duplicon in one genome (both parents)</text>')
    for gx, name in ((0, "D1"), (190, "D3"), (400, "D2")):
        col, n = dup[name]
        s.append(f'<text x="{gx}" y="40" font-size="12" font-weight="650" fill="var(--fg)">{name}</text>')
        for k in range(n):
            s.append(f'<rect x="{gx + 24 + k * 15}" y="28" width="12" height="14" rx="2" fill="{col}" fill-opacity=".8"/>')
        s.append(f'<text x="{gx + 24 + n * 15 + 4}" y="40" font-size="12" fill="var(--muted)">× {n}</text>')
    s.append('<text x="280" y="76" font-size="11" fill="var(--muted)" text-anchor="middle">read depth along the gene (each copy counts)</text>')
    s.append('<text x="458" y="76" font-size="11" fill="var(--muted)">famCN = average depth</text>')
    s.append('<text x="612" y="76" font-size="11" fill="var(--muted)">Soto\'s rule</text>')
    rows = [("gene X", "D1"), ("gene Z", "D1"), ("gene Y", "D2")]
    RH, Y0, sc = 78, 86, 3.2
    mids = []
    for r, (name, other) in enumerate(rows):
        y0 = Y0 + r * RH
        yb = y0 + 40
        d3, d2 = dup["D3"][1], dup[other][1]
        fam = (d3 + d2) / 2
        s.append(f'<text x="0" y="{y0 + 30}" font-size="12.5" font-weight="650" fill="var(--fg)">{name}</text>')
        s.append(f'<text x="0" y="{y0 + 45}" font-size="10.5" fill="var(--muted)">on D3 + {other}</text>')
        for (a, b, depth) in ((120, 280, d3), (280, 440, d2)):
            h = depth * sc
            s.append(f'<rect x="{a}" y="{yb - h}" width="{b - a}" height="{h}" fill="var(--faint)" fill-opacity=".35"/>')
            s.append(f'<text x="{(a + b) / 2}" y="{yb - h - 3}" font-size="10.5" fill="var(--muted)" text-anchor="middle">{depth}</text>')
        s.append(f'<line x1="120" x2="440" y1="{yb + 10}" y2="{yb + 10}" stroke="var(--igv)" stroke-width="1.5"/>')
        for a, b in ((132, 160), (205, 235), (300, 330), (380, 420)):
            s.append(f'<rect x="{a}" y="{yb + 5}" width="{b - a}" height="10" fill="var(--igv)"/>')
        s.append(f'<rect x="120" y="{yb + 19}" width="159" height="6" fill="{dup["D3"][0]}" fill-opacity=".8"/>')
        s.append(f'<rect x="281" y="{yb + 19}" width="159" height="6" fill="{dup[other][0]}" fill-opacity=".8"/>')
        s.append(f'<text x="458" y="{yb + 2}" font-size="11.5" fill="var(--muted)">({d3} + {d2}) ÷ 2 =</text>')
        s.append(f'<text x="560" y="{yb + 4}" font-size="17" font-weight="700" fill="var(--fg)">{fam:g}</text>')
        mids.append(yb)
    ytop, ybot = mids[0] + 15, mids[-1] + 5
    s.append(f'<line x1="220" x2="220" y1="{ytop}" y2="{ybot}" stroke="var(--cross)" stroke-width="1.5" stroke-dasharray="4 3"/>')
    s.append(f'<text x="220" y="{mids[-1] + 42}" font-size="11" fill="var(--cross)" text-anchor="middle">same exon (≥ 98%) in all three: one sequence family</text>')
    s.append(f'<path d="M600 {mids[0] - 8} H606 V{mids[1] + 4} H600" fill="none" stroke="var(--fg)" stroke-width="1.5"/>')
    yt = (mids[0] + mids[1]) / 2
    s.append(f'<text x="614" y="{yt - 10}" font-size="12.5" font-weight="650" fill="var(--fg)">Soto family 1</text>')
    s.append(f'<text x="614" y="{yt + 5}" font-size="11" fill="var(--muted)">famCN 6 and 6</text>')
    s.append(f'<text x="614" y="{yt + 19}" font-size="11" fill="var(--muted)">differ by &lt; 2: together</text>')
    s.append(f'<text x="614" y="{mids[2] - 10}" font-size="12.5" font-weight="650" fill="var(--fg)">Soto family 2</text>')
    s.append(f'<text x="614" y="{mids[2] + 5}" font-size="11" fill="var(--muted)">famCN 9 vs 6</text>')
    s.append(f'<text x="614" y="{mids[2] + 19}" font-size="11" fill="var(--cross)">differ by 3: cut here</text>')
    H = mids[-1] + 52
    return (f'<svg viewBox="0 0 760 {H}" width="100%" role="img" aria-label="Three genes share an exon on duplicon D3; their read '
            f'depth also covers D1 or D2, which have different copy numbers, so their famCN differ and Soto\'s cut separates the '
            f'genes on D2 from the genes on D1">' + "".join(s) + '</svg>')


def dotplot():
    rows = sorted(CLUSTERS, key=lambda r: -r[1])
    W, L, R, T, rh = 760, 260, 20, 26, 15
    H = T + rh * len(rows) + 40
    lo, hi = -0.3, 0.9
    X = lambda v: L + (v - lo) / (hi - lo) * (W - L - R)
    s = []
    for v in (-0.2, 0, 0.2, 0.4, 0.6, 0.8):
        s.append(f'<line x1="{X(v):.1f}" x2="{X(v):.1f}" y1="{T - 6}" y2="{H - 34}" stroke="{"var(--fg)" if v == 0 else "var(--line)"}"/>')
        s.append(f'<text x="{X(v):.1f}" y="{H - 20}" font-size="11" fill="var(--muted)" text-anchor="middle">{v:+.1f}</text>')
    s.append(f'<text x="{(L + W - R) / 2}" y="{H - 4}" font-size="11.5" fill="var(--fg)" text-anchor="middle">duplicon sharing inside Soto\'s families minus across them (0 = no relation)</text>')
    s.append(f'<text x="{X(0) + 6}" y="{T - 10}" font-size="10.5" fill="var(--muted)">no relation</text>')
    for k, (half, dl, fams, ex) in enumerate(rows):
        y = T + rh * k + rh / 2
        names = ", ".join(ex.split(", ")[:2])
        bad = dl <= 0.05
        s.append(f'<text x="{L - 10}" y="{y + 4}" font-size="10.5" fill="{"var(--cross)" if bad else "var(--fg)"}" text-anchor="end">{names}</text>')
        fill = "var(--accent)" if half == "heldout" else "var(--surface)"
        s.append(f'<circle cx="{X(dl):.1f}" cy="{y}" r="5" fill="{fill}" stroke="var(--accent)" stroke-width="1.8"><title>{fams}: {dl:+.3f} ({"held-out" if half == "heldout" else "dev"})</title></circle>')
    return (f'<svg viewBox="0 0 {W} {H}" width="100%" role="img" aria-label="Per sequence family, how much more duplicon '
            f'content genes share inside Soto\'s families than across them">' + "".join(s) + '</svg>')


def pane():
    rows = "".join(f'<tr><th scope="row">{a}</th><td>{b}</td><td>{c}</td></tr>' for a, b, c in ROWS)
    quotes = "".join(f'<blockquote>"{q}"<cite>Soto et al. 2025, Cell 188:5363, {w}</cite></blockquote>' for q, w in QUOTES)
    return f'''<div role="tabpanel" id="pane-def" aria-labelledby="tab-def" class="pane" hidden>
<section class="how-intro">
  <div class="eyebrow">Two meanings of "gene family" · Soto 2025 and the evolutionary definition</div>
  <h2>Same word, different scope</h2>
  <p class="note">Soto define gene families operationally for their question, and state the definition explicitly. It is
  narrower than the evolutionary definition: it sees only recent duplications and splits homologs by copy number. Where both
  apply, they agree.</p>
</section>
<section>
  <div class="def-grid">
    <div class="chart">{nested()}</div>
    <div class="cn-facts">
      <div class="tile"><div class="lab">Soto families inside one sequence family</div><div class="big">440 / 444</div><div class="sub">99.1%; the exceptions include two families Soto merged by hand</div></div>
      <div class="tile"><div class="lab">Sequence families cut by copy number</div><div class="big">33</div><div class="sub">holding 87 Soto families (NPIP: 4)</div></div>
      <div class="tile"><div class="lab">Paralog pairs Soto's table keeps apart</div><div class="big">83 pairs</div><div class="sub">median identity 0.858: 65 have no ≥ 98% exon link (outside the SD98 scope), 18 are linked but split by copy number</div></div>
    </div>
  </div>
</section>
<section>
  <h2>Which duplications each definition sees</h2>
  <div class="chart">{axis()}</div>
</section>
<section>
  <h2>Side by side</h2>
  <div class="tbl"><table class="deftbl">
    <thead><tr><th></th><th>Soto 2025 (operational)</th><th>Evolutionary (this thesis)</th></tr></thead>
    <tbody>{rows}</tbody>
  </table></div>
</section>
<section>
  <h2>Why Soto's copy-number cut is not arbitrary: it follows duplicons</h2>
  <p class="note">Segmental duplications are mosaics of duplicons, ancestral duplication units. Three objects, not one:</p>
  <div class="chart">{three_objects()}</div>
  <h3 class="sub-h">Why copy number tracks duplicons</h3>
  <div class="chart">{famcn_diagram()}</div>
  <p class="note">famCN is read depth with every read counted at each copy it matches, so it counts how many copies of a gene's
  sequence one genome carries. Genes on the same duplicons get the same famCN; a gene on a more-copied duplicon gets a higher one.
  Soto keep two genes together only if their famCN differ by less than 2, so the cut falls where the duplicon make-up changes
  (toy numbers above). That follows from how famCN is measured; the test below checks that it held in Soto's table.</p>
  <p class="note">Test, pre-registered before any number was computed: inside each of the 33 sequence families Soto cuts by copy
  number, do genes in the same Soto family share more duplicon content than genes in different ones? Duplicons from Vollger et al.
  2022 (the DupMasker annotation Soto cite).</p>
  <div class="def-grid">
    <div class="chart">{dotplot()}<div class="legend" style="margin-top:6px"><span><svg width="12" height="12" aria-hidden="true"><circle cx="6" cy="6" r="5" fill="var(--accent)"/></svg> held-out half (decides)</span><span><svg width="12" height="12" aria-hidden="true"><circle cx="6" cy="6" r="4.5" fill="none" stroke="var(--accent)" stroke-width="1.8"/></svg> development half</span><span style="color:var(--cross)">red = cut inside one duplicon make-up</span></div></div>
    <div class="cn-facts">
      <div class="tile"><div class="lab">Soto's pieces sit on different duplicons</div><div class="big">14 / 16</div><div class="sub">held-out sequence families (p = 0.0001); development half 15 / 17</div></div>
      <div class="tile"><div class="lab">Copy-number gaps follow duplicon differences</div><div class="big">ρ = 0.58</div><div class="sub">8,565 gene pairs; 0.31 across Soto's cuts alone</div></div>
      <div class="tile"><div class="lab">Exceptions</div><div class="sub">6 cuts inside one duplicon make-up, among them TBC1D3 and CICP (red in the plot). Not tested: whether every duplicon boundary is a family boundary.</div></div>
    </div>
  </div>
  <p class="note">Reading: different duplicons carry different copy numbers, so Soto's copy-number gate largely separates genes riding on different duplicons. Duplicons are the shared unit: a family is a set of duplicon-aligned pieces, an SD block can carry several families, and a fusion gene crosses a duplicon boundary.</p>
</section>
<section>
  <h2>In Soto's own words</h2>
  <div class="quotes">{quotes}</div>
</section>
<section>
  <h2>One example: NPIPB3, NPIPB4, NPIPB5</h2>
  <div class="tiles">
    <div class="tile"><div class="lab">Homology family</div><div class="big">one</div><div class="sub">all three are NPIP copies</div></div>
    <div class="tile"><div class="lab">Sequence family</div><div class="big">one</div><div class="sub">NPIPB3 and NPIPB4 are 98.6% identical</div></div>
    <div class="tile"><div class="lab">Soto families</div><div class="big">three</div><div class="sub">copy numbers 26.5, 29.1 and 22.2 differ by 2 or more</div></div>
  </div>
  <div class="bar"><button class="chip" data-real="ID_151">Open NPIPB3 in the browser</button><button class="chip" data-real="ID_154">Open the NPIP A-clade family</button></div>
</section>
</div>'''


if __name__ == "__main__":
    print(pane())
