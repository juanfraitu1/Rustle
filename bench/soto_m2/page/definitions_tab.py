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
