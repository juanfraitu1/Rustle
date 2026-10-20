#!/usr/bin/env python3
"""The gorilla tab of the NPIP Read Pools page: STON1 on KB3781's mother's chromosome 12 and not his father's, from
ston1_fig.py's data (user 2026-10-20: "can we make another tab in this with those data from gorilla?").

    gorilla_tab_render.py STON1.json OUT.html
"""
import json
import sys


def glyph_exons(exons, x0, sc, y, h, strand, cls):
    """A gene as its exons on a thin intron line, with a strand point on the last exon."""
    o = [f'<line class="intr" x1="{x0 + exons[0][0] * sc:.1f}" y1="{y + h / 2:.1f}" x2="{x0 + exons[-1][1] * sc:.1f}" y2="{y + h / 2:.1f}"/>']
    for a, b in exons:
        o.append(f'<rect class="{cls}" x="{x0 + a * sc:.1f}" y="{y}" width="{max(1.2, (b - a) * sc):.1f}" height="{h}"/>')
    tipx = x0 + (exons[-1][1] if strand == "+" else exons[0][0]) * sc
    d = 5 if strand == "+" else -5
    o.append(f'<polygon class="{cls}" points="{tipx:.1f},{y - 2} {tipx + d:.1f},{y + h / 2:.1f} {tipx:.1f},{y + h + 2}"/>')
    return "".join(o)


def map_svg(d):
    W, X0, SPAN = 960, 128, 800
    L = 600000
    sc = SPAN / L
    ym, yp = 70, 168
    m, p = d["mother"], d["father"]
    o = [f'<svg viewBox="0 0 {W} 214" role="img" aria-label="KB3781 mother against father around STON1">']
    for r in d["ribbons"]:
        o.append(f'<polygon class="rib gk" points="{X0 + r["q0"] * sc:.1f},{ym + 12} {X0 + r["q1"] * sc:.1f},{ym + 12} '
                 f'{X0 + r["t1"] * sc:.1f},{yp - 4} {X0 + r["t0"] * sc:.1f},{yp - 4}"/>')
    gap0 = max(r["q1"] for r in d["ribbons"] if r["q0"] == 0)
    gap1 = min(r["q0"] for r in d["ribbons"] if r["q0"] > 0)
    o.append(f'<rect class="only" x="{X0 + gap0 * sc:.1f}" y="{ym - 34}" width="{(gap1 - gap0) * sc:.1f}" height="{ym - 18 + 34 - (ym - 34) + 16}"/>')
    o.append(f'<text class="onl" x="{X0 + (gap0 + gap1) / 2 * sc:.1f}" y="{ym - 38}" text-anchor="middle">'
             f'{(gap1 - gap0) / 1000:.1f} kb only the mother has</text>')
    # coverage of the 113 reads on the mother (per 2 kb)
    bins = [0] * (L // 2000 + 1)
    for s, e, _ in d["reads"]["mother"]["in_window"]:
        for b in range(s // 2000, e // 2000 + 1):
            bins[b] += 1
    mx = max(bins) or 1
    pts = [f"{X0:.1f},{ym - 4}"] + [f"{X0 + (i * 2000 + 1000) * sc:.1f},{ym - 4 - 24 * v / mx:.1f}" for i, v in enumerate(bins)] + [f"{X0 + SPAN:.1f},{ym - 4}"]
    o.append(f'<polygon class="cov" points="{" ".join(pts)}"/>')
    o.append(f'<text class="tm" x="{X0 + 330000 * sc:.1f}" y="{ym - 20}">113 reads</text>')
    o.append(f'<line class="axisl" x1="{X0}" y1="{ym + 6}" x2="{X0 + SPAN}" y2="{ym + 6}"/>')
    o.append(f'<line class="axisl" x1="{X0}" y1="{yp + 6}" x2="{X0 + SPAN}" y2="{yp + 6}"/>')
    for side, y in ((m, ym), (p, yp)):
        for name, g in side["genes"].items():
            if g["qcov"] < 0.2:
                continue                                   # a stray piece (STON1 on the father: 5% of the mRNA)
            cls = "gs" if name == "STON1" else "gg"
            o.append(glyph_exons(g["exons"], X0, sc, y, 12, g["strand"], cls))
            if g["qcov"] > 0.9:
                o.append(f'<text class="gl" x="{X0 + (g["exons"][0][0] + g["exons"][-1][1]) / 2 * sc:.1f}" y="{y + 30}" text-anchor="middle">{name}</text>')
            else:                                          # a partial gene: label to its right, clear of its neighbour's
                o.append(f'<text class="gl" x="{X0 + g["exons"][-1][1] * sc + 8:.1f}" y="{y + 10}">{name}, {100 * g["qcov"]:.0f}% of the mRNA</text>')
    o.append(f'<text class="tn" x="8" y="{ym + 10}">Mother</text><text class="tm" x="8" y="{ym + 24}">chr12 (CM054594.2)</text>')
    o.append(f'<text class="tn" x="8" y="{yp + 10}">Father</text><text class="tm" x="8" y="{yp + 24}">chr12 = reference</text>')
    for t in range(0, L + 1, 100000):
        o.append(f'<text class="tk" x="{X0 + t * sc:.1f}" y="{yp + 44}" text-anchor="middle">{t // 1000} kb</text>')
    o.append("</svg>")
    return "".join(o)


def dot_svg(d, thr=150, full=500, binsz=2000):
    L = 600000
    PW = 300
    sc = PW / L
    Lm, T, B = 40, 10, 40
    W, H = Lm + PW + 10, T + PW + B
    o = [f'<svg viewBox="0 0 {W} {H}" role="img" aria-label="Dot plot, mother (x) against father (y), around STON1">',
         f'<rect class="plotbg" x="{Lm}" y="{T}" width="{PW}" height="{PW}"/>']
    cell = binsz * sc
    for key, (fw, rv) in d["dots"].items():
        xb, yb = map(int, key.split(","))
        for n, cls in ((fw, "df"), (rv, "dr")):
            if n >= thr:
                o.append(f'<rect class="{cls}" x="{Lm + xb * cell:.1f}" y="{T + yb * cell:.1f}" width="{cell + 0.4:.2f}" '
                         f'height="{cell + 0.4:.2f}" opacity="{min(1.0, 0.35 + 0.65 * n / full):.2f}"/>')
    g = d["mother"]["genes"]["STON1"]
    o.append(f'<rect class="gs" x="{Lm + g["exons"][0][0] * sc:.1f}" y="{T + PW + 3}" width="{(g["exons"][-1][1] - g["exons"][0][0]) * sc:.1f}" height="5"/>')
    for t in range(0, L + 1, 200000):
        o.append(f'<text class="tk" x="{Lm + t * sc:.1f}" y="{T + PW + 20}" text-anchor="middle">{t // 1000}</text>'
                 f'<text class="tk" x="{Lm - 5}" y="{T + t * sc + 3:.1f}" text-anchor="end">{t // 1000}</text>')
    o.append(f'<text class="tk" x="{Lm + PW / 2}" y="{H - 6}" text-anchor="middle">mother (kb); father on the y axis</text></svg>')
    return "".join(o)


STYLE = """
<style>
/* gorilla tab: the page's tokens plus a STON1 accent */
:root{--gs:#c2410c;--gk:#4B4F58;--onlybg:#FCEED4}
@media (prefers-color-scheme: dark){:root:not([data-theme="light"]){--gs:#fb923c;--gk:#C2C5CE;--onlybg:#3A2C14}}
:root[data-theme="dark"]{--gs:#fb923c;--gk:#C2C5CE;--onlybg:#3A2C14}
#gorilla svg text{font-family:var(--mono)}
#gorilla .gs{fill:var(--gs)} #gorilla .gg{fill:var(--ink-soft)} #gorilla .intr{stroke:var(--ink-soft);stroke-width:1}
#gorilla .rib.gk{fill:var(--gk);opacity:.16}
#gorilla .only{fill:var(--onlybg)} #gorilla .onl{font-size:11px;fill:var(--warn);font-weight:600}
#gorilla .cov{fill:var(--gs);opacity:.55}
#gorilla .axisl{stroke:var(--rule-strong);stroke-width:1.2}
#gorilla .gl{font-size:10.5px;fill:var(--ink-soft)} #gorilla .tn{font-size:12px;font-weight:600;fill:var(--ink)}
#gorilla .tm{font-size:11px;fill:var(--muted)} #gorilla .tk{font-size:9.5px;fill:var(--muted)}
#gorilla .plotbg{fill:var(--plotbg);stroke:var(--rule)} #gorilla .df{fill:var(--dotf)} #gorilla .dr{fill:var(--dotr)}
#gorilla .two{display:grid;grid-template-columns:minmax(0,360px) minmax(0,1fr);gap:16px;align-items:start}
@media (max-width:760px){#gorilla .two{grid-template-columns:1fr}}
#gorilla .units{display:grid;grid-template-columns:repeat(auto-fit,minmax(170px,1fr));gap:10px;margin:14px 0 4px}
#gorilla .unit{background:var(--surface);border:1px solid var(--rule);border-radius:10px;padding:10px 12px;min-width:0}
#gorilla .unit .n{font-family:var(--mono);font-size:11px;letter-spacing:.05em;text-transform:uppercase;color:var(--muted)}
#gorilla .unit .v{font-family:var(--mono);font-size:20px;font-weight:600;font-variant-numeric:tabular-nums}
#gorilla .unit .s{font-size:12.5px;color:var(--ink-soft)}
#gorilla table.fate td,#gorilla table.fate th{white-space:normal;text-align:left}
</style>
"""


def main(json_p, out_p):
    d = json.load(open(json_p))
    rm, rp = d["reads"]["mother"], d["reads"]["father"]
    gap0 = max(r["q1"] for r in d["ribbons"] if r["q0"] == 0)
    gap1 = min(r["q0"] for r in d["ribbons"] if r["q0"] > 0)
    t0 = max(r["t1"] for r in d["ribbons"] if r["q0"] == 0)
    t1 = min(r["t0"] for r in d["ribbons"] if r["q0"] > 0)
    flank_id = min(r["identity"] for r in d["ribbons"])
    s1 = d["mother"]["genes"]["STON1"]
    html = f"""{STYLE}
<div id="gorilla">
  <p class="eyebrow">Rustle · O3 · gorilla chr12 · KB3781 ("Jim") · his own two parental assemblies (mGorGor1 v2.0) and his own Iso-Seq and DNA</p>
  <h1>STON1 is on Jim's mother's chromosome 12 and not on his father's</h1>
  <p class="lede">KB3781 is the gorilla the T2T reference mGorGor1 was assembled from, so both of his parental chromosomes are
  known base by base. The reference sequence of chromosome 12 (NC_073236.2) is his father's. His mother's chromosome carries
  STON1; at the same place his father's carries a 3 kb join. His STON1 reads therefore have no place on the reference.</p>
  <div class="units">
    <div class="unit"><div class="n">Mother, between the shared flanks</div><div class="v">{(gap1 - gap0) / 1000:.1f} kb</div><div class="s">LHCGR … GTF2A1L, STON1</div></div>
    <div class="unit"><div class="n">Father, same place</div><div class="v">{(t1 - t0) / 1000:.1f} kb</div><div class="s">STON1 and the 5′ end of GTF2A1L are gone</div></div>
    <div class="unit"><div class="n">STON1 mRNA aligned</div><div class="v">{100 * s1["qcov"]:.1f}% / none</div><div class="s">human NM_006873 on the mother / on the father</div></div>
    <div class="unit"><div class="n">Jim's STON1 Iso-Seq reads</div><div class="v">{rm["placed"]} / {rp["placed"]}</div><div class="s">align to the mother's genome / to the father's</div></div>
    <div class="unit"><div class="n">Jim's DNA (Illumina WGS)</div><div class="v">97.1%</div><div class="s">of STON1 covered, median depth 6.0 = one chromosome (λ 4.98)</div></div>
  </div>

  <section>
    <h2>The mother's and father's chromosome 12 around STON1</h2>
    <p class="sub">600 kb of each, aligned to each other (minimap2 -x asm5): the grey ribbons are the shared sequence, {100 * flank_id:.2f}%
    identical or more on both sides. Genes are human mRNAs spliced onto each window. The orange area is the coverage of Jim's 113
    STON1 reads on his mother's chromosome.</p>
    <div class="win wscroll">{map_svg(d)}</div>
  </section>

  <section>
    <h2>Dot plot: the diagonal breaks where only the mother has sequence</h2>
    <div class="two">
      <div class="card">{dot_svg(d)}</div>
      <div class="card"><h3>How to read it</h3><ul>
        <li>Every point is an exact 32-bp match between the mother (x) and the father (y). Where the two chromosomes are the same,
        the points form one diagonal.</li>
        <li>Around STON1 (orange bar under the x axis) the diagonal stops: about {(gap1 - gap0) / 1000:.0f} kb of the mother match
        nothing in the father. It resumes only {(t1 - t0) / 1000:.1f} kb further along the father, so the two segments sit side by
        side with a horizontal gap between them.</li>
        <li>That is what a deletion on one chromosome looks like. A misassembly would not leave both sides
        colinear at {100 * flank_id:.2f}% identity.</li>
      </ul></div>
    </div>
  </section>

  <section>
    <h2>Where Jim's STON1 reads go</h2>
    <div class="scroll"><table class="fate">
      <tr><th>aligned to (minimap2 -x splice:hq, whole genome)</th><th>reads placed</th><th>where</th></tr>
      <tr><td class="lab">the reference (= his father's chromosomes)</td><td class="hi">0 of 113</td><td class="lab">unaligned: these reads were found among the 959 reads that align nowhere on it</td></tr>
      <tr><td class="lab">his father's whole assembly</td><td class="hi">{rp["placed"]} of {rp["placed"] + rp["unmapped"]}</td><td class="lab">unaligned</td></tr>
      <tr><td class="lab">his mother's whole assembly</td><td class="hi">{rm["placed"]} of {rm["placed"] + rm["unmapped"]}</td><td class="lab">all on chromosome 12, on STON1</td></tr>
    </table></div>
    <p class="sub">STON1 has no close relative in the gorilla genome, so on the father's sequence its reads do not even tie
    elsewhere: they go unaligned. The human tab's LOC124907808 / LOC124907807 are the other case, a copy with near-identical
    relatives, whose reads tie to NPIPB15 instead. Same event, a chromosome lacking a gene its partner has, two outcomes.</p>
  </section>

  <section>
    <h2>The 16q22 NPIP repeat in Jim: the same count on both parents</h2>
    <p class="sub">The orthologous region holds five NPIPB15-like copies on each of Jim's chromosomes 16 (human NPIPB15 as the probe,
    95-96% across species), two inverted and three forward, in the same order. His mother's and father's 7 Mb align colinearly
    at 99.4-99.75% identity. This gorilla does not vary there, so his NPIP reads behave as on a reference that matches him.</p>
  </section>

  <div class="card checks" style="margin-top:20px"><h3>Why this is not an assembly or alignment error</h3><ul>
    <li><b>One animal throughout.</b> The two chromosomes, the RNA and the DNA are all Jim's (KB3781): the T2T assemblies of his
    parental haplotypes, his fibroblast Iso-Seq, his Illumina whole-genome reads.</li>
    <li><b>His DNA confirms it.</b> 97.1% of the STON1 sequence his reads build is covered by his own whole-genome k-mers, at a
    median depth of 6.0, the depth of one chromosome (λ = 4.98 per haploid genome): one copy, as the assemblies say.</li>
    <li><b>Anchored both sides.</b> Outside the deleted stretch the two chromosomes are colinear at {100 * flank_id:.2f}% or more.</li>
    <li><b>Whole-genome, not windowed.</b> The reads were aligned to each parent's entire genome: {rm["placed"]} of 113 on the mother,
    {rp["placed"]} on the father, so nothing on the father's chromosomes resembles them closely enough to align.</li>
  </ul></div>
  <p class="sub" style="margin-top:14px">Data: <code>bench/o3_graph/ston1_fig.py</code>, <code>gorilla_tab_render.py</code>; mother
  CM054594.2:95,700,000-96,300,000 (GCA_028885495.2), father CM054570.2:90,900,000-91,500,000 (GCA_028885475.2; reference
  NC_073236.2 has the same length and sequence); the 113 reads of O3 cluster <code>clSRR27438212.9740719</code>
  (<code>docs/PREREG_unmapped_rescue_2026-10-08.md</code>, Amendments 28-34, DNA check §T1). Gorilla numbers are never pooled
  with the human tab's. Built 20 Oct 2026.</p>
</div>
"""
    open(out_p, "w").write(html)
    print(f"[gorilla_tab] -> {out_p} ({len(html):,} bytes); mother-only {gap1 - gap0} bp, father {t1 - t0} bp, flank id {flank_id}")
    return 0


if __name__ == "__main__":
    sys.exit(main(*sys.argv[1:3]))
