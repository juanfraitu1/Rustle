#!/usr/bin/env python3
"""Render npip_cnv_fig.py's data.json as the copy-number section of the NPIP Read Pools page (inline SVG, the page's tokens).

    npip_cnv_render.py DATA.json READS.json OUT.html
READS.json: the A119b primaries on LOC124907808 / LOC124907807 by margin to the best other NPIP copy (exact / near / clear).
"""
import json
import sys

UNIT_CLASS = {0: "uf", 1: "u1", 2: "u2", 3: "u3"}


def kb(x):
    return f"{x / 1000:.1f}"


def glyph(x0, x1, y, h, strand, cls):
    """A gene box with a strand point."""
    w = max(2.0, x1 - x0)
    tip = min(6.0, w / 3)
    if strand == "+":
        pts = [(x0, y), (x0 + w - tip, y), (x0 + w, y + h / 2), (x0 + w - tip, y + h), (x0, y + h)]
    else:
        pts = [(x0 + tip, y), (x0 + w, y), (x0 + w, y + h), (x0 + tip, y + h), (x0, y + h / 2)]
    return f'<polygon class="g {cls}" points="{" ".join(f"{a:.1f},{b:.1f}" for a, b in pts)}"/>'


def map_panel(d, h):
    W, X0, SPAN = 960, 128, 800
    sc = SPAN / (d["region"]["end"] - d["region"]["start"])
    yt, yb = 34, 112
    o = [f'<svg viewBox="0 0 {W} 160" role="img" aria-label="CHM13 against {h["name"]}: each CHM13 unit and where it aligns">']
    for r in h["ribbons"]:
        a, b, c_, e = X0 + r["q0"] * sc, X0 + r["q1"] * sc, X0 + r["t0"] * sc, X0 + r["t1"] * sc
        cls = UNIT_CLASS[r["unit"]] + (" sec" if not r["primary"] else "")
        o.append(f'<polygon class="rib {cls}" points="{a:.1f},{yt + 9} {b:.1f},{yt + 9} {e:.1f},{yb - 3} {c_:.1f},{yb - 3}"/>')
    o.append(f'<line class="axisl" x1="{X0}" y1="{yt + 4}" x2="{X0 + SPAN}" y2="{yt + 4}"/>')
    for g in d["genes"]:
        o.append(glyph(X0 + g["start"] * sc, X0 + g["end"] * sc, yt - 2, 12, g["strand"], UNIT_CLASS[g["unit"]]))
        if g["cls"] == "NPIP":
            o.append(f'<text class="gl" x="{X0 + (g["start"] + g["end"]) / 2 * sc:.1f}" y="{yt - 8}" text-anchor="middle">{g["name"]}</text>')
    hl = h["length"] * sc
    o.append(f'<line class="axisl" x1="{X0}" y1="{yb + 4}" x2="{X0 + hl:.1f}" y2="{yb + 4}"/>')
    gstrand = {g["name"]: g["strand"] for g in d["genes"]}
    for c in h["copies"]:                    # the gene's own strand: CHM13's, flipped when its alignment is reverse
        st = gstrand[c["best_gene"]] if c["strand"] == "+" else ("-" if gstrand[c["best_gene"]] == "+" else "+")
        cls = UNIT_CLASS[c["unit"]] if c["cls"] == "NPIP" else "uf"     # only the NPIP copy names its unit (see the checks)
        o.append(glyph(X0 + c["ts"] * sc, X0 + c["te"] * sc, yb - 2, 12, st, cls))
        if c["cls"] == "NPIP":
            nxt = sorted(c["others"].values(), reverse=True)
            cx = X0 + (c["ts"] + c["te"]) / 2 * sc
            o.append(f'<text class="gl" x="{cx:.1f}" y="{yb + 24}" text-anchor="middle">like {c["best_gene"]} {100 * c["identity"]:.2f}%</text>')
            if len(nxt) > 1:
                o.append(f'<text class="tm" x="{cx:.1f}" y="{yb + 37}" text-anchor="middle">next best {100 * nxt[1]:.2f}%</text>')
    units = sum(1 for c in h["copies"] if c["cls"] == "NPIP")
    o.append(f'<text class="tn" x="8" y="{yt + 8}">CHM13</text><text class="tm" x="8" y="{yt + 22}">450.0 kb · 3 units</text>')
    o.append(f'<text class="tn" x="8" y="{yb + 8}">{h["name"]}</text>'
             f'<text class="tm" x="8" y="{yb + 22}">{kb(h["syntenic_len"])} kb · {units} unit{"s" if units > 1 else ""}</text>')
    o.append("</svg>")
    return "".join(o)


def dot_panel(d, title, sub, ylen, dots, binsz=2000, thr=150, full=500):
    X = d["region"]["end"] - d["region"]["start"]
    PW = 300
    sc = PW / X
    PH = ylen * sc
    L, T, B = 40, 10, 44
    W, H = L + PW + 10, T + PH + B
    o = [f'<svg viewBox="0 0 {W:.0f} {H:.0f}" role="img" aria-label="Dot plot: CHM13 (x) against {title} (y)">',
         f'<rect class="plotbg" x="{L}" y="{T}" width="{PW}" height="{PH:.1f}"/>']
    for g in d["genes"]:                                        # unit bands under the x axis
        o.append(f'<rect class="band {UNIT_CLASS[g["unit"]]}" x="{L + g["start"] * sc:.1f}" y="{T + PH + 3:.1f}" '
                 f'width="{max(1.0, (g["end"] - g["start"]) * sc):.1f}" height="5"/>')
    cell = binsz * sc
    for key, (fw, rv) in dots.items():
        xb, yb = map(int, key.split(","))
        for n, cls in ((fw, "df"), (rv, "dr")):
            if n >= thr:
                op = min(1.0, 0.35 + 0.65 * n / full)
                o.append(f'<rect class="{cls}" x="{L + xb * cell:.1f}" y="{T + yb * cell:.1f}" width="{cell + 0.4:.2f}" '
                         f'height="{cell + 0.4:.2f}" opacity="{op:.2f}"/>')
    for t in range(0, X + 1, 100000):
        x = L + t * sc
        o.append(f'<line class="tick" x1="{x:.1f}" y1="{T + PH + 9:.1f}" x2="{x:.1f}" y2="{T + PH + 12:.1f}"/>'
                 f'<text class="tk" x="{x:.1f}" y="{T + PH + 22:.1f}" text-anchor="middle">{t // 1000}</text>')
    for t in range(0, int(ylen) + 1, 100000):
        y = T + t * sc
        o.append(f'<text class="tk" x="{L - 5}" y="{y + 3:.1f}" text-anchor="end">{t // 1000}</text>')
    o.append(f'<text class="tk" x="{L + PW / 2:.1f}" y="{H - 6:.1f}" text-anchor="middle">CHM13 (kb)</text>')
    o.append("</svg>")
    return f'<figure class="dp"><figcaption><b>{title}</b><span>{sub}</span></figcaption>{"".join(o)}</figure>'


def reads_panel(rd):
    W, X0, SPAN = 960, 150, 700
    tot = max(v["total"] for v in rd.values())
    sc = SPAN / tot
    o = [f'<svg viewBox="0 0 {W} 112" role="img" aria-label="A119b primaries on the two copies by AS margin to another NPIP copy">']
    for i, (name, v) in enumerate(rd.items()):
        y = 14 + i * 46
        x = X0
        o.append(f'<text class="tn" x="8" y="{y + 13}">{name}</text><text class="tm" x="8" y="{y + 27}">{v["total"]} primaries</text>')
        for k, cls, lab in (("exact", "rx", "exact tie"), ("near", "rn", "near tie"), ("clear", "rc", "placed")):
            w = v[k] * sc
            o.append(f'<rect class="{cls}" x="{x:.1f}" y="{y}" width="{w:.1f}" height="20"/>')
            if w > 60:
                o.append(f'<text class="bl" x="{x + 6:.1f}" y="{y + 14}">{v[k]} {lab}</text>')
            x += w
        o.append(f'<text class="tm" x="{X0}" y="{y + 34}">also aligned within 10 AS to: {v["partners"]}</text>')
    o.append("</svg>")
    return "".join(o)


STYLE = """
<style>
/* copy-number section: unit colours and dot-plot marks, same token pattern as the page */
:root{--u1:#2f6fbf;--u2:#b8741a;--u3:#7d4fbf;--uf:#a3a8b0;--dotf:#23262c;--dotr:#c2410c;--plotbg:#F6F7F4}
@media (prefers-color-scheme: dark){:root:not([data-theme="light"]){--u1:#5b95e0;--u2:#e0a04a;--u3:#a982e6;--uf:#6d727b;--dotf:#e6e8ee;--dotr:#fb923c;--plotbg:#202228}}
:root[data-theme="dark"]{--u1:#5b95e0;--u2:#e0a04a;--u3:#a982e6;--uf:#6d727b;--dotf:#e6e8ee;--dotr:#fb923c;--plotbg:#202228}
#cnv .lede2{color:var(--ink-soft);max-width:76ch}
#cnv .units{display:grid;grid-template-columns:repeat(auto-fit,minmax(150px,1fr));gap:10px;margin:14px 0 4px}
#cnv .unit{background:var(--surface);border:1px solid var(--rule);border-radius:10px;padding:10px 12px;min-width:0}
#cnv .unit .n{font-family:var(--mono);font-size:11px;letter-spacing:.05em;text-transform:uppercase;color:var(--muted)}
#cnv .unit .v{font-family:var(--mono);font-size:20px;font-weight:600;font-variant-numeric:tabular-nums}
#cnv .unit .s{font-size:12.5px;color:var(--ink-soft)}
#cnv svg text{font-family:var(--mono)}
#cnv .g{stroke:var(--surface);stroke-width:.6}
#cnv .u1{fill:var(--u1)} #cnv .u2{fill:var(--u2)} #cnv .u3{fill:var(--u3)} #cnv .uf{fill:var(--uf)}
#cnv .rib{opacity:.30;stroke:none} #cnv .rib.sec{opacity:.10}
#cnv .axisl{stroke:var(--rule-strong);stroke-width:1.2}
#cnv .gl{font-size:10.5px;fill:var(--ink-soft)}
#cnv .tn{font-size:12px;font-weight:600;fill:var(--ink)} #cnv .tm{font-size:11px;fill:var(--muted)}
#cnv .tk{font-size:9.5px;fill:var(--muted)} #cnv .tick{stroke:var(--rule-strong)}
#cnv .plotbg{fill:var(--plotbg);stroke:var(--rule)}
#cnv .band{opacity:.9}
#cnv .df{fill:var(--dotf)} #cnv .dr{fill:var(--dotr)}
#cnv .dps{display:grid;grid-template-columns:repeat(auto-fit,minmax(280px,1fr));gap:14px;align-items:start}
#cnv figure.dp{margin:0;background:var(--surface);border:1px solid var(--rule);border-radius:12px;padding:12px 12px 6px;min-width:0}
#cnv figure.dp figcaption{display:flex;flex-direction:column;gap:2px;font-size:12.5px;color:var(--ink-soft);margin-bottom:6px}
#cnv figure.dp figcaption b{font-family:var(--mono);font-size:12.5px;color:var(--ink)}
#cnv figure.dp svg{display:block;width:100%;height:auto}
#cnv .rx{fill:var(--tied)} #cnv .rn{fill:var(--warn)} #cnv .rc{fill:var(--good)}
#cnv .bl{font-size:11px;fill:var(--surface)}
#cnv .checks li{margin-bottom:7px}
#cnv .ulegend{display:flex;flex-wrap:wrap;gap:14px;font-family:var(--mono);font-size:11.5px;color:var(--muted);margin:8px 0 0}
#cnv .ulegend span{display:inline-flex;align-items:center;gap:6px}
#cnv .sw.u1,#cnv .sw.u2,#cnv .sw.u3,#cnv .sw.uf,#cnv .sw.df,#cnv .sw.dr{display:inline-block;width:14px;height:9px;border-radius:2px}
#cnv .sw.u1{background:var(--u1)} #cnv .sw.u2{background:var(--u2)} #cnv .sw.u3{background:var(--u3)} #cnv .sw.uf{background:var(--uf)}
#cnv .sw.df{background:var(--dotf)} #cnv .sw.dr{background:var(--dotr)}
#cnv .sw.rxs,#cnv .sw.rns,#cnv .sw.rcs{display:inline-block;width:14px;height:9px;border-radius:2px}
#cnv .sw.rxs{background:var(--tied)} #cnv .sw.rns{background:var(--warn)} #cnv .sw.rcs{background:var(--good)}
</style>
"""


def main(data_p, reads_p, out_p):
    d = json.load(open(data_p))
    rd = json.load(open(reads_p))
    period = sum(d["unit_period"]) / 2
    tiles = [('<div class="unit"><div class="n">CHM13 (T2T)</div><div class="v">450.0 kb</div>'
              '<div class="s">3 units: NPIPB15, LOC124907808, LOC124907807</div></div>')]
    for h in d["haps"]:
        miss = 450000 - h["syntenic_len"]
        units = sum(1 for c in h["copies"] if c["cls"] == "NPIP")
        tiles.append(f'<div class="unit"><div class="n">{h["name"]}</div><div class="v">{kb(h["syntenic_len"])} kb</div>'
                     f'<div class="s">{units} unit{"s" if units > 1 else ""}; {kb(miss)} kb shorter = {miss / period:.2f} units</div></div>')
    maps = "".join(f'<div class="win wscroll">{map_panel(d, h)}</div>' for h in d["haps"])
    dps = [dot_panel(d, "CHM13 against itself", "the diagonal plus parallel lines 114 and 229 kb off it: one unit, three times",
                     450000, d["self_dots"])]
    for h in d["haps"]:
        units = sum(1 for c in h["copies"] if c["cls"] == "NPIP")
        sub = ("all three CHM13 units fall on the same 100 kb" if units == 1
               else "units 2 and 3 fall on the same second unit")
        dps.append(dot_panel(d, h["name"], sub, h["length"], h["dots"]))
    html = f"""{STYLE}
  <section id="cnv">
    <h2>Why LOC124907808 and LOC124907807 rest on tied reads: people carry fewer copies of this repeat than T2T</h2>
    <p class="lede2">T2T (CHM13) has three tandem copies of a {kb(period)} kb unit at 16q22. Each unit holds a PDPR-like
    piece, one NPIPB15 copy and one CLEC18 gene. Four independently assembled human haplotypes (HG002, the T2T consortium's
    Q100 v1.1; HG00733, HPRC year 1) carry one or two units. A gene of unit 2 or 3 can therefore sit on one parent's
    chromosome and not on the other's.</p>
    <div class="units">{"".join(tiles)}</div>
    <div class="ulegend"><span><i class="sw u1"></i>CHM13 unit 1 (NPIPB15)</span><span><i class="sw u2"></i>unit 2 (LOC124907808)</span>
      <span><i class="sw u3"></i>unit 3 (LOC124907807)</span><span><i class="sw uf"></i>unique flank</span>
      <span>boxes: PDPR-like (−) · NPIP (+) · CLEC18 (−) per unit</span></div>

    <div class="win"><h3>Where each CHM13 unit aligns in each haplotype</h3>
      <p class="wsub">Each CHM13 segment aligned on its own (minimap2 -x asm10; every hit covering half the segment, faint = secondary).
      A haplotype NPIP copy is coloured by the CHM13 copy it matches best (identity, and the runner-up's); its PDPR-like and
      CLEC18 genes are grey. Ribbons that meet on one unit are CHM13 units that copy the same sequence: the haplotype has one
      copy where CHM13 has two or three.</p>
      {maps}
    </div>

    <div class="win"><h3>Dot plots: every exact 32-bp match between CHM13 (x) and each sequence (y)</h3>
      <p class="wsub">Same bp scale on both axes. A diagonal is the same sequence in the same order; a parallel diagonal is a repeat.
      Coloured bars under each plot mark CHM13's units. <span class="ulegend" style="display:inline-flex;margin:0">
      <span><i class="sw df"></i>same strand</span><span><i class="sw dr"></i>reverse strand</span></span></p>
      <div class="dps">{"".join(dps)}</div>
    </div>

    <div class="win"><h3>What the reads do on T2T (A119b testis, minimap2 -N 50 -p 0.1)</h3>
      <p class="wsub">A read comes from a copy the person carries. Aligned to T2T's three near-identical units, it scores almost
      the same on each, so most reads on LOC124907808 and LOC124907807 are ties with NPIPB15: a person with one or two units
      spreads their reads over three. Margins are AS to the best alignment on another NPIP copy, the 10-AS rule O2 uses
      (the 0.98 ratio in the table above counts all of these as tied).</p>
      <div class="ulegend"><span><i class="sw rxs"></i>exact tie: the same AS on another NPIP copy</span>
        <span><i class="sw rns"></i>near tie: within 1-9 AS</span><span><i class="sw rcs"></i>placed: every other copy ≥ 10 AS lower</span></div>
      <div class="wscroll">{reads_panel(rd)}</div>
      <p class="wsub">HG002's RNA barely expresses this cluster (6 reads on it, none on the two LOC copies), so its reads cannot be traced.
      A119b has no assembly; its 40 clearly placed reads suggest it carries at least one of the LOC copies, which only its DNA can confirm.</p>
    </div>

    <div class="card checks"><h3>Why this is not an assembly or alignment error</h3><ul>
      <li><b>Two people, two assembly pipelines, the same pattern.</b> HG002 (Q100 v1.1, curated diploid T2T assembly) and HG00733
      (HPRC hifiasm): in both, one haplotype carries one unit and the other two.</li>
      <li><b>Whole units only.</b> Each haplotype is shorter than CHM13 by 0.98-2.00 units of {kb(period)} kb (tiles above). A
      misassembly or a misplaced read pile would not come in exact multiples of the repeat.</li>
      <li><b>Anchored on both sides.</b> The unique flanks of every haplotype align at MAPQ 60 and 99.6-99.97% identity, so the
      region is in its true place.</li>
      <li><b>Gene order kept, unit count changed.</b> Every haplotype unit keeps PDPR-like, NPIP and CLEC18 in CHM13's order and
      orientation; only the number of units differs, which is what an unequal crossover between tandem units produces.</li>
      <li><b>Which CHM13 copy a haplotype copy "is" is close to a tie.</b> Each haplotype NPIP copy matches all three CHM13
      copies at 99.3-99.9%, its best by 0.04-0.43 points. The labels name the best; the same closeness is why reads tie.</li>
      <li><b>The direct test left to run:</b> DNA reads (HG002's public HiFi or ultra-long ONT) spanning a unit junction would show
      each haplotype's unit count read by read.</li>
    </ul></div>
    <p class="sub">Data: <code>bench/o3_graph/npip_cnv_fig.py</code> (alignments, gene hits, 32-mer dot plots; region CHM13
    chr16:80,100,000-80,550,000), <code>npip_cnv_render.py</code>; HG002 <code>hg002v1.1.fa</code> chr16_MATERNAL /
    chr16_PATERNAL; HG00733 HPRC year-1 maternal / paternal (paternal shown reverse-complemented). Built 20 Oct 2026.</p>
  </section>
"""
    open(out_p, "w").write(html)
    print(f"[npip_cnv] section -> {out_p} ({len(html):,} bytes)", flush=True)
    return 0


if __name__ == "__main__":
    sys.exit(main(*sys.argv[1:4]))
