"""The "Gene ends" tab: per NPIP / TBC1D3 copy (human, A119b), the RefSeq gene, Soto's gene model (CAT v4) and our best model, each
end judged by exons against the RefSeq gene (data: bench/soto_m2/soto_m2_gene_ends.py)."""
import html
import json
from collections import Counter

LONG, SHORT = "var(--cross)", "var(--s1)"


def fmt(n):
    return f"{n:,}"


def call_txt(k, bp):
    if k > 0:
        return f"+{k} exon{'s' if k > 1 else ''}", LONG
    if k < 0:
        return f"−{-k} exon{'s' if k < -1 else ''}", SHORT
    if bp:
        return f"same ({'+' if bp > 0 else '−'}{fmt(abs(bp))} bp)", "var(--muted)"
    return "same", "var(--muted)"


def lane_svg(exons, y, X, r_lo, r_hi, tip):
    """One gene model at height y: exons and intron line; parts beyond the RefSeq gene in red; the part of the RefSeq gene it
    lacks as a pale blue band."""
    s = [f'<g data-tip="{html.escape(tip, quote=True)}">']
    lo, hi = exons[0][0], exons[-1][1]
    s.append(f'<rect x="{X(lo) - 2:.1f}" y="{y - 7}" width="{X(hi) - X(lo) + 4:.1f}" height="14" fill="transparent"/>')
    if lo > r_lo:
        s.append(f'<rect x="{X(r_lo):.1f}" y="{y - 5}" width="{X(min(lo, r_hi)) - X(r_lo):.1f}" height="10" fill="{SHORT}" fill-opacity=".16"/>')
    if hi < r_hi:
        s.append(f'<rect x="{X(max(hi, r_lo)):.1f}" y="{y - 5}" width="{X(r_hi) - X(max(hi, r_lo)):.1f}" height="10" fill="{SHORT}" fill-opacity=".16"/>')

    def seg(a, b, inside_col, kind):
        out = []
        cuts = sorted({a, b} | ({r_lo} if a < r_lo < b else set()) | ({r_hi} if a < r_hi < b else set()))
        for p, q in zip(cuts, cuts[1:]):
            col = inside_col if r_lo <= p and q <= r_hi else LONG
            if kind == "line":
                out.append(f'<line x1="{X(p):.1f}" x2="{X(q):.1f}" y1="{y}" y2="{y}" stroke="{col}" stroke-width="1.2"/>')
            else:
                out.append(f'<rect x="{X(p):.1f}" y="{y - 5}" width="{max(X(q) - X(p), 1.2):.1f}" height="10" fill="{col}"/>')
        return out
    s += seg(lo, hi, "var(--igv)", "line")
    for a, b in exons:
        s += seg(a, b, "var(--igv)", "box")
    s.append("</g>")
    return "".join(s)


def copy_svg(r):
    ref = r["ref"]
    r_lo, r_hi = ref[0][0], ref[-1][1]
    lanes = [("RefSeq", ref, None, f"RefSeq {r['name']} ({r['strand']}): {r['chrom']}:{fmt(r_lo + 1)}-{fmt(r_hi)}")]
    so, ou = r.get("soto"), r.get("ours")
    if so:
        lanes.append(("Soto (CAT)", so["exons"], so["ends"], f"Soto's gene model: CAT v4 {so['name']}, in Soto family "
                      f"{', '.join(so['fams']) or 'none'}"))
    else:
        lanes.append(("Soto (CAT)", None, None, ""))
    if ou:
        fz = f"; fused with {', '.join(ou['fused_with'])}" if ou["fused_with"] else ""
        lanes.append(("ours", ou["exons"], ou["ends"], f"our best model {ou['tid']} (class {ou['cls']}, {ou['support']:g} reads; "
                      f"{ou['n_models']} models at this copy{fz})"))
    else:
        lanes.append(("ours", None, None, ""))
    lo = min([r_lo] + [e[0][0] for _, e, _, _ in lanes if e])
    hi = max([r_hi] + [e[-1][1] for _, e, _, _ in lanes if e])
    pad = (hi - lo) * 0.02 + 50
    lo, hi = lo - pad, hi + pad
    G0, G1 = 120, 620
    X = lambda p: G0 + (p - lo) / (hi - lo) * (G1 - G0)
    H = 22 + 18 * len(lanes) + 4
    s = [f'<svg viewBox="0 0 900 {H}" width="100%" role="img" aria-label="{html.escape(r["name"])}: RefSeq gene, Soto gene model and our model">']
    arrow = "→" if r["strand"] == "+" else "←"
    s.append(f'<text x="0" y="13" font-size="12.5" font-weight="650" fill="var(--fg)">{html.escape(r["name"])} {arrow}</text>')
    s.append(f'<text x="{max(G0, (len(r["name"]) + 2) * 7.8 + 12):.0f}" y="13" font-size="10.5" fill="var(--muted)">{r["chrom"]}:{fmt(r_lo + 1)}-{fmt(r_hi)} · '
             f'{fmt(r_hi - r_lo)} bp</text>')
    s.append(f'<line x1="{X(r_lo):.1f}" x2="{X(r_lo):.1f}" y1="20" y2="{H - 2}" stroke="var(--faint)" stroke-dasharray="3 3"/>')
    s.append(f'<line x1="{X(r_hi):.1f}" x2="{X(r_hi):.1f}" y1="20" y2="{H - 2}" stroke="var(--faint)" stroke-dasharray="3 3"/>')
    for i, (lab, ex, ends, tip) in enumerate(lanes):
        y = 30 + 18 * i
        s.append(f'<text x="{G0 - 8}" y="{y + 4}" font-size="10.5" fill="var(--muted)" text-anchor="end">{lab}</text>')
        if not ex:
            s.append(f'<text x="{G0}" y="{y + 4}" font-size="10.5" fill="var(--faint)">no {"CAT gene on this strand" if lab.startswith("Soto") else "model"}</text>')
            continue
        s.append(lane_svg(ex, y, X, r_lo, r_hi, tip))
        if ends:
            (k5, b5), (k3, b3) = ends
            t5, c5 = call_txt(k5, b5)
            t3, c3 = call_txt(k3, b3)
            s.append(f'<text x="{G1 + 14}" y="{y + 4}" font-size="10.5" xml:space="preserve"><tspan fill="var(--muted)">5′\u00a0</tspan>'
                     f'<tspan fill="{c5}" font-weight="{650 if k5 else 400}">{t5}</tspan><tspan fill="var(--muted)">\u00a0·\u00a03′\u00a0</tspan>'
                     f'<tspan fill="{c3}" font-weight="{650 if k3 else 400}">{t3}</tspan></text>')
    s.append("</svg>")
    return "".join(s)


def pane(data):
    rows = data
    tiles = []
    for lane, lab in (("soto", "Soto's gene models (CAT v4)"), ("ours", "Our best model (A119b reads)")):
        n = sum(1 for r in rows if lane in r)
        c = Counter()
        for r in rows:
            if lane in r:
                for side, (k, _b) in zip(("5", "3"), r[lane]["ends"]):
                    c[(side, (k > 0) - (k < 0))] += 1
        cells = "".join(
            f'<tr><th scope="row">{side}′ end</th><td style="color:{LONG}">{c[(side, 1)]}</td><td>{c[(side, 0)]}</td>'
            f'<td style="color:{SHORT}">{c[(side, -1)]}</td></tr>' for side in ("5", "3"))
        tiles.append(f'<div class="tile"><div class="lab">{lab} · {n} copies</div><table class="endtbl"><thead><tr><th></th>'
                     f'<th style="color:{LONG}">longer</th><th>same exon</th><th style="color:{SHORT}">shorter</th></tr></thead>'
                     f'<tbody>{cells}</tbody></table></div>')
    pk = [r for r in rows if r.get("soto") and r["soto"]["ends"][0][0] > 0 and "ID_149" in r["soto"]["fams"]]
    pk_txt = ", ".join(f'{html.escape(r["name"])} (+{r["soto"]["ends"][0][0]} exons, CAT {html.escape(r["soto"]["name"])})' for r in pk)
    body = []
    for fam in ("NPIP", "TBC1D3"):
        fr = [r for r in rows if r["family"] == fam]
        body.append(f'<h3 class="sub-h">{fam} · {len(fr)} copies</h3>')
        body += [f'<div class="endrow">{copy_svg(r)}</div>' for r in fr]
    return f'''<div role="tabpanel" id="pane-ends" aria-labelledby="tab-ends" class="pane" hidden>
<section class="how-intro">
  <div class="eyebrow">Gene ends · human NPIP and TBC1D3 copies, A119b Iso-Seq reads</div>
  <h2>Where each prediction starts and stops, against the RefSeq gene</h2>
  <p class="note">One row per annotated copy. <b>RefSeq</b> is the reference gene (all its transcripts collapsed); <b>Soto (CAT)</b> is the
  gene model Soto's families are made of (CAT v4, matched by exon overlap, so its name can differ); <b>ours</b> is our best model at the
  copy (exact chain first, then the best match class, then read support). Each end is judged by exons, with no bp threshold: <span
  style="color:{LONG};font-weight:650">+k exons</span> = k exons beyond the RefSeq gene's first or last exon, <span
  style="color:{SHORT};font-weight:650">−k exons</span> = k of its exons missing; "same" means the end falls in the same exon, with how far it moves.</p>
</section>
<section>
  <div class="cn-facts" style="grid-template-columns:repeat(auto-fit,minmax(280px,1fr))">{"".join(tiles)}</div>
  <p class="note">Soto's models that run longer at the 5′ end and sit in Soto's PKD1P family (ID_149): {pk_txt or "none"}. These CAT models
  carry the upstream PKD1 pieces, so the NPIP copy they cover is grouped with PKD1P, not with NPIP — the fusion architecture of the
  duplication block, seen in the annotation.</p>
  <div class="legend"><span><svg width="18" height="10" aria-hidden="true"><rect y="1" width="18" height="8" fill="var(--igv)"/></svg> exon inside the RefSeq gene</span>
  <span><svg width="18" height="10" aria-hidden="true"><rect y="1" width="18" height="8" fill="{LONG}"/></svg> beyond the RefSeq gene</span>
  <span><svg width="18" height="10" aria-hidden="true"><rect y="1" width="18" height="8" fill="{SHORT}" fill-opacity=".16"/></svg> part of the RefSeq gene the model lacks</span>
  <span>dashed lines = RefSeq gene start and end · coordinates CHM13 v2.0 (CAT moved from v1.0: chr16 −5 bp, chr17 −291 bp)</span></div>
</section>
<section class="ends">{"".join(body)}</section>
<p class="note">Five TBC1D3 pseudogenes (TBC1D3P1, P3, P4, P7, LOC100420311) are annotated without exons and are not drawn. Models from the
2026-09-29 copy-recovery comparison (our shipped defaults on A119b); Soto's genes from Table S1C's CAT v4 models.</p>
</div>'''


if __name__ == "__main__":
    import sys
    print(pane(json.load(open(sys.argv[1]))))
