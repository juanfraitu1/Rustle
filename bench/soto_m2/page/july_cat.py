"""CAT/Liftoff v2.0 labels in the Detection tab (2026-10-01; user decision: CAT is the default human annotation). Applied after
july_dna_only.det. Data from bench/soto_m2/soto_m2_cat_labels.py: every member's class, annotation and suggested exclusion are
replaced by the CAT labels (same classifier as the rebuilt RefSeq labels; the July RefSeq label stays in each row's tooltip and,
where the flag changed, beside the row). Each DNA-mode extra copy gets the CAT gene(s) it overlaps and whether a family-sized CAT gene
explains its size outlier. Each edit must match exactly once."""
import html as H
import json


def _rep(s, a, b):
    assert s.count(a) == 1, (s.count(a), a[:80])
    return s.replace(a, b)


def _swap_json(js, head, fn):
    """Replace the JSON object that follows `head` in js with fn(object)."""
    i = js.index(head) + len(head)
    i = js.index("{", i)
    obj, n = json.JSONDecoder().raw_decode(js[i:])
    out = json.dumps(fn(obj), separators=(",", ":"), ensure_ascii=False)
    assert "</" not in out
    return js[:i] + out + js[i + n:]


def panel(s):
    c = s["counts"]
    ex = s["excluded"]
    e = s["extras"]
    g = lambda k: e.get(k, 0)
    big, small = g("big|block") + g("big|nofit") + g("big|nogene"), g("small|piece") + g("small|nofit") + g("small|nogene")
    rows = [("unannotated (no gene overlaps)", "unannotated"), ("fragment (&lt; 800 bp)", "fragment"),
            ("piece of another gene", "piece-of-other")]
    tr = "".join(f'<tr><th scope="row">{lab}</th><td>{c["registered"].get(k, 0)}</td><td>{c["cat"].get(k, 0)}</td></tr>'
                 for lab, k in rows)
    return f'''<style>
  .catcmp{{display:grid;grid-template-columns:repeat(auto-fit,minmax(300px,1fr));gap:14px;margin:4px 0 16px}}
  .catcmp .box{{border:1px solid var(--line);border-radius:12px;background:var(--panel);padding:12px 14px;min-width:0}}
  .catcmp h3{{font-size:13px;margin:0 0 8px;color:var(--ink)}}
  .catcmp table{{border-collapse:collapse;font-size:12.5px;width:100%;font-variant-numeric:tabular-nums}}
  .catcmp th,.catcmp td{{padding:4px 8px;border-bottom:1px solid var(--line);text-align:right}}
  .catcmp th[scope=row],.catcmp thead th:first-child{{text-align:left;font-weight:500;color:var(--muted)}}
  .catcmp tfoot td,.catcmp tfoot th{{font-weight:650;color:var(--ink);border-bottom:0}}
  .catcmp p{{font-size:12.5px;color:var(--muted);margin:8px 0 0;line-height:1.5}}
  .catcmp b.k{{color:var(--ink)}}
  .was{{font-family:var(--mono);font-size:10px;color:var(--faint);margin-left:6px;white-space:nowrap}}
  .det.szfit{{background:var(--excl-soft);color:var(--excl);white-space:nowrap}}
</style>
<div class="catcmp">
  <div class="box"><h3>Soto's 362 members: suggested "segdup piece" exclusions</h3>
    <table><thead><tr><th></th><th>RefSeq (July)</th><th>CAT v2.0</th></tr></thead><tbody>{tr}</tbody>
    <tfoot><tr><th scope="row">suggested exclusions</th><td>{ex["registered"]}</td><td>{ex["cat"]}</td></tr></tfoot></table>
    <p>Same rules for both. The July script was not kept; the rebuilt rules give the July label for {s["rebuilt_agreement"]} of
    the 362 members and nearly its counts. Every Soto member is exactly a CAT gene, so under CAT no member is unannotated or a
    piece of another gene: the <b class="k">{c["registered"].get("piece-of-other", 0) + c["registered"].get("unannotated", 0)}</b>
    RefSeq flags of those two kinds came from RefSeq lacking or renaming the gene. What is left is size: <b class="k">{c["cat"].get("fragment", 0)}</b>
    members shorter than 800 bp.</p></div>
  <div class="box"><h3>{s["n_extras"]} extra copies (DNA mode): size against the family median</h3>
    <table><thead><tr><th></th><th>&gt; 2× median</th><th>&lt; 0.5× median</th></tr></thead><tbody>
    <tr><th scope="row">contains, or lies inside, a family-sized CAT gene</th><td>{g("big|block")}</td><td>{g("small|piece")}</td></tr>
    <tr><th scope="row">CAT genes, none of family size</th><td>{g("big|nofit")}</td><td>{g("small|nofit")}</td></tr>
    <tr><th scope="row">no CAT gene at all</th><td>{g("big|nogene")}</td><td>{g("small|nogene")}</td></tr></tbody>
    <tfoot><tr><th scope="row">size outliers</th><td>{big}</td><td>{small}</td></tr></tfoot></table>
    <p>The size flag uses no annotation, so CAT does not change it. A family-sized CAT gene is 0.5–2× the family median and
    overlaps the copy by at least half of the shorter of the two. <b class="k">{g("big|block")}</b> oversized copies contain such a
    gene and <b class="k">{g("small|piece")}</b> undersized copies lie inside one: measured by that gene instead of by the aligned
    block, these <b class="k">{g("big|block") + g("small|piece")}</b> copies are family-sized. The other
    <b class="k">{big + small - g("big|block") - g("small|piece")}</b> outliers stay outliers. Gene content is not homology: the gene
    is named in each row, not proven to be the paralog, and a large block can contain a small family-sized gene.</p></div>
</div>
'''


def det(markup, js, labels):
    L = labels
    mem, ext = L["members"], L["extras"]

    def fix_d(D):
        for f in D["fams"]:
            for m in f["members"]:
                v = mem[f"{f['fam']}:{m['start']}"]
                m.update(cls=v["cls"], ann=v["ann"], excl=v["excl"], nested=v["nested"], cls_rs=v["cls_rs"], ann_rs=v["ann_rs"])
        cc = {}
        for f in D["fams"]:
            for m in f["members"]:
                cc[m["cls"]] = cc.get(m["cls"], 0) + 1
        D["cls_counts"] = cc
        return D

    def fix_x(X):
        for fam, xs in X.items():
            assert len(xs) == len(ext[fam]), fam
            for x, v in zip(xs, ext[fam]):
                x.update(cs=v["status"], cg=[dict(n=g["name"], b=g["bp"], t=g["bt"]) for g in v["genes"]], ng=v["n_genes"],
                         cf=dict(n=v["fit"]["name"], b=v["fit"]["bp"], t=v["fit"]["bt"]) if v["fit"] else None, cr=v["fit_ratio"])
        return X

    js = _swap_json(js, "const D = ", fix_d)
    js = _swap_json(js, "const DEXTRA=", fix_x)
    # member rows: CAT label, nesting as information, the July RefSeq label in the tooltip and beside rows whose flag changed
    js = _rep(js, "const annTxt = m.ann? `${CLS_LABEL[m.cls]} · ${m.ann}` : CLS_LABEL[m.cls];",
              "const annTxt = (m.ann? `${CLS_LABEL[m.cls]} · ${m.ann}` : CLS_LABEL[m.cls]) + (m.nested? ` · inside ${m.nested}` : '');\n"
              "        const annTip = `CAT v2.0: ${annTxt} | RefSeq (July): ${CLS_LABEL[m.cls_rs]||m.cls_rs}${m.ann_rs? ' · '+m.ann_rs : ''}`;\n"
              "        const was = (EXCL_CLS.has(m.cls_rs) && !EXCL_CLS.has(m.cls))? `<span class=\"was\">RefSeq: ${CLS_LABEL[m.cls_rs]}</span>` : '';")
    js = _rep(js, '<span class="annb ${flag?\'ex\':\'ok\'}" title="${annTxt}">${annTxt}</span>',
              '<span style="min-width:0;display:flex;align-items:center"><span class="annb ${flag?\'ex\':\'ok\'}" title="${annTip}">'
              '${annTxt}</span>${was}</span>')
    # extra-copy rows: CAT gene content and whether a family-sized CAT gene explains the size
    js = _rep(js, "const X = c.extra ? `unlisted paralog · ${((c.end-c.start)/1000).toFixed(1)} kb` :",
              "const catTxt = c.cs===undefined? '' : c.cs==='nogene'? ' · no CAT gene' : c.cf? ` · CAT ${c.cf.n} ${(c.cf.b/1000).toFixed(1)} kb "
              "(${c.cr}× median)` : ` · CAT ${c.cg[0].n} ${(c.cg[0].b/1000).toFixed(1)} kb${c.ng>1?` +${c.ng-1} more`:''}`;\n"
              "          const X = c.extra ? `unlisted paralog · ${((c.end-c.start)/1000).toFixed(1)} kb${catTxt}` :")
    js = _rep(js, "const status = out? `<span class=\"det szout\" title=\"size ${r.toFixed(2)}× the family's median member — likely "
                  "artifact\">⚠ ${r.toFixed(1)}× size</span>` : `<span class=\"det c\">candidate</span>`;",
              "const fits = c.cs==='block'||c.cs==='piece';\n"
              "          const why = c.cs==='block'? 'a block around a family-sized CAT gene' : c.cs==='piece'? 'inside a family-sized CAT gene' : "
              "c.cs==='nogene'? 'no CAT gene here' : 'CAT genes here, none of family size';\n"
              "          const status = out? `<span class=\"det ${fits?'szfit':'szout'}\" title=\"size ${r.toFixed(2)}× the family's median member; "
              "${why}\">⚠ ${r.toFixed(1)}× · ${fits?(c.cs==='block'?'family-sized gene inside':'inside a family-sized gene'):(c.cs==='nogene'?'no CAT gene':'no gene fits')}"
              "</span>` : `<span class=\"det c\">candidate</span>`;")
    # one more lever: size outliers that no family-sized CAT gene explains
    js = _rep(js, "const candHtml = `<button class=\"clsbtn cand\" data-cls=\"__candout\" aria-pressed=\"${outOn}\">Remove size-outlier "
                  "candidates <span class=\"ct\">${outs.length}</span></button>`",
              "const nofit=curCands().filter(c=>isOutlier(c)&&!(c.cs==='block'||c.cs==='piece')), nofitOn = nofit.length && nofit.every(c=>excl.has(ckey(c)));\n"
              "  const candHtml = `<button class=\"clsbtn cand\" data-cls=\"__candout\" aria-pressed=\"${outOn}\">Remove size-outlier "
              "candidates <span class=\"ct\">${outs.length}</span></button>`\n"
              "    + `<button class=\"clsbtn cand\" data-cls=\"__candnofit\" aria-pressed=\"${nofitOn}\">Remove size outliers with no "
              "family-sized CAT gene <span class=\"ct\">${nofit.length}</span></button>`")
    js = _rep(js, "else if(cls==='__cand'){",
              "else if(cls==='__candnofit'){ const o=curCands().filter(c=>isOutlier(c)&&!(c.cs==='block'||c.cs==='piece')); const on=o.every(c=>excl.has(ckey(c))); "
              "o.forEach(c=> on?excl.delete(ckey(c)):excl.add(ckey(c))); }\n    else if(cls==='__cand'){")
    js = _rep(js, "'Nothing excluded — full Soto truth set + all 37 candidates as potential FPs.'",
              "`Nothing excluded — full Soto truth set + all ${curCands().length} ${mode==='dna'?'extra copies':'candidates'} as potential FPs.`")
    # "With candidates" filter: use the current mode's candidates (DNA mode = the extra copies)
    js = _rep(js, "if(filter==='cand') return (f.cands||[]).length>0;", "if(filter==='cand') return curFamCands(f).length>0;")
    markup = _rep(markup, '<div class="clsbtns" id="clsbtns"></div>', panel(L["summary"]) + '<div class="clsbtns" id="clsbtns"></div>')
    markup = _rep(markup, "Member annotation from T2T-CHM13 v2.0 RefSeq (GCF_009914755): <code>unannotated</code> = no gene overlaps the "
                          "locus; <code>fragment</code> = &lt; 800 bp; <code>piece of another gene</code> = mostly inside a "
                          "differently-named gene.",
                  "Member annotation from T2T-CHM13 v2.0 CAT/Liftoff (2026-10-01; the July labels used RefSeq GCF_009914755 and "
                  "stay in each row's tooltip): <code>unannotated</code> = no gene overlaps the locus; <code>fragment</code> = "
                  "&lt; 800 bp; <code>piece of another gene</code> = no gene of the member's own identity covers 75% of it and "
                  "75% of it lies inside one other gene at least twice its length. <code>inside …</code> = the member is its own "
                  "CAT gene but also lies inside a longer CAT gene (information, not a flag). Extra copies list the CAT genes they "
                  "overlap; <code>✓</code> = a family-sized CAT gene (0.5–2× the family median) explains the size outlier.")
    markup = _rep(markup, '<button class="filterbtn" data-f="cand" aria-pressed="false">With candidates</button>',
                  '<button class="filterbtn" data-f="cand" aria-pressed="false">With extra copies</button>')
    return markup, js


def load(path):
    return json.load(open(path))


if __name__ == "__main__":
    import sys
    print(H.escape(panel(load(sys.argv[1])["summary"]))[:400])
