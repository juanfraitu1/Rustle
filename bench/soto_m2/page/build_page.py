#!/usr/bin/env python3
"""Builds the Soto-vs-ours meeting page (claude.ai artifact J12TB7ebNq3uxZ7aErd7q7) from:
  template.html                    the page (tabs, styles, viewer script) with placeholders
  method_tab.py / misses_panel.py  the "How Soto builds families" tab and the RNA miss-reason cards
  definitions_tab.py, mcl_tab.py   moved to the companion page (build_definitions_page.py, artifact 9yXmXZd2kRrQffQVzEJGti)
  july/*.html                      the two July pages, embedded as tabs in shadow roots
  families.json                    from bench/soto_m2/soto_m2_families.py --out-json
  sd_regions.json                  from bench/soto_m2/soto_m2_sd_regions.py (SD98 regions, their genes and duplicons)
  gene_ends.json                   from bench/soto_m2/soto_m2_gene_ends.py (the "Gene ends" tab; frozen copy
                                   docs/SOTO_GENE_ENDS_2026-09-30.json)
  cat_labels.json                  from bench/soto_m2/soto_m2_cat_labels.py: CAT/Liftoff v2.0 labels for the Detection tab
                                   (--cat-labels; without it the tab keeps the July RefSeq labels)

    python3 bench/soto_m2/page/build_page.py --data families.json --sd sd_regions.json --out soto_vs_ours.html

Two fixes are applied to the July member page's copy only: its "RNA missed" filter kept the RNA-found rows (313) instead
of the missed ones (49), and its BED button attempted a file download the viewer blocks (it already shows and copies the
text, so the download is dropped)."""
import argparse
import json
import os
import re
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import definitions_tab  # noqa: E402
import gene_ends_tab  # noqa: E402
import july_cat  # noqa: E402
import july_dna_only  # noqa: E402
import mcl_tab  # noqa: E402
import method_tab  # noqa: E402
import misses_panel  # noqa: E402


def split_old(path):
    t = open(path).read()
    t = t[t.index("<title>"):]
    t = re.sub(r"<title>.*?</title>", "", t, count=1, flags=re.S)
    t = re.sub(r"</body>\s*</html>\s*$", "", t.strip())
    m = list(re.finditer(r"<script>(.*?)</script>", t, re.S))
    assert len(m) == 1, len(m)
    js = m[0].group(1)
    rest = t[:m[0].start()] + t[m[0].end():]
    css = "\n".join(re.findall(r"<style>(.*?)</style>", rest, re.S))
    markup = re.sub(r"<style>.*?</style>", "", rest, flags=re.S).strip()
    css = re.sub(r':root\[data-theme="?dark"?\]', ':host([data-theme="dark"])', css)
    css = re.sub(r':root\[data-theme="?light"?\]', ':host([data-theme="light"])', css)
    css = css.replace(":root", ":host")
    css = re.sub(r"(^|[}\s])body\s*\{", r"\1.sbody{", css)
    css += "\n.sbody .wrap{padding-left:0!important;padding-right:0!important;padding-top:8px!important;max-width:none!important}"
    js = js.replace("if(filter==='rnamiss' && d.rna!=='yes') return false;", "if(filter==='rnamiss' && d.rna!=='no') return false;")
    js = re.sub(r"\n\s*try\{ const a=document\.createElement\('a'\);.*?catch\(e\)\{\}", "", js, count=1, flags=re.S)
    markup = markup.replace("⬇ Get BED — the included loci", "Show BED — the included loci")
    for bad in ("</template>", "</script>"):
        assert bad not in markup + css + (js if bad == "</script>" else ""), bad
    return f"<style>{css}</style><div class=\"sbody\">{markup}</div>", js


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--data", required=True, help="families.json from soto_m2_families.py")
    ap.add_argument("--sd", required=True, help="sd_regions.json from soto_m2_sd_regions.py")
    ap.add_argument("--ends", help="gene_ends.json from soto_m2_gene_ends.py (only if the template has the Gene ends tab)")
    ap.add_argument("--cat-labels", help="cat_labels.json from soto_m2_cat_labels.py (Detection tab on CAT/Liftoff v2.0)")
    ap.add_argument("--out", required=True)
    a = ap.parse_args(argv)
    t = open(os.path.join(HERE, "template.html")).read()
    d = open(a.data).read()
    assert "</" not in d
    for key, name in (("DET", "detection_2026-07-28.html"), ("MEM", "members_2026-07-27.html")):
        html, js = split_old(os.path.join(HERE, "july", name))
        html, js = (july_dna_only.det if key == "DET" else july_dna_only.mem)(html, js)   # DNA to DNA only (2026-09-30)
        if key == "DET" and a.cat_labels:
            html, js = july_cat.det(html, js, july_cat.load(a.cat_labels))           # CAT/Liftoff v2.0 labels (2026-10-01)
        t = t.replace(f"__{key}_HTML__", html, 1).replace(f"__{key}_JS__", js, 1)
    t = t.replace("__HOW_PANE__", method_tab.pane() + "\n", 1)
    t = t.replace("__DEF_PANE__", definitions_tab.pane(), 1)
    t = t.replace("__MCL_PANE__", mcl_tab.pane(), 1)
    if "__ENDS_PANE__" in t:
        t = t.replace("__ENDS_PANE__", gene_ends_tab.pane(json.load(open(a.ends))), 1)
    t = t.replace("__MISSES_DET__", misses_panel.misses_html("d"), 1).replace("__MISSES_MEM__", misses_panel.misses_html("m"), 1)
    t = t.replace("__DATA__", d, 1)
    sd = open(a.sd).read()
    assert "</" not in sd
    t = t.replace("__SD__", sd, 1)
    assert not re.search(r"__[A-Z_]+__", t), re.findall(r"__[A-Z_]+__", t)[:3]
    open(a.out, "w").write(t)
    print(f"built {a.out}: {len(t):,} bytes")


if __name__ == "__main__":
    main()
