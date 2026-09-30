#!/usr/bin/env python3
"""Builds the Soto-vs-ours meeting page (claude.ai artifact J12TB7ebNq3uxZ7aErd7q7) from:
  template.html                    the page (tabs, styles, viewer script) with placeholders
  method_tab.py / misses_panel.py  the "How Soto builds families" tab and the RNA miss-reason cards
  july/*.html                      the two July pages, embedded as tabs in shadow roots
  families.json                    from bench/soto_m2/soto_m2_families.py --out-json

    python3 bench/soto_m2/page/build_page.py --data families.json --out soto_vs_ours.html

Two fixes are applied to the July member page's copy only: its "RNA missed" filter kept the RNA-found rows (313) instead
of the missed ones (49), and its BED button attempted a file download the viewer blocks (it already shows and copies the
text, so the download is dropped)."""
import argparse
import os
import re
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
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
    css += "\n.sbody .wrap{padding-left:0!important;padding-right:0!important;padding-top:8px!important}"
    js = js.replace("if(filter==='rnamiss' && d.rna!=='yes') return false;", "if(filter==='rnamiss' && d.rna!=='no') return false;")
    js = re.sub(r"\n\s*try\{ const a=document\.createElement\('a'\);.*?catch\(e\)\{\}", "", js, count=1, flags=re.S)
    markup = markup.replace("⬇ Get BED — the included loci", "Show BED — the included loci")
    for bad in ("</template>", "</script>"):
        assert bad not in markup + css + (js if bad == "</script>" else ""), bad
    return f"<style>{css}</style><div class=\"sbody\">{markup}</div>", js


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--data", required=True, help="families.json from soto_m2_families.py")
    ap.add_argument("--out", required=True)
    a = ap.parse_args(argv)
    t = open(os.path.join(HERE, "template.html")).read()
    d = open(a.data).read()
    assert "</" not in d
    for key, name in (("DET", "detection_2026-07-28.html"), ("MEM", "members_2026-07-27.html")):
        html, js = split_old(os.path.join(HERE, "july", name))
        t = t.replace(f"__{key}_HTML__", html, 1).replace(f"__{key}_JS__", js, 1)
    t = t.replace("__HOW_PANE__", method_tab.pane() + "\n", 1)
    t = t.replace("__MISSES_DET__", misses_panel.misses_html("d"), 1).replace("__MISSES_MEM__", misses_panel.misses_html("m"), 1)
    t = t.replace("__DATA__", d, 1)
    assert not re.search(r"__[A-Z_]+__", t), re.findall(r"__[A-Z_]+__", t)[:3]
    open(a.out, "w").write(t)
    print(f"built {a.out}: {len(t):,} bytes")


if __name__ == "__main__":
    main()
