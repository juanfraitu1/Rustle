#!/usr/bin/env python3
"""Add the "DNA-mode seeds" tab to the LRPAP1 artifact page (docs/archive/2026-10/merged/PREREG_lrpap1_seed_reach_2026-10-11.md).

    add_seed_tab.py <live page.html> <seed_reach.json> <tab dir> <out.html>

<live page.html> = the artifact as read back (publish skeleton included or not); <tab dir> holds seed_tab.{html,css,js}. The existing content
becomes the "Copy browser" tab unchanged; the new tab reads SEEDS (a subset of seed_reach.json) and the page's own DATA and helpers."""
import json
import re
import sys


def strip_skeleton(page):
    i = page.find("<body>")
    if i >= 0:
        page = page[i + len("<body>"):]
    return re.sub(r"\s*</body>\s*</html>\s*$", "\n", page).lstrip()


def inject(page, css, tab_html, js, seeds):
    if 'id="tab-seeds"' in page:
        raise SystemExit("the page already has the seeds tab")
    page = strip_skeleton(page)
    lede = re.search(r'<p class="lede">.*?</p>', page, re.S)
    assert lede, "no lede"
    bar = ('\n  <div class="tabs" role="tablist" aria-label="views">'
           '<button id="tb-browser" role="tab" aria-controls="tab-browser" aria-selected="true">Copy browser</button>'
           '<button id="tb-seeds" role="tab" aria-controls="tab-seeds" aria-selected="false" tabindex="-1">DNA-mode seeds</button></div>\n'
           '  <div id="tab-browser" role="tabpanel" aria-labelledby="tb-browser">')
    page = page[:lede.end()] + bar + page[lede.end():]
    f = page.index("</footer>") + len("</footer>")
    page = page[:f] + "\n  </div>\n" + tab_html + page[f:]
    w = page.index('<div class="wrap">')
    page = page[:w] + "<style>\n" + css + "</style>\n" + page[w:]
    s = page.rindex("</script>") + len("</script>")
    page = page[:s] + "\n<script>\nconst SEEDS=" + json.dumps(seeds, separators=(",", ":")) + ";\n" + js + "</script>\n" + page[s:]
    return page


def seeds_subset(r):
    keep = ("ids", "finds", "finds_transcript", "finds_body", "extra_loci", "identity", "min_seeds_one_pass", "min_seeds_iterate")
    out = {k: r[k] for k in keep}
    out["sweep"] = [{k: x[k] for k in ("threshold", "min_seeds_one_pass", "min_seeds_iterate", "components")} for x in r["sweep"]]
    out["curves"] = {k: {"one_pass": r["curves"][k]["one_pass"]} for k in ("0.800", "0.980", "0.990", "0.995")}
    return out


def main(live, res, tabdir, out):
    r = json.load(open(res))
    page = inject(open(live).read(), open(f"{tabdir}/seed_tab.css").read(), open(f"{tabdir}/seed_tab.html").read(),
                  open(f"{tabdir}/seed_tab.js").read(), seeds_subset(r))
    open(out, "w").write(page)
    print(f"wrote {out} ({len(page):,} bytes)")


if __name__ == "__main__":
    main(*sys.argv[1:5])
