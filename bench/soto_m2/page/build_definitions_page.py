#!/usr/bin/env python3
"""Builds the companion page "Gene Family Definitions" (the "Two definitions" and "How we build families" tabs, moved out of the
Soto-vs-ours meeting page on 2026-09-30). Styles come from the meeting page's template, so the two pages always look alike.

    python3 build_definitions_page.py --template template.html --out definitions.html
"""
import argparse
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import definitions_tab  # noqa: E402
import mcl_tab  # noqa: E402

JS = r"""<script>
(function () {
  const tip = document.getElementById("tip");
  const esc = s => String(s).replace(/[&<>"]/g, c => ({ "&": "&amp;", "<": "&lt;", ">": "&gt;", '"': "&quot;" }[c]));
  document.addEventListener("mousemove", e => {
    const dt = e.target.closest && e.target.closest("[data-tip]");
    if (!dt) { tip.hidden = true; return; }
    tip.innerHTML = esc(dt.dataset.tip); tip.hidden = false;
    const w = tip.offsetWidth, h = tip.offsetHeight;
    let x = e.clientX + 14, y = e.clientY + 14;
    if (x + w > innerWidth - 8) x = e.clientX - w - 14;
    if (y + h > innerHeight - 8) y = e.clientY - h - 14;
    tip.style.left = x + "px"; tip.style.top = y + "px";
  });
  const tabs = [...document.querySelectorAll('.tabs [role="tab"]')];
  function select(p, focus) {
    tabs.forEach(b => { const on = b.dataset.p === p; b.setAttribute("aria-selected", on); b.tabIndex = on ? 0 : -1;
      document.getElementById("pane-" + b.dataset.p).hidden = !on; if (on && focus) b.focus(); });
    try { localStorage.setItem("defs-tab", p); } catch (e) {}
  }
  tabs.forEach((b, k) => {
    b.addEventListener("click", () => select(b.dataset.p));
    b.addEventListener("keydown", e => { if (e.key === "ArrowRight" || e.key === "ArrowLeft") {
      const n = tabs[(k + (e.key === "ArrowRight" ? 1 : tabs.length - 1)) % tabs.length]; select(n.dataset.p, true); } });
  });
  const h = (location.hash || "").slice(1);
  let start = h === "mcl" ? "mcl" : h === "definitions" ? "def" : null;
  if (!start) { try { start = localStorage.getItem("defs-tab"); } catch (e) {} }
  select(start === "mcl" ? "mcl" : "def");
})();
</script>"""


def build(template_path):
    t = open(template_path).read()
    head = t[:t.index("</style>") + len("</style>")].replace("<title>Soto vs Our Families</title>", "<title>Gene Family Definitions</title>", 1)
    defs = definitions_tab.pane()
    # the chips that opened a family in the meeting page's viewer have no viewer here
    a = defs.index('<div class="bar"><button class="chip" data-real=')
    b = defs.index("</div>", a) + len("</div>")
    defs = defs[:a] + ('<p class="note">In the "Soto vs Our Families" page, open NPIPB3 as ID_151 and the NPIP A-clade as ID_154 in '
                       'the family viewer.</p>') + defs[b:]
    body = f"""
<div class="wrap">
<header>
  <div class="eyebrow">Companion to "Soto vs Our Families" · human CHM13 · 30 Sep 2026</div>
  <h1>Gene Family Definitions</h1>
</header>
<nav class="tabs" role="tablist" aria-label="Views">
  <button role="tab" id="tab-def" aria-controls="pane-def" aria-selected="true" data-p="def">Two definitions<span>what "gene family" means</span></button>
  <button role="tab" id="tab-mcl" aria-controls="pane-mcl" aria-selected="false" data-p="mcl">How we build families<span>MCL, explained plainly</span></button>
</nav>
{defs}
{mcl_tab.pane()}
</div>
<div id="tip" hidden></div>
"""
    return head + body + JS + "\n"


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--template", default=os.path.join(HERE, "template.html"))
    ap.add_argument("--out", required=True)
    a = ap.parse_args(argv)
    page = build(a.template)
    open(a.out, "w").write(page)
    print(f"built {a.out}: {len(page):,} bytes")


if __name__ == "__main__":
    main()
