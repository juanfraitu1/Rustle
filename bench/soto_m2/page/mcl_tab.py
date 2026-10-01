"""The "How MCL builds ours" tab: Markov clustering, step by step, on the same toy graph as the unit test that guards it
in our code (annotation_families.rs, mcl_splits_two_cliques_joined_by_a_single_bridge), computed live in the page with
the same rule set (self-loops 1, column-normalised, expand = M x M, inflate = entries^r then prune < 1e-9 and
renormalise, stop when no entry moves by 1e-7, attractor = each column's heaviest row, ties to the lowest index).
Prints the pane HTML (with its own script)."""

STEPS = [
    ("Draw the copies and their likeness",
     "Each gene copy is a dot. When two copies look alike, we draw a line between them; the more alike, the thicker the line. "
     "Each dot also gets a small loop to itself, so a walker is allowed to stay put."),
    ("Let the walkers wander (\"expand\")",
     "Picture a crowd standing on every dot. At each move, people step to a neighbouring dot, preferring thick lines. Inside a "
     "real family there are many lines, so people keep circulating inside it. Between two families there is only one stray "
     "line, so only a few people cross."),
    ("Favour the busy paths (\"inflate\")",
     "Like a path across a lawn that gets worn in because people use it, busy routes are made busier and quiet routes "
     "quieter. The knob r says how strongly; our pipeline uses 2.8. After this step the few people crossing the stray line "
     "become even fewer."),
    ("Repeat until nothing changes, then read the families",
     "Wander, favour, wander, favour... until the picture stops changing. In the end every dot sends almost all of its people "
     "to one place. Dots whose people end up in the same place are one family."),
]

PARAMS = [
    ("A dot", "one gene copy (a locus): found from the RNA reads, or taken from the annotation in guided mode"),
    ("A line is drawn when", "the two copies line up letter by letter with at least 70% of letters matching over at least 300 "
                             "letters, the match covers at least 30% of the longer copy (or 70% of the shorter one), and it "
                             "includes at least 60% of the smaller copy's exons"),
    ("Line thickness", "how alike the two copies are × how much of the copy the match covers"),
    ("How strongly busy paths win", "r = 2.8; anything below one in a billion is dropped; at most 100 rounds"),
    ("Smallest family reported", "2 copies"),
]

def pane():
    steps = "".join(
        f'<li class="step"><div class="sh"><span class="sn">{k}</span><div><h3>{t}</h3><p>{d}</p></div></div></li>'
        for k, (t, d) in enumerate(STEPS, 1))
    params = "".join(f'<tr><th scope="row">{a}</th><td>{b}</td></tr>' for a, b in PARAMS)
    return f'''<div role="tabpanel" id="pane-mcl" aria-labelledby="tab-mcl" class="pane" hidden>
<section class="how-intro">
  <div class="eyebrow">How our families are built · MCL (Markov clustering), explained without math</div>
  <h2>Finding families by letting walkers wander</h2>
  <p class="note">Below are six gene copies. A1, A2 and A3 are copies of one gene and look very alike (thick lines, 0.95);
  B1, B2 and B3 are copies of another. One stray line (0.72, red) joins A3 and B1: the kind of accidental link a shared
  repeat or a fused transcript creates. If we simply called "anything connected" one family, that single line would glue
  the two families together. MCL is how we avoid that. Press the buttons in order to watch it happen.
  (This is the same example our own code is tested on.)</p>
</section>
<section>
  <div class="bar mcl-ctl">
    <button class="chip" id="mclStart">1 · Start over</button>
    <button class="chip" id="mclExp">2 · Let them wander</button>
    <button class="chip" id="mclInf">3 · Favour busy paths</button>
    <button class="chip" id="mclEnd">Repeat to the end</button>
    <label class="mcl-r" for="mclR">how strongly busy paths win: r = <b id="mclRv">2.8</b></label>
    <input type="range" id="mclR" min="1.2" max="6" step="0.1" value="2.8" aria-label="Inflation">
  </div>
  <div class="mcl-state" id="mclState"></div>
  <div class="mcl-grid">
    <div class="chart" id="mclGraph"></div>
    <div class="chart" id="mclMat"></div>
  </div>
  <div class="tiles">
    <div class="tile"><div class="lab">Families MCL sees at this moment</div><div class="big" id="mclNow"></div><div class="sub" id="mclNowS"></div></div>
    <div class="tile"><div class="lab">If "any line" meant "same family"</div><div class="big">1 family</div><div class="sub">all six glued together by the one stray line</div></div>
    <div class="tile"><div class="lab">Other settings of the knob r</div><div class="sub" id="mclSweep"></div><div class="sub">too gentle and families stay glued; stronger keeps them apart</div></div>
  </div>
</section>
<section>
  <h2>What is happening, step by step</h2>
  <ol class="steps">{steps}</ol>
</section>
<section>
  <h2>The same idea on real gene copies</h2>
  <div class="tbl"><table class="deftbl"><tbody>{params}</tbody></table></div>
  <p class="note">Why not simply "connected = same family": on our data, that rule glued families into giant groups of
  145 and 114 genes through stray links (register 474). The value r = 2.8 was not tuned to get a result: our anchored
  families come out the same anywhere from 2.0 to 4.0 (NPIP keeps all 44 copies). Soto's recipe does not use MCL: it calls
  "connected" one family, after first removing links between genes whose copy numbers differ by 2 or more.</p>
</section>
<script>
(function () {{
  const names = ["A1", "A2", "A3", "B1", "B2", "B3"];
  const E = [[0, 1, 0.95], [0, 2, 0.95], [1, 2, 0.95], [3, 4, 0.95], [3, 5, 0.95], [4, 5, 0.95], [2, 3, 0.72]];
  const n = names.length;
  const POS = [[70, 60], [70, 210], [190, 135], [330, 135], [450, 60], [450, 210]];
  const COLS = ["var(--s1)", "var(--s2)", "var(--s3)", "var(--s4)", "var(--s5)", "var(--s7)"];
  const norm = M => {{ for (let j = 0; j < n; j++) {{ let s = 0; for (let i = 0; i < n; i++) s += M[i][j]; if (s > 0) for (let i = 0; i < n; i++) M[i][j] /= s; }} return M; }};
  const base = () => {{ const M = Array.from({{ length: n }}, () => Array(n).fill(0)); E.forEach(([a, b, w]) => {{ M[a][b] = w; M[b][a] = w; }}); for (let i = 0; i < n; i++) M[i][i] = 1; return norm(M); }};
  const expand = M => {{ const R = Array.from({{ length: n }}, () => Array(n).fill(0)); for (let i = 0; i < n; i++) for (let j = 0; j < n; j++) {{ let s = 0; for (let k = 0; k < n; k++) s += M[i][k] * M[k][j]; R[i][j] = s; }} return R; }};
  const inflate = (M, r) => norm(M.map(row => row.map(v => {{ const x = Math.pow(v, r); return x < 1e-9 ? 0 : x; }})));
  const converged = (A, B) => A.every((row, i) => row.every((v, j) => Math.abs(v - B[i][j]) < 1e-7));
  const groups = M => {{
    const p = [...Array(n).keys()], f = x => {{ while (p[x] !== x) {{ p[x] = p[p[x]]; x = p[x]; }} return x; }};
    for (let j = 0; j < n; j++) {{ let r = 0; for (let i = 1; i < n; i++) if (M[i][j] > M[r][j]) r = i; const a = f(j), b = f(r); if (a !== b) p[Math.max(a, b)] = Math.min(a, b); }}
    const g = new Map(); for (let i = 0; i < n; i++) {{ const r = f(i); if (!g.has(r)) g.set(r, []); g.get(r).push(i); }}
    return [...g.values()];
  }};
  const runToEnd = r => {{ let M = base(), it = 0; while (it < 100) {{ const N = inflate(expand(M), r); it++; const done = converged(M, N); M = N; if (done) break; }} return {{ M, it }}; }};
  let M = base(), phase = "start", round = 0;
  const r = () => +document.getElementById("mclR").value;
  const fmtG = gs => gs.map(g => "{{" + g.map(i => names[i]).join(", ") + "}}").join("  ");
  function draw(label) {{
    const gs = groups(M), colOf = new Map();
    gs.forEach((g, k) => g.forEach(i => colOf.set(i, COLS[k % COLS.length])));
    let s = "";
    E.forEach(([a, b, w]) => {{
      const f = (M[a][b] + M[b][a]) / 2, [x1, y1] = POS[a], [x2, y2] = POS[b];
      s += `<line x1="${{x1}}" y1="${{y1}}" x2="${{x2}}" y2="${{y2}}" stroke="${{a === 2 && b === 3 ? "var(--cross)" : "var(--muted)"}}" stroke-width="${{(1 + 14 * f).toFixed(1)}}" stroke-linecap="round" opacity=".8"/>`;
      s += x1 === x2 ? `<text x="${{x1 + (x1 < 260 ? -14 : 14)}}" y="${{(y1 + y2) / 2 + 4}}" font-size="11" fill="var(--muted)" text-anchor="${{x1 < 260 ? "end" : "start"}}">${{w}}</text>` : `<text x="${{(x1 + x2) / 2}}" y="${{(y1 + y2) / 2 - 8}}" font-size="11" fill="var(--muted)" text-anchor="middle">${{w}}</text>`;
    }});
    names.forEach((nm, i) => {{
      const [x, y] = POS[i];
      s += `<circle cx="${{x}}" cy="${{y}}" r="20" fill="${{colOf.get(i)}}" stroke="var(--surface)" stroke-width="3"/><text x="${{x}}" y="${{y + 4.5}}" font-size="13" font-weight="700" fill="#fff" text-anchor="middle">${{nm}}</text>`;
    }});
    s += `<text x="260" y="262" font-size="11.5" fill="var(--muted)" text-anchor="middle">thicker line = more walkers using it now · same colour = same family so far</text>`;
    document.getElementById("mclGraph").innerHTML = `<svg viewBox="0 0 520 272" width="100%" role="img" aria-label="Toy homology graph with current flow">${{s}}</svg>`;
    const cs = 44, x0 = 46, y0 = 30;
    let h = "";
    names.forEach((nm, j) => {{ h += `<text x="${{x0 + j * cs + cs / 2}}" y="${{y0 - 10}}" font-size="12" font-weight="600" fill="var(--fg)" text-anchor="middle">${{nm}}</text>`; }});
    names.forEach((nm, i) => {{
      h += `<text x="${{x0 - 10}}" y="${{y0 + i * cs + cs / 2 + 4}}" font-size="12" font-weight="600" fill="var(--fg)" text-anchor="end">${{nm}}</text>`;
      for (let j = 0; j < n; j++) {{
        const v = M[i][j];
        h += `<rect x="${{x0 + j * cs + 1}}" y="${{y0 + i * cs + 1}}" width="${{cs - 2}}" height="${{cs - 2}}" rx="3" fill="var(--s1)" fill-opacity="${{(0.06 + 0.94 * v).toFixed(3)}}"/>`;
        h += `<text x="${{x0 + j * cs + cs / 2}}" y="${{y0 + i * cs + cs / 2 + 4}}" font-size="10.5" fill="${{v > 0.5 ? "#fff" : "var(--fg)"}}" text-anchor="middle">${{v < 0.005 ? "0" : Math.round(100 * v) + "%"}}</text>`;
      }}
    }});
    h += `<text x="${{x0 + 3 * cs}}" y="${{y0 + n * cs + 18}}" font-size="11.5" fill="var(--muted)" text-anchor="middle">each column: where the walkers that started</text>`;
    h += `<text x="${{x0 + 3 * cs}}" y="${{y0 + n * cs + 33}}" font-size="11.5" fill="var(--muted)" text-anchor="middle">on that copy are now (adds up to 100%)</text>`;
    document.getElementById("mclMat").innerHTML = `<svg viewBox="0 0 330 ${{y0 + n * cs + 42}}" width="100%" role="img" aria-label="Where the walkers are">${{h}}</svg>`;
    document.getElementById("mclState").textContent = label;
    document.getElementById("mclNow").textContent = gs.length + (gs.length === 1 ? " family" : " families");
    document.getElementById("mclNowS").textContent = fmtG(gs);
    document.getElementById("mclExp").disabled = phase === "expanded";
    document.getElementById("mclInf").disabled = phase !== "expanded";
  }}
  function start() {{ M = base(); phase = "start"; round = 0; draw("Start: walkers stand on every copy, ready to step along the lines (thicker lines are more tempting)"); }}
  document.getElementById("mclStart").addEventListener("click", start);
  document.getElementById("mclExp").addEventListener("click", () => {{ M = expand(M); round++; phase = "expanded"; draw(`Round ${{round}}: the walkers wandered; inside each trio they mix, and a few crossed the stray line`); }});
  document.getElementById("mclInf").addEventListener("click", () => {{
    const N = inflate(M, r()); M = N; phase = "inflated";
    draw(`Round ${{round}}: busy paths made busier (strength r = ${{r()}}); the trickle across the stray line shrinks`);
  }});
  document.getElementById("mclEnd").addEventListener("click", () => {{ const o = runToEnd(r()); M = o.M; phase = "done"; round = o.it; draw(`Finished after ${{o.it}} rounds (r = ${{r()}}): each copy sends nearly all its walkers to one place; copies sharing that place form a family`); }});
  document.getElementById("mclR").addEventListener("input", () => {{ document.getElementById("mclRv").textContent = r(); start(); }});
  document.getElementById("mclSweep").innerHTML = [1.2, 1.4, 2, 2.8, 4, 6].map(x => {{
    const o = runToEnd(x), gs = groups(o.M);
    return `r = ${{x}}: <b>${{gs.length}} famil${{gs.length === 1 ? "y" : "ies"}}</b>`;
  }}).join("<br>");
  start();
}})();
</script>
</div>'''


if __name__ == "__main__":
    print(pane())
