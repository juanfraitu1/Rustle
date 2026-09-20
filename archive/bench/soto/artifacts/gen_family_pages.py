#!/usr/bin/env python3
"""Generate 3 IGV-style family artifact pages from family_full.json."""
import json
FAM = json.load(open("bench/soto/artifacts/family_full.json"))

TEMPLATE = r"""<title>__TITLE__</title>
<style>
  :root{
    --ground:#f5f7f9; --panel:#ffffff; --panel-2:#eef2f5; --track:#e8edf1;
    --ink:#17222e; --muted:#5a6a79; --faint:#8496a5;
    --line:#dde4ea; --line-strong:#c6d0d9;
    --accent:#0e6e78; --accent-soft:#0e6e7818;
    --recovered:#157f5a; --recovered-soft:#157f5a1c;
    --missed:#b1462d; --missed-soft:#b1462d18;
    --discovery:#a9760f; --discovery-soft:#a9760f18;
    --cov:#4a6b86;
    --mono:ui-monospace,"SF Mono","JetBrains Mono","Cascadia Code",Menlo,Consolas,monospace;
    --sans:-apple-system,"Segoe UI","Helvetica Neue",system-ui,Roboto,sans-serif;
    --shadow:0 1px 2px #17222e0a,0 8px 24px #17222e0c;
  }
  @media (prefers-color-scheme:dark){:root{
    --ground:#0c1219; --panel:#121b24; --panel-2:#0f171f; --track:#0b1219;
    --ink:#e7eef4; --muted:#93a4b4; --faint:#67788a;
    --line:#213040; --line-strong:#2c3f52;
    --accent:#3cb7c4; --accent-soft:#3cb7c422;
    --recovered:#37c98c; --recovered-soft:#37c98c22;
    --missed:#e37a5d; --missed-soft:#e37a5d22;
    --discovery:#e0b24d; --discovery-soft:#e0b24d22;
    --cov:#7aa0c0;
    --shadow:0 1px 2px #0006,0 10px 30px #0004;
  }}
  :root[data-theme="dark"]{
    --ground:#0c1219; --panel:#121b24; --panel-2:#0f171f; --track:#0b1219;
    --ink:#e7eef4; --muted:#93a4b4; --faint:#67788a;
    --line:#213040; --line-strong:#2c3f52;
    --accent:#3cb7c4; --accent-soft:#3cb7c422;
    --recovered:#37c98c; --recovered-soft:#37c98c22;
    --missed:#e37a5d; --missed-soft:#e37a5d22;
    --discovery:#e0b24d; --discovery-soft:#e0b24d22;
    --cov:#7aa0c0; --shadow:0 1px 2px #0006,0 10px 30px #0004;
  }
  :root[data-theme="light"]{
    --ground:#f5f7f9; --panel:#ffffff; --panel-2:#eef2f5; --track:#e8edf1;
    --ink:#17222e; --muted:#5a6a79; --faint:#8496a5;
    --line:#dde4ea; --line-strong:#c6d0d9;
    --accent:#0e6e78; --accent-soft:#0e6e7818;
    --recovered:#157f5a; --recovered-soft:#157f5a1c;
    --missed:#b1462d; --missed-soft:#b1462d18;
    --discovery:#a9760f; --discovery-soft:#a9760f18;
    --cov:#4a6b86; --shadow:0 1px 2px #17222e0a,0 8px 24px #17222e0c;
  }
  *{box-sizing:border-box}
  body{margin:0;background:var(--ground);color:var(--ink);font-family:var(--sans);line-height:1.55;font-size:16px;-webkit-font-smoothing:antialiased}
  .wrap{max-width:1000px;margin:0 auto;padding:clamp(20px,4vw,48px) clamp(14px,4vw,36px) 80px}
  .eyebrow{font-family:var(--mono);font-size:12px;letter-spacing:.14em;text-transform:uppercase;color:var(--accent);font-weight:600;margin:0 0 10px;display:flex;gap:10px;align-items:center;flex-wrap:wrap}
  .regime{font-family:var(--mono);font-size:11px;letter-spacing:.06em;text-transform:uppercase;padding:3px 10px;border-radius:11px;font-weight:600}
  .r-all{background:var(--recovered-soft);color:var(--recovered)}
  .r-some{background:var(--discovery-soft);color:var(--discovery)}
  .r-none{background:var(--missed-soft);color:var(--missed)}
  h1{font-size:clamp(30px,5vw,46px);line-height:1.05;letter-spacing:-.02em;margin:0 0 14px;font-weight:700}
  h1 .sub{font-family:var(--mono);font-size:.42em;font-weight:500;color:var(--faint);letter-spacing:0}
  .lede{font-size:clamp(15.5px,1.9vw,18px);color:var(--muted);max-width:66ch;margin:0}
  .lede b{color:var(--ink);font-weight:600}

  .stats{display:flex;flex-wrap:wrap;gap:10px;margin:26px 0 6px}
  .stat{background:var(--panel);border:1px solid var(--line);border-radius:11px;padding:12px 18px;box-shadow:var(--shadow);min-width:118px}
  .stat .k{font-family:var(--mono);font-size:10.5px;letter-spacing:.1em;text-transform:uppercase;color:var(--faint);margin:0 0 5px}
  .stat .v{font-family:var(--mono);font-size:22px;font-weight:700;letter-spacing:-.01em;font-variant-numeric:tabular-nums}

  .section-h{font-family:var(--mono);font-size:12px;letter-spacing:.13em;text-transform:uppercase;color:var(--faint);margin:38px 0 14px;display:flex;align-items:center;gap:12px}
  .section-h::after{content:"";flex:1;height:1px;background:var(--line)}

  .legend{display:flex;flex-wrap:wrap;gap:16px;font-size:12.5px;color:var(--muted);font-family:var(--mono);margin-bottom:14px}
  .legend span{display:inline-flex;align-items:center;gap:6px}
  .sw{width:13px;height:13px;border-radius:3px;display:inline-block}
  .sw.d{background:var(--discovery);border:1.5px dashed var(--discovery)}

  .browser{background:var(--panel);border:1px solid var(--line);border-radius:13px;box-shadow:var(--shadow);overflow:hidden}
  .locus{border-bottom:1px solid var(--line)}
  .locus:last-child{border-bottom:none}
  .locus-head{font-family:var(--mono);font-size:12px;color:var(--muted);background:var(--panel-2);padding:8px 16px;border-bottom:1px solid var(--line);display:flex;justify-content:space-between;font-variant-numeric:tabular-nums}
  .locus-head b{color:var(--ink)}
  .igv{padding:12px 16px 4px}
  .cov{position:relative;height:58px;border-bottom:1px dashed var(--line)}
  .covbar{position:absolute;bottom:0;background:var(--cov);border-radius:2px 2px 0 0;min-width:5px;transform:translateX(-50%)}
  .covnum{position:absolute;bottom:0;transform:translateX(-50%);font-family:var(--mono);font-size:9.5px;color:var(--faint);white-space:nowrap}
  .covlabel{position:absolute;left:0;top:2px;font-family:var(--mono);font-size:9.5px;letter-spacing:.08em;text-transform:uppercase;color:var(--faint)}
  .feats{position:relative;height:26px;margin-top:8px}
  .feat{position:absolute;top:0;height:22px;border-radius:4px;min-width:11px;display:flex;align-items:center;justify-content:center;overflow:hidden}
  .feat.det{background:var(--recovered)}
  .feat.miss{background:var(--missed)}
  .feat.disc{background:var(--discovery-soft);border:1.5px dashed var(--discovery)}
  .feat .fin{font-family:var(--mono);font-size:9px;color:#fff;font-weight:600;white-space:nowrap;padding:0 3px}
  .labels{position:relative;height:70px;margin-top:2px}
  .lab{position:absolute;top:0;font-family:var(--mono);font-size:10.5px;transform-origin:left top;transform:rotate(34deg);white-space:nowrap;color:var(--muted)}
  .lab.det{color:var(--recovered)} .lab.miss{color:var(--missed)} .lab.disc{color:var(--discovery);font-weight:600}
  .lab .rr{color:var(--faint)}
  .ruler{position:relative;height:22px;border-top:1px solid var(--line-strong)}
  .tick{position:absolute;top:0;transform:translateX(-50%);font-family:var(--mono);font-size:10px;color:var(--faint);font-variant-numeric:tabular-nums}
  .tick::before{content:"";position:absolute;top:-1px;left:50%;width:1px;height:5px;background:var(--line-strong)}

  table{width:100%;border-collapse:collapse;font-size:14px;min-width:560px}
  .tablewrap{border:1px solid var(--line);border-radius:12px;overflow:hidden;box-shadow:var(--shadow);background:var(--panel)}
  .tscroll{overflow-x:auto}
  thead th{background:var(--panel-2);text-align:left;font-family:var(--mono);font-size:10.5px;letter-spacing:.07em;text-transform:uppercase;color:var(--faint);font-weight:600;padding:10px 14px;border-bottom:1px solid var(--line-strong);white-space:nowrap}
  thead th.num{text-align:right}
  tbody td{padding:9px 14px;border-bottom:1px solid var(--line);vertical-align:middle}
  tbody tr:last-child td{border-bottom:none}
  .g{font-family:var(--mono);font-size:13px;font-weight:600}
  .loc{font-family:var(--mono);font-size:12px;color:var(--muted);font-variant-numeric:tabular-nums}
  .num{text-align:right;font-family:var(--mono);font-variant-numeric:tabular-nums}
  .chip{display:inline-flex;align-items:center;gap:5px;font-family:var(--mono);font-size:11px;padding:2px 9px;border-radius:11px;font-weight:600}
  .chip.det{background:var(--recovered-soft);color:var(--recovered)}
  .chip.miss{background:var(--missed-soft);color:var(--missed)}
  .why{font-size:13px;color:var(--muted)}

  .disc-grid{display:grid;grid-template-columns:1fr;gap:12px}
  .disc{background:var(--panel);border:1px solid var(--line);border-left:3px solid var(--discovery);border-radius:12px;padding:16px 18px;box-shadow:var(--shadow)}
  .disc .top{display:flex;flex-wrap:wrap;gap:14px;align-items:baseline}
  .disc .co{font-family:var(--mono);font-size:14px;font-weight:600;font-variant-numeric:tabular-nums}
  .disc .ev{display:flex;flex-wrap:wrap;gap:8px;margin:10px 0 8px}
  .evchip{font-family:var(--mono);font-size:12px;padding:4px 11px;border-radius:8px;border:1px solid var(--line);background:var(--panel-2);font-variant-numeric:tabular-nums}
  .evchip b{color:var(--ink)} .evchip .lbl{color:var(--faint)}
  .evchip.reads b{color:var(--cov)} .evchip.homol b{color:var(--recovered)}
  .evchip.szout{border-color:var(--missed)} .evchip.szout b{color:var(--missed)}
  .disc .note{font-size:13.5px;color:var(--muted);margin:0}
  .verdict{margin-left:auto;font-family:var(--mono);font-size:11px;letter-spacing:.05em;text-transform:uppercase;padding:3px 10px;border-radius:11px;font-weight:600;background:var(--panel-2);color:var(--muted);border:1px solid var(--line)}

  .fpbox{margin:14px 0 0;background:var(--accent-soft);border:1px solid var(--accent-soft);border-left:3px solid var(--accent);border-radius:11px;padding:14px 18px;font-size:14px;color:var(--muted)}
  .fpbox b{color:var(--ink)}
  footer{margin-top:38px;padding-top:18px;border-top:1px solid var(--line);color:var(--faint);font-size:12.5px;line-height:1.5}
  footer code{font-family:var(--mono);font-size:11.5px;color:var(--muted)}
</style>

<div class="wrap">
  <p class="eyebrow"><span id="famtag"></span><span class="regime" id="regtag"></span></p>
  <h1 id="title"></h1>
  <p class="lede" id="lede"></p>

  <div class="stats" id="stats"></div>

  <div class="section-h">Genome browser view · read coverage over the locus</div>
  <div class="legend">
    <span><span class="sw" style="background:var(--recovered)"></span>recovered member</span>
    <span><span class="sw" style="background:var(--missed)"></span>missed member</span>
    <span><span class="sw d"></span>discovered copy</span>
    <span><span class="sw" style="background:var(--cov)"></span>read coverage (log)</span>
  </div>
  <div class="browser" id="browser"></div>

  <div class="section-h">Members</div>
  <div class="tablewrap"><div class="tscroll"><table>
    <thead><tr><th>Gene</th><th>Locus</th><th class="num">bp</th><th class="num">reads</th><th>Status</th><th>Why</th></tr></thead>
    <tbody id="mtb"></tbody>
  </table></div></div>

  <div class="section-h" id="disc-h">Discovered copies of this family · is each one real, or a false positive?</div>
  <div class="fpbox" id="disc-fp">A discovery is a <b>false positive only if it has no evidence</b> — no reads and no homology. Every copy below is a candidate <b>member of this family</b>, judged on exactly two things you can read straight off the browser: <b>read coverage</b> (is it expressed?) and <b>full-length identity to a family member</b> (is it a real copy of <i>this</i> family?). A phantom would show neither.</div>
  <div class="disc-grid" id="disc" style="margin-top:14px"></div>

  <footer>
    Reads = primary alignments (<code>samtools view -F 2308</code>) over each locus in gorilla A119b IsoSeq on T2T-CHM13 v2.0.
    Discovery identity = best <code>minimap2 asm20</code> hit to any Soto member. Source:
    <code>soto_member_detection.tsv</code>, <code>candidate_verification.tsv</code>.
  </footer>
</div>

<script>
const FAM = /*__FAM__*/;
const CLUSTER_GAP = 600000;
const fmt = n => n.toLocaleString('en-US');
const mb = n => (n/1e6).toFixed(3);

// build feature list
const feats = [];
FAM.members.forEach(m=>feats.push({...m, type:m.det?'det':'miss'}));
FAM.discoveries.forEach(d=>feats.push({gene:'new copy', start:d.start, end:d.end, reads:d.reads, det:false, type:'disc', homolog:d.homolog, id:d.id}));
feats.sort((a,b)=>a.start-b.start);
const maxReads = Math.max(...feats.map(f=>f.reads), 1);

// cluster
const clusters=[]; let cur=null;
feats.forEach(f=>{
  if(!cur || f.start - cur.hiRaw > CLUSTER_GAP){ cur={feats:[f], loRaw:f.start, hiRaw:f.end}; clusters.push(cur); }
  else { cur.feats.push(f); cur.hiRaw=Math.max(cur.hiRaw,f.end); }
});
clusters.forEach(c=>{
  const span=c.hiRaw-c.loRaw;
  const pad = c.feats.length===1 ? 14000 : Math.max(span*0.07, 2500);
  c.lo=c.loRaw-pad; c.hi=c.hiRaw+pad; c.span=c.hi-c.lo;
});

function pos(c,x){ return (x-c.lo)/c.span*100; }
function covH(r){ return 6 + 46*Math.log(r+1)/Math.log(maxReads+1); }

// render browser
const bro=document.getElementById('browser');
bro.innerHTML = clusters.map(c=>{
  const cov = c.feats.map(f=>{
    const cx=pos(c,(f.start+f.end)/2), w=Math.max(pos(c,f.end)-pos(c,f.start),0.8), h=covH(f.reads);
    return `<div class="covbar" style="left:${cx}%;width:${w}%;height:${h}px;${f.type==='disc'?'background:var(--discovery);':''}"></div>
            <div class="covnum" style="left:${cx}%;bottom:${h+2}px">${fmt(f.reads)}</div>`;
  }).join('');
  const blocks = c.feats.map(f=>{
    const l=pos(c,f.start), w=Math.max(pos(c,f.end)-pos(c,f.start),0);
    return `<div class="feat ${f.type}" style="left:${l}%;width:${w}%" title="${f.gene} · ${fmt(f.reads)} reads"></div>`;
  }).join('');
  const labs = c.feats.map(f=>{
    const l=pos(c,(f.start+f.end)/2);
    const nm = f.type==='disc' ? '◇ new copy' : f.gene;
    return `<div class="lab ${f.type}" style="left:${l}%">${nm} <span class="rr">${fmt(f.reads)}r</span></div>`;
  }).join('');
  const nt=5; let ticks='';
  for(let i=0;i<nt;i++){ const p=i/(nt-1); const gpos=c.lo+p*c.span; ticks+=`<div class="tick" style="left:${p*100}%">${mb(gpos)}</div>`; }
  return `<div class="locus">
    <div class="locus-head"><span><b>${FAM.chrom}</b>:${fmt(Math.round(c.lo))}–${fmt(Math.round(c.hi))}</span><span>${c.feats.length} feature${c.feats.length>1?'s':''} · ${(c.span/1000).toFixed(0)} kb</span></div>
    <div class="igv">
      <div class="cov"><span class="covlabel">coverage</span>${cov}</div>
      <div class="feats">${blocks}</div>
      <div class="labels">${labs}</div>
      <div class="ruler">${ticks}<div class="tick" style="left:100%;transform:translateX(-100%)">Mb</div></div>
    </div>
  </div>`;
}).join('');

// header
document.getElementById('title').innerHTML = `${FAM.gene} <span class="sub">${FAM.id} · ${FAM.chrom}</span>`;
document.getElementById('famtag').textContent = `Soto family · ${FAM.chrom}`;
const reg=document.getElementById('regtag'); reg.textContent = {all:'Discovered all',some:'Discovered some',none:'Discovered none'}[FAM.regime];
reg.className='regime r-'+FAM.regime;
document.getElementById('lede').innerHTML = FAM.lede;
const totReads=FAM.members.reduce((s,m)=>s+m.reads,0);
document.getElementById('stats').innerHTML = [
  ['Sensitivity',FAM.sens],['Precision',FAM.prec],['Member reads',fmt(totReads)],['Discoveries','+'+FAM.discoveries.length]
].map(([k,v])=>`<div class="stat"><p class="k">${k}</p><div class="v">${v}</div></div>`).join('');

// member table
document.getElementById('mtb').innerHTML = FAM.members.map(m=>`<tr>
  <td class="g">${m.gene}</td>
  <td class="loc">${FAM.chrom}:${fmt(m.start)}-${fmt(m.end)}</td>
  <td class="num">${fmt(m.end-m.start)}</td>
  <td class="num">${fmt(m.reads)}</td>
  <td><span class="chip ${m.det?'det':'miss'}">${m.det?'recovered':'missed'}</span></td>
  <td class="why">${m.why}</td></tr>`).join('');

// discovery cards — only genuine candidate members OF THIS FAMILY
if(!FAM.discoveries.length){
  document.getElementById('disc-h').textContent = 'Discovered copies of this family';
  document.getElementById('disc-fp').style.display = 'none';
  document.getElementById('disc').innerHTML =
    `<div class="disc" style="border-left-color:var(--faint)"><p class="note" style="margin:0">No new copies of <b>this</b> family were discovered — every recovered call is one of its catalogued members. An off-catalogue locus that sits nearby but is <b>not homologous to this family</b> is not a member and is not shown here.</p></div>`;
} else
document.getElementById('disc').innerHTML = FAM.discoveries.map(d=>{
  const L=d.end-d.start;
  const sizeChip = `<span class="evchip ${d.outlier?'szout':''}"><span class="lbl">size</span> <b>${(L/1000).toFixed(1)} kb</b> · <b>${d.ratio.toFixed(1)}×</b> <span class="lbl">members</span>${d.outlier?' ⚠':''}</span>`;
  return `<div class="disc" ${d.outlier?'style="border-left-color:var(--missed)"':''}>
  <div class="top">
    <span class="co">${FAM.chrom}:${fmt(d.start)}-${fmt(d.end)}</span>
    <span class="verdict">candidate for ${FAM.gene}</span>
  </div>
  <div class="ev">
    <span class="evchip reads"><span class="lbl">reads</span> <b>${fmt(d.reads)}</b></span>
    <span class="evchip homol"><span class="lbl">identity to a member</span> <b>${(d.id*100).toFixed(1)}%</b></span>
    ${sizeChip}
  </div>
  <p class="note">${d.note}</p>
</div>`;}).join('');
</script>
"""

TITLES = {
 "ID_8":"PMS2P — discovered all (Rustle × Soto)",
 "ID_131":"Amylase — discovered some (Rustle × Soto)",
 "ID_402":"NCF1 — discovered none (Rustle × Soto)",
}
FN = {"ID_8":"family_pms2p.html","ID_131":"family_amylase.html","ID_402":"family_ncf1.html"}
for fid,f in FAM.items():
    html = TEMPLATE.replace("__TITLE__", TITLES[fid]).replace("/*__FAM__*/", json.dumps(f))
    open("bench/soto/artifacts/"+FN[fid],"w").write(html)
    print("wrote", FN[fid])
