/* ---------- tabs ---------- */
(function(){
  const TABS=[['tb-browser','tab-browser','browser'],['tb-seeds','tab-seeds','seeds']];
  function show(key,focus){
    TABS.forEach(([b,p,k])=>{const on=k===key,btn=document.getElementById(b);btn.setAttribute('aria-selected',on);btn.tabIndex=on?0:-1;document.getElementById(p).hidden=!on;if(on&&focus)btn.focus();});
    try{localStorage.setItem('lrpap1-tab',key);}catch(e){}
    try{history.replaceState(null,'',key==='seeds'?'#seeds':location.pathname+location.search);}catch(e){}
  }
  TABS.forEach(([b,p,k],i)=>{const btn=document.getElementById(b);btn.addEventListener('click',()=>show(k));
    btn.addEventListener('keydown',e=>{if(e.key==='ArrowRight'||e.key==='ArrowLeft'){e.preventDefault();show(TABS[(i+(e.key==='ArrowRight'?1:TABS.length-1))%TABS.length][2],true);}});});
  let start='browser';try{start=localStorage.getItem('lrpap1-tab')||'browser';}catch(e){}
  if(location.hash==='#seeds')start='seeds';
  show(start==='seeds'?'seeds':'browser');
})();

/* ---------- seeds tab ---------- */
(function(){
  const S=SEEDS,ids=S.ids,n=ids.length;
  const byId=Object.fromEntries(DATA.map(c=>[c.cid,c]));
  const lab=c=>`${c.chrom} ${(c.lo/1e6).toFixed(2)}`;
  const idv=(i,j)=>i===j?1:(S.identity[i+'|'+j]||S.identity[j+'|'+i]||null);
  const pair=[];ids.forEach(i=>ids.forEach(j=>{if(i!==j&&idv(i,j))pair.push(idv(i,j));}));
  const lo=Math.min(...pair),hi=Math.max(...pair);
  const allFound=ids.every(i=>S.finds[i].length===n);
  const extra=ids.reduce((a,i)=>a+S.extra_loci[i],0);
  const firstTwo=S.sweep.find(r=>r.min_seeds_iterate>1);
  const t=[[`${S.min_seeds_one_pass}`,'seed needed: any one of the eleven'],[`${allFound?n:'—'} / ${n}`,'copies found from every seed'],[`${extra}`,'hits outside the eleven copies'],
    [`${(lo*100).toFixed(1)}–${(hi*100).toFixed(1)}%`,'identity between copies'],[`${(firstTwo.threshold*100).toFixed(2)}%`,'identity bar where a 2nd seed is needed (shipped: 80%)']];
  document.getElementById('seed-tiles').innerHTML=t.map(([v,l])=>`<div class="tile"><div class="n">${v}</div><div class="l">${l}</div></div>`).join('');

  /* 4 · heatmap */
  const D0=0.955,shade=v=>Math.round(18+82*Math.max(0,Math.min(1,(v-D0)/(1-D0))));
  const lg=document.getElementById('hm-legend');
  lg.innerHTML=`<span>gene-body identity</span>`+[0.96,0.97,0.98,0.99,1].map(v=>`<span><i class="sw" style="background:color-mix(in oklab,var(--fib) ${shade(v)}%,var(--surface))"></i>${(v*100).toFixed(0)}%</span>`).join('')+`<span><i class="sw self"></i>the seed itself</span><span><i class="sw corner"></i>found by the gene-body search only</span>`;
  {const host=document.getElementById('hm');const CW=44,CH=30,X0=138,Y0=104,W=X0+n*CW+48,H=Y0+n*CH+28;const s=svg(W,H,host);
    el('text',{x:X0+n*CW/2,y:12,class:'rowhd','text-anchor':'middle'},s,'copies found');
    el('text',{x:0,y:Y0-10,class:'rowhd'},s,'seed copy');
    ids.forEach((j,c)=>{const x=X0+c*CW+CW/2,cp=byId[j];const tx=el('text',{x:x,y:Y0-8,class:'tk',transform:`rotate(-50 ${x} ${Y0-8})`},s,lab(cp));bindTip(tx,cp.name);
      if(cp.cls==='fragment')el('rect',{x:X0+c*CW+1,y:Y0-2,width:CW-2,height:n*CH+4,rx:3,style:'fill:var(--tes-soft);opacity:.35'},s);});
    ids.forEach((i,r)=>{const y=Y0+r*CH,ci=byId[i];const rt=el('text',{x:0,y:y+CH/2+4,class:'lab'},s,lab(ci));bindTip(rt,ci.name+(ci.cls==='fragment'?' · 5′ fragment':''));
      const body=new Set(S.finds_body[i]),tr=new Set(S.finds_transcript[i]);
      ids.forEach((j,c)=>{const found=S.finds[i].includes(j),v=idv(i,j),x=X0+c*CW,cj=byId[j];
        const r2=el('rect',{x:x+2,y:y+2,width:CW-4,height:CH-4,rx:3,class:'hcell'+(i===j?' self':''),style:found?`fill:color-mix(in oklab,var(--fib) ${shade(v||D0)}%,var(--surface))`:'fill:none'},s);
        if(found&&body.has(j)&&!tr.has(j))el('path',{d:`M${x+CW-12},${y+2} L${x+CW-2},${y+2} L${x+CW-2},${y+12} Z`,class:'corner'},s);
        bindTip(r2,`seed ${ci.name} (${lab(ci)}) → ${cj.name} (${lab(cj)})\n${i===j?'the seed itself':found?'found':'not found'}${v&&i!==j?` · gene-body identity ${(v*100).toFixed(2)}%`:''}\nby ${[tr.has(j)?'transcript':null,body.has(j)?'gene body':null].filter(Boolean).join(' + ')||'—'}`);});
      el('text',{x:X0+n*CW+10,y:y+CH/2+4,class:'num'},s,S.finds[i].length);});
    el('text',{x:X0+n*CW+10,y:Y0-8,class:'tk'},s,'found');
  }

  /* 5 · small multiples: copies found vs number of seeds */
  const PANELS=[['0.800','identity ≥ 80% · shipped'],['0.980','identity ≥ 98%'],['0.990','identity ≥ 99%'],['0.995','identity ≥ 99.5%']];
  const sm=document.getElementById('sm');
  PANELS.forEach(([key,title])=>{const cur=S.curves[key].one_pass;const box=document.createElement('div');box.className='fig';sm.appendChild(box);
    const W=300,H=212,L=34,R=10,T=26,B=34,pw=W-L-R,ph=H-T-B;const s=svg(W,H,box);
    const X=k=>L+(k-1)/(n-1)*pw,Y=v=>T+ph-(v/n)*ph;
    el('text',{x:L,y:14,class:'num'},s,title);
    [0,3,6,9,11].forEach(v=>{el('line',{x1:L,x2:L+pw,y1:Y(v),y2:Y(v),class:v===n?'allline':'grid'},s);el('text',{x:L-6,y:Y(v)+4,class:'tk','text-anchor':'end'},s,v);});
    [1,3,5,7,9,11].forEach(k=>el('text',{x:X(k),y:T+ph+14,class:'tk','text-anchor':'middle'},s,k));
    el('text',{x:L+pw/2,y:H-4,class:'tk','text-anchor':'middle'},s,'seeds used (k)');
    let band='M'+cur.map(r=>`${X(r.k)},${Y(r.max)}`).join('L')+'L'+[...cur].reverse().map(r=>`${X(r.k)},${Y(r.min)}`).join('L')+'Z';
    el('path',{d:band,class:'band'},s);
    el('path',{d:'M'+cur.map(r=>`${X(r.k)},${Y(r.median)}`).join('L'),class:'medline'},s);
    cur.forEach(r=>{const hit=el('rect',{x:X(r.k)-pw/(2*(n-1)),y:T,width:pw/(n-1),height:ph,style:'fill:transparent'},s);
      bindTip(hit,`${title}\n${r.k} seed${r.k>1?'s':''} (${fmt(r.n_sets)} sets): ${r.min}–${r.max} copies found, median ${r.median}\n${Math.round(r.all_found*100)}% of the sets find all 11`);});
    const k1=cur[0];el('text',{x:W-R,y:14,class:'tk','text-anchor':'end'},s,k1.max===n&&k1.min===n?'1 seed: all 11':`1 seed: ${k1.min}–${k1.max} copies`);});

  /* 6 · seeds needed vs identity bar (piecewise x: 80–95% compressed, 95–100% wide) */
  {const host=document.getElementById('sw');const W=820,H=300,L=46,R=118,T=18,B=40,pw=W-L-R,ph=H-T-B,BRK=0.95,F=0.22;const s=svg(W,H,host);
    const X=v=>v<=BRK?L+(v-0.8)/(BRK-0.8)*pw*F:L+pw*F+(v-BRK)/(1-BRK)*pw*(1-F);const Y=k=>T+ph-(k-1)/(n-1)*ph;
    el('rect',{x:X(lo),y:T,width:X(hi)-X(lo),height:ph,class:'idband'},s);
    [1,3,5,7,9,11].forEach(k=>{el('line',{x1:L,x2:L+pw,y1:Y(k),y2:Y(k),class:'grid'},s);el('text',{x:L-8,y:Y(k)+4,class:'tk','text-anchor':'end'},s,k);});
    [0.80,0.85,0.90,0.95,0.96,0.97,0.98,0.99,1.00].forEach(v=>el('text',{x:X(v),y:T+ph+15,class:'tk','text-anchor':'middle'},s,(v*100).toFixed(0)+'%'));
    el('path',{d:`M${X(BRK)-4},${T+ph-5} l4,10 M${X(BRK)},${T+ph-5} l4,10`,class:'brk'},s);
    el('text',{x:L+pw/2,y:H-4,class:'tk','text-anchor':'middle'},s,"finder identity bar (axis compressed below 95%)");
    el('text',{x:12,y:T+ph/2,class:'tk','text-anchor':'middle',transform:`rotate(-90 12 ${T+ph/2})`},s,'seeds needed for all 11');
    el('line',{x1:X(0.8),x2:X(0.8),y1:T,y2:T+ph,class:'ship'},s);el('text',{x:X(0.8)+5,y:T+12,class:'tk'},s,'shipped bar 80%');
    const step=key=>{let d='',prev=null;S.sweep.forEach(r=>{const x=X(r.threshold),y=Y(r[key]);d+=prev===null?`M${x},${y}`:`L${x},${Y(prev)}L${x},${y}`;prev=r[key];});return d;};
    el('path',{d:step('min_seeds_one_pass'),class:'step op'},s);el('path',{d:step('min_seeds_iterate'),class:'step it'},s);
    const last=S.sweep[S.sweep.length-1];
    el('text',{x:L+pw+8,y:Y(last.min_seeds_iterate)+4,class:'num'},s,'11 at 100%');
    const a=S.sweep.find(r=>r.min_seeds_iterate>1);el('line',{x1:X(a.threshold),x2:X(a.threshold),y1:Y(1),y2:Y(4.2),class:'ann'},s);
    el('text',{x:X(a.threshold)-6,y:Y(4.6),class:'tk','text-anchor':'end'},s,`2nd seed needed from ${(a.threshold*100).toFixed(2)}%`);
    S.sweep.forEach(r=>{const hit=el('rect',{x:X(r.threshold)-3,y:T,width:6,height:ph,style:'fill:transparent'},s);
      bindTip(hit,`identity bar ${(r.threshold*100).toFixed(2)}%\nseeds needed: ${r.min_seeds_one_pass} in one pass · ${r.min_seeds_iterate} with iteration\n${r.components.length} subcluster${r.components.length>1?'s':''} need their own seed`);});
  }
})();
