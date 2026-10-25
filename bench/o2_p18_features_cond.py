import sys
src=open(__file__.replace('_cond', '')).read().split("W = collections.defaultdict")[0]
exec(src)
c=collections.Counter()
for name, rs in recs.items():
    t = truth[name]; p = next(r for r in rs if not r.is_secondary)
    band = [r for r in rs if r.get_tag('ms') >= p.get_tag('ms') - BAND]
    if len(band) < 2 or right(p,t): continue
    for r in band:
        if r.query_sequence is None:
            r.query_sequence = p.query_sequence if r.is_reverse == p.is_reverse else rc(p.query_sequence)
    good=[r for r in band if right(r,t)]
    if not good: c[('primary wrong, truth not in band',)]+=1; continue
    g=max(good,key=lambda r:r.get_tag('ms')); fp,fg=feats(p),feats(g)
    grp='tie' if abs(fp['ll']-fg['ll'])<1e-9 else 'decidable'
    c[(grp,'reads')]+=1; c[(grp,'allele',t[6])]+=1
    for k in ('ms','NM','ll','annot_junctions','noncanon_junctions','downstream_A','cm','s1','introns','de'):
        a,b=(fg[k],fp[k]) if k not in ('NM','de') else (fp[k],fg[k])
        c[(grp,k,'favours truth' if a>b else 'favours primary' if a<b else 'equal')]+=1
for k,v in sorted(c.items(), key=str): print(k,v)
