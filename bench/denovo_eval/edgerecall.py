"""EDGE-HUMAN pre-MCL edge recall: Compara Simiiformes pairs (expressed, judgeable) connected by the shipped families vs families + U1 lift edges."""
import sys,itertools,collections,csv,json
sys.path.insert(0,'/mnt/linuxdisk/home/juanfraitu/rustle_m2_soto/bench/denovo_eval'); import famtruth as F
s=sys.argv[1]; R=f'/mnt/linuxdisk/tmp/rustle_figures/runs/{s}/{s}'
G=F.load_cat(); idx=F.exon_index(G)
# loci -> gene via gff3 exons
gene_of_key={}; ex=collections.defaultdict(list); keyinfo={}
for l in open(R+'.fam.loci.gff3'):
    f=l.rstrip('\n').split('\t')
    if len(f)<9: continue
    a=dict(x.split('=',1) for x in f[8].split(';') if '=' in x)
    if f[2]=='gene': keyinfo[a['Name']]=(f'{f[0]}:{f[3]}-{f[4]}',f[0],f[6])
    elif f[2]=='exon': ex[a.get('gene') or a['Parent'].replace('gene-','')].append((int(f[3])-1,int(f[4])))
for n,(k,c,st) in keyinfo.items():
    g=F.map_copy(idx,c,st,sorted(ex.get(n,[])) or [(int(k.split(':')[1].split('-')[0])-1,int(k.split('-')[1]))])
    if g: gene_of_key[k]=g
fams,_,_=F.load_catalog(R+'.fam.copies.tsv',G,idx)
base_pairs={(a,b) for v in fams.values() for a,b in itertools.combinations(sorted(v),2)}
u1=set()
for l in open(f'{s}_u1/edges.keyed.tsv'):
    a,b,w=l.split('\t'); ga,gb=gene_of_key.get(a),gene_of_key.get(b)
    if ga and gb and ga!=gb: u1.add((min(ga,gb),max(ga,gb)))
expressed=F.expressed_from_gtf(R+'.families.gtf',idx) | set().union(*fams.values())
sim,universe=F.truth_pairs(G,F.PARALOGS,F.LEVELS['simiiformes']); sim={tuple(sorted(p)) for p in sim}
ue=universe & expressed; truth={p for p in sim if p[0] in ue and p[1] in ue}
# connectivity: components of base-family pairs (families are cliques) ∪ U1 edges, restricted to truth evaluation
def comp_pairs(edges):
    par={}
    def f(x):
        while par.setdefault(x,x)!=x: par[x]=par[par[x]]; x=par[x]
        return x
    for a,b in edges: par[f(a)]=f(b)
    return f
fb=comp_pairs(base_pairs); fu=comp_pairs(base_pairs|u1)
rb=sum(1 for a,b in truth if fb(a)==fb(b)); ru=sum(1 for a,b in truth if fu(a)==fu(b)); direct=sum(1 for p in truth if p in u1)
# judge the NEW gene pairs U1 adds (direct edges not already in a base family)
anyp=set()
for r in csv.DictReader(open(F.PARALOGS),delimiter='\t'):
    a,b=r['Gene stable ID'].split('.')[0],r['Human paralogue gene stable ID'].split('.')[0]
    if a and b: anyp.add((min(a,b),max(a,b)))
groups=collections.defaultdict(set)
for r in csv.DictReader(open('/mnt/linuxdisk/tmp/hgnc/hgnc_complete_set.txt'),delimiter='\t'):
    for g in r['gene_group_id'].split('|'):
        if g and r['ensembl_gene_id']: groups[r['ensembl_gene_id']].add(g)
def related(a,b):
    if G[a]['name']==G[b]['name']: return 'same-symbol'
    ea,eb={e.split('.')[0] for e in G[a]['ens']},{e.split('.')[0] for e in G[b]['ens']}
    if any((min(x,y),max(x,y)) in anyp for x in ea for y in eb): return 'older'
    if any(groups[x]&groups[y] for x in ea for y in eb): return 'hgnc'
    return None
new=u1-base_pairs; c=collections.Counter(); wrong=[]
for a,b in new:
    if a in universe and b in universe:
        if (a,b) in sim: c['true']+=1
        else:
            k=related(a,b); c[k or 'wrong']+=1
            if not k: wrong.append((G[a]['name'],G[b]['name']))
    else: c['unjudged']+=1
print(f"{s}: truth pairs (expressed, judgeable) {len(truth)} | connected by shipped families {rb} ({rb/len(truth):.3f}) | + U1 lift edges {ru} ({ru/len(truth):.3f}); U1 direct hits {direct}")
print(f"  U1 new gene pairs {len(new)}: true {c['true']} older {c['older']} hgnc {c['hgnc']} same-symbol {c['same-symbol']} WRONG {c['wrong']} unjudged {c['unjudged']} | judged precision (true+older+hgnc+same)/(judged) {(c['true']+c['older']+c['hgnc']+c['same-symbol'])/max(1,len(new)-c['unjudged']):.3f}")
print('  wrong examples', sorted(wrong)[:12])
json.dump(dict(truth=len(truth),base=rb,u1=ru,direct=direct,new=dict(c)),open(f'{s}_u1/edgerecall.json','w'))
