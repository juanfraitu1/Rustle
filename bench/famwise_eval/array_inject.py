"""Add the read classes seen in A119b arrays to a famsim scenario before alignment (appends to reads.fq):
antisense (reverse complement of a unit's genomic span +-flanks), read-through (unit i transcript + unit i+1 transcript),
run-on (sense genomic span from the transcript start to 1.5 kb past its end). Names carry the class and the unit."""
import sys, random, collections, re
d=sys.argv[1]; rnd=random.Random(5)
seqs={}; name=None
for l in open(f'{d}/genome.truth.fa'):
    l=l.strip()
    if l.startswith('>'): name=l[1:].split()[0]; seqs[name]=[]
    else: seqs[name].append(l)
seqs={k:''.join(v) for k,v in seqs.items()}
ex=collections.defaultdict(list); strand={}; ctg={}
for l in open(f'{d}/truth.gtf'):
    f=l.rstrip('\n').split('\t')
    if len(f)>8 and f[2]=='exon':
        t=re.search(r'transcript_id "([^"]+)"',f[8]).group(1); ex[t].append((int(f[3])-1,int(f[4]))); strand[t]=f[6]; ctg[t]=f[0]
expr=[l.split('\t')[0] for l in open(f'{d}/copies.tsv') if l.startswith('U') and int(l.split('\t')[8])>0]
iso0={t.split('.')[0]:sorted(ex[t]) for t in ex if t.endswith('.iso0') and t.split('.')[0] in expr}
rc=lambda s: s[::-1].translate(str.maketrans('ACGTacgtN','TGCAtgcaN'))
def mut(s,e=0.002): return ''.join(c if rnd.random()>e else rnd.choice('ACGT') for c in s)
out=open(f'{d}/reads.fq','a'); n=collections.Counter()
units=sorted(iso0, key=lambda u: iso0[u][0][0])
allu=sorted({t.split('.')[0] for t in ex if t.startswith('U')}, key=lambda u: min(e[0][0] for t,e in ex.items() if t.startswith(u+'.')))
for u in units:
    E=iso0[u]; c='sim'; s0,e0=E[0][0],E[-1][1]
    for k in range(int(0.2*100)):          # antisense, 20 per expressed unit
        a=s0-rnd.randint(0,500); b=e0+rnd.randint(300,900)
        q=mut(rc(seqs[c][a:b])); out.write(f"@anti_{u}|x|{k}\n{q}\n+\n{'I'*len(q)}\n"); n['antisense']+=1
    for k in range(10):                     # run-on: sense genomic span past the end
        a=s0+rnd.randint(0,100); b=e0+rnd.randint(800,1500)
        q=mut(seqs[c][a:b]); out.write(f"@runon_{u}|x|{k}\n{q}\n+\n{'I'*len(q)}\n"); n['run_on']+=1
    j=allu.index(u)
    if j+1<len(allu):                       # read-through into the next unit (its canonical transcript)
        nxt=[t for t in ex if t==allu[j+1]+'.iso0'][0]
        tx=''.join(seqs[c][a:b] for a,b in E)+''.join(seqs[c][a:b] for a,b in sorted(ex[nxt]))
        for k in range(8):
            q=mut(tx); out.write(f"@thru_{u}|x|{k}\n{q}\n+\n{'I'*len(q)}\n"); n['read_through']+=1
out.close(); print(dict(n))
