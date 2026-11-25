# audit round 60: internal consistency of the O1 and O2 tables of one run (PREFIX given). Prints one line per check.
import csv, sys, collections
P = sys.argv[1]
rd = lambda f: list(csv.DictReader(open(P + f), delimiter='\t'))
O2 = __import__('os').path.exists(P + '.assign.assignments.tsv')   # since register 1903 `all` runs the family report, not `assign`
cop, clu = rd('.fam.copies.tsv'), rd('.fam.clusters.tsv')
asg, q, fam = (rd('.assign.assignments.tsv'), rd('.assign.quant.tsv'), rd('.assign.families.tsv')) if O2 else ([], [], [])
skip = set(l.strip() for l in open(P + '.assign.excluded.txt')) if __import__('os').path.exists(P + '.assign.excluded.txt') else set()
def chk(name, bad, n): print(f"{'ok  ' if not bad else 'FAIL'} {name}: {len(bad)} of {n}" + (f"  e.g. {list(bad)[:3]}" if bad else ''))
# O1
csize = collections.Counter(r['cluster_id'] for r in clu); decl = {r['cluster_id']: int(r['size']) for r in clu}
chk('clusters: declared size = member rows', [c for c in csize if csize[c] != decl[c]], len(csize))
chk('clusters: every cluster has >= 2 members', [c for c in csize if csize[c] < 2], len(csize))
mem = collections.Counter((r['chrom'], r['start'], r['end']) for r in clu)
chk('clusters: a member locus is in one cluster', [m for m in mem if mem[m] > 1], len(mem))
tid = collections.Counter(r['tid'] for r in cop)
chk('copies: transcript id unique', [t for t in tid if tid[t] > 1], len(tid))
key = collections.Counter((r['family_id'], r['copy_idx']) for r in cop)
chk('copies: (family, copy_idx) unique', [k for k in key if key[k] > 1], len(key))
byf = collections.defaultdict(list)
for r in cop: byf[r['family_id']].append(int(r['copy_idx']))
chk('copies: copy_idx is 0..n-1 per family', [f for f, v in byf.items() if sorted(v) != list(range(len(v)))], len(byf))
chk('copies: every family has >= 2 copies', [f for f, v in byf.items() if len(v) < 2], len(byf))
chk('copies: family is a cluster', [f for f in byf if f not in csize], len(byf))
def exbad(r):
    ex = [tuple(map(int, x.split('-'))) for x in r['exons'].split(',')]
    return len(ex) != int(r['n_exon']) or ex[0][0] != int(r['start']) or ex[-1][1] != int(r['end']) or any(a >= b for a, b in ex) or any(ex[i][1] > ex[i + 1][0] for i in range(len(ex) - 1))
chk('copies: exons sorted, disjoint, span = start..end, count = n_exon', [r['tid'] for r in cop if exbad(r)], len(cop))
chk('copies: strand is + or -', [r['tid'] for r in cop if r['strand'] not in '+-' or not r['strand']], len(cop))
hdr = [l[1:].strip() for l in open(P + '.fam.copies.fa') if l.startswith('>')]
want = [f"{r['family_id']}|{r['copy_idx']}|{r['chrom']}:{r['start']}-{r['end']}|{r['strand']}|nexon={r['n_exon']}" for r in cop]
chk('copies.fa: one record per copy row, same order', [] if hdr == want else ['order or content differs'], len(cop))
# O2
if not O2:
    print(f"info {P}.assign.assignments.tsv is absent (the legacy `assign` stage was not run): the family tables were checked, the assignment checks skipped"); sys.exit(0)
ckey = {(r['family_id'], r['copy_idx']) for r in cop}
qkey = collections.Counter((r['family_id'], r['copy_index']) for r in q)
chk('quant: copy exists in the copy table', [k for k in qkey if k not in ckey], len(qkey))
chk('quant: one row per copy', [k for k in qkey if qkey[k] > 1], len(qkey))
chk('quant: no excluded family', sorted({r['family_id'] for r in q} & skip), len(qkey))
a_fam = {r['family_id'] for r in asg}
chk('assignments: family was assigned (in quant)', sorted(a_fam - {r['family_id'] for r in q}), len(a_fam))
rf = collections.Counter((r['read_name'], r['family_id']) for r in asg)
chk('assignments: one row per read x family', [k for k in rf if rf[k] > 1], len(rf))
ass = [r for r in asg if r['status'] == 'assigned']
chk('assigned rows: copy exists', [(r['read_name'], r['family_id'], r['assigned_copy']) for r in ass if (r['family_id'], r['assigned_copy']) not in ckey], len(ass))
chk('assigned rows: AS-tied (as_margin 0), at least one candidate, never an orphan', [r['read_name'] for r in ass if r['as_margin'] != '0' or int(r['n_candidates']) < 1 or r['origin_rejected'] == '1'], len(ass))
chk('sole_candidate = assigned with one candidate and contested', [r['read_name'] for r in asg if (r['sole_candidate'] == '1') != (r['status'] == 'assigned' and r['n_candidates'] == '1' and r['contested'] == '1')], len(asg))
# the hard count is by best copy over every status, minus reads with no candidate at all (copy_assign_bin.rs, quant write)
hard = collections.Counter((r['family_id'], r['assigned_copy']) for r in asg if not (r['origin_rejected'] == '1' and r['n_candidates'] == '0'))
qh = {(r['family_id'], r['copy_index']): int(r['n_reads_hard']) for r in q}
chk('quant: n_reads_hard = rows whose best copy it is, orphans out', [k for k in qh if qh[k] != hard.get(k, 0)], len(qh))
nf = collections.Counter(r['family_id'] for r in asg)
chk('families: n_reads = assignment rows of the family', [r['family_id'] for r in fam if int(r['n_reads']) != nf.get(r['family_id'], 0)], len(fam))
multi = collections.Counter(r['read_name'] for r in ass)
print(f"info assigned molecules {len(multi)}, assigned in more than one family {sum(1 for v in multi.values() if v > 1)}; statuses {dict(collections.Counter(r['status'] for r in asg))}")
