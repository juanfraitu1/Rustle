#!/usr/bin/env python3
"""Full evaluation: compute copy-exclusive node sets, assign copy from GAF
(graph) and from minimap2 SAM (linear), compare against ground truth derived
from read names."""
import re, sys, collections

# ---------- ground truth from read name ----------
def truth_of(read):
    if read.startswith("uniq_A"):   return "A"
    if read.startswith("uniq_B"):   return "B"
    if read.startswith("multi_r_"): return "B"   # multi_r = TRUE copy B
    if read.startswith("multi_"):   return "A"   # multi   = TRUE copy A
    return "?"

def group_of(read):
    if read.startswith("uniq_A1"): return "uniq_A1"
    if read.startswith("uniq_A2"): return "uniq_A2"
    if read.startswith("uniq_B1"): return "uniq_B1"
    if read.startswith("multi_r_"): return "multi_r"
    if read.startswith("multi_"):  return "multi"
    return "?"

# ---------- parse GFA for exclusive node sets ----------
paths = {}
node_seq_len = {}
with open("graph.gfa") as f:
    for line in f:
        if line.startswith("S"):
            t = line.rstrip("\n").split("\t")
            node_seq_len[t[1]] = len(t[2])
        elif line.startswith("P"):
            t = line.rstrip("\n").split("\t")
            paths[t[1]] = [x[:-1] for x in t[2].split(",")]

A_nodes = set(paths["A1"]) | set(paths["A2"])
B_nodes = set(paths["B1"])
A_excl = A_nodes - B_nodes
B_excl = B_nodes - A_nodes
shared = A_nodes & B_nodes

graph_stats = dict(
    nodes=len(node_seq_len),
    A_path=len(A_nodes), B_path=len(B_nodes), shared=len(shared),
    A_excl=len(A_excl), B_excl=len(B_excl),
    A_excl_bp=sum(node_seq_len[n] for n in A_excl),
    B_excl_bp=sum(node_seq_len[n] for n in B_excl),
)

# ---------- GRAPH assignment from GAF ----------
node_re = re.compile(r'[<>](\d+)')
def parse_gaf(path):
    best = {}
    with open(path) as f:
        for line in f:
            t = line.rstrip("\n").split("\t")
            if len(t) < 12: continue
            read = t[0]
            try: matches = int(t[9])
            except ValueError: matches = 0
            if read not in best or matches > best[read][0]:
                best[read] = (matches, t[5])
    calls = {}
    detail = {}
    for read, (m, pathstr) in best.items():
        nodes = node_re.findall(pathstr)
        a = sum(1 for n in nodes if n in A_excl)
        b = sum(1 for n in nodes if n in B_excl)
        call = "A" if a > b else "B" if b > a else "ambiguous"
        calls[read] = call
        detail[read] = (a, b)
    return calls, detail

# ---------- LINEAR assignment from minimap2 SAM vs A1+B1 ----------
def parse_sam_linear(path):
    # best AS per (read, target-copy). copy = first char-group of ref name A/B
    best_score = collections.defaultdict(dict)  # read -> {copy: AS}
    with open(path) as f:
        for line in f:
            if line.startswith("@"): continue
            t = line.rstrip("\n").split("\t")
            read, flag, ref = t[0], int(t[1]), t[2]
            if flag & 0x4: continue          # unmapped
            if ref == "*": continue
            copy = "A" if ref.startswith("A") else "B" if ref.startswith("B") else "?"
            AS = None
            for f_ in t[11:]:
                if f_.startswith("AS:i:"): AS = int(f_[5:])
            if AS is None: continue
            d = best_score[read]
            if copy not in d or AS > d[copy]:
                d[copy] = AS
    calls = {}
    detail = {}
    for read, d in best_score.items():
        a = d.get("A", None); b = d.get("B", None)
        if a is None and b is None: continue
        if a is None: call = "B"
        elif b is None: call = "A"
        elif a > b: call = "A"
        elif b > a: call = "B"
        else: call = "ambiguous"
        calls[read] = call
        detail[read] = (a, b)
    return calls, detail

graph_calls, graph_detail = parse_gaf("aln.gaf")
linear_calls, linear_detail = parse_sam_linear("linear.sam")

all_reads = sorted(set(graph_calls) | set(linear_calls))

# ---------- accuracy on uniq_* (known, unambiguous truth) ----------
def accuracy(calls, read_filter):
    n=0; correct=0; amb=0
    for r, c in calls.items():
        if not read_filter(r): continue
        n += 1
        if c == "ambiguous": amb += 1; continue
        if c == truth_of(r): correct += 1
    return correct, n, amb

is_uniq = lambda r: r.startswith("uniq_")
ga_c, ga_n, ga_amb = accuracy(graph_calls, is_uniq)
li_c, li_n, li_amb = accuracy(linear_calls, is_uniq)

print("="*64)
print("GRAPH STATS")
for k,v in graph_stats.items(): print(f"  {k:12s}: {v}")
print()
print("="*64)
print("ACCURACY ON uniq_* READS (n correct / n total ; ambiguous)")
print(f"  {'method':14s} {'correct':>8s} {'total':>6s} {'amb':>5s} {'acc':>7s}")
print(f"  {'graph (GAF)':14s} {ga_c:>8d} {ga_n:>6d} {ga_amb:>5d} {ga_c/ga_n*100:>6.1f}%")
print(f"  {'linear (mm2)':14s} {li_c:>8d} {li_n:>6d} {li_amb:>5d} {li_c/li_n*100:>6.1f}%")
print()

# also break out by group A vs B for uniq
for label, filt in [("uniq_A* (truth A)", lambda r: r.startswith("uniq_A")),
                    ("uniq_B* (truth B)", lambda r: r.startswith("uniq_B"))]:
    gc,gn,gamb = accuracy(graph_calls, filt)
    lc,ln,lamb = accuracy(linear_calls, filt)
    print(f"  {label:18s}  graph {gc}/{gn} (amb {gamb})   linear {lc}/{ln} (amb {lamb})")
print()

# ---------- multi reads: truth IS known from names ----------
print("="*64)
print("MULTI-MAPPER SPLIT (truth derivable: multi=A, multi_r=B)")
for label, filt in [("multi   (truth A)", lambda r: r.startswith("multi_") and not r.startswith("multi_r")),
                    ("multi_r (truth B)", lambda r: r.startswith("multi_r"))]:
    g = collections.Counter(graph_calls[r] for r in graph_calls if filt(r))
    l = collections.Counter(linear_calls[r] for r in linear_calls if filt(r))
    gc,gn,gamb = accuracy(graph_calls, filt)
    lc,ln,lamb = accuracy(linear_calls, filt)
    print(f"  {label}")
    print(f"     graph : A={g['A']} B={g['B']} amb={g['ambiguous']}   correct {gc}/{gn}")
    print(f"     linear: A={l['A']} B={l['B']} amb={l['ambiguous']}   correct {lc}/{ln}")
print()

# overall accuracy including multi (truth from names)
all_filt = lambda r: True
gc,gn,gamb = accuracy(graph_calls, all_filt)
lc,ln,lamb = accuracy(linear_calls, all_filt)
print("="*64)
print("OVERALL ACCURACY (all 103 reads, truth from names)")
print(f"  graph : {gc}/{gn}  ({gc/gn*100:.1f}%)  ambiguous={gamb}")
print(f"  linear: {lc}/{ln}  ({lc/ln*100:.1f}%)  ambiguous={lamb}")
print()

# ---------- per-read disagreement / error dump ----------
print("="*64)
print("ERRORS / DISAGREEMENTS (read  truth  graph(a,b)  linear(a,b))")
for r in all_reads:
    g = graph_calls.get(r,"-"); l = linear_calls.get(r,"-")
    tr = truth_of(r)
    gerr = (g!=tr and g!="ambiguous"); lerr = (l!=tr and l!="ambiguous")
    if gerr or lerr or g=="ambiguous" or l=="ambiguous":
        print(f"  {r:14s} truth={tr}  graph={g}{graph_detail.get(r,'')}  linear={l}{linear_detail.get(r,'')}")
print("  (none beyond those listed)" )
