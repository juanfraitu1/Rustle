#!/usr/bin/env python3
"""Evaluate paralog copy assignment for one structural regime.

Usage: python3 evaluate.py <regime_dir>

Methods compared against ground truth (from read name uniqA_/uniqB_):
  GRAPH        : GraphAligner GAF; assign by exclusive-node vote + (inversion)
                 orientation-through-flipped-segment vote.
  LINEAR-AS    : minimap2 to copyA and copyB separately; assign by best AS.
  LINEAR-NAIVE : minimap2 to concat(copyA,copyB); assign by where the PRIMARY
                 alignment lands (what a coordinate-ingesting tool sees).

Reports per-class accuracy AND decision margins (AS gaps, node-vote gaps),
plus split/soft-clip diagnostics for the naive primary.
"""
import sys, os, re, collections

D = sys.argv[1]
REGIME = os.path.basename(os.path.abspath(D))


def truth_of(rd):
    if rd.startswith("uniqA_"):
        return "A"
    if rd.startswith("uniqB_"):
        return "B"
    return "?"


def class_of(rd):
    # uniqA_<class>_<idx>_...
    m = re.match(r"uniq[AB]_([a-z_]+?)_\d", rd)
    return m.group(1) if m else "?"


# ---------------------------------------------------------------- meta
meta = {}
for line in open(os.path.join(D, "meta.txt")):
    k, v = line.strip().split("=")
    meta[k] = v

# ---------------------------------------------------------------- GFA
paths = {}
nodelen = {}
for line in open(os.path.join(D, "graph.gfa")):
    t = line.rstrip("\n").split("\t")
    if t[0] == "S":
        nodelen[t[1]] = len(t[2])
    elif t[0] == "P":
        paths[t[1]] = [(x[:-1], x[-1]) for x in t[2].split(",")]

A = paths["copyA"]
B = paths["copyB"]
Anodes = {n for n, o in A}
Bnodes = {n for n, o in B}
Aor = {n: o for n, o in A}
Bor = {n: o for n, o in B}
shared = Anodes & Bnodes
A_excl = Anodes - Bnodes
B_excl = Bnodes - Anodes
# inversion signature: nodes shared but with FLIPPED orientation between paths
flip_nodes = {n for n in shared if Aor[n] != Bor[n]}
# For flipped nodes: copyA's orientation is the "A vote" orientation.
flipA_or = {n: Aor[n] for n in flip_nodes}  # orientation that means copy A

# ---------------------------------------------------------------- GRAPH (GAF)
# GAF col6 path string e.g. >12>13<30<31...  We extract (node, orient).
gaf_token = re.compile(r'([<>])(\d+)')
graph_assign = {}
graph_margin = {}
graph_detail = {}
for line in open(os.path.join(D, "aln.gaf")):
    t = line.rstrip("\n").split("\t")
    if len(t) < 12:
        continue
    rd = t[0]
    pathstr = t[5]
    toks = gaf_token.findall(pathstr)  # list of ('>'/'<', node)
    # node-orientation as traversed by the read
    read_or = {}  # node -> '+' if '>', '-' if '<'
    visited = []
    for sym, node in toks:
        read_or[node] = '+' if sym == '>' else '-'
        visited.append(node)
    vset = set(visited)
    # vote 1: exclusive nodes
    a_excl_hit = len(vset & A_excl)
    b_excl_hit = len(vset & B_excl)
    # vote 2: orientation through flipped (inversion) nodes
    a_or_hit = 0
    b_or_hit = 0
    for n in vset & flip_nodes:
        # copyA traverses node n in orientation flipA_or[n].
        # The read may traverse the whole path in reverse; to be orientation
        # invariant we compare the RELATIVE orientation pattern. Simpler &
        # robust: count how the read's orientation at flipped nodes matches A's
        # vs B's, but normalized by the read's global strand. We approximate the
        # read's global strand by majority orientation on NON-flipped shared
        # nodes (which are same-orient in both copies).
        pass
    # determine read global strand from non-flipped shared nodes
    nonflip_shared = (vset & shared) - flip_nodes
    fwd_votes = sum(1 for n in nonflip_shared if read_or[n] == Aor[n])
    rev_votes = sum(1 for n in nonflip_shared if read_or[n] != Aor[n])
    read_strand_fwd = fwd_votes >= rev_votes  # True if read traverses graph in copyA's forward sense
    for n in vset & flip_nodes:
        # copyA orientation at n (in graph-forward sense) = flipA_or[n]
        # copyB orientation at n = the flip of that.
        ro = read_or[n]
        if not read_strand_fwd:
            ro = '+' if ro == '-' else '-'  # normalize read to graph-forward
        if ro == flipA_or[n]:
            a_or_hit += 1
        else:
            b_or_hit += 1
    a_total = a_excl_hit + a_or_hit
    b_total = b_excl_hit + b_or_hit
    if a_total > b_total:
        graph_assign[rd] = "A"
    elif b_total > a_total:
        graph_assign[rd] = "B"
    else:
        graph_assign[rd] = "?"
    graph_margin[rd] = abs(a_total - b_total)
    graph_detail[rd] = (a_excl_hit, b_excl_hit, a_or_hit, b_or_hit)


# ---------------------------------------------------------------- LINEAR-AS
def parse_sam_as(path):
    """read -> best AS among PRIMARY/secondary records (max)."""
    best = {}
    for line in open(path):
        if line.startswith("@"):
            continue
        t = line.rstrip("\n").split("\t")
        rd = t[0]
        flag = int(t[1])
        if flag & 4:
            continue  # unmapped
        as_v = None
        for f in t[11:]:
            if f.startswith("AS:i:"):
                as_v = int(f[5:])
        if as_v is None:
            continue
        if rd not in best or as_v > best[rd]:
            best[rd] = as_v
    return best


asA = parse_sam_as(os.path.join(D, "linA.sam"))
asB = parse_sam_as(os.path.join(D, "linB.sam"))
all_reads = set()
for line in open(os.path.join(D, "reads.fa")):
    if line.startswith(">"):
        all_reads.add(line[1:].strip())

linear_as_assign = {}
linear_as_margin = {}
for rd in all_reads:
    a = asA.get(rd, None)
    b = asB.get(rd, None)
    if a is None and b is None:
        linear_as_assign[rd] = "?"
        linear_as_margin[rd] = 0
    elif b is None:
        linear_as_assign[rd] = "A"
        linear_as_margin[rd] = a  # only A hit
    elif a is None:
        linear_as_assign[rd] = "B"
        linear_as_margin[rd] = b
    else:
        if a > b:
            linear_as_assign[rd] = "A"
        elif b > a:
            linear_as_assign[rd] = "B"
        else:
            linear_as_assign[rd] = "?"
        linear_as_margin[rd] = abs(a - b)


# ---------------------------------------------------------------- LINEAR-NAIVE
# concat reference; assign by where the PRIMARY (not secondary/suppl) lands.
# Also collect split diagnostics: does the read get supplementary/secondary or
# heavy soft-clip on its primary?
naive_assign = {}
naive_primary_clip = {}      # read -> soft-clip bp on primary
naive_has_suppl = collections.defaultdict(bool)
naive_has_secondary = collections.defaultdict(bool)
naive_primary_as = {}


def softclip_bp(cigar):
    return sum(int(n) for n, op in re.findall(r'(\d+)([MIDNSHP=X])', cigar) if op == 'S')


for line in open(os.path.join(D, "naive.sam")):
    if line.startswith("@"):
        continue
    t = line.rstrip("\n").split("\t")
    rd = t[0]
    flag = int(t[1])
    if flag & 4:
        naive_assign.setdefault(rd, "?")
        continue
    rname = t[2]  # copyA or copyB
    cig = t[5]
    if flag & 2048:
        naive_has_suppl[rd] = True
        continue
    if flag & 256:
        naive_has_secondary[rd] = True
        continue
    # PRIMARY
    naive_assign[rd] = "A" if rname == "copyA" else "B"
    naive_primary_clip[rd] = softclip_bp(cig)
    for f in t[11:]:
        if f.startswith("AS:i:"):
            naive_primary_as[rd] = int(f[5:])

for rd in all_reads:
    naive_assign.setdefault(rd, "?")

# ---------------------------------------------------------------- score
classes = collections.OrderedDict()
for rd in sorted(all_reads):
    c = class_of(rd)
    classes.setdefault(c, []).append(rd)

methods = [
    ("GRAPH", graph_assign, graph_margin),
    ("LINEAR-AS", linear_as_assign, linear_as_margin),
    ("LINEAR-NAIVE-PRIMARY", naive_assign, None),
]


def acc_table(read_subset):
    rows = {}
    for name, assign, margin in methods:
        cor = wrong = amb = 0
        margins = []
        for rd in read_subset:
            tr = truth_of(rd)
            pr = assign.get(rd, "?")
            if pr == "?":
                amb += 1
            elif pr == tr:
                cor += 1
                if margin:
                    margins.append(margin[rd])
            else:
                wrong += 1
        rows[name] = (cor, wrong, amb, len(read_subset), margins)
    return rows


print("=" * 78)
print("REGIME: %s    (%d reads total)" % (REGIME, len(all_reads)))
print("=" * 78)
print("\nGraph structure:")
print("  nodes=%d  A-excl=%d (%dbp)  B-excl=%d (%dbp)  flipped/inversion nodes=%d (%dbp)"
      % (len(nodelen), len(A_excl), sum(nodelen[n] for n in A_excl),
         len(B_excl), sum(nodelen[n] for n in B_excl),
         len(flip_nodes), sum(nodelen[n] for n in flip_nodes)))

# Decide which classes are "resolvable" (exclude nospan_nosnp from headline)
resolvable = [rd for rd in all_reads if class_of(rd) != "nospan_nosnp"]
print("\n--- ACCURACY by read class ---")
hdr = "%-14s | %-22s | %-22s | %-22s" % ("class (n)", "GRAPH", "LINEAR-AS", "LINEAR-NAIVE-PRIM")
print(hdr)
print("-" * len(hdr))
for c, reads in classes.items():
    rows = acc_table(reads)
    cells = []
    for name, _, _ in methods:
        cor, wr, amb, tot, _ = rows[name]
        cells.append("%d/%d cor (%dwr %damb)" % (cor, tot, wr, amb))
    tag = c + (" [UNRESOLVABLE]" if c == "nospan_nosnp" else "")
    print("%-14s | %-22s | %-22s | %-22s" % ("%s (%d)" % (tag, len(reads)),
                                             cells[0], cells[1], cells[2]))

print("\n--- HEADLINE (resolvable reads only, n=%d) ---" % len(resolvable))
rows = acc_table(resolvable)
for name, _, _ in methods:
    cor, wr, amb, tot, margins = rows[name]
    mtxt = ""
    if margins:
        mtxt = "  margin min=%g mean=%.1f max=%g" % (min(margins), sum(margins) / len(margins), max(margins))
    print("  %-22s %d/%d correct  %d wrong  %d ambiguous%s"
          % (name, cor, tot, wr, amb, mtxt))

# ---------------------------------------------------------------- split diagnostics
print("\n--- NAIVE split / soft-clip diagnostics (inversion focus) ---")
span_reads = [rd for rd in all_reads if class_of(rd) == "span"]
n_suppl = sum(1 for rd in span_reads if naive_has_suppl[rd])
n_sec = sum(1 for rd in span_reads if naive_has_secondary[rd])
clips = [naive_primary_clip.get(rd, 0) for rd in span_reads]
heavy_clip = sum(1 for c in clips if c >= 50)
print("  span reads n=%d:  with supplementary=%d  with secondary=%d  primary softclip>=50bp=%d"
      % (len(span_reads), n_suppl, n_sec, heavy_clip))
if clips:
    print("  span primary soft-clip bp: min=%d mean=%.1f max=%d"
          % (min(clips), sum(clips) / len(clips), max(clips)))
# Did naive primary land on the WRONG copy for any span read?
wrong_span = [(rd, naive_assign.get(rd)) for rd in span_reads if naive_assign.get(rd) != truth_of(rd)]
print("  span reads where naive PRIMARY mis-assigned: %d" % len(wrong_span))
for rd, pr in wrong_span[:10]:
    print("     %s  truth=%s  naive=%s  softclip=%dbp"
          % (rd, truth_of(rd), pr, naive_primary_clip.get(rd, -1)))

# Per-class margin detail for graph (node votes) and linear (AS)
print("\n--- DECISION MARGINS by class ---")
print("%-16s | %-26s | %-26s" % ("class", "GRAPH node-vote gap", "LINEAR-AS gap"))
for c, reads in classes.items():
    gm = [graph_margin[rd] for rd in reads if graph_assign.get(rd) == truth_of(rd)]
    lm = [linear_as_margin[rd] for rd in reads if linear_as_assign.get(rd) == truth_of(rd)]
    def stat(xs):
        return "min=%g mean=%.1f max=%g" % (min(xs), sum(xs) / len(xs), max(xs)) if xs else "n/a"
    print("%-16s | %-26s | %-26s" % (c, stat(gm), stat(lm)))
