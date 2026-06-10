#!/usr/bin/env python3
"""
§6L — Default rustle-only (187) first-divergence stage distribution.

Classifies each of the 187 rustle-only multi-intron chains by the EARLIEST
pipeline stage where rustle has it and StringTie does not.

Stages (earliest → latest):
  junction    — ST didn't accept ≥1 of the chain's junctions (rustle did)
  graph       — all junctions accepted by both, but ST didn't build graph nodes
                spanning the chain (rustle did) [approximate: bundlenode_list vs graphnode_list]
  collapse    — graphed by both, but not an ST transfrag_collapse rep (rustle keeps it)
  seed        — ST has it as collapse rep but never seeds it; rustle seeds it
  flow        — ST seeds it but never extracts it (path_extracted); rustle extracts it
  post_flow   — ST extracts it but not in ST's final GTF; rustle keeps it
  unknown     — couldn't classify

Usage: python3 bench/stage_divergence_187.py [ru.gtf] [st.gtf] [ru_ev.jsonl] [st_ev.jsonl]
"""

import re, sys, json, collections, itertools

RU_GTF   = sys.argv[1] if len(sys.argv) > 1 else "/tmp/ru.gtf"
ST_GTF   = sys.argv[2] if len(sys.argv) > 2 else "/tmp/st.gtf"
RU_EV    = sys.argv[3] if len(sys.argv) > 3 else "/tmp/ru_ev.jsonl"
ST_EV    = sys.argv[4] if len(sys.argv) > 4 else "/tmp/st_ev.jsonl"

# ── helpers ────────────────────────────────────────────────────────────────────

def load_gtf_chains(path):
    """Return set of (strand, intron-chain-tuple) for all multi-intron transcripts."""
    tx = collections.defaultdict(list)
    strand = {}
    for line in open(path):
        if line.startswith('#'): continue
        f = line.rstrip('\n').split('\t')
        if len(f) < 9 or f[2] != 'exon': continue
        m = re.search(r'transcript_id "([^"]+)"', f[8])
        if not m: continue
        tx[m.group(1)].append((int(f[3]), int(f[4])))
        strand[m.group(1)] = f[6]
    chains = set()
    for t, ex in tx.items():
        ex.sort()
        if len(ex) < 2: continue
        introns = tuple((ex[i][1], ex[i+1][0]) for i in range(len(ex)-1))
        chains.add((strand[t], introns))
    return chains

def chain_junctions(introns):
    """Convert intron-chain tuple to set of (donor, acceptor) junction keys.
    introns: tuple of (exon_end, next_exon_start) in 1-based inclusive coords.
    junction_accept convention: start = exon_end (donor), end = next_exon_start (acceptor).
    These are the SAME values stored in the GTF intron tuple — no offset needed.
    """
    return frozenset(introns)

def parse_intron_str(s):
    """Parse 'A-B,C-D,...' intron string (1-based intron coords from path_extracted/transfrag_collapse)
    into sorted tuple of (exon_end, next_exon_start) pairs — matching GTF chain format.
    path_extracted intron A-B: A=intron_start=exon_end+1, B=intron_end=next_exon_start-1
    → GTF chain coord: (exon_end, next_exon_start) = (A-1, B+1)
    junction_accept also uses (exon_end, next_exon_start).
    """
    if not s or s == "SE" or s == "":
        return ()
    parts = []
    for seg in s.split(','):
        seg = seg.strip()
        if not seg: continue
        if '-' not in seg: continue
        a, b = seg.split('-')
        # Convert from (intron_start, intron_end) to (exon_end, next_exon_start)
        parts.append((int(a) - 1, int(b) + 1))
    return tuple(sorted(parts))

_STEP_RE      = re.compile(r'"step":"([^"]+)"')
_TOOL_RE      = re.compile(r'"tool":"([^"]+)"')
_CHROM_RE     = re.compile(r'"chrom":"([^"]+)"')
_START_RE     = re.compile(r'"start":(\d+)')
_END_RE       = re.compile(r'"end":(\d+)')
_STRAND_RE    = re.compile(r'"strand":"([^"]+)"')
_REP_INT_RE   = re.compile(r'"rep_introns":"([^"]*)"')
_INTRONS_RE   = re.compile(r'"introns":"([^"]*)"')
_ACCEPTED_RE  = re.compile(r'"accepted":(true|false)')
_NODES_RE     = re.compile(r'"nodes":"([^"]*)"')
_N_NODES_RE   = re.compile(r'"n_nodes":(\d+)')
_NACCEPTED_RE = re.compile(r'"n_introns":(\d+)')

def _fallback_parse(line):
    """Regex-based fallback parser for truncated JSON lines."""
    d = {}
    m = _STEP_RE.search(line); d['step'] = m.group(1) if m else ''
    m = _TOOL_RE.search(line); d['tool'] = m.group(1) if m else ''
    m = _CHROM_RE.search(line); d['chrom'] = m.group(1) if m else ''
    m = _START_RE.search(line); d['start'] = int(m.group(1)) if m else 0
    m = _END_RE.search(line); d['end'] = int(m.group(1)) if m else 0
    m = _STRAND_RE.search(line)
    # "strand" appears twice: once at top level, once maybe in payload
    if m: d['strand'] = m.group(1)
    else: d['strand'] = '.'
    payload = {}
    m = _REP_INT_RE.search(line)
    if m: payload['rep_introns'] = m.group(1)
    m = _INTRONS_RE.search(line)
    if m: payload['introns'] = m.group(1)
    m = _ACCEPTED_RE.search(line)
    if m: payload['accepted'] = (m.group(1) == 'true')
    m = _NODES_RE.search(line)
    if m: payload['nodes'] = m.group(1)
    m = _N_NODES_RE.search(line)
    if m: payload['n_nodes'] = int(m.group(1))
    m = _NACCEPTED_RE.search(line)
    if m: payload['n_introns'] = int(m.group(1))
    d['payload'] = payload
    return d

def load_events(path):
    """Parse JSONL events, return dict step -> list of event dicts."""
    events = collections.defaultdict(list)
    for line in open(path):
        line = line.strip()
        if not line: continue
        try:
            d = json.loads(line)
        except json.JSONDecodeError:
            # Truncated line: use regex fallback
            d = _fallback_parse(line)
        step = d.get('step', '')
        if not step or step.startswith('_'): continue
        events[step].append(d)
    return events

# ── load data ──────────────────────────────────────────────────────────────────

print("Loading GTFs...", flush=True)
ru_chains = load_gtf_chains(RU_GTF)
st_chains = load_gtf_chains(ST_GTF)
ru_only = ru_chains - st_chains
print(f"  Rustle-only chains: {len(ru_only)}")

print("Loading parity events...", flush=True)
ru_ev = load_events(RU_EV)
st_ev = load_events(ST_EV)

# Print event counts per step
for step in ['junction_accept','graphnode_list','transfrag_collapse','parse_trflong_seed','path_extracted']:
    print(f"  {step}: rustle={len(ru_ev[step])} st={len(st_ev[step])}")

# ── build lookup sets for each stage ──────────────────────────────────────────

# Stage 1: junction_accept — accepted junctions per tool
# Key: (start, end, strand) — strand may be '.' for rustle events (no chrom filter used here)
# rustle junction_accept has no chrom; we match by (start,end,strand)
def build_junc_accepted(events):
    """Returns set of (start, end, strand) for accepted junctions."""
    s = set()
    for ev in events:
        p = ev.get('payload', {})
        if p.get('accepted', False):
            strand = ev.get('strand', '.')
            s.add((ev['start'], ev['end'], strand))
    return s

ru_junc_accepted = build_junc_accepted(ru_ev['junction_accept'])
st_junc_accepted = build_junc_accepted(st_ev['junction_accept'])

print(f"\n  ru_junc_accepted: {len(ru_junc_accepted)}, st_junc_accepted: {len(st_junc_accepted)}")

# Stage 2: graphnode_list — graph nodes per tool
# We collect all graphnode intervals per (chrom, strand) bundle so we can check
# if all exon segments of a chain fall within graph-node coverage.
# graphnode_list payload: {"n_nodes":N,"nodes":"start-end:cov,..."}
# We'll collect all node intervals per bundle (bundle keyed by chrom+strand).
# A chain is "graphed" if ALL its exon intervals are covered by some node in any
# single bundle that spans the chain region.
# Approximation: we check if each exon has any overlapping node in the same (chrom,strand) pool.
# Note: rustle graphnode_list has chrom; ST's bundlenode_list is different from graphnode_list.
# Both tools emit "graphnode_list" for the flow-graph nodes.

def build_graph_node_intervals(events):
    """Build dict (chrom, strand) -> sorted list of (node_start, node_end)."""
    ivals = collections.defaultdict(list)
    for ev in events:
        chrom = ev.get('chrom', '')
        strand = ev.get('strand', '.')
        p = ev.get('payload', {})
        nodes_str = p.get('nodes', '')
        for seg in nodes_str.split(','):
            seg = seg.strip()
            if not seg: continue
            m = re.match(r'^(\d+)-(\d+)', seg)
            if m:
                ivals[(chrom, strand)].append((int(m.group(1)), int(m.group(2))))
    return ivals

print("Building graph node intervals...", flush=True)
ru_graph_ivals = build_graph_node_intervals(ru_ev['graphnode_list'])
st_graph_ivals = build_graph_node_intervals(st_ev['graphnode_list'])

def exon_in_nodes(exon_start, exon_end, node_list):
    """Return True if exon [exon_start, exon_end] overlaps any node interval."""
    for ns, ne in node_list:
        if ns <= exon_end and ne >= exon_start:
            return True
    return False

# Stage 3: transfrag_collapse — rep_introns per tool
# Key: (strand, intron-chain-tuple)
def build_collapse_chains(events):
    s = set()
    for ev in events:
        p = ev.get('payload', {})
        ri = p.get('rep_introns', '')
        if ri == 'SE' or ri == '': continue
        # Handle truncated rep_introns (ST has a 480-char line buffer):
        # drop the last segment if it doesn't look like "N-N"
        parts = ri.split(',')
        if parts and '-' not in parts[-1]:
            parts = parts[:-1]
        if not parts: continue
        valid = True
        introns_list = []
        for seg in parts:
            seg = seg.strip()
            if not seg: continue
            m = re.match(r'^(\d+)-(\d+)$', seg)
            if not m: valid = False; break
            introns_list.append((int(m.group(1)), int(m.group(2))))
        if not valid or not introns_list: continue
        introns = tuple(sorted(introns_list))
        strand = ev.get('strand', '.')
        s.add((strand, introns))
    return s

ru_collapse = build_collapse_chains(ru_ev['transfrag_collapse'])
st_collapse = build_collapse_chains(st_ev['transfrag_collapse'])
print(f"  collapse chains: ru={len(ru_collapse)}, st={len(st_collapse)}")

# Stage 4: parse_trflong_seed — seeded chains per tool
def build_seed_chains(events):
    s = set()
    for ev in events:
        p = ev.get('payload', {})
        introns_str = p.get('introns', '')
        if not introns_str or p.get('n_introns', 0) == 0: continue
        introns = parse_intron_str(introns_str)
        if not introns: continue
        strand = ev.get('strand', '.')
        s.add((strand, introns))
    return s

ru_seeds = build_seed_chains(ru_ev['parse_trflong_seed'])
st_seeds = build_seed_chains(st_ev['parse_trflong_seed'])
print(f"  seed chains: ru={len(ru_seeds)}, st={len(st_seeds)}")

# Stage 5: path_extracted — extracted chains per tool
def build_extracted_chains(events):
    s = set()
    for ev in events:
        p = ev.get('payload', {})
        introns_str = p.get('introns', '')
        if not introns_str: continue
        introns = parse_intron_str(introns_str)
        if not introns: continue
        strand = ev.get('strand', '.')
        s.add((strand, introns))
    return s

ru_extracted = build_extracted_chains(ru_ev['path_extracted'])
st_extracted = build_extracted_chains(st_ev['path_extracted'])
print(f"  extracted chains: ru={len(ru_extracted)}, st={len(st_extracted)}")

# ── need chrom info for rustle-only chains ─────────────────────────────────────
# Load chrom + strand info from rustle's path_extracted events for each chain.
# Also need to know which chrom each chain lives on for graph-node lookup.

def build_chain_to_chrom(events):
    """Map (strand, introns-tuple) -> chrom from path_extracted events."""
    m = {}
    for ev in events:
        p = ev.get('payload', {})
        introns_str = p.get('introns', '')
        if not introns_str: continue
        introns = parse_intron_str(introns_str)
        if not introns: continue
        strand = ev.get('strand', '.')
        chrom = ev.get('chrom', '')
        m[(strand, introns)] = chrom
    return m

# Also load from GTF for chains not in path_extracted (shouldn't happen for final GTF chains)
def build_chain_to_chrom_from_gtf(path):
    tx = collections.defaultdict(list)
    strand = {}
    chrom_map = {}
    for line in open(path):
        if line.startswith('#'): continue
        f = line.rstrip('\n').split('\t')
        if len(f) < 9 or f[2] != 'exon': continue
        m = re.search(r'transcript_id "([^"]+)"', f[8])
        if not m: continue
        tid = m.group(1)
        tx[tid].append((int(f[3]), int(f[4])))
        strand[tid] = f[6]
        chrom_map[tid] = f[0]
    result = {}
    for t, ex in tx.items():
        ex.sort()
        if len(ex) < 2: continue
        introns = tuple((ex[i][1], ex[i+1][0]) for i in range(len(ex)-1))
        result[(strand[t], introns)] = chrom_map[t]
    return result

print("Building chain->chrom maps...", flush=True)
chain_to_chrom_ru_ev = build_chain_to_chrom(ru_ev['path_extracted'])
chain_to_chrom_gtf   = build_chain_to_chrom_from_gtf(RU_GTF)

def get_chrom(strand_chain_key):
    return chain_to_chrom_ru_ev.get(strand_chain_key) or chain_to_chrom_gtf.get(strand_chain_key, '')

# ── stage 2: check if chain's exons are covered by graph nodes ─────────────────

def chain_graphed_by(strand_chain, graph_ivals_by_cs):
    """True if all exons of chain are covered by graph nodes for (chrom, strand)."""
    strand, introns = strand_chain
    chrom = get_chrom(strand_chain)
    if not chrom: return False
    key = (chrom, strand)
    nodes = graph_ivals_by_cs.get(key, [])
    if not nodes: return False
    # Reconstruct exon intervals from introns
    # introns: [(exon1_end, exon2_start), ...]  (1-based inclusive exon coords)
    # Exon spans: exon1=[introns[0][0] ... first junction], etc.
    # We check each intron endpoint (which bounds an exon boundary)
    # For simplicity: check if each intron-boundary coordinate has a covering node.
    # A chain with introns [(a1,b1),(a2,b2),...] has exon boundaries at a1, b1, a2, b2, ...
    # We'll check that at least one node covers b1-1 (last base of exon 1) and a1+1 etc.
    # Simplified check: ALL intron-boundary coords must lie within some node interval.
    all_coords = set()
    for a, b in introns:
        all_coords.add(a - 1)  # last base of left exon
        all_coords.add(b)      # first base of right exon
    for coord in all_coords:
        if not any(ns <= coord <= ne for ns, ne in nodes):
            return False
    return True

# ── classify each rustle-only chain ───────────────────────────────────────────

stage_counts = collections.Counter()
examples = collections.defaultdict(list)  # stage -> list of chain info

# For strand matching in junction lookup: rustle emits strand correctly
# ST emits strand for junctions. For chains without strand info, we try both.

print("\nClassifying 187 rustle-only chains...", flush=True)

for strand_chain in ru_only:
    strand, introns = strand_chain

    # ── Stage 1: junction ──────────────────────────────────────────────────
    junctions = chain_junctions(introns)
    st_missing_junc = False
    missing_junc_example = None
    for jstart, jend in junctions:
        # Try exact strand match first, then '.' (rustle emits strand but ST may differ)
        found = (
            (jstart, jend, strand) in st_junc_accepted or
            (jstart, jend, '+') in st_junc_accepted or
            (jstart, jend, '-') in st_junc_accepted or
            (jstart, jend, '.') in st_junc_accepted
        )
        if not found:
            st_missing_junc = True
            missing_junc_example = (jstart, jend)
            break

    if st_missing_junc:
        stage_counts['junction'] += 1
        if len(examples['junction']) < 5:
            examples['junction'].append({
                'strand': strand, 'introns': introns,
                'missing_junc': missing_junc_example
            })
        continue

    # ── Stage 2: graph ──────────────────────────────────────────────────────
    ru_graphed = chain_graphed_by(strand_chain, ru_graph_ivals)
    st_graphed = chain_graphed_by(strand_chain, st_graph_ivals)

    if ru_graphed and not st_graphed:
        stage_counts['graph'] += 1
        if len(examples['graph']) < 5:
            examples['graph'].append({'strand': strand, 'introns': introns})
        continue

    # ── Stage 3: collapse ──────────────────────────────────────────────────
    # Check if chain is a transfrag_collapse rep in each tool.
    in_ru_collapse = (strand, introns) in ru_collapse
    in_st_collapse = (strand, introns) in st_collapse

    if in_ru_collapse and not in_st_collapse:
        stage_counts['collapse'] += 1
        if len(examples['collapse']) < 5:
            examples['collapse'].append({'strand': strand, 'introns': introns})
        continue

    # ── Stage 4: seed ──────────────────────────────────────────────────────
    in_ru_seed = (strand, introns) in ru_seeds
    in_st_seed = (strand, introns) in st_seeds

    if in_ru_seed and not in_st_seed:
        stage_counts['seed'] += 1
        if len(examples['seed']) < 5:
            examples['seed'].append({'strand': strand, 'introns': introns})
        continue

    # ── Stage 5: flow ──────────────────────────────────────────────────────
    in_ru_extracted = (strand, introns) in ru_extracted
    in_st_extracted = (strand, introns) in st_extracted

    if in_ru_extracted and not in_st_extracted:
        stage_counts['flow'] += 1
        if len(examples['flow']) < 5:
            examples['flow'].append({'strand': strand, 'introns': introns})
        continue

    # ── Stage 6: post_flow ─────────────────────────────────────────────────
    # ST extracted it (or we couldn't tell), but it's not in ST's final GTF
    if in_st_extracted:
        stage_counts['post_flow'] += 1
        if len(examples['post_flow']) < 5:
            examples['post_flow'].append({'strand': strand, 'introns': introns})
        continue

    # ── Unknown ────────────────────────────────────────────────────────────
    stage_counts['unknown'] += 1
    if len(examples['unknown']) < 5:
        examples['unknown'].append({
            'strand': strand, 'introns': introns,
            'ru_graphed': ru_graphed, 'st_graphed': st_graphed,
            'in_ru_collapse': in_ru_collapse, 'in_st_collapse': in_st_collapse,
            'in_ru_seed': in_ru_seed, 'in_st_seed': in_st_seed,
            'in_ru_extracted': in_ru_extracted, 'in_st_extracted': in_st_extracted,
        })

# ── print results ──────────────────────────────────────────────────────────────

total = len(ru_only)
print("\n=== STAGE DISTRIBUTION ===")
print(f"{'Stage':<12} {'Count':>6} {'%':>7}")
print("-" * 28)
stage_order = ['junction', 'graph', 'collapse', 'seed', 'flow', 'post_flow', 'unknown']
for stage in stage_order:
    n = stage_counts.get(stage, 0)
    pct = 100.0 * n / total if total else 0
    print(f"{stage:<12} {n:>6} {pct:>6.1f}%")
print(f"{'TOTAL':<12} {total:>6}")

print("\n=== EXAMPLES (top stages) ===")
for stage in stage_order:
    exs = examples.get(stage, [])
    if not exs: continue
    n = stage_counts.get(stage, 0)
    if n == 0: continue
    print(f"\n--- {stage} ({n} chains) ---")
    for e in exs[:3]:
        chain_str = ",".join(f"{a}-{b}" for a,b in e['introns'])
        print(f"  strand={e['strand']} introns={chain_str}")
        if 'missing_junc' in e and e['missing_junc']:
            j = e['missing_junc']
            print(f"    first missing junction: ({j[0]},{j[1]})")
        if 'ru_graphed' in e:
            print(f"    ru_graphed={e['ru_graphed']} st_graphed={e['st_graphed']} "
                  f"collapse={e['in_ru_collapse']}/{e['in_st_collapse']} "
                  f"seed={e['in_ru_seed']}/{e['in_st_seed']} "
                  f"extracted={e['in_ru_extracted']}/{e['in_st_extracted']}")
