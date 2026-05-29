#!/usr/bin/env python3
"""Phase-0 diagnostic harness for the PATHPAT flow-parity project.

For each over-enumeration extra in Rustle's output (gffcompare class != '='),
attribute it to: its originating seed (seed_tf) + source (flow/checktrf), and
classify the Rustle-vs-ST divergence by comparing its intron chain to ST's
extracted-path set.

Inputs (all already generated):
  RU_GTF   final Rustle GTF
  TMAP     gffcompare tmap (qry class codes)
  RU_PE    Rustle path_extracted JSONL (seed_tf, source, introns)
  ST_PE    ST path_extracted JSONL (introns)
"""
import json, re, sys
from collections import defaultdict, Counter

RU_GTF = "/tmp/full_default.gtf"
TMAP   = "/tmp/cmp_default.full_default.gtf.tmap"
RU_PE  = "/tmp/ru_pe.jsonl"
ST_PE  = "/tmp/st_pe.jsonl"

def introns_from_exons(exons):
    """exons: sorted [(s,e)] 1-based GTF inclusive -> intron tuples (e_i+1, s_{i+1}-1)."""
    ex = sorted(exons)
    return tuple((ex[i][1] + 1, ex[i + 1][0] - 1) for i in range(len(ex) - 1))

def parse_introns_field(s):
    if not s:
        return tuple()
    out = []
    for tok in s.split(","):
        a, b = tok.split("-")
        out.append((int(a), int(b)))
    return tuple(out)

# --- load Rustle GTF: tid -> (chrom, strand, exons, source) ---
gtf = {}
exons_by_tid = defaultdict(list)
for line in open(RU_GTF):
    f = line.rstrip("\n").split("\t")
    if len(f) < 9:
        continue
    m = re.search(r'transcript_id "([^"]+)"', f[8])
    if not m:
        continue
    tid = m.group(1)
    if f[2] == "transcript":
        src = re.search(r'source "([^"]+)"', f[8])
        gtf[tid] = {"chrom": f[0], "strand": f[6], "source": src.group(1) if src else "?"}
    elif f[2] == "exon":
        exons_by_tid[tid].append((int(f[3]), int(f[4])))
for tid, ex in exons_by_tid.items():
    if tid in gtf:
        gtf[tid]["introns"] = introns_from_exons(ex)
        gtf[tid]["nexons"] = len(ex)

# --- load tmap class codes: qry_id -> class ---
cls = {}
with open(TMAP) as fh:
    next(fh, None)
    for line in fh:
        f = line.rstrip("\n").split("\t")
        if len(f) >= 5:
            cls[f[4]] = f[2]

# --- ST extracted intron-chain set (by strand) ---
st_chains = defaultdict(set)        # strand -> set(intron tuple)
st_chain_list = defaultdict(list)   # strand -> list(intron tuple) for nearest search
for line in open(ST_PE):
    try:
        e = json.loads(line)
    except Exception:
        continue
    if e.get("step") != "path_extracted":
        continue
    p = e.get("payload", e)
    ich = parse_introns_field(p.get("introns", ""))
    if ich:
        strand = e.get("strand", ".")
        st_chains[strand].add(ich)
        st_chain_list[strand].append(ich)

# --- Rustle path_extracted: intron chain -> list of (seed_tf, source, cov) ---
ru_pe = defaultdict(list)
for line in open(RU_PE):
    try:
        e = json.loads(line)
    except Exception:
        continue
    if e.get("step") != "path_extracted":
        continue
    p = e.get("payload", e)
    ich = parse_introns_field(p.get("introns", ""))
    ru_pe[ich].append((p.get("seed_tf"), p.get("source"), p.get("cov")))

def classify(extra_ich, strand):
    """Compare extra chain to ST chains on same strand; return (type, best_shared, st_nintrons)."""
    eset = set(extra_ich)
    best = None
    best_shared = -1
    for st in st_chain_list[strand]:
        shared = len(eset & set(st))
        if shared > best_shared:
            best_shared = shared
            best = st
    if best is None:
        return ("no_st_overlap", 0, 0)
    bset = set(best)
    shared = eset & bset
    extra_only = eset - bset       # introns Rustle has that ST doesn't
    st_only = bset - eset          # introns ST has that Rustle doesn't
    if not extra_only and not st_only:
        return ("exact_in_st(!?)", len(shared), len(best))
    if extra_only and not st_only:
        return ("superset_extra_junction", len(shared), len(best))   # Rustle adds junctions
    if st_only and not extra_only:
        return ("subset_contained", len(shared), len(best))           # Rustle drops junctions
    # both differ: alt-splice (shifted donor/acceptor) vs skip
    return ("alt_splice_diverge", len(shared), len(best))

# --- attribute each non-'=' extra ---
EXTRA_CLASSES = set("jcmknoxiep s".split())
rows = []
for tid, info in gtf.items():
    c = cls.get(tid, "?")
    if c == "=":
        continue
    ich = info.get("introns", tuple())
    if not ich:           # single-exon: separate concern, skip here
        continue
    strand = info["strand"]
    in_st = ich in st_chains[strand]
    pe = ru_pe.get(ich, [])
    seed = pe[0][0] if pe else None
    dtype, shared, st_n = classify(ich, strand)
    rows.append((tid, c, info["source"], info["nexons"], seed, in_st, dtype, shared, st_n))

# --- report ---
print(f"Total non-'=' multi-exon extras: {len(rows)}")
print("\n=== by gffcompare class ===")
for k, v in Counter(r[1] for r in rows).most_common():
    print(f"  {v:4d}  {k}")
print("\n=== by Rustle source ===")
for k, v in Counter(r[2] for r in rows).most_common():
    print(f"  {v:4d}  {k}")
print("\n=== by divergence type (extra chain vs nearest ST chain) ===")
for k, v in Counter(r[6] for r in rows).most_common():
    print(f"  {v:4d}  {k}")
print("\n=== was the extra's exact chain ALSO extracted by ST? ===")
for k, v in Counter(('in_ST_extracted' if r[5] else 'NOT_in_ST') for r in rows).most_common():
    print(f"  {v:4d}  {k}")
print("\n=== cross-tab: source x divergence_type ===")
ct = Counter((r[2], r[6]) for r in rows)
for (src, dt), v in ct.most_common():
    print(f"  {v:4d}  {src:18s} {dt}")
# dump full table
with open("/tmp/phase0_attribution.tsv", "w") as out:
    out.write("tid\tclass\tsource\tnexons\tseed_tf\tin_ST_extracted\tdivergence\tshared_introns\tst_chain_nintrons\n")
    for r in rows:
        out.write("\t".join(str(x) for x in r) + "\n")
print("\nFull table: /tmp/phase0_attribution.tsv")
