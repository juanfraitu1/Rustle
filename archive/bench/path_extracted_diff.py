#!/usr/bin/env python3
"""Cov-parity gate for the read->transfrag accumulation port (Task 3 / Phase 2a).

Joins `path_extracted` events from Rustle and StringTie by (strand, intron-chain)
and reports the WITHIN-LOCUS coverage ratio for each chain in each tool. ST `cov`
is tlen-scaled and not directly comparable to Rustle's flow cov, so we normalize
each chain's cov by the locus-dominant cov (max cov among chains sharing the locus,
per tool). The within-locus ratio is therefore tool-invariant and comparable.

A "contested minority" is a chain that is:
  - a reference TP (intron chain present in the reference GTF),
  - extracted by Rustle (present in ru path_extracted),
  - killed under RUSTLE_PREDCLUSTER_ST by isofrac or retained_intron (ru pred_kill),
  - kept / extracted by ST (present in st path_extracted),
  - where ST within-locus cov ratio is materially HIGHER than Rustle's.

Usage:
  path_extracted_diff.py RU_PE.jsonl ST_PE.jsonl REF.gtf [--killed-only]

  RU_PE.jsonl  Rustle parity log with steps path_extracted,pred_kill
  ST_PE.jsonl  StringTie parity log with step path_extracted
  REF.gtf      reference annotation (defines TP chains)

With --killed-only, restrict the contested-minority table to chains Rustle killed.
Otherwise also print all-chain cov-ratio parity (movement under RUSTLE_FLOW_ST).
"""
import sys, json, collections

def parse_pe(path):
    """Return list of (strand, chain_tuple, cov, start, end, killed_reason or None)."""
    extracted = {}   # (strand, chain) -> dict(cov, start, end)
    killed = {}      # (strand, chain) -> reason  (only isofrac/retained_intron)
    killed_any = {}  # (strand, chain) -> reason  (any reason)
    for line in open(path):
        line = line.strip()
        if not line:
            continue
        try:
            e = json.loads(line)
        except Exception:
            continue
        step = e.get("step")
        p = e.get("payload", {})
        strand = e.get("strand", ".")
        introns = p.get("introns", "")
        chain = tuple(introns.split(",")) if introns else ()
        key = (strand, chain)
        if step == "path_extracted":
            cov = float(p.get("cov", 0.0))
            # keep the max-cov extraction per chain (chain may be extracted twice)
            prev = extracted.get(key)
            if prev is None or cov > prev["cov"]:
                extracted[key] = {"cov": cov, "start": e.get("start"), "end": e.get("end")}
        elif step == "pred_kill":
            reason = str(p.get("reason", "?"))
            killed_any[key] = reason
            if reason in ("isofrac", "retained_intron"):
                # pred_kill carries no introns in the rustle payload; derive chain
                # from coordinates by matching against extracted set later.
                killed[key] = reason
    return extracted, killed, killed_any

def parse_killed_by_coord(path):
    """Rustle pred_kill payload has no introns; capture (strand,start,end,reason)."""
    rows = []
    for line in open(path):
        line = line.strip()
        if not line:
            continue
        try:
            e = json.loads(line)
        except Exception:
            continue
        if e.get("step") != "pred_kill":
            continue
        p = e.get("payload", {})
        reason = str(p.get("reason", "?"))
        rows.append((e.get("strand", "."), e.get("start"), e.get("end"), reason,
                     float(p.get("cov", 0.0)), p.get("nexons")))
    return rows

def parse_ref_chains(path):
    """Return set of (strand, chain_tuple) for reference transcripts (>=2 exons)."""
    tx_exons = collections.defaultdict(list)
    tx_strand = {}
    for line in open(path):
        if line.startswith("#"):
            continue
        f = line.rstrip("\n").split("\t")
        if len(f) < 9 or f[2] != "exon":
            continue
        attrs = f[8]
        tid = None
        for part in attrs.split(";"):
            part = part.strip()
            if part.startswith("transcript_id"):
                tid = part.split('"')[1]
                break
        if tid is None:
            continue
        tx_exons[tid].append((int(f[3]), int(f[4])))
        tx_strand[tid] = f[6]
    ref = set()
    for tid, exons in tx_exons.items():
        exons.sort()
        if len(exons) < 2:
            continue
        introns = tuple(f"{exons[i][1]+1}-{exons[i+1][0]-1}" for i in range(len(exons)-1))
        ref.add((tx_strand[tid], introns))
    return ref

def assign_loci(keys, starts_ends):
    """Group (strand,chain) keys into loci by overlapping spans per strand."""
    # build (strand) -> sorted list of (start,end,key)
    by_strand = collections.defaultdict(list)
    for k in keys:
        se = starts_ends.get(k)
        if se is None:
            continue
        by_strand[k[0]].append((se["start"], se["end"], k))
    locus_of = {}
    for strand, items in by_strand.items():
        items.sort()
        cur_end = None
        locus_id = -1
        for s, e, k in items:
            if cur_end is None or s > cur_end:
                locus_id += 1
                cur_end = e
            else:
                cur_end = max(cur_end, e)
            locus_of[k] = (strand, locus_id)
    return locus_of

def within_locus_ratio(extracted):
    """For each extracted key, cov / locus-dominant cov."""
    starts_ends = {k: {"start": v["start"], "end": v["end"]} for k, v in extracted.items()}
    locus_of = assign_loci(extracted.keys(), starts_ends)
    dom = {}  # locus -> max cov
    for k, v in extracted.items():
        loc = locus_of.get(k)
        if loc is None:
            continue
        dom[loc] = max(dom.get(loc, 0.0), v["cov"])
    ratio = {}
    for k, v in extracted.items():
        loc = locus_of.get(k)
        if loc is None or dom.get(loc, 0.0) <= 0:
            ratio[k] = 0.0
        else:
            ratio[k] = v["cov"] / dom[loc]
    return ratio, locus_of

def main():
    args = [a for a in sys.argv[1:] if not a.startswith("--")]
    killed_only = "--killed-only" in sys.argv
    if len(args) < 3:
        print(__doc__)
        sys.exit(1)
    ru_path, st_path, ref_path = args[0], args[1], args[2]

    ru_ext, _, ru_killed_any = parse_pe(ru_path)
    st_ext, _, _ = parse_pe(st_path)
    ref = parse_ref_chains(ref_path)
    ru_killed_rows = parse_killed_by_coord(ru_path)

    ru_ratio, _ = within_locus_ratio(ru_ext)
    st_ratio, _ = within_locus_ratio(st_ext)

    # Map killed coordinates -> chains by matching to ru extracted spans.
    # A pred_kill with isofrac/retained_intron near an extracted chain's span marks it killed.
    killed_chains = {}  # (strand,chain) -> reason
    ext_by_span = collections.defaultdict(list)
    for k, v in ru_ext.items():
        ext_by_span[(k[0], v["start"], v["end"])].append(k)
    for strand, s, e, reason, cov, nex in ru_killed_rows:
        if reason not in ("isofrac", "retained_intron"):
            continue
        # match exact span first
        for k in ext_by_span.get((strand, s, e), []):
            killed_chains[k] = reason
        # also tolerant: same strand, overlapping span, cov close
        for k, v in ru_ext.items():
            if k[0] != strand:
                continue
            if abs(v["start"] - s) <= 5 and abs(v["end"] - e) <= 5:
                killed_chains.setdefault(k, reason)

    # Contested minorities.
    print("=== Cov-parity (path_extracted within-locus ratio) ===")
    print(f"ref TP chains: {len(ref)}  ru extracted: {len(ru_ext)}  st extracted: {len(st_ext)}")
    print(f"ru killed (isofrac/RI) chains matched: {len(killed_chains)}")
    contested = []
    for k, reason in killed_chains.items():
        if k not in ref:
            continue
        if k not in st_ext:
            continue
        rr = ru_ratio.get(k, 0.0)
        sr = st_ratio.get(k, 0.0)
        if sr > rr + 1e-9:
            contested.append((k, reason, rr, sr, sr - rr))
    contested.sort(key=lambda x: -x[4])
    print(f"\nCONTESTED MINORITIES (ref TP, ru-killed by isofrac/RI, kept by ST, ST ratio > RU ratio): {len(contested)}")
    print(f"{'strand':6} {'reason':16} {'ru_ratio':>9} {'st_ratio':>9} {'delta':>8}  chain[:1]")
    for k, reason, rr, sr, d in contested[:60]:
        ch = k[1][0] if k[1] else "-"
        print(f"{k[0]:6} {reason:16} {rr:9.4f} {sr:9.4f} {d:8.4f}  {ch}")
    if contested:
        import statistics
        print(f"\nsummary: n={len(contested)} mean ru_ratio={statistics.mean(c[2] for c in contested):.4f} "
              f"mean st_ratio={statistics.mean(c[3] for c in contested):.4f} "
              f"mean gap={statistics.mean(c[4] for c in contested):.4f}")
    return contested

if __name__ == "__main__":
    main()
