#!/usr/bin/env python3
"""POSITIONAL & CLIPPING channel analysis on the 40 tied DAZ reads.

Channel signals interrogated, per read, for the DAZ1 placement vs the DAZ3 placement:
  (A) genomic 5'/3' ends of the alignment vs each copy's annotated TSS/TES (GFF)
  (B) soft/hard-clip lengths at the leading/trailing end of each placement's CIGAR
  (C) mapping-extent / coverage-shape asymmetry: aligned ref length, # blocks,
      fraction of read consumed, distance of each end to the copy boundary

The two genes are INVERTED paralogs, so DAZ1(-) and DAZ3(+) placements are reverse
complements of each other. We therefore compare each placement to ITS OWN copy's
annotated start/end, and ask whether one placement fits its copy's gene model
better than the other (which would be a copy-assigning signal).

Reads from /tmp/daz_tied_records.json (produced by tied_set_records.py).
"""
import json, re, collections

DAZ1 = (42783133, 42859657)   # - strand gene span
DAZ3 = (42879918, 42945552)   # + strand gene span

# annotated transcript extremes (union over all isoforms), from GFF:
# DAZ1 mRNAs all start 42783133; ends 42859562..42859657 -> use gene span.
# DAZ3 mRNAs 42879918..42945552.
# 5' TSS / 3' TES depend on strand:
#  DAZ1 (-) : TSS = high coord 42859657 ; TES = low coord 42783133
#  DAZ3 (+) : TSS = low  coord 42879918 ; TES = high coord 42945552
ANN = {
    'DAZ1': dict(strand='-', lo=DAZ1[0], hi=DAZ1[1], tss=DAZ1[1], tes=DAZ1[0]),
    'DAZ3': dict(strand='+', lo=DAZ3[0], hi=DAZ3[1], tss=DAZ3[0], tes=DAZ3[1]),
}

def nblocks(cig):
    return sum(1 for n, o in re.findall(r'(\d+)([MIDNSH=X])', cig) if o in 'M=X')

def n_introns(cig):
    return sum(1 for n, o in re.findall(r'(\d+)([MIDNSH=X])', cig) if o == 'N')

def matched_bases(cig):
    """matched/aligned READ bases only (M/=/X) -- excludes intron N (the bug fix)."""
    return sum(int(n) for n, o in re.findall(r'(\d+)([MIDNSH=X])', cig) if o in 'M=X')

def total_intron(cig):
    return sum(int(n) for n, o in re.findall(r'(\d+)([MIDNSH=X])', cig) if o == 'N')

recs = json.load(open('/tmp/daz_tied_records.json'))
qs = sorted(recs.keys())
print(f"analysing {len(qs)} tied reads\n")

# -------- per-read derived positional/clip features --------
rows = []
for q in qs:
    c = recs[q]
    a1, a3 = c['DAZ1'], c['DAZ3']
    row = dict(q=q)
    for tag, a in (('1', a1), ('3', a3)):
        ann = ANN['DAZ1' if tag == '1' else 'DAZ3']
        # alignment genomic ends (ref coords); pos is 1-based leftmost, endp exclusive-ish
        gstart, gend = a['pos'], a['endp']
        # The read's own 5' / 3' in genomic terms depends on the alignment strand.
        # On '+' alignment: read 5' = gstart, read 3' = gend.
        # On '-' alignment: read 5' = gend,  read 3' = gstart.
        if a['strand'] == '+':
            r5, r3 = gstart, gend
        else:
            r5, r3 = gend, gstart
        # distance of alignment ends to the copy's annotated TSS / TES
        d_tss = abs(r5 - ann['tss'])   # read-5' to annotated transcript start
        d_tes = abs(r3 - ann['tes'])   # read-3' to annotated transcript end
        # distance of alignment extremes to gene boundaries (extent fit)
        d_lo = gstart - ann['lo']      # +inside / -outside-left
        d_hi = ann['hi'] - gend        # +inside / -outside-right
        row[f'A{tag}'] = a['AS']
        row[f'strand{tag}'] = a['strand']
        row[f'lead{tag}'] = a['lead_clip']
        row[f'trail{tag}'] = a['trail_clip']
        row[f'alen{tag}'] = a['alen']                 # ref span incl. introns (raw)
        row[f'mb{tag}'] = matched_bases(a['cig'])     # matched READ bases (M/=/X) -- the real extent
        row[f'tin{tag}'] = total_intron(a['cig'])     # total intron length
        row[f'nblk{tag}'] = nblocks(a['cig'])
        row[f'nint{tag}'] = n_introns(a['cig'])
        row[f'd_tss{tag}'] = d_tss
        row[f'd_tes{tag}'] = d_tes
        row[f'd_lo{tag}'] = d_lo
        row[f'd_hi{tag}'] = d_hi
        row[f'sec{tag}'] = a['sec']
        row[f'mapq{tag}'] = a['mapq']
        row[f'seqlen{tag}'] = a['seqlen']
    rows.append(row)

# ================= CHANNEL A: ends vs annotated TSS/TES =================
print("=" * 78)
print("CHANNEL A: alignment ends vs each copy's annotated TSS / TES")
print("=" * 78)
print("Idea: if a read's 5'/3' ends snap to ONE copy's annotated start/end much")
print("better than the other's, that copy is favoured. Compare d_tss & d_tes.\n")
# signed: positive => DAZ1 placement fits its TSS/TES tighter (smaller distance)
sepA = 0
print(f"{'read':<16}{'5end d_tss D1/D3':>22}{'3end d_tes D1/D3':>22}{'call':>8}")
better_tss = collections.Counter(); better_tes = collections.Counter()
THRESH = 50   # bp: a meaningful end-snap difference for IsoSeq full-length
for r in rows:
    d1t, d3t = r['d_tss1'], r['d_tss3']
    d1e, d3e = r['d_tes1'], r['d_tes3']
    # which copy does the 5' end fit better?
    if abs(d1t - d3t) >= THRESH:
        better_tss['DAZ1' if d1t < d3t else 'DAZ3'] += 1
    if abs(d1e - d3e) >= THRESH:
        better_tes['DAZ1' if d1e < d3e else 'DAZ3'] += 1
    # combined call only if BOTH ends agree and exceed threshold
    call = '--'
    fav = set()
    if abs(d1t - d3t) >= THRESH: fav.add('DAZ1' if d1t < d3t else 'DAZ3')
    if abs(d1e - d3e) >= THRESH: fav.add('DAZ1' if d1e < d3e else 'DAZ3')
    if len(fav) == 1:
        call = fav.pop(); sepA += 1
    elif len(fav) == 2:
        call = 'CONFLICT'
for r in rows[:10]:
    print(f"{r['q'][-14:]:<16}"
          f"{str(r['d_tss1'])+'/'+str(r['d_tss3']):>22}"
          f"{str(r['d_tes1'])+'/'+str(r['d_tes3']):>22}")
print(f"\n  5'-end snaps better to a copy (|diff|>={THRESH}bp): {dict(better_tss)}")
print(f"  3'-end snaps better to a copy (|diff|>={THRESH}bp): {dict(better_tes)}")
print(f"  reads with a NON-CONFLICTING end-snap call: {sepA}/{len(rows)}")
print("  NOTE: because the copies are inverted near-identical, an end that snaps")
print("  to DAZ1's TSS necessarily lands the same number of bp from DAZ3's TES.")
print("  This is a SHARED geometry, not copy-specific signal -- see below.\n")

# Show the raw symmetry: for each read, is d_tss1 ~ d_tes3 and d_tes1 ~ d_tss3?
sym5 = sum(1 for r in rows if abs(r['d_tss1'] - r['d_tes3']) <= 5)
sym3 = sum(1 for r in rows if abs(r['d_tes1'] - r['d_tss3']) <= 5)
print(f"  inversion symmetry check: d_tss(DAZ1)==d_tes(DAZ3) within 5bp: {sym5}/{len(rows)}")
print(f"                            d_tes(DAZ1)==d_tss(DAZ3) within 5bp: {sym3}/{len(rows)}")

# ================= CHANNEL B: soft/hard clips =================
print("\n" + "=" * 78)
print("CHANNEL B: leading/trailing clip lengths, DAZ1 vs DAZ3 placement")
print("=" * 78)
print("Idea: if one placement clips far fewer read bases, the read 'fits' that")
print("copy better -> assignable. Compare total clip and per-end clip.\n")
sepB = 0
CLIP_THRESH = 20   # bp difference in total clipping to call a winner
print(f"{'read':<16}{'D1 lead/trail':>16}{'D3 lead/trail':>16}{'totD1':>7}{'totD3':>7}{'call':>8}")
clip_winner = collections.Counter()
for r in rows:
    t1 = r['lead1'] + r['trail1']
    t3 = r['lead3'] + r['trail3']
    call = '--'
    if abs(t1 - t3) >= CLIP_THRESH:
        call = 'DAZ1' if t1 < t3 else 'DAZ3'
        clip_winner[call] += 1
        sepB += 1
    r['_clipcall'] = call; r['_t1'] = t1; r['_t3'] = t3
for r in rows[:12]:
    print(f"{r['q'][-14:]:<16}"
          f"{str(r['lead1'])+'/'+str(r['trail1']):>16}"
          f"{str(r['lead3'])+'/'+str(r['trail3']):>16}"
          f"{r['_t1']:>7}{r['_t3']:>7}{r['_clipcall']:>8}")
print(f"\n  total-clip winner (|diff|>={CLIP_THRESH}bp): {dict(clip_winner)}")
print(f"  reads separated by clip asymmetry: {sepB}/{len(rows)}")
# distribution of |t1 - t3|
diffs = sorted(abs(r['_t1'] - r['_t3']) for r in rows)
print(f"  |clipD1 - clipD3| distribution: min={diffs[0]} median={diffs[len(diffs)//2]} max={diffs[-1]}")
print(f"  exact ties (clipD1==clipD3): {sum(1 for r in rows if r['_t1']==r['_t3'])}/{len(rows)}")

# ================= CHANNEL C: mapping extent / coverage shape =================
print("\n" + "=" * 78)
print("CHANNEL C: mapping-extent / coverage-shape asymmetry")
print("=" * 78)
print("Idea: aligned ref length, #exon-blocks, #introns -- if a read aligns with")
print("more matched bases / a cleaner structure to one copy, favour it.\n")
sepC = 0
print("EXTENT = MATCHED READ BASES (M/=/X) only -- intron N excluded (artifact fix).\n")
print(f"{'read':<16}{'matchD1/D3':>16}{'nblkD1/D3':>11}{'nintD1/D3':>11}{'intronD1/D3':>16}{'call':>8}")
ext_winner = collections.Counter()
MB_THRESH = 20  # bp difference in MATCHED read bases to call a winner
for r in rows:
    dm = r['mb1'] - r['mb3']
    call = '--'
    if abs(dm) >= MB_THRESH:
        call = 'DAZ1' if dm > 0 else 'DAZ3'
        ext_winner[call] += 1; sepC += 1
    r['_extcall'] = call
for r in rows[:12]:
    print(f"{r['q'][-14:]:<16}"
          f"{str(r['mb1'])+'/'+str(r['mb3']):>16}"
          f"{str(r['nblk1'])+'/'+str(r['nblk3']):>11}"
          f"{str(r['nint1'])+'/'+str(r['nint3']):>11}"
          f"{str(r['tin1'])+'/'+str(r['tin3']):>16}{r['_extcall']:>8}")
print(f"\n  MATCHED-BASE winner (|diff|>={MB_THRESH}bp): {dict(ext_winner)}")
mdiffs = sorted(abs(r['mb1'] - r['mb3']) for r in rows)
print(f"  |matchedD1 - matchedD3| dist: min={mdiffs[0]} median={mdiffs[len(mdiffs)//2]} max={mdiffs[-1]}")
print(f"  exact matched-base ties: {sum(1 for r in rows if r['mb1']==r['mb3'])}/{len(rows)}")
# the raw ref-span (incl introns) for contrast -- shows the earlier artifact
araw = sorted(abs(r['alen1'] - r['alen3']) for r in rows)
print(f"  [contrast] raw ref-span incl introns |diff|: median={araw[len(araw)//2]} max={araw[-1]} (intron-size diff, NOT extent)")
nblk_tie = sum(1 for r in rows if r['nblk1']==r['nblk3'])
nint_tie = sum(1 for r in rows if r['nint1']==r['nint3'])
print(f"  #block ties: {nblk_tie}/{len(rows)}; #intron ties: {nint_tie}/{len(rows)}")
print(f"  reads separated by matched-base extent: {sepC}/{len(rows)}")

# ================= AGGREGATE =================
print("\n" + "=" * 78)
print("AGGREGATE: any tied read assigned by ANY positional/clip channel?")
print("=" * 78)
assigned = []
for r in rows:
    chans = []
    # channel A only counts NON-conflicting end-snap (already symmetric => exclude)
    if r['_clipcall'] != '--': chans.append(('clip', r['_clipcall']))
    if r['_extcall'] != '--': chans.append(('extent', r['_extcall']))
    if chans:
        assigned.append((r['q'], chans))
print(f"  reads with a clip/extent channel call: {len(assigned)}")
for q, ch in assigned:
    print(f"    {q[-16:]}: {ch}")

json.dump([dict(q=r['q'], t1=r['_t1'], t3=r['_t3'], alen1=r['alen1'], alen3=r['alen3'],
                nblk1=r['nblk1'], nblk3=r['nblk3'], strand1=r['strand1'], strand3=r['strand3'],
                d_tss1=r['d_tss1'], d_tss3=r['d_tss3'], d_tes1=r['d_tes1'], d_tes3=r['d_tes3'],
                sec1=r['sec1'], sec3=r['sec3'], mapq1=r['mapq1'], mapq3=r['mapq3'])
           for r in rows], open('/tmp/daz_channel_rows.json', 'w'))
