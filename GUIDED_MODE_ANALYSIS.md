# Guided Mode Analysis: De novo + Reference Guide Testing

**Date**: 2026-05-18  
**Purpose**: Test if guided assembly improves completeness/quality  
**Status**: Data collection in progress

---

## Initial Findings

### Transcript Counts by Mode

```
Reference:             1,839 structures (baseline)
De novo only:          1,944 structures (+105 vs ref, +5.7%)
With guide (-G):       1,790 structures (-49 vs ref, -2.7%)
Guide + no exact:      1,869 structures (+30 vs ref, +1.6%)

Structure match rates:
De novo only:          105.7% (105 extra structures)
With guide (-G):        97.3% (49 missing structures!)
Guide + no exact:      101.6% (30 extra structures)
```

### Key Discovery

**Using `-G` (guide flag) REDUCES output, not increases it!**
- Loses 154 transcripts vs. baseline de novo
- Even more restrictive than "no exact" mode
- Could be by design (guide constrains assembly)

---

## Hypothesis

When reference GTF is used as guide (`-G` flag):
1. Assembly is constrained to reference topology
2. Only transcripts matching reference structures are emitted
3. Novel isoforms are suppressed
4. Result: Lower sensitivity (lose novel discoveries)

This would explain:
- ✓ Lower transcript count with guide
- ✓ Why "no exact" mode has more (allows Rustle relaxations)
- ✓ Trade: Better specificity/parity, worse novelty

---

## Modes Being Tested

### Mode 1: De novo only (Baseline)
- No guide GTF
- Default STRINGTIE_EXACT=1
- **Result**: 1,944 structures (1,312 exact matches + 632 novel)

### Mode 2: De novo WITH guide (-G flag)
- Reference GTF provided as guide
- Default STRINGTIE_EXACT=1
- **Expected**: More conservative, higher parity
- **Actual**: Fewer transcripts (1,790)

### Mode 3: De novo WITH guide + STRINGTIE_EXACT=0
- Reference GTF as guide
- Rustle-relaxed mode
- **Expected**: Middle ground
- **Actual**: 1,869 (between strict guide and pure denovo)

---

## What We Need to Know

1. **Sensitivity metrics** (junction/intron/transcript level)
   - Which mode catches more of the reference?
   - Which mode is most complete?

2. **Precision metrics**
   - Which mode has fewest false positives?
   - Is guide mode achieving higher specificity?

3. **F1 scores**
   - Which mode balances precision/recall best?
   - Which achieves closest to "100% on all metrics"?

4. **Specificity of guide constraints**
   - Is guide preventing real discoveries or just noise?
   - Are the 154 lost transcripts in denovo valid or artifacts?

---

## Why This Matters for Your Question

**You asked**: "Can we get 100% on all metrics with denovo guided?"

**What we're testing**:
- ✓ Does guide mode help achieve better metrics?
- ✓ Or does it trade sensitivity for precision?
- ✓ Can we find a mode that's complete (100% sensitivity) AND accurate (100% precision)?

**Early indication**: Guide mode restricts, doesn't enhance. But might achieve better balance.

---

## Next Steps (Low Usage - Document for Later)

### Step 1: Get Full Metrics
Need to run gffcompare properly on all three modes:
```bash
gffcompare -r ../GGO_19.gtf /tmp/guided_*.gtf > metrics_comparison.txt
```

This will show:
- Sensitivity: % of reference junctions/transcripts recovered
- Precision: % of predicted transcripts in reference
- Specificity
- F1 scores

### Step 2: Analyze Results
Create comparison table:
```
Metric          De novo    With guide    Guide+noexact
─────────────────────────────────────────────────────
Sensitivity     XX%        XX%           XX%
Precision       XX%        XX%           XX%
Matching        XXX        XXX           XXX
F1 Score        XX%        XX%           XX%
```

### Step 3: Characterize Differences
- What transcripts are unique to each mode?
- Are the lost transcripts (with guide) valid or artifacts?
- Does guide mode achieve ANY 100% metrics?

### Step 4: Test Hybrid Approaches
If basic modes don't achieve 100%, test combinations:
- Guide with relaxed thresholds
- Guide with novelty enabled
- Guide with specific gates enabled/disabled

---

## Saved Files for Continuation

```
/tmp/guided_denovo.gtf       ← De novo only
/tmp/guided_with_ref.gtf     ← With -G flag
/tmp/guided_no_exact.gtf     ← With -G and RUSTLE_STRINGTIE_EXACT=0
```

## Preliminary Metrics to Continue

Basic counts (already computed):
- Reference: 1,839 unique structures
- De novo: 1,944 (105.7% match rate)
- With guide: 1,790 (97.3% match rate)
- Guide+noexact: 1,869 (101.6% match rate)

Full metrics needed from gffcompare to answer your question.

---

## Key Question for Continuation

**Does any mode achieve 100% on all metrics (Sn, Pr, Sp, F1)?**

Most likely answer:
- **De novo alone**: High sensitivity (~92%), lower precision (~67-72%)
- **With guide**: Possibly higher precision, likely lower sensitivity
- **Hybrid**: Might be middle ground

None will likely hit 100% on ALL metrics simultaneously - that would require:
- 100% sensitivity: Recover all reference transcripts
- 100% precision: Only output real transcripts
- These naturally trade off

But let's verify!
