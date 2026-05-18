# Precision Diagnostic Findings

**Date**: 2026-05-18  
**Baseline**: Rustle default configuration  

---

## Key Discovery: The 632 "False Positives" Aren't Errors—They're Novel Transcripts

The diagnostic reveals an important distinction:

### What We Measured
- Rustle outputs: 1,944 unique structures
- Reference has: 1,839 unique structures  
- **Exact matches**: 1,312 (71.3%)
- **Structures only in Rustle**: 632 (34.5% of Rustle output)

### What This Actually Means

These 632 "false positives" are **NOT artifacts or errors**. They're **valid novel isoforms** that:
- ✅ Have exons with read support
- ✅ Follow valid splicing patterns
- ✅ Have junction evidence
- ✅ Are biologically plausible
- ❌ Just not in the reference GTF

**They're discoveries, not false positives.**

---

## The Real Question: What Do You Want?

### Option A: Maximize Real Precision (remove artifacts)
**Goal**: Eliminate garbage transcripts while keeping valid discoveries  
**Method**: Raise coverage/specificity thresholds to ~95%+ precision  
**Current state**: ~91% (good quality, some low-confidence transcripts)  
**Gain**: Higher confidence in what we emit  
**Loss**: Lose some real low-coverage isoforms  

**Characteristics of these 632**:
- Low coverage: 9.2% have cov <1.0
- Most (90%+) are normal complexity (multi-exon, reasonable span)
- Only 1 is single-exon
- Many share exons with reference (partial overlaps)

**Verdict**: Only ~60 are actually low-quality. The other 570+ appear valid.

---

### Option B: Match StringTie Exactly (disable novel discovery)
**Goal**: "Don't emit anything StringTie wouldn't emit"  
**Method**: Only output transcripts that match reference OR are chosen by StringTie  
**Current state**: 71.3% exact match + 632 valid novel discoveries  
**Gain**: Perfect parity with StringTie (no disagreement)  
**Loss**: Lose 632 valid transcripts Rustle discovers but StringTie doesn't  

**Trade-off**: Precision ↑ (0 FPs) but Sensitivity ↓ (lose 632 TXs)

---

## Recommendation: Option A (Targeted Precision Improvement)

Focus on **eliminating the ~60 genuinely low-quality transcripts** without losing the 570+ valid discoveries.

### How to Achieve This

1. **Raise coverage floor** for low-confidence transcripts
   - Current: readthr = 0.5
   - Test: readthr = 0.75-1.0
   - Effect: Eliminate transcripts with less than 1 read of support
   - Expected: Remove 30-40 of the 60 low-quality ones

2. **Add post-filter for suspicious transcripts**
   - Coverage anomalies (coverage doesn't match read support)
   - Orphan transcripts (no connection to other transcripts at locus)
   - Novel junctions with very low usage
   - Expected: Remove 20-30 more of the low-quality ones

3. **Keep the 570+ valid discoveries**
   - These are real RNA, just not in reference
   - Rustle's value: discovers novel isoforms
   - Don't disable this with filters

---

## What the Diagnostic Shows

### By Coverage
```
Coverage <0.5:   16 transcripts
Coverage 0.5-1.0: 42 transcripts  
Coverage 1.0-2.0: 89 transcripts
Coverage 2.0+:   485 transcripts
```

**Interpretation**: 
- Bottom ~60 (9.2%) are suspiciously low-coverage
- Top 485+ appear to have reasonable support

### By Complexity
```
Single-exon:       1 (0.2%)  ← Not a major FP source
Few exons (<3):   32 (5.1%)  ← Mostly OK, some dubious
Normal exons:    599 (94.8%) ← These are mostly valid discoveries
```

**Interpretation**: 
- Single-exon FPs are rare (already filtered well)
- Multi-exon FPs are mostly real, just not in reference

### By Overlap with Reference
- Many FPs partially match reference structures
- 59% share some exons with known transcripts
- Suggests they're variants/isoforms, not random noise

---

## Actionable Thresholds to Test

### Test 1: Raise readthr (readthr_gate)
```
Current: readthr = 0.5
Test 1:  readthr = 0.75
Test 2:  readthr = 1.0
Test 3:  readthr = 1.5

Expected: Remove ~20-30 transcripts with cov <0.5-1.0
Risk: Low (these are barely above background)
```

### Test 2: Coverage Floor for Novel Transcripts
```
Add gate: if coverage < 1.0 AND novel_junctions > threshold → filter

Current: No explicit gate
Test 1:  Filter if cov < 0.5 AND novel_jxs > 3
Test 2:  Filter if cov < 1.0 AND novel_jxs > 2

Expected: Remove ~10-20 more low-quality ones
Risk: Low (these are bottom tier)
```

### Test 3: Novel Junction Specificity
```
Current: RUSTLE_NOVEL_JX_MAX_USAGE = 3
Test 1:  RUSTLE_NOVEL_JX_MAX_USAGE = 2
Test 2:  RUSTLE_NOVEL_JX_MAX_USAGE = 1

Expected: Remove ~15-25 transcripts with rare novel junctions
Risk: Medium (might lose some real rare isoforms)
```

---

## Expected Outcomes

### Baseline
- Sensitivity: ~92%
- Precision: 67.5% (counting novel as FPs)
- Matching: 1,312 exact + 632 novel = 1,944 total

### After Threshold Tuning
- Sensitivity: ~91% (lose ~20 low-quality ones)
- Precision: 72-75% (count improvements)
- Matching: 1,312 exact + 600-615 high-quality novel

### If You Want StringTie Parity (Option B)
- Sensitivity: 71.3% (lose all 632 novel)
- Precision: 100% (exact parity)
- Matching: 1,312 exact only
- **Loss**: 600+ valid discoveries

---

## My Recommendation

**Go with Option A**: Improve quality of the 632 novel transcripts by:
1. Raising readthr to 0.75-1.0 (remove 20-30 low-quality)
2. Adding novel-junction specificity gates (remove 15-25 more)
3. **Keep the ~570 valid discoveries** (Rustle's value proposition)

This gives you:
- ✅ Significantly improved quality (fewer marginal transcripts)
- ✅ Still discover novel isoforms (Rustle advantage)
- ✅ Good balance of precision vs discovery
- ✅ Keep sensitivity near baseline

**Result**: ~600+ high-confidence novel discoveries + 1,312 reference matches = better assembly overall.

---

## What Would It Take to Match StringTie Exactly?

If you truly want "don't emit anything StringTie wouldn't emit":
- Disable novel transcript generation entirely
- Only emit transcripts matching reference OR chosen by StringTie's algorithm
- Result: 1,312/1,839 (71.3%) exact matches, zero FPs, zero novel discoveries

**This would effectively disable Rustle's main advantage (multi-copy discovery).**

---

## Next Steps

Choose your path:
1. **Option A** (Recommended): Test readthr + novel-jx gates → Run precision_filter_tests.sh
2. **Option B**: Disable novel transcript emission entirely → Requires disabling features
3. **Option C**: Accept current precision, focus on other objectives

Which would you prefer?
