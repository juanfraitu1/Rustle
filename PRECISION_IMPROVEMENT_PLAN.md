# Precision Improvement Plan: Eliminate False Positives

**Goal**: "Don't emit anything StringTie wouldn't emit"  
**Current baseline**: ~91% precision  
**Target**: 95%+ precision (or match StringTie exactly)

---

## The Challenge

Improving precision is **harder than improving sensitivity** because:
- Any filter we add might also kill true positives
- We lose some valid but low-confidence transcripts
- Tradeoff: better precision = lower sensitivity

But the goal is clear: **eliminate false positives without hurting true positives**.

---

## Types of False Positives We Can Target

Based on the codebase, likely FP categories:

### 1. **Low-Coverage Artifacts** (easiest to eliminate)
- Coverage too low for true transcript (readthr_gate issue)
- Solution: Raise `readthr` threshold or add post-filter

### 2. **Single-Exon False Discoveries**
- Novel SE genes that look good to Rustle but aren't real
- Solution: Stricter SE-specific gates (higher cov requirement)

### 3. **Chimeric/Bridge Transcripts**
- Already filtered by chimera filter (commit 53e992f)
- Could be further tuned

### 4. **Short Transcripts**
- May be noise/incomplete assemblies
- Solution: Raise min_transcript_length

### 5. **High-Novelty, Low-Support Transcripts**
- Novel junctions with low read support
- Solution: RUSTLE_NOVEL_JX_MAX_USAGE already targets this

### 6. **Coverage Calculation Issues**
- Rustle's coverage might be inflated compared to StringTie
- Solution: Align coverage calculation or scale thresholds

---

## Three-Phase Approach

### Phase 1: Diagnostic (30 min)
1. Run `precision_analysis.py` to categorize FPs
2. Identify top 3 FP categories by frequency
3. Measure what % of FPs each category represents

### Phase 2: Targeted Filtering (1-2 hours)
For each top FP category:
1. Design a gate/filter that catches it
2. Test with gffcompare
3. Measure precision gain vs sensitivity loss
4. Keep only beneficial gates (Pr gain > Sn loss)

### Phase 3: Integration (30 min)
1. Combine beneficial gates
2. Test final configuration
3. Document threshold changes
4. Verify no regressions

---

## Candidate Filters to Test

### Filter 1: Raise readthr_gate
```rust
// Current: config.readthr = 0.5
// Test: config.readthr = 1.0 or higher
// Effect: Eliminate low-coverage transcripts
// Risk: Lose some true low-abundance isoforms
```

### Filter 2: Single-Exon Coverage Threshold
```rust
// Add: Min coverage for single-exon transcripts
// Current gate: config.singlethr = 0.5
// Test: Raise to 1.0, 2.0, or higher
// Effect: More selective about single-exon discoveries
```

### Filter 3: Short Transcript Length
```rust
// Current: min_transcript_length = 0 (no filter)
// Test: 100bp, 200bp, 500bp minimum
// Effect: Eliminate very short noise
// Risk: Some real short exons are valid
```

### Filter 4: Coverage Calculation
```rust
// Review: How is coverage calculated?
// Current: flow-based (nodeflux * noderate) / txlen
// Compare: StringTie's coverage calculation
// Align if different
```

### Filter 5: Novel Junction Specificity
```rust
// Current: RUSTLE_NOVEL_JX_MAX_USAGE = 3
// Test: Raise to 2, 1, or 0 (no novel junctions)
// Effect: Only use reference or highly-supported junctions
```

---

## Implementation Plan

### Step 1: Diagnostic Run
```bash
# Generate baseline if needed
./target/release/rustle -L ../GGO_19.bam -o /tmp/rustle_baseline.gtf

# Run analysis
python3 precision_analysis.py > precision_diagnostic.txt

# Review: What are the top 3 FP categories?
```

### Step 2: For Each Candidate Filter
```bash
# Test variant
RUSTLE_READTHR=1.0 ./target/release/rustle -L ../GGO_19.bam -o /tmp/test.gtf

# Compare metrics
gffcompare -r ../GGO_19.gtf /tmp/test.gtf | tail -1
# Note: Sensitivity, Precision, # matching

# Decision: Keep if (Pr_gain > 1%) AND (Sn_loss < 0.5%)
```

### Step 3: Combined Configuration
```bash
# Test all beneficial gates together
RUSTLE_READTHR=1.0 \
RUSTLE_SINGLETHR=1.0 \
RUSTLE_MIN_TRANSCRIPT_LENGTH=200 \
./target/release/rustle -L ../GGO_19.bam -o /tmp/final.gtf

# Validate
gffcompare -r ../GGO_19.gtf /tmp/final.gtf
```

---

## Expected Outcomes

### Optimistic (High Probability)
- Find 2-3 simple threshold relaxations that improve precision
- Precision: 91% → 93-95%
- Sensitivity loss: <1pp (acceptable)
- Configuration: 2-3 parameter changes

### Realistic (Medium Probability)
- Precision improves to 92-93%
- Requires tuning multiple gates
- Some sensitivity loss (1-2pp)
- Configuration: 3-4 parameter changes

### Conservative (Lower Probability)
- Precision plateaus at 91-92%
- Large sensitivity tradeoff needed for further improvement
- Suggests 91% is good balance point

---

## Key Metrics to Track

For each test:
```
Baseline (current):
  Sensitivity: XX%
  Precision:   XX%
  Matching:    XXX transcripts
  
Test variant:
  Sensitivity: XX% (Δ=±Y%)
  Precision:   XX% (Δ=±Y%)
  Matching:    XXX transcripts
  
Decision: Keep if Pr gain > Sn loss AND improves overall F1
```

---

## Timeline

| Phase | Estimated |
|-------|-----------|
| Diagnostic | 30 min |
| Filter testing | 90 min |
| Integration | 30 min |
| **Total** | **2.5 hours** |

---

## Success Criteria

- ✓ Identify top 3 FP categories (diagnostic)
- ✓ Test 3-5 candidate filters
- ✓ Achieve 93%+ precision if possible
- ✓ Document final configuration
- ✓ Verify no major sensitivity regression (<1pp loss acceptable)

---

## Related Infrastructure

Already in codebase:
- RUSTLE_CHIMERA_FILTER_GTF (chimeric transcripts)
- RUSTLE_NOVEL_JX_MAX_USAGE (novel junction specificity)
- readthr_gate (coverage threshold)
- singlethr (single-exon threshold)
- min_transcript_length (length threshold)

New thresholds can use existing env vars or add new ones.

---

## Next: Run Diagnostic

Ready to:
1. Generate baseline GTF
2. Run precision_analysis.py
3. Identify top FP categories
4. Test filters in priority order

Let's start!
