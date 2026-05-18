#!/bin/bash

echo "=========================================="
echo "Precision Diagnostic: Identify False Positives"
echo "=========================================="
echo ""

# Step 1: Generate baseline
if [ ! -f /tmp/rustle_precision_baseline.gtf ]; then
    echo "Step 1: Generating Rustle baseline GTF..."
    ./target/release/rustle -L ../GGO_19.bam -o /tmp/rustle_precision_baseline.gtf 2>&1 | tail -3
    echo "✓ Baseline generated"
else
    echo "✓ Baseline already exists"
fi

echo ""
echo "Step 2: Running precision analysis..."
python3 precision_analysis.py > precision_diagnostic.txt 2>&1

echo "Step 3: Summary of findings"
echo "=========================================="
head -50 precision_diagnostic.txt

echo ""
echo "=========================================="
echo "Full analysis saved to: precision_diagnostic.txt"
echo ""
echo "Next steps:"
echo "1. Review FP categories in precision_diagnostic.txt"
echo "2. Test filters in PRECISION_IMPROVEMENT_PLAN.md"
echo "3. Run: ./test_precision_filters.sh"
