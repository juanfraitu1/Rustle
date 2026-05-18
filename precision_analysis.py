#!/usr/bin/env python3
"""
Precision Analysis: Identify false positives (transcripts Rustle has that StringTie doesn't)
and characterize them to find patterns for elimination.
"""

def extract_transcripts(gtf_path):
    """Extract all transcripts with their exon structures."""
    transcripts = {}

    with open(gtf_path) as f:
        for line in f:
            if line.startswith('#') or not line.strip():
                continue
            parts = line.split('\t')
            if len(parts) < 9 or parts[2] != 'exon':
                continue

            attrs = parts[8]
            tx_id = None
            if 'transcript_id "' in attrs:
                tx_id = attrs.split('transcript_id "')[1].split('"')[0]

            if not tx_id:
                continue

            chrom = parts[0]
            strand = parts[6]
            start = int(parts[3])
            end = int(parts[4])

            if tx_id not in transcripts:
                transcripts[tx_id] = {
                    'chrom': chrom,
                    'strand': strand,
                    'exons': [],
                    'coverage': 0,
                }

            transcripts[tx_id]['exons'].append((start, end))

            # Try to extract coverage if available
            if 'cov "' in attrs:
                try:
                    cov_str = attrs.split('cov "')[1].split('"')[0]
                    transcripts[tx_id]['coverage'] = float(cov_str)
                except:
                    pass

    # Normalize exons
    for tx_id in transcripts:
        transcripts[tx_id]['exons'] = tuple(sorted(transcripts[tx_id]['exons']))

    return transcripts

def analyze_precision(ref_txs, rustle_txs):
    """Analyze precision: identify FPs and characterize them."""

    print("=" * 80)
    print("PRECISION ANALYSIS: False Positive Characterization")
    print("=" * 80)

    # Convert to exon-structure sets for comparison
    ref_structures = set()
    ref_by_structure = {}
    for tx_id, tx in ref_txs.items():
        struct = tx['exons']
        ref_structures.add(struct)
        if struct not in ref_by_structure:
            ref_by_structure[struct] = []
        ref_by_structure[struct].append(tx_id)

    rustle_structures = set()
    rustle_by_structure = {}
    for tx_id, tx in rustle_txs.items():
        struct = tx['exons']
        rustle_structures.add(struct)
        if struct not in rustle_by_structure:
            rustle_by_structure[struct] = []
        rustle_by_structure[struct].append(tx_id)

    # Find false positives (in Rustle but not in reference)
    fp_structures = rustle_structures - ref_structures

    print(f"\nReference structures: {len(ref_structures)}")
    print(f"Rustle structures: {len(rustle_structures)}")
    print(f"Exact matches: {len(ref_structures & rustle_structures)}")
    print(f"False positive structures: {len(fp_structures)}")
    print(f"Precision: {100 * len(ref_structures & rustle_structures) / len(rustle_structures):.1f}%")

    if not fp_structures:
        print("\n✅ Perfect precision! No false positives.")
        return

    # Characterize FPs
    print(f"\nAnalyzing {len(fp_structures)} false positive structures...")

    fp_txs = []
    for struct in fp_structures:
        for tx_id in rustle_by_structure[struct]:
            tx = rustle_txs[tx_id]
            fp_txs.append({
                'tx_id': tx_id,
                'exons': len(struct),
                'coverage': tx['coverage'],
                'span': struct[-1][1] - struct[0][0] if struct else 0,
                'structure': struct,
                'chrom': tx['chrom'],
                'strand': tx['strand'],
            })

    # Sort by coverage (likely artifacts have low coverage)
    fp_txs.sort(key=lambda x: x['coverage'])

    print(f"\nFalse positives by coverage (lowest first):")
    print(f"{'Rank':>5} {'TX_ID':20} {'Exons':>6} {'Coverage':>10} {'Span':>10}")
    print("-" * 60)

    for i, fp in enumerate(fp_txs[:20]):
        print(f"{i+1:5d} {fp['tx_id']:20} {fp['exons']:6d} {fp['coverage']:10.2f} {fp['span']:10d}")

    # Categorize FPs
    print(f"\nFalse positive categorization:")

    single_exon = sum(1 for fp in fp_txs if fp['exons'] == 1)
    low_cov = sum(1 for fp in fp_txs if fp['coverage'] < 1.0)
    low_exons = sum(1 for fp in fp_txs if fp['exons'] < 3)
    short_span = sum(1 for fp in fp_txs if fp['span'] < 1000)

    print(f"  Single-exon: {single_exon}/{len(fp_txs)} ({100*single_exon/len(fp_txs):.1f}%)")
    print(f"  Low coverage (<1.0): {low_cov}/{len(fp_txs)} ({100*low_cov/len(fp_txs):.1f}%)")
    print(f"  Few exons (<3): {low_exons}/{len(fp_txs)} ({100*low_exons/len(fp_txs):.1f}%)")
    print(f"  Short span (<1kb): {short_span}/{len(fp_txs)} ({100*short_span/len(fp_txs):.1f}%)")

    # Overlapping FPs with ref structures (partial matches)
    print(f"\nPartial overlaps with reference:")
    partial_overlaps = 0
    for fp_struct in list(fp_structures)[:10]:
        for ref_struct in ref_structures:
            shared = len(set(fp_struct) & set(ref_struct))
            if shared > 0 and shared < len(fp_struct):
                partial_overlaps += 1
                print(f"  FP {len(fp_struct)} exons shares {shared} with ref {len(ref_struct)} exons")
                break

    print(f"\nRecommendations:")
    print(f"1. Lower min_transcript_length threshold")
    print(f"2. Raise readthr_gate (min coverage requirement)")
    print(f"3. Add single-exon specificity filters")
    print(f"4. Review coverage calculation (may be inflated)")

if __name__ == '__main__':
    print("Loading reference GTF...")
    ref_txs = extract_transcripts('../GGO_19.gtf')

    print("Loading Rustle GTF...")
    rustle_txs = extract_transcripts('/tmp/rustle_baseline.gtf')

    analyze_precision(ref_txs, rustle_txs)
