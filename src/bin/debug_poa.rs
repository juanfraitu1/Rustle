//! Debug poasta poa_msa on real GOLGA6L7 exon-0 sequences.
use rustle::genome::GenomeIndex;

use poasta::graphs::poa::POAGraph;
use poasta::aligner::PoastaAligner;
use poasta::aligner::config::AffineDijkstra;
use poasta::aligner::scoring::{GapAffine, AlignmentType};
use poasta::io::fasta::poa_graph_to_fasta;

fn main() -> anyhow::Result<()> {
    let genome = GenomeIndex::from_fasta("/scratch/jxi21/Assembler/GGO.fasta")?;
    let chrom = "NC_073243.2";
    let s0 = genome.fetch_sequence(chrom, 104789646, 104789918).unwrap();
    let s1 = genome.fetch_sequence(chrom, 104830535, 104830737).unwrap();
    let s2 = genome.fetch_sequence(chrom, 104871355, 104871545).unwrap();
    println!("len0={} len1={} len2={}", s0.len(), s1.len(), s2.len());

    let gap_costs = GapAffine::new(1, 1, 2);
    let aligner: PoastaAligner<'_, AffineDijkstra> =
        PoastaAligner::new(AffineDijkstra(gap_costs), AlignmentType::Global);

    let mut graph: POAGraph<u32> = POAGraph::new();
    let unit_weights: Vec<usize> = vec![1; s0.len().max(s1.len()).max(s2.len())];

    let inputs = [&s0, &s1, &s2];
    for (i, seq) in inputs.iter().enumerate() {
        let w = &unit_weights[..seq.len()];
        if graph.is_empty() {
            graph.add_alignment_with_weights(&format!("seq{i}"), seq, None, w).unwrap();
        } else {
            let aln_result = aligner.align::<u32, _>(&graph, seq);
            let alignment = if aln_result.alignment.is_empty() { None }
                else { Some(&aln_result.alignment as &poasta::aligner::Alignment<poasta::graphs::poa::POANodeIndex<u32>>) };
            graph.add_alignment_with_weights(&format!("seq{i}"), seq, alignment, w).unwrap();
        }
        println!("after seq{}: graph empty={}", i, graph.is_empty());
    }

    // Dump FASTA buffer to disk for inspection.
    let mut buf: Vec<u8> = Vec::new();
    poa_graph_to_fasta(&graph, &mut buf).unwrap();
    std::fs::write("/tmp/poa_debug.fa", &buf).unwrap();
    println!("wrote {} bytes to /tmp/poa_debug.fa", buf.len());

    // Re-parse.
    let mut rows: Vec<(String, Vec<u8>)> = Vec::new();
    let mut cur_name = String::new();
    let mut cur_seq: Vec<u8> = Vec::new();
    for line in buf.split(|&b| b == b'\n') {
        if line.is_empty() { continue; }
        if line[0] == b'>' {
            if !cur_name.is_empty() {
                rows.push((std::mem::take(&mut cur_name), std::mem::take(&mut cur_seq)));
            }
            cur_name = String::from_utf8_lossy(&line[1..]).to_string();
        } else {
            cur_seq.extend_from_slice(line);
        }
    }
    if !cur_name.is_empty() { rows.push((cur_name, cur_seq)); }
    for (n, s) in &rows {
        let nongap = s.iter().filter(|&&b| b != b'-').count();
        println!("row {}: total_len={} non_gap={} (head: {:?})",
            n, s.len(), nongap, &s[..s.len().min(60)]);
    }
    Ok(())
}
