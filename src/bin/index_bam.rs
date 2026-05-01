//! Build a .bai index for a coordinate-sorted BAM.
//!
//! Usage: index_bam <in.bam> [out.bam.bai]

use anyhow::{anyhow, Context, Result};
use noodles_bam as bam;
use noodles_bgzf as bgzf;
use noodles_core::Position;
use noodles_csi::binning_index::{
    index::reference_sequence::bin::Chunk,
    Indexer,
};
use noodles_sam::{
    self as sam,
    alignment::Record as _,
    header::record::value::map::header::{sort_order::COORDINATE, tag::SORT_ORDER},
};
use std::env;
use std::fs::File;
use std::path::PathBuf;

fn is_coordinate_sorted(header: &sam::Header) -> bool {
    header
        .header()
        .and_then(|hdr| hdr.other_fields().get(&SORT_ORDER))
        .map(|sort_order| sort_order == COORDINATE)
        .unwrap_or_default()
}

fn alignment_context(
    record: &bam::Record,
) -> std::io::Result<(Option<usize>, Option<Position>, Option<Position>)> {
    Ok((
        record.reference_sequence_id().transpose()?,
        record.alignment_start().transpose()?,
        record.alignment_end().transpose()?,
    ))
}

fn main() -> Result<()> {
    let args: Vec<String> = env::args().collect();
    if args.len() < 2 {
        eprintln!("usage: index_bam <in.bam> [out.bam.bai]");
        std::process::exit(2);
    }
    let in_path = PathBuf::from(&args[1]);
    let out_path = if args.len() >= 3 {
        PathBuf::from(&args[2])
    } else {
        let mut p = in_path.clone();
        let s = p.file_name().unwrap().to_string_lossy().to_string() + ".bai";
        p.set_file_name(s);
        p
    };
    eprintln!("[index_bam] in={} out={}", in_path.display(), out_path.display());

    let mut reader = bam::io::reader::Builder::default()
        .build_from_path(&in_path)
        .with_context(|| format!("opening {}", in_path.display()))?;
    let header = reader.read_header()?;

    if !is_coordinate_sorted(&header) {
        return Err(anyhow!(
            "BAM is not marked as coordinate-sorted (SO:coordinate). Refusing to index."
        ));
    }

    let mut record = bam::Record::default();
    let mut builder = Indexer::default();
    let mut start_position = reader.get_ref().virtual_position();
    let mut n = 0u64;
    while reader.read_record(&mut record)? != 0 {
        n += 1;
        let end_position = reader.get_ref().virtual_position();
        let chunk = Chunk::new(start_position, end_position);

        let context = match alignment_context(&record)? {
            (Some(id), Some(start), Some(end)) => {
                let is_mapped = !record.flags().is_unmapped();
                Some((id, start, end, is_mapped))
            }
            _ => None,
        };
        builder.add_record(context, chunk)?;
        start_position = end_position;

        if n % 500000 == 0 {
            eprintln!("[index_bam] indexed {} records", n);
        }
    }

    let index = builder.build(header.reference_sequences().len());
    eprintln!("[index_bam] built index over {} records, {} reference sequences", n, header.reference_sequences().len());

    let mut writer = bam::bai::io::Writer::new(File::create(&out_path)?);
    writer.write_index(&index)?;
    eprintln!("[index_bam] wrote {}", out_path.display());
    Ok(())
}
