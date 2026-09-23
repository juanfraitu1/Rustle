//! §6zb reader probe: how fast can noodles iterate a BAM's records when NOTHING is extracted?
//! Two paths: (a) indexed `query` over a region (what `reads_in_region` does, one `Record` allocated per
//! iteration); (b) index-seek to the region's first chunk, then sequential `read_record` into ONE reused
//! `Record` until the reference id changes. Prints record counts and wall seconds.
use noodles_csi::BinningIndex as _;
use std::io::BufReader;
use std::time::Instant;

fn main() -> anyhow::Result<()> {
    let args: Vec<String> = std::env::args().collect();
    let (bam, chrom) = (&args[1], &args[2]);
    let mode = args.get(3).map(|s| s.as_str()).unwrap_or("both");

    if mode == "both" || mode == "query" {
        let t = Instant::now();
        let file = std::fs::File::open(bam)?;
        let mut reader = noodles_bam::io::Reader::from(noodles_bgzf::Reader::new(BufReader::with_capacity(1 << 20, file)));
        let header = reader.read_header()?;
        let index = noodles_bam::bai::read(format!("{bam}.bai"))?;
        let len = header.reference_sequences().get(chrom.as_bytes()).map(|r| usize::from(r.length())).unwrap();
        let region: noodles_core::Region = format!("{chrom}:1-{len}").parse()?;
        let mut n = 0usize;
        for r in reader.query(&header, &index, &region)? {
            let _ = r?;
            n += 1;
        }
        eprintln!("query      : {n} records in {:.1} s", t.elapsed().as_secs_f64());
    }
    if mode == "recordbuf" {
        // the query path PLUS the RecordBuf conversion `reads_in_region` performs on every record
        let t = Instant::now();
        let file = std::fs::File::open(bam)?;
        let mut reader = noodles_bam::io::Reader::from(noodles_bgzf::Reader::new(BufReader::with_capacity(1 << 20, file)));
        let header = reader.read_header()?;
        let index = noodles_bam::bai::read(format!("{bam}.bai"))?;
        let len = header.reference_sequences().get(chrom.as_bytes()).map(|r| usize::from(r.length())).unwrap();
        let region: noodles_core::Region = format!("{chrom}:1-{len}").parse()?;
        let mut n = 0usize;
        let mut cig = 0usize;
        for r in reader.query(&header, &index, &region)? {
            let r = r?;
            let rb = noodles_sam::alignment::RecordBuf::try_from_alignment_record(&header, &r)?;
            cig += rb.cigar().as_ref().len();
            n += 1;
        }
        eprintln!("query+RecordBuf: {n} records ({cig} cigar ops) in {:.1} s", t.elapsed().as_secs_f64());
    }
    if mode == "mt" {
        // the query path through the MultithreadedReader (4 workers) `reads_in_region_indexed` uses
        let t = Instant::now();
        let file = std::fs::File::open(bam)?;
        let worker = std::num::NonZeroUsize::new(4).unwrap();
        let bgzf = noodles_bgzf::MultithreadedReader::with_worker_count(worker, BufReader::with_capacity(1 << 20, file));
        let mut reader = noodles_bam::io::Reader::from(bgzf);
        let header = reader.read_header()?;
        let index = noodles_bam::bai::read(format!("{bam}.bai"))?;
        let len = header.reference_sequences().get(chrom.as_bytes()).map(|r| usize::from(r.length())).unwrap();
        let region: noodles_core::Region = format!("{chrom}:1-{len}").parse()?;
        let mut n = 0usize;
        for r in reader.query(&header, &index, &region)? {
            let _ = r?;
            n += 1;
        }
        eprintln!("query via MultithreadedReader(4): {n} records in {:.1} s", t.elapsed().as_secs_f64());
    }
    if mode == "both" || mode == "seq" {
        let t = Instant::now();
        let file = std::fs::File::open(bam)?;
        let mut reader = noodles_bam::io::Reader::from(noodles_bgzf::Reader::new(BufReader::with_capacity(1 << 20, file)));
        let header = reader.read_header()?;
        let index = noodles_bam::bai::read(format!("{bam}.bai"))?;
        let tid = header.reference_sequences().get_index_of(chrom.as_bytes()).unwrap();
        // first chunk start of this reference: the minimum chunk begin over its bins
        let rs = &index.reference_sequences()[tid];
        let start = rs.bins().values().flat_map(|b| b.chunks().iter().map(|c| c.start())).min().unwrap();
        reader.get_mut().seek(start)?;
        let mut rec = noodles_bam::Record::default();
        let mut n = 0usize;
        while reader.read_record(&mut rec)? > 0 {
            match rec.reference_sequence_id() {
                Some(Ok(id)) if id == tid => n += 1,
                Some(Ok(id)) if id > tid => break,
                _ => {}
            }
        }
        eprintln!("seq+reuse  : {n} records in {:.1} s", t.elapsed().as_secs_f64());
    }
    Ok(())
}
