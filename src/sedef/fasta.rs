//! Minimal streaming FASTA reader that keeps full header lines and
//! transparently handles gzip-compressed input. needletail would be
//! faster, but it only exposes the first whitespace-delimited token of
//! a header; we need the description (e.g. "chromosome 1" in NCBI
//! files) for `--chr` filtering.

use std::error::Error;
use std::fs::File;
use std::io::{BufRead, BufReader};
use std::path::Path;

pub fn open_maybe_gz(path: &Path) -> std::io::Result<Box<dyn BufRead>> {
    let f = File::open(path)?;
    let inner = BufReader::with_capacity(1 << 20, f);
    if path.extension().is_some_and(|e| e == "gz") {
        Ok(Box::new(BufReader::new(flate2::read::MultiGzDecoder::new(
            inner,
        ))))
    } else {
        Ok(Box::new(inner))
    }
}

/// Reads all records (or those whose full header contains `chr_filter`,
/// case-insensitively). Headers are returned without the leading '>'.
pub fn read_fasta(
    path: &Path,
    chr_filter: Option<&str>,
) -> Result<Vec<(String, Vec<u8>)>, Box<dyn Error>> {
    let mut reader = open_maybe_gz(path)?;
    let needle = chr_filter.map(|s| s.to_lowercase());

    let mut seqs = Vec::new();
    let mut header: Option<String> = None;
    let mut seq: Vec<u8> = Vec::new();
    let mut line = Vec::new();

    let flush = |header: Option<String>, seq: &mut Vec<u8>, seqs: &mut Vec<(String, Vec<u8>)>| {
        if let Some(h) = header {
            let keep = match &needle {
                Some(n) => h.to_lowercase().contains(n.as_str()),
                None => true,
            };
            if keep {
                seqs.push((h, std::mem::take(seq)));
            } else {
                seq.clear();
            }
        }
    };

    loop {
        line.clear();
        if reader.read_until(b'\n', &mut line)? == 0 {
            break;
        }
        if line[0] == b'>' {
            flush(header.take(), &mut seq, &mut seqs);
            header = Some(
                String::from_utf8_lossy(&line[1..])
                    .trim_end()
                    .to_string(),
            );
        } else {
            seq.extend(line.iter().copied().filter(|b| !b.is_ascii_whitespace()));
        }
    }
    flush(header, &mut seq, &mut seqs);
    Ok(seqs)
}
