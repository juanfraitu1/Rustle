//! Intermediary anchor store shared by stages 3 and 4.
//!
//! `refine` computes the shared-k-mer anchor pairs of every block, and
//! `decompose` needs exactly the same lists. Persisting them avoids a
//! full re-anchoring (the dominant cost of both stages) whenever
//! decomposition parameters change.
//!
//! Format v1 (gzip-compressed, little-endian):
//!
//! ```text
//! [8]   magic "SDEFANC1"
//! [u32] k
//! [u64] n_blocks
//! per block:
//!   [u8]  status          (1 = ok, 0 = weak, 2 = self)
//!   [u32] chr_a len + bytes
//!   [u64] a0, [u64] a1
//!   [u32] chr_b len + bytes
//!   [u64] b0, [u64] b1
//!   [u64] n_anchors
//!   n_anchors x ([u32] pb_delta, [u32] pa_delta_wrapping)
//! ```
//!
//! Anchors are sorted by (pb, pa); deltas keep 16-mer lists in the tens
//! of MB even for repeat-saturated chromosomes.

use std::error::Error;
use std::io::{Read, Write};
use std::path::Path;

use flate2::read::GzDecoder;
use flate2::write::GzEncoder;
use flate2::Compression;

const MAGIC: &[u8; 8] = b"SDEFANC1";

/// Block metadata plus its shared-k-mer anchors, in original (pa, pb)
/// coordinates.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct BlockAnchors {
    pub chr_a: String,
    pub a0: u64,
    pub a1: u64,
    pub chr_b: String,
    pub b0: u64,
    pub b1: u64,
    /// 1 = ok, 0 = weak, 2 = self (same encoding as the refine TSV).
    pub status: u8,
    pub anchors: Vec<(u32, u32)>,
}

fn w_u32(buf: &mut Vec<u8>, v: u32) {
    buf.extend_from_slice(&v.to_le_bytes());
}
fn w_u64(buf: &mut Vec<u8>, v: u64) {
    buf.extend_from_slice(&v.to_le_bytes());
}

fn r_u32(buf: &[u8], pos: &mut usize) -> Result<u32, Box<dyn Error>> {
    if *pos + 4 > buf.len() {
        return Err("anchors file truncated (u32)".into());
    }
    let v = u32::from_le_bytes(buf[*pos..*pos + 4].try_into().unwrap());
    *pos += 4;
    Ok(v)
}
fn r_u64(buf: &[u8], pos: &mut usize) -> Result<u64, Box<dyn Error>> {
    if *pos + 8 > buf.len() {
        return Err("anchors file truncated (u64)".into());
    }
    let v = u64::from_le_bytes(buf[*pos..*pos + 8].try_into().unwrap());
    *pos += 8;
    Ok(v)
}

pub fn write(path: &Path, k: usize, blocks: &[BlockAnchors]) -> Result<(), Box<dyn Error>> {
    let f = std::fs::File::create(path)?;
    let mut enc = GzEncoder::new(std::io::BufWriter::new(f), Compression::default());

    let mut head = Vec::with_capacity(24);
    head.extend_from_slice(MAGIC);
    w_u32(&mut head, k as u32);
    w_u64(&mut head, blocks.len() as u64);
    enc.write_all(&head)?;

    for b in blocks {
        let mut rec = Vec::with_capacity(64 + b.anchors.len() * 8);
        rec.push(b.status);
        w_u32(&mut rec, b.chr_a.len() as u32);
        rec.extend_from_slice(b.chr_a.as_bytes());
        w_u64(&mut rec, b.a0);
        w_u64(&mut rec, b.a1);
        w_u32(&mut rec, b.chr_b.len() as u32);
        rec.extend_from_slice(b.chr_b.as_bytes());
        w_u64(&mut rec, b.b0);
        w_u64(&mut rec, b.b1);

        let mut sorted = b.anchors.clone();
        sorted.sort_unstable();
        w_u64(&mut rec, sorted.len() as u64);
        let (mut pb_prev, mut pa_prev) = (0u32, 0u32);
        for &(pa, pb) in &sorted {
            w_u32(&mut rec, pb.wrapping_sub(pb_prev));
            w_u32(&mut rec, pa.wrapping_sub(pa_prev));
            pb_prev = pb;
            pa_prev = pa;
        }
        enc.write_all(&rec)?;
    }
    enc.finish()?;
    Ok(())
}

pub fn read(path: &Path) -> Result<(usize, Vec<BlockAnchors>), Box<dyn Error>> {
    let f = std::fs::File::open(path)?;
    let mut dec = GzDecoder::new(f);
    let mut buf = Vec::new();
    dec.read_to_end(&mut buf)?;
    let mut pos = 0usize;

    if buf.len() < 20 || &buf[0..8] != MAGIC {
        return Err("not an anchors file (bad magic)".into());
    }
    pos += 8;
    let k = r_u32(&buf, &mut pos)? as usize;
    let n = r_u64(&buf, &mut pos)? as usize;
    let mut out = Vec::with_capacity(n);
    for _ in 0..n {
        if pos >= buf.len() {
            return Err("anchors file truncated (record)".into());
        }
        let status = buf[pos];
        pos += 1;
        let la = r_u32(&buf, &mut pos)? as usize;
        if pos + la > buf.len() {
            return Err("anchors file truncated (chr_a)".into());
        }
        let chr_a = String::from_utf8_lossy(&buf[pos..pos + la]).into_owned();
        pos += la;
        let a0 = r_u64(&buf, &mut pos)?;
        let a1 = r_u64(&buf, &mut pos)?;
        let lb = r_u32(&buf, &mut pos)? as usize;
        if pos + lb > buf.len() {
            return Err("anchors file truncated (chr_b)".into());
        }
        let chr_b = String::from_utf8_lossy(&buf[pos..pos + lb]).into_owned();
        pos += lb;
        let b0 = r_u64(&buf, &mut pos)?;
        let b1 = r_u64(&buf, &mut pos)?;
        let na = r_u64(&buf, &mut pos)? as usize;
        let mut anchors = Vec::with_capacity(na);
        let (mut pb, mut pa) = (0u32, 0u32);
        for _ in 0..na {
            pb = pb.wrapping_add(r_u32(&buf, &mut pos)?);
            pa = pa.wrapping_add(r_u32(&buf, &mut pos)?);
            anchors.push((pa, pb));
        }
        out.push(BlockAnchors {
            chr_a,
            a0,
            a1,
            chr_b,
            b0,
            b1,
            status,
            anchors,
        });
    }
    Ok((k, out))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn roundtrip() {
        let blocks = vec![
            BlockAnchors {
                chr_a: "chr1 test".into(),
                a0: 100,
                a1: 200,
                chr_b: "chr1 test".into(),
                b0: 500,
                b1: 700,
                status: 1,
                anchors: vec![(3, 10), (1, 5), (7, 10), (2, 40), (9, 12)],
            },
            BlockAnchors {
                chr_a: "c2".into(),
                a0: 0,
                a1: 50,
                chr_b: "c3".into(),
                b0: 9,
                b1: 99,
                status: 0,
                anchors: vec![],
            },
        ];
        let path = std::env::temp_dir().join(format!(
            "sedef-anchors-test-{}-{}.bin",
            std::process::id(),
            "x"
        ));
        write(&path, 16, &blocks).unwrap();
        let (k, back) = read(&path).unwrap();
        std::fs::remove_file(&path).ok();
        assert_eq!(k, 16);
        assert_eq!(blocks.len(), back.len());
        for (a, b) in blocks.iter().zip(back.iter()) {
            let mut sa = a.anchors.clone();
            sa.sort_unstable();
            assert_eq!(sa, b.anchors);
            assert_eq!(a.chr_a, b.chr_a);
            assert_eq!(a.status, b.status);
        }
    }
}
