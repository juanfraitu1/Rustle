//! Quick utility — break down read counts in each family region by category.
fn main() -> anyhow::Result<()> {
    let path = std::env::args().nth(1)
        .unwrap_or_else(|| "/scratch/jxi21/Assembler/GGO_19.bam".to_string());
    let mut reader = noodles_bam::io::indexed_reader::Builder::default()
        .build_from_path(&path)?;
    let header = reader.read_header()?;
    println!("BAM: {}", path);
    let regions: &[(&str, &str)] = &[
        ("GOLGA6L7",        "NC_073243.2:104789647-104877901"),
        ("  COPY0",         "NC_073243.2:104789647-104796276"),
        ("  COPY1",         "NC_073243.2:104830535-104837094"),
        ("  COPY2",         "NC_073243.2:104871355-104877901"),
        ("KRAB-ZNF",        "NC_073244.2:65299948-65342879"),
        ("AMY2A/B",         "NC_073224.2:136262907-136334077"),
        ("TBC1D3 LOC1",     "NC_073228.2:44678473-44689131"),
        ("OR_chr1 LOC1",    "NC_073224.2:6032726-6033683"),
        ("NBPF1",           "NC_073224.2:124449205-124474484"),
        ("GOLGA8 LOC1",     "NC_073240.2:31359988-31374377"),
    ];
    println!("  {:<18} {:<40} {:>9} {:>9} {:>9} {:>9}",
        "name", "region", "primary", "secondary", "suppl.", "unmapped");
    for (name, region_str) in regions {
        let region = match region_str.parse() { Ok(r) => r, Err(_) => continue };
        let q = match reader.query(&header, &region) { Ok(it) => it, Err(_) => continue };
        let (mut p, mut s, mut x, mut u) = (0u64, 0u64, 0u64, 0u64);
        for r in q {
            if let Ok(rec) = r {
                let f = rec.flags();
                if f.is_unmapped() { u += 1; }
                else if f.is_secondary() { s += 1; }
                else if f.is_supplementary() { x += 1; }
                else { p += 1; }
            }
        }
        println!("  {:<18} {:<40} {:>9} {:>9} {:>9} {:>9}", name, region_str, p, s, x, u);
    }
    Ok(())
}
