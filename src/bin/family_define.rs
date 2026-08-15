//! family_define -- native Rust O1 family-definition DRIVER (interest I / O1).
//!
//! CAPSTONE of the O1 migration: reproduces the shipped RNA multi-copy family catalog
//! `bench/family_rna_refine.tsv` byte-for-byte (default flags -> md5 dca64cbd), replacing
//! the Python driver `bench/family_rna_refine.py`. All the building blocks are ported +
//! byte-parity tested in `rustle::vg_family`; this binary wires the orchestration in
//! `rustle::vg_family::driver` behind the same CLI + env flags as the Python driver.
//!
//! DEFAULT (no flag) = the shipped family definition: core+aln edge gate + repeat-hub gate
//! + gamma=0.20 refinement + recombinant-split gate + multi-repeat-bridge gate + allele
//! demote. Opt out with --legacy; ablate individual gates with --no-*; --high-precision
//! swaps ONLY gamma (0.20 -> 0.40).

use std::io::Write;

use clap::Parser;

use rustle::vg_family::driver::{self, Inputs, Options, GAMMA, HIGH_PRECISION_GAMMA};

#[derive(Parser, Debug)]
#[command(about = "LEGACY parity fixture: reproduces the frozen Python catalog from precomputed TSVs. \
                   The shipped genome-wide catalog is `gw_family_catalog` (exon-sum E_r edge). \
                   Kept for regression parity only; not the live family definition. See rustle_mechanism.html.")]
struct Args {
    /// Output catalog TSV (Python default: bench/family_rna_refine.tsv).
    #[arg(long, default_value = "bench/family_rna_refine.tsv")]
    out: String,

    /// Opt OUT of the RNA-only refinement: write nothing and print the legacy line
    /// (also RUSTLE_RNA_ORACLE=0).
    #[arg(long, default_value_t = false)]
    legacy: bool,

    /// Deprecated no-op (the RNA-only refinement is now the default).
    #[arg(long = "rna-oracle", default_value_t = false)]
    rna_oracle: bool,

    /// Ablation: DISABLE just the repeat-hub gate (min_shared_mult>=20 cut)
    /// (also RUSTLE_NO_REPEAT_GATE=1).
    #[arg(long = "no-repeat-gate", default_value_t = false)]
    no_repeat_gate: bool,

    /// Ablation: DISABLE just the recombinant-split gate (also RUSTLE_NO_SPLIT_RECOMBINANTS=1).
    #[arg(long = "no-split-recombinants", default_value_t = false)]
    no_split_recombinants: bool,

    /// Ablation: DISABLE just the multi-repeat-bridge gate (also RUSTLE_NO_REPEAT_BRIDGE_GATE=1).
    #[arg(long = "no-repeat-bridge-gate", default_value_t = false)]
    no_repeat_bridge_gate: bool,

    /// HIGH-PRECISION operating point: swap ONLY gamma 0.20 -> 0.40
    /// (PRECISION_RECALL_FRONTIER.md). Also RUSTLE_HIGH_PRECISION=1.
    #[arg(long = "high-precision", default_value_t = false)]
    high_precision: bool,

    /// Run + report but do NOT write outputs (used by the self-check).
    #[arg(long = "no-write", default_value_t = false)]
    no_write: bool,

    // ---- input paths ----
    /// Denovo transcript meta TSV. REQUIRED — no hardcoded winloci_scratch fallback.
    #[arg(long)]
    meta: Option<String>,
    /// Annotation intervals TSV. REQUIRED — no hardcoded winloci_scratch fallback.
    #[arg(long)]
    annot: Option<String>,
    #[arg(long)]
    raw_fams: Option<String>,
    #[arg(long)]
    edges: Option<String>,
    #[arg(long)]
    univ_aln: Option<String>,
    #[arg(long)]
    allele: Option<String>,
    #[arg(long)]
    vg_repeat: Option<String>,
    /// De-novo skeletons TSV. REQUIRED — no hardcoded winloci_scratch fallback.
    #[arg(long)]
    skeletons: Option<String>,
    /// Genome FASTA (soft-masked GGO.fasta). REQUIRED — no hardcoded winloci_scratch fallback.
    #[arg(long)]
    genome: Option<String>,
}

fn env_on(k: &str, want: &str) -> bool {
    std::env::var(k).ok().as_deref() == Some(want)
}

fn main() {
    let args = Args::parse();

    // DEFAULT-ON: the RNA-only refinement IS the family definition unless the legacy opt-out.
    if args.legacy || env_on("RUSTLE_RNA_ORACLE", "0") {
        println!(
            "legacy: RNA-only refinement DISABLED \
             (core_recip>=0.13 shipped catalog; run bench/denovo_families.py for the legacy path)"
        );
        return;
    }

    let meta = args
        .meta
        .expect("--meta is required (no hardcoded winloci_scratch fallback)");
    let annot = args
        .annot
        .expect("--annot is required (no hardcoded winloci_scratch fallback)");
    let skeletons = args
        .skeletons
        .expect("--skeletons is required (no hardcoded winloci_scratch fallback)");
    let genome = args
        .genome
        .expect("--genome is required (no hardcoded winloci_scratch fallback)");
    let mut inp = Inputs::new(meta, annot, skeletons, genome);
    if let Some(v) = args.raw_fams {
        inp.raw_fams = v;
    }
    if let Some(v) = args.edges {
        inp.edges = v;
    }
    if let Some(v) = args.univ_aln {
        inp.univ_aln = v;
    }
    if let Some(v) = args.allele {
        inp.allele = v;
    }
    if let Some(v) = args.vg_repeat {
        inp.vg_repeat = v;
    }

    let repeat_gate = !(args.no_repeat_gate || env_on("RUSTLE_NO_REPEAT_GATE", "1"));
    let split_recombinants =
        !(args.no_split_recombinants || env_on("RUSTLE_NO_SPLIT_RECOMBINANTS", "1"));
    let repeat_bridge_gate =
        !(args.no_repeat_bridge_gate || env_on("RUSTLE_NO_REPEAT_BRIDGE_GATE", "1"));
    let high_precision = args.high_precision || env_on("RUSTLE_HIGH_PRECISION", "1");
    let gamma = if high_precision {
        HIGH_PRECISION_GAMMA
    } else {
        GAMMA
    };

    let opt = Options {
        repeat_gate,
        split_recombinants,
        repeat_bridge_gate,
        gamma,
    };

    let built = match driver::build_catalog(&inp, opt) {
        Ok(b) => b,
        Err(e) => {
            eprintln!("family_define: build failed: {e}");
            std::process::exit(1);
        }
    };

    let hp = (gamma - GAMMA).abs() > 1e-9;
    eprintln!(
        "==================== RNA-ONLY FAMILY DEFINITION ({}) ====================",
        if hp {
            format!("HIGH-PRECISION gamma={gamma:.2}")
        } else {
            "default".to_string()
        }
    );
    eprintln!(
        "gates            : repeat-hub {}  recombinant-split {}  repeat-bridge {}  gamma={}",
        on_off(repeat_gate),
        on_off(split_recombinants),
        on_off(repeat_bridge_gate),
        gamma
    );
    eprintln!(
        "DN edges         : within-gene kept={}  cross-gene kept={}  cross-gene CUT={} (repeat-hub={})",
        built.n_dn_within, built.n_dn_cross_kept, built.n_dn_cross_cut, built.n_dn_cross_cut_repeat
    );
    eprintln!(
        "recombinant split: families split={}   multi-repeat-bridge: families cut={}",
        built.n_families_split, built.n_families_repeat_bridge_split
    );
    eprintln!("alleles demoted  : {}", built.n_demoted);
    eprintln!("n_families       : {}", built.catalog.len());

    if !args.no_write {
        let mut buf: Vec<u8> = Vec::new();
        if let Err(e) = driver::write_outputs(&built, &mut buf) {
            eprintln!("family_define: write failed: {e}");
            std::process::exit(1);
        }
        if let Err(e) = std::fs::write(&args.out, &buf) {
            eprintln!("family_define: could not write {}: {e}", args.out);
            std::process::exit(1);
        }
        // md5-free stdout receipt (the byte-parity check hashes the file itself).
        let mut stdout = std::io::stdout();
        let _ = writeln!(stdout, "wrote {} ({} families)", args.out, built.catalog.len());
    }
}

fn on_off(b: bool) -> &'static str {
    if b {
        "ON"
    } else {
        "OFF"
    }
}
