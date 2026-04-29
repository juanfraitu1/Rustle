use rustle::vg_hmm::profile::{ProfileHmm, poa_msa};

#[test]
fn empty_profile_hmm_constructs() {
    let p = ProfileHmm::empty(0);
    assert_eq!(p.n_columns, 0);
    assert!(p.match_emit.is_empty());
}

#[test]
fn singleton_class_profile_has_one_hot_match_emissions() {
    let p = ProfileHmm::from_singleton(b"ACGT");
    assert_eq!(p.n_columns, 4);
    // P(A | col 0) > 0.9 with smoothing floor.
    assert!(p.match_emit[0][0] > 0.9_f64.ln(), "match emit[col 0]['A'] = {}", p.match_emit[0][0]);
}

#[test]
fn poa_msa_returns_equal_length_rows() {
    let inputs: Vec<Vec<u8>> = vec![
        b"ACGTACGT".to_vec(),
        b"ACGAACGT".to_vec(),
        b"ACGTACGT".to_vec(),
    ];
    let rows = poa_msa(&inputs).expect("POA failed");
    let l = rows[0].len();
    assert!(rows.iter().all(|r| r.len() == l), "rows must be equal length, got {:?}", rows.iter().map(|r| r.len()).collect::<Vec<_>>());
}

#[test]
fn fully_conserved_column_has_peaked_emission() {
    let aligned: Vec<Vec<u8>> = vec![b"ACGT".to_vec(), b"ACGT".to_vec(), b"ACGT".to_vec()];
    let p = ProfileHmm::from_msa(&aligned).expect("fit failed");
    // Column 0 is fully A; with n_seq=3 and mixing alpha=0.25, log P(A) ≈ −0.21.
    let log_pa = p.match_emit[0][0];
    assert!(log_pa > -0.30, "expected fully conserved log P(A) > -0.30, got {}", log_pa);
    // Conserved A column should dominate C column by more than 1 nat.
    let log_pc = p.match_emit[0][1];
    assert!(log_pa > log_pc + 1.0, "conserved A column should beat C column by > 1 nat: log_pa={} log_pc={}", log_pa, log_pc);
}

#[test]
fn divergent_column_is_flatter_than_conserved() {
    let aligned: Vec<Vec<u8>> = vec![b"AC".to_vec(), b"AG".to_vec(), b"AT".to_vec()];
    let p = ProfileHmm::from_msa(&aligned).expect("fit failed");
    // Col 0: fully A; col 1: 1/3 each of C/G/T.
    let log_p_a_col0 = p.match_emit[0][0];
    let log_p_c_col1 = p.match_emit[1][1];
    assert!(log_p_a_col0 > log_p_c_col1, "conserved must beat divergent: {} vs {}", log_p_a_col0, log_p_c_col1);
}

#[test]
fn no_gaps_msa_makes_mm_dominant() {
    let rows: Vec<Vec<u8>> = vec![b"ACGT".to_vec(), b"ACGT".to_vec()];
    let p = ProfileHmm::from_msa(&rows).unwrap();
    // M->M should be the dominant transition out of every column.
    for tr in &p.trans {
        assert!(tr.mm > tr.mi && tr.mm > tr.md, "expected M->M dominant, got {:?}", tr);
    }
}

#[test]
fn gappy_msa_increases_md_or_mi() {
    let rows: Vec<Vec<u8>> = vec![b"AC-GT".to_vec(), b"ACGGT".to_vec(), b"AC-GT".to_vec()];
    let p = ProfileHmm::from_msa(&rows).unwrap();
    // Transition out of column 1 (after the C) should show non-trivial M->D
    // because two of three rows have a deletion at column 2.
    assert!(p.trans[1].md > (1e-3_f64).ln(), "expected M->D > 0 at col 1, got {}", p.trans[1].md);
}
