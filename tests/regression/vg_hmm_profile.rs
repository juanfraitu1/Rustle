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
    // Column 0 is fully A; expect P(A) close to 1.
    let log_pa = p.match_emit[0][0];
    assert!(log_pa > (-0.05_f64), "expected log P(A | col 0) ≈ 0, got {}", log_pa);
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
