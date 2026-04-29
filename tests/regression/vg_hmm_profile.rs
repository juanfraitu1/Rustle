use rustle::vg_hmm::profile::ProfileHmm;

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
