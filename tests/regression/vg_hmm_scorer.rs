use rustle::vg_hmm::profile::ProfileHmm;
use rustle::vg_hmm::scorer::forward_against_profile;

#[test]
fn perfectly_matching_read_scores_higher_than_random() {
    let p = ProfileHmm::from_singleton(b"ACGTACGT");
    let good = forward_against_profile(&p, b"ACGTACGT");
    let bad  = forward_against_profile(&p, b"TTTTTTTT");
    assert!(good > bad, "good={} bad={}", good, bad);
}
