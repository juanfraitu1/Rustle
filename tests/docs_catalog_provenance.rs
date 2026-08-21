//! Lint: no stale catalog reference may sit in the docs without its provenance.
//!
//! WHY THIS TEST EXISTS. The shipped 494-family / 1,415-copy catalog was built with `refine` on — a
//! default removed on 2026-08-20 — so **no invocation of the current binary reproduces it**; the
//! current default emits 627 families / 2,019 copies. Figures quoted per `/494` or `/1415` are
//! therefore historical, and a reader who meets one without a label will take it for current. That
//! happened: `24.88%` of copies ≤ 2 kb (352/1415) was quoted repeatedly as a live fact when the
//! current catalog reads 21.40% (432/2019), and two provenance labels were written wrong on the same
//! document in two days.
//!
//! A promise not to do it again is not a mechanism. This is the mechanism.
//!
//! A document passes if EITHER it carries the standard `CATALOG PROVENANCE.` banner, OR every line
//! mentioning `/494` or `/1415` carries an inline provenance token. Adding a new stale reference to
//! an un-bannered document fails the suite.

use std::fs;

const BANNER: &str = "CATALOG PROVENANCE.";

/// Inline tokens that make a single line self-describing without the document-level banner.
const INLINE: &[&str] = &[
    "494-family", "494 catalog", "superseded", "Jul-17", "2026-07-17", "HISTORICAL",
    "1415 copies", "1,415 copies", "pre-guard", "guard OFF", "shipped catalog",
];

#[test]
fn no_unlabelled_stale_catalog_references_in_docs() {
    let mut offenders: Vec<String> = Vec::new();
    let dir = std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("docs");
    let mut checked = 0usize;

    for entry in fs::read_dir(&dir).expect("docs/ must exist") {
        let path = entry.expect("readable dir entry").path();
        if path.extension().and_then(|e| e.to_str()) != Some("md") {
            continue;
        }
        let text = match fs::read_to_string(&path) {
            Ok(t) => t,
            Err(_) => continue,
        };
        // A stale reference is `/494` or `/1415` not followed by another digit (so `/4940` is fine).
        let stale = |l: &str| {
            ["/494", "/1415"].iter().any(|pat| {
                l.match_indices(pat).any(|(i, _)| {
                    !l[i + pat.len()..].chars().next().is_some_and(|c| c.is_ascii_digit())
                })
            })
        };
        if !text.lines().any(stale) {
            continue;
        }
        checked += 1;
        if text.contains(BANNER) {
            continue;
        }
        for (n, line) in text.lines().enumerate() {
            if stale(line) && !INLINE.iter().any(|t| line.contains(t)) {
                offenders.push(format!(
                    "{}:{}",
                    path.file_name().and_then(|s| s.to_str()).unwrap_or("?"),
                    n + 1
                ));
            }
        }
    }

    assert!(
        checked > 0,
        "the lint found no /494 or /1415 references at all — it has stopped testing anything, \
         which is worse than failing. Check the stale-reference patterns."
    );
    assert!(
        offenders.is_empty(),
        "stale catalog references without provenance ({} site(s)): {:?}\n\
         The 494-family / 1,415-copy catalog is SUPERSEDED — the current default emits 627 families \
         / 2,019 copies and no invocation reproduces the old one. Either add the standard \
         `CATALOG PROVENANCE.` banner to the document, or label the line inline (e.g. \"494 catalog\", \
         \"superseded\"). Do not simply delete the marker to make this pass.",
        offenders.len(),
        offenders
    );
}
