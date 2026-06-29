"""Pure scoring logic for the VG copy-recovery protocol (U3, U4, U6).
No file or subprocess I/O in these functions — they take parsed data and
return plain dict rows, so they are unit-testable in isolation."""

def _union_find_families(genes, pairs):
    parent = {g: g for g in genes}
    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x
    for a, b in pairs:
        if a in parent and b in parent:
            ra, rb = find(a), find(b)
            if ra != rb:
                parent[ra] = rb
    fam = {}
    for g in genes:
        fam.setdefault(find(g), []).append(g)
    return fam  # root_gene -> [genes]


def build_universe(gene_tx, paralog_pairs):
    """gene_tx: {gene_id: [transcript_id,...]}. paralog_pairs: [(geneA,geneB),...]
    that already passed identity/coverage thresholds.
    Returns rows: {transcript_id, gene_id, family_id, n_family_copies} for genes
    in families of >=2 genes. family_id is the deterministic min gene_id in the family."""
    families = _union_find_families(list(gene_tx.keys()), paralog_pairs)
    rows = []
    for _root, genes in families.items():
        if len(genes) < 2:
            continue
        family_id = min(genes)
        n_copies = len(genes)
        for g in sorted(genes):
            for tx in gene_tx[g]:
                rows.append({
                    "transcript_id": tx,
                    "gene_id": g,
                    "family_id": family_id,
                    "n_family_copies": n_copies,
                })
    rows.sort(key=lambda r: (r["family_id"], r["gene_id"], r["transcript_id"]))
    return rows


FSM = "full-splice_match"
ISM = "incomplete-splice_match"


def norm_tx_id(s):
    """Normalize a transcript id to its bare NCBI accession so IDs join across
    sources. The annotation's `transcript_id=` attribute is the clean accession
    (e.g. `XM_055380435.2`), but gffread -T (which feeds SQANTI3's reference)
    emits the GFF3 `ID=` form with an `rna-`/`gene-` prefix (e.g.
    `rna-XR_008678095.2`). SQANTI3's `associated_transcript` therefore carries
    that prefix while `universe.tsv` does not — stripping the prefix on both
    sides makes the join work. Idempotent."""
    if s is None:
        return ""
    for pre in ("rna-", "gene-"):
        if s.startswith(pre):
            return s[len(pre):]
    return s


def recovery_set(classification_rows, universe_tx):
    """classification_rows: list of dicts with 'structural_category' and
    'associated_transcript'. universe_tx: set of reference transcript ids (U).
    Returns {ref_transcript: {'fsm': bool, 'ism': bool}} restricted to U.
    IDs are normalized via norm_tx_id on both sides before matching, so a
    prefixed SQANTI3 associated_transcript joins the clean universe accession.
    The returned keys are the normalized (clean) ids."""
    uni = {norm_tx_id(t) for t in universe_tx}
    rec = {}
    for row in classification_rows:
        cat = row.get("structural_category", "")
        ref = norm_tx_id(row.get("associated_transcript", ""))
        if ref not in uni:
            continue
        if cat not in (FSM, ISM):
            continue
        cur = rec.setdefault(ref, {"fsm": False, "ism": False})
        if cat == FSM:
            cur["fsm"] = True
        elif cat == ISM:
            cur["ism"] = True
    return rec


AUTHENTIC = "authentic"
PHANTOM = "phantom"
UNRESOLVABLE = "unresolvable"


def classify_authenticity(n_decisive, n_tied, k_decisive=2, k_tied=1):
    """Three-bucket authenticity from own-copy read evidence.

    - n_decisive: reads whose best alignment is *strictly better* on this copy than
      on any sister copy (or that map uniquely to it) -> the copy genuinely owns them.
    - n_tied: reads that align *equally well* to this copy and a sister -> genuinely
      undecidable (the DAZ inverted-backbone case).

    Buckets:
      authentic    -- n_decisive >= k_decisive (proven recovery).
      unresolvable -- not authentic but n_tied >= k_tied (real-or-not undecidable;
                      do NOT call it a phantom, that would lose true matches).
      phantom      -- neither: whatever coverage exists is sister spillover that the
                      sister owns decisively (the fabrication case, e.g. prim=0 copies).

    A *count* of primary alignments is deliberately NOT used: minimap2's primary flag
    among MAPQ-0 multimappers is a coin-flip, so it both over-credits phantoms (a copy
    that wins the flip) and can starve real copies whose reads all tie.
    """
    if n_decisive >= k_decisive:
        return AUTHENTIC
    if n_tied >= k_tied:
        return UNRESOLVABLE
    return PHANTOM


def _coerce_status(v):
    """Accept either a status string or a legacy bool (True->authentic, False->phantom)."""
    if isinstance(v, bool):
        return AUTHENTIC if v else PHANTOM
    return v


def head_to_head(rustle_rec, stringtie_rec, status, family_of):
    """rustle_rec / stringtie_rec: {ref_tx: {'fsm','ism'}}.
    status: {ref_tx: 'authentic'|'phantom'|'unresolvable'} (legacy bools accepted)
    for rustle recoveries. family_of: {ref_tx: family_id}.
    Returns the headline split (authentic wins / phantom / unresolvable) + counts."""
    def fsm_set(rec):
        return {tx for tx, v in rec.items() if v["fsm"]}
    r_fsm = fsm_set(rustle_rec)
    s_fsm = fsm_set(stringtie_rec)
    rustle_only = r_fsm - s_fsm
    st = {tx: _coerce_status(status.get(tx, PHANTOM)) for tx in rustle_only}
    auth = sorted(tx for tx in rustle_only if st[tx] == AUTHENTIC)
    phantom = sorted(tx for tx in rustle_only if st[tx] == PHANTOM)
    unresolvable = sorted(tx for tx in rustle_only if st[tx] == UNRESOLVABLE)
    return {
        "rustle_only_fsm_authentic": auth,
        "rustle_only_fsm_phantom": phantom,
        "rustle_only_fsm_unresolvable": unresolvable,
        "n_win": len(auth),
        "n_phantom": len(phantom),
        "n_unresolvable": len(unresolvable),
        "families_won": sorted({family_of.get(tx, "NA") for tx in auth}),
    }
