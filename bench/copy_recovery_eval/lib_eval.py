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


def head_to_head(rustle_rec, stringtie_rec, authentic, family_of):
    """rustle_rec / stringtie_rec: {ref_tx: {'fsm','ism'}}. authentic: {ref_tx: bool}
    (only meaningful for rustle recoveries). family_of: {ref_tx: family_id}.
    Returns the headline split as sorted lists + counts."""
    def fsm_set(rec):
        return {tx for tx, v in rec.items() if v["fsm"]}
    r_fsm = fsm_set(rustle_rec)
    s_fsm = fsm_set(stringtie_rec)
    rustle_only = r_fsm - s_fsm
    auth = sorted(tx for tx in rustle_only if authentic.get(tx, False))
    phantom = sorted(tx for tx in rustle_only if not authentic.get(tx, False))
    return {
        "rustle_only_fsm_authentic": auth,
        "rustle_only_fsm_phantom": phantom,
        "n_win": len(auth),
        "n_phantom": len(phantom),
        "families_won": sorted({family_of.get(tx, "NA") for tx in auth}),
    }
