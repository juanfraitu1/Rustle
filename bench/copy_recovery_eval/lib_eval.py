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
