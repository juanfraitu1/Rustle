#!/usr/bin/env python3
"""gen_o2_assemble_fixture.py -- golden fixture for the Rust port of
`o2_vg_visualization.materialize_family`'s ASSEMBLY core (assemble_vg).

For each fixture family it runs the REAL Python recompute path
(`materialize_family(fam, use_cache=False)`) against a PRIVATE temp cache (so the shipped
`/home/juanfra/winloci_scratch/o2vg` cache is NOT touched), then recovers the exact
assembly INPUTS from the temp working dir the run leaves behind:

  regions (dedup+sorted-by-length-desc), backbone (bc,bs,be), refseq (uppercased),
  names, copy_cols = column_alleles(copies.bam, copy0), read_cols = column_alleles(
  reads.bam, copy0), best_names (= sorted(best), recovered from the pre-abstain reads
  order), gene_lbl (subset for the regions).

copy_cols / read_cols are TRIMMED to the CANDIDATE PSV positions only -- positions where
every copy is defined AND >=2 distinct alleles. This is a faithful reduction: the Python
assembly consults `copy_cols`/`read_cols` at NO other positions (the bubble filter
`len(hap)==len(names) and set>=2 and read-support>=MIN` short-circuits before touching
`read_cols` off-candidate, and a dropped `copy_cols` entry only ever turns a would-be
non-candidate into a non-candidate). Rust `assemble_vg` on the trimmed maps therefore
produces byte-identical output while keeping the fixture small.

The fixture stores BOTH the pre-abstain VG (the cache dump) and the post-abstain VG (the
`materialize_family` return), so the Rust test checks the assembly AND the recombinant-
abstain merge.

Run: PYTHONHASHSEED=0 /home/juanfra/miniforge3/bin/python bench/gen_o2_assemble_fixture.py
"""
import json
import os
import sys
import tempfile

os.environ.setdefault("PYTHONHASHSEED", "0")
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import psv_graph_genomewide as pg   # noqa: E402
import o2_vg_visualization as o2    # noqa: E402
import pysam                        # noqa: E402

FIXTURE_FAMS = [39, 96, 90, 11]     # resolvable(+68 flips) / no_psv floor(tier3) / small(+flips) / partial(tier3,K<n)
OUT = os.path.join(HERE, "..", "src", "rustle", "vg_family", "testdata",
                   "o2_assemble_fixture.json")


def candidate_positions(refseq, names, copy_cols):
    """Positions where every copy is defined AND >=2 distinct alleles (read_cols/copy_cols
    are consulted ONLY here in the assembly)."""
    cand = []
    for p in range(len(refseq)):
        hap = {names[0]: refseq[p]}
        for nm in names[1:]:
            if nm in copy_cols.get(p, {}):
                hap[nm] = copy_cols[p][nm]
        if len(hap) == len(names) and len(set(hap.values())) >= 2:
            cand.append(p)
    return cand


def cols_to_json(cols, positions):
    """{pos -> {name -> base}} restricted to `positions`, as {str(pos): {name: base}}."""
    out = {}
    for p in positions:
        inner = cols.get(p, {})
        out[str(p)] = {nm: b for nm, b in inner.items()}
    return out


def main():
    genome = pysam.FastaFile(pg.GENOME)
    contig_len = dict(zip(genome.references, genome.lengths))
    gene_lbl = o2.load_gene_labels()
    fams = o2.load_families()

    fam_entries = []
    real_cache = o2.CACHE
    for fam in FIXTURE_FAMS:
        with tempfile.TemporaryDirectory() as td:
            o2.CACHE = td
            post = o2.materialize_family(fam, use_cache=False)
            pre = json.load(open(os.path.join(td, f"fam{fam}.json")))
            tmp = os.path.join(td, f"tmp_fam{fam}")

            regions = pg.dedup_copies(fams[fam])
            regions = sorted(regions, key=lambda r: -(r[2] - r[1]))
            names = [f"copy{i}" for i in range(len(regions))]
            bc, bs, be, _ = regions[0]
            bs, be = max(0, bs), min(contig_len.get(bc, be), be)
            refseq = genome.fetch(bc, bs, be).upper()

            copy_cols = pg.column_alleles(f"{tmp}/copies.bam", names[0])
            read_cols = pg.column_alleles(f"{tmp}/reads.bam", names[0])
            cand = candidate_positions(refseq, names, copy_cols)

            best_names = [r["name"] for r in pre["reads"]]   # sorted(best) order
            gl = {}
            for (c, s, e, nr) in regions:
                if (c, s) in gene_lbl:
                    gl[f"{c}\t{s}"] = gene_lbl[(c, s)]

            fam_entries.append(dict(
                fam_id=fam,
                regions=[[c, s, e, nr] for (c, s, e, nr) in regions],
                backbone=[bc, bs, be],
                refseq=refseq,
                names=names,
                copy_cols=cols_to_json(copy_cols, cand),
                read_cols=cols_to_json(read_cols, cand),
                best_names=best_names,
                gene_lbl=gl,
                n_candidate_positions=len(cand),
                n_bubbles=pre["n_bubbles"],
                n_flipped=sum(1 for r in post["reads"] if r.get("gate_status_pre_abstain")),
                vg_pre=pre,
                vg_post=post,
            ))
            o2.CACHE = real_cache
        print(f"  fam{fam}: {len(regions)} copies, {len(refseq)} bp backbone, "
              f"{len(cand)} candidate cols, bubbles={fam_entries[-1]['n_bubbles']}, "
              f"reads={len(best_names)}, flipped={fam_entries[-1]['n_flipped']}",
              file=sys.stderr)

    fixture = dict(
        note="assemble_vg golden: inputs (columns trimmed to candidate positions) + "
             "pre/post-abstain VG for each family",
        min_reads_psv=pg.MIN_READS_PSV,
        families=fam_entries,
    )
    os.makedirs(os.path.dirname(OUT), exist_ok=True)
    with open(OUT, "w") as f:
        json.dump(fixture, f)
    sz = os.path.getsize(OUT)
    print(f"wrote {OUT}  ({sz/1e6:.2f} MB, {len(fam_entries)} families)", file=sys.stderr)


if __name__ == "__main__":
    main()
