import importlib.util
import os
import random
import subprocess

HERE = os.path.dirname(os.path.abspath(__file__))


def _load(name):
    spec = importlib.util.spec_from_file_location(f"o3_graph_{name}", os.path.join(HERE, f"{name}.py"))
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


L = _load("localref")


def test_regions_are_padded_clamped_and_merged():
    hits = [("c1", 100, 200), ("c1", 250, 300), ("c1", 5000, 5100), ("c2", 10, 20)]
    out = L.merge_regions(hits, pad=100, lengths={"c1": 6000, "c2": 50})
    assert out == [("c1", 0, 400), ("c1", 4900, 5200), ("c2", 0, 50)]


def test_contig_names_round_trip():
    assert L.parse_contig(L.contig_name("chr1", 10, 20)) == ("chr1", 10, 20)


def test_a_record_and_its_supplementary_tag_are_lifted_to_genome_coordinates():
    line = ("r\t0\tchr2:1000-5000\t11\t60\t50M\t*\t0\t0\t" + "A" * 50 + "\t*\tNM:i:0\t"
            "SA:Z:chr2:1000-5000,201,+,20S30M,60,0;chr1:0-900,5,-,30M20S,3,1;\n")
    f = L.lift_line(line).rstrip("\n").split("\t")
    assert f[2] == "chr2" and f[3] == "1011" and f[11] == "NM:i:0"
    assert f[12] == "SA:Z:chr2,1201,+,20S30M,60,0;chr1,5,-,30M20S,3,1;"


def test_an_unmapped_record_is_left_unmapped():
    line = "u\t4\t*\t0\t0\t*\t*\t0\t0\tACGT\t*\n"
    assert L.lift_line(line) == line


def _mut(s, rate, rng):
    return "".join(rng.choice([b for b in "ACGT" if b != c]) if rng.random() < rate else c for c in s)


def test_reads_mapped_to_the_local_reference_land_where_whole_genome_mapping_puts_them(tmp_path):
    import pysam
    rng = random.Random(7)
    rand = lambda n: "".join(rng.choice("ACGT") for _ in range(n))  # noqa: E731
    gene = rand(2000)
    chroms = {"chr1": rand(5000) + gene + rand(13000),                                   # the copy at 5000-7000
              "chr2": rand(10000) + _mut(gene, 0.02, rng) + rand(18000),                 # a paralog at 10000-12000
              "chr3": rand(30000)}                                                       # unrelated
    g = tmp_path / "genome.fa"
    g.write_text("".join(f">{c}\n{s}\n" for c, s in chroms.items()))
    pysam.faidx(str(g))
    regions = L.closure_regions({"t0": gene}, str(g), pad=1000, lengths={c: len(s) for c, s in chroms.items()})
    assert [r[0] for r in regions] == ["chr1", "chr2"]                                    # both homologous loci, not chr3
    local = tmp_path / "local"
    L.write_local(str(g), regions, str(local))
    fq = tmp_path / "reads.fq"
    with open(fq, "w") as fh:
        for i in range(20):
            src = chroms["chr1"] if i % 2 == 0 else chroms["chr2"]
            lo = (5000 if i % 2 == 0 else 10000) + rng.randint(0, 400)
            s = _mut(src[lo:lo + 1200], 0.002, rng)
            fh.write(f"@r{i}\n{s}\n+\n{'~' * len(s)}\n")
    whole = tmp_path / "whole.bam"
    subprocess.run(f"{L.MM2} -t 1 {g} {fq} 2>/dev/null | samtools sort -o {whole} - && samtools index {whole}",
                   shell=True, check=True, executable="/bin/bash")
    lifted = tmp_path / "lifted.bam"
    L.map_reads(str(fq), str(local), str(g) + ".fai", str(lifted), threads=1)

    def prim(p):
        with pysam.AlignmentFile(str(p)) as b:
            return {r.query_name: (r.reference_name, r.reference_start, r.cigarstring, r.mapping_quality)
                    for r in b if not (r.is_secondary or r.is_supplementary)}
    a, b = prim(whole), prim(lifted)
    assert len(a) == 20 and a == b
