"""bench/crosscopy/tool_table.py: which arms emit the junction of a copy-specific call at the copy that uses it, and at the one that does not."""
import os, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from tool_table import table
from test_calls import write
import calls as K

def sample(tmp_path, name, genome, units, arms, reads=30, first=0, s_hi=0.0, miskeyed=()):
    """`units` LARGE calls (family Fi: U on c1, S on c2, junction at 1000 i + 101 .. 1000 i + 400), copy-specific unless s_hi > 0.05;
    arms = {name: (n at U, n at S)}: the arm's GTF holds the junction of the first n calls at U, and of the first n at S.
    Member names carry the sample's name, as locus spans do; the sites table marks the calls in `miskeyed`."""
    d = tmp_path / name; d.mkdir()
    iv = lambda i: (1000 * i + 101, 1000 * i + 400)
    rows = []
    for i in range(first, first + units):
        rows += [(f"F{i}", f"c1:{name}", "donor", f"C{i}", f"D{i}", "c1", "+", *iv(i), reads // 2, reads, 0.5, 0.5, 0.5, 0),
                 (f"F{i}", f"c2:{name}", "donor", f"C{i}", f"D{i}", "c2", "+", *iv(i), 0, reads, 0.0, s_hi, 0.0, int(not s_hi))]
    psi = write(d, rows)
    (d / "sites.tsv").write_text("family\tanchor\tcol_anchor\tcol_other\tverdict_1\tunit\n" + "".join(
        f"{c.family}\t{c.anchor}\t{c.col_anchor}\t{c.col_other}\t{'MISKEYED' if c.family in miskeyed else 'AGREES'}\t{K.unit_name(c)}\n" for c in K.load_calls(psi)))
    spec = [genome, psi, str(d / "sites.tsv")]
    for arm, (nu, ns) in arms.items():
        ex = lambda ch, i: "".join(f'{ch}\tx\texon\t{s}\t{e}\t.\t+\t.\tgene_id "g{ch}{i}"; transcript_id "t{ch}{i}";\n' for s, e in ((iv(i)[0] - 50, iv(i)[0] - 1), (iv(i)[1] + 1, iv(i)[1] + 50)))
        (d / f"{arm}.gtf").write_text("".join(ex("c1", i) for i in range(first, first + nu)) + "".join(ex("c2", i) for i in range(first, first + ns)))
        spec.append(f"{arm}={d / f'{arm}.gtf'}")
    return ",".join(spec)

def test_shares_and_reading(tmp_path):
    t = table(None, [sample(tmp_path, "a", "G", 40, {"x": (36, 0), "y": (30, 4), "v": (21, 0), "w": (19, 0)})])
    assert t["x"] == (40, 0.9, 0.0) and t["y"] == (40, 0.75, 0.1) and t.reading("y", "x") == "behind" and t.reading("x", "y") == "level"
    assert t.reading("w", "v") == "level"                                        # 19 against 21 of 40: lower by 0.05, not by more
    t = table(None, [sample(tmp_path, "b", "G", 29, {"x": (29, 0), "y": (0, 0)})])
    assert t.reading("y", "x") == "level"                                        # under 30 units

def test_units_once_per_genome_and_without_the_development_sample(tmp_path):
    dev = sample(tmp_path, "dev", "G", 5, {})                                    # units 0-4 of genome G were seen in development
    a = sample(tmp_path, "a", "G", 20, {"x": (20, 0)}, reads=30)                 # units 0-19, x emits all at U
    b = sample(tmp_path, "b", "G", 20, {"x": (0, 0), "z": (20, 20)}, reads=90)   # the same 20 units, more reads: this sample speaks for them
    c = sample(tmp_path, "c", "H", 20, {"x": (20, 0)})                           # another genome: other units, none in development
    t = table("G=" + dev.split(",")[2], [a, b, c])
    assert t["x"] == (35, 20 / 35, 0.0)                                          # 15 units of G (from b: not emitted) + 20 of H
    assert t["z"] == (15, 1.0, 1.0) and t.reading("x", "z") == "level"           # x and z share sample b only: 15 units

def test_the_speaking_sample_decides_copy_specific_and_miskeyed_units_are_out(tmp_path):
    a = sample(tmp_path, "a", "G", 10, {"x": (10, 0)}, reads=30)                 # copy-specific here ...
    b = sample(tmp_path, "b", "G", 10, {"x": (10, 0)}, reads=90, s_hi=0.2)       # ... but not in the sample that speaks for the units
    assert table(None, [a, b])["x"][0] == 0
    c = sample(tmp_path, "c", "H", 10, {"x": (10, 0)}, miskeyed={"F3", "F4"})
    assert table(None, [c])["x"] == (8, 1.0, 0.0)
