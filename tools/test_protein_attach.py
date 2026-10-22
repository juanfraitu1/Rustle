"""tools/protein_attach.py: the exon-gate guard (a families PAF written under the gate holds only the gated pairs)."""
import importlib.util
import os
import tempfile

HERE = os.path.dirname(os.path.abspath(__file__))
spec = importlib.util.spec_from_file_location("protein_attach", os.path.join(HERE, "protein_attach.py"))
PA = importlib.util.module_from_spec(spec)
spec.loader.exec_module(PA)


def test_gated_paf_is_read_from_the_families_params():
    with tempfile.TemporaryDirectory() as d:
        fam = os.path.join(d, "x.fam")
        open(fam + ".params.tsv", "w").write("paf\tx\nexon_gate\ttrue\n")
        assert PA.gated_paf(fam + ".loci.paf")
        open(fam + ".params.tsv", "w").write("exon_gate\tfalse\n")
        assert not PA.gated_paf(fam + ".loci.paf")
        os.remove(fam + ".params.tsv")
        assert not PA.gated_paf(fam + ".loci.paf")


def test_a_gated_run_is_recognised_by_its_own_paf_name():
    with tempfile.TemporaryDirectory() as d:
        fam = os.path.join(d, "x.fam")
        open(fam + ".loci.gated.paf", "w").close()
        assert PA.gated_paf(fam + ".loci.gated.paf")
        assert PA.gated_paf(fam + ".loci.paf"), "no full PAF beside a gated one: the run was gated"
        open(fam + ".loci.paf", "w").close()
        assert not PA.gated_paf(fam + ".loci.paf"), "a full PAF (no exon_gate row in params)"
