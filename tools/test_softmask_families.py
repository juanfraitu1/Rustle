import os, sys
sys.path.insert(0, os.path.dirname(__file__))
from softmask_families import listed

def test_majority_of_copies_lowercase(tmp_path):
    fa = tmp_path / 'g.fa'; seq = 'ACGTACGTAC' + 'acgtacgtac' + 'ACGTAacgta'  # 10 upper, 10 lower, 5 + 5; lines of 7 bases
    fa.write_text('>c1\n' + '\n'.join(seq[i:i + 7] for i in range(0, 30, 7)) + '\n')
    (tmp_path / 'g.fa.fai').write_text('c1\t30\t4\t7\t8\n')
    cp = tmp_path / 'c.tsv'
    cp.write_text('family_id\tchrom\texons\n'
                  'SAT\tc1\t11-20\nSAT\tc1\t12-15,17-20\nSAT\tc1\t1-10\n'  # 2 of 3 copies lowercase -> listed
                  'GENE\tc1\t1-10\nGENE\tc1\t11-20\n'                      # 1 of 2: not a majority
                  'HALF\tc1\t21-30\n')                                     # 0.5 lowercase < 0.9
    assert listed(str(cp), str(fa), 0.9) == (['SAT'], 3)
    assert listed(str(cp), str(fa), 0.5) == (['HALF', 'SAT'], 3)
