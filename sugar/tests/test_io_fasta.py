# (C) 2024, Tom Eulenfeld, MIT license

from sugar import read, iter_
from sugar._io.fasta import _id_from_header

def test_fasta_in():
    seqs = read('!data/example.fasta')
    assert len(seqs) == 2
    for seq in iter_('!data/example.fasta'):
        assert str(seq[0]) in 'ML'


def test_fasta_id_from_header():
    assert _id_from_header('id42X') == 'id42X'
    assert _id_from_header('id42X dd    ') == 'id42X'
    assert _id_from_header('id42X,dd    ') == 'id42X'
    assert _id_from_header('id42X; dd   ') == 'id42X'
    assert _id_from_header('id42X| dd   ') == 'id42X'
    assert _id_from_header('|||gb|id42X|') == 'id42X'
    assert _id_from_header('|||gb:id42X|') == 'id42X'
