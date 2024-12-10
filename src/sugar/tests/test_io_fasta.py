# (C) 2024, Tom Eulenfeld, MIT license
from io import StringIO

from sugar import read, iter_
from sugar._io.fasta import _id_from_header

def test_fasta_in():
    seqs = read('!data/example.fasta')
    assert len(seqs) == 2
    assert seqs[0].id == 'MCHU'
    assert seqs[1].id == 'AAD44166.1'
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
    assert _id_from_header('id42X gb:id5') == 'id42X'
    assert _id_from_header('  gb:5') is None


def test_fasta_header():
    seqs = read('!data/example.fasta')
    header = seqs[0].meta._fasta.header
    with StringIO() as f:
        seqs.write(f, 'fasta')
        f.seek(0)
        seqs2 = read(f)
    assert seqs2[0].meta._fasta.header == header
