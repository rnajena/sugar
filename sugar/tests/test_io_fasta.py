# (C) 2024, Tom Eulenfeld, MIT license

from sugar import read, iter_

def test_in():
    seqs = read('!data/example.fasta')
    assert len(seqs) == 2
    for seq in iter_('!data/example.fasta'):
        assert str(seq[0]) in 'ML'
