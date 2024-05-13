# (C) 2024, Tom Eulenfeld, MIT license

from sugar import read, read_fts


def test_gff_fts():
    fts = read_fts('!data/fts_example.gff')
    assert fts[0].meta._fmt == 'gff'

def test_gff_seq():
    seqs = read('!data/example.gff')
    assert len(seqs) == 1
    assert str(seqs[0]['cds']) == 'CAT'
    
    # TODO test write and load again
    # TODO test cds -