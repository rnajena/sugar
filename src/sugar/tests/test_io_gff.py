# (C) 2024, Tom Eulenfeld, MIT license

from sugar import read, read_fts
from sugar.tests.util import _clean_seqs, _clean_fts, tempfilename


def test_gff_fts():
    fts = read_fts('!data/fts_example.gff')
    assert fts[0].meta._fmt == 'gff'


def test_gff_seq():
    seqs = read('!data/example.gff')
    assert len(seqs) == 1
    assert str(seqs[0]['cds']) == 'CAT'


def test_gtf():
    fts = _clean_fts(read_fts('!data/example.gff'))
    with tempfilename(suffix='.gtf') as fname:
        fts.write(fname)
        fts2 = _clean_fts(read_fts(fname))
    assert fts == fts2


    # TODO test write and load again
    # TODO test cds -

