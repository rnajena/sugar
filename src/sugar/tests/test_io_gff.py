# (C) 2024, Tom Eulenfeld, MIT license

import pytest
from sugar import read, read_fts
from sugar.tests.util import _clean_fts, tempfilename


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


def test_gff_quote():
    fts = read_fts('!data/fts_example.gff')
    fts[0].seqid = '%,;&'
    fts[0].name = ';-), :-)'
    with tempfilename(suffix='.gff') as fname:
        fts.write(fname)
        fts2 = read_fts(fname)
    assert fts2[0].seqid == fts[0].seqid
    assert fts2[0].name == fts[0].name
    assert ':' in fts.tofmtstr('gff')


def test_gff_parsing_error():
    with pytest.raises(ValueError, match='Parsing error'):
        fts = read_fts('!data/io_gff_parsing_error.gff')
    with pytest.raises(ValueError, match='Parsing error'):
        fts = read_fts('!data/io_gff_parsing_error.gff', 'gtf')
    with pytest.raises(ValueError, match='Parsing error'):
        fts = read_fts('!data/io_gff_parsing_error.gff', 'gff', gff_version='3')


    # TODO test write and load again
    # TODO test cds -

