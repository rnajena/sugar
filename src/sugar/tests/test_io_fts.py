# (C) 2024, Tom Eulenfeld, MIT license

import sugar
from sugar import read_fts
import tempfile
from sugar.tests.util import _clean_fts, tempfilename


def test_gff():
    fts = read_fts()
    with tempfilename(suffix='.gff') as fname:
        fts.write(fname)
        fts2 = read_fts(fname)
        assert isinstance(fts2, sugar.FeatureList)
        assert fts2 == fts


def test_ft_gff():
    ft = read_fts()[0]
    with tempfilename(suffix='.gff') as fname:
        ft.write(fname)
        ft2 = read_fts(fname)[0]
        assert isinstance(ft2, sugar.Feature)
        assert ft2 == ft


def test_genbank2gff():
    fts = read_fts('!data/example.gb')
    with tempfilename(suffix='.gff') as fname:
        fts.write(fname)
        fts2 = read_fts(fname)
        assert isinstance(fts2, sugar.FeatureList)
        assert _clean_fts(fts2) == _clean_fts(fts)
        assert 'seqid' in fts2[0].meta


def test_blast2gff():
    fts = read_fts('!data/example.gb')
    with tempfilename(suffix='.gff') as fname:
        fts.write(fname)
        fts2 = read_fts(fname)
        assert isinstance(fts2, sugar.FeatureList)
        assert _clean_fts(fts2) == _clean_fts(fts)
        assert 'seqid' in fts2[0].meta
