# (C) 2024, Tom Eulenfeld, MIT license

import sugar
from sugar import read_fts
import tempfile
from sugar.tests.util import _clean_fts


def test_gff():
    fts = read_fts()
    with tempfile.NamedTemporaryFile(suffix='.gff') as f:
        fts.write(f.name)
        fts2 = read_fts(f.name)
        assert isinstance(fts2, sugar.FeatureList)
        assert fts2 == fts


def test_ft_gff():
    ft = read_fts()[0]
    with tempfile.NamedTemporaryFile(suffix='.gff') as f:
        ft.write(f.name)
        ft2 = read_fts(f.name)[0]
        assert isinstance(ft2, sugar.Feature)
        assert ft2 == ft


def test_genbank2gff():
    fts = read_fts('!data/example.gb')
    with tempfile.NamedTemporaryFile(suffix='.gff') as f:
        fts.write(f.name)
        fts2 = read_fts(f.name)
        assert isinstance(fts2, sugar.FeatureList)
        assert _clean_fts(fts2) == _clean_fts(fts)
        assert 'seqid' in fts2[0].meta


def test_blast2gff():
    fts = read_fts('!data/example.gb')
    with tempfile.NamedTemporaryFile(suffix='.gff') as f:
        fts.write(f.name)
        fts2 = read_fts(f.name)
        assert isinstance(fts2, sugar.FeatureList)
        assert _clean_fts(fts2) == _clean_fts(fts)
        assert 'seqid' in fts2[0].meta
