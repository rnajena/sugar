# (C) 2024, Tom Eulenfeld, MIT license

import pytest
from sugar import read, read_fts


def test_fts_defect():
    from sugar.core.fts import Defect
    assert Defect(Defect.NONE)._reverse() == Defect(Defect.NONE)
    assert Defect(Defect.MISS_LEFT)._reverse() == Defect(Defect.MISS_RIGHT)
    defect = Defect(Defect.MISS_LEFT | Defect.MISS_RIGHT)
    assert defect._reverse() == defect
    assert Defect(Defect.BEYOND_LEFT)._reverse() == Defect(Defect.BEYOND_RIGHT)
    defect = Defect(Defect.BEYOND_LEFT | Defect.BEYOND_RIGHT)
    assert defect._reverse() == defect
    assert Defect(Defect.UNKNOWN_LEFT)._reverse() == Defect(Defect.UNKNOWN_RIGHT)
    defect = Defect(Defect.UNKNOWN_LEFT | Defect.UNKNOWN_RIGHT)
    assert defect._reverse() == defect


def test_fts_get_select():
    seq = read()[0]
    ft = seq.fts.get('cds')
    ft2 = seq.fts.get(['cds'])
    assert ft.type == 'CDS'
    assert ft == ft2
    assert len(seq[ft]) == ft.loc.stop - ft.loc.start
    fts = seq.fts.select('cds')
    fts2 = seq.fts.select(['cds'])
    assert fts[0] == ft
    assert fts2 == fts


def test_fts_str():
    fts = read_fts()
    assert 'CDS' in fts.tostr()


def test_fts_d():
    fts = read_fts()
    assert fts.d[fts[-1].id] == fts[-1]


def test_fts_rc():
    fts = read_fts()
    fts2 = fts.copy().rc(seqlen=fts.loc_range[1]).rc(seqlen=fts.loc_range[1])
    assert fts2 == fts


def test_fts_magic_methods():
    fts = read_fts()
    fts2 = fts.copy()
    fts2 += fts
    assert len(fts2) == len(fts + fts) == len(fts.data + fts) == 2 * len(fts)
    fts2 = fts.copy()
    fts2 &= fts[:1]
    assert len(fts2) == len(fts & fts[:1]) == len(fts.data & fts[:1]) == 1
    fts2 = fts.copy()
    fts2 |= fts[:1]
    assert len(fts2) == len(fts | fts[:1]) == len(fts.data | fts[:1]) == len(fts)
    fts2 = fts.copy()
    fts2 -= fts[:1]
    assert len(fts2) == len(fts - fts[:1]) == len(fts.data - fts[:1]) == len(fts) - 1
    fts2 = fts[1:].copy()
    fts2 ^= fts[:1]
    assert len(fts2) == len(fts[1:] ^ fts[:1]) == len(fts.data[1:] ^ fts[:1]) == len(fts)
    assert fts.data == fts


def test_fts_slice():
    fts = read_fts()
    assert len(fts.slice(0, 1000)) < len(fts)
    assert fts.slice(None, None) == fts
    fts2 = fts.slice(1, fts.loc_range[1]-1)
    assert len(fts2) == len(fts)
    assert fts[0].loc.defect == fts2[0].loc.defect.NONE
    assert fts2[0].loc.defect == fts2[0].loc.defect.MISS_LEFT | fts2[0].loc.defect.MISS_RIGHT


def test_ft_overlaps():
    fts = read_fts()
    assert fts.get('cds').overlaps(fts.get('region'))
    assert fts.get('cds').overlaps(fts.get('region').locs)


def test_ft_shortcuts():
    ft = read_fts()[0]
    assert ft.id == ft.meta.id
    assert ft.type == ft.meta.type
    assert ft.seqid == ft.meta.seqid
    ft.id = 'XXX'
    ft.type = 'CDS'
    ft.seqid = 'xxx'
    assert ft.id == ft.meta.id
    assert ft.type == ft.meta.type
    assert ft.seqid == ft.meta.seqid


def test_locs_magic_methods():
    fts = read_fts()
    assert fts.get('cds').locs < fts.get('cDNA_match').locs
    assert fts.get('cds').locs <= fts.get('cDNA_match').locs
    assert fts.get('cDNA_match').locs > fts.get('cds').locs
    assert fts.get('cDNA_match').locs >= fts.get('cds').locs
    assert fts.get('cds').locs - fts.get('cds').locs == 0


def test_locs_new():
    from sugar.core.fts import Location, LocationTuple
    with pytest.raises(ValueError, match='No location'):
        LocationTuple()
    with pytest.raises(ValueError, match='One of'):
        LocationTuple([Location(0, 1)], start=0, stop=1)
    with pytest.raises(TypeError, match='LocationTuple'):
        LocationTuple([(0,)])
    with pytest.raises(ValueError, match='Found multiple'):
        LocationTuple([Location(0, 1, '+'), Location(1, 2, '-')])
    with pytest.raises(ValueError, match='.*at least one location'):
        LocationTuple([])


def test_loc_hash():
    locs = read_fts()[0].locs
    assert locs in {locs}


def test_fts_tolists():
    fts = read_fts().select('CDS')
    for type_, start, stop, strand in fts.tolists():
        assert type_ == 'CDS'
        assert start == 61943
        assert strand == '+'


def test_fts_topandas():
    pandas = pytest.importorskip('pandas', reason='require pandas module')
    fts = read_fts().select('CDS')
    df = fts.topandas()
    assert isinstance(df, pandas.DataFrame)
    assert df['type'][0] == 'CDS'
    assert df['start'][0] == 61943
    assert df['strand'][0] == '+'
