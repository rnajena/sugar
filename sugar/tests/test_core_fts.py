# (C) 2024, Tom Eulenfeld, MIT license

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


def test_fts_slice():
    fts = read_fts()
    assert len(fts.slice(0, 1000)) < len(fts)


def test_fts_d():
    fts = read_fts()
    assert len(fts.d[fts[0].meta.seqid]) > 0


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
    assert fts.slice(None, None) == fts
    fts2 = fts.slice(1, fts.loc_range[1]-1)
    assert len(fts2) == len(fts)
    assert fts[0].loc.defect == fts2[0].loc.Defect.NONE
    assert fts2[0].loc.defect == fts2[0].loc.Defect.MISS_LEFT | fts2[0].loc.Defect.MISS_RIGHT
