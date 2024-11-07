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
    asser
    t 'CDS' in fts.tostr()

def test_fts_slice():
    fts = read_fts()
    assert len(fts.slice(0, 1000)) < len(fts)


def test_fts_d():
    fts = read_fts()
    assert len(fts.d[fts[0].meta.seqid]) > 0
