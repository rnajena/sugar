# (C) 2024, Tom Eulenfeld, MIT license

from sugar import read, read_fts


def test_features():
    seq = read()[0]
    ft = seq.fts.get('cds')
    assert ft.type == 'CDS'
    #assert len(seq[ft]) == ft.stop - ft.start
    fts = seq.fts.select('cds')
    assert fts[0] == ft


def test_features_special():
    fts = read_fts()
    assert 'CDS' in fts.tostr()
    assert fts.select('CDS')[0].type == 'CDS'
    assert fts.get('CDS').type == 'CDS'
    assert len(fts.slice(0, 1000)) < len(fts)
    assert len(fts.d[fts[0].meta.seqid]) > 0
    fts.sort()
    assert fts.loc_range[0] == fts[0].locs[0].start
