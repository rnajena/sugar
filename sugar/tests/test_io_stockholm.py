# (C) 2024, Tom Eulenfeld, MIT license
import pytest

from sugar import read, iter_
from sugar.tests.util import normalize_content, tempfilename


def test_stk():
    fname = '!data/example.stk'
    seqs = read(fname)
    assert seqs.meta._stockholm.GF.ID == 'UPSK'
    assert seqs.meta._stockholm.GC.SS_cons == '.AAA....<<<<aaa....>>>>'
    assert 'virus RNA' in seqs.meta._stockholm.GF.RT
    for i, seq in enumerate(iter_(fname)):
        assert seq == seqs[i]
    with tempfilename() as fn2:
        seqs.write(fn2, 'stockholm')
        ignore=['#=GF RT']
        assert normalize_content(fname) != normalize_content(fn2)
        assert normalize_content(fname, ignore=ignore) == normalize_content(fn2, ignore=ignore)


def test_stk_more_metadata():
    fname = '!data/example_stockholm_more_metadata.stk'
    seqs = read(fname)
    assert seqs.todict()['O83071/192-246'].meta._stockholm.GR.SA.startswith('999')
    with tempfilename() as fn2:
        seqs.write(fn2, 'stockholm')
        ignore=['#=GF CC']
        assert normalize_content(fname) != normalize_content(fn2)
        assert normalize_content(fname, ignore=ignore) == normalize_content(fn2, ignore=ignore)
    fname = '!data/example_stockholm_more_metadata_linebreak.stk'
    seqs2 = read(fname)
    assert seqs2 == seqs


def test_row2fts2row():
    from sugar._io.stockholm import row2fts, fts2row
    row1 = '.....1....||....2....|3|...4...|'
    row2 = '.....1....||....2....|3|...4....'
    row3 = '|' + '.' * 1000 + '1|'
    fts1 = row2fts(row1)
    fts2 = row2fts(row2)
    fts3 = row2fts(row3)
    assert len(fts1) == 4
    assert len(fts2) == 4
    assert len(fts3) == 1
    assert fts1[0].name == '1'
    assert fts2[0].name == '1'
    Defect = fts1[0].loc.Defect
    assert fts1[0].loc.defect == Defect.MISS_LEFT
    assert fts1[1].loc.defect == Defect.NONE
    assert fts2[-1].loc.defect == Defect.MISS_RIGHT
    assert len(fts3[0]) == 1003
    row1b = fts2row(fts1)
    row2b = fts2row(fts2)
    row3b = fts2row(fts3)
    assert row1 == row1b
    assert row2 == row2b
    assert len(row3b) == 1003
    assert row3b.count('1') > 1
    assert row3b.startswith('|' + 10 * '.')
    assert row3b.endswith(10 * '.' + '|')
    fts3[0].loc.stop = 6
    fts3[0].name = 'long_name'
    with pytest.warns(match='too long'):
        row3b = fts2row(fts3)
    assert len(row3b) == 6
    assert row3b == '|long|'
