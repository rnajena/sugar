# (C) 2024, Tom Eulenfeld, MIT license
from importlib.resources import files
import pytest

from sugar import read, iter_
from sugar.tests.util import normalize_content, tempfilename


def test_stockholm():
    fname = '!data/example.stk'
    seqs = read(fname)
    assert seqs.meta._stockholm.GF.ID == 'UPSK'
    assert seqs.meta._stockholm.GC.SS_cons == '.AAA....<<<<aaa....>>>>'
    assert 'virus RNA' in seqs.meta._stockholm.GF.RT
    for i, seq in enumerate(iter_(fname)):
        assert seq == seqs[i]
    with tempfilename() as fn2:
        seqs.write(fn2, 'stockholm')
        ignore=['#=GF RT', '#=GF CC']
        assert normalize_content(fname) != normalize_content(fn2)
        assert normalize_content(fname, ignore=ignore) == normalize_content(fn2, ignore=ignore)


def test_stockholm_more_metadata():
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


def test_stockholm_row2fts2row():
    from sugar._io.stockholm import row2fts, fts2row
    row1 = '-----1----||----2----|3|---4---|'
    row2 = '-----1----||----2----|3|---4----'
    row3 = '|' + '-' * 1000 + '1|'
    row4 = '...|---1---|....|2|...'

    fts1 = row2fts(row1)
    fts2 = row2fts(row2)
    fts3 = row2fts(row3)
    fts4 = row2fts(row4)
    assert len(fts1) == 4
    assert len(fts2) == 4
    assert len(fts3) == 1
    assert len(fts4) == 2
    assert fts1[0].name == '1'
    assert fts2[0].name == '1'
    defect = fts1[0].loc.defect
    assert fts1[0].loc.defect == defect.MISS_LEFT
    assert fts1[1].loc.defect == defect.NONE
    assert fts2[0].loc.defect == defect.MISS_LEFT
    assert fts2[-1].loc.defect == defect.MISS_RIGHT
    assert len(fts3[0]) == 1003

    row1b = fts2row(fts1)
    row2b = fts2row(fts2)
    row3b = fts2row(fts3)
    row4b = fts2row(fts4, lensec=len(row4))
    assert row1 == row1b
    assert row2 == row2b
    assert len(row3b) == 1003
    assert row3b.count('1') > 1
    assert row3b.startswith('|' + 10 * '-')
    assert row3b.endswith(10 * '-' + '|')
    assert row4b == row4

    fts1b = row2fts(row1b)
    fts2b = row2fts(row2b)
    fts3b = row2fts(row3b)
    fts4b = row2fts(row4b)
    assert fts1b == fts1
    assert fts2b == fts2
    assert fts3b == fts3
    assert fts4b == fts4

    # test long name
    fts3[0].loc.stop = 6
    fts3[0].name = 'long_name'
    with pytest.warns(match='too long'):
        row3b = fts2row(fts3)
    assert len(row3b) == 6
    assert row3b == '|long|'

    # test different chars
    row4c = row4b.replace('|', ':').replace('.', ' ').replace('-', ' ')
    fts4c = row2fts(row4c, splitchar=':', fillchars=' ')
    assert fts4c == fts4
    row4d = fts2row(fts4c, splitchar=':', inchar=' ', outchar=' ', lensec=len(row4c))
    assert row4d == row4c


def test_stockholm_multiline():
    fname = '!data/example_stockholm_multi.stk'
    seqs = read(fname)
    assert len(seqs) == 3
    N = len(seqs[0])
    assert len(seqs[1]) == N
    assert len(seqs[2]) == N
    assert len(seqs.meta._stockholm.GC.SS_cons) == N
    assert len(seqs[2].meta._stockholm.GR.AS) == N


def test_stockholm_multialignment():
    fname = str(files('sugar.tests.data').joinpath('example_stockholm_multi.stk'))
    with open(fname) as f:
        seqs1 = read(f, 'stockholm')  # read 1st alignment
        seqs2 = read(f, 'stockholm')  # read 2nd alignment
    assert len(seqs1) == 3
    assert len(seqs2) == 5
