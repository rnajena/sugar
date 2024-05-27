# (C) 2024, Tom Eulenfeld, MIT license

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

def test_stk2():
    fname = '!data/example2.stk'
    seqs = read(fname)
    assert seqs.todict()['O83071/192-246'].meta._stockholm.GR.SA.startswith('999')
    with tempfilename() as fn2:
        seqs.write(fn2, 'stockholm')
        ignore=['#=GF CC']
        assert normalize_content(fname) != normalize_content(fn2)
        assert normalize_content(fname, ignore=ignore) == normalize_content(fn2, ignore=ignore)
    fname = '!data/example2b.stk'
    seqs2 = read(fname)
    assert seqs2 == seqs
