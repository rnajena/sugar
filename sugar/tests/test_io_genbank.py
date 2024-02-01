# (C) 2024, Tom Eulenfeld, MIT license

from sugar.io.main import read, iter_

def test_gb():
    fname = '!data/example.gb'
    seqs = read(fname)
    for i, seq in enumerate(iter_(fname)):
        assert seq == seqs[i]
    seqs2 = read(fname, translation=False)
    assert sum('translation' in f for f in seqs[0].meta.features) == 1
    assert sum('translation' in f for f in seqs2[0].meta.features) == 0
    assert seqs[0].fts.get('source').misc == ['test', 'test2']
