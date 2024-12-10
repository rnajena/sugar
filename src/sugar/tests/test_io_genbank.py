# (C) 2024, Tom Eulenfeld, MIT license

from sugar import read, iter_, read_fts

def test_gb():
    fname = '!data/example.gb'
    seqs = read(fname)
    for i, seq in enumerate(iter_(fname)):
        assert seq == seqs[i]
    seqs2 = read(fname, exclude=('translation',))
    assert sum('translation' in f.meta._genbank for f in seqs[0].fts) == 1
    assert sum('translation' in f.meta._genbank for f in seqs2[0].fts) == 0
    assert seqs[0].fts.get('source').meta._genbank.misc == ['test', 'test2']


def test_gb_fts():
    fname = '!data/example.gb'
    fts = read_fts(fname)
    fts2 = read_fts(fname, exclude=('translation',))
    assert sum('translation' in f.meta._genbank for f in fts) == 2
    assert sum('translation' in f.meta._genbank for f in fts2) == 0
    assert fts.get('source').meta._genbank.misc == ['test', 'test2']
