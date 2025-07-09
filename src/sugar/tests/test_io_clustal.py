# (C) 2024, Tom Eulenfeld, MIT license
import pytest
from sugar import read


def test_clustal_io(tmpfname):
    seqs = read('!data/example.clustal')
    assert len(seqs) == 3
    seqs.write(tmpfname, 'clustal')
    seqs2 = read(tmpfname)
    assert seqs2 == seqs
    seqs[1].data = seqs[1].data[:len(seqs[1])-1]
    with pytest.warns(match='different'):
        seqs.write(tmpfname, 'clustal')
    seqs3 = read(tmpfname)
    assert seqs3 == seqs
