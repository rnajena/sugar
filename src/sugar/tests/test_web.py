# (C) 2023, Tom Eulenfeld, MIT license

import pytest
from os.path import join
import tempfile
from time import perf_counter

from sugar import read
from sugar.web import Entrez


@pytest.mark.webtest
def test_entrez():
    seqid = 'AB047639'
    with tempfile.TemporaryDirectory() as tempdir:
        client = Entrez(path=tempdir)
        fname = client.fetch_seq(seqid)
        assert fname.endswith(f'{seqid}.gb')
        seqs1 = read(fname)
        fname2 = client.fetch_seq(seqid)
        seq2 = client.get_seq(seqid)
        assert fname2 == fname
        assert seq2 == seqs1[0]
        fnames = client.fetch_basket([seqid])
        assert fnames[0] == fname
        seqs3 = client.get_basket([seqid])
        assert seqs3 == seqs1
        client2 = Entrez(path=join(tempdir, 'not_existing_dir'))
        fname3 = client2.fetch_seq(seqid)
        assert fname3.endswith(f'{seqid}.gb')
    assert fname == fname2
    assert len(client._request_times) == 1  # only single request


@pytest.mark.webtest
def test_entrez_time_constraint():
    seqid = 'AB047639'
    client = Entrez(path=None)
    client._requests = 1
    client._sleep = 1
    time0 = perf_counter()
    client.fetch_seq(seqid)
    client.fetch_seq(seqid)
    assert perf_counter() - time0 > 1
