# (C) 2024, Tom Eulenfeld, MIT license

import glob
from importlib.resources import files
from io import StringIO
import os.path
import tempfile

import pytest
import sugar
from sugar import read, iter_, read_fts
from sugar._io.main import detect, detect_ext
from sugar.tests.util import tempfilename


GLOBEXPR = str(files('sugar.tests.data').joinpath('example*.*'))
GLOBEXPR_FTS = str(files('sugar.tests.data').joinpath('fts_example*.*'))
FNAMES = glob.glob(GLOBEXPR)
FNAMES_FTS = glob.glob(GLOBEXPR_FTS)
# formats with read and write support in sugar._io module
TESTIOFMTS = ('fasta', 'sjson', 'stockholm')


def test_detect():
    for fname in FNAMES:
        assert (detect(fname) == detect_ext(fname) != None) or detect_ext(fname) is None


def test_detect_fts():
    for fname in FNAMES_FTS:
        assert (detect(fname, 'fts') == detect_ext(fname, 'fts') != None) or detect_ext(fname, 'fts') is None


def test_read_fts():
    for fname in FNAMES_FTS:
        try:
            fts = read_fts(fname)
        except ImportError:  # ignore pandas ImportError for csv and tsv
            continue
        assert isinstance(fts, sugar.FeatureList)
        assert len(fts) > 0



def test_read():
    for fname in FNAMES:
        seqs = read(fname)
        assert isinstance(seqs, sugar.BioBasket)
        assert len(seqs) > 0


def test_iter():
    for fname in FNAMES:
        for seq in iter_(fname):
            assert isinstance(seq, sugar.BioSeq)
            assert len(seq) > 0


def test_io_file():
    seqs = read()
    for fmt in TESTIOFMTS:
        with tempfilename() as fname:
            seqs.write(fname, fmt)
            # test read
            seqs2 = read(fname)
            assert isinstance(seqs2, sugar.BioBasket)
            for seq2, seq1 in zip(seqs2, seqs):
                assert str(seq2) == str(seq1)
            # test iter_
            for seq2, seq1 in zip(iter_(fname), seqs):
                assert str(seq2) == str(seq1)


def test_write_fmtstr_seq():
    seqs = read()
    with tempfilename() as fname:
        seqs[0].write(fname, 'fasta')
        seqs2 = read(fname)
        assert str(seqs2[0]) == str(seqs[0])
    s = seqs[0].tofmtstr('fasta')
    seqs2 = seqs.fromfmtstr(s)
    assert str(seqs2[0]) == str(seqs[0])


def test_io_fmtstr():
    seqs = read()
    for fmt in TESTIOFMTS:
        s = seqs.tofmtstr(fmt)
        # test read via fromfmtstr
        seqs2 = seqs.fromfmtstr(s)

        # seqs2 = read(fh)
        assert isinstance(seqs2, sugar.BioBasket)
        for seq2, seq1 in zip(seqs2, seqs):
            assert str(seq2) == str(seq1)
        # test iter_
        fh = StringIO(s)
        for seq2, seq1 in zip(iter_(fh), seqs):
            assert str(seq2) == str(seq1)


def test_read_glob():
    seqs = read(GLOBEXPR)
    assert isinstance(seqs, sugar.BioBasket)
    assert len(seqs) >= 2
    for seq in iter_(GLOBEXPR):
        assert isinstance(seq, sugar.BioSeq)


def test_uncompress():
    assert read('!data/io_test.gz', archive='gz') == read()
    assert read('!data/io_test.gz') == read()


def test_archive():
    assert read()[0] in read('!data/io_test.zip')


def test_write_archive():
    with tempfile.TemporaryDirectory() as tmpdir:
        fname = os.path.join(tmpdir, 'test.fasta')
        seqs = read()
        seqs.write(fname, archive='zip')
        seqs2 = read(fname + '.zip')
        assert len(seqs) == len(seqs2)
        assert str(seqs[0]) == str(seqs2[0])
        seqs.write(fname + '.gff', 'GFF', archive=True)
        seqs.write(fname + '.zip')


@pytest.mark.webtest
def test_download():
    url = 'https://raw.githubusercontent.com/rnajena/sugar/master/src/sugar/tests/data/example.gb'
    assert read(url) == read()


@pytest.mark.webtest
def test_download_uncompress():
    url = 'https://raw.githubusercontent.com/rnajena/sugar/master/src/sugar/tests/data/io_test.gz'
    assert read(url) == read()


@pytest.mark.webtest
def test_download_zip():
    url = 'https://raw.githubusercontent.com/rnajena/sugar/master/src/sugar/tests/data/io_test.zip'
    assert read()[0] in read(url)


def test_empty_file():
    with tempfilename() as fname:
        with open(fname, 'wb') as f:
            f.write(b'   \n')
        with pytest.raises(match='empty file'):
            read(fname)
        with pytest.raises(match='empty file'):
            read_fts(fname)
        with pytest.raises(match='empty file'):
            for _ in iter_(fname):
                pass
