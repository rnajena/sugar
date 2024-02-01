# (C) 2024, Tom Eulenfeld, MIT license

import glob
from importlib.resources import files
import tempfile
from io import StringIO

import pytest
import sugar
from sugar.io.main import detect, detect_ext, read, iter_


GLOBEXPR = str(files('sugar.tests.data').joinpath('*.*'))
FNAMES = glob.glob(GLOBEXPR)
# formats with read and write support in sugar.io module
TESTIOFMTS = ('fasta', 'sjson', 'stockholm')


def test_detect():
    for fname in FNAMES:
        assert (detect(fname) == detect_ext(fname) != None) or detect_ext(fname) is None


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
        with tempfile.NamedTemporaryFile() as f:
            seqs.write(f.name, fmt)
            # test read
            f.seek(0)
            seqs2 = read(f.name)
            assert isinstance(seqs2, sugar.BioBasket)
            for seq2, seq1 in zip(seqs2, seqs):
                assert str(seq2) == str(seq1)
            # test iter_
            f.seek(0)
            for seq2, seq1 in zip(iter_(f.name), seqs):
                assert str(seq2) == str(seq1)


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


def test_io_tool():
    try:
        import Bio
    except ImportError:
        pytest.skip('test needs biopython')
    seqs = read()
    with tempfile.NamedTemporaryFile() as f:
        seqs.write(f.name, 'fasta', tool='biopython')
        f.seek(0)
        seqs2 = read(f.name, tool='biopython')
        assert isinstance(seqs2, sugar.BioBasket)
        for seq2, seq1 in zip(seqs2, seqs):
            assert str(seq2) == str(seq1)
