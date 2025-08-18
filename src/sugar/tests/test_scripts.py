# (C) 2023, Tom Eulenfeld, MIT license
from subprocess import check_output
import sys
import unittest.mock

import pytest

from sugar import read, read_fts
from sugar.scripts import run
from sugar.tests.util import tempfilename


def test_print():
    with tempfilename() as fname:
        read().write(fname, 'fasta')
        assert b'ACCTG' in check_output(f'sugar print {fname}'.split())
        assert check_output(f'sugar print {fname} --raw'.split()).startswith(b'ACCTG')


def test_printf():
    with tempfilename() as fname:
        fts = read_fts()
        fts.write(fname, 'gff')
        assert b'region' in check_output(f'sugar printf {fname}'.split())
        assert check_output(f'sugar printf {fname} --raw'.split()).startswith(b'region')


def test_cat():
    with tempfilename() as fname:
        read().write(fname, 'fasta')
        assert check_output(f'sugar cat {fname}'.split()).startswith(b'>')
        with tempfilename() as fname2:
            assert b'' == check_output(f'sugar cat {fname} -o {fname2} -fo sjson'.split())


def test_catf():
    with tempfilename() as fname:
        read_fts().write(fname, 'gff')
        assert b'region' in check_output(f'sugar catf {fname}'.split())
        with tempfilename() as fname2:
            assert b'' == check_output(f'sugar catf {fname} -o {fname2} -fo gff'.split())


def test_merge():
    with tempfilename() as fname:
        read().write(fname, 'fasta')
        assert check_output(f'sugar merge {fname}'.split()).startswith(b'>')
        with tempfilename() as fname2:
            assert b'' == check_output(f'sugar merge {fname} -o {fname2} -fo sjson'.split())


def test_load(capsys):
    sys.modules['IPython'] = unittest.mock.MagicMock()
    with tempfilename() as fname:
        read().write(fname, 'fasta')
        run('load', fname=fname)
    captured = capsys.readouterr()
    assert 'Bye' in captured.out
    with tempfilename() as fname:
        read_fts().write(fname, 'gff')
        run('loadf', fname=fname)
    captured = capsys.readouterr()
    assert 'Bye' in captured.out


def test_rc():
    assert b'CAT' == check_output('sugar rc ATG'.split()).strip()
    with tempfilename() as fname:
        read().write(fname, 'fasta')
        assert check_output(f'sugar rc {fname}'.split()).startswith(b'>')
        with tempfilename() as fname2:
            assert b'' == check_output(f'sugar rc {fname} -o {fname2} -fo sjson'.split())


def test_translate():
    assert b'S*R' == check_output('sugar translate --complete TCTTGAAGG'.split()).strip()
    assert b'MS' == check_output('sugar translate ATGTCTTGAAGG'.split()).strip()
    with tempfilename() as fname:
        read().write(fname, 'fasta')
        assert check_output(f'sugar translate --complete {fname}'.split()).startswith(b'>')
        with tempfilename() as fname2:
            assert b'' == check_output(f'sugar translate --complete {fname} -o {fname2} -fo sjson'.split())


def test_tutorial(tmp_path_cd):
    assert b'' == check_output('sugar tutorial'.split())
    fnames = [f.name for f in tmp_path_cd.iterdir()]
    assert 'hits.blastn' in fnames
    assert len(fnames) == 5


def test_logo(tmp_path_cd):
    pytest.importorskip('matplotlib')
    assert b'' == check_output('sugar logo sugar_logo.png'.split())
    fnames = [f.name for f in tmp_path_cd.iterdir()]
    assert 'sugar_logo.png' in fnames


def test_index():
    # TODO: write test for index script
    pass
