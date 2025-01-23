# (C) 2024, Tom Eulenfeld, MIT license

from contextlib import redirect_stderr, redirect_stdout
from importlib.resources import files
import io
import os
import pytest
import tempfile
from sugar import FastaIndex
from sugar.scripts import cli
import sys

pytest.importorskip('binarysearchfile', reason='require binarysearchfile module')

try:
    import platformdirs
except ImportError:
    platformdirs = None


@pytest.mark.xfail(sys.platform == 'darwin', reason='db module not working as expected for MacOS')
def test_fastaindex():
    fastafiles = str(files('sugar.tests.data').joinpath('example*.fasta'))
    with tempfile.TemporaryDirectory() as tmpdir:
        fname = os.path.join(tmpdir, 'test_bsf.sugarindex')
        fname2 = os.path.join(tmpdir, 'test_db.sugarindex')
        f1 = FastaIndex(fname, create=True, mode='binary')
        with FastaIndex(fname2, create=True, mode='db') as f2:
            f1.add(fastafiles, silent=True)
            f2.add(fastafiles, silent=True)
            assert f1.totalsize < f2.totalsize
            id_ = 'BTBSCRYR'
            assert f1.get_basket(id_)[0].str.startswith('tgcaccaaacatgtcta'.upper())
            assert f1.get_basket(id_)[0] == f1.get_seq(id_)
            assert id_ in f1.get_fasta(id_)
            assert id_ in f1.get_fastaheader(id_)
            assert f1.get_basket(id_) == f2.get_basket(id_)
            assert f1.get_fasta(id_) == f2.get_fasta(id_)
            assert f1.get_fastaheader(id_) == f2.get_fastaheader(id_)
            assert('dbname' in str(f1))
            assert('dbname' in str(f2))
            seqids = (id_, 'MCHU')
            assert len(list(f1.iter(seqids))) == 2
            assert len(list(f1.iter_fasta(seqids))) == 2
            assert len(list(f1.iter_fastaheader(seqids))) == 2
            assert list(f1.iter(seqids)) == list(f2.iter(seqids))
            assert list(f1.iter_fasta(seqids)) == list(f2.iter_fasta(seqids))
            assert list(f1.iter_fastaheader(seqids)) == list(f2.iter_fastaheader(seqids))
            assert len(f1) == len(f2)

        # reload files
        # extracting only part of sequence
        f1 = FastaIndex(fname)
        with FastaIndex(fname2) as f2:
            id1, id2 = 'BTBSCRYR', 'MCHU'
            assert f1.get_basket((id1, None, 10)) == f1.get_basket(id1)[:, :10]
            assert f1.get_basket([(id1, None, 10), (id2, None, 10)]) == f1.get_basket([id1, id2])[:, :10]
            assert f1.get_basket([(id1, 5, 9), (id2, 5, 9)]) == f1.get_basket([id1, id2])[:, 5:9]
            assert f1.get_basket([(id1, None, None), (id2, None, None)]) == f1.get_basket([id1, id2])
            with pytest.warns(UserWarning, match='Start index'):
                assert f1.get_basket([(id1, -10, None), (id2, -10, None)]) == f1.get_basket([id1, id2])
            with pytest.warns(UserWarning, match='End index'):
                assert f1.get_basket([(id1, None, 1000), (id2,  None, 1000)]) == f1.get_basket([id1, id2])

        # recreate files
        # check readding stuff
        os.remove(fname)
        try:
            os.remove(fname2)
        except FileNotFoundError:
            # work-around for dbm.dumb, which creates 3 files
            from glob import glob
            for _fname in glob(fname2 + '.*'):
                os.remove(_fname)
        with pytest.warns(UserWarning, match='new file'):
            f1 = FastaIndex(fname, mode='binary')
        f1.add(fastafiles, seek=10, silent=True)
        assert len(f1) == 1
        with pytest.raises(ValueError, match='performance'):
            f1.add(fastafiles)
        f1.add(fastafiles, seek=10, force=True, silent=True)
        assert len(f1) == 2
        with FastaIndex(fname2, create=True, mode='db') as f2:
            f2.add(fastafiles, seek=10, silent=True)
            assert len(f2) == 1
            f2.add(fastafiles, seek=10, silent=True)
            assert len(f2) == 1


@pytest.mark.xfail(sys.platform == 'darwin', reason='db module not working as expected for MacOS')
def test_fastaindex_script():
    fastafiles = str(files('sugar.tests.data').joinpath('example*.fasta'))
    with tempfile.TemporaryDirectory() as tmpdir:
        fname = os.path.join(tmpdir, 'test.sugarindex')
        cli(f'index create -m binary {fname}'.split())
        odb = f'-d {fname}' if platformdirs is None else ''
        with redirect_stdout(None):
            cli(f'index info {odb}'.split())
        with redirect_stderr(io.StringIO()):
            cli(f'index add {odb} {fastafiles}'.split())
        with redirect_stdout(None):
            cli(f'index print {odb} BTBSCRYR'.split())
            cli(f'index print {odb} BTBSCRYR,1,5'.split())
        data = """BTBSCRYR, 5, 10
        MCHU, 5, 10"""
        fnameinp = os.path.join(tmpdir, 'inp.txt')
        with open(fnameinp, 'w') as f:
            f.write(data)
        with redirect_stdout(None):
            cli(f'index fetch {odb} {fnameinp}'.split())
        # os.system(f'cat {fnameinp} | sugar index fetch -')
