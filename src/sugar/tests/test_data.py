# (C) 2024, Tom Eulenfeld, MIT license

from glob import glob
from importlib.resources import files
import pytest
from sugar.data import submat, gcode


def test_submat():
    sm = submat('blosum62')
    assert sm['Q']['E'] == sm['E']['Q'] == 2
    with pytest.raises(FileNotFoundError):
        submat('unknown_matrix')


@pytest.mark.filterwarnings('ignore:Letter.*not found in table header')
def test_submat_all():
    globexpr = str(files('sugar.data.data_submat').joinpath('*'))
    for fname in glob(globexpr):
        if 'README' in fname:
            continue
        sm = submat(fname)
        assert len(sm.keys()) >= 4


def test_gcode():
    gc = gcode(1)
    assert gc.id == 1
    assert gc.name == 'Standard'
    assert gc.tt['ATG'] == 'M'
    assert gc.tt['GAY'] == 'D'
    assert gc.tt['CTN'] == 'L'
    assert gc.tt['CTY'] == 'L'
    assert gc.ttinv['D'] == {'GAT', 'GAC'}
    assert 'ATG' in gc.starts
    assert 'ATG' not in gc.astarts
    assert 'ATN' in gc.astarts
    assert 'NNN' in gc.astarts
    assert 'TAG' in gc.stops
    assert 'TAG' not in gc.astops
    assert 'TAN' in gc.astops
    assert 'NNN' in gc.astops
