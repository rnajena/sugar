# (C) 2024, Tom Eulenfeld, MIT license
import io
import pytest
from sugar import read_fts


def test_tsv():
    pytest.importorskip('pandas', reason='require pandas module')
    fts = read_fts()
    out = fts.write(None, 'tsv')
    fts2 = read_fts(io.StringIO(out))
    fts3 = read_fts('!data/fts_example.tsv')
    assert fts2 == fts3


def test_csv():
    pytest.importorskip('pandas', reason='require pandas module')
    fts = read_fts()
    out = fts.write(None, 'csv')
    fts2 = read_fts(io.StringIO(out))
    fts3 = read_fts('!data/fts_example.csv')
    assert fts2 == fts3


def test_xsv_edgecases():
    pytest.importorskip('pandas', reason='require pandas module')
    fts = read_fts()
    out = fts.write(None, 'csv', sep=' ', keys='type start len strand defect')
    fts2 = read_fts(io.StringIO(out), sep=None, engine='python')
    fts3 = read_fts('!data/fts_example.tsv')
    assert fts2 == fts3
    out = fts.write(None, 'csv', sep='|', keys='type stop len strand defect')
    fts2 = read_fts(io.StringIO(out), sep='|')
    assert fts2 == fts3
