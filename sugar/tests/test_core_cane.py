# (C) 2024, Tom Eulenfeld, MIT license
import pytest

from sugar import read, read_fts
from  sugar.core.cane import translate


def test_translate():
    s = 'CCC-AT-GAT-NCC--CCCCT---ANT-A--GGGN'
    aas = 'P-MX-PP-X-*G'

    assert translate(s, complete=True) == aas
    with pytest.warns(UserWarning, match='might be a stop'):
        assert translate(s[3:-3], warn=True) == aas[1:-2]
    assert translate(s[3:-3], warn=False) == aas[1:-2]
    with pytest.warns(UserWarning, match='not a start|might be a stop'):
        assert translate(s[8:-3], warn=True) == aas[3:-2]
    with pytest.raises(ValueError, match='not a start'):
        translate(s)


def test_translate_seq():
    seqs = read()['cds'].translate()
    assert str(seqs[0]) == seqs[0].fts.get('cds').meta._genbank.translation


def test_translate_final_stop():
    seq = read()['cds'][0]
    assert str(seq.copy().translate(complete=True))[-1] == '*'
    assert str(seq.copy().translate(complete=True, final_stop=True))[-1] == '*'
    assert str(seq.copy().translate(complete=True, final_stop=False))[-1] != '*'
    seq[3:6] = 'TAG'
    assert str(seq.copy().translate())[-1] != '*'
    assert len(seq.copy().translate()) == 1
    assert str(seq.copy().translate(final_stop=False))[-1] != '*'
    assert str(seq.copy().translate(final_stop=True))[-1] == '*'
    assert len(seq.copy().translate(final_stop=True)) == 2

# TODO more translation tests


def test_filter_fts():
    fts = read_fts()
    fts2 = fts.copy().filter(len_gt=30_000)
    assert len(fts2) == 5
    fts2 = fts.copy().filter(len_le=3000, type_eq='cDNA_match')
    assert len(fts2) == 1
    fts2 = fts.filter(len_min=1000, type_in=('CDS', 'exon'), inplace=False)
    assert len(fts2) == 2
    assert len(fts) > len(fts2)


def test_filter_seqs():
    seqs = read()
    seqs.filter(len_gt=9500)
    assert len(seqs) == 1


def test_groupby_fts_nested():
    fts = read_fts()
    d = fts.groupby(('type', len))
    assert 'CDS' in d
    keys = list(d['CDS'].keys())
    assert all(isinstance(k, int) for k in keys)
    assert isinstance(d['CDS'][keys[0]], fts.__class__)


def test_groupby_seqs():
    seqs = read()
    grouped = seqs.groupby()
    assert len(grouped) == 2
    assert set(grouped.keys()) == set(seqs.ids)


def test_sort_fts_nested():
    fts = read_fts()
    fts.sort(('type', len))
    types = [ft.type for ft in fts]
    assert types == sorted(types)
    assert len(fts[1]) <= len(fts[2])

def test_sort_seqs():
    seqs = read()
    seqs.sort(len)
    lens = [len(seq) for seq in seqs]
    assert lens == sorted(lens)
