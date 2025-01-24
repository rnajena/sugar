# (C) 2024, Tom Eulenfeld, MIT license
import pytest

from sugar import read, read_fts, BioSeq, BioBasket
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

def test_match():
    seq = BioSeq('NNNUAGDDDUAGAUG')
    seqs = BioBasket([seq])
    seq2 = BioSeq('-UU-U-AG')
    assert seq.match('stop').start() == 3
    assert seq.match('start').end() == len(seq)
    matches = seq.matchall('stop')
    assert matches[0].span() == seq.match('stop').span()
    assert len(matches) == 2
    assert seqs.match('stop')[0].start() == 3
    matches = seq2.matchall('stop', gap=None)
    assert seqs.matchall('stop')[0].start() == 3
    assert len(matches) == 0
    match = seq2.match('stop', gap='-')
    assert match.group() == 'U-AG'
    assert seq2.match('stop', gap='-', rf=1) == None
    assert seq2.match('stop', gap='-', rf=2).group() == 'U-AG'
    assert seq2.match('stop', gap='-', rf=(1, 2)).group() == 'U-AG'
    assert seq2.match('stop', gap='-', rf=(0, 1)) == None
    seq3 = seq2.copy().rc()
    match3 = seq3.match('stop', gap='-', rf='bwd')
    assert match.span() == match3._match.span()
    assert match.span() != match3.span()


def test_orf():
    seqs = read()
    orfs = seqs[0].find_orfs()
    assert len(orfs) > 0
    longest_orf = orfs.sort(len)[-1]
    assert seqs[0][longest_orf] == seqs[0]['cds']

    orfs2 = seqs[0].find_orfs(rf='both')
    assert len(orfs2) > len(orfs)

    orfs = seqs.find_orfs()
    for id_ in seqs.ids:
        assert seqs.d[id_][orfs.groupby('seqid')[id_].sort(len)[-1]] == seqs.d[id_]['cds']


def test_select_fts_second_modus():
    fts = read_fts()
    fts2 = fts.copy().select(len_gt=30_000)
    assert len(fts2) == 5
    fts2 = fts.copy().select(len_le=3000, type_eq='cDNA_match')
    assert len(fts2) == 1
    fts2 = fts.select(len_min=1000, type_in=('CDS', 'exon'), inplace=False)
    assert len(fts2) == 2
    assert len(fts) > len(fts2)


def test_select_seqs():
    seqs = read()
    seqs.select(len_gt=9500, inplace=True)
    assert len(seqs) == 1
    seqs = read()
    seqs2 = seqs.select(len_gt=9500, inplace=False)
    assert len(seqs2) == 1
    assert len(seqs2) < len(seqs)


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
    fts.sort()
    assert fts.loc_range[0] == fts[0].locs[0].start

def test_sort_seqs():
    seqs = read()
    seqs.sort(len)
    lens = [len(seq) for seq in seqs]
    assert lens == sorted(lens)
