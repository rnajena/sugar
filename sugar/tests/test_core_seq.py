# (C) 2024, Tom Eulenfeld, MIT license

import pytest
from sugar import read, Attr, BioSeq, BioBasket, Feature, FeatureList
from sugar.core.seq import MutableMetaString as MMS
from sugar.tests.util import tempfilename


def test_attr():
    assert Attr(a=1) == Attr(a=1)
    assert Attr(a=1) != Attr(a=2)

def test_mutablemetastring():
    s1 = MMS('bla', id='5')
    s2 = MMS('bla', id='5')
    s3 = MMS('bla', id='6')
    assert s1 == s2
    assert s1 != s3
    s4 = MMS(s1)
    assert s1 == s4


def test_bioseq_to_from_biopython():
    SeqRecord = pytest.importorskip('Bio.SeqRecord', reason='need biopython')
    seqs = read()
    seq = seqs[0]
    obj = seq.toobj('biopython')
    seq2 = seq.fromobj(obj)
    assert isinstance(obj, SeqRecord.SeqRecord)
    assert seq2.id == obj.id == seq.id
    assert str(seq2) == str(obj.seq) == str(seq)


def test_biobasket_to_from_biopython():
    SeqRecord = pytest.importorskip('Bio.SeqRecord', reason='need biopython')
    seqs = read()
    obj = seqs.toobj('biopython')
    seqs2 = seqs.fromobj(obj)
    assert isinstance(obj[0], SeqRecord.SeqRecord)
    assert seqs2[0].id == obj[0].id == seqs[0].id
    assert str(seqs2[0]) == str(obj[0].seq) == str(seqs[0])


def test_todict():
    seqs = read()
    d = seqs.todict()
    assert set(d.keys()) == {seq.id for seq in seqs}
    assert sorted(d.values()) == sorted(seqs)
    assert seqs.d == d


def test_complement():
    seq = BioSeq('ACGT')
    seq2 = seq.copy().complement()
    assert str(seq2) != str(seq)
    seq2.complement()
    assert str(seq2) == str(seq)

    seq = BioSeq('ACGU')
    seq2 = seq.copy().complement()
    assert str(seq2) != str(seq)
    seq2.complement()
    assert str(seq2) == str(seq)


def test_ids():
    seqs = read()
    assert seqs.ids[0] == seqs[0].id


def test_copy():
    seq = read()[0]
    n = len(seq)
    assert seq.copy()[10:] != seq
    assert len(seq.copy()[10:]) == n - 10
    assert seq.copy() == seq


def test_countall():
    seqs = read()
    assert 'T' in seqs[0].countall()
    assert abs(sum(seqs[0].countall(rtype='prob').values()) - 1) < 1e-8
    try:
        import pandas
    except ImportError:
        pass
    else:
        df = seqs.countall(rtype='df')
        assert abs(df['prob'].sum() - 2) < 1e-8


def test_countplot():
    pytest.importorskip('pandas', reason='need pandas')
    pytest.importorskip('seaborn', reason='need seaborn')
    seqs = read()
    with tempfilename() as fname:
        seqs[0].countplot(plot=fname)


def test_meta_str():
    meta = read()[0].meta
    assert 'id' in str(meta)
    assert 'CDS' in str(meta)


def test_shortcuts():
    seq = read()[0]
    assert seq.id == seq.meta.id
    assert seq.fts == seq.meta.features


def test_getitem():
    seqs = read()
    olen = len(seqs[0])
    assert isinstance(seqs[0], BioSeq)
    assert isinstance(seqs[0:1], BioBasket)
    assert len(seqs[0:1]) == 1
    assert seqs[:] == seqs[::-1][::-1]
    assert seqs is not seqs[:]
    assert seqs[:][0] is seqs[0]
    assert len(seqs[0, 1:5]) == 4
    assert seqs[0:1, 1:5][0] == seqs[0, 1:5]
    assert seqs[0:1, 1:5][0] is not seqs[0, 1:5]

    assert len(seqs[0]) == olen
    for seq in seqs:
        seq.meta.features = FeatureList(
            [Feature(type='cds', start=1, stop=5)])
    seq2 = seqs[0]['cds']
    seqs2 = seqs[:, 'cds']
    seqs3 = seqs['cds']
    assert seq2 == seqs2[0]
    assert seqs2 == seqs3
    assert len(seqs2[0]) == 4
    # features do not change
    assert seqs2[0].meta.features[0].loc.start == 1
    assert seqs[0][1:10].meta.features[0].loc.start == 1
    assert seqs[:, ::-1][:, 'cds'] != seqs[:, 'cds']
    # assert 'orig_len' not in seqs[0][1:10].meta.features[0]
    # assert seqs[0][3:6].meta.features[0].stop == 3
    # assert seqs[0][3:6].meta.features[0].orig_len == 4
    # assert len(seqs[0][10:20].meta.features) == 0


def test_setitem():
    seqs = read()
    seqs[:, 2] = 'X'
    assert seqs[0, 2] == 'X'
    seqs[:, 2:4] = 'AB'
    assert seqs[0, 2:4] == 'AB'
    seqs[0] = 'ABC'
    assert isinstance(seqs[0], BioSeq)
    assert seqs[0] == 'ABC'


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
    matches = seq2.matchall('stop')
    assert len(matches) == 0
    match = seq2.match('stop', gap='-')
    assert match.group() == 'U-AG'
