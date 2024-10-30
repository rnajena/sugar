# (C) 2024, Tom Eulenfeld, MIT license

import pytest
from sugar import read, Attr, BioSeq, BioBasket, Feature, FeatureList
from sugar.tests.util import tempfilename


def test_siformat():
    from sugar.core.seq import _si_format
    assert _si_format(10000) == '10k'
    assert _si_format(0) == '0'


def test_attr():
    assert Attr(a=1) == Attr(a=1)
    assert Attr(a=1) != Attr(a=2)


def test_bioseq_equal():
    s1 = BioSeq('bla', id='5')
    s2 = BioSeq('bla', id='5')
    s3 = BioSeq('bla', id='6')
    assert s1 == s2
    assert s1 != s3
    s4 = BioSeq(s1)
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
    seqs = read()
    seqs2 = seqs.copy()
    assert seqs2 == seqs
    seqs2[0].data = 'NNN'
    assert seqs2 != seqs


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


def test_biobasket_str():
    seqs = read()
    seqs2 = seqs.copy()
    seqs2.data = []
    assert str(seqs2).startswith('0 seq')
    seqs2 = seqs.copy()
    seqs.data = 10 * seqs.data
    assert '...' in str(seqs2)


def test_shortcuts():
    seq = read()[0]
    assert seq.id == seq.meta.id
    assert seq.fts == seq.meta.fts
    seq.id = 'XXX'
    assert seq.id == seq.meta.id


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
        seq.meta.fts = FeatureList(
            [Feature(type='cds', start=1, stop=5)])
    seq2 = seqs[0]['cds']
    seqs2 = seqs[:, 'cds']
    seqs3 = seqs['cds']
    assert seq2 == seqs2[0]
    assert seqs2 == seqs3
    assert len(seqs2[0]) == 4
    # features do not change
    assert seqs2[0].meta.fts[0].loc.start == 1
    assert seqs[0][1:10].meta.fts[0].loc.start == 1
    assert seqs[:, ::-1][:, 'cds'] != seqs[:, 'cds']
    # assert 'orig_len' not in seqs[0][1:10].meta.features[0]
    # assert seqs[0][3:6].meta.features[0].stop == 3
    # assert seqs[0][3:6].meta.features[0].orig_len == 4
    # assert len(seqs[0][10:20].meta.features) == 0

    ## TODO!!!


def test_sl_slicable_inplace():
    seqs = read()
    assert seqs.sl()[:1] == seqs[:1]


def test_setitem():
    seqs = read()
    seqs[:, 2] = 'X'
    assert seqs[0, 2] == 'X'
    seqs[:, 2:4] = 'AB'
    assert seqs[0, 2:4] == 'AB'
    seqs[0] = 'ABC'
    assert isinstance(seqs[0], BioSeq)
    assert seqs[0] == 'ABC'
    seqs = read()
    seqs[:2] = ['AGT', 'TGA']
    assert str(seqs[0]) == 'AGT'


def test_add_fts():
    seqs = read()
    nfts = len(seqs.fts)
    seqs.add_fts([seqs.fts[0]])
    assert len(seqs.fts) == nfts + 1
    seqs = read()
    seq = seqs[0]
    nfts = len(seq.fts)
    ft = seq.fts[0]
    seq.add_fts([ft])
    assert len(seq.fts) == nfts + 1
    assert seq.fts[0] == ft
    assert seq.fts[1] == ft
    assert seq.fts[-1] != ft

    ft = seqs.fts[0]
    ft.seqid = 'unknown'
    with pytest.warns(UserWarning, match='.*unknown'):
        seqs.add_fts([ft])
    with pytest.warns(UserWarning, match='.*unknown'):
        seqs.fts = [ft]
    with pytest.warns(UserWarning, match='.*mismatch'):
        seqs[0].add_fts([ft])


def test_biobasket_rc():
    seqs = read()
    seqs2 = seqs.copy().rc()
    assert seqs[0].rc() == seqs2[0]


def test_repr():
    from sugar import Location, Meta
    seqs = read()
    assert eval(repr(seqs[0])) == seqs[0]
    assert eval(repr(seqs)) == seqs


def test_magic_methods():
    # TODO
    pass


def test_str_methods():
    seqs = read()
    seq = seqs[0]
    assert seq.copy().str.lower().str.islower()
    assert seq.copy().str.upper().str.isupper()
    assert seq.copy().str.swapcase().str.swapcase() == seq
    seq2 = seq.copy().str.center(len(seq)+4)
    assert seq2.str.startswith('  ') and seq2.str.endswith('  ')
    seq2.str.strip(' ')
    assert seq2 == seq
    seq2 = seq.copy().str.ljust(len(seq)+2)
    assert seq2.str.endswith('  ')
    seq2.str.rstrip(' ')
    assert seq2 == seq
    seq2 = seq.copy().str.rjust(len(seq)+2)
    assert seq2.str.startswith('  ')
    seq2.str.lstrip(' ')
    assert seq2 == seq
    assert seq.copy().str.removeprefix('ACCT') == seq[4:]
    assert seq.copy().str.removesuffix('TGT') == seq[:-3]
    assert seq.str.encode().decode() == seq.data
    assert seq.str.index('ACCT') == seq.str.find('ACCT') == 0
    assert seq.str.rindex('TGT') == seq.str.rfind('TGT') == len(seq) - 3
    assert seq.str.isalpha() and seq.str.isascii()
    seq.str.maketrans('A', 'T')
    assert len(seq.str.split('ACCT', 1)) == len(seq.str.rsplit('ACCT', 1)) == 2
    assert len(seq.str.splitlines()) == 1
    assert len(seqs.str.find('A')) == 2
    assert all(seqs.copy().str.lower().str.islower())
