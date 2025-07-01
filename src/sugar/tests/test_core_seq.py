# (C) 2024, Tom Eulenfeld, MIT license

import pytest
from sugar import read, Attr, BioSeq, BioBasket, Feature, FeatureList


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


def test_countall_pandas():
    pytest.importorskip('pandas', reason='require pandas module')
    seqs = read()
    df = seqs.countall(rtype='df')
    assert abs(df['prob'].sum() - 2) < 1e-8


def test_countplot(tmpfname):
    pytest.importorskip('pandas', reason='require pandas module')
    pytest.importorskip('seaborn', reason='require seaborn module')
    seqs = read()
    seqs[0].countplot(plot=tmpfname)


def test_meta_str():
    meta = read()[0].meta
    assert 'id' in str(meta)
    assert 'CDS' in str(meta)


def test_seqs_str():
    seqs = read()
    seqs2 = seqs.copy()
    seqs2.data = []
    assert str(seqs2).startswith('0 seq')
    seqs2 = seqs.copy()
    seqs.data = 10 * seqs.data
    assert '...' in str(seqs2)


def test_seqs_shortcuts():
    seq = read()[0]
    assert seq.id == seq.meta.id
    assert seq.fts == seq.meta.fts
    seq.id = 'XXX'
    assert seq.id == seq.meta.id


def test_seqs_getitem():
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


def test_seqs_getitem_special():
    seq = read()[1]
    seq[4:8] = '----'
    seq2 = seq.sl(update_fts=True)[1:100]
    seq3 = seq.sl(update_fts=True, gap='-')[1:100-4]
    assert seq2.data == seq3.data
    assert len(seq2.fts) == 2
    assert all(len(ft) == 99 for ft in seq2.fts)
    assert all(len(ft) == 99-4 for ft in seq3.fts)
    assert all(ft.loc.defect == ft.loc.defect.MISS_LEFT | ft.loc.defect.MISS_RIGHT for ft in seq2.fts)
    seq2 = seq.sl(update_fts=True)['cds']
    assert seq2['cds'].data == seq['cds'].data
    assert len(seq2) == len(seq2['cds'])
    defect = seq2.fts.get('cds').loc.defect
    assert defect == defect.NONE

    seq.fts.get('cds').locs = [(0, 10), (15, 25), (30, 40)]
    assert len(seq['cds']) == 30
    with pytest.raises(ValueError, match='.*Sorry'):
        seq2 = seq.sl(update_fts=True)['cds']
    seq2 = seq['cds']
    assert len(seq2) < 40
    seq2 = seq.sl(filler='-')['cds']
    assert len(seq2) == 40
    seq3 = seq.sl(splitter='-----')['cds']
    assert seq3 == seq2
    seq3 = seq.sl(splitter='X', filler='-')['cds']
    assert len(seq3) > 40
    assert seq3.str.count('X') == 2
    assert seq3.str.replace('X', '') == seq2


def test_seqs_getitem_special_edgecases():
    seq = read()[1]
    seq.fts[1].loc.strand = '-'
    seq.fts[2].loc.strand = '-'
    seq2 = seq.sl(update_fts=True)[seq.fts[1]]
    seq3 = seq.sl(update_fts=True)[seq.fts[2]]
    assert seq2 == seq2[seq2.fts[1]]
    assert seq3 == seq3[seq3.fts[1]]


def test_seqs_getitem_typeerror():
    seqs = read()
    fts = seqs.fts
    with pytest.raises(TypeError):
        seqs[fts]

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


def test_bioseq_magic_methods():
    s1, s2 = read()
    with pytest.warns(UserWarning, match='Join'):
        assert len(s1+s2) == len(s1) + len(s2)
        assert str(s1 + s2) == str(s1 + s2.data) == str(s1.data + s2)
    s3 = s1.copy()
    s3 += s2.data
    assert len(s3) == len(s1) + len(s2)
    s1, s2 = read().sort()
    assert s1 < s2
    assert s1.id < s2.id
    with pytest.raises(TypeError):
        s1 < ''


def test_biobasket_magic_methods():
    seqs = read()
    seqs2 = seqs.copy()
    seqs2 += seqs
    assert len(seqs2) == len(seqs + seqs) == len(seqs.data + seqs) == 2 * len(seqs)
    seqs2 = seqs.copy()
    seqs2 &= seqs[:1]
    assert len(seqs2) == len(seqs & seqs[:1]) == len(seqs.data & seqs[:1]) == 1
    seqs2 = seqs.copy()
    seqs2 |= seqs[:1]
    assert len(seqs2) == len(seqs | seqs[:1]) == len(seqs.data | seqs[:1]) == 2
    seqs2 = seqs.copy()
    seqs2 -= seqs[:1]
    assert len(seqs2) == len(seqs - seqs[:1]) == len(seqs.data - seqs[:1]) == 1
    seqs2 = seqs[1:].copy()
    seqs2 ^= seqs[:1]
    assert len(seqs2) == len(seqs[1:] ^ seqs[:1]) == len(seqs.data[1:] ^ seqs[:1]) == 2
    assert seqs.data == seqs


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


def test_seqs_tostr():
    assert read()[0].tostr().strip() == read().tostr(add_header=False).splitlines()[0].strip()


def test_seqs_reverse_complement():
    seqs = read()
    assert seqs.copy().reverse().complement() == seqs.rc()


def test_seqs_init():
    seqs1 = BioBasket(['ATG'])
    with pytest.warns(match='initialized with a seq'):
        seqs2 = BioBasket(seqs1[0])
    with pytest.warns(match='initialized with a seq'):
        seqs3 = BioBasket('ATG')
    assert seqs2 == seqs1
    assert seqs3 == seqs1
