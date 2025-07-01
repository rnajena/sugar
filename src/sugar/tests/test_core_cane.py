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


def test_find_orfs():
    seqs = read()
    orfs = seqs[0].find_orfs(rf='fwd')
    assert len(orfs) > 0
    longest_orf = orfs.sort(len)[-1]
    assert seqs[0][longest_orf] == seqs[0]['cds']

    orfs2 = seqs[0].find_orfs(rf='all')
    assert len(orfs2) > len(orfs)

    orfsb = seqs[0].find_orfs(rf='fwd', nested='all')
    assert len(orfsb) > len(orfs)

    orfsc = seqs[0].find_orfs(rf='fwd', need_stop=False)
    assert len(orfsc) == len(orfs) + 2

    orfs4 = seqs[0].find_orfs(rf='all', nested='all')
    assert len(orfs4) > len(orfs2)

    orfs = seqs.find_orfs()
    for id_ in seqs.ids:
        assert seqs.d[id_][orfs.groupby('seqid')[id_].sort(len)[-1]] == seqs.d[id_]['cds']


def test_find_orfs_extended():
    from functools import reduce
    from operator import add
    teststr = 'CCC_CAT_GCT_GAC_CTA_CCC_CAT_CAT_GTG_ACC_CCC_CTG_ACC_CAT_GCC_CTG_ACC_ATG_CCC_C'
    seq = BioSeq(teststr.replace('_', ''))
    assert (seq.find_orfs(rf='fwd') + seq.find_orfs(rf='bwd')) == seq.find_orfs(rf='all')
    assert reduce(add, [seq.find_orfs(rf=rf) for rf in (0, 1, 2, -1, -2, -3)]) == seq.find_orfs(rf='all')
    orfs1 = seq.find_orfs()
    orfs2 = seq.find_orfs(need_start='never')
    orfs3 = seq.find_orfs(need_start='once')
    orfs1b = seq.find_orfs(nested='all')
    orfs2b = seq.find_orfs(need_start='never', nested='all')
    orfs3b = seq.find_orfs(need_start='once', nested='all')
    orfs1c = seq.find_orfs(nested='no')
    orfs2c = seq.find_orfs(need_start='never', nested='no')
    orfs3c = seq.find_orfs(need_start='once', nested='no')
    assert len(orfs1c) < len(orfs1) < len(orfs1b)
    assert len(orfs1c) == 3
    assert len(orfs1) == 4
    assert len(orfs1b) == 7
    assert orfs2 == orfs2b
    assert len(orfs2) == 6
    assert len(orfs2c) == 2
    assert orfs3 == orfs3b
    assert len(orfs3) == 6
    assert len(orfs3c) == 3
    assert orfs1c == orfs1c.copy().sort(lambda ft: [0, 1, 2, -1, -2, -3].index(ft.meta.rf))
    orfs3d = seq.copy().rc().find_orfs(need_start='once', nested='all')
    for orf in orfs3d:
        orf.meta.rf = - 1 - orf.meta.rf
    assert orfs3d.copy().rc(seqlen=len(seq)).sort() == orfs3b.copy().sort()
    teststr3 = teststr[:30] + '-----' + teststr[30:]
    seq = BioSeq(teststr3.replace('_', ''))
    assert seq.find_orfs(rf='fwd', gap='-') + seq.find_orfs(rf='bwd', gap='-') == seq.find_orfs(rf='all', gap='-')
    assert len(orfs1) == len(seq.find_orfs(gap='-'))
    assert len(orfs2) == len(seq.find_orfs(need_start='never', gap='-'))
    assert len(orfs3) == len(seq.find_orfs(need_start='once', gap='-'))
    assert len(orfs1b) == len(seq.find_orfs(nested='all', gap='-'))
    assert len(orfs2b) == len(seq.find_orfs(need_start='never', nested='all', gap='-'))
    assert len(orfs3b) == len(seq.find_orfs(need_start='once', nested='all', gap='-'))
    assert len(orfs1c) == len(seq.find_orfs(nested='no', gap='-'))
    assert len(orfs2c) == len(seq.find_orfs(need_start='never', nested='no', gap='-'))
    assert len(orfs3c) == len(seq.find_orfs(need_start='once', nested='no', gap='-'))


def test_find_orfs_gaps():
    teststr = 'CCC_CAT_GC_____T_GAC_CTA_CCC_CAT_CAT_GTG_ACC_CCC_CTG_ACC_CAT_GCC_CTG_ACC_ATG_CCC_C'
    seq1 = BioSeq(teststr.replace('_', ''))
    seq2 = BioSeq(teststr)
    orfs1 = seq1.find_orfs()
    orfs2 = seq2.find_orfs(gap='_')
    assert len(orfs2) == len(orfs1)
    for o1, o2 in zip(orfs1, orfs2):
        assert o2.meta.rf == o1.meta.rf
        assert seq1[o1] == seq2[o2].str.replace('_', '')
    orfs1 = seq1.find_orfs(need_start='never', need_stop=False)
    with pytest.warns(match='stop position'):
        orfs2 = seq2.find_orfs(gap='_', need_start='never', need_stop=False)
    assert len(orfs2) == len(orfs1)
    for o1, o2 in zip(orfs1, orfs2):
        assert o2.meta.rf == o1.meta.rf
        if len(seq2[o2].str.replace('_', '')) % 3 == 0:
            assert seq1[o1] == seq2[o2].str.replace('_', '')


def test_find_orfs_vs_orffinder():
    teststr = 'CCC_CAT_GCT_GAC_CTA_CCC_CAT_CAT_GTG_ACC_CCC_CTG_ACC_CAT_GCC_CTG_ACC_ATG_CCC_C'
    seq = BioSeq(teststr.replace('_', ''))
    # https://www.bioinformatics.org/sms2/orf_find.html
    # ATG, fwd, 1, 2, 3
    # >ORF number 1 in reading frame 1 on the direct strand extends from base 52 to base 57.
    # ATGCCC
    # >ORF number 1 in reading frame 2 on the direct strand extends from base 5 to base 28.
    # ATGCTGACCTACCCCATCATGTGA
    # >ORF number 2 in reading frame 2 on the direct strand extends from base 41 to base 49.
    # ATGCCCTGA
    # No ORFs were found in reading frame 3.
    orfs = seq.find_orfs(rf='fwd', need_stop=False)
    assert len(orfs) == 3
    assert orfs[0].locs.range == (51, 57)
    assert orfs[0].meta.rf == 0
    assert orfs[1].locs.range == (4, 28)
    assert orfs[1].meta.rf == 1
    assert orfs[2].locs.range == (40, 49)
    assert orfs[2].meta.rf == 1
    assert all(orf.meta.has_start for orf in orfs)
    assert not orfs[0].meta.has_stop
    assert orfs[1].meta.has_stop
    assert orfs[2].meta.has_stop
    orfs2 = seq.find_orfs(rf='fwd', need_stop=True)
    assert len(orfs2) == 2
    assert orfs2 == orfs[1:]
    # same with any codon:
    # >ORF number 1 in reading frame 1 on the direct strand extends from base 1 to base 57.
    # >ORF number 1 in reading frame 2 on the direct strand extends from base 2 to base 28.
    # >ORF number 2 in reading frame 2 on the direct strand extends from base 29 to base 37.
    # >ORF number 3 in reading frame 2 on the direct strand extends from base 38 to base 49.
    # >ORF number 4 in reading frame 2 on the direct strand extends from base 50 to base 58.
    # >ORF number 1 in reading frame 3 on the direct strand extends from base 3 to base 11.
    # >ORF number 2 in reading frame 3 on the direct strand extends from base 12 to base 56.
    orfs_orffinder = [
        [0, 0, 57],
        [1, 1, 28],
        [1, 28, 37],
        [1, 37, 49],
        [1, 49, 58],
        [2, 2, 11],
        [2, 11, 56],
    ]
    orfs = seq.find_orfs(rf='fwd', need_stop=False, need_start='never')
    orfs.sort('rf')
    assert len(orfs) == 7
    for orf, (rf, start, end) in zip(orfs, orfs_orffinder):
        assert orf.locs.range == (start, end)
        assert orf.meta.rf == rf
    # ATG bwd strand
    # No ORFs were found in reading frame 1.
    # >ORF number 1 in reading frame 2 on the reverse strand extends from base 17 to base 46.
    # >ORF number 2 in reading frame 2 on the reverse strand extends from base 53 to base 58.
    # >ORF number 1 in reading frame 3 on the reverse strand extends from base 6 to base 38.
    orfs_orffinder = [
        [-2, 16, 46],
        [-2, 52, 58],
        [-3, 5, 38],
    ]
    orfs = seq.find_orfs(rf='bwd', need_stop=False)
    assert len(orfs) == 3
    for orf, (rf, start, end) in zip(orfs, orfs_orffinder):
        assert orf.meta.rf == rf
        assert orf.loc.start == len(seq) - end
        assert orf.loc.stop == len(seq) - start
    orfs2 = seq.find_orfs(rf='bwd', need_stop=True)
    assert len(orfs2) == 2
    assert orfs2 == orfs[:1] + orfs[-1:]
    # any codon bwd strand
    # >ORF number 1 in reading frame 1 on the reverse strand extends from base 1 to base 57.
    # >ORF number 1 in reading frame 2 on the reverse strand extends from base 2 to base 46.
    # >ORF number 2 in reading frame 2 on the reverse strand extends from base 47 to base 58.
    # >ORF number 1 in reading frame 3 on the reverse strand extends from base 3 to base 38.
    # >ORF number 2 in reading frame 3 on the reverse strand extends from base 39 to base 56.
    orfs_orffinder = [
        [-1, 0, 57],
        [-2, 1, 46],
        [-2, 46, 58],
        [-3, 2, 38],
        [-3, 38, 56],
    ]
    orfs = seq.find_orfs(rf='bwd', need_stop=False, need_start='never')
    orfs.sort('rf', reverse=True)
    assert len(orfs) == 5
    for orf, (rf, start, end) in zip(orfs, orfs_orffinder):
        assert orf.meta.rf == rf
        assert orf.loc.start == len(seq) - end
        assert orf.loc.stop == len(seq) - start


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


def test_select_seqs_warn_not():
    seqs = read()
    with pytest.warns(match='Attribute'):
        seqs2 = seqs.select(notpresentattr=5)
        seqs3 = seqs.select(notpresentattr_noteq=5)
        seqs4 = seqs.select(notpresentattr_ne=5)
    assert len(seqs2) == 0
    assert len(seqs3) == 2
    assert len(seqs4) == 2


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
