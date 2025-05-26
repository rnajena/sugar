# (C) 2024, Tom Eulenfeld, MIT license
import pytest
from sugar import read, iter_, read_fts
from sugar._io.genbank import _parse_locs


def test_gb():
    fname = '!data/example.gb'
    seqs = read(fname)
    for i, seq in enumerate(iter_(fname)):
        assert seq == seqs[i]
    seqs2 = read(fname, exclude=('translation',))
    assert sum('translation' in f.meta._genbank for f in seqs[0].fts) == 1
    assert sum('translation' in f.meta._genbank for f in seqs2[0].fts) == 0
    assert seqs[0].fts.get('source').meta._genbank.misc == ['test', 'test2']


def test_gb_fts():
    fname = '!data/example.gb'
    fts = read_fts(fname)
    fts2 = read_fts(fname, exclude=('translation',))
    assert sum('translation' in f.meta._genbank for f in fts) == 2
    assert sum('translation' in f.meta._genbank for f in fts2) == 0
    assert fts.get('source').meta._genbank.misc == ['test', 'test2']


def test_gb_loc_parsing():
    """
    Examples from https://www.insdc.org/submitting-standards/feature-table/#3.4
    """
    locs = _parse_locs('467')
    assert len(locs) == 1
    assert len(locs[0]) == 1
    assert locs[0].start == 466
    assert locs[0].strand == '+'
    locs = _parse_locs('340..565')
    assert len(locs) == 1
    assert locs[0].start == 339
    assert locs[0].strand == '+'
    locs = _parse_locs('<345..500')
    assert locs[0].defect == locs[0].defect.BEYOND_LEFT
    locs = _parse_locs('<1..888')
    assert locs[0].defect == locs[0].defect.BEYOND_LEFT
    locs = _parse_locs('1..>888')
    assert locs[0].defect == locs[0].defect.BEYOND_RIGHT
    locs = _parse_locs('102.110')
    assert locs[0].defect == locs[0].defect.UNKNOWN_SINGLE_BETWEEN
    locs = _parse_locs('123^124')
    assert locs[0].defect == locs[0].defect.BETWEEN_CONSECUTIVE
    locs = _parse_locs('join(12..78,134..202)')
    assert len(locs) == 2
    assert locs[0].start == 11
    assert locs[1].start == 133
    locs = _parse_locs('complement(34..126)')
    assert locs[0].start == 33
    assert locs[0].strand == '-'
    locs = _parse_locs('complement(join(2691..4571,4918..5163))')
    assert len(locs) == 2
    assert locs[0].start == 4917
    assert locs[0].strand == '-'
    assert locs[1].start == 2690
    locs2 = _parse_locs('join(complement(4918..5163),complement(2691..4571))')
    assert locs2 == locs
    locs2 = _parse_locs('order(complement(4918..5163),complement(2691..4571))')
    assert locs2 == locs
    with pytest.warns():
        locs = _parse_locs('J00194.1:100..202')
        locs[0].start == 99
        assert locs[0].meta._genbank.seqid == 'J00194.1'
        locs2 = _parse_locs('join(1..100,J00194.1:100..202)')
        locs2[0].start == 0
        locs2[1].start == 0
        assert locs2[1] == locs[0]
