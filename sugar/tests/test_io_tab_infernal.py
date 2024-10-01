# (C) 2024, Tom Eulenfeld, MIT license
from importlib.resources import files
from sugar import read_fts


def test_infernal():
    fts = read_fts('!data/fts_example.infernal', ftype='testq')
    assert fts[0].loc.start == 361955 - 1
    assert fts[0].loc.strand == '-'
    assert fts[0].type == 'testq'
    assert fts[0].meta.seqid == 'NC_013790.1'
    assert fts[2].loc.strand == '+'
    assert len(fts) == 56
    assert fts[0].meta._fmt == 'infernal'


def test_infernal_different_fmt_cmsearch():
    fts1 = read_fts('!data/fts_example.infernal')  # fmt 1
    fts2 = read_fts('!data/fts_example_tab_infernal_cmsearch_fmt3.txt')
    assert len(fts2) == len(fts1)


def test_infernal_different_fmt_cmscan():
    fts1 = read_fts('!data/fts_example_tab_infernal_cmscan_fmt1.txt')
    fts2 = read_fts('!data/fts_example_tab_infernal_cmscan_fmt2.txt')
    fts3 = read_fts('!data/fts_example_tab_infernal_cmscan_fmt3.txt')
    assert len(fts1) == 4
    assert len(fts2) == 4
    assert len(fts3) == 4


def test_infernal_different_fmt_cmscan_overlap():
    fts1 = read_fts('!data/fts_example_tab_infernal_cmscan_fmt2_mrum.txt')
    assert len(fts1) == 25
