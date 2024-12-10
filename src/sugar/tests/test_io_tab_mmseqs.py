# (C) 2024, Tom Eulenfeld, MIT license
from importlib.resources import files
from sugar import read_fts
from sugar._io.tab.core import _CONVERTH, _DEFAULT_OUTFMT


def test_mmseqs():
    fts = read_fts('!data/fts_example.mmseqs2', ftype='testq')
    assert fts[0].loc.start == 39922089 - 1
    assert fts[0].loc.strand == '-'
    assert fts[0].type == 'testq'
    assert fts[0].meta.seqid == 'NC_081844.1'
    assert fts[1].loc.strand == '+'
    assert len(fts) == 4
    assert fts[0].meta._fmt == 'mmseqs'


def test_mmseqs_outfmt():
    fname = str(files('sugar.tests.data').joinpath('fts_example.mmseqs2'))
    with open(fname) as f:
        f.readline()
        fts = read_fts(f, outfmt=' '.join(_DEFAULT_OUTFMT['mmseqs']), fmt='mmseqs')
    assert len(fts) == 4


def test_mmseqs_different_outfmt():
    fts1 = read_fts('!data/fts_example.mmseqs2')  # fmtmode 4
    fts2 = read_fts('!data/fts_example_blast_mmseqs2_fmtmode0.txt')
    assert len(fts1) == 4
    assert len(fts2) == 4
    assert fts2 == fts1


def test_mmseqs_different_outfmt_all_headers():
    outfmt = ('query,target,evalue,gapopen,pident,fident,nident,qstart,qend,qlen,tstart,tend,tlen,alnlen,raw,'
              'bits,cigar,qseq,qheader,theader,qaln,taln,qframe,tframe,mismatch,qcov,tcov,qset,qsetid,tset,'
              'tsetid,qorfstart,qorfend,torfstart,torfend,ppos').replace(',', ' ')
    fts1 = read_fts('!data/fts_example_blast_mmseqs2_fmtmode4_all.txt')
    fts2 = read_fts('!data/io_blast_mmseqs2_fmtmode0_all.txt', 'mmseqs', outfmt=outfmt)
    assert len(fts1) == 4
    assert len(fts2) == 4
    assert fts2 == fts1


def test_mmseqs_vs_blast():
    fts1 = read_fts('!data/fts_example.blastn')
    fts2 = read_fts('!data/fts_example.mmseqs2')
    assert len(fts2) == len(fts1)
    assert fts2[0].locs == fts1[0].locs
    for key in fts1[0].meta._blast:
        if key not in ('qstart', 'qend', 'sstart', 'send', 'bitscore'):
            assert fts2[0].meta._mmseqs[_CONVERTH['mmseqs'].get(key, key)] == fts1[0].meta._blast[key]
