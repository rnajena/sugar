# (C) 2024, Tom Eulenfeld, MIT license
from importlib.resources import files
from sugar import read_fts
from sugar._io.tab.core import _DEFAULT_OUTFMT


def test_blast():
    fts = read_fts('!data/fts_example.blastn', ftype='testq')
    assert fts[0].loc.start == 39922089 - 1
    assert fts[0].loc.strand == '-'
    assert fts[0].type == 'testq'
    assert fts[0].meta.seqid == 'NC_081844.1'
    assert fts[1].loc.strand == '+'
    assert len(fts) == 4
    assert fts[0].meta._fmt == 'blast'


def test_blast_outfmt():
    fname = str(files('sugar.tests.data').joinpath('fts_example.blastn'))
    with open(fname) as f:
        while f.readline().startswith('#'):
            pass
        fts = read_fts(f, outfmt=' '.join(_DEFAULT_OUTFMT['blast']), fmt='blast')
    assert len(fts) == 3


def test_blast_different_outfmt():
    fts1 = read_fts('!data/fts_example.blastn')  # outfmt 7
    fts2 = read_fts('!data/fts_example_blast_outfmt6.blastn')
    fts3 = read_fts('!data/io_blast_outfmt10.blastn', sep=',')
    assert len(fts1) == 4
    assert len(fts2) == 4
    assert len(fts3) == 4
    assert fts2 == fts1
    assert fts3 == fts1


def test_blast_different_outfmt_all_headers():
    outfmt = ('qseqid qgi qacc qaccver qlen sseqid sallseqid sgi sallgi sacc saccver '
              'sallacc slen qstart qend sstart send qseq evalue bitscore score length '
              'pident nident mismatch positive gapopen gaps ppos frames qframe sframe '
              'btop staxid ssciname scomname sblastname sskingdom staxids sscinames '
              'scomnames sskingdoms stitle salltitles sstrand qcovs qcovhsp qcovus')
    fts1 = read_fts('!data/fts_example_blast_outfmt7_all.blastn')
    fts2 = read_fts('!data/io_blast_outfmt6_all.blastn', 'blast', outfmt=outfmt)
    fts3 = read_fts('!data/io_blast_outfmt10_all.blastn', 'blast', sep=',', outfmt=outfmt)
    assert len(fts1) == 4
    assert len(fts2) == 4
    assert len(fts3) == 4
    assert fts2 == fts1
    assert fts3 == fts1
