# (C) 2024, Tom Eulenfeld, MIT license
from importlib.resources import files
from sugar import read_fts


def test_blast():
    fts = read_fts('!data/fts_example.blastn', ftype='testq')
    assert fts[0].loc.start == 332
    assert fts[0].loc.strand == '+'
    assert fts[0].type == 'testq'
    assert fts[0].meta.seqid == 'tests1'
    assert fts[1].loc.strand == '-'


def test_blast_outfmt():
    fname = str(files('sugar.tests.data').joinpath('fts_example.blastn'))
    with open(fname) as f:
        while f.readline().startswith('#'):
            pass
        outfmt = ('qseqid sseqid bitscore evalue pident length mismatch '
                  'gapopen qstart qend qlen sstart send sstrand slen')
        fts = read_fts(f, outfmt=outfmt)
    assert len(fts) == 1
