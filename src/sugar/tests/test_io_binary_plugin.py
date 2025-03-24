# (C) 2024, Tom Eulenfeld, MIT license

import sugar


### define a simple FASTA-like binary plugin
from sugar.core.seq import BioBasket, BioSeq
binary_fmt = True
filename_extensions_testbin = ['bintest']
MAGIC = b'\xfe\x8abintestformat'
def is_testbin(f, **kw):
    return f.read(len(MAGIC)) == MAGIC
def read_testbin(f):
    return BioBasket(list(iter_testbin(f)))
def iter_testbin(f):
    assert is_testbin(f)
    while line := f.readline():
        if line.strip() != b'':
            if line.startswith(b'>'):
                id_ = line[1:].decode('latin1').strip()
            else:
                data = line.decode('latin1').strip()
                yield BioSeq(data, id=id_)
def append_testbin(seq, f, magic=True):
    if magic:
        f.write(MAGIC)
    f.write(b'>' + seq.id.encode('latin1'))
    f.write(str(seq).encode('latin1'))
def write_testbin(seqs, f):
    f.write(MAGIC)
    for seq in seqs:
        append_testbin(seq, f, magic=False)
###

### define another fake Feature binary plugin based on gff
binary_fmt_fts = True
filename_extensions_fts_testbin = ['bintest']
MAGIC2 = b'\xfe\x8aftsbintestformat'
def is_fts_testbin(f, **kw):
    return f.read(len(MAGIC2)) == MAGIC2
def read_fts_testbin(f):
    assert is_fts_testbin(f)
    return sugar.read_fts(f)
def write_fts_testbin(fts, f):
    f.write(MAGIC2)
    fts.write(f, fmt='gff')
###


def test_create_binary_seq_test_file(tmpfname):
    # from importlib.resources import files
    # fname = str(files('sugar.tests.data').joinpath('io_binary_plugin_atg.bintest'))
    seqs = sugar.read()
    seqs.write(tmpfname.with_suffix('.bintest'))


def test_create_binary_fts_test_file(tmpfname):
    # from importlib.resources import files
    # fname = str(files('sugar.tests.data').joinpath('io_binary_plugin_fake_fts.bintest'))
    fts = sugar.read_fts()
    fts.write(tmpfname.with_suffix('.bintest'))


def test_binary_seq_plugin():
    seq1 = sugar.read('!data/io_binary_plugin_atg.bintest')[0]
    seq2 = sugar.read()[0]
    assert seq1.id == seq2.id
    assert str(seq1) == str(seq2)


def test_binary_fts_plugin():
    ft1 = sugar.read_fts('!data/io_binary_plugin_fake_fts.bintest')[0]
    ft2 = sugar.read_fts()[0]
    del ft1.meta._fmt
    del ft2.meta._fmt
    assert ft1 == ft2
