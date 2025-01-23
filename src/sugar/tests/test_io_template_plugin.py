# (C) 2024, Tom Eulenfeld, MIT license

import sugar
from sugar.tests.util import tempfilename


### The following lines are copied from the sugar wiki
### How to create a new sequence or feature plugin
"""
My Fancy Seq Plugin

This is an example sequence file format. The layout is as follows:
#MyFancySeqFormat
seq1 AAATTGGGCCC
seq2 ATGGCT
"""

# To add read support you must define either an iter_fancy or read_fancy function, or both.
# To add write support you must define either an append_fancy or write_fancy function, or both.

from sugar import BioBasket, BioSeq
from sugar._io.util import _add_fmt_doc

# Use the following flag to indicate, that your file format is binary rather than text-based,
# the passed file handlers will be opened in binary mode.
#binary_fmt = True

# optional, filename extensions for automatic detection of file format
# when writing
filename_extensions_fancy = ['fancy']

def is_fancy(f, **kw):
    """
    Function is optional, used for auto-detection of format when reading

    It should return True if the format is detected,
    otherwise it may raise any exception or return False.
    """
    content = f.read(50)
    return content.strip().lower().startswith('#myfancyseqformat')

# The function decorators are used to automatically add a warning
# to the docstring, that this function should be called via the main
# iter_ or read functions.
@_add_fmt_doc('read')
def iter_fancy(f, optional_argument=None):
    """
    The iter_fancy function expects a file handler and has to yield BioSeq objects.

    You can define optional arguments.
    """
    for line in f:
        if line.strip() != '' and not line.startswith('#'):
            seqid, data = line.split()
            yield BioSeq(data, id=seqid)

@_add_fmt_doc('read')
def read_fancy(f, **kw):
    """
    The read_fancy function expects a file handler and has to return a BioBasket object
    """
    # We are lazy here and reuse iter_fancy
    return BioBasket(list(iter_fancy(f, **kw)))

@_add_fmt_doc('write')
def append_fancy(seq, f, **kw):
    """
    Write a single seq to file handler
    """
    f.write(f'{seq.id} {seq.data}\n')

@_add_fmt_doc('write')
def write_fancy(seqs, f, **kw):
    """
    Write a BioBasket object to file handler
    """
    f.write('#MyFancySeqFormat 3.14159\n')
    for seq in seqs:
        # be lazy again
        append_fancy(seq, f, **kw)
### End of copied lines


def test_create_template_seq_test_file():
    # from importlib.resources import files
    # fname = str(files('sugar.tests.data').joinpath('io_template_plugin.fancy'))
    seqs = sugar.read()
    with tempfilename(suffix='.fancy') as fname:
        seqs.write(fname)


def test_template_seq_plugin():
    seq1 = sugar.read('!data/io_template_plugin.fancy')[0]
    seq2 = sugar.read()[0]
    assert seq1.id == seq2.id
    assert str(seq1) == str(seq2)
