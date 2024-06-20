"""
Welcome
=======

Welcome to sugar's API documentation!

The sugar project aims to provide a lightweight framework
to facilitate rapid application development for bioinformatics.

It thus provides classes and functions to deal with DNA and RNA sequences and
related annotations and it provides parsers for common file formats using a
plugin interface.

The core of sugar are the sequence handling classes `.BioSeq` and `.BioBasket`.
The `.BioSeq` class behaves like a string with useful bioinformatics methods attached.
The `.BioBasket` class is a container for several `.BioSeq` objects and
behaves like a list with useful methods attached.
An example of such a useful method is certainly the `~.BioBasket.translate()` method.

To read sequences, use the powerful `~._io.main.read()` routine.
It can handle glob expressions, web resources, archives,
and automatically detects file formats by inspecting the file contents.
To write sequences to files, use the `.BioBasket.write()` method.

Since sugar uses a plugin system, it is easy to add support for new file formats on the fly.
The following sequence file formats are supported out of the box:

{format_table_seqs}

The table links to the used modules and function.

Sequence features, respective annotations,
can be handled with the `.Feature` and `.FeatureList` classes.
To read features use the `~._io.main.read_fts()` routine and
to write features use the `.FeatureList.write()` method.
The following feature formats are supported out of the box:

{format_table_fts}

sugar provides an indexing tool to quickly retrieve
sequences or subsequences from large FASTA files.
Use it via the `.FastaIndex` class or the ``sugar index`` command.

The `.Entrez` class can be used to fetch sequences
from the NCBI online database.

sugar also provides access to codon translation tables and substitution
matrices (i.e. BLOSUM62) with the `sugar.data` package.

The command line interface of sugar can be used for common tasks.
Call ``sugar -h`` for an overview of available commands.
The test suite can be run with ``sugar test``.
"""

__version__ = '0.2.0'

from sugar.core.meta import Attr, Meta
from sugar.core.fts import Feature, FeatureList, Location
from sugar.core.seq import BioBasket, BioSeq
from sugar._io import read, iter_, read_fts
from sugar.index.fastaindex import FastaIndex
from sugar import data, web


import os as __os
from sugar._io.util import _insert_format_plugin_table as __insert_table
from sugar._io.util import _create_format_plugin_table as __create_table
__sphinx = __os.environ.get('SPHINX_BUILD') == '1'
__insert_table('seqs', 'in', fancy=__sphinx)(read)
__insert_table('seqs', 'out', fancy=__sphinx)(BioBasket.write)
__insert_table('fts', 'in', fancy=__sphinx)(read_fts)
__insert_table('fts', 'out', fancy=__sphinx)(FeatureList.write)
if __doc__ is not None:
    _original_doc = __doc__
    __doc__ = __doc__.format(
        format_table_seqs=__create_table('seqs', 'io', fancy=__sphinx),
        format_table_fts=__create_table('fts', 'io', fancy=__sphinx)
        )
