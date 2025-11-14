"""
Sugar -- A Python framework for bioinformatics
"""

__version__ = '1.0.0'

from sugar.core.meta import Attr, Meta
from sugar.core.fts import Feature, FeatureList, Location
from sugar.core.seq import BioBasket, BioSeq
from sugar._io import read, iter_, read_fts, write as __write, write_fts as __write_fts
from sugar.index.fastaindex import FastaIndex
from sugar import data, web


import os as __os
from sugar._io.util import _insert_format_plugin_table as __insert_table
__sphinx = __os.environ.get('SPHINX_BUILD') == '1'
__insert_table('seqs', 'in', fancy=__sphinx)(read)
__insert_table('seqs', 'out', fancy=__sphinx)(__write)
__insert_table('fts', 'in', fancy=__sphinx)(read_fts)
__insert_table('fts', 'out', fancy=__sphinx)(__write_fts)
