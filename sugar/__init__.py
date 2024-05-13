"""
sugar

Simple containers for RNA and DNA sequence analysis
"""

__version__ = '0.0.1-dev'

from sugar.core.meta import Attr, Meta
from sugar.core.fts import Feature, FeatureList, Location
from sugar.core.seq import BioBasket, BioSeq
from sugar.io import read, iter_, read_fts
from sugar.io.fastaindex import FastaIndex
from sugar import data, web
