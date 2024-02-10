"""
sugar

Simple containers for RNA and DNA sequence analysis
"""
from sugar.core.meta import Attr, Meta
from sugar.core.fts import Feature, FeatureList, Location
from sugar.core.seq import BioBasket, BioSeq
from sugar.io import read, iter_, read_fts
from sugar import data, web

__version__ = '0.0.1-dev'
