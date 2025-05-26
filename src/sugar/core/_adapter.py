# (C) 2025, Tom Eulenfeld, MIT license
# The functions in this module are documented in the corresponding classes
"""
Adapters for other bio libraries

The functions in this module are not meant to be called directly.
Rather, they should be called via the corresponding sugar class methods.
"""

from copy import copy
from sugar.core.fts import Feature, FeatureList, Location
from sugar.core.meta import Attr
from sugar.core.seq import BioSeq, BioBasket


### BioPython sequences

def seq2biopython(seq):
    """
    See `.BioSeq.tobiopython()`
    """
    from Bio.Seq import Seq
    from Bio.SeqRecord import SeqRecord
    features = fts2biopython(seq.fts) if len(seq.fts) > 0 else None
    kw = {k: dict(v) if isinstance(v, Attr) else v for k, v in seq.meta.get('_biopython', {}).items()}
    if seq.id:
        kw['id'] = seq.id
    return SeqRecord(Seq(seq.data), features=features, **kw)


def seqs2biopython(seqs, msa=False):
    """
    See `.BioBasket.tobiopython()`
    """
    seqs = [seq2biopython(seq) for seq in seqs]
    if msa:
        from Bio.Align import MultipleSeqAlignment
        seqs = MultipleSeqAlignment(seqs)
    return seqs


def biopython2seq(obj, cls=BioSeq):
    """
    See `.BioSeq.frombiopython()`
    """
    if hasattr(obj, 'seq'):  # SeqRecord
        data = str(obj.seq)
        id_ = obj.id if obj.id != '<unknown id>' else None
        biopy = {}
        if obj.name != '<unknown name>':
            biopy['name'] = obj.name
        if obj.description != '<unknown description>':
            biopy['description'] = obj.description
        if obj.annotations:
            biopy['annotations'] = obj.annotations
        if obj.annotations:
            biopy['annotations'] = obj.annotations
        if obj.dbxrefs:
            biopy['dbxrefs'] = obj.dbxrefs
        if obj.letter_annotations:
            biopy['letter_annotations'] = obj.letter_annotations
        meta = {'id': id_}
        if biopy:
            meta['_biopython'] = biopy
        if obj.features:
            meta['fts'] = biopython2fts(obj.features)
    else:  # Seq
        data = str(obj)
        meta = None
    return cls(data, meta=meta)


def biopython2seqs(obj, cls=BioBasket):
    """
    See `.BioBasket.frombiopython()`
    """
    seqs = [biopython2seq(seq) for seq in obj]
    return cls(seqs)


### BioPython features

_biopy_strand = {'+': 1, '-': -1, '.': None, '?': None}
_biopy_strand_r = {1: '+', -1: '-', None: '.'}


def ft2biopython(ft):
    """
    See `.Feature.tobiopython()`
    """
    from Bio.SeqFeature import CompoundLocation, SeqFeature, SimpleLocation
    # we ignore any defects, for now
    locs = [SimpleLocation(loc.start, loc.stop, _biopy_strand[loc.strand], ref=ft.seqid) for loc in ft.locs]
    loc = locs[0] if len(locs) == 1 else CompoundLocation(locs)
    kw = {attr: getattr(ft, attr) for attr in 'id type'.split() if getattr(ft, attr)}
    if '_bioypthon' in ft.meta or 'name' in ft.meta:
        kw['qualifiers'] = dict(ft.meta.get('_biopython', {}))
        if 'name' in ft.meta:
            kw['qualifiers']['name'] = ft.meta.name
    biopyft = SeqFeature(location=loc, **kw)
    return biopyft


def fts2biopython(fts):
    """
    See `.FeatureList.tobiopython()`
    """
    return [ft2biopython(ft) for ft in fts]


def biopython2ft(obj, cls=Feature):
    """
    See `.Feature.frombiopython()`
    """
    biopylocs = obj.location.parts
    locs = [Location(int(loc.start), int(loc.end), _biopy_strand_r[loc.strand]) for loc in biopylocs]
    meta = {}
    if obj.id != '<unknown id>':
        meta['id'] = obj.id
    if obj.qualifiers:
        meta['_biopython'] = obj.qualifiers
        if 'name' in meta['_biopython']:
            meta['name'] = meta['_biopython']['name']
    if biopylocs[0].ref:
        meta['seqid'] = biopylocs[0].ref
    ft = cls(obj.type or None, locs, meta=meta)
    return ft


def biopython2fts(obj, cls=FeatureList):
    """
    See `.FeatureList.frombiopython()`
    """
    return cls([biopython2ft(biopyft) for biopyft in obj])


### Biotite sequences

def seq2biotite(seq, type=None, gap='-', warn=True):
    """
    See `.BioSeq.tobiotite()`
    """
    from biotite.sequence import NucleotideSequence, ProteinSequence
    data = seq.data
    if gap:
        for g in gap:
            if g in data:
                if warn:
                    from warnings import warn
                    warn(f'Remove gap characters {gap} for the conversion to biotite')
                for g in gap:
                    data = data.replace(g, '')
                break
    type = type or seq.type
    cls = {'nt': NucleotideSequence, 'aa': ProteinSequence}[type]
    return cls(data)


def seqs2biotite(oseqs, type=None, msa=False, gap='-', warn=True):
    """
    See `.BioBasket.tobiotite()`
    """
    seqs = [seq.tobiotite(type=type, gap=gap, warn=not msa and warn) for seq in oseqs]
    if msa:
        from biotite.sequence.align import Alignment
        trace = Alignment.trace_from_strings([seq.data for seq in oseqs])
        seqs = Alignment(seqs, trace)
    return seqs


def biotite2seq(obj, cls=BioSeq):
    """
    See `.BioSeq.frombiotite()`
    """
    from biotite.sequence import NucleotideSequence, ProteinSequence
    data = ''.join(obj.symbols)
    type_ = None
    if isinstance(obj, NucleotideSequence):
        type_ = 'nt'
    elif isinstance(obj, ProteinSequence):
        type_ = 'aa'
    return cls(data, type=type_)


def biotite2seqs(obj, cls=BioBasket):
    """
    See `.BioBasket.frombiotite()`
    """
    if hasattr(obj, 'sequences'):  # Alignment object
        types = [BioSeq.frombiotite(seq).type for seq in obj.sequences]
        ali = [BioSeq(data, type=type_) for data, type_ in zip(obj.get_gapped_sequences(), types)]
        return cls(ali)
    else:
        seqs = [BioSeq.frombiotite(seq) for seq in obj]
        return cls(seqs)
