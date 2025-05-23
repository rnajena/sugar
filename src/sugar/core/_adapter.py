# (C) 2025, Tom Eulenfeld, MIT license
# The functions in this module are documented in the corresponding classes
"""
Adapters for other bio libraries

The functions in this module are not meant to be called directly.
Rather, they should be called via the corresponding sugar class methods or their instances.
"""


from sugar.core.seq import BioSeq, BioBasket
from sugar.core.util import _add_doc_from_other

### BioPython sequences

@_add_doc_from_other(BioSeq.tobiopython)
def seq2biopython(seq):
    from Bio.Seq import Seq
    from Bio.SeqRecord import SeqRecord
    return SeqRecord(Seq(seq.data), id=seq.id)


@_add_doc_from_other(BioBasket.tobiopython)
def seqs2biopython(seqs, msa=False):
    seqs = [seq.tobiopython() for seq in seqs]
    if msa:
        from Bio.Align import MultipleSeqAlignment
        seqs = MultipleSeqAlignment(seqs)
    return seqs


@_add_doc_from_other(BioSeq.frombiopython)
def biopython2seq(obj, cls=BioSeq):
    if hasattr(obj, 'seq'):  # SeqRecord
        data = str(obj.seq)
        id_ = obj.id
    else:  # Seq
        data = str(obj)
        id_ = None
    return cls(data, id=id_)


@_add_doc_from_other(BioBasket.frombiopython)
def biopython2seqs(obj, cls=BioBasket):
    seqs = [BioSeq.frombiopython(seq) for seq in obj]
    return cls(seqs)


### Biotite sequences

@_add_doc_from_other(BioSeq.tobiotite)
def seq2biotite(seq, type=None, gap='-', warn=True):
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


@_add_doc_from_other(BioBasket.tobiotite)
def seqs2biotite(oseqs, type=None, msa=False, gap='-', warn=True):
    seqs = [seq.tobiotite(type=type, gap=gap, warn=not msa and warn) for seq in oseqs]
    if msa:
        from biotite.sequence.align import Alignment
        trace = Alignment.trace_from_strings([seq.data for seq in oseqs])
        seqs = Alignment(seqs, trace)
    return seqs


@_add_doc_from_other(BioSeq.frombiotite)
def biotite2seq(obj, cls=BioSeq):
    from biotite.sequence import NucleotideSequence, ProteinSequence
    data = ''.join(obj.symbols)
    type_ = None
    if isinstance(obj, NucleotideSequence):
        type_ = 'nt'
    elif isinstance(obj, ProteinSequence):
        type_ = 'aa'
    return cls(data, type=type_)


@_add_doc_from_other(BioBasket.frombiotite)
def biotite2seqs(obj, cls=BioBasket):
    if hasattr(obj, 'sequences'):  # Alignment object
        types = [BioSeq.frombiotite(seq).type for seq in obj.sequences]
        ali = [BioSeq(data, type=type_) for data, type_ in zip(obj.get_gapped_sequences(), types)]
        return cls(ali)
    else:
        seqs = [BioSeq.frombiotite(seq) for seq in obj]
        return cls(seqs)
