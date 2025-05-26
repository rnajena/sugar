# (C) 2025, Tom Eulenfeld, MIT license

from importlib.resources import files
import pytest
from sugar import BioBasket, BioSeq, Feature, read, read_fts
from sugar.tests.util import _clean_fts, _clean_seqs


def test_seq2biopython2seq():
    SeqRecord = pytest.importorskip('Bio.SeqRecord', reason='require biopython module')
    seqs = read()
    seq = seqs[0]
    obj = seq.tobiopython()
    seq2 = seq.frombiopython(obj)
    assert isinstance(obj, SeqRecord.SeqRecord)
    assert seq2.id == obj.id == seq.id
    assert str(seq2) == str(obj.seq) == str(seq)
    assert _clean_seqs([seq])[0] == seq2


def test_seqs2biopython2seqs():
    SeqRecord = pytest.importorskip('Bio.SeqRecord', reason='require biopython module')
    seqs = read()
    obj = seqs.tobiopython()
    seqs2 = seqs.frombiopython(obj)
    assert isinstance(obj[0], SeqRecord.SeqRecord)
    assert seqs2[0].id == obj[0].id == seqs[0].id
    assert str(seqs2[0]) == str(obj[0].seq) == str(seqs[0])
    assert _clean_seqs(seqs) == seqs2


def test_biopython2seqs2biopython():
    SeqIO = pytest.importorskip('Bio.SeqIO', reason='require biopython module')
    fname = files('sugar.tests.data').joinpath('example.gb')
    with open(fname) as f:
        records = list(SeqIO.parse(f, 'genbank'))
    seqs = BioBasket.frombiopython(records)
    assert [len(seq.fts) for seq in seqs] == [len(seq.features) for seq in records]
    records2 = seqs.tobiopython()
    for attr in 'seq id name description dbxrefs annotations'.split():
        assert getattr(records2[0], attr) == getattr(records[0], attr)
    from Bio.SeqRecord import Seq, SeqRecord
    record = SeqRecord(Seq('ATG'))
    seq = BioSeq.frombiopython(record)
    assert seq.id is None
    assert '_biopython' not in seq.meta
    assert seq.data == 'ATG'
    record2 = seq.tobiopython()
    for attr in 'seq id name description dbxrefs annotations'.split():
        assert getattr(record2, attr) == getattr(record, attr)


def test_seqs2biopython2seqs_msa():
    SeqRecord = pytest.importorskip('Bio.SeqRecord', reason='require biopython module')
    seqs = read().str.ljust(10_000, '-')
    obj = seqs.tobiopython(msa=True)
    seqs2 = seqs.frombiopython(obj)
    assert isinstance(obj[0], SeqRecord.SeqRecord)
    assert seqs2[0].id == obj[0].id == seqs[0].id
    assert str(seqs2[0]) == str(obj[0].seq) == str(seqs[0])
    assert _clean_seqs(seqs) == seqs2


def test_fts2biopython2fts():
    SeqFeature = pytest.importorskip('Bio.SeqFeature', reason='require biopython module')
    from sugar.core._adapter import fts2biopython, biopython2fts
    fts = read_fts()
    obj = fts2biopython(fts)
    fts2 = biopython2fts(obj)
    assert isinstance(obj[0], SeqFeature.SeqFeature)
    assert _clean_fts(fts) == _clean_fts(fts2)


def test_ft2biopython2ft():
    SeqFeature = pytest.importorskip('Bio.SeqFeature', reason='require biopython module')
    pass
    from sugar.core._adapter import fts2biopython, biopython2fts
    fts = read_fts()
    fts2 = fts2biopython(fts)
    fts3 = biopython2fts(fts2)
    assert _clean_fts(fts) == _clean_fts(fts3)


def test_biopython2ft2biopython():
    SeqFeature = pytest.importorskip('Bio.SeqFeature', reason='require biopython module')
    bioft = SeqFeature.SeqFeature(SeqFeature.SimpleLocation(1, 10, -1))
    ft = Feature.frombiopython(bioft)
    assert ft.id is None
    assert ft.seqid is None
    assert '_biopython' not in ft.meta
    bioft2 = ft.tobiopython()
    assert bioft2 == bioft


def test_bioseq_to_from_biotite():
    pytest.importorskip('biotite', reason='require biotite module')
    seqs = read()
    seq = seqs[0]
    obj = seq.tobiotite()
    seq2 = seq.frombiotite(obj)
    assert str(seq2) == ''.join(obj.symbols) == str(seq)


def test_bioseq_toftsviewer():
    featview = pytest.importorskip('dna_features_viewer', reason='require dna_features_viewer module')
    seqs = read()
    seq = seqs[0]
    obj = seq.toftsviewer()
    assert isinstance(obj, featview.GraphicRecord)


def test_seqs2biotite2seqs():
    pytest.importorskip('biotite', reason='require biotite module')
    seqs = read()
    obj = seqs.tobiotite()
    seqs2 = seqs.frombiotite(obj)
    assert str(seqs2[0]) == ''.join(obj[0].symbols) == str(seqs[0])
    seqs.str.rjust(10_000, '-')
    with pytest.warns():
        obj = seqs.tobiotite()
    seqs2 = seqs.frombiotite(obj)
    assert str(seqs2[0]) == ''.join(obj[0].symbols) == str(seqs[0].str.replace('-', ''))


def test_seqs2biotite2seqs_msa():
    pytest.importorskip('biotite', reason='require biotite module')
    seqs = read().str.rjust(10_000, '-')
    obj = seqs.tobiotite(msa=True)
    seqs2 = seqs.frombiotite(obj)
    assert str(seqs2[0]) == str(seqs[0])
