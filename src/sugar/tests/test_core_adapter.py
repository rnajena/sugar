# (C) 2024, Tom Eulenfeld, MIT license

import pytest
from sugar import read


def test_bioseq_to_from_biopython():
    SeqRecord = pytest.importorskip('Bio.SeqRecord', reason='require biopython module')
    seqs = read()
    seq = seqs[0]
    obj = seq.tobiopython()
    seq2 = seq.frombiopython(obj)
    assert isinstance(obj, SeqRecord.SeqRecord)
    assert seq2.id == obj.id == seq.id
    assert str(seq2) == str(obj.seq) == str(seq)


def test_biobasket_to_from_biopython():
    SeqRecord = pytest.importorskip('Bio.SeqRecord', reason='require biopython module')
    seqs = read()
    obj = seqs.tobiopython()
    seqs2 = seqs.frombiopython(obj)
    assert isinstance(obj[0], SeqRecord.SeqRecord)
    assert seqs2[0].id == obj[0].id == seqs[0].id
    assert str(seqs2[0]) == str(obj[0].seq) == str(seqs[0])


def test_biobasket_to_from_biopython_msa():
    SeqRecord = pytest.importorskip('Bio.SeqRecord', reason='require biopython module')
    seqs = read().str.ljust(10_000, '-')
    obj = seqs.tobiopython(msa=True)
    seqs2 = seqs.frombiopython(obj)
    assert isinstance(obj[0], SeqRecord.SeqRecord)
    assert seqs2[0].id == obj[0].id == seqs[0].id
    assert str(seqs2[0]) == str(obj[0].seq) == str(seqs[0])


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


def test_biobasket_to_from_biotite():
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


def test_biobasket_to_from_biotite_msa():
    pytest.importorskip('biotite', reason='require biotite module')
    seqs = read().str.rjust(10_000, '-')
    obj = seqs.tobiotite(msa=True)
    seqs2 = seqs.frombiotite(obj)
    assert str(seqs2[0]) == str(seqs[0])
