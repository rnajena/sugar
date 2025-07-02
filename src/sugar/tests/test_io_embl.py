# (C) 2025, Tom Eulenfeld, MIT license
from sugar import read, read_fts


def test_embl():
    seqs1 = read()
    seqs2 = read('!data/example_ena.embl')
    seqs3 = read('!data/example_uniprot.embl')
    assert seqs2[0].data == seqs1[0].data
    assert seqs2[0].id == seqs1[0].id
    assert len(set(seqs2[0].meta._embl.keys()) & set(seqs1[0].meta._genbank.keys())) > 0
    assert seqs2[0].meta._embl.accession.rstrip(';') == seqs1[0].meta._genbank.accession
    assert seqs3[0].data == seqs1[0].fts.get('CDS').meta._genbank.translation
    assert seqs3[0].data == seqs2[0].fts.get('CDS').meta._embl.translation


def test_embl_fts():
    fts1 = read_fts()
    fts2 = read_fts('!data/example_ena.embl')
    fts3 = read_fts('!data/example_uniprot.embl')
    assert len(fts2) == 2
    assert len(fts3) == 255
