# (C) 2026, Tom Eulenfeld, MIT license
from sugar import BioBasket, read_fts

# Sequences used to generate test files with MEME suite at https://meme-suite.org
# mode anr, nmotifs 3, minw 3
# command: meme sequences.fa -dna -oc . -nostatus -time 14400 -mod anr -nmotifs 3 -minw 3 -maxw 50 -objfun classic -minsites 2 -revcomp -markov_order 0
DATA = """
>seq1
ACGTTTGGTTATTTTACCAACTTTAAC
>seq2
TCGTTTCGGTTACCAACTTTAAC
>seq3
GGTAACCAAACGC
"""

MOTIFS = """
CGTTTGGT
CGTTTGGT
AAGTTGGT
AAGTTGGT
CGTTTCGG
AAC
AAC
TTT
TAT
"""


def test_meme_txt():
    fts = read_fts('!data/fts_example_meme.txt', ftype='motif')
    assert len(fts) == 9
    assert fts[0].loc.start == 4
    assert fts[0].loc.strand == '-'
    assert fts[0].seqid == 'seq3'
    assert fts[0].meta._fmt == 'meme_txt'
    fts2 = read_fts('!data/fts_example_meme.txt', engine='txt', ftype='motif')
    assert fts2 == fts

    seqs = BioBasket.fromfmtstr(DATA)
    assert set(map(str, seqs[fts])) == set(MOTIFS.split())
    seqs.fts = fts
    assert len(seqs.fts) == len(fts)
