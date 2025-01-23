# (C) 2024, Tom Eulenfeld, MIT license

import tempfile

from sugar import read, read_fts
from sugar.tests.util import _clean_seqs, _clean_fts, tempfilename


def test_sjson():
    seqs = _clean_seqs(read())
    oseqs = seqs.copy()
    with tempfilename(suffix='.sjson') as fname:
        seqs.write(fname)
        seqs2 = _clean_seqs(read(fname))
    assert seqs == oseqs
    assert seqs2 == seqs
    assert not hasattr(seqs, '_fmtcomment')
    assert not hasattr(seqs2, '_fmtcomment')
    seqs2 = _clean_seqs(seqs.fromfmtstr(seqs.tofmtstr('sjson')))
    assert seqs2 == seqs
    assert not hasattr(seqs, '_fmtcomment')
    assert not hasattr(seqs2, '_fmtcomment')


def test_fts_sjson():
    fts = _clean_fts(read_fts())
    ofts = fts.copy()
    with tempfilename(suffix='.sjson') as fname:
        fts.write(fname)
        fts2 = _clean_fts(read_fts(fname))
    assert fts == ofts
    assert fts2 == fts
    assert not hasattr(fts, '_fmtcomment')
    assert not hasattr(fts2, '_fmtcomment')
    fts_out = fts.tofmtstr('sjson', indent=2)
    assert '"_cls": "Feature"' in fts_out


def test_sjson_private():
    seqs = read()
    for seq in seqs:
        del seq.meta._fmt
    # we still have the _genbank metadata
    oseqs = seqs.copy()
    with tempfilename(suffix='.sjson') as fname:
        seqs.write(fname, private=True)
        seqs2 = read(fname)
    for seq in seqs2:
        del seq.meta._fmt
    assert seqs == oseqs
    assert seqs2 == seqs

def test_sjson_write_arbritrary_object():
    from sugar._io.sjson import read_sjson, write_sjson
    seqs = _clean_seqs(read())
    with tempfilename(suffix='.sjson') as fname:
        with open(fname, 'w') as f:
            write_sjson([1, seqs], f)
        with open(fname) as f:
            obj = read_sjson(f)
    assert obj == [1, seqs]
