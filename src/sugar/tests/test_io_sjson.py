# (C) 2024, Tom Eulenfeld, MIT license

import tempfile

from sugar import read
from sugar.tests.util import _clean_seqs, tempfilename


def test_sjson():
    seqs = _clean_seqs(read())
    with tempfilename(suffix='.sjson') as fname:
        seqs.write(fname)
        seqs2 = _clean_seqs(read(fname))
    assert seqs2 == seqs
    assert not hasattr(seqs, '_fmtcomment')
    assert not hasattr(seqs2, '_fmtcomment')
    seqs2 = _clean_seqs(seqs.fromfmtstr(seqs.tofmtstr('sjson')))
    assert seqs2 == seqs
    assert not hasattr(seqs, '_fmtcomment')
    assert not hasattr(seqs2, '_fmtcomment')
