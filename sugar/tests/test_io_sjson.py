# (C) 2024, Tom Eulenfeld, MIT license

from sugar.io.main import read
import tempfile

def _clean_meta(seqs):
    for seq in seqs:
        for key in list(seq.meta):
            if key.startswith('_'):
                delattr(seq.meta, key)
    return seqs


def test_sjson():
    seqs = _clean_meta(read())
    with tempfile.NamedTemporaryFile(suffix='.sjson', delete=False) as f:
        f.close()
        seqs.write(f.name)
        seqs2 = _clean_meta(read(f.name))
    assert seqs2 == seqs
    assert not hasattr(seqs, '_fmtcomment')
    assert not hasattr(seqs2, '_fmtcomment')
    seqs2 = _clean_meta(seqs.fromfmtstr(seqs.tofmtstr('sjson')))
    assert seqs2 == seqs
    assert not hasattr(seqs, '_fmtcomment')
    assert not hasattr(seqs2, '_fmtcomment')
