# (C) 2023, Tom Eulenfeld, MIT license

from contextlib import contextmanager
from importlib.resources import files
from pathlib import Path
import tempfile
import os


@contextmanager
def tempfilename(suffix=''):
    """
    A context manager yielding the name of a temporary file,

    use this or the pytest fixture ``tmpfname``
    """
    with tempfile.TemporaryDirectory(prefix='sugar_') as tmp:
        yield os.path.join(tmp, 'tmp' + suffix)


def normalize_content(fname, ignore=()):
    if fname.startswith('!data/'):
        fname = fname.removeprefix('!data/')
        fname = str(files('sugar.tests.data').joinpath(fname))
    with open(fname) as f:
        return ' '.join([l for line in f.read().strip().split('\n')
                           if not any(line.startswith(pre) for pre in ignore)
                           for l in line.strip().split()])


def _clean_fts(fts):
    for ft in fts:
        for key in list(ft.meta):
            if key.startswith('_'):
                delattr(ft.meta, key)
        for loc in ft.locs:
            for key in list(loc.meta):
                if key.startswith('_'):
                    delattr(loc.meta, key)
    return fts


def _clean_seqs(seqs):
    for seq in seqs:
        for key in list(seq.meta):
            if key.startswith('_'):
                delattr(seq.meta, key)
        if fts := seq.meta.get('fts'):
           _clean_fts(fts)
    return seqs


# @contextmanager
# def _changetmpdir(path=None):
#     if path is None:
#         origin = Path().resolve()
#         with tempfile.TemporaryDirectory() as tmpdir:
#             try:
#                 os.chdir(tmpdir)
#                 yield Path(tmpdir)
#             finally:
#                 os.chdir(origin)
#     else:
#         yield Path(path)
