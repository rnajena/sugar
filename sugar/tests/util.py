# (C) 2023, Tom Eulenfeld, MIT license

from contextlib import contextmanager
from importlib.resources import files
import tempfile


@contextmanager
def tempfilename(**kw):
    with tempfile.NamedTemporaryFile(delete=False, **kw) as fp:
        fp.close()
        yield fp.name


def normalize_content(fname, ignore=()):
    if fname.startswith('!data/'):
        fname = fname.removeprefix('!data/')
        fname = str(files('sugar.tests.data').joinpath(fname))
    with open(fname) as f:
        return ' '.join([l for line in f.read().strip().split('\n')
                           if not any(line.startswith(pre) for pre in ignore)
                           for l in line.strip().split()])
