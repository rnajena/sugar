# (C) 2024, Tom Eulenfeld, MIT license

from contextlib import contextmanager
from functools import reduce, wraps
import glob
from importlib.metadata import entry_points
from importlib.resources import files
import io
import itertools
import operator
import os.path
from pathlib import PurePath
import shutil
import sys
import tempfile
from urllib.parse import urlparse

from sugar import BioBasket, BioSeq


def _epsname_key(epsname):
    try:
        return FMTS.index(epsname)
    except ValueError:
        return len(FMTS)


EPS = entry_points(group='sugar.io')
FMTS = ['fasta', 'genbank', 'stockholm', 'sjson']
FMTS_ALL = sorted(EPS.names, key=_epsname_key)

ARCHIVE_EXTS = [ext.removeprefix('.') for fmt in shutil.get_unpack_formats()
                for ext in fmt[1]]


@contextmanager
def file_opener(f, mode='r'):
    if isinstance(f, PurePath):
        f = str(f)
    if not isinstance(f, str):
        # not a string - we assume a file-like object
        yield f
    else:
        with open(f, mode=mode) as fh:
            yield fh


def detect(fname):
    with file_opener(fname) as f:
        fpos = f.tell()
        for fmt in FMTS_ALL:
            module = EPS[fmt].load()
            if hasattr(module, 'is_format'):
                f.seek(0)
                try:
                    if module.is_format(f):
                        return fmt
                except Exception:
                    pass
                finally:
                    f.seek(fpos)


def detect_ext(fname):
    try:
        _, ext = os.path.splitext(fname)
        ext = ext.removeprefix('.')
    except Exception:
        return
    for fmt in FMTS_ALL:
        module = EPS[fmt].load()
        if hasattr(module, 'EXT'):
            if ext in module.EXT:
                return fmt

def __get_ext(fname):
    return os.path.split(fname)[1].split('.', 1)[-1]


def resolve_fname(reader):
    @wraps(reader)
    def new_reader(fname=None, *args, archive=None, **kw):
        if fname is None:
            fname = '!data/example.gb'
        elif isinstance(fname, PurePath):
            fname = str(fname)
        if isinstance(fname, str):
            if fname.startswith('!data/'):
                fname = fname.removeprefix('!data/')
                fname = str(files('sugar.tests.data').joinpath(fname))
            if fname == '-':
                fname = io.StringIO(sys.stdin.read())
            elif '://' in fname[:10]:
                import requests
                r = requests.get(fname)
                r.raise_for_status()
                bname = os.path.basename(urlparse(fname).path)
                if (__get_ext(bname) in ARCHIVE_EXTS or
                        archive is not None):
                    with tempfile.NamedTemporaryFile(suffix=bname, delete=False) as f:
                        f.write(r.content)
                        f.close()
                        return new_reader(f.name, *args, archive=archive, **kw)
                else:
                    fname = io.StringIO(r.text)
            elif glob.has_magic(fname):
                fnames = glob.glob(fname, recursive=True)
                if not fnames:
                    raise IOError(f'No file matching glob pattern {fname}')
                seqs = [new_reader(fname, *args, archive=archive, **kw) for fname in fnames]
                if isinstance(seqs[0], BioBasket):  # read was wrapped
                    return reduce(operator.add, seqs)
                else:  # iter_ was wrapped
                    return reduce(itertools.chain, seqs)
            elif (__get_ext(fname) in ARCHIVE_EXTS or
                  archive is not None):
                if archive is True:
                    archive = None
                with tempfile.TemporaryDirectory() as tmpdir:
                    shutil.unpack_archive(fname, tmpdir, archive)
                    globexpr = os.path.join(tmpdir, '**/*.*')
                    return new_reader(globexpr, *args, **kw)
        return reader(fname, *args, **kw)
    return new_reader


@resolve_fname
def iter_(fname, fmt=None, tool=None, mode='r', **kw):
    if fmt is None:
        fmt = detect(fname)
    if fmt is None:
        raise IOError('Format cannot be auto-detected')
    if tool == 'biopython':
        from Bio import SeqIO
        for seq in SeqIO.parse(fname, fmt):
            yield BioSeq.fromobj(seq, 'biopython')
    elif tool:
        raise ValueError(f'Unrecognized tool: {tool}')
    else:
        module = EPS[fmt].load()
        if hasattr(module, 'iter_'):
            with file_opener(fname, mode=mode) as f:
                seqs = module.iter_(f, **kw)
                for seq in seqs:
                    seq.meta._fmt = fmt
                    yield seq
        elif hasattr(module, 'read'):
            with file_opener(fname, mode=mode) as f:
                seqs = module.read(f, **kw)
            for seq in seqs:
                seq.meta._fmt = fmt
                yield seq
        else:
            raise RuntimeError(f'No read support for format {fmt}')


@resolve_fname
def read(fname, fmt=None, mode='r', tool=None, **kw):
    if fmt is None:
        fmt = detect(fname)
    if fmt is None:
        raise IOError('Format cannot be auto-detected')
    if tool == 'biopython':
        from Bio import SeqIO
        seqs = SeqIO.parse(fname, fmt, **kw)
        return BioBasket.fromobj(seqs, 'biopython')
    elif tool:
        raise ValueError(f'Unrecognized tool: {tool}')
    module = EPS[fmt].load()
    if hasattr(module, 'read'):
        with file_opener(fname, mode=mode) as f:
            seqs = module.read(f, **kw)
    elif hasattr(module, 'iter_'):
        with file_opener(fname, mode=mode) as f:
            seqs = list(module.iter_(f, **kw))
    else:
        raise RuntimeError(f'No read support for format {fmt}')
    seqs = BioBasket(seqs)
    for seq in seqs:
        seq.meta._fmt = fmt
    return seqs


def write(seqs, fname, fmt=None, mode='w', tool=None, **kw):
    if fmt is None:
        fmt = detect_ext(fname)
    if fmt is None:
        raise IOError('Format cannot be auto-detected')
    if tool == 'biopython':
        from Bio import SeqIO
        return SeqIO.write(seqs.toobj('biopython'), fname, fmt)
    elif tool:
        raise ValueError(f'Unrecognized tool: {tool}')
    module = EPS[fmt].load()
    if hasattr(module, 'append') and 'a' in mode:
        with file_opener(fname, mode=mode) as f:
            for seq in seqs:
                module.append(seq, f, **kw)
    elif hasattr(module, 'write'):
        with file_opener(fname, mode=mode) as f:
            module.write(seqs, f, **kw)
    elif hasattr(module, 'append') and 'w' in mode:
        with file_opener(fname, mode=mode) as f:
            for seq in seqs:
                module.append(seq, f, **kw)
    else:
        raise RuntimeError(f'No write support for format {fmt}')
