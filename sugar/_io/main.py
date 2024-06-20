# (C) 2024, Tom Eulenfeld, MIT license

from contextlib import contextmanager
from functools import reduce, wraps
import glob
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

from sugar.core.fts import FeatureList
from sugar.core.seq import BioBasket, BioSeq
from sugar._io.util import ARCHIVE_EXTS, EPS, FMTS_ALL


def _binary(module, what='seqs'):
    assert what in ('seqs', 'fts')
    prop = 'binary_fmt' if what == 'seqs' else 'binary_fmt_fts'
    return hasattr(module, prop) and getattr(module, prop)


@contextmanager
def _file_opener(f, mode='r', binary=False, encoding=None):
    if isinstance(f, str):
        if binary and 'b' not in mode:
            mode = mode + 'b'
        with open(f, mode=mode, encoding=encoding) as fh:
            yield fh
    else:
        # not a string - we assume a file-like object
        if not binary and isinstance(f, io.BufferedIOBase):
            f = _NonClosingTextIOWrapper(f, encoding=encoding)
        yield f

class _NonClosingTextIOWrapper(io.TextIOWrapper):
    def __del__(self):
        try:
            self.detach()
        except Exception:
            pass

def detect(fname, what='seqs', *, encoding=None, **kw):
    """
    Try to detect file format from contents

    :param what: ``'seqs'`` or ``'fts'``
    """
    assert what in ('seqs', 'fts')
    suf = '' if what == 'seqs' else '_fts'
    with _file_opener(fname, binary=True) as f:
        fpos = f.tell()
        for fmt in FMTS_ALL[what]:
            module = EPS[what][fmt].load()
            if hasattr(module, 'is_format' + suf):
                if _binary(module, what) and not isinstance(f, io.BufferedIOBase):
                    continue
                if not _binary(module, what) and isinstance(f, io.BufferedIOBase):
                    f_text_or_binary = _NonClosingTextIOWrapper(f, encoding=encoding)
                else:
                    f_text_or_binary = f
                try:
                    if getattr(module, 'is_format' + suf)(f_text_or_binary, **kw):
                        return fmt
                except Exception:
                    pass
                finally:
                    f.seek(fpos)


# def detect_fts(fname, **kw):
#     with file_opener(fname) as f:
#         fpos = f.tell()
#         for fmt in FMTS_FTS_ALL:
#             module = EPS_FTS[fmt].load()
#             if hasattr(module, 'is_format_fts'):
#                 f.seek(fpos)
#                 try:
#                     if module.is_format_fts(f, **kw):
#                         return fmt
#                 except Exception:
#                     pass
#                 finally:
#                     f.seek(fpos)


def detect_ext(fname, what='seqs'):
    """
    Try to detect file format for writing from extension

    :param what: ``'seqs'`` or ``'fts'``
    """
    assert what in ('seqs', 'fts')
    suf = '' if what == 'seqs' else '_fts'
    try:
        _, ext = os.path.splitext(fname)
        ext = ext.removeprefix('.')
    except Exception:
        return
    for fmt in FMTS_ALL[what]:
        module = EPS[what][fmt].load()
        if hasattr(module, 'filename_extensions' + suf):
            if ext in getattr(module, 'filename_extensions' + suf):
                return fmt


# def detect_ext_fts(fname):
#     try:
#         _, ext = os.path.splitext(fname)
#         ext = ext.removeprefix('.')
#     except Exception:
#         return
#     for fmt in FMTS_FTS_ALL:
#         module = EPS_FTS[fmt].load()
#         if hasattr(module, 'filename_extensions_fts'):
#             if ext in module.filename_extensions_fts:
#                 return fmt



def _resolve_archive(writer):
    @wraps(writer)
    def new_writer(objs, fname, *args, archive=None, **kw):
        if isinstance(fname, PurePath):
            fname = str(fname)
        if archive is None:
            return  writer(objs, fname, *args, **kw)
        elif not isinstance(fname, str):
            msg = 'archive option is only allowed for file names, not file-like objects'
            raise ValueError(msg)
        else:
            if archive is True:
                archive = shutil.get_archive_formats()[0][0]
            with tempfile.TemporaryDirectory() as tmpdir:
                tmpfname = os.path.join(tmpdir, os.path.basename(fname))
                writer(objs, tmpfname, *args, **kw)
                shutil.make_archive(fname, archive, tmpdir)
    return new_writer


def __get_ext(fname):
    return os.path.split(fname)[1].split('.', 1)[-1]



def _resolve_fname(example_fname='!data/example.gb'):
    """
    Decorator, takes filename as string and resolved the filename

    Can also deal with online resources, glob expressions,
    BytesIO and TextIO objects are just passed through.
    """
    def wrapper(reader):
        @wraps(reader)
        def new_reader(fname=None, *args, archive=None, **kw):
            # if isinstance(fname, io.StringIO):
            #     msg = 'fname must be string, Path object or io.BytesIO object'
            #     raise ValueError(msg)
            if isinstance(fname, bytes):
                msg = ('The read function cannot take a bytes object, '
                       'but you can wrap it in an instance of io.BytesIO.')
                raise ValueError(msg)
            if fname is None:  # load example file
                fname = example_fname
            elif isinstance(fname, PurePath):  # comnvert Path object to string
                fname = str(fname)
            if not isinstance(fname, str):  # it is a file-like object, pass through
                fl = fname
            else:  # it's a string
                if fname.startswith('!data/'):  # files from example folder
                    fname = fname.removeprefix('!data/')
                    fname = str(files('sugar.tests.data').joinpath(fname))
                if fname == '-':  # pipe from stdin
                    fl = io.BytesIO(sys.stdin.buffer.read())
                elif '://' in fname[:10]:  # it's a urL
                    import requests
                    r = requests.get(fname)
                    r.raise_for_status()
                    bname = os.path.basename(urlparse(fname).path)
                    if (archive is not None and archive != 'gz' or any(
                            bname.endswith('.' + ext) for ext in ARCHIVE_EXTS)):
                        # download archive and run function again
                        with tempfile.NamedTemporaryFile(suffix=bname, delete=False) as f:
                            f.write(r.content)
                            f.close()
                            return new_reader(f.name, *args, archive=archive, **kw)
                    elif archive == 'gz' or bname.endswith('.gz'):  # uncompress download
                        import gzip
                        fl = io.BytesIO(gzip.decompress(r.content))
                    else:
                        # fl = io.StringIO(r.text)  # download is just data
                        fl = io.BytesIO(r.content)  # download is just data
                elif glob.has_magic(fname):  # it's a glob expression
                    fnames = glob.glob(fname, recursive=True)
                    if not fnames:
                        raise IOError(f'No file matching glob pattern {fname}')
                    # run function with all individual files
                    objs = [new_reader(fname, *args, archive=archive, **kw) for fname in fnames]
                    if isinstance(objs[0], (BioBasket, FeatureList)):  # read or read_fts was wrapped
                        return reduce(operator.add, objs)
                    else:  # iter_ was wrapped
                        return reduce(itertools.chain, objs)
                elif (archive is not None and archive != 'gz' or any(
                        fname.endswith('.' + ext) for ext in ARCHIVE_EXTS)):
                    # unpack archive and run again with glob expr
                    if archive is True:
                        archive = None
                    with tempfile.TemporaryDirectory() as tmpdir:
                        shutil.unpack_archive(fname, tmpdir, archive)
                        globexpr = os.path.join(tmpdir, '**/*.*')
                        return new_reader(globexpr, *args, **kw)
                elif archive == 'gz' or fname.endswith('.gz'):  # uncompress file
                    import gzip
                    with gzip.open(fname) as f:
                        fl = io.BytesIO(f.read())
                else:  # fname is just a simple filename, nothing to see, go further
                    fl = fname
            return reader(fl, *args, **kw)
        return new_reader
    return wrapper


@_resolve_fname()
def iter_(fname, fmt=None, *, mode='r', encoding=None, tool=None, **kw):
    """
    Iterate over a file and yield `.BioSeq` objects of each sequence

    See `read()` function.

    Example
    -------

    >>> from sugar import iter_
    >>> for seq in iter_():  # use the example file
    ...     print(f'GC content of seq {seq.id} is {100*seq.gc:.0f}%.')
    GC content of seq AB047639 is 58%.
    GC content of seq AB677533 is 57%.
    """
    if fmt is None:
        fmt = detect(fname, encoding=encoding, **kw)
    if fmt is None:
        raise IOError('Format cannot be auto-detected')
    fmt = fmt.lower()
    if tool == 'biopython':
        from Bio import SeqIO
        for seq in SeqIO.parse(fname, fmt):
            yield BioSeq.fromobj(seq, 'biopython')
    elif tool:
        raise ValueError(f'Unrecognized tool: {tool}')
    else:
        module = EPS['seqs'][fmt].load()
        with _file_opener(fname, mode=mode, binary=_binary(module), encoding=encoding) as f:
            if hasattr(module, 'iter_'):
                seqs = module.iter_(f, **kw)
                for seq in seqs:
                    seq.meta._fmt = fmt
                    yield seq
            elif hasattr(module, 'read'):
                seqs = module.read(f, **kw)
            else:
                raise RuntimeError(f'No read support for format {fmt}')
        for seq in seqs:
            seq.meta._fmt = fmt
            yield seq


@_resolve_fname()
def read(fname, fmt=None, *, mode='r', encoding=None, tool=None, **kw):
    """
    Read a file or file-like object with sequences into `.BioBasket`

    :param fname: Filename, can also be a glob expression,
        a web resource,
        an archive, gzipped file,
        or a file-like object (e.g. `~io.BytesIO`, `~io.StringIO`)
    :param fmt: format of the file (defaul: auto-detect)
    :param mode: mode for opening the file, change this only if you know what
        you do
    :param encoding: encoding of the file
    :param tool: use alternative tool for reading the file,
        supported tools are: ``'biopython'``
    :param archive: Explicitely request reading an archive, type may be specified
       (default: auto-detected by file extension)

    All other kwargs are passed to the underlaying reader routine.

    The following formats are supported, for documentation of supported kwargs
    follow the provided links.

    {format_table}

    Example
    -------

    >>> from sugar import read
    >>> seqs = read('crazy_virus.fasta')  # read a local file  # doctest: +SKIP
    >>> seqs = read()  # load example file
    >>> print(seqs)  # doctest: +SKIP
    2 seqs in basket
    AB047639  9678  ACCTGCCCCTAATAGGGGCGACACTCCGCCATGAATCACTCCCCTGTGA...  GC:58.26%
    AB677533  9471  GCCCGCCCCCTGATGGGGGCGACACTCCGCCATGAATCACTCCCCTGTG...  GC:57.46%
      customize output with BioBasket.tostr() method

    >>> url = 'https://raw.githubusercontent.com/rnajena/sugar/master/sugar/tests/data/io_test.zip'
    >>> seqs = read(url)  # load an archive from the web  # doctest: +SKIP
    >>> print(seqs)  # doctest: +SKIP
    5 seqs in basket
    MCHU         150  MADQLTEEQIAEFKEAFSLFDKDGDGTITTKELGTVMRSLGQNPTEA...
    AAD44166.1   284  LCLYTHIGRNIYYGSYLYSETWNTGIMLLLITMATAFMGYVLPWGQM...
    BTBSCRYR     620  TGCACCAAACATGTCTAAAGCTGGAACCAAAATTACTTTCTTTGAAG...  GC:52.58%
    AB047639    9678  ACCTGCCCCTAATAGGGGCGACACTCCGCCATGAATCACTCCCCTGT...  GC:58.26%
    AB677533    9471  GCCCGCCCCCTGATGGGGGCGACACTCCGCCATGAATCACTCCCCTG...  GC:57.46%
      customize output with BioBasket.tostr() method
    """
    if fmt is None:
        fmt = detect(fname, **kw)
    if fmt is None:
        raise IOError('Format cannot be auto-detected')
    fmt = fmt.lower()
    if tool == 'biopython':
        from Bio import SeqIO
        seqs = SeqIO.parse(fname, fmt, **kw)
        return BioBasket.fromobj(seqs, 'biopython')
    elif tool:
        raise ValueError(f'Unrecognized tool: {tool}')
    module = EPS['seqs'][fmt].load()
    with _file_opener(fname, mode=mode, binary=_binary(module), encoding=encoding) as f:
        if hasattr(module, 'read'):
            seqs = module.read(f, **kw)
        elif hasattr(module, 'iter_'):
            seqs = list(module.iter_(f, **kw))
        else:
            raise RuntimeError(f'No read support for format {fmt}')
    seqs = BioBasket(seqs)
    for seq in seqs:
        seq.meta._fmt = fmt
    return seqs


@_resolve_fname(example_fname='!data/fts_example.gff')
def read_fts(fname, fmt=None, *, mode='r', encoding=None, **kw):
    """
    Read a file or file-like object with features into `.FeatureList`

    :param fname: Filename, can also be a glob expression,
        a web resource,
        an archive, gzipped file,
        or a file-like object (e.g. `~io.BytesIO`, `~io.StringIO`)
    :param fmt: format of the file (defaul: auto-detect)
    :param mode: mode for opening the file, change this only if you know what
        you do
    :param encoding: encoding of the file
    :param archive: Explicitely request reading an archive, type may be specified
       (default: auto-detected by file extension)

    All other kwargs are passed to the underlaying reader routine.

    The following formats are supported, for documentation of supported kwargs
    follow the provided links.

    {format_table}
    """
    if fmt is None:
        fmt = detect(fname, what='fts', **kw)
    if fmt is None:
        raise IOError('Format cannot be auto-detected')
    fmt = fmt.lower()
    module = EPS['fts'][fmt].load()
    with _file_opener(fname, mode=mode, binary=_binary(module, 'fts'), encoding=encoding) as f:
        if hasattr(module, 'read_fts'):
            fts = module.read_fts(f, **kw)
        else:
            raise RuntimeError(f'No fts read support for format {fmt}')
    for ft in fts:
        ft.meta._fmt = fmt
    return FeatureList(fts)


@_resolve_archive
def write(seqs, fname, fmt=None, *, mode='w', tool=None, encoding=None, **kw):
    """
    Write sequences to file, see `.BioBasket.write()`
    """
    if fmt is None:
        fmt = detect_ext(fname)
    if fmt is None:
        raise IOError('Format cannot be auto-detected')
    fmt = fmt.lower()
    if tool == 'biopython':
        from Bio import SeqIO
        return SeqIO.write(seqs.toobj('biopython'), fname, fmt)
    elif tool:
        raise ValueError(f'Unrecognized tool: {tool}')
    module = EPS['seqs'][fmt].load()
    with _file_opener(fname, mode=mode, binary=_binary(module), encoding=encoding) as f:
        if hasattr(module, 'append') and 'a' in mode:
            for seq in seqs:
                module.append(seq, f, **kw)
        elif hasattr(module, 'write'):
            module.write(seqs, f, **kw)
        elif hasattr(module, 'append') and 'w' in mode:
            for seq in seqs:
                module.append(seq, f, **kw)
        else:
            raise RuntimeError(f'No write support for format {fmt}')


@_resolve_archive
def write_fts(fts, fname, fmt=None, *, mode='w', **kw):
    """
    Write features to file, see `.FeatureList.write()`
    """
    if fmt is None:
        fmt = detect_ext(fname, 'fts')
    if fmt is None:
        raise IOError('Format cannot be auto-detected')
    fmt = fmt.lower()
    module = EPS['fts'][fmt].load()
    with _file_opener(fname, mode=mode, binary=_binary(module, 'fts')) as f:
        if hasattr(module, 'binary_fmt_fts') and module.binary_fmt_fts and 'b' not in mode:
            mode = 'b' + mode
        if hasattr(module, 'write_fts'):
            module.write_fts(fts, f, **kw)
        else:
            raise RuntimeError(f'No fts write support for format {fmt}')
