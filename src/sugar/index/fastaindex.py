# (C) 2024, Tom Eulenfeld, MIT license
"""
FASTA indexer
"""

# Comparison of binary db with esl-sfetch
# * allow fasta data to be located in multiple files
# * faster lookup for range queries in wrapped fasta files
# * ~10% higher file sizeq

import dbm
from glob import glob
import mmap
import sys
from warnings import warn
try:
    from tqdm import tqdm
except ImportError:
    tqdm = None
import os
try:
    from binarysearchfile import BinarySearchFile
except ImportError:
    BinarySearchFile = object
import sugar


SUGARFASTAINDEX_VERSION = '0.1.0'


class FastaBinarySearchFile(BinarySearchFile):
    """
    Custom `BinarySeachFile`_ used to store the index
    """
    magic = b'\xfe\x8a\x01\x01'
    headerstart = f'SugarFASTAindex v{SUGARFASTAINDEX_VERSION}, sugar v{sugar.__version__}\n'.encode('latin1')
    # seqid, file number, line length, offset in fasta file
    record = (0, 50, 50, 50)


def _iter_fasta_index(fname, seek=None, silent=False):
    """Read fasta file and yield (seqid, linelen, startbyte) triplets"""
    total = os.path.getsize(fname)
    if tqdm and not silent:
        desc = 'Index ' + os.path.split(fname)[1] + ' {:>11_d}'
        pbar = tqdm(desc=desc.format(0), total=total, unit='Byte', unit_scale=True)
        numseqs = 0
    else:
        pbar = None
    with open(fname, 'r+b') as file:
        with mmap.mmap(file.fileno(), 0) as f:
            if seek:
                f.seek(seek)
            if hasattr(f, 'madvise'):  # not available on Windows
                if seek:
                    ps = mmap.PAGESIZE
                    f.madvise(mmap.MADV_SEQUENTIAL, seek // ps * ps)
                else:
                    f.madvise(mmap.MADV_SEQUENTIAL)
            start = f.find(b'>')
            while start != -1:
                oldstart = start
                f.seek(start + 1)
                try:
                    seqid = f.readline().split()[0]
                except IndexError as ex:
                    raise ValueError(f'Pos {start}, empty id') from ex
                startline = f.tell()
                f.readline()
                endline = f.tell()
                start = f.find(b'>', endline)
                if endline == start or start == -1 and endline == total:
                    linelen = 0
                else:
                    linelen = endline - startline
                yield seqid, linelen, oldstart
                if pbar:
                    numseqs += 1
                    len_ = start - oldstart if start != -1 else total - oldstart
                    if pbar.update(len_):
                        pbar.set_description(desc.format(numseqs))
    if pbar:
        pbar.set_description(desc.format(numseqs))
        pbar.close()


def _extract_seqdata(fname, linelen, start, index=(None, None),
                     onlyheader=False):
    with open(fname, 'r+b') as file:
        with mmap.mmap(file.fileno(), 0) as f:
            f.seek(start)
            if onlyheader:
                return f.readline()
            elif index == (None, None):
                # full header + sequence
                k = f.find(b'>', start + 1)
                numbytes = None if k == -1 else k - start
                return f.read(numbytes)
            else:
                # we extract only part o the full sequence
                header = f.readline()
                offset = start + len(header)
                i, j = index
                if linelen != 0:
                    # fasta data might stretch over more than one line
                    if header[-2] in (b'\n', b'\r'):
                        nlec = 2  # Windows file ending, two chars
                    else:
                        nlec = 1  # Linux file ending, one char
                if i is None:
                    i = 0
                elif i < 0:
                    warn('Start index is not allowed to be smaller than 0, set to 0')
                    i = 0
                elif linelen != 0:
                    i += (i // (linelen - nlec)) * nlec
                if j is None:
                    # need to find end of record
                    k = f.find(b'>', offset + i)
                    j = None if k == -1 else k - offset  # read to end of file, if no other record afterwards
                elif j < 0:
                    raise ValueError('End index is not allowed to be smaller than 0')
                else:
                    if linelen != 0:
                        j += (j // (linelen - nlec)) * nlec
                    if (k := f.find(b'>', offset + i, offset + j)) != -1:
                        warn('End index is too large, set to sequence end')
                        j = k - offset
                f.seek(i, 1)
                numbytes = None if j is None else j - i
                return header + f.read(numbytes)


def _humanb(size,p=2):
    suf = ('','k', 'M', 'G', 'T')
    i = 0
    while size > 1024:
        i += 1
        size = size / 1024
    return f'{size:.{p}f} {suf[i]}Byte'


def _int(b):
    return int.from_bytes(b)

def _pack(fn, linelen, start):
    return fn.to_bytes(2) + linelen.to_bytes(2) + start.to_bytes(start.bit_length())

def _unpack(b):
    return _int(b[:2]), _int(b[2:4]), _int(b[4:])


def _cache_dbname(dbname):
    try:
        from platformdirs import user_cache_dir
    except ImportError:
        pass
    else:
        cachedir = user_cache_dir('sugar', 'rnajena')
        fname = os.path.join(cachedir, 'last_used_fastaindexfile.txt')
        if dbname is None and (os.path.exists(fname) or dbm.whichdb(fname)):
                with open(fname) as f:
                    dbname = f.read().strip()
        elif dbname is not None:
            os.makedirs(cachedir, exist_ok=True)
            with open(fname, 'w') as f:
                f.write(os.path.abspath(dbname))
    if dbname is None:
        raise ValueError('Dbname not specified')
    return os.path.abspath(dbname)


class FastaIndex():

    """
    Index FASTA files and query (sub-) sequences

    :param str dbname: Name of index file, default: last used index file
    :param str create: Create a new index file
    :param str path: Common path for FASTA file names, default ``'{dbpath}'``, i.e.
        filenames are relative to path of index file, only needed for new index file
    :param str mode: ``'binary'`` uses a binary search file,
        ``'db'`` uses a database via Python's `dbm` module,
        only needed for new index file
    """

    def __init__(self, dbname=None, create=False, path='{dbpath}', mode=None):
        self.dbname = _cache_dbname(dbname)
        exists = self._exists()
        if mode not in (None, 'binary', 'db'):
            raise ValueError(f"Mode '{mode}' is not a valid mode")
        if create or not exists:  # new index file
            if create and exists:
                raise ValueError(f'File {dbname} exists. Delete the file before creating a new database.')
            elif not create and not exists:
                warn(f'File {dbname} does not exist. Creating a new file.')
            if mode is None:
                raise ValueError('Please specify mode, binary (binary search index) or db (database)')
            if path != '' and path[-1] != '/':
                path = path + '/'
            self.mode = mode
            self.path = path
            self.files = []
            if self.mode == 'db':
                self.db = dbm.open(self.dbname, flag='c')
            else:
                self.db = FastaBinarySearchFile(self.dbname)
        else:  # not create and exists, load existing index file
            if mode is None:
                if dbm.whichdb(self.dbname) not in ('', None):
                    mode = 'db'
                elif FastaBinarySearchFile.check_magic(self.dbname, n=2):
                    mode = 'binary'
                else:
                    raise ValueError(f'File {dbname} is unknown index file')
            self.mode = mode
            self._read_header()
        dbpath = os.path.dirname(self.dbname)
        self.pathf = self.path.format(dbpath=dbpath)

    def __enter__(self):
        return self

    def __exit__(self, type, value, traceback):
        if self.mode == 'db':
            self.db.close()

    def _read_header(self):
        if self.mode == 'db':
            self.db = dbm.open(self.dbname)
            self.path, *self.files = self.db['header'].decode('latin1').split(',')
        else:
            self.db = FastaBinarySearchFile(self.dbname)
            self.path, *self.files = map(str.strip, self.db.read_header().decode('latin1').split('\n')[1].split(','))

    def add(self, fnameexpr, seek=None, force=False, silent=False):
        if self.mode == 'binary' and self._exists() and len(self)>0 and not force:
            raise ValueError('For best performance of binary search index creation '
                             'use the add command only once. '
                             'Force update of the file with the force flag/keyword.')
        if not isinstance(fnameexpr, str):
            fnames = sorted(fname for expr in fnameexpr for fname in glob(expr))
        else:
            fnames = sorted(glob(fnameexpr))
        adddata = []
        for fname in fnames:
            fnamerel = os.path.relpath(fname, self.pathf)
            try:
                fn = int(self.files.index(fnamerel))
            except ValueError:
                fn = len(self.files)
                self.files.append(fnamerel)
            for seqid, *data in _iter_fasta_index(fname, seek=seek, silent=silent):
                if self.mode == 'db':
                    self.db[seqid] = _pack(fn, *data)
                else:
                    adddata.append((seqid, fn) + tuple(data))
        # add some whitespace at path variable
        # in this way it can be easily changed without having to rewrite the
        # whole binary search file
        header = (','.join([self.path + (' ' * 50) * (self.mode == 'binary')] +
                           self.files)).encode('latin1')
        if self.mode == 'db':
            self.db['header'] = header
        elif force:
            self.db.update(adddata, header=header)
        else:
            self.db.write(adddata, header=header)

    def _search(self, seqids, **kw):
        if (isinstance(seqids, (bytes, str)) or
                len(seqids) == 3 and not isinstance(seqids[1], (bytes, str))):
            seqids = [seqids]
        for seqid in seqids:
            if not isinstance(seqid, (bytes, str)):
                seqid, start, stop = seqid
                index = (start, stop)
            else:
                index = (None, None)
            if self.mode == 'db':
                fn, *args = _unpack(self.db[seqid])
            else:
                if isinstance(seqid, str):
                    seqid = seqid.encode('latin1')
                seqid_, fn, *args = self.db[seqid]
                assert seqid_ == seqid
            yield _extract_seqdata(self.pathf + self.files[fn], *args, index=index, **kw).decode('latin1')

    def iter(self, seqids):
        """
        Yield `.BioSeq` sequences from given sequence ids

        :param seqids: List of sequence ids, might also be a single seqid,
           a 3-len tuple with start and stop indices ``(seqid, start, stop)``
           (start or stop can be None), or a list with 3-len tuples
        """
        from sugar import BioBasket
        for seqdata in self._search(seqids):
            yield BioBasket.fromfmtstr(seqdata, fmt='fasta')[0]

    def iter_fasta(self, seqids):
        """
        Yield FASTA strings from given sequence ids

        See `iter()` for documentation of seqids argument.
        """
        for seqdata in self._search(seqids):
            yield seqdata

    def iter_fastaheader(self, seqids):
        """
        Yield FASTA headers from given sequence ids

        See `iter()` for documentation of seqids argument.
        """
        for seqdata in self._search(seqids, onlyheader=True):
            yield seqdata

    def get_seq(self, seqid):
        """
        Return a `.BioSeq` given its id

        See `iter()` for documentation of seqid argument.
        """
        return list(self.iter([seqid]))[0]

    def get_basket(self, seqids):
        """
        Return `.BioBasket` with all corresponding sequences

        See `iter()` for documentation of seqids argument.
        """
        from sugar import BioBasket
        return BioBasket(list(self.iter(seqids)))

    def get_fasta(self, seqids):
        """
        Return FASTA string with all sequences

        See `iter()` for documentation of seqid argument.
        """
        return ''.join(self.iter_fasta(seqids))

    def get_fastaheader(self, seqids):
        """
        Return FASTA header from all sequences

        See `iter()` for documentation of seqid argument.
        """
        return ''.join(self.iter_fastaheader(seqids))

    def _exists(self):
        # work-around for dbm.dumb, which creates 3 files
        return os.path.exists(self.dbname) or dbm.whichdb(self.dbname)


    def __len__(self):
        if self.mode == 'db':
            return len(self.db) - 1
        else:
            return len(self.db)

    def __str__(self):
        if not self._exists(): #self.files:
            return f'File {self.dbname} does not exist yet'
        if self.mode == 'binary':
            try:
                sizeinfo = f'   recsize: {self.db.attr.reclen} Byte  {tuple(self.db.size)}\n'
            except TypeError:
                sizeinfo = '    recsize: no records\n'
        else:
            sizeinfo = ''
        return ('Sugar index\n'
                f'    dbname: {self.dbname}\n'
                f'      type: {self.dbtype}\n'
                f'  num seqs: {len(self):_d}\n'
                f'    dbsize: {_humanb(self.totalsize)}\n' + sizeinfo +
                f'      path: {self.path}\n'
                '     files:\n         ' +
                '        ' + ' '.join(self.files) + '\n')

        return (f'{type(self)}\n'
                f'     fname: {self.attr.fname}\n'
                f'   records: {len(self):_d}\n'
                f'      size: {_humanb(self.attr.totalsize)}\n'
                f'   recsize: {self.attr.recsize} Byte  {tuple(self.size)}\n')

    @property
    def dbtype(self):
        return type(self.db).__name__

    @property
    def totalsize(self):
        try:
            return os.path.getsize(self.dbname)
        except FileNotFoundError:
            if self._exists():
                # work-around for dbm.dumb, which creates 3 files
                return os.path.getsize(self.dbname+'.dat')
            raise


def _parse_cmd_seqids(seqids):
    if ',' in seqids:  # seqid, start, stop
        seqid, start, stop = map(str.strip, seqids.split(','))
        return seqid, int(start), None if stop in '0-' else int(stop)
    if not os.path.exists(seqids) and seqids != '-':
        return seqids
    if seqids == '-':
        data = sys.stdin.read()
    else:
        with open(seqids) as f:
            data = f.read()
    query = []
    for line in data.splitlines():
        id_, i, j = map(str.strip, line.split(',' if ',' in data else None))
        query.append([id_, int(i), None if j in '0-' else int(j)])
    return query


def _fastaindex_cmd(idxcommand, dbname, mode=None, path=None, seqids=None,
          fnames=None, out='-', force=False):
    match idxcommand:
        case 'info':
            print(FastaIndex(dbname))
        case 'create':
            index = FastaIndex(dbname, create=True, mode=mode, path=path)
            index.add([])
        case 'add':
            index = FastaIndex(dbname)
            index.add(fnames, force=force)
        case 'print':
            index = FastaIndex(dbname)
            print(index.get_basket(_parse_cmd_seqids(seqids)))
        case 'load':
            from sugar.scripts import _start_ipy
            index = FastaIndex(dbname)
            seqs = index.get_basket(_parse_cmd_seqids(seqids))
            _start_ipy(seqs)
        case 'fetch':
            index = FastaIndex(dbname)
            print(index.get_fasta(_parse_cmd_seqids(seqids)))
