# (C) 2024, Tom Eulenfeld, MIT license
"""
Sequence related classes, `.BioSeq`, `.BioBasket`
"""

import collections
import collections.abc
import copy
from functools import reduce
import io
import operator
import sys
from warnings import warn

from sugar.data import CODES
from sugar.core.fts import Feature, FeatureList, Location
from sugar.core.meta import Attr, Meta


CODES_INV = {frozenset(v): k for k, v in CODES.items()}
COMPLEMENT = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G', '.': '.', '-': '-'}
COMPLEMENT_ALL = {c: CODES_INV[frozenset(COMPLEMENT[nt] for nt in nts)] for c, nts in CODES.items()}
COMPLEMENT_TRANS = str.maketrans(COMPLEMENT_ALL)


class _Slicable_GetItem():
    def __init__(self, obj, **kw):
        self.obj = obj
        self.kw = kw

    def __getitem__(self, i):
        return self.obj.getitem(i, **self.kw)


class _BioSeqStr():
    """
    Helper class to hold all string methods in the `BioSeq.str` namespace

    :meta public:
    """
    def __init__(self, parent):
        self.__parent = parent

    def __deepcopy__(self, orig):
        # TODO test
        return self


    def center(self, width, *args):
        self.__parent.data = self.__parent.data.center(width, *args)
        return self.__parent

    def count(self, sub, start=0, end=sys.maxsize):
        return self.__parent.data.count(str(sub), start, end)

    def removeprefix(self, prefix, /):
        self.__parent.data = self.__parent.data.removeprefix(str(prefix))
        return self.__parent

    def removesuffix(self, suffix, /):
        self.__parent.data = self.__parent.data.removesuffix(str(suffix))
        return self.__parent

    def encode(self, encoding='utf-8', errors='strict'):
        encoding = 'utf-8' if encoding is None else encoding
        errors = 'strict' if errors is None else errors
        return self.__parent.data.encode(encoding, errors)

    def endswith(self, suffix, start=0, end=sys.maxsize):
        return self.__parent.data.endswith(str(suffix), start, end)

    # def expandtabs(self, tabsize=8):
    #     return self.__class__(self.data.expandtabs(tabsize))

    def find(self, sub, start=0, end=sys.maxsize):
        return self.__parent.data.find(str(sub), start, end)

    def format(self, /, *args, **kwds):
        return self.__parent.data.format(*args, **kwds)

    def format_map(self, mapping):
        return self.__parent.data.format_map(mapping)

    def index(self, isub, start=0, end=sys.maxsize):
        return self.__parent.data.index(str(sub), start, end)

    def isalpha(self):
        return self.__parent.data.isalpha()

    # def isalnum(self):
    #     return self.__parent.data.isalnum()

    def isascii(self):
        return self.__parent.data.isascii()

    # def isdecimal(self):
    #     return self.data.isdecimal()

    # def isdigit(self):
    #     return self.data.isdigit()

    # def isidentifier(self):
    #     return self.data.isidentifier()

    def islower(self):
        return self.__parent.data.islower()

    # def isnumeric(self):
    #     return self.data.isnumeric()

    # def isprintable(self):
    #     return self.data.isprintable()

    # def isspace(self):
    #     return self.data.isspace()

    # def istitle(self):
    #     return self.data.istitle()

    def isupper(self):
        return self.__parent.data.isupper()

    # def join(self, seq):
    #     return self.data.join(seq)

    def ljust(self, width, *args):
        self.__parent.data = self.__parent.data.ljust(width, *args)
        return self.__parent

    def lower(self):
        self.__parent.data = self.__parent.data.lower()
        return self.__parent

    def lstrip(self, chars=None):
        self.__parentdata = self.__parent.data.lstrip(chars)
        return self.__parent

    @staticmethod
    def maketrans(*args):
        return str.maketrans(*args)

    # def partition(self, sep):
    #     return self.data.partition(sep)

    def replace(self, old, new, maxsplit=-1):
        self.__parent.data = self.__parent.data.replace(
            str(old), str(new), maxsplit)
        return self.__parent

    def rfind(self, sub, start=0, end=sys.maxsize):
        return self.__parent.data.rfind(str(sub), start, end)

    def rindex(self, sub, start=0, end=sys.maxsize):
        return self.__parent.data.rindex(str(sub), start, end)

    def rjust(self, width, *args):
        self.__parent.data = self.__parent.data.rjust(width, *args)
        return self.__parent

    # def rpartition(self, sep):
    #     return self.data.rpartition(sep)

    def rstrip(self, chars=None):
        if chars is not None:
            chars = str(chars)
        self.__parent.data = self.__parent.data.rstrip(chars)
        return self.__parent

    def split(self, sep=None, maxsplit=-1):
        if sep is not None:
            sep = str(sep)
        return self.__parent.data.split(sep, maxsplit)

    def rsplit(self, sep=None, maxsplit=-1):
        if sep is not None:
            sep = str(sep)
        return self.__parent.data.rsplit(sep, maxsplit)

    def splitlines(self, keepends=False):
        return self.__parent.data.splitlines(keepends)

    def startswith(self, prefix, start=0, end=sys.maxsize):
        return self.__parent.data.startswith(str(prefix), start, end)

    def strip(self, chars=None):
        if chars is not None:
            chars = str(chars)
        self.__parent.data = self.__parent.data.strip(chars)
        return self.__parent

    def swapcase(self):
        self.__parent.data = self.__parent.data.swapcase()
        return self.__parent

    # def title(self):
    #     return self.__class__(self.data.title())

    def translate(self, *args):
        self.__parent.data = self.__parent.data.translate(*args)
        return self.__parent

    def upper(self):
        self.__parent.data = self.__parent.data.upper()
        return self.__parent

    # def zfill(self, width):
    #     return self.__class__(self.data.zfill(width))

class _BioBasketStr():
    """
    Helper class to move all string methods into the `BioBasket.str` namespace

    It calls the corresponding `BioSeq.str` method under the hood and returns
    either the altered `BioBasket` object or a list of results.

    :meta public:
    """
    def __init__(self, parent):
        self.__parent = parent

    def __getattr__(self, name):
        def method(*args, **kw):
            results = [
                getattr(seq.str, name)(*args, **kw)
                for seq in self.__parent
            ]
            if name in (
                'center', 'remove_prefix', 'ljust', 'lower',
                'lstrip', 'replace', 'rjust', 'rstrip',
                'strip', 'swapcase', 'translate', 'upper'
                ):
                return self.__parent
            else:
                return results
        return method


def _detect_tool(obj):
    try:
        from Bio.Seq import Seq
        from Bio.SeqRecord import SeqRecord
    except ImportError:
        pass
    else:
        if isinstance(obj, (Seq, SeqRecord)):
            return 'biopython'


class BioSeq():
    """
    Class holding sequence data and metadata, exposing bioinformatics methods.

    Most methods operate by default in-place, but return the BioSeq object again.
    Therefore, method chaining can be used.
    """

    def __init__(self, data, id='', meta=None, type=None):
        #: Namespace holding all available string methods,
        #: see `_BioSeqStr` for available methods
        #: and `python:str` for documentation of the methods
        #:
        #: .. rubric:: Example:
        #:
        #: >>> seq = read()[0]
        #: >>> seq.str.find('ATG')  # Use string method
        #: 30
        self.str = _BioSeqStr(self)
        #: Property holding the data string
        self.data = str(data).upper()
        if hasattr(data, 'meta'):
            meta = data.meta
        elif 'meta' in data:
            meta = data['meta']
        elif meta is None:
            meta = {}
        #: Property holding metadata
        self.meta = Meta(meta)
        if id or 'id' not in self.meta:
            self.meta.id = id
        assert type in (None, 'nt', 'aa')
        if type is None:
            codes = set(CODES) | {'U'}
            type = 'nt' if all(nb in codes for nb in self.data) else 'aa'
        if type == 'nt' and 'U' in data:
            from warnings import warn
            warn('Found U in nucleotide sequence, '
                 'some methods will not work as expected')
        #: type of the sequence, either ``'nt'`` or ``'aa'``
        self.type = type

    def _repr_pretty_(self, p, cycle):
        if cycle:
            p.text('...')
        else:
            p.text(str(self))

    def __str__(self):
        return str(self.data)

    def __repr__(self):
        metastr = ', '.join(f'{prop}={repr(val)}' for prop, val in vars(self.meta).items())
        return f'{type(self).__name__}([{repr(self.data)}, meta=dict({metastr}))'

    def __eq__(self, string):
        if isinstance(string, BioSeq):
            return self.data == string.data and self.meta == string.meta
        return self.data == string

    def __lt__(self, string):
        if isinstance(string, BioSeq):
            return self.id < string.id
        self.id < ''

    def __le__(self, string):
        if isinstance(string, BioSeq):
            return self.id <= string.id
        self.id <= ''

    def __gt__(self, string):
        if isinstance(string, BioSeq):
            return self.id > string.id
        self.id > ''

    def __ge__(self, string):
        if isinstance(string, BioSeq):
            return self.id >= string.id
        self.id >= ''

    def __contains__(self, char):
        return str(char) in self.data

    def __len__(self):
        return len(self.data)

    def __setitem__(self, index, value):
        l = list(self.data)
        l[index] = value
        self.data = ''.join(l)

    def __add__(self, other):
        if isinstance(other, BioSeq) and self.meta != other.meta:
            warn('Try to add two BioSeq objects with different meta data')
        return self.__class__(self.data + str(other), meta=self.meta)

    def __iadd__(self, other):
        self = self + other
        return self

    def __radd__(self, other):
        return self.__class__(str(other) + self.data, meta=self.meta)

    @property
    def id(self):
        """Alias for ``BioSeq.meta.id``"""
        return self.meta.id

    @id.setter
    def id(self, value):
        self.meta.id = value

    @property
    def fts(self):
        """
        Alias for ``BioSeq.meta.fts``

        The fts object holds all feature metadata.
        It is an instance of `.FeatureList`.
        """
        return self.meta.setdefault('fts', FeatureList())

    @fts.setter
    def fts(self, value):
        self.meta.fts = FeatureList(value)
        for ft in self.meta.fts:
            if ft.seqid != self.id:
                warn('Feature seqid and sequence id mismatch')

    def add_fts(self, fts):
        """
        Add some features to the feature list.

        If you want set all features use the `BioSeq.fts` attribute.

        :param fts: features to add
        """
        self.fts = self.fts + fts
        self.fts.sort()

    @property
    def gc(self):
        """
        GC content of the sequence
        """
        GC = self.str.count('G') + self.str.count('C')
        AT = self.str.count('A') + self.str.count('T') + self.str.count('U')
        if GC + AT > 0:
            return GC / (GC + AT)
        else:
            return 0

    @property
    def i(self):
        """Return slicable object to support in-place slicing

        Deprecated: Use getitem() or sl attribute.
        """
        msg = 'BioSeq.i is deprecated, use geitem() method or sl attribute'
        warnings.warn(msg, DeprecationWarning, stacklevel=2)
        return _Slicable_GetItem(self, inplace=True)

    def rc(self):
        """
        Reverse complement, alias for ``BioSeq.reverse().complement()``
        """
        return self.reverse().complement()


    def __getitem__(self, index):
        return self.getitem(index)

    def getitem(self, index, inplace=False, gap=None):
        """
        Slice the sequence and return a subsequence

        This is the method which is called if you slice with ``BioSeq[]`` syntax.
        If you want to use non-default options call this method directly,
        or by the `BioSeq.sl` attribute.

        .. rubric:: Example:

        >>> from sugar import read
        >>> seq = read()[0]
        >>> print(seq[5:10])
        CCCCT
        >>> print(seq[5])
        C
        >>> print(seq['cds'][:3])
        ATG

        :param index: Specifies which part of the sequence is returned.
           The following types are supported.

           int,slice
             location is specified by int or slice
           `.Location`
             specified by location
           `.Feature`
             specified by feature
           str
             position of first feature of given type, e.g. ``'cds'``
             will return sequence with first coding sequence
        :param bool inplace:
            The subsequence is not only returned, but the original
            sequence is modified in-place (default: False)
        :param str gap:
            gaps of the given characters will be accounted for when
            slicing the sequence (default: gaps will not be accounted for)
        """
        # TODO: add correct_fts kwargs
        try:
            if gap is not None:
                # from bisect import bisect
                nogaps = [i for i, nt in enumerate(self.data) if nt not in gap]
                adj = lambda i: nogaps[i] if i is not None and i < len(nogaps) else None
                if isinstance(index, int):
                    index = adj(index)
                elif isinstance(index, slice):
                    index = slice(adj(index.start), adj(index.stop), index.step)
            subseq = self.__class__(self.data[index], meta=self.meta)
            # subseq = super().getitem(index, gap=gap)
        except:
            if isinstance(index, str):
                index = self.fts.get(index)
                if index is None:
                    raise TypeError('Feature not found')
            if isinstance(index, Location):
                from sugar.core.fts import _slice_locs
                subseq = _slice_locs(self, [index], gap=gap)
            elif isinstance(index, Feature):
                from sugar.core.fts import _slice_locs
                subseq = _slice_locs(self, index.locs, gap=gap)
                # index = index._slice()
                # if index is None:
                #     msg = f'Feature {index.type} of seq {self.id} has no location'
                #     raise TypeError(msg)
                # subseq = super().__getitem__(index)
            else:
                raise TypeError('Index not supported')
        if inplace:
            self.data = subseq.data
        return subseq
        # it follows a lot of code to keep the feature indices intact
        # this is not really necessary I guess, but I try anyway
        # if 'features' in subseq.meta:
        #     subseq.meta.features = FeatureList([])
        #     if isinstance(index, int):
        #         index = slice(index, index+1, 1)
        #     (_, _, istep) = index.indices(len(self))
        #     if istep in (1, -1):
        #         if istep == 1:
        #             for ft in self.meta.features:
        #                 if (getattr(ft, 'start') is None or
        #                         istep == 1 and getattr(ft, 'stop', -1) == -1):
        #                     continue
        #                 stride = getattr(ft, 'stride', 1)
        #                 if stride == 1:
        #                     assert ft.start < ft.stop
        #                     start = max(0, ft.start - index.start)  # index.start might be None
        #                     stop = ft.stop - index.start if index.stop is None else min(ft.stop - index.start, index.stop)
        #                     if start >= stop:
        #                         continue
        #                 else:## todo wahaha
        #                     assert ft.start > ft.stop or ft.stop is None
        #                     start = max(ft.start - index.stop, -1)
        #                     stop = min(ft.stop, index.start-1)
        #                     if stop < 0:
        #                         stop = None
        #                     if start <= stop:
        #                         continue
        #                 ft2 = ft.copy()
        #                 ft2.start = start
        #                 ft2.stop = stop
        #                 ft2.stride = stride
        #                 if abs(stop - start) < abs(ft.start - ft.stop):
        #                     ft2.orig_len = getattr(ft, 'orig_len', abs(ft.start - ft.stop))
        #                 subseq.meta.features.append(ft2)
        #         else:
        #             for ft in self.meta.features[::-1]:
        #                 if (getattr(ft, 'start') is None or
        #                         istep == 1 and getattr(ft, 'stop', -1) == -1):
        #                     continue
        #                 stride = getattr(ft, 'stride', 1)
        #                 if stride == 1:
        #                     assert ft.start < ft.stop
        #                     start = ft.stop-1 if index.start is None else max(ft.stop-1, index.start)
        #                     stop = min(ft.start-1, index.stop)  # index.stop might be None
        #                     if stop < 0:
        #                         stop = None
        #                     if start <= stop:
        #                         continue
        #                 else:
        #                     assert ft.start > ft.stop or ft.stop is None
        #                     start = 0 if ft.stop is None and index.stop is None else max(ft.stop+1, index.stop+1)
        #                     stop = min(ft.start-1, index.start-1)
        #                     if start >= stop:
        #                         continue
        #                 print(ft)
        #                 print(index, istep)
        #                 print(start, stop)
        #                 ft2 = ft.copy()
        #                 ft2.start = start
        #                 ft2.stop = stop
        #                 ft2.stride = -stride
        #                 if abs(stop-start) < abs(ft.start - ft.stop):
        #                     ft2.orig_len = getattr(ft, 'orig_len', abs(ft.start - ft.stop))
        #                 subseq.meta.features.append(ft2)
        # return subseq


    def sl(self, **kw):
        """
        Method allowing to call `BioSeq.getitem()` with non-default options and extended indexing syntax

        Returns a slicable object. Use the ``BioSeq[]`` notation directly if you use default arguments.

        .. rubric:: Example:

        >>> from sugar import read
        >>> seq = read()[0]
        >>> print(seq[:5])
        ACCTG
        >>> print(seq.sl(inplace=True, gap='-')[:5:2])
        ACG
        >>> print(seq)  # was modified in-place
        ACG
        """
        return _Slicable_GetItem(self, **kw)


    def biotranslate(self, *args, **kw):
        from sugar.core.translate import translate
        import warnings
        msg = 'BioSeq.biotranslate() is deprecated, use translate() method'
        warnings.warn(msg, DeprecationWarning, stacklevel=2)
        self.data = translate(self.data, *args, **kw)
        self.type = 'aa'
        return self

    def complement(self):
        """
        Complementary sequence, i.e. transcription
        """
        if 'U' in self.data:
            self.str.replace('U', 'T').str.translate(COMPLEMENT_TRANS).str.replace('T', 'U')
        else:
            self.str.translate(COMPLEMENT_TRANS)
        return self

    def countall(self, **kw):
        return BioBasket([self]).countall(**kw)

    def countplot(self, hue=None, **kw):
        return BioBasket([self]).countplot(hue=hue, **kw)

    def copy(self):
        """
        Return a deep copy of the object
        """
        return copy.deepcopy(self)

    @classmethod
    def fromobj(cls, obj, tool=None):
        """
        Create a `BioSeq` object from a Python object created by another tool.

        This method is *WIP* - work in progress.

        :param obj: The object to convert.
        :param tool: the used tool (default: autodetect,
            allowed: ``'biopython'``)
        """
        if tool is None:
            tool = _detect_tool(obj)
        if tool is None:
            raise ValueError('Cannot determine origin of object')
        if tool == 'biopython':
            if hasattr(obj, 'seq'):  # SeqRecord
                seq = str(obj.seq)
                id_ = obj.id
            else:  # Seq
                seq = str(obj)
                id_ = None
            return cls(seq, id=id_)
        else:
            raise ValueError(f'Unsupported tool: {tool}')

    def match(self, *args, **kw):
        """
        Search regex and return match, see ``~.cane.match()``
        """
        from sugar.core.cane import match as _match
        return _match(self, *args, **kw)

    def matchall(self, *args, **kw):
        """
        Search regex and return `.BioMatchList` with all matches, see `~.cane.match()`
        """
        kw['matchall'] = True
        return self.match(*args, **kw)

    def find_orfs(self, *args, **kw):
        """
        Find ORFS in the sequence, see `~.cane.find_orfs()`
        """
        from sugar.core.cane import find_orfs
        return find_orfs(self, *args, **kw)

    def tofmtstr(self, fmt, **kw):
        """
        Write object to a string of specified format, see `~.main.write()`
        """
        return BioBasket([self]).tofmtstr(fmt, **kw)

    def toobj(self, tool=None):
        """
        Convert the object to an object of another bioinformatics Python tool.

        *WIP* - work in progress.

        :param tool: Only ``'biopython'`` is supported.
        """
        if tool == 'biopython':
            from Bio.Seq import Seq
            from Bio.SeqRecord import SeqRecord
            return SeqRecord(Seq(self.data), id=self.id)
        else:
            raise ValueError(f'Unsupported tool: {tool}')

    def reverse(self):
        """
        Reverse the sequence
        """
        self.data = self.data[::-1]
        return self

    def translate(self, *args, **kw):
        """
        Translate nucleotide sequence to amino acid sequence, see `~.cane.translate()`.

        The original translate method of the str class can be used via ``BioBasket.str.translate()``.
        """
        from sugar.core.cane import translate
        self.data = translate(self.data, *args, **kw)
        self.type = 'aa'
        return self

    def write(self, fname, fmt=None, **kw):
        """
        Write sequence to file, see `~.main.write()`
        """
        BioBasket([self]).write(fname, fmt, **kw)


import math


def _si_format(v, l=4):
    if v == 0:
        return '0'
    l2 = int(math.log10(v)) + 1
    if l2 > l:
        si = {1: 'k', 2: 'M', 3: 'G', 4: 'T'}
        j = (3 + l2 - l) // 3
        s = si.get(j, '?')
        v = v / 10 ** (3*j)
        # p = min(max(0, l - (l2 - 3*j) - 1), 1)
        p = 1 if l2 % 3 == l % 3 else 0
        return f'{v:.{p}f}{s}'
    else:
        return str(v)


class BioBasket(collections.UserList):
    """
    Class holding a list of `BioSeq` objects

    The BioBasket object can be used like a list.
    It has useful bioinformatics methods attached.

    The list itself is stored in the ``data`` property.
    The BioBasket object may also have an metadata
    attribute.
    """
    def __init__(self, data=None, meta=None):
        # Documentation for str attribute:
        #: Namespace holding all available string methods,
        #:
        #: The `BioBasket.str` methods call the corresponding `BioSeq.str` methods under the hood
        #: and return either the altered `BioBasket` object or a list with results.
        #: See `_BioSeqStr` for available methods
        #: and `python:str` for documentation of the methods
        #:
        #: .. rubric:: Example:
        #:
        #: >>> seqs = read()
        #: >>> seqs.str.find('ATG')  # Use string method
        #: [30, 12]
        self.str = _BioBasketStr(self)
        if data is None:
            data = []
        if hasattr(data, 'meta'):
            meta = data.meta
        elif 'meta' in data:
            meta = data['meta']
        elif meta is None:
            meta = {}
        super().__init__(data)
        # FAKE, just for documenting tha data property
        #: Property holding the list of sequences
        self.data = self.data
        #: Property holding metadata
        self.meta = Meta(meta)

    def __eq__(self, other):
        if isinstance(other, BioBasket):
            return self.data == other.data and self.meta == other.meta
        return self.data == other

    @property
    def ids(self):
        """List of sequence ids"""
        return [seq.meta.id for seq in self]

    @property
    def fts(self):
        """
        `.FeatureList` of containing features of all sequences

        Can also be used as setter.
        Code example: ``seqs.fts = new_fts``.
        """
        fts = [ft for seq in self for ft in seq.fts]
        return FeatureList(fts)

    @fts.setter
    def fts(self, value):
        fts = FeatureList(value).todict()
        for seq in self:
            if seq.id in fts:
                seq.fts = fts.pop(seq.id)
        if len(fts) > 0:
            missing_ids = ', '.join(fts.keys())
            warn(f'Features for seqids {missing_ids} could not be '
                 'attached to any sequence')

    def add_fts(self, fts):
        """
        Add some features to the feature list of corresponding sequences.

        If you want set all features use the `BioBasket.fts` attribute.

        :param fts: features to add
        """
        fts = FeatureList(fts).todict()
        for seq in self:
            if seq.id in fts:
                seq.fts = seq.fts + fts.pop(seq.id)
                seq.meta.fts.sort()
        if len(fts) > 0:
            missing_ids = ', '.join(fts.keys())
            warn(f'Features for seqids {missing_ids} could not be '
                 'attached to any sequence')

    def rc(self):
        """
        Reverse complement, alias for ``BioSeq.reverse().complement()``
        """
        return self.reverse().complement()

    def __str__(self):
        return self.tostr(add_hint=True)

    def _repr_pretty_(self, p, cycle):
        if cycle:
            p.text('...')
        else:
            p.text(str(self))

    def __repr__(self):
        metastr = ', '.join(f'{prop}={repr(val)}' for prop, val in vars(self.meta).items())
        return f'{type(self).__name__}({super().__repr__()}, meta=dict({metastr}))'

    def __getitem__(self, i):
        return self.getitem(i)

    def getitem(self, i, **kw):
        """
        Slice sequences

        This is the method which is called if you slice with ``BioBasket[]`` syntax.
        If you want to use non-default options call this method directly,
        or by the `BioBasket.sl` attribute.

        .. rubric:: Example:

        >>> from sugar import read
        >>> seqs = read()
        >>> print(seqs[:2, 5:10])
        2 seqs in basket
        AB047639  5  CCCCT  ...
        AB677533  5  CCCCC  ...
        >>> print(seqs[:2, 'cds'][:, 0:3])
        2 seqs in basket
        AB047639  3  ATG  ...
        AB677533  3  ATG  ...

        :param index:
            Specifies which part of the sequences or which sequences are returned.

            int
                Returns a `BioSeq` from the basket
            slice
                Returns a new `BioBasket` object with a subset of the sequences
            str,feature,location
                Updates all sequences inisde the basket, see `BioSeq.getitem()`
            (int, object)
                Returns a `BioSeq` from the basket and slices it with the object, see `BioSeq.getitem()`
            (slice, object)
                Returns a new `BioBasket` object with a subset of the sequences which are replaced
                by subsequences according to `BioSeq.getitem()`
        :param \*\*kw:
            Aditional kwargs are passed to `BioSeq.getitem()`.
        """
        if isinstance(i, int):
            return self.data[i]
        elif isinstance(i, slice):
            seqs = self.__class__(self.data[i], meta=self.meta)
        elif isinstance(i, (str, Feature, Location)):
            seqs = self.__class__(self.data, meta=self.meta)
            seqs.data = [seq.getitem(i, **kw) for seq in seqs.data]
        elif len(i) == 2:
            i, j = i
            if isinstance(i, int):
                return self.data[i].getitem(j, **kw)
            elif isinstance(i, slice):
                seqs = self.__class__(self.data[i], meta=self.meta)
                seqs.data = [seq.getitem(j, **kw) for seq in seqs.data]
            else:
                raise TypeError('Index not supported')
        else:
            raise TypeError('Index not supported')
        return seqs

    def sl(self, **kw):
        """
        Method allowing to call `BioBasket.getitem()` with non-default options and extended indexing syntax

        Returns a slicable object. Use the ``BioBasket[]`` notation directly if you use default arguments.
        """
        return _Slicable_GetItem(self, **kw)

    def __setitem__(self, i, value):
        if isinstance(i, (int, slice)):
            if isinstance(value, (str, BioSeq)):
                value = BioSeq(value)
            else:  # assume its list or BioBasket or similar
                value = [BioSeq(v) for v in value]
            self.data[i] = value
        elif len(i) == 2:
            i, j = i
            for seq in self[i]:
                seq[j] = value
        else:
            raise TypeError('Index not supported')

    def biotranslate(self, *args, **kw):
        import warnings
        msg = 'BioBasket.biotranslate() is deprecated, use translate() method'
        warnings.warn(msg, DeprecationWarning, stacklevel=2)
        for seq in self:
            seq.translate(*args, **kw)
        return self

    def complement(self):
        """
        Complementary sequences, i.e. transcription
        """
        for seq in self:
            seq.complement()
        return self

    def translate(self, *args, **kw):
        """
        Translate nucleotide sequences to amino acid sequences, see `~.cane.translate()`.

        The original translate method of the str class can be used via ``BioBasket.str.translate()``.
        """
        for seq in self:
            seq.translate(*args, **kw)
        return self

    def reverse(self, *args, **kw):
        """
        Reverse sequences
        """
        for seq in self:
            seq.reverse(*args, **kw)
        return self

    def copy(self):
        """
        Return a deep copy of the BioBasket object.
        """
        return copy.deepcopy(self)

    def countall(self, rtype='counter'):
        """
        Count letters in sequences

        This method might undergo disrupting changes or it might be removed in a later version.

        :param rtype:
          * ``'counter'`` Return `~collections.Counter` object
          * ``'prob'`` Return dictionary with normalized counts
          * ``'df'`` Return pandas DataFrame object with count, prob and tprob (total prob) fields
        """
        if rtype == 'df':
            import pandas as pd
            records = [{'id': seq.id, 'letter': letter, 'count': count}
                       for seq in self for letter, count in collections.Counter(seq.data).items()]
            df = pd.DataFrame.from_records(records)
            df['prob'] =  df.groupby('id', group_keys=False)['count'].apply(lambda c: c/c.sum())
            df['tprob']= df['count'] / df['count'].sum()
            return df
        else:
            counters = [collections.Counter(seq.data) for seq in self]
            counter = reduce(operator.add, counters)
            if rtype == 'counter':
                return counter
            elif rtype == 'prob':
                s = counter.total()
                return {k: v / s for k, v in counter.items()}

    def countplot(self, y='letter', x='count', hue='id', order=None, plot='show',
                  figsize=None, ax=None, savefigkw={}, **kw):
        """
        Create a plot of letter counts

        This method might undergo disrupting changes or it might be removed in a later version.

        Under the hood this method uses pandas and seaborn libraries.
        For a help on most arguments, see
        `seaborn.barplot() <https://seaborn.pydata.org/generated/seaborn.barplot.html#seaborn.barplot>`_.
        """
        import matplotlib.pyplot as plt
        import seaborn as sns
        df = self.countall(rtype='df')
        if order in df:
            df = df.sort_values(order, ascending=False)
            order=None
        if ax is None:
            fig = plt.figure(figsize=figsize)
            ax = fig.add_subplot(111)
        ax = sns.barplot(df, y=y, x=x, hue=hue, order=order, ax=ax, **kw)
        if plot == 'show':
            plt.show()
        elif plot == 'return':
            return ax
        else:
            import matplotlib.pyplot as plt
            fig = ax.get_figure()
            fig.savefig(plot, **savefigkw)
            plt.close(fig)

    @classmethod
    def fromobj(cls, obj, tool=None):
        """
        Create a `BioBasket` object from a Python object created by another tool.

        This method is *WIP* - work in progress.

        :param obj: The object to convert.
        :param tool: the used tool (default: autodetect,
            allowed: ``'biopython'``)
        """
        if tool is None and len(obj)>0:
            tool = _detect_tool(obj[0])
        if tool is None:
            raise ValueError('Cannot determine type of object')
        if tool == 'biopython':
            seqs = [BioSeq.fromobj(seq) for seq in obj]
            return cls(seqs)
        else:
            raise ValueError(f'Unsupported tool: {tool}')

    @staticmethod
    def fromfmtstr(in_, **kw):
        """
        Read sequences from a string
        """
        from sugar import read
        if not isinstance(in_, bytes):
            in_ = in_.encode('latin1')
        return read(io.BytesIO(in_), **kw)

    def todict(self):
        """
        Return a dictionary with sequence ids as keys and sequences as values

        .. note::
            This method is different from the `BioBasket.groupby()` method.
            Each value of the dict returned by ``todict()`` is a sequence,
            whereas each value of the dict returned by ``groupby()`` is a
            BioBasket.
        """
        return {seq.id: seq for seq in self}

    @property
    def d(self):
        """
        Alias for `BioBasket.todict()`
        """
        return self.todict()

    def match(self, *args, **kw):
        """
        Search regex and return `.BioMatchList` of matches, see `~.cane.match()`
        """
        from sugar.core.cane import BioMatchList
        matches = BioMatchList()
        for seq in self:
            m = seq.match(*args, **kw)
            if kw.get('matchall'):
                matches.extend(m)
            else:
                matches.append(m)
        return matches

    def matchall(self, *args, **kw):
        """
        Search regex and return `.BioMatchList` of all matches, see `~.cane.match()`
        """
        kw['matchall'] = True
        return self.match(*args, **kw)

    def find_orfs(self, *args, **kw):
        """
        Find ORFS in sequences, see `~.cane.find_orfs()`
        """
        return reduce(lambda orfs1, orfs2: orfs1 + orfs2,
                      [seq.find_orfs(*args, **kw) for seq in self])

    def sort(self, keys=('id',), reverse=False):
        """
        Sort sequences in-place

        :param keys: Tuple of meta keys or functions to use for sorting.
            May also be a single string or callable.
            Defaults to sorting by id.
        :param reverse: Use reversed order (default: False)

        :return: Sorted sequences

        .. rubric:: Example:

        >>> from sugar import read
        >>> seqs = read()
        >>> seqs.sort(len)  # doctest: +SKIP
        """
        from sugar.core.cane import _sorted
        self.data = _sorted(self.data, keys=keys, reverse=reverse, attr='meta')
        return self

    def groupby(self, keys=('id',)):
        """
        Group sequences

        :param keys: Tuple of meta keys or functions to use for grouping.
            May also be a single string or callable.
            By default the method groups by only id.
        :return: Nested dict structure

        .. rubric:: Example:

        >>> from sugar import read
        >>> seqs = read()
        >>> grouped = seqs.groupby()
        """
        from sugar.core.cane import _groupby
        return _groupby(self, keys, attr='meta')

    def filter(self, inplace=True, **kw):
        """
        Filter sequences

        :param \*\*kw: All kwargs need to be of the form
            ``key_op=value``, where op is one of
            the operators from the `python:operator` module.
            Additionally, the operators ``'in'`` (membership),
            ``'max'`` (alias for le)
            ``'min'`` (alias for ge) are supported.
            The different filter conditions are combined with
            the *and* operator.
        :param inplace: Whether to modify the original object.
        :return: Filtered sequences

        .. rubric:: Example:

        >>> from sugar import read
        >>> seqs = read()
        >>> seqs.filter(len_gt=9500)  # doctest: +SKIP
        """
        from sugar.core.cane import _filter
        filtered = _filter(self.data, **kw)
        if inplace:
            self.data = filtered
            return self
        else:
            return self.__class__(filtered)

    def tofmtstr(self, fmt, **kw):
        """
        Write sequences to a string of specified format, see `~.main.write()`
        """
        out = io.StringIO()
        self.write(out, fmt=fmt, **kw)
        return out.getvalue()

    def tostr(self, h=19, w=80, wid=19, wlen=4, showgc=True, add_hint=False, raw=False):
        """
        Return string with information about sequences, used by ``__str__()`` method
        """
        if raw:
            return '\n'.join(str(seq) for seq in self)
        if len(self) == 0:
            return '0 seqs in basket!'
        out = [f'{len(self)} seqs in basket']
        lenid = min(max([len(seq.id) for seq in self if seq.id], default=-2), wid)+2
        maxwlen = max(len(str(len(seq))) for seq in self)
        if 0 < wlen < 4:
            wlen = max(wlen, min(4, maxwlen))
        elif wlen in (0, None) or wlen > maxwlen:
            wlen = maxwlen
        sgclen = 11 if showgc and any(getattr(seq, 'type') == 'nt' for seq in self) else 0
        for i, seq in enumerate(self):
            if (h in (None, 0) or (h2:=(n:=len(self))-h) <= 0 or
                    not (n) // 2 - (h2+1) // 2 <= i <= (n) // 2 + (h2) // 2):
                id_ = ('' if not seq.id else
                       seq.id if wid in (None, 0) or len(seq.id) <= wid else
                       seq.id[:wid-3] + '...').ljust(lenid)
                data = seq.data if w in (None, 0) or len(seq) <= w-lenid-wlen-3-sgclen else seq.data[:w-lenid-wlen-6-sgclen] + '...'
                if showgc and getattr(seq, 'type') == 'nt':
                    data = data + f'  GC:{100*seq.gc:.2f}%'
                out.append(id_ + _si_format(len(seq), max(4, wlen)).rjust(wlen) + '  ' + data)
            elif i == len(self) // 2:
                out.append('...')
        if add_hint:
            out.append('  customize output with BioBasket.tostr() method')
        return '\n'.join(out)

    def toobj(self, tool=None):
        """
        Convert the object to an object of another bioinformatics Python tool.

        *WIP* - work in progress.

        :param tool: Only ``'biopython'`` is supported.
        """
        if tool == 'biopython':
            return [seq.toobj(tool) for seq in self]
        else:
            raise ValueError(f'Unsupported tool: {tool}')

    def write(self, fname, fmt=None, **kw):
        """
        Write sequences to file, see `~.main.write()`
        """
        from sugar._io import write
        write(self, fname, fmt=fmt, **kw)

    # def consensus(self, gap='-'):
    #     n = len(self)
    #     data = [seq.data for seq in self]
    #     cons = []
    #     perc = []
    #     percnongaps = []
    #     for nt in zip(*data):
    #         max_count = 0
    #         nt = str(nt)
    #         num_gaps = nt.count(gap)
    #         for letter in list(set(nt) - {gap}):
    #             count = nt.count(letter)
