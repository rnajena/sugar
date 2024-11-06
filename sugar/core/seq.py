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
from sugar.core.fts import Feature, FeatureList, Location, LocationTuple
from sugar.core.meta import Attr, Meta


CODES_INV = {frozenset(v): k for k, v in CODES.items()}
COMPLEMENT = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G', '.': '.', '-': '-'}
COMPLEMENT_ALL = {c: CODES_INV[frozenset(COMPLEMENT[nt] for nt in nts)] for c, nts in CODES.items()}
COMPLEMENT_TRANS = str.maketrans(COMPLEMENT_ALL)


class _Sliceable_GetItem():
    def __init__(self, obj, **kw):
        self.obj = obj
        self.kw = kw

    def __getitem__(self, i):
        return self.obj._getitem(i, **self.kw)


class _BioSeqStr():
    """
    Helper class to hold all string methods in the `BioSeq.str` namespace

    :meta public:
    """
    def __init__(self, parent):
        self.__parent = parent

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

    def find(self, sub, start=0, end=sys.maxsize):
        return self.__parent.data.find(str(sub), start, end)

    def index(self, sub, start=0, end=sys.maxsize):
        return self.__parent.data.index(str(sub), start, end)

    def isalpha(self):
        return self.__parent.data.isalpha()

    def isascii(self):
        return self.__parent.data.isascii()

    def islower(self):
        return self.__parent.data.islower()

    def isupper(self):
        return self.__parent.data.isupper()

    def ljust(self, width, *args):
        self.__parent.data = self.__parent.data.ljust(width, *args)
        return self.__parent

    def lower(self):
        self.__parent.data = self.__parent.data.lower()
        return self.__parent

    def lstrip(self, chars=None):
        if chars is not None:
            chars = str(chars)
        self.__parent.data = self.__parent.data.lstrip(chars)
        return self.__parent

    @staticmethod
    def maketrans(*args):
        return str.maketrans(*args)

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

    def translate(self, *args):
        self.__parent.data = self.__parent.data.translate(*args)
        return self.__parent

    def upper(self):
        self.__parent.data = self.__parent.data.upper()
        return self.__parent


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
        return f'{type(self).__name__}({repr(self.data)}, meta=dict({metastr}))'

    def __eq__(self, string):
        if isinstance(string, BioSeq):
            return self.data == string.data and self.meta == string.meta
        return self.data == string

    def __lt__(self, other):
        if isinstance(other, BioSeq):
            return self.meta.get('id', '') < other.meta.get('id', '')
        msg = f"'<' not supported between instances of '{type(self).__name__}' and '{type(other).__name__}'"
        raise TypeError(msg)

    def __len__(self):
        return len(self.data)

    def __setitem__(self, index, value):
        l = list(self.data)
        l[index] = str(value)
        self.data = ''.join(l)

    def __add__(self, other):
        if isinstance(other, BioSeq) and self.meta != other.meta:
            warn('Join two BioSeq objects with different meta data')
        return self.__class__(self.data + str(other), meta=self.meta)

    def __iadd__(self, other):
        if isinstance(other, BioSeq) and self.meta != other.meta:
            warn('Join two BioSeq objects with different meta data')
        self.data = self.data + str(other)
        return self

    def __radd__(self, other):
        return self.__class__(str(other) + self.data, meta=self.meta)

    @property
    def str(self):
        """
        Namespace holding all available string methods,
        see `_BioSeqStr` for available methods
        and `python:str` for documentation of the methods

        .. rubric:: Example:

        >>> seq = read()[0]
        >>> seq.str.find('ATG')  # Use string method
        30
        """
        return _BioSeqStr(self)

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
        self.fts = self.fts + FeatureList(fts)
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
        """Return sliceable object to support in-place slicing

        Deprecated: Use getitem() or sl attribute.
        """
        msg = 'BioSeq.i is deprecated, use geitem() method or sl attribute'
        warnings.warn(msg, DeprecationWarning, stacklevel=2)
        return _Sliceable_GetItem(self, inplace=True)

    def rc(self, update_fts=False):
        """
        Reverse complement, alias for ``BioSeq.reverse().complement()``
        """
        self.reverse().complement()
        if update_fts:
            self.fts = self.fts.reverse(lenseq=len(self))
        return self

    def __getitem__(self, index):
        return self._getitem(index)

    def sl(self, **kw):
        """
        Method allowing to slice the `BioSeq` object with non-default options.

        If you want to use the default options, you can slice the BioSeq object directly.
        For non-default options, slice the sliceable object returned by this method.

        :param bool inplace:
            The subsequence is not only returned, but the original
            sequence is modified in-place (default: False)
        :param str gap:
            gaps of the given characters will be accounted for when
            slicing the sequence (default: gaps will not be accounted for)

        .. rubric:: Slicing options:

        The slice specifies which part of the sequence is returned and
        is defined inside the square brackets ``[]``
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

        .. rubric:: Example:

        >>> from sugar import read
        >>> seq = read()[0]
        >>> print(seq[:5])  # use direct slicing for default options
        ACCTG
        >>> print(seq[4])
        G
        >>> print(seq['cds'][:3])
        ATG
        >>> print(seq.sl(inplace=True, gap='-')[:5:2])  # non-default options
        ACG
        >>> print(seq)  # was modified in-place
        ACG
        """
        return _Sliceable_GetItem(self, **kw)

    def _slice_locs(self, locs, splitter=None, filler=None, gap=None, update_fts=False):
        # TODO document
        # Merge the sequences corresponding to the ordered locations
        sub_seqs = []
        prev_loc = None
        for i, loc in enumerate(locs):
            add_seq = self.sl(gap=gap)[loc.start:loc.stop]
            if loc.strand == '-':
                add_seq = add_seq.rc()
            if filler is not None and prev_loc is not None:
                if loc.strand == '-':
                    num = prev_loc.start - loc.stop
                else:
                    num = loc.start - prev_loc.stop
                if num > 0:
                    sub_seqs.append(num * filler)
            sub_seqs.append(add_seq.data)
            prev_loc = loc
        sub_seq = self[:0]
        if splitter is None:
            splitter = ''
        sub_seq.data = splitter.join(sub_seqs)
        if update_fts:
            start, stop = locs.range[0]
            fts = FeatureList()
            for loc in locs:
                fts.extend(self.fts.slice(loc.start, loc.stop, rel=start))
            if locs[0].strand == '-':
                fts = fts.reverse(lenseq=stop)
            sub_seq.fts = fts
        return sub_seq

    def _getitem(self, index, inplace=False, gap=None, update_fts=False, **kw):
        """
        Slice the sequence and return a subsequence

        This is the method which is called if you slice with ``BioSeq[]`` syntax
        or with ``BioSeq.sl()[]`` syntax
        If you want to use non-default options call this method directly,
        or by the `BioSeq.sl` attribute.
        """
        # TODO: add update_fts tests
        if not isinstance(index, (int, slice)):
            if isinstance(index, str):
                index = self.fts.get(index)
                if index is None:
                    raise TypeError('Feature not found')
            if isinstance(index, Location):
                index = LocationTuple([Location])
            elif isinstance(index, Feature):
                index = index.locs
            if isinstance(index, LocationTuple):
                subseq = self._slice_locs(index, gap=gap, update_fts=update_fts, **kw)
            else:
                raise TypeError('Index not supported')
        else:
            if gap is not None:
                # from bisect import bisect
                nogaps = [i for i, nt in enumerate(self.data) if nt not in gap]
                adj = lambda i: nogaps[i] if i is not None and i < len(nogaps) else None
                if isinstance(index, int):
                    index = adj(index)
                elif isinstance(index, slice):
                    index = slice(adj(index.start), adj(index.stop), index.step)
            subseq = self.__class__(self.data[index], meta=self.meta)
            if update_fts:
                if isinstance(index, int):
                    start, stop = index, index + 1
                elif isinstance(index, slice):
                    start, stop = index.start, index.stop
                    if index.step != 1:
                        raise ValueError('update_fts for slices only supported with step==1')
                subseq.fts = self.fts.slice(start, stop, rel=start)
        if inplace:
            self.data = subseq.data
        return subseq

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

    def tostr(self, **kw):
        """
        Return nice string, see `BioBasket.tostr()`
        """
        kw.setdefault('add_header', False)
        return BioBasket([self]).tostr(**kw)

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

    # Implement all variants of &, |, -, ^
    def __and__(self, other):
        return self.__class__([seq for seq in self if seq in other])

    def __rand__(self, other):
        return self & other

    def __iand__(self, other):
        self.data = [seq for seq in self if seq in other]
        return self

    def __or__(self, other):
        return self + [seq for seq in other if seq not in self]

    def __ror__(self, other):
        return self | other

    def __ior__(self, other):
        self.data += [seq for seq in other if seq not in self]
        return self

    def __sub__(self, other):
        return self.__class__([seq for seq in self if seq not in other])

    def __rsub__(self, other):
        return self.__class__(other) - self

    def __isub__(self, other):
        self.data = [seq for seq in self if seq not in other]
        return self

    def __xor__(self, other):
        return (self | other) - (self & other)

    def __rxor__(self, other):
        return self ^ other

    def __ixor__(self, other):
        self.data = (self ^ other).data
        return self

    @property
    def str(self):
        """
        Namespace holding all available string methods.

        The `BioBasket.str` methods call the corresponding `BioSeq.str` methods under the hood
        and return either the altered `BioBasket` object or a list with results.
        See `_BioSeqStr` for available methods
        and `python:str` for documentation of the methods.

        .. rubric:: Example:

        >>> seqs = read()
        >>> seqs.str.find('ATG')  # Use string method
        [30, 12]
        """
        return _BioBasketStr(self)

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
                seq.fts.sort()
        if len(fts) > 0:
            missing_ids = ', '.join(fts.keys())
            warn(f'Features for seqids {missing_ids} could not be '
                 'attached to any sequence')

    def rc(self, **kw):
        """
        Reverse complement, alias for ``BioBasket.reverse().complement()``
        """
        for seq in self:
            seq.rc(**kw)
        return self

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
        return self._getitem(i)

    def sl(self, **kw):
        """
        Method allowing to slice the `BioBasket` object with non-default options.

        If you want to use the default options, you can slice the BioBasket object directly.
        For non-default options, slice the sliceable object returned by this method.

        :param \*\*kw:
            All kwargs are documented in `BioSeq.sl()`.

        .. rubric:: Slice options:

        The slice specifies which part of the sequence(s) are returned and
        is defined inside the square brackets ``[]``
        The following options are supported.

        int
            Returns a `BioSeq` from the basket
        slice
            Returns a new `BioBasket` object with a subset of the sequences
        str,feature,location
            Returns a new `BioBasket` object with updated sequences inside, see `BioSeq.sl()`
        (int, object)
            Returns a `BioSeq` from the basket and slices it with the object, see `BioSeq.sl()`
        (slice, object)
            Returns a new `BioBasket` object with a subset of the sequences which are replaced
            by subsequences according to `BioSeq.sl()`

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
        """
        return _Sliceable_GetItem(self, **kw)

    def _getitem(self, i, **kw):
        """
        Slice sequences

        This is the method which is called if you slice with ``BioBasket[]`` syntax.
        If you want to use non-default options call this method directly,
        or by the `BioBasket.sl` attribute.
        """
        if isinstance(i, int):
            return self.data[i]
        elif isinstance(i, slice):
            seqs = self.__class__(self.data[i], meta=self.meta)
        elif isinstance(i, (str, Feature, Location)):
            seqs = self.__class__(self.data, meta=self.meta)
            seqs.data = [seq._getitem(i, **kw) for seq in seqs.data]
        elif len(i) == 2:
            i, j = i
            if isinstance(i, int):
                return self.data[i]._getitem(j, **kw)
            elif isinstance(i, slice):
                seqs = self.__class__(self.data[i], meta=self.meta)
                seqs.data = [seq._getitem(j, **kw) for seq in seqs.data]
            else:
                raise TypeError('Index not supported')
        else:
            raise TypeError('Index not supported')
        return seqs

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
            the *and* operator. If you need *or*, call filter twice
            and combine the results with ``|`` operator, e.g.
            ``seqs.filter(...) | seqs.filter(...)``
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

    def tostr(self, h=19, w=80, wid=19, wlen=4, showgc=True,
              add_hint=False, raw=False, add_header=True):
        """
        Return string with information about sequences, used by ``__str__()`` method
        """
        if raw:
            return '\n'.join(str(seq) for seq in self)
        if len(self) == 0:
            return '0 seqs in basket!' * add_header
        out = [f'{len(self)} seqs in basket'] * add_header
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
