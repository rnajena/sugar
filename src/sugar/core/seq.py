# (C) 2024, Tom Eulenfeld, MIT license
"""
Sequence related classes, `.BioSeq`, `.BioBasket`
"""

import collections
import collections.abc
import copy
from functools import reduce
import io
import math
import sys
from warnings import warn

from sugar.data import CODES
from sugar.core.fts import Feature, FeatureList, Location, LocationTuple
from sugar.core.meta import Attr, Meta
from sugar.core.util import _add_inplace_doc


CODES_INV = {frozenset(v): k for k, v in CODES.items()}
COMPLEMENT = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G', '.': '.', '-': '-'}
COMPLEMENT_ALL = {c: CODES_INV[frozenset(COMPLEMENT[nt] for nt in nts)] for c, nts in CODES.items()}
COMPLEMENT_TRANS = str.maketrans(COMPLEMENT_ALL)


class _Sliceable_GetItem():
    def __init__(self, obj, method='_getitem', **kw):
        self.kw = kw
        self.call = getattr(obj, method)

    def __getitem__(self, i):
        return self.call(i, **self.kw)


class _BioSeqStr():
    """
    Helper class to hold all string methods in the `BioSeq.str` namespace.

    The methods modify the data in-place, if applicable,
    which is different from the behavior of the original string methods.

    See `str` for documentation of the methods.

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
        """
        The ``ljust()`` and ``rjust()`` methods can be used to fill up an alignment with gaps.

        Example: ``seqs.ljust(500, '-')``
        """
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
    either the modified `BioBasket` object or a list of results.

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


class BioSeq():
    """
    Class holding sequence data and metadata, exposing bioinformatics methods.

    Most methods work in-place by default, but return the BioSeq object again.
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
        and `str` for documentation of the methods

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

        If you want to set all features, use the `BioSeq.fts` attribute.

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

    _add_inplace_doc
    def rc(self, update_fts=False):
        """
        Reverse complement, alias for ``BioSeq.reverse().complement()``
        """
        self.reverse().complement()
        if update_fts:
            self.fts = self.fts.rc(seqlen=len(self))
        return self

    def __getitem__(self, index):
        return self._getitem(index)

    def sl(self, **kw):
        """
        Method that allows you to slice the `BioSeq` object with non-default options.

        If you want to use the default options, you can slice the BioSeq object directly.
        For non-default options, slice the sliceable object returned by this method.

        :param bool inplace:
            Not only will the subsequence be returned,
            but the original sequence will be modified in-place
            (default: False).
        :param str gap:
            Gaps of the given characters are taken into account
            when slicing the sequence
            (default: gaps are not taken into account)

        .. rubric:: Slicing options:

        The slice specifies which part of the sequence is returned, and
        is defined inside the square brackets ``[]``
        The following types are supported.

        int,slice
            The location is specified by int or slice
        `.Location`
            specified by location
        `.Feature`
            specified by feature
        str
            Position of the first feature of the given type, e.g. ``'cds'``
            will return the sequence with the first coding sequence.

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
        >>> print(seq)  # has been modified in-place
        ACG
        """
        return _Sliceable_GetItem(self, **kw)

    def slindex(self, gap=None):
        """
        Method that translates an index to account for gaps

        .. rubric:: Example:

        >>> from sugar import BioSeq
        >>> seq = BioSeq('ATG---GGA')
        >>> print(seq)
        ATG---GGA
        >>> print(seq[1:5])
        TG--
        >>> print(seq.sl(gap='-')[1:5])
        TG---GG
        >>> print(seq.slindex(gap='-')[1:5])
        slice(1, 8, None)
        >>> print(seq[seq.slindex(gap='-')[1:5]])
        TG---GG
        """
        return _Sliceable_GetItem(self, method='_getindex', gap=gap)

    def _slice_locs(self, locs, splitter=None, filler=None, gap=None, update_fts=False):
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
            if splitter is not None and prev_loc is not None:
                sub_seqs.append(splitter)
            sub_seqs.append(add_seq.data)
            prev_loc = loc
        sub_seq = BioSeq(''.join(sub_seqs), meta=self.meta.copy())
        if update_fts:
            if len(locs) > 1:
                raise ValueError('update_fts is only allowed for locs wit a single Location. Sorry.')
            start, stop = locs.range
            fts = FeatureList()
            for loc in locs:
                fts.extend(self.fts.slice(loc.start, loc.stop, rel=start))
            if locs[0].strand == '-':
                fts = fts.rc(seqlen=stop-start)
            sub_seq.fts = fts
        return sub_seq

    def _getindex(self, index, gap=None):
        if not isinstance(index, (int, slice)):
            if isinstance(index, str):
                ft = self.fts.get(index)
                if ft is None:
                    raise ValueError(f'Feature of type {index} not found')
                index = ft
            if isinstance(index, Location):
                index = LocationTuple([Location])
            elif isinstance(index, Feature):
                index = index.locs
            if isinstance(index, LocationTuple):
                start, stop = index.range
                return self.slindex(gap=gap)[start:stop]
            else:
                raise TypeError(f"Index of type '{type(index).__name__}' not supported")
        # from bisect import bisect
        nogaps = [i for i, nt in enumerate(self.data) if nt not in gap]
        adj = lambda i: nogaps[i] if i is not None and i < len(nogaps) else None
        adj_stop = lambda i: nogaps[i-1]+1 if i is not None and 0 < i < len(nogaps) else 0 if i == 0 else None
        if isinstance(index, int):
            index = adj(index)
        elif isinstance(index, slice):
            index = slice(adj(index.start), adj_stop(index.stop), index.step)
        return index

    def _getitem(self, index, inplace=False, gap=None, update_fts=False, **kw):
        """
        Slice the sequence and return a subsequence

        This is the method which is called if you slice with ``BioSeq[]`` syntax
        or with ``BioSeq.sl()[]`` syntax
        If you want to use non-default options call this method directly,
        or by the `BioSeq.sl` attribute.
        """
        # TODO: doc for splitter and filler
        if not isinstance(index, (int, slice)):
            if isinstance(index, str):
                ft = self.fts.get(index)
                if ft is None:
                    raise ValueError(f'Feature of type {index} not found')
                index = ft
            if isinstance(index, Location):
                index = LocationTuple([index])
            elif isinstance(index, Feature):
                index = index.locs
            if isinstance(index, LocationTuple):
                subseq = self._slice_locs(index, gap=gap, update_fts=update_fts, **kw)
            else:
                raise TypeError(f"Index of type '{type(index).__name__}' not supported")
        else:
            subseq = self.__class__(
                self.data[index if gap is None else self._getindex(index, gap=gap)],
                meta=self.meta
                )
            if update_fts:
                if isinstance(index, int):
                    start, stop = index, index + 1
                elif isinstance(index, slice):
                    start, stop = index.start, index.stop
                    if index.step not in (None, 1):
                        raise ValueError('update_fts for slices only supported with step==1')
                subseq.fts = self.fts.slice(start, stop, rel=start or 0)
        if inplace:
            self.data = subseq.data
        return subseq

    @_add_inplace_doc
    def complement(self):
        """
        Complementary sequence, i.e. transcription
        """
        if 'U' in self.data:
            self.str.replace('U', 'T').str.translate(COMPLEMENT_TRANS).str.replace('T', 'U')
        else:
            self.str.translate(COMPLEMENT_TRANS)
        return self

    def countall(self, *args, **kw):
        return BioBasket([self]).countall(*args, **kw)

    def countplot(self, *args, hue=None, **kw):
        return BioBasket([self]).countplot(*args, hue=hue, **kw)

    def copy(self):
        """
        Return a deep copy of the object
        """
        return copy.deepcopy(self)

    @classmethod
    def frombiopython(cls, obj):
        """
        Create a `.BioSeq` object from a biopython_ `~Bio.SeqRecord.SeqRecord` or `~Bio.Seq.Seq` object.

        :param obj: The object to convert.

        .. note::
            BioPython Features in the ``SeqRecord.features`` attribute are automatically converted.
        """
        from sugar.core._adapter import biopython2seq
        return biopython2seq(obj, cls=cls)

    @classmethod
    def frombiotite(cls, obj):
        """
        Create a `.BioSeq` object from a biotite_ sequence object.

        :param obj: The object to convert.
        """
        from sugar.core._adapter import biotite2seq
        return biotite2seq(obj, cls=cls)

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

    def plot_ftsviewer(self, *args, **kw):
        """
        Plot features of the sequence using DNAFeaturesViewer_, see `~.imaging.ftsviewer.plot_ftsviewer()`

        .. note::
            Using `BioSeq <.BioSeq.plot_ftsviewer>` or `.BioBasket.plot_ftsviewer()`
            over `.FeatureList.plot_ftsviewer()` has the advantage,
            that sequence lengths are used automatically.
        """
        return BioBasket([self]).plot_ftsviewer(*args, **kw)

    def tostr(self, **kw):
        """
        Return a nice string, see `BioBasket.tostr()`
        """
        kw.setdefault('add_header', False)
        return BioBasket([self]).tostr(**kw)

    def tofmtstr(self, fmt, **kw):
        """
        Write object to a string of given format, see `~.main.write()`
        """
        return BioBasket([self]).tofmtstr(fmt, **kw)

    def tobiopython(self):
        """
        Convert BioSeq to biopython_ `~Bio.SeqRecord.SeqRecord` instance

        Attached ``BioSeq.fts`` features are automatically converted.
        """
        from sugar.core._adapter import seq2biopython
        return seq2biopython(self)

    def tobiotite(self, **kw):
        """
        Convert BioSeq to biotite_ `~biotite.sequence.NucleotideSequence` or `~biotite.sequence.ProteinSequence` instance

        :param str type: ``'nt'`` creates a `~biotite.sequence.NucleotideSequence` instance,
            ``'aa'`` creates a `~biotite.sequence.ProteinSequence` instance,
            by default the class is inferred from the sequence itself.
        :param str gap: Gap characters that must be removed from the sequence string.
        :param bool warn: Whether to warn if gap characters have been removed, default is True.
        """
        from sugar.core._adapter import seq2biotite
        return seq2biotite(self, **kw)

    def toftsviewer(self, **kw):
        r"""
        Convert features of this sequence to DNAFeaturesViewer_ `~dna_features_viewer.GraphicRecord`

        See `.FeatureList.toftsviewer`.
        """
        return self.fts.toftsviewer(seq=self, **kw)

    @_add_inplace_doc
    def reverse(self):
        """
        Reverse the sequence
        """
        self.data = self.data[::-1]
        return self

    @_add_inplace_doc
    def translate(self, *args, update_fts=False, **kw):
        """
        Translate nucleotide sequence to amino acid sequence, see `~.cane.translate()`.

        The original translate method of the str class can be used via ``BioBasket.str.translate()``.
        """
        from sugar.core.cane import translate
        self.data = translate(self.data, *args, **kw)
        self.type = 'aa'
        if update_fts:
            fts = self.fts.slice(0, len(self) * 3)
            for ft in fts:
                for loc in ft.locs:
                    loc.start = loc.start // 3
                    loc.stop = loc.stop // 3
            self.fts = fts
        return self

    def write(self, fname=None, fmt=None, **kw):
        """
        Write sequence to file, see `~.main.write()`
        """
        return BioBasket([self]).write(fname=fname, fmt=fmt, **kw)


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
    It has useful bioinformatics methods attached to it.

    The list itself is stored in the ``data`` property.
    The BioBasket object may also have a metadata
    attribute.
    """
    def __init__(self, data=None, meta=None):
        if isinstance(data, (str, BioSeq)):
            warn(f'{self.__class__.__name__} object initialized with a sequence. '
                 'To hide this warning initialize with a list of sequences.')
            data = [data]
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
        and return either the modified `BioBasket` object or a list of results.
        See `_BioSeqStr` for available methods
        and `str` for method documentation.

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
        fts = FeatureList(value).groupby('seqid')
        for seq in self:
            if seq.id in fts:
                seq.fts = fts.pop(seq.id)
            else:
                try:
                    del seq.meta.fts
                except KeyError:
                    pass
        if len(fts) > 0:
            missing_ids = ', '.join(fts.keys())
            warn(f'Features for seqids {missing_ids} could not be '
                 'attached to any sequence')

    def add_fts(self, fts):
        """
        Add some features to the feature list of the corresponding sequences.

        If you want to set all features, use the `BioBasket.fts` attribute.

        :param fts: features to add
        """
        fts = FeatureList(fts).groupby('seqid')
        for seq in self:
            if seq.id in fts:
                seq.fts = seq.fts + fts.pop(seq.id)
                seq.fts.sort()
        if len(fts) > 0:
            missing_ids = ', '.join(fts.keys())
            warn(f'Features for seqids {missing_ids} could not be '
                 'attached to any sequence')

    @_add_inplace_doc
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
        r"""
        Method that allows you to slice the `BioBasket` object with non-default options.

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
        elif isinstance(i, (str, Feature, LocationTuple, Location)):
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
                raise TypeError('First dimension of 2D index must be of type int or slice')
        else:
            raise TypeError(f"Index of type '{type(i).__name__}' not supported")
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
            raise TypeError(f"Index of type '{type(i).__name__}' not supported")

    @_add_inplace_doc
    def complement(self):
        """
        Complementary sequences, i.e. transcription
        """
        for seq in self:
            seq.complement()
        return self

    @_add_inplace_doc
    def translate(self, *args, **kw):
        """
        Translate nucleotide sequences to amino acid sequences, see `~.cane.translate()`.

        The original translate method of the str class can be used via ``BioBasket.str.translate()``.
        """
        for seq in self:
            seq.translate(*args, **kw)
        return self

    @_add_inplace_doc
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

    def countall(self, rtype='counter', k=1):
        """
        Count letters in sequences

        This method may undergo disrupting changes or it may be removed in a later release.

        :param rtype:
          * ``'counter'`` Return `~collections.Counter` object
          * ``'prob'`` Return dictionary with normalized counts
          * ``'df'`` Return pandas DataFrame object with count, prob and tprob (total prob) fields
        """
        from operator import add
        if rtype == 'df':
            import pandas as pd
            records = [{'id': seq.id, 'word': word, 'count': count}
                       for seq in self for word, count in collections.Counter(
                           seq.data if k == 1 else
                           [reduce(add, kmer) for kmer in zip(*[seq.data[i:] for i in range(k)])]
                           ).items()]
            df = pd.DataFrame.from_records(records)
            df['prob'] =  df.groupby('id', group_keys=False)['count'].apply(lambda c: c/c.sum())
            df['tprob']= df['count'] / df['count'].sum()
            return df
        else:
            counters = [
                collections.Counter(
                    seq.data if k == 1 else
                    [reduce(add, kmer) for kmer in zip(*[seq.data[i:] for i in range(k)])]
                    ) for seq in self]
            counter = reduce(add, counters)
            if rtype == 'counter':
                return counter
            elif rtype == 'prob':
                s = counter.total()
                return {k: v / s for k, v in counter.items()}

    def countplot(self, y='word', x='count', hue='id', order=None, plot='show',
                  figsize=None, ax=None, savefigkw={}, **kw):
        """
        Create a plot of letter counts

        This method may undergo disruptive changes, or it may be removed in a later release.

        Under the hood this method uses the pandas and seaborn libraries.
        For a help on most of the arguments, see
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
    def frombiopython(cls, obj):
        """
        Create a `.BioBasket` object from a list of biopython_ `~Bio.SeqRecord.SeqRecord` or `~Bio.Seq.Seq` objects.

        :param obj: The object to convert, can also be a `~Bio.Align.MultipleSeqAlignment` object.

        .. note::
            BioPython Features in the ``SeqRecord.features`` attribute are automatically converted.
        """
        from sugar.core._adapter import biopython2seqs
        return biopython2seqs(obj, cls=cls)

    @classmethod
    def frombiotite(cls, obj):
        """
        Create a `.BioBasket` object from a list of biotite_ sequence objects.

        :param obj: The object to convert, can also be a biotite ``Alignment`` object.
        """
        from sugar.core._adapter import biotite2seqs
        return biotite2seqs(obj, cls=cls)

    @staticmethod
    def fromfmtstr(in_, fmt=None, **kw):
        """
        Read sequences from a string
        """
        from sugar import read
        if not isinstance(in_, bytes):
            in_ = in_.encode('latin1')
        return read(io.BytesIO(in_), fmt=fmt, **kw)

    def todict(self):
        """
        Return a dictionary with sequence ids as keys and sequences as values

        .. note::
            This method is different from the `BioBasket.groupby()` method.
            Each value of the dict returned by ``todict()`` is a sequence,
            while each value of the dict returned by ``groupby()`` is a
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

    @_add_inplace_doc
    def sort(self, keys=('id',), reverse=False):
        """
        Sort sequences in-place

        :param keys: Tuple of meta keys or functions to use for sorting.
            Can also be a single string or a callable.
            Defaults to sorting by id.
        :param reverse: Use reverse order (default: False)

        :return: Sorted sequences

        .. rubric:: Example:

        >>> from sugar import read
        >>> seqs = read()
        >>> seqs.sort(len)  # doctest: +SKIP
        """
        from sugar.core.cane import _sorted
        self.data = _sorted(self.data, keys=keys, reverse=reverse, attr='meta')
        return self

    def groupby(self, keys=('id',), flatten=False):
        """
        Group sequences

        :param keys: Tuple of meta keys or functions to use for grouping.
            Can also be a single string or a callable.
            By default, the method groups only by id.
        :return: Nested dict structure

        .. rubric:: Example:

        >>> from sugar import read
        >>> seqs = read()
        >>> grouped = seqs.groupby()
        """
        from sugar.core.cane import _groupby
        return _groupby(self, keys, attr='meta', flatten=flatten)

    def select(self, inplace=False, **kw):
        r"""
        Select sequences

        :param \*\*kw: All kwargs must be of the form
            ``key_op=value``, where op is one of
            the operators from the `operator` module.
            Additionally, the operator ``'in'`` (membership) is supported.
            The different select conditions are combined with
            the *and* operator. If you need *or*, call select twice
            and combine the results with the ``|`` operator, e.g.
            ``seqs.select(...) | seqs.select(...)``
        :param inplace: Whether to modify the original object (default: False)
        :return: Selected sequences

        .. rubric:: Example:

        >>> from sugar import read
        >>> seqs = read()
        >>> seqs2 = seqs.select(len_gt=9500)  # Select sequences with length > 9500
        """
        from sugar.core.cane import _select
        selected = _select(self.data, **kw)
        if inplace:
            self.data = selected
            return self
        else:
            return self.__class__(selected)

    def tofmtstr(self, fmt, **kw):
        """
        Write sequences to a string of the specified format, see `~.main.write()`
        """
        return self.write(None, fmt)

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

    def tobiopython(self, *, msa=False):
        """
        Convert the BioBasket to a list of biopython_ `~Bio.SeqRecord.SeqRecord` objects

        :param bool msa: Return a biopython `~Bio.Align.MultipleSeqAlignment` object instead of a list.

        Attached ``BioSeq.fts`` features are not converted.
        """
        from sugar.core._adapter import seqs2biopython
        return seqs2biopython(self, msa=msa)

    def tobiotite(self, **kw):
        """
        Convert BioBasket to a list of biotite_ `~biotite.sequence.NucleotideSequence` or `~biotite.sequence.ProteinSequence` instance

        :param str type: ``'nt'`` creates a `~biotite.sequence.NucleotideSequence` instance,
            ``'aa'`` creates a `~biotite.sequence.ProteinSequence` instance,
            by default the class is inferred from the sequence itself.
        :param bool msa: Return a biotite ``Alignment`` object instead of a list, default is False
        :param str gap: Gap characters that must be removed from the sequence strings.
        :param bool warn: Wether to warn if gap characters have been removed, default is True,
            not used with ``msa=True``
        """
        from sugar.core._adapter import seqs2biotite
        return seqs2biotite(self, **kw)

    def write(self, fname=None, fmt=None, **kw):
        """
        Write sequences to file, see `~.main.write()`
        """
        from sugar._io import write
        return write(self, fname=fname, fmt=fmt, **kw)

    def plot_alignment(self, *args, **kw):
        """
        Plot an alignment, see `~.imaging.alignment.plot_alignment()`
        """
        from sugar.imaging import plot_alignment
        return plot_alignment(self, *args, **kw)

    def plot_ftsviewer(self, *args, **kw):
        """
        Plot features of the sequences using DNAFeaturesViewer_, see `~.imaging.ftsviewer.plot_ftsviewer()`

        .. note::
            Using `BioSeq <.BioSeq.plot_ftsviewer>` or `.BioBasket.plot_ftsviewer()`
            over `.FeatureList.plot_ftsviewer()` has the advantage,
            that sequence lengths are used automatically.
        """
        from sugar.imaging import plot_ftsviewer
        return plot_ftsviewer(self.fts, *args, seqs=self, **kw)

    def merge(self, spacer='', update_fts=False, keys=('id',)):
        data = []
        for group in self.groupby(keys, flatten=True).values():
            seq = group[0]
            for oseq in group[1:]:
                seq.data = seq.data + spacer + oseq.data
                if update_fts:
                    nfts = []
                    for ft in oseq.fts:
                        rel = len(seq) - len(oseq)
                        locs = [Location(l.start + rel, l.stop + rel, l.strand, l.defect) for l in ft.locs]
                        nfts.append(Feature(locs=locs, meta=ft.meta))
                    seq.add_fts(nfts)
            data.append(seq)
        self.data = data
        return self


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
