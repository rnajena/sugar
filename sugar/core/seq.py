# (C) 2024, Tom Eulenfeld, MIT license

import collections
import collections.abc
import copy
from functools import reduce
import operator
import sys
import io

from sugar.data import CODES
from sugar.core.fts import Feature, FeatureList
from sugar.core.meta import Attr, Meta


CODES_INV = {frozenset(v): k for k, v in CODES.items()}
COMPLEMENT = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G', '.': '.', '-': '-'}
COMPLEMENT_ALL = {c: CODES_INV[frozenset(COMPLEMENT[nt] for nt in nts)] for c, nts in CODES.items()}
COMPLEMENT_TRANS = str.maketrans(COMPLEMENT_ALL)


# class Feature(Attr):
#     def _slice(self):
#         stride = getattr(self, 'stride', 1)
#         return slice(self.start, self.stop, stride)


# class FeatureList(collections.UserList):
#     def __init__(self, data=None):
#         super().__init__(data)

#     def __str__(self):
#         return self.tostr()

#     def _repr_pretty_(self, p, cycle):
#         if cycle:
#             p.text('...')
#         else:
#             p.text(str(self))

#     def tostr(self, w=80, wt=12, wl=12, wle=6, exclude_features=()):
#         out = []
#         for ft in self:
#             t = getattr(ft, 'type', '')
#             if t in exclude_features:
#                 continue
#             exclude_keys = ('start', 'stop', 'stride', 'type', 'loc', 'translation')
#             ftstr = ';'.join(f'{k}={v}' for k, v in vars(ft).items()
#                              if k not in exclude_keys)
#             l = getattr(ft, 'loc', '')
#             le = ft.stop - ft.start if getattr(ft, 'start', None) is not None else '?'
#             le = f'({le})'
#             ftstr = f'{t:>{wt}} {l:<{wl}} {le:<{wle}}  {ftstr}'
#             if w and len(ftstr) > w:
#                 ftstr = ftstr[:w-3] + '...'
#             out.append(ftstr)
#         return '\n'.join(out)

#     def get(self, type_):
#         for ft in self.data:
#             if ft.type == type_:
#                 return ft

#     def geta(self, type_):
#         features = []
#         for ft in self.data:
#             if ft.type == type_:
#                 features.append(ft)
#         return FeatureList(features)




class _Slicable():
    def __init__(self, seq):
        self.seq = seq

    def __getitem__(self, i):
        self.seq.data = self.seq[i].data
        return self.seq


class MutableMetaString(collections.abc.Sequence):
    def __init__(self, data, id='', meta=None, type=None):
        self.data = str(data).upper()
        if hasattr(data, 'meta'):
            meta = data.meta
        elif 'meta' in data:
            meta = data['meta']
        elif meta is None:
            meta = {}
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
        self.type = type

    @property
    def id(self):
        return self.meta.id

    @id.setter
    def id(self, value):
        self.meta.id = value

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

    ## adapted from Pythons UserString

    # def __int__(self):
    #     return int(self.data)

    # def __float__(self):
    #     return float(self.data)

    # def __complex__(self):
    #     return complex(self.data)

    #def __hash__(self):
    #    return hash(self.data)

    # def __getnewargs__(self):
    #     return (self.data[:],)

    def __eq__(self, string):
        if isinstance(string, MutableMetaString):
            return self.data == string.data and self.meta == string.meta
        return self.data == string

    def __lt__(self, string):
        if isinstance(string, MutableMetaString):
            return self.id < string.id
        self.id < ''

    def __le__(self, string):
        if isinstance(string, MutableMetaString):
            return self.id <= string.id
        self.id <= ''

    def __gt__(self, string):
        if isinstance(string, MutableMetaString):
            return self.id > string.id
        self.id > ''

    def __ge__(self, string):
        if isinstance(string, MutableMetaString):
            return self.id >= string.id
        self.id >= ''

    def __contains__(self, char):
        if isinstance(char, MutableMetaString):
            char = char.data
        return char in self.data

    def __len__(self):
        return len(self.data)

    def __getitem__(self, index):
        return self.__class__(self.data[index], meta=self.meta)

    def __setitem__(self, index, value):
        l = list(self.data)
        l[index] = value
        self.data = ''.join(l)

    def __add__(self, other):
        if isinstance(other, MutableMetaString):
            meta = self.meta if self.meta == other.meta else None
            id_ = self.id if self.id == other.id else None
            return self.__class__(self.data + other.data, id=id_, meta=meta)
        elif isinstance(other, str):
            return self.__class__(self.data + other, meta=self.meta)
        return self.__class__(self.data + str(other), meta=self.meta)

    def __iadd__(self, other):
        if isinstance(other, MutableMetaString):
            self.data = self.data + other.data
        elif isinstance(other, str):
            self.data = self.data + other
        else:
            self.data = self.data + str(other)
        return self

    def __radd__(self, other):
        if isinstance(other, str):
            return self.__class__(other + self.data, meta=self.meta)
        return self.__class__(str(other) + self.data, meta=self.meta)

    # def __mul__(self, n):
    #     return self.__class__(self.data * n, meta=self.meta)

    # def __imul__(self, n):
    #     self.data = self.data * n

    # __rmul__ = __mul__

    # def __mod__(self, args):
    #     return self.__class__(self.data % args, meta=self.meta)

    # def __rmod__(self, template):
    #     return self.__class__(str(template) % self, meta=self.meta)

    # the following methods are defined in alphabetical order:
    # def capitalize(self):
    #     self.data.capitalize()
    #     return self

    # def casefold(self):
    #     self.data = self.data.casefold()
    #     return self

    def center(self, width, *args):
        self.data = self.data.center(width, *args)
        return self

    def count(self, sub, start=0, end=sys.maxsize):
        if isinstance(sub, MutableMetaString):
            sub = sub.data
        return self.data.count(sub, start, end)

    def removeprefix(self, prefix, /):
        if isinstance(prefix, MutableMetaString):
            prefix = prefix.data
        self.data = self.data.removeprefix(prefix)
        return self

    def removesuffix(self, suffix, /):
        if isinstance(suffix, MutableMetaString):
            suffix = suffix.data
        self.data = self.data.removesuffix(suffix)
        return self

    def encode(self, encoding='utf-8', errors='strict'):
        encoding = 'utf-8' if encoding is None else encoding
        errors = 'strict' if errors is None else errors
        return self.data.encode(encoding, errors)

    def endswith(self, suffix, start=0, end=sys.maxsize):
        return self.data.endswith(suffix, start, end)

    # def expandtabs(self, tabsize=8):
    #     return self.__class__(self.data.expandtabs(tabsize))

    def find(self, sub, start=0, end=sys.maxsize):
        if isinstance(sub, MutableMetaString):
            sub = sub.data
        return self.data.find(sub, start, end)

    def format(self, /, *args, **kwds):
        return self.data.format(*args, **kwds)

    def format_map(self, mapping):
        return self.data.format_map(mapping)

    def index(self, sub, start=0, end=sys.maxsize):
        return self.data.index(sub, start, end)

    def isalpha(self):
        return self.data.isalpha()

    # def isalnum(self):
    #     return self.data.isalnum()

    def isascii(self):
        return self.data.isascii()

    # def isdecimal(self):
    #     return self.data.isdecimal()

    # def isdigit(self):
    #     return self.data.isdigit()

    # def isidentifier(self):
    #     return self.data.isidentifier()

    def islower(self):
        return self.data.islower()

    # def isnumeric(self):
    #     return self.data.isnumeric()

    # def isprintable(self):
    #     return self.data.isprintable()

    # def isspace(self):
    #     return self.data.isspace()

    # def istitle(self):
    #     return self.data.istitle()

    def isupper(self):
        return self.data.isupper()

    # def join(self, seq):
    #     return self.data.join(seq)

    def ljust(self, width, *args):
        self.data = self.data.ljust(width, *args)
        return self

    def lower(self):
        self.data = self.data.lower()
        return self

    def lstrip(self, chars=None):
        self.data = self.data.lstrip(chars)
        return self

    maketrans = str.maketrans

    # def partition(self, sep):
    #     return self.data.partition(sep)

    def replace(self, old, new, maxsplit=-1):
        if isinstance(old, MutableMetaString):
            old = old.data
        if isinstance(new, MutableMetaString):
            new = new.data
        self.data = self.data.replace(old, new, maxsplit)
        return self

    def rfind(self, sub, start=0, end=sys.maxsize):
        if isinstance(sub, MutableMetaString):
            sub = sub.data
        return self.data.rfind(sub, start, end)

    def rindex(self, sub, start=0, end=sys.maxsize):
        return self.data.rindex(sub, start, end)

    def rjust(self, width, *args):
        self.data = self.data.rjust(width, *args)
        return self

    # def rpartition(self, sep):
    #     return self.data.rpartition(sep)

    def rstrip(self, chars=None):
        self.data = self.data.rstrip(chars)
        return self

    def split(self, sep=None, maxsplit=-1):
        return self.data.split(sep, maxsplit)

    def rsplit(self, sep=None, maxsplit=-1):
        return self.data.rsplit(sep, maxsplit)

    def splitlines(self, keepends=False):
        return self.data.splitlines(keepends)

    def startswith(self, prefix, start=0, end=sys.maxsize):
        return self.data.startswith(prefix, start, end)

    def strip(self, chars=None):
        self.data = self.data.strip(chars)
        return self

    def swapcase(self):
        self.data = self.data.swapcase()
        return self

    # def title(self):
    #     return self.__class__(self.data.title())

    def strtranslate(self, *args):
        self.data = self.data.translate(*args)
        return self

    def upper(self):
        self.data = self.data.upper()
        return self

    # def zfill(self, width):
    #     return self.__class__(self.data.zfill(width))

def _detect_tool(obj):
    try:
        from Bio.Seq import Seq
        from Bio.SeqRecord import SeqRecord
    except ImportError:
        pass
    else:
        if isinstance(obj, (Seq, SeqRecord)):
            return 'biopython'


class BioSeq(MutableMetaString):

    @property
    def fts(self):
        return self.meta.features

    @fts.setter
    def fts(self, value):
        self.meta.features = FeatureList(value)

    @property
    def gc(self):
        GC = self.count('G') + self.count('C')
        AT = self.count('A') + self.count('T') + self.count('U')
        if GC + AT > 0:
            return GC / (GC + AT)
        else:
            return 0

    @property
    def i(self):
        return _Slicable(self)

    @property
    def rc(self):
        return self.reverse().complement()

    def __getitem__(self, index):
        try:
            subseq = super().__getitem__(index)
        except:
            if isinstance(index, str):
                index = self.meta.features.get(index)
                if index is None:
                    raise TypeError('Feature not found')
            if isinstance(index, Location):
                from sugar.core.fts import _slice_locs
                subseq = _slice_locs(self, [index])
            elif isinstance(index, Feature):
                from sugar.core.fts import _slice_locs
                subseq = _slice_locs(self, index.locs)
                # index = index._slice()
                # if index is None:
                #     msg = f'Feature {index.type} of seq {self.id} has no location'
                #     raise TypeError(msg)
                # subseq = super().__getitem__(index)
            else:
                raise TypeError('Index not supported')
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


    def biotranslate(self, *args, **kw):
        from sugar.core.translate import translate
        import warnings
        msg = 'BioSeq.biotranslate() is deprecated, use translate() method'
        warnings.warn(msg, DeprecationWarning, stacklevel=2)
        self.data = translate(self.data, *args, **kw)
        self.type = 'aa'
        return self

    def complement(self):
        if 'U' in self.data:
            self.replace('U', 'T').strtranslate(COMPLEMENT_TRANS).replace('T', 'U')
        else:
            self.strtranslate(COMPLEMENT_TRANS)
        return self

    def countall(self, **kw):
        return BioBasket([self]).countall(**kw)

    def countplot(self, hue=None, **kw):
        return BioBasket([self]).countplot(hue=hue, **kw)

    def copy(self):
        return copy.deepcopy(self)

    @classmethod
    def fromobj(cls, obj, tool=None):
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

    def match(self, *args, **kwargs):
        """
        Return match object for first found occurence of regex sub, None if not found

        Args:
            sub (str): regex or ``'start'`` or ``'stop'`` to find start/stop codon
            orf (int): May be set to an integer between 0 and 2
                inclusive to respect the corresponding open reading frame.
                Defaults to None to use all 3 ORFs.
            start (int): Index of nucleobase to start matching. Defaults to 0.
            gap (str): Consider gaps of given character, Defaults to None.
            findall (bool): False will return first match, True will
                return all matches. Defaults to False.

        Returns:
            match (match or None)
        """
        return self._match(*args, **kwargs)

    def matchall(self, *args, **kwargs):
        kwargs['matchall'] = True
        return self._match(*args, **kwargs)


    def _match(self, sub, *, orf=None,
               start=0, gap=None, matchall=False):
        from bisect import bisect
        import re

        if isinstance(sub, MutableMetaString):
            sub = sub.data
        if sub == 'start':
            sub = '(AUG|ATG)'
            if gap is not None:
                sub = f'(A[{gap}]*U[{gap}]*G|A[{gap}]*T[{gap}]*G)'
        elif sub == 'stop':
            sub = '(UAG|UAA|UGA|TAG|TAA|TGA)'
            if gap is not None:
                gapstr = f'[{gap}]*'
                sub = ''.join(ch + (gapstr if ch in 'UTAG' and sub[i+1] in 'UATG' else '')
                              for i, ch in enumerate(sub))
        if gap is None:
            gaps = None
        else:
            gaps = [i for i, nt in enumerate(str(self)) if nt == gap if i >= start]
        matches = []
        for m in re.finditer(sub, str(self)):
            if (start is None or (i := m.start()) >= start and (orf is None or (
                    (i - start - bisect(gaps, i)) % 3 == orf
                    if gaps else (i - start) % 3 == orf))):
                     # bisect(gaps, i) gives number of gaps before index i
                if matchall:
                    matches.append(m)
                else:
                    return m
        if matchall:
            return matches

    def tofmtstr(self, *args, **kw):
        return BioBasket([self]).tofmtstr(*args, **kw)

    def toobj(self, tool=None):
        if tool == 'biopython':
            from Bio.Seq import Seq
            from Bio.SeqRecord import SeqRecord
            return SeqRecord(Seq(self.data), id=self.id)
        else:
            raise ValueError(f'Unsupported tool: {tool}')

    def reverse(self):
        self.data = self.data[::-1]
        return self

    def translate(self, *args, **kw):
        from sugar.core.translate import translate
        self.data = translate(self.data, *args, **kw)
        self.type = 'aa'
        return self

    def write(self, *args, **kw):
        BioBasket([self]).write(*args, **kw)


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
        self.meta = Meta(meta)

    def __eq__(self, other):
        if isinstance(other, BioBasket):
            return self.data == other.data and self.meta == other.meta
        return self.data == other

    @property
    def ids(self):
        return [seq.meta.id for seq in self]

    @property
    def rc(self):
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
        if isinstance(i, int):
            seqs = self.data[i]
        elif isinstance(i, slice):
            seqs = self.__class__(self.data[i], meta=self.meta)
        elif isinstance(i, (str, Feature, Location)):
            seqs = self.__class__(self.data, meta=self.meta)
            seqs.data = [seq[i] for seq in seqs.data]
        elif len(i) == 2:
            i, j = i
            if isinstance(i, int):
                seqs = self.data[i][j]
            elif isinstance(i, slice):
                seqs = self.__class__(self.data[i], meta=self.meta)
                seqs.data = [seq[j] for seq in seqs.data]
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
        for seq in self:
            seq.complement()
        return self

    def strtranslate(self, *args, **kw):
        for seq in self:
            seq.strtranslate(*args, **kw)
        return self

    def translate(self, *args, **kw):
        for seq in self:
            seq.translate(*args, **kw)
        return self

    def replace(self, *args, **kw):
        for seq in self:
            seq.replace(*args, **kw)
        return self

    def reverse(self, *args, **kw):
        for seq in self:
            seq.reverse(*args, **kw)
        return self

    def copy(self):
        return copy.deepcopy(self)

    def countall(self, rtype='counter'):
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
        if tool is None and len(obj)>0:
            tool = _detect_tool(obj[0])
        if tool is None:
            raise ValueError('Cannot determine type of object')
        if tool == 'biopython':
            seqs = [BioSeq.fromobj(seq) for seq in obj]
            return cls(seqs)
        else:
            raise ValueError(f'Unsupported tool: {tool}')

    @classmethod
    def fromfmtstr(cls, in_, **kw):
        from sugar import read
        if not isinstance(in_, bytes):
            in_ = in_.encode('latin1')
        return read(io.BytesIO(in_), **kw)

    def todict(self):
        return {seq.id: seq for seq in self}

    @property
    def d(self):
        return self.todict()

    def tofmtstr(self, fmt, **kw):
        out = io.StringIO()
        self.write(out, fmt=fmt, **kw)
        return out.getvalue()

    def tostr(self, h=19, w=80, wid=19, wlen=4, showgc=True, add_hint=False, raw=False):
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
        if tool == 'biopython':
            return [seq.toobj(tool) for seq in self]
        else:
            raise ValueError(f'Unsupported tool: {tool}')

    def write(self, fname, fmt=None, **kw):
        """
        Write sequences to file

        This method calls the underlaying writer routines via `~.main.write()`

        :param fname: filename or file-like object
        :param fmt: format of the file (defaul: auto-detect with file extension)
        :param mode: mode for opening the file, change this only if you know what
            you do
        :param encoding: encoding of the file
        :param tool: use alternative tool for writing the file,
            supported tools are: ``'biopython'``
        :param archive: Explicitely request writing an archive, type may be specified
           (default: auto-detected by file extension)

        All other kwargs are passed to the underlaying writer routine.

        The following formats are supported, for documentation of supported kwargs
        follow the provided links.

        {format_table}
        """
        from sugar._io import write
        write(self, fname, fmt=fmt, **kw)

    def match(self, *args, **kw):
        return [seq.match(*args, **kw) for seq in self]

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




from sugar.core.fts import Location, Defect, Strand
SUGAR = (Attr, BioBasket, BioSeq, Feature, FeatureList, Meta, Location, Defect, Strand)
