# (C) 2024, Tom Eulenfeld, MIT license

"""
Feature related classes `.Feature`, `.FeatureList`, `.Location`, `.Strand`, `.Defect`
"""

# Originally, this file was based on the annotation module from the biotite package.
# Later, it was rewritten. Constants in Defect amd Strand classes mostly retained their names.

from copy import deepcopy
import collections
import io
import sys
from enum import IntFlag, StrEnum, auto
from sugar.core.meta import Meta
from sugar.core.util import _add_inplace_doc
from sugar.core.util import deprecated


class Defect(IntFlag):
    """
    Types of location defects

    A location has a defect,
    when the feature is not exactly located between start and stop base
    """
    #: No location defect
    NONE = 0
    #: Part of the feature has been truncated
    #: before the start base
    #: (e.g. by slicing with `FeatureList.slice()`)
    MISS_LEFT = auto()
    #: Part of the feature has been truncated
    #: after or at the stop base
    #: (e.g. by slicing with `FeatureList.slice()`)
    MISS_RIGHT = auto()
    #: The feature starts at an unknown position
    #: before the start base
    BEYOND_LEFT = auto()
    #: The feature stops at an unknown position
    #: after or at the stop base
    BEYOND_RIGHT = auto()
    #: The feature starts at an unknown position
    UNKNOWN_LEFT = auto()
    #: The feature stops at an unknown position
    UNKNOWN_RIGHT = auto()
    #: The position is between two consecutive bases
    BETWEEN_CONSECUTIVE = auto()
    #: The exact position is unknown, but it is at a
    #: single base between the start and stop residue
    UNKNOWN_SINGLE_BETWEEN = auto()

    def _reverse(self):
        """
        Return reversed defect

        i.e. describe the same defect from the view of the reverse strand.
        """
        defect = Defect(self)
        if len((self.MISS_LEFT | self.MISS_RIGHT) & self) == 1:
            defect ^= self.MISS_LEFT | self.MISS_RIGHT
        if len((self.BEYOND_LEFT | self.BEYOND_RIGHT) & self) == 1:
            defect ^= self.BEYOND_LEFT | self.BEYOND_RIGHT
        if len((self.UNKNOWN_LEFT | self.UNKNOWN_RIGHT) & self) == 1:
            defect ^= self.UNKNOWN_LEFT | self.UNKNOWN_RIGHT
        return defect


class Strand(StrEnum):
    """
    Types of strand of feature location
    """
    #: The feature is located on the forward strand
    FORWARD = '+'
    #: The feature is located on the reverse strand
    REVERSE = '-'
    #: The feature is not associated with any strand
    NONE = '.'
    #: The strandness of the feature is unknown
    UNKNOWN = '?'

    def _reverse(self):
        """
        Return reversed strand
        """
        return Strand({'+': '-', '-': '+'}.get(self, self))


class Location():
    """
    Class describing the contiguous position of a feature
    """
    def __init__(self, start, stop, strand='+', defect=0, meta=None):
        if start >= stop:
            raise ValueError('start must be lower than stop')
        #: Start location (zero-based numbering)
        self.start = start
        #: Stop location (zero-based numbering)
        self.stop = stop
        self.strand = strand
        self.defect = defect
        self.meta = meta

    def __repr__(self):
        return (f"Location({self.start}, {self.stop}, strand='{self.strand}', "
                f"defect={self.defect.value})")

    def __eq__(self, other):
        if not isinstance(other, Location):
            return False
        return (    self.start  == other.start
                and self.stop   == other.stop
                and self.strand == other.strand
                and self.defect == other.defect
                and self.meta == other.meta
                )

    def __hash__(self):
        return hash((self.start, self.stop, self.strand, self.defect))

    @property
    def meta(self):
        """
        Location can optionally have metadata
        """
        if self._meta is None:
            self._meta = Meta()
        return self._meta

    @meta.setter
    def meta(self, v):
        self._meta = None if v is None else Meta(v)

    @property
    def _stride(self):
        """
        Stride is -1 for the reverse strand, else +1
        """
        return -1 if self.strand == '-' else 1

    def __len__(self):
        return self.stop - self.start

    @property
    def strand(self):
        """Strand of the location"""
        return self._strand

    @strand.setter
    def strand(self, v):
        self._strand = Strand(v)

    @property
    def defect(self):
        """Defect of the location"""
        return self._defect

    @defect.setter
    def defect(self, v):
        self._defect = Defect(v)

    def _reverse(self, seqlen=0):
        """
        Return reversed location

        :param seqlen: Length of the sequence which the location belongs to,
        the default 0 wil return negative start and stop base locations.
        """
        loc = self
        start, stop = seqlen - loc.stop, seqlen - loc.start
        strand = loc.strand._reverse()
        defect = loc.defect._reverse()
        return Location(start, stop, strand, defect, meta=loc.meta)


class LocationTuple(tuple):
    """
    Tuple of contiguous locations, describing the position of a feature
    """
    def __new__(cls, locs=None, start=None, stop=None, strand='+'):
        if start is not None or stop is not None:
            if locs is not None:
                raise ValueError('One of locs or start/stop can be given')
            locs = (Location(start, stop, strand),)
        if locs is None:
            raise ValueError('No location specified')
        if len(locs) == 0:
            raise ValueError('LocationTuple must include at least one location')
        for i, loc in enumerate(locs):
            if not isinstance(loc, Location):
                try:
                    locs[i] = Location(*loc)
                except Exception as ex:
                    msg = 'LocationTuple needs ot be initialized with a tuple of Locations or tuples'
                    raise TypeError(msg) from ex
        locs = tuple(locs)
        if len(locs) > 0:
            strands = set(loc.strand for loc in locs)
            if len(strands) > 1:
                msg = f'Found multiple strand values in Locations: {" ".join(strands)}'
                raise ValueError(msg)
            if locs[0].strand == '-':
                locs = sorted(locs, key=lambda loc: loc.stop, reverse=True)
            else:
                locs = sorted(locs, key=lambda loc: loc.start)
        return super().__new__(cls, locs)

    @property
    def range(self):
        """
        Get the range of locations of this feature

        :returns:
            tuple ``start, stop`` with start and stop location
            (zero-based numbering)
        """
        start = min(loc.start for loc in self)
        stop = max(loc.stop for loc in self)
        return start, stop

    def __lt__(self, other):
        if isinstance(other, LocationTuple):
            start, stop = self.range
            start2, stop2 = other.range
            return start < start2 or (start==start2 and stop < stop2)
        msg = f"'<' not supported between instances of '{type(self).__name__}' and '{type(other).__name__}'"
        raise TypeError(msg)

    def __le__(self, other):
        if isinstance(other, LocationTuple):
            start, stop = self.range
            start2, stop2 = other.range
            return start < start2 or (start==start2 and stop <= stop2)
        msg = f"'<=' not supported between instances of '{type(self).__name__}' and '{type(other).__name__}'"
        raise TypeError(msg)

    def __gt__(self, other):
        if isinstance(other, LocationTuple):
            start, stop = self.range
            start2, stop2 = other.range
            return start > start2 or (start==start2 and stop > stop2)
        msg = f"'>' not supported between instances of '{type(self).__name__}' and '{type(other).__name__}'"
        raise TypeError(msg)

    def __ge__(self, other):
        if isinstance(other, LocationTuple):
            start, stop = self.range
            start2, stop2 = other.range
            return start > start2 or (start==start2 and stop >= stop2)
        msg = f"'>=' not supported between instances of '{type(self).__name__}' and '{type(other).__name__}'"
        raise TypeError(msg)

    def __sub__(self, other):
        """Return position in the middle of both ranges"""
        if isinstance(other, LocationTuple):
            lr1 = self.range
            lr2 = other.range
            return (sum(lr1) - sum(lr2)) // 2
        msg = f"'-' not supported between instances of '{type(self).__name__}' and '{type(other).__name__}'"
        raise TypeError(msg)

    def overlaps(self, other):
        """
        Whether the location ranges overlap with the other location range
        """
        if isinstance(other, LocationTuple):
            lr1 = self.range
            lr2 = other.range
            return lr1[0] < lr2[1] and lr1[1] > lr2[0]
        msg = f"LocationTuple.overlaps() not supported for instances of '{type(other).__name__}'"
        raise TypeError(msg)

    def _reverse(self, seqlen=0):
        """Return reversed LocationTuple"""
        return LocationTuple([loc._reverse(seqlen=seqlen) for loc in self])


class Feature():
    """
    A single feature/annotation

    :param str type: The name of the feature class, e.g. *gene* or *CDS*
    :param list locs:
        A list of feature locations. In most cases this list will
        contain only one location
        but multiple locations are possible,
        for example in virus genomes (due to frame shifts).
    :param start,stop,strand:
        Instead of specifying the locations, a single location can be given
        by start and stop indices and optionally strand.
    :param dict meta:
        The metadata describing the feature.

    .. note::
        The following metadata attributes are directly accessible as
        attributes of Feature: *type*, *name*, *id* and *seqid*.
        For example, the feature id can be obtained by both `Feature.id`
        and ``Feature.meta.id``.
    """

    def __init__(self, type=None, locs=None, meta=None, **kw):
        if meta is None:
            meta = {}
        self.meta = Meta(meta)
        if type is not None:
            self.meta.type = type
        self._locs = LocationTuple(locs=locs, **kw)

    @property
    def locs(self):
        """
        `LocationTuple` of feature locations
        """
        return self._locs

    @locs.setter
    def locs(self, value):
        self._locs = LocationTuple(value)

    @property
    def type(self):
        """
        Alias for ``Feature.meta.type``
        """
        return self.meta.get('type')

    @type.setter
    def type(self, value):
        self.meta.type = value

    @property
    def id(self):
        """
        Alias for ``Feature.meta.id``
        """
        return self.meta.get('id')

    @id.setter
    def id(self, value):
        self.meta.id = value

    @property
    def seqid(self):
        """
        Alias for ``Feature.meta.seqid``
        """
        return self.meta.get('seqid')

    @seqid.setter
    def seqid(self, value):
        self.meta.seqid = value

    @property
    def name(self):
        """
        Alias for ``Feature.meta.name``
        """
        return self.meta.get('name')

    @name.setter
    def name(self, value):
        self.meta.name = value

    def __repr__(self):
        meta = self.meta.copy()
        meta.pop('type', None)
        return f'Feature("{self.type}", [{", ".join([loc.__repr__() for loc in self.locs])}], meta={meta!r})'

    @property
    def loc(self):
        """
        Access first location
        """
        l, *_ = self.locs
        return l

    def __eq__(self, other):
        if not isinstance(other, Feature):
            return False
        return (self.type  == other.type
                and self.locs == other.locs
                and self.meta == other.meta)

    def __lt__(self, other):
        if isinstance(other, Feature):
            if self.seqid != other.seqid:
                return self.seqid < other.seqid
            return self.locs < other.locs
        if isinstance(other, LocationTuple):
            return self.locs < other
        msg = f"'<' not supported between instances of '{type(self).__name__}' and '{type(other).__name__}'"
        raise TypeError(msg)

    def __len__(self):
        lr = self.locs.range
        return lr[1] - lr[0]

    def overlaps(self, other):
        """
        Whether the location ranges overlap with the other feature
        """
        if isinstance(other, Feature):
            return self.locs.overlaps(other.locs)
        if isinstance(other, LocationTuple):
            return self.locs.overlaps(other)
        msg = f"Feature.overlaps() not supported for instances of '{type(other).__name__}'"
        raise TypeError(msg)

    def rc(self, seqlen=0):
        """
        Reverse complement the feature.

        After the in-place operation, the feature will be described by the reverse complement strand.

        :param int seqlen: The sequence length, the default of 0 will result in negative
            location indices.
        """
        self.locs = self.locs._reverse(seqlen=seqlen)
        return self

    def write(self, fname=None, fmt=None, **kw):
        """
        Write feature to file, see `~.main.write_fts()`
        """
        return FeatureList([self]).write(fname=fname, fmt=fmt, **kw)

    # def __hash__(self):
    #     return hash((self.type, self.locs, frozenset(self.meta.items())))


class FeatureList(collections.UserList):
    def __init__(self, data=None):
        """
        A `FeatureList` is a list of features belonging to a single sequence or to different sequences.

        :param list data: the features
        """
        if hasattr(data, 'data'):
            data = data.data
        super().__init__(data)

    @classmethod
    def frompandas(cls, df, ftype=None, one_based=False):
        """
        Convert `pandas.DataFrame` object to `FeatureList`

        :param df: Dataframe with at least start and stop columns.
            The following columns can be used: type, start, stop, len, strand, defect.
            Other columns are stored as metadata.
        :param ftype: If the dataframe has no type column,
            the ``ftype`` column is used instead,
            if it does not exist, ``ftype`` is used directly as type.
        :param one_based: Whether the data uses one-based numbering.
            It will be converted to the zero-based numbering used by sugar.

        :return: created `FeatureList` instance
        """
        if ftype is not None and 'type' not in df:
            df = df.copy()
            if ftype in df:
                df['type'] = df[ftype]
            else:
                df['type'] = ftype
        if 'len' in df:
            df = df.copy()
            if 'start' in df and 'stop' not in df:
                df['stop'] = df['start'] + df['len']
            elif 'start' not in df and 'stop' in df:
                df['start'] = df['stop'] - df['len']
            del df['len']
        fts = []
        for rec in df.to_dict('records'):
            loc = Location(rec.pop('start') - one_based,
                           rec.pop('stop'),
                           strand=rec.pop('strand', '?'),
                           defect=rec.pop('defect', Defect.NONE))
            ft = Feature(locs=[loc], meta=rec)
            fts.append(ft)
        return cls(fts)

    def __str__(self):
        return self.tostr()

    def _repr_pretty_(self, p, cycle):
        if cycle:
            p.text('...')
        else:
            p.text(str(self))

    # Implement all variants of &, |, -, ^
    def __and__(self, other):
        return self.__class__([ft for ft in self if ft in other])

    def __rand__(self, other):
        return self & other

    def __iand__(self, other):
        self.data = [ft for ft in self if ft in other]
        return self

    def __or__(self, other):
        return self + [ft for ft in other if ft not in self]

    def __ror__(self, other):
        return self | other

    def __ior__(self, other):
        self.data += [ft for ft in other if ft not in self]
        return self

    def __sub__(self, other):
        return self.__class__([ft for ft in self if ft not in other])

    def __rsub__(self, other):
        return self.__class__(other) - self

    def __isub__(self, other):
        self.data = [ft for ft in self if ft not in other]
        return self

    def __xor__(self, other):
        return (self | other) - (self & other)

    def __rxor__(self, other):
        return self ^ other

    def __ixor__(self, other):
        self.data = (self ^ other).data
        return self

    def tostr(self, raw=False, w=80, wt=12, wl=20, h=80, exclude_fts=()):
        """
        Return string with information about features, used by ``__str__()`` method
        """
        def _sort_meta_key(m):
            order = ['name', 'gene']
            try:
                return order.index(m[0])
            except ValueError:
                return len(order) + m[0].startswith('_')
        if raw:
            out = []
            for ft in self:
                if str(getattr(ft, 'type', None)) in exclude_fts:
                    continue
                for l in ft.locs:
                    ftstr = (f'{getattr(ft, "type", ".")} {l.start} {l.stop} {l.strand}'
                             f' {ft.meta.get("name", ".")} {ft.meta.get("id", ".")} {ft.meta.get("seqid", ".")}')
                    out.append(ftstr)
            return '\n'.join(out)
        out = []
        wt, wtmax = 0, wt
        wl, wlmax = 0, wl
        wlstart = 0
        wllen = 0
        for ift, ft in enumerate(self):
            if h and ift+1 == h:
                break
            wt = min(max(len(str(ft.type)), wt), wtmax)
            wlstart = max(max(len(f'{l.start:_}') for l in ft.locs), wlstart)
            wllen = max(max(len(f'{len(l):_}') for l in ft.locs), wllen)
        wllen = min(wllen, max(wlmax-wlstart-3, 5))
        for ift, ft in enumerate(self):
            if h and ift+1 == h:
                out.append(f'... and {len(self)-h+1:_} more')
                break
            t = str(getattr(ft, 'type', None))
            if t in exclude_fts:
                continue
            exclude_types = ('translation', 'type')
            metastr = ';'.join(f'{k}={v}' for k, v in
                               sorted(vars(ft.meta).items(), key=_sort_meta_key)
                               if k not in exclude_types)
            for i, l in enumerate(ft.locs):
                locstr = f'{l.start:>{wlstart}_}{l.strand} {len(l):>{wllen}_}'
                ftstr = f'{t if i == 0 else "":>{wt}} {locstr}'
                if i == 0:
                    ftstr = ftstr + f'  {metastr}'
                elif l.meta is not None:
                    locmetastr = ';'.join(f'{k}={v}' for k, v in
                                          sorted(l.meta.items(), key=_sort_meta_key)
                                          if k not in exclude_types)
                    ftstr = ftstr + f'  {locmetastr}'

                if w and len(ftstr) > w:
                    ftstr = ftstr[:w-3] + '...'
                out.append(ftstr)
        return '\n'.join(out)

    def tofmtstr(self, fmt, **kw):
        """
        Write features to a string of the given format, see `~.main.write_fts()`
        """
        return self.write(None, fmt, **kw)

    def tolists(self, keys='type start stop strand'):
        """
        Return a generator yielding a list for each feature

        :param keys: Parameters from the metadata or location to return,
            ``'len'`` is also allowed,
            can be a string or tuple, defaults to ``'type start stop strand'``

        .. rubric:: Example:

        >>> from sugar import read_fts
        >>> fts = read_fts().select('cDNA_match')
        >>> for record in fts.tolists('type start strand len'):
        ...     print(*record)
        cDNA_match 101888622 - 4245
        cDNA_match 103140200 - 30745
        cDNA_match 103944892 - 7136
        cDNA_match 107859806 - 2392
        """
        if isinstance(keys, str):
            keys = keys.split()
        for ft in self:
            yield [
                ft.loc.strand if k == 'strand' else
                ft.loc.defect if k == 'defect' else
                ft.locs.range[0] if k == 'start' else
                ft.locs.range[1] if k == 'stop' else
                len(ft) if k == 'len' else
                ft.meta.get(k)
                for k in keys
            ]

    def topandas(self, keys='type start stop strand', **kw):
        """
        Return a `pandas.DataFrame` of the features

        :param keys: Parameters from the metadata or location to return,
            ``'len'`` is also allowed,
            can be a string or tuple, defaults to ``'type start stop strand'``.

        .. rubric:: Example:

        >>> from sugar import read_fts
        >>> fts = read_fts().select('cDNA_match')
        >>> df = fts.topandas()  # doctest: +SKIP
        >>> print(df)  # doctest: +SKIP
                type      start      stop   strand
        0  cDNA_match  101888622  101892867      -
        1  cDNA_match  103140200  103170945      -
        2  cDNA_match  103944892  103952028      -
        3  cDNA_match  107859806  107862198      -
        """
        import pandas
        if isinstance(keys, str):
            keys = keys.split()
        kw.setdefault('columns', keys)
        return pandas.DataFrame(self.tolists(keys=keys), **kw)

    def get(self, type):
        """
        Return the first feature of given feature type, e.g. ``'cds'``

        :param type: String or list of multiple strings
        """
        type_ = type
        if not isinstance(type_, str):
            type_ = tuple(t.lower() for t in type_)
        for ft in self.data:
            if (isinstance(type_, str) and ft.type.lower() == type_.lower() or
                    isinstance(type_, tuple) and ft.type.lower() in type_):
                return ft

    @deprecated('The filter method is deprecated, use `select()` instead')
    def filter(self, **kw):
        """
        Filter features
        """
        return self.select(**kw)

    def select(self, type=None, *, inplace=False, **kw):
        r"""
        Select features

        Two different operating modi can be used, or both.
        Use the ``type`` argument to select features of one type (use a string)
        or of different types (use a list).

        All other kwargs must be of the form
        ``key_op=value``, where op is one of
        the operators from the `python:operator` module.
        Additionally, the operator ``'in'`` (membership) is supported.
        The different selection criteria are combined with
        the *and* operator. If you need *or*, call select twice
        and combine the results with ``|`` operator, e.g.
        ``fts.select(...) | fts.select(...)``

        :param type: String or list of multiple strings
        :param inplace: Whether to modify the original object (default: False)
        :param \*\*kw: Selection criteria
        :return: Selected features

        .. rubric:: Example:

        >>> from sugar import read_fts
        >>> fts = read_fts()
        >>> fts2 = fts.select('CDS')  # select all CDS fts
        >>> fts3 = fts.select(len_gt=100_000)  # select all fts with length > 100 kB
        """
        from sugar.core.cane import _select
        if type is None:
            selected = self.data
        else:
            if not isinstance(type, str):
                type = tuple(t.lower() for t in type)
            selected = []
            for ft in self.data:
                if (isinstance(type, str) and ft.type.lower() == type.lower() or
                        isinstance(type, tuple) and ft.type.lower() in type):
                    selected.append(ft)
        selected = _select(selected, **kw)
        if inplace:
            self.data = selected
            return self
        else:
            return self.__class__(selected)

    def todict(self):
        """
        Return a dictionary with feature ids as keys and features as values

        .. note::
            This method is different from the `FeatureList.groupby()` method.
            Each value of the dict returned by ``todict()`` is a feature,
            while each value of the dict returned by ``groupby()`` is a
            FeatureList.
        """
        return {ft.id: ft for ft in self}

    def groupby(self, keys=('seqid',)):
        """
        Group features

        :param keys: Tuple of meta keys or functions to use for grouping.
            Can also be a single string or a callable.
            By default, the method groups by seqid only.
        :return: Nested dict structure

        .. rubric:: Example:

        >>> from sugar import read_fts
        >>> fts = read_fts()
        >>> fts.groupby('type')  # doctest: +SKIP
        """
        from sugar.core.cane import _groupby
        return _groupby(self, keys, attr='meta')

    @property
    def d(self):
        """
        Alias for ``FeatureList.todict()``
        """
        return self.todict()

    @property
    def loc_range(self):
        """
        Get the range of locations over all features

        :returns:
            tuple ``start, stop`` with start and stop locations
            (zero-based numbering)
        """
        if len(self) == 0:
            return None
        mins, maxs = zip(*[ft.locs.range for ft in self])
        return min(mins), max(maxs)

    def write(self, fname=None, fmt=None, **kw):
        """
        Write features to file, see `~.main.write_fts()`
        """
        from sugar._io import write_fts
        return write_fts(self, fname=fname, fmt=fmt, **kw)

    def slice(self, start, stop, *, rel=0):
        """
        Return a sub-feature between start and stop

        :param start,stop: start and stop locations
        :param int rel: Subtracts the value ``rel`` from each location position.
        """
        if start is None:
            start = -sys.maxsize
        if stop is None:
            stop = sys.maxsize
        sub_annot = []
        for ft in self:
            sublocs = []
            for loc in ft.locs:
                if loc.start < stop and loc.stop > start:
                    defect = loc.defect
                    if loc.start < start:
                        defect |= Defect.MISS_LEFT
                    if loc.stop > stop:
                        defect |= Defect.MISS_RIGHT
                    lstart = max(start, loc.start) - rel
                    lstop = min(stop, loc.stop) - rel
                    sublocs.append(Location(
                        lstart, lstop, loc.strand, defect, meta=loc.meta
                    ))
            if len(sublocs) > 0:
                new_ft = Feature(locs=sublocs, meta=ft.meta)
                sub_annot.append(new_ft)
        return self.__class__(sub_annot)

    @_add_inplace_doc
    def rc(self, seqlen=0):
        """
        Reverse complement all features, see `Feature.rc()`

        :param int seqlen: The sequence length, the default 0 will result in negative
            location positions.
        """
        for ft in self:
            ft.rc(seqlen=seqlen)
        return self

    @_add_inplace_doc
    def sort(self, keys=None, reverse=False):
        """
        Sort features in-place

        :param keys: Tuple of meta keys or functions to use for sorting.
            None can be used as a single value or in the tuple
            to apply the default sorting by position.
            Can also be a single string or a callable.
        :param reverse: Use reverse order (default: False)

        :return: Sorted features

        .. rubric:: Example:

        >>> from sugar import read_fts
        >>> fts = read_fts()
        >>> fts.sort(('type', len))  # doctest: +SKIP
        """
        from sugar.core.cane import _sorted
        self.data = _sorted(self.data, keys=keys, reverse=reverse, attr='meta')
        return self

    def copy(self):
        """
        Return a deep copy of the object
        """
        return deepcopy(self)
