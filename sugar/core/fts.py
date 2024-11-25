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
from enum import Flag, StrEnum, auto
from sugar.core.meta import Meta
from sugar.core.util import _add_inplace_doc


class Defect(Flag):
    """
    This enum type describes location defects.

    A location has a defect, when the feature itself is not directly
    located in the range of the start to the stop base.
    """
    #: No location defect
    NONE         = 0
    #: Part of the feature has been truncated
    #: before the start base
    #: (e.g. by slicing with `FeatureList.slice()`)
    MISS_LEFT    = auto()
    #: Part of the feature has been truncated
    #: after the stop base, inclusive
    #: (e.g. by slicing with `FeatureList.slice()`)
    MISS_RIGHT   = auto()
    #: The feature starts at an unknown position
    #: before the start base
    BEYOND_LEFT  = auto()
    #: The feature stops at an unknown position
    #: after the stop base, inclusive
    BEYOND_RIGHT = auto()
    #: The feature starts at an unknown position
    UNKNOWN_LEFT  = auto()
    #: The feature stops at an unknown position
    UNKNOWN_RIGHT = auto()
    #: The position is between two consecutive bases
    BETWEEN_CONSECUTIVE = auto()
    #: The exact position is unknown, but it is at a
    #: single base between the start and stop residue
    UNKNOWN_SINGLE_BETWEEN = auto()

    def _reverse(self):
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
    This enum type describes the strand of the feature location.
    """
    #: The feature is located on the forward strand
    FORWARD = '+'
    #: The feature is located on the reverse strand
    REVERSE = '-'
    #: The feature is not associated with any strand
    NONE = '.'
    #: The strandness of the feature is unknown
    UNKNOWN = '?'


class Location():

    def __init__(self, start, stop, strand='+', defect=0, meta=None):
        """
        Class describing the contiguous position of a feature
        """
        if start >= stop:
            raise ValueError('start must be lower than stop')
        #:
        self.start = start
        #:
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
        Optionally location may have metadata
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
        return self._strand

    @strand.setter
    def strand(self, v):
        self._strand = Strand(v)

    @property
    def defect(self):
        return self._defect

    @defect.setter
    def defect(self, v):
        self._defect = Defect(v)

    def _reverse(self, seqlen=0):
        loc = self
        start, stop = seqlen-loc.stop, seqlen-loc.start
        strand = loc.strand
        if strand == '+':
            strand = '-'
        elif strand == '-':
            strand = '+'
        defect = loc.defect._reverse()
        return Location(start, stop, strand, defect, meta=loc.meta)


class LocationTuple(tuple):
    """
    Tuple of locations
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
        Get the minimum start base and maximum stop base of all locations
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
        if isinstance(other, LocationTuple):
            lr1 = self.range
            lr2 = other.range
            return (sum(lr1) - sum(lr2)) / 2
        msg = f"'-' not supported between instances of '{type(self).__name__}' and '{type(other).__name__}'"
        raise TypeError(msg)

    def overlaps(self, other):
        """
        Weather the location ranges overlaps with other feature
        """
        if isinstance(other, LocationTuple):
            lr1 = self.range
            lr2 = other.range
            return lr1[0] < lr2[1] and lr1[1] > lr2[0]
        msg = f"LocationTuple.overlaps() not supported for instances of '{type(other).__name__}'"
        raise TypeError(msg)

    def _reverse(self, seqlen=0):
        return LocationTuple([loc._reverse(seqlen=seqlen) for loc in self])


class Feature():
    """
    A single feature/annotation

    :param str type: The name of the feature class, e.g. *gene* or *CDS*
    :param list locs:
        A list of feature locations. In most cases this list will only
        contain one location, but multiple ones are also possible for
        example in virus genomes (due to frame shifts).
    :param start,stop,strand:
        Instead of specifying the locations, a single location can be given
        by start and stop indices and optionally strand.
    :param dict meta:
        The metadata describing the feature.

    .. note::
        The following metadata attributes can be accessed directly as an
        attribute of Feature: *type*, *name*, *id* and *seqid*.
        For example the feature id can be obtained by both `Feature.id`
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
        List of feature locations
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
        """Represent Feature as a string for debugging."""
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
        Weather the location ranges overlaps with other feature
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

        After the operation the feature will be located on the reverse complement strand.

        :param int seqlen: The sequence length, the default 0 will result in negative
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
        A `FeatureList` is a list of features belonging to one or several sequences.

        :param list data: the features
        """
        if hasattr(data, 'data'):
            data = data.data
        super().__init__(data)

    @classmethod
    def frompandas(cls, df, ftype=None):
        if ftype is not None and 'type' not in df:
            if ftype in df:
                df['type'] = df[ftype]
            else:
                df['type'] = ftype
        fts = []
        for rec in df.to_dict('records'):
            loc = Location(rec.pop('start'),
                           rec.pop('stop'),
                           strand=rec.pop('strand', '+'),
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
                    # print(getattr(ft, "type"))
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
        Write features to a string of specified format, see `~.main.write_fts()`
        """
        return self.write(None, fmt)

    def tolist(self, vals='type start stop strand'):
        """
        Return a generator yielding a list for each feature

        :param vals: Parameters from the metadata or location to return,
            might be a string or tuple, defaults to ``'type start stop strand'``

        .. rubric:: Example:

        >>> from sugar import read_fts
        >>> fts = read_fts().select('CDS')
        >>> for type_, start, stop, strand in fts.tolist():
        ...     print(type_, start, strand,)
        CDS 61943 +
        """
        if isinstance(vals, str):
            vals = vals.split()
        for ft in self:
            yield [
                ft.loc.strand if v == 'strand' else
                ft.locs.range[0] if v == 'start' else
                ft.locs.range[1] if v == 'stop' else
                ft.meta.get(v)
                for v in vals
            ]

    def topandas(self, vals='type start stop strand', **kw):
        """
        Return a pandas DataFrame of the features

        :param vals: Parameters from the metadata or location to return,
            might be a string or tuple, defaults to ``'type start stop strand'``

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
        if isinstance(vals, str):
            vals = vals.split()
        kw.setdefault('columns', vals)
        return pandas.DataFrame(self.tolist(vals=vals), **kw)


    def get(self, type):
        """
        Return the first feature of specified feature type, e.g. ``'cds'``

        :param type: String or list of multiple strings
        """
        type_ = type
        if not isinstance(type_, str):
            type_ = tuple(t.lower() for t in type_)
        for ft in self.data:
            if (isinstance(type_, str) and ft.type.lower() == type_.lower() or
                    isinstance(type_, tuple) and ft.type.lower() in type_):
                return ft

    def select(self, type):
        """
        Return features of specified feature type, e.g. ``'cds'``

        For a more powerful selection method, use `filter()`.

        :param type: String or list of multiple strings
        """
        type_ = type
        if not isinstance(type_, str):
            type_ = tuple(t.lower() for t in type_)
        fts = []
        for ft in self.data:
            if (isinstance(type_, str) and ft.type.lower() == type_.lower() or
                    isinstance(type_, tuple) and ft.type.lower() in type_):
                fts.append(ft)
        return self.__class__(fts)

    def todict(self):
        """
        Return a dictionary with sequence ids as keys and FeatureLists as values
        """
        d = {}
        for ft in self:
            seqid = ft.meta.get('seqid', '')
            d.setdefault(seqid, self.__class__()).append(ft)
        return d

    def groupby(self, keys=('seqid',)):
        """
        Group features

        :param keys: Tuple of meta keys or functions to use for grouping.
            May also be a single string or callable.
            By default the method groups by only seqid.
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
        Alias for ``FeatureList.groupby('seqid')``
        """
        return self.groupby('seqid')

    @property
    def loc_range(self):
        """
        Get the range of feature locations,
        i.e. the start and exclusive stop base/residue.

        Returns
        -------
        int : start
            Start location.
        int : stop
            Exclusive stop location.
        """
        start = sys.maxsize
        stop = -sys.maxsize
        for ft in self:
            for loc in ft.locs:
                if loc.start < start:
                    start = loc.start
                if loc.stop > stop:
                    stop = loc.stop
        return start, stop


    def write(self, fname=None, fmt=None, **kw):
        """
        Write features to file, see `~.main.write_fts()`
        """
        from sugar._io import write_fts
        return write_fts(self, fname=fname, fmt=fmt, **kw)


    def slice(self, start, stop, rel=0):
        """
        Return a sub-annotation between start and stop
        """
        if start is None:
            start = -sys.maxsize
        if stop is None:
            stop = sys.maxsize
        # else:
        #     i_stop = stop - 1
        sub_annot = []
        for ft in self:
            locs_in_scope = []
            for loc in ft.locs:
                # Always true for maxsize values
                # in case no start or stop index is given
                if loc.start < stop and loc.stop > start:
                    # The location is at least partly in the
                    # given location range
                    defect = loc.defect
                    if loc.start < start:
                        defect |= Defect.MISS_LEFT
                    if loc.stop > stop:
                        defect |= Defect.MISS_RIGHT
                    lstart = max(start, loc.start) - rel
                    lstop = min(stop, loc.stop) - rel
                    locs_in_scope.append(Location(
                        lstart, lstop, loc.strand, defect, meta=loc.meta
                    ))
            if len(locs_in_scope) > 0:
                # The feature is present in the new annotation
                # if any of the original locations is in the new
                # scope
                new_ft = Feature(locs=locs_in_scope, meta=ft.meta)
                sub_annot.append(new_ft)
        return self.__class__(sub_annot)

    @_add_inplace_doc
    def rc(self, seqlen=0):
        """
        Reverse complement all features, see `Feature.rc()`

        :param int seqlen: The sequence length, the default 0 will result in negative
            location indices.
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
            May also be a single string or callable.
        :param reverse: Use reversed order (default: False)

        :return: Sorted features

        .. rubric:: Example:

        >>> from sugar import read_fts
        >>> fts = read_fts()
        >>> fts.sort(('type', len))  # doctest: +SKIP
        """
        from sugar.core.cane import _sorted
        self.data = _sorted(self.data, keys=keys, reverse=reverse, attr='meta')
        return self

    def filter(self, inplace=False, **kw):
        r"""
        Filter features

        :param \*\*kw: All kwargs need to be of the form
            ``key_op=value``, where op is one of
            the operators from the `python:operator` module.
            Additionally, the operators ``'in'`` (membership),
            ``'max'`` (alias for le)
            ``'min'`` (alias for ge) are supported.
            The different filter conditions are combined with
            the *and* operator. If you need *or*, call filter twice
            and combine the results with ``|`` operator, e.g.
            ``fts.filter(...) | fts.filter(...)``
        :param inplace: Whether to modify the original object (default: False)
        :return: Filtered features

        .. rubric:: Example:

        >>> from sugar import read_fts
        >>> fts = read_fts()
        >>> fts.filter(len_gt=100_000)  # doctest: +SKIP
        """
        from sugar.core.cane import _filter
        filtered = _filter(self.data, **kw)
        if inplace:
            self.data = filtered
            return self
        else:
            return self.__class__(filtered)

    def copy(self):
        """
        Return a deep copy of the object
        """
        return deepcopy(self)
