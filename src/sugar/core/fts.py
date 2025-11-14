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

    @property
    def range(self):
        """
        Get the range of the location or location tuple

        :returns:
            tuple ``start, stop`` with start and stop location
            (zero-based numbering)
        """
        return self.start, self.stop

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

    def __lt__(self, other):
        if isinstance(other, (Location, LocationTuple)):
            start, stop = self.range
            start2, stop2 = other.range
            return start < start2 or (start==start2 and stop < stop2)
        msg = f"'<' not supported between instances of '{type(self).__name__}' and '{type(other).__name__}'"
        raise TypeError(msg)

    def __le__(self, other):
        if isinstance(other, (Location, LocationTuple)):
            start, stop = self.range
            start2, stop2 = other.range
            return start < start2 or (start==start2 and stop <= stop2)
        msg = f"'<=' not supported between instances of '{type(self).__name__}' and '{type(other).__name__}'"
        raise TypeError(msg)

    def __gt__(self, other):
        if isinstance(other, (Location, LocationTuple)):
            start, stop = self.range
            start2, stop2 = other.range
            return start > start2 or (start==start2 and stop > stop2)
        msg = f"'>' not supported between instances of '{type(self).__name__}' and '{type(other).__name__}'"
        raise TypeError(msg)

    def __ge__(self, other):
        if isinstance(other, (Location, LocationTuple)):
            start, stop = self.range
            start2, stop2 = other.range
            return start > start2 or (start==start2 and stop >= stop2)
        msg = f"'>=' not supported between instances of '{type(self).__name__}' and '{type(other).__name__}'"
        raise TypeError(msg)

    @property
    def mid(self):
        """
        Return the middle position of the location or location tuple
        """
        return sum(self.range) // 2

    def shift(self, offset):
        """
        Shift the location by the given offset in-place

        :param int offset: The offset to shift the location
        :return: The shifted location
        :rtype: `.Location`
        """
        self.start += offset
        self.stop  += offset
        return self

    def overlaps(self, other):
        """
        Whether the location/locations overlap with the other location/locations
        """
        if isinstance(other, (Location, LocationTuple)):
            lr1 = self.range
            lr2 = other.range
            return lr1[0] < lr2[1] and lr1[1] > lr2[0]
        msg = f"{type(self).__name__}.overlaps() not supported for instances of '{type(other).__name__}'"
        raise TypeError(msg)

    def overlaplen(self, other):
        """
        Return overlap length with the other location or location tuple
        """
        if isinstance(other, (Location, LocationTuple)):
            lr1 = self.range
            lr2 = other.range
            return max(0, min(lr1[1], lr2[1]) - max(lr1[0], lr2[0]))
        msg = f"{type(self).__name__}.overlaplen() not supported for instances of '{type(other).__name__}'"
        raise TypeError(msg)

    def contains(self, other):
        """
        Whether the location range contains the other location range
        """
        if isinstance(other, (Location, LocationTuple)):
            lr1 = self.range
            lr2 = other.range
            return lr1[0] <= lr2[0] and lr2[1] <= lr1[1]
        msg = f"{type(self).__name__}.overlaps() not supported for instances of '{type(other).__name__}'"
        raise TypeError(msg)

    def distance(self, other, *, pos='inner', sign=False):
        """
        Distance to other location or location tuple

        :param str pos: ``'inner'`` returns the shortest distance between the locations,
            ``'middle'`` returns the distance between the mid locations
        :param bool sign: If set to True, the returned distance will have a negative sign
            if the other location has a smaller position.
            Otherwise the distance will always be larger equal than zero.
        """
        assert pos in ('inner', 'middle')
        if isinstance(other, (Location, LocationTuple)):
            if self.overlaps(other):
                return 0
            if pos == 'middle':
                dist = other.mid - self.mid
            elif self > other:
                dist = other.stop - self.start
            else:
                dist = other.start - self.stop
            if not sign:
                dist = abs(dist)
            return dist
        msg = f"{type(self).__name__}.distance() not supported for instances of '{type(other).__name__}'"
        raise TypeError(msg)

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
    def start(self):
        """
        Get the start position of location tuple
        """
        return min(loc.start for loc in self)

    @property
    def stop(self):
        """
        Get the stop position of location tuple
        """
        return max(loc.stop for loc in self)

    def shift(self, offset):
        """
        Shift the locations by the given offset in-place

        :param int offset: The offset to shift the locations
        :return: The shifted location tuple
        :rtype: `.LocationTuple`
        """
        for loc in self:
            loc.shift(offset)
        return self

    __lt__ = Location.__lt__
    __le__ = Location.__le__
    __gt__ = Location.__gt__
    __ge__ = Location.__ge__
    range = Location.range
    mid = Location.mid
    contains = Location.contains
    distance = Location.distance
    overlaplen = Location.overlaplen
    overlaps = Location.overlaps

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

    def contains(self, other):
        """
        Whether the feature location range contains other
        """
        if isinstance(other, Feature):
            return self.locs.contains(other.locs)
        if isinstance(other, (Location, LocationTuple)):
            return self.locs.contains(other)
        msg = f"Feature.contains() not supported for instances of '{type(other).__name__}'"
        raise TypeError(msg)

    def distance(self, other, **kw):
        """
        Distance to other location or location tuple, see `LocationTuple.distance()`
        """
        if isinstance(other, Feature):
            return self.locs.distance(other.locs, **kw)
        if isinstance(other, (Location, LocationTuple)):
            return self.locs.distance(other, **kw)
        msg = f"Feature.distance() not supported for instances of '{type(other).__name__}'"
        raise TypeError(msg)

    def overlaps(self, other):
        """
        Whether the feature location overlaps with the other
        """
        if isinstance(other, Feature):
            return self.locs.overlaps(other.locs)
        if isinstance(other, (Location, LocationTuple)):
            return self.locs.overlaps(other)
        msg = f"Feature.overlaps() not supported for instances of '{type(other).__name__}'"
        raise TypeError(msg)

    def overlaplen(self, other):
        """
        Return overlap length with the other location or location tuple
        """
        if isinstance(other, Feature):
            return self.locs.overlaplen(other.locs)
        if isinstance(other, (Location, LocationTuple)):
            return self.locs.overlaplen(other)
        msg = f"Feature.overlaplen() not supported for instances of '{type(other).__name__}'"
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

    @classmethod
    def frombiopython(cls, obj):
        """
        Create a `.Feature` object from a biopython_ `~Bio.SeqFeature.SeqFeature` object.

        :param obj: The object to convert.

        Location defects are ignored.
        """
        from sugar.core._adapter import biopython2ft
        return biopython2ft(obj, cls=cls)

    def tobiopython(self):
        """
        Convert Feature to biopython_ `~Bio.SeqFeature.SeqFeature` instance
        """
        from sugar.core._adapter import ft2biopython
        return ft2biopython(self)

    def toftsviewer(self, *, label='default', **kw):
        r"""
        Convert feature to DNAFeaturesViewer_ `~dna_features_viewer.GraphicFeature`

        :param label: The label of the feature,
            may be a str key of the meta dictionary,
            or a function taking the feature and returning the label,
            or the str label itself,
            defaults to ``'name'`` and if that is not present in the metadata, ``'type'``.
        :param \*\*kw: All other kwargs are passed to `~dna_features_viewer.GraphicFeature`.

        Instead of passing label, color and hatch to this function, corresponding values can also be passed via
        the ``Feature.meta`` attribute with the keys ``'_ftsviewer_label'``, ``'_ftsviewer_color'`` and ``'_ftsviewer_hatch'``.
        """
        try:
            from dna_features_viewer import GraphicFeature
        except ImportError as ex:
            raise ImportError('Please install dna_features_viewer to use ftsviewer functionality') from ex
        start, stop = self.locs.range
        strand = {'+': 1, '-': -1, '.': 0, '?': 0}[self.loc.strand]
        if '_ftsviewer_label' in self.meta:
            label = self.meta['_ftsviewer_label']
        elif label == 'default':
            label = self.meta.get('name') or self.type
        elif isinstance(label, str):
            label = self.meta.get(label, label)
        elif label is not None:
            label = str(label(self))
        if '_ftsviewer_color' in self.meta:
            kw['color'] = self.meta['_ftsviewer_color']
        if '_ftsviewer_hatch' in self.meta:
            kw['hatch'] = self.meta['_ftsviewer_hatch']
        return GraphicFeature(
            start=start, end=stop, strand=strand,
            open_left=Defect.MISS_LEFT in self.loc.defect,
            open_right=Defect.MISS_RIGHT in self.locs[-1].defect,
            label=label,
            **kw)


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

    @classmethod
    def frombiopython(cls, obj):
        """
        Create a `.FeatureList` object from a list of biopython_ `~Bio.SeqFeature.SeqFeature` objects.

        :param obj: The object to convert.

        Location defetcs are ignored.
        """
        from sugar.core._adapter import biopython2fts
        return biopython2fts(obj, cls=cls)

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

    def select(self, type=None, *, inplace=False, strand=None, **kw):
        r"""
        Select features

        Two different operating modi can be used, or both.
        Use the ``type`` argument to select features of one type (use a string)
        or of different types (use a list).

        All other kwargs must be of the form
        ``key_op=value``, where op is one of
        the operators from the `operator` module.
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
        selected = self.data
        if type is not None:
            if not isinstance(type, str):
                type = tuple(t.lower() for t in type)
            selected = [
                ft for ft in selected
                if (isinstance(type, str) and ft.type.lower() == type.lower() or
                    isinstance(type, tuple) and ft.type.lower() in type)]
        if strand is not None:
            selected = [ft for ft in selected if ft.loc.strand == strand]
        selected = _select(selected, **kw)
        if inplace:
            self.data = selected
            return self
        else:
            return self.__class__(selected)

    def tobiopython(self):
        """
        Convert the FeatureList to a list of biopython_ `~Bio.SeqFeature.SeqFeature` objects
        """
        from sugar.core._adapter import fts2biopython
        return fts2biopython(self)

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

    def toftsviewer(self, *, label='default', colorby='type', color=None,
                            circular=False,
                            seqlen=None, seq=None,
                            first_index=0,
                            **kw):
        r"""
        Convert features to DNAFeaturesViewer_ `~dna_features_viewer.GraphicRecord`

        :param label: The label of the feature,
            may be a str key of the meta dictionary,
            or a function taking the feature and returning the label,
            or the str label itself,
            defaults to ``'name'`` and if that is not present in the metadata, ``'type'``.
        :param colorby: How to define the color of the features, might be any key in the metadata,
            defaults to ``'type'``, but can also be a function taking a Feature and returning an identifier
        :param color: The color of the features,
            this might be a constant color,
            a list of colors, or
            None for the default matplotlib color cycle (the default), or
            a dictionary mapping the feature identifiers to colors.
        :param circular: If True return an instance of `~dna_features_viewer.CircularGraphicRecord` instead
        :param seq: sequence or sequence data
        :param seqlen: length of sequence, defaults to the length of ``seq`` or the stop location of the last feature.
        :param \*\*kw: All other kwargs are passed to
            `~dna_features_viewer.GraphicFeature` or
            `~dna_features_viewer.GraphicRecord` or
            `~dna_features_viewer.CircularGraphicRecord`, respectively.
        """
        from sugar.core.util import _pop_kws_for_func
        from sugar.imaging.alignment import _get_fts_colordict
        try:
            if circular:
                from dna_features_viewer import CircularGraphicRecord as GR
            else:
                from dna_features_viewer import GraphicRecord as GR
        except ImportError as ex:
            raise ImportError('Please install dna_features_viewer to use ftsviewer functionality') from ex
        kw2 = _pop_kws_for_func(kw, GR)
        color, colorby = _get_fts_colordict(self, color, colorby)
        gfts = [ft.toftsviewer(label=label, color=color[colorby(ft)], **kw) for ft in self]
        if seqlen is None:
            try:
                seqlen = len(seq)
            except TypeError:
                seqlen = self.loc_range[1] - first_index
        return GR(sequence_length=seqlen, sequence=str(seq), features=gfts, first_index=first_index, **kw2)

    def plot_ftsviewer(self, *args, **kw):
        """
        Plot features using DNAFeaturesViewer_, see `~.imaging.ftsviewer.plot_ftsviewer()`
        """
        from sugar.imaging import plot_ftsviewer
        return plot_ftsviewer(self, *args, **kw)

    def groupby(self, keys=('seqid',), flatten=False):
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
        return _groupby(self, keys, attr='meta', flatten=flatten)

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

    def remove_duplicates(self):
        """
        Remove duplicate features
        """
        self.reverse()
        self.data = [ft for i, ft in enumerate(self) if not ft in self[i+1:]]
        self.reverse()
        return self

    def remove_nested(self):
        """
        Remove nested features, i.e. features contained within others
        """
        self.remove_duplicates()
        fts = sorted(self.data, key=len, reverse=True)
        remove = []
        for i, ft in enumerate(fts):
            if ft in remove:
                continue
            for ft2 in fts[i+1:]:
                if ft2 in remove:
                    continue
                if ft.contains(ft2):
                    remove.append(ft2)
        self.data = [ft for ft in self.data if ft not in remove]
        return self

    def remove_overlapping(self):
        """
        Remove overlapping features

        Features on earlier positions in the list are preferred.
        For example, to keep longer features, sort the list beforehand with
        ``fts.sort(len, reverse=True)``.
        """
        self.remove_duplicates()
        remove = []
        for i, ft in enumerate(self):
            if ft in remove:
                continue
            for ft2 in self[i+1:]:
                if ft2 in remove:
                    continue
                if ft.overlaps(ft2):
                    remove.append(ft2)
        self.data = [ft for ft in self.data if ft not in remove]
        return self
