# BSD 3-Clause License
# --------------------

# Copyright 2017, The Biotite contributors
# All rights reserved.

# Redistribution and use in source and binary forms, with or without modification,
# are permitted provided that the following conditions are met:

# 1. Redistributions of source code must retain the above copyright notice, this
# list of conditions and the following disclaimer.

# 2. Redistributions in binary form must reproduce the above copyright notice,
# this list of conditions and the following disclaimer in the documentation and/or
# other materials provided with the distribution.

# 3. Neither the name of the copyright holder nor the names of its contributors
# may be used to endorse or promote products derived from this software without
# specific prior written permission.

# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
# ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
# WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
# DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR
# ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
# (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
# LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON
# ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
# (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
# SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.


# This code is adapted from the biotite package, 2024 Tom Eulenfeld
"""
Feature related classes `.Feature`, `.FeatureList`, `.Location`, `.Location.Strand`, `.Location.Defect`

Some code is modified from bitotite package.
"""

__all__ = ["Location", "Feature", "FeatureList"]

from copy import deepcopy
import collections
import io
import sys
from enum import Flag, StrEnum, auto
from sugar.core.meta import Meta


class Defect(Flag):
    """
    This enum type describes location defects.

    A location has a defect, when the feature itself is not directly
    located in the range of the start to the stop base.
    """
    #: No location defect
    NONE         = auto()
    #: A part of the feature has been truncated
    #: before the start base/residue of the :class:`Location`
    #: (probably by indexing an :class:`FeatureList` object)
    MISS_LEFT    = auto()
    #: A part of the feature has been truncated
    #: after the stop base/residue of the :class:`Location`
    #: (probably by indexing an :class:`FeatureList` object)
    MISS_RIGHT   = auto()
    #: The feature starts at an unknown position
    #: before the start base/residue of the :class:`Location`
    BEYOND_LEFT  = auto()
    #: The feature starts at an unknown position
    #: before the start base/residue of the :class:`Location`
    BEYOND_RIGHT = auto()
    #: The exact position is unknown, but it is at a
    #: single base/residue between the start and stop residue of
    #: the :class:`Location`, inclusive
    UNK_LOC      = auto()
    #: The position is between to consecutive
    #: bases/residues.
    BETWEEN      = auto()

    def reverse(self):
        defect = Defect(self)
        if (self.MISS_LEFT | MISS_RIGHT) & self:
            defect ^= self.MISS_LEFT | MISS_RIGHT
        if (self.BEYOND_LEFT | BEYOND_RIGHT) & self:
            defect ^= self.BEYOND_LEFT | BEYOND_RIGHT
        return defect

class Strand(StrEnum):
    """
    This enum type describes the strand of the feature location.
    """
    #: The feature is located on the forward strand
    FORWARD = '+'
    #: The feature is located on the reverse strand
    REVERSE = '-'
    #: The feature is not associated with any strand, e.g. for proteins
    NONE = '.'
    #: The strandness of the feature is unknown
    UNKNOWN = '?'


class Location():
    Strand = Strand
    Defect = Defect

    def __init__(self, start, stop, strand=Strand.FORWARD,
                 defect=Defect.NONE):
        """
        A :class:`Location` defines at which base(s)/residue(s) a feature is
        located.

        A feature can have multiple :class:`Location` instances if multiple
        locations are joined.

        Objects of this class are immutable.

        Attributes
        ----------
        start : int
            Starting base or residue position of the feature.
        stop : int
            Inclusive ending base or residue position of the feature.
        strand : Strand
            The strand direction.
            Always `Strand.FORWARD` for peptide features.
        defect : Defect
            A possible defect of the location.
        """
        if start >= stop:
            raise ValueError(
                "The start position must be lower than the stop position"
            )
        #:
        self.start = start
        #:
        self.stop = stop
        self.strand = strand
        self.defect = defect

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
                and getattr(self, '_gff', None) == getattr(other, '_gff', None)
                )

    def __hash__(self):
        return hash((self.start, self.stop, self.strand, self.defect))

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
        defect = loc.defect.reverse()
        return Location(start, stop, strand, defect)

    # def _slice(self):
    #     if self._stride == 1:
    #         return slice(self.start, self.stop, self._stride)
    #     else:
    #         start = self.start - 1
    #         if start == -1:
    #             start = None
    #         return slice(self.stop-1, start, self._stride)


# def _slice_loc(seq, loc):
#     nseq = seq[loc.start:loc.stop]
#     if loc.strand == '-':
#         nseq = nseq.reverse().complement()
#     return nseq


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
            raise ValueError('A Locations must include at least one location')
        for loc in locs:
            if not isinstance(loc, Location):
                msg = 'Locations needs ot be initialized with a tuple of Locations'
                raise TypeError(msg)
        strands = set(loc.strand for loc in locs)
        if len(strands) > 1:
            msg = f'Found multiple strand values in Locations: {" ".join(strands)}'
            raise ValueError(msg)
        locs = tuple(locs)
        if locs[0].strand == '-':
            locs = sorted(locs, key=lambda loc: loc.stop, reverse=True)
        else:
            locs = sorted(locs, key=lambda loc: loc.start)
        return super().__new__(cls, locs)

    @property
    def range(self):
        """
        Get the minimum start base/residue and maximum stop base/residue
        of all feature locations.

        This can be used to create a location, that spans all of the
        feature's locations.

        Returns
        -------
        start : int
            The minimum start base/residue of all locations.
        stop : int
            The maximum stop base/residue of all locations.
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
    This class represents a single sequence feature, for example from a
    GenBank feature table.
    A feature describes a functional part of a sequence.
    It consists of a feature type, describing the general class of the
    feature, at least one location, describing its position on the
    reference, and metadata, describing the feature in detail.

    Objects of this class are immutable.

    :param str type: The name of the feature class, e.g. *gene*, *CDS* or
        *regulatory*.
    :param list locs:
        A list of feature locations. In most cases this list will only
        contain one location, but multiple ones are also possible for
        example in eukaryotic CDS (due to splicing) or in virus genomes
        (due to frame shifts).
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

    def reverse(self, seqlen=0):
        """
        Reverse the feature.

        After the operation the feature will be located on the reverse complement strand.
        """
        self.locs = self.locs._reverse(seqlen=seqlen)
        return self

    def write(self, *args, **kw):
        """
        Write feature to file, see `~.main.write_fts()`
        """
        FeatureList([self]).write(self, *args, **kw)

    # def __hash__(self):
    #     return hash((self.type, self.locs, frozenset(self.meta.items())))


class FeatureList(collections.UserList):
    def __init__(self, data=None):
        """
        A `FeatureList` is a set of features belonging to one or several sequences.

        Its advantage over a simple list is the base/residue position based
        indexing:
        When using the `slice()` method call, a subannotation is
        created, containing copies of all :class:`Feature` objects whose
        start and stop base/residue are in range of the slice.
        If the slice starts after the start base/residue or/and the slice
        ends before the stop residue, the position out of range is set to
        the boundaries of the slice (the `Feature` is truncated).
        In this case the :class:`Feature` obtains the
        `Location.Defect.MISS_LEFT` and/or
        `Location.Defect.MISS_RIGHT` defect.
        The third case occurs when a `Feature` starts after the slice
        ends or a `Feature` ends before the slice starts.
        In this case the`Feature` will not appear in the
        subannotation.

        The start or stop position in the slice indices can be omitted, then
        the subannotation will include all features from the start or up to
        the stop, respectively. Step values are ignored.
        The stop values are still exclusive, i.e. the subannotation will
        contain a not truncated :class:`Feature` only if its stop
        base/residue is smaller than the stop value of the slice.

        Multiple :class:`FeatureList` objects can be concatenated to one
        :class:`FeatureList` object using the '+' operator.
        If a feature is present in both :class:`FeatureList` objects, the
        resulting :class:`FeatureList` will contain this feature twice.

        Parameters
        ----------
        data : list
            The features to create the :class:`FeatureList` from. if not
            provided, an empty :class:`FeatureList` is created.

        """
        # TODO
        """
        Examples
        --------
        Creating an annotation from a feature list:

        >>> feature1 = Feature("CDS", [Location(-10, 30 )], meta={"gene" : "test1"})
        >>> feature2 = Feature("CDS", [Location(20,  50 )], meta={"gene" : "test2"})
        >>> annotation = FeatureList([feature1, feature2])
        >>> for f in sorted(list(annotation)):
        ...     print(f.meta["gene"], "".join([str(loc) for loc in f.locs]))
        test1 -10-30 >
        test2 20-50 >

        Merging two annotations and a feature:

        >>> feature3 = Feature("CDS", [Location(100, 130 )], meta={"gene" : "test3"})
        >>> feature4 = Feature("CDS", [Location(150, 250 )], meta={"gene" : "test4"})
        >>> annotation2 = FeatureList([feature3, feature4])
        >>> feature5 = Feature("CDS", [Location(-50, 200 )], meta={"gene" : "test5"})
        >>> annotation = annotation + annotation2 + feature5
        >>> for f in sorted(list(annotation)):
        ...     print(f.meta["gene"], "".join([str(loc) for loc in f.locs]))
        test5 -50-200 >
        test1 -10-30 >
        test2 20-50 >
        test3 100-130 >
        test4 150-250 >

        Location based indexing, note the defects:

        >>> annotation = annotation[40:150]
        >>> for f in sorted(list(annotation)):
        ...     gene = f.meta["gene"]
        ...     loc_str = "".join([f"{loc}    {loc.defect}" for loc in f.locs])
        ...     print(gene, loc_str)
        test5 40-149 >    Defect.MISS_RIGHT|MISS_LEFT
        test2 40-50 >    Defect.MISS_LEFT
        test3 100-130 >    Defect.NONE
        """
        if hasattr(data, 'data'):
            data = data.data
        super().__init__(data)

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
                elif hasattr(l, '_gff'):
                    locmetastr = ';'.join(f'{k}={v}' for k, v in
                                          sorted(l._gff.items(), key=_sort_meta_key)
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
        out = io.StringIO()
        self.write(out, fmt=fmt, **kw)
        return out.getvalue()

    def get(self, type):
        """
        Return the first feature of specified feature type, e.g. ``'cds'``

        :param type: String or list of multiple strings
        """
        type_ = type
        if isinstance(type_, tuple):
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
        if isinstance(type_, tuple):
            type_ = tuple(t.lower() for t in type_)
        fts = []
        for ft in self.data:
            if (isinstance(type_, str) and ft.type.lower() == type_.lower() or
                    isinstance(type_, tuple) and ft.type.lower() in type_):
                fts.append(ft)
        return FeatureList(fts)

    def todict(self):
        """
        Return a dictionary with sequence ids as keys and FeatureLists as values
        """
        d = {}
        for ft in self:
            seqid = ft.meta.get('seqid', '')
            d.setdefault(seqid, FeatureList()).append(ft)
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


    def write(self, fname, fmt=None, **kw):
        """
        Write features to file, see `~.main.write_fts()`
        """
        from sugar._io import write_fts
        write_fts(self, fname, fmt=fmt, **kw)


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
                if loc.start <= stop and loc.stop >= start:
                    # The location is at least partly in the
                    # given location range
                    # Handle defects
                    lstart = loc.start - rel
                    lstop = loc.stop - rel
                    defect = loc.defect
                    if loc.start < start:
                        defect |= Location.Defect.MISS_LEFT
                        lstart = start - rel
                    if loc.stop > stop:
                        defect |= Location.Defect.MISS_RIGHT
                        lstop = stop - rel
                    locs_in_scope.append(Location(
                        lstart, lstop, loc.strand, defect
                    ))
            if len(locs_in_scope) > 0:
                # The feature is present in the new annotation
                # if any of the original locations is in the new
                # scope
                new_ft = Feature(locs=locs_in_scope, meta=ft.meta)
                sub_annot.append(new_ft)
        return FeatureList(sub_annot)

    def reverse(self, seqlen=0):
        """
        Reverse all features, see `Feature.reverse()`

        :return: Reversed features
        """
        for ft in self:
            ft.reverse(seqlen=seqlen)
        return self

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

    def filter(self, inplace=True, **kw):
        """
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
        :param inplace: Whether to modify the original object.
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
