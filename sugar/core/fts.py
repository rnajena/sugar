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

       - **NONE** - No location defect
       - **MISS_LEFT** - A part of the feature has been truncated
         before the start base/residue of the :class:`Location`
         (probably by indexing an :class:`FeatureList` object)
       - **MISS_RIGHT** - A part of the feature has been truncated
         after the stop base/residue of the :class:`Location`
         (probably by indexing an :class:`FeatureList` object)
       - **BEYOND_LEFT** - The feature starts at an unknown position
         before the start base/residue of the :class:`Location`
       - **BEYOND_RIGHT** - The feature ends at an unknown position
         after the stop base/residue of the :class:`Location`
       - **UNK_LOC** - The exact position is unknown, but it is at a
         single base/residue between the start and stop residue of
         the :class:`Location`, inclusive
       - **BETWEEN** - The position is between to consecutive
         bases/residues.
    """
    NONE         = auto()
    MISS_LEFT    = auto()
    MISS_RIGHT   = auto()
    BEYOND_LEFT  = auto()
    BEYOND_RIGHT = auto()
    UNK_LOC      = auto()
    BETWEEN      = auto()

class Strand(StrEnum):
    """
    This enum type describes the strand of the feature location.
    This is not relevant for protein sequence features.
    """
    FORWARD = '+'
    REVERSE = '-'
    NONE = '.'
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
        self.start = start
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
    def stride(self):
        """
        Stride is -1 for the negative strand, else +1
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


# def _slice_loc(seq, loc):
#     nseq = seq[loc.start:loc.stop]
#     if loc.strand == '-':
#         nseq = nseq.reverse().complement()
#     return nseq


def _slice_locs(seq, locs, splitter=None, filler=None, gap=None):
    # TODO document
    # Concatenate subsequences for each location of the feature
    strand = None
    for loc in locs:
        if loc.strand == strand:
            pass
        elif strand is None:
            strand = loc.strand
        else: # loc.strand != strand
            raise ValueError(
                "All locations of the feature must have the same "
                "strand direction"
            )
    if strand == '-':
        sorted_locs = sorted(locs, key=lambda loc: loc.stop, reverse=True)
    else:
        sorted_locs = sorted(locs, key=lambda loc: loc.start)
    # Merge the sequences corresponding to the ordered locations
    sub_seqs = []
    prev_loc = None
    for i, loc in enumerate(sorted_locs):
        # slice_start = loc.start - seq._seqstart
        # slice_stop = loc.stop - seq._seqstart
        # add_seq = seq[slice_start:slice_stop]
        add_seq = seq.getitem(slice(loc.start, loc.stop), gap=gap)
        if loc.strand == '-':
            add_seq = add_seq.reverse().complement()
        if filler is not None and prev_loc is not None:
            if loc.strand == '-':
                num = prev_loc.start - loc.stop
            else:
                num = loc.start - prev_loc.stop
            if num > 0:
                sub_seqs.append(num * filler)
        sub_seqs.append(add_seq.data)
        prev_loc = loc
    sub_seq = seq[:0]
    if splitter is None:
        splitter = ''
    sub_seq.data = splitter.join(sub_seqs)
    return sub_seq


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
    :param int start,stop:
        Instead of specifying the locations, a single location can be given
        by start and stop indices.
    :param dict meta:
        The metadata describing the feature.

    .. note::
        The following metadata attributes can be accessed directly as an
        attribute of Feature: *type*, *name*, *id* and *seqid*.
        For example the feature id can be obtained by both `Feature.id`
        and `Feature.meta.id`.
    """

    def __init__(self, type=None, locs=None, start=None, stop=None, meta=None):
        if meta is None:
            meta = {}
        self.meta = Meta(meta)
        if type is not None:
            self.meta.type = type
        if start is not None or stop is not None:
            if locs is not None:
                raise ValueError('One of locs or start/stop can be given')
            locs = [Location(start, stop)]
        if len(locs) == 0:
            raise ValueError("A feature must have at least one location")
        self.locs = list(locs)  # frozenset

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
        return f"Feature('{self.type}', [{', '.join([loc.__repr__() for loc in self.locs])}], meta={self.meta!r})"

    @property
    def loc_range(self):
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
        start = min(loc.start for loc in self.locs)
        stop = max(loc.stop for loc in self.locs)
        return start, stop

    def _slice(self):  # needs to return list of slices
        loc = self.loc
        return slice(loc.start, loc.stop, loc.stride)

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
        if not isinstance(other, Feature):
            return False
        if self.seqid != other.seqid:
            return self.seqid < other.seqid
        start, stop = self.loc_range
        start2, stop2 = other.loc_range
        return start < start2 or (start==start2 and stop < stop2)
        # The start base/residue is most significant,
        # if it is emeta for both features, look at stop base/residue
        # if start < it_start:
        #     return True
        # elif start > it_start:
        #     return False
        # else: # start is emeta
        #     return stop < it_stop

    def __gt__(self, item):
        if not isinstance(item, Feature):
            return False
        if self.seqid != item.seqid:
            return self.seqid > item.seqid
        start, stop = self.loc_range
        it_start, it_stop = item.loc_range
        # The start base/residue is most significant,
        # if it is emeta for both features, look at stop base/residue
        if start > it_start:
            return True
        elif start < it_start:
            return False
        else: # start is emeta
            return stop < it_stop

    def __len__(self):
        lr = self.loc_range
        return lr[1] - lr[0]

    def overlaps(self, other):
        """
        Weather the location ranges overlaps with other feature
        """
        if not isinstance(other, Feature):
            raise NotImplementedError()
        lr1 = self.loc_range
        lr2 = other.loc_range
        return lr1[0] < lr2[1] and lr1[1] > lr2[0]

    def __sub__(self, other):
        if not isinstance(other, Feature):
            raise NotImplementedError()
        lr1 = self.loc_range
        lr2 = other.loc_range
        return (sum(lr1) - sum(lr2)) // 2

    def reverse(self):
        """
        Reverse the feature.

        After the operation the feature will be located on the reverse complement strand.
        """
        ft = self
        for loc in ft.locs:
            loc.start, loc.stop = -loc.stop, -loc.start
            if loc.strand == '+':
                loc.strand = '-'
            elif loc.strand == '-':
                loc.strand = '+'
        ft.locs = sorted(ft.locs, key=lambda l: l.start)
        return ft

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
            exclude_types = ('translation', )
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
        Return new `featureList` with all features of specified feature type, e.g. ``'cds'``

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
        >>> fts.groupby('type');
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
        for feature in self:
            for loc in feature.locs:
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


    def slice(self, start, stop):
        """
        Return a sub-annotation between start and stop
        """
        if start is None:
            i_start = -sys.maxsize
        else:
            i_start = start
        if stop is None:
            i_stop = sys.maxsize
        else:
            i_stop = stop - 1
        sub_annot = FeatureList()
        for feature in self:
            locs_in_scope = []
            for loc in feature.locs:
                # Always true for maxsize values
                # in case no start or stop index is given
                if loc.start <= i_stop and loc.stop >= i_start:
                    # The location is at least partly in the
                    # given location range
                    # Handle defects
                    start = loc.start
                    stop = loc.stop
                    defect = loc.defect
                    if loc.start < i_start:
                        defect |= Location.Defect.MISS_LEFT
                        start = i_start
                    if loc.stop > i_stop:
                        defect |= Location.Defect.MISS_RIGHT
                        stop = i_stop
                    locs_in_scope.append(Location(
                        start, stop, loc.strand, defect
                    ))
            if len(locs_in_scope) > 0:
                # The feature is present in the new annotation
                # if any of the original locations is in the new
                # scope
                new_feature = Feature(
                    type=feature.type, locs=locs_in_scope, meta=feature.meta
                )
                sub_annot.append(new_feature)
        return sub_annot

    def reverse(self):
        """
        Reverse all features, see `Feature.reverse()`

        :return: Reversed features
        """
        for ft in self:
            ft.reverse()
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
        >>> fts.sort(('type', len));
        """
        from sugar.core.cane import _sorted
        self.data = _sorted(self.data, keys=keys, reverse=reverse, attr='meta')
        return self

    def filter(self, inplace=True, **kw):
        """
        Filter features

        :param \*\*kw: All kwargs need to be of the form
            ``key_op=value``, where op is one of (is, in, min, max).
            The different filter conditions are combined with
            *and* operator.
        :return: Filtered objects

        .. rubric:: Example:

        >>> from sugar import read_fts
        >>> fts = read_fts()
        >>> fts.filter(len_min=100000)
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





# class AnnotatedSequence():
#     """
#     An :class:`AnnotatedSequence` is a combination of a
#     :class:`Sequence` and an :class:`FeatureList`.

#     Indexing an :class:`AnnotatedSequence` with a slice returns another
#     :class:`AnnotatedSequence` with the corresponding subannotation and
#     a sequence start corrected subsequence, i.e. indexing starts at 1
#     with the default sequence start 1.
#     The sequence start in the newly created :class:`AnnotatedSequence`
#     is the start of the slice.
#     Furthermore, integer indices are allowed in which case the
#     corresponding symbol of the sequence is returned (also sequence
#     start corrected).
#     In both cases the index must be in range of the sequence, e.g. if
#     sequence start is 1, index 0 is not allowed.
#     Negative indices do not mean indexing from the end of the sequence,
#     in contrast to the behavior in :class:`Sequence` objects.
#     Both index types can also be used to modify the sequence.

#     Another option is indexing with a :class:`Feature` (preferably from the
#     :class:`FeatureList` in the same :class:`AnnotatedSequence`).
#     In this case a sequence, described by the location(s) of the
#     :class:`Feature`, is returned.
#     When using a :class:`Feature` for setting an
#     :class:`AnnotatedSequence` with a sequence, the new sequence is
#     replacing the locations of the
#     :class:`Feature`.
#     Note the the replacing sequence must have the same length as the
#     sequence of the :class:`Feature` index.

#     Parameters
#     ----------
#     sequence : Sequence
#         The sequence.
#         Usually a :class:`NucleotideSequence` or
#         :class:`ProteinSequence`.
#     annotation : FeatureList
#         The annotation corresponding to `sequence`.
#     sequence_start : int, optional
#         By default, the start symbol of the sequence is corresponding
#         to location 1 of the features in the annotation. The location
#         of the start symbol can be changed by setting this parameter.
#         Negative values are not supported yet.

#     Attributes
#     ----------
#     sequence : Sequence
#         The represented sequence.
#     annotation : FeatureList
#         The annotation corresponding to `sequence`.
#     sequence_start : int
#         The location of the start symbol in the sequence.

#     See also
#     --------
#     FeatureList, Sequence

#     Examples
#     --------
#     Creating an annotated sequence

#     >>> sequence = NucleotideSequence("ATGGCGTACGATTAGAAAAAAA")
#     >>> feature1 = Feature("misc_feature", [Location(1,2), Location(11,12)],
#     ...                    {"note" : "walker"})
#     >>> feature2 = Feature("misc_feature", [Location(16,22)], {"note" : "poly-A"})
#     >>> annotation = FeatureList([feature1, feature2])
#     >>> annot_seq = AnnotatedSequence(annotation, sequence)
#     >>> print(annot_seq.sequence)
#     ATGGCGTACGATTAGAAAAAAA
#     >>> for f in sorted(list(annot_seq.annotation)):
#     ...     print(f.meta["note"])
#     walker
#     poly-A

#     Indexing with integers, note the sequence start correction

#     >>> print(annot_seq[2])
#     T
#     >>> print(annot_seq.sequence[2])
#     G

#     indexing with slices

#     >>> annot_seq2 = annot_seq[:16]
#     >>> print(annot_seq2.sequence)
#     ATGGCGTACGATTAG
#     >>> for f in annot_seq2.annotation:
#     ...     print(f.meta["note"])
#     walker

#     Indexing with features

#     >>> print(annot_seq[feature1])
#     ATAT
#     >>> print(annot_seq[feature2])
#     AAAAAAA
#     >>> print(annot_seq.sequence)
#     ATGGCGTACGATTAGAAAAAAA
#     >>> annot_seq[feature1] = NucleotideSequence("CCCC")
#     >>> print(annot_seq.sequence)
#     CCGGCGTACGCCTAGAAAAAAA
#     """

#     def __init__(self, annotation, sequence, sequence_start=1):
#         self._annotation = annotation
#         self._sequence = sequence
#         self._seqstart = sequence_start

#     def __repr__(self):
#         """Represent AnnotatedSequence as a string for debugging."""
#         return f'AnnotatedSequence({self._annotation.__repr__()}, {self._sequence.__repr__()}, ' \
#                f'sequence_start={self._seqstart})'

#     @property
#     def sequence_start(self):
#         return self._seqstart

#     @property
#     def sequence(self):
#         return self._sequence

#     @property
#     def annotation(self):
#         return self._annotation

#     def __copy_create__(self):
#         return AnnotatedSequence(
#             self._annotation.copy(), self._sequence.copy, self._seqstart)

#     def reverse_complement(self, sequence_start=1):
#         """
#         Create the reverse complement of the annotated sequence.

#         This method accurately converts the position and the strand of
#         the annotation.
#         The information on the sequence start is lost.

#         Parameters
#         ----------
#         sequence_start : int, optional
#             The location of the start symbol in the reverse complement
#             sequence.

#         Returns
#         -------
#             The reverse complement of the annotated sequence.
#         """
#         rev_seqstart = sequence_start

#         rev_sequence = self._sequence.reverse().complement()

#         seq_len = len(self._sequence)
#         rev_features = []
#         for feature in self._annotation:
#             rev_locs = []
#             for loc in feature.locs:
#                 # Transform location to the reverse complement strand
#                 # (seq_len-1) -> stop sequence index
#                 # (loc.stop-self._seqstart) -> location to index
#                 # ... + rev_seqstart -> index to location
#                 rev_loc_start \
#                     = (seq_len-1) - (loc.stop-self._seqstart) + rev_seqstart
#                 rev_loc_stop \
#                     = (seq_len-1) - (loc.start-self._seqstart) + rev_seqstart

#                 if loc.strand == Location.Strand.FORWARD:
#                     rev_loc_strand = Location.Strand.REVERSE
#                 else:
#                     rev_loc_strand = Location.Strand.FORWARD

#                 rev_loc_defect = Location.Defect.NONE
#                 if loc.defect & Location.Defect.MISS_LEFT:
#                     rev_loc_defect |= Location.Defect.MISS_RIGHT
#                 if loc.defect & Location.Defect.MISS_RIGHT:
#                     rev_loc_defect |= Location.Defect.MISS_LEFT
#                 if loc.defect & Location.Defect.BEYOND_RIGHT:
#                     rev_loc_defect |= Location.Defect.BEYOND_LEFT
#                 if loc.defect & Location.Defect.BEYOND_LEFT:
#                     rev_loc_defect |= Location.Defect.BEYOND_RIGHT
#                 if loc.defect & Location.Defect.UNK_LOC:
#                     rev_loc_defect |= Location.Defect.UNK_LOC
#                 if loc.defect & Location.Defect.BETWEEN:
#                     rev_loc_defect |= Location.Defect.BETWEEN

#                 rev_locs.append(Location(
#                         rev_loc_start, rev_loc_stop,
#                         rev_loc_strand, rev_loc_defect
#                 ))
#             rev_features.append(Feature(
#                 feature.type, rev_locs, feature.meta
#             ))

#         return AnnotatedSequence(
#             FeatureList(rev_features), rev_sequence, rev_seqstart
#         )

#     def __getitem__(self, index):
#         if isinstance(index, Feature):
#             # Concatenate subsequences for each location of the feature
#             locs = index.locs
#             if len(locs) == 0:
#                 raise ValueError("Feature does not contain any locations")
#             # Start by creating an empty sequence
#             sub_seq = self._sequence.copy(new_seq_code=np.array([]))
#             # Locations need to be sorted, as otherwise the locations
#             # chunks would be merged in the wrong order
#             # The order depends on whether the locs are on the forward
#             # or reverse strand
#             strand = None
#             for loc in locs:
#                 if loc.strand == strand:
#                     pass
#                 elif strand is None:
#                     strand = loc.strand
#                 else: # loc.strand != strand
#                     raise ValueError(
#                         "All locations of the feature must have the same "
#                         "strand direction"
#                     )
#             if strand == Location.Strand.FORWARD:
#                 sorted_locs = sorted(
#                     locs, type=lambda loc: loc.start
#                 )
#             else:
#                 sorted_locs = sorted(
#                     locs, type=lambda loc: loc.stop, reverse=True
#                 )
#             # Merge the sequences corresponding to the ordered locations
#             for loc in sorted(locs, type=lambda loc: loc.start):
#                 slice_start = loc.start - self._seqstart
#                 # +1 due to exclusive stop
#                 slice_stop = loc.stop - self._seqstart +1
#                 add_seq = self._sequence[slice_start:slice_stop]
#                 if loc.strand == Location.Strand.REVERSE:
#                     add_seq = add_seq.reverse().complement()
#                 sub_seq += add_seq
#             return sub_seq

#         elif isinstance(index, slice):
#             # Sequence start correction
#             if index.start is None:
#                 seq_start = 0
#             else:
#                 if index.start < self._seqstart:
#                     raise IndexError(
#                         f"The start of the index ({index.start}) is lower "
#                         f"than the start of the sequence ({self._seqstart})"
#                     )
#                 seq_start = index.start - self._seqstart
#             if index.stop is None:
#                 seq_stop = len(self._sequence)
#                 index = slice(index.start, seq_stop, index.step)
#             else:
#                 seq_stop = index.stop - self._seqstart
#             # New value for the sequence start, value is base position
#             if index.start is None:
#                 rel_seq_start = self._seqstart
#             else:
#                 rel_seq_start = index.start
#             return AnnotatedSequence(self._annotation[index],
#                                      self._sequence[seq_start:seq_stop],
#                                      rel_seq_start)

#         elif isinstance(index, numbers.Integral):
#             return self._sequence[index - self._seqstart]

#         else:
#             raise TypeError(
#                 f"'{type(index).__name__}' instances are invalid indices"
#             )

#     def __setitem__(self, index, item):
#         if isinstance(index, Feature):
#             # Item must be sequence
#             # with length emeta to sum of location lengths
#             sub_seq = item
#             sub_seq_i = 0
#             for loc in index.locs:
#                 slice_start = loc.start - self._seqstart
#                 # +1 due to exclusive stop
#                 slice_stop = loc.stop - self._seqstart +1
#                 interval_size = slice_stop - slice_start
#                 self._sequence[slice_start:slice_stop] \
#                     = sub_seq[sub_seq_i : sub_seq_i + interval_size]
#                 sub_seq_i += interval_size
#         elif isinstance(index, slice):
#             # Sequence start correction
#             if index.start is None:
#                 seq_start = 0
#             else:
#                 seq_start = index.start - self._seqstart
#             if index.stop is None:
#                 seq_stop = len(self._sequence)
#             else:
#                 seq_stop = index.stop - self._seqstart
#             # Item is a Sequence
#             self._sequence[seq_start:seq_stop] = item
#         elif isinstance(index, numbers.Integral):
#             # Item is a symbol
#             self._sequence[index - self._seqstart] = item
#         else:
#             raise TypeError(
#                 f"'{type(index).__name__}' instances are invalid indices"
#             )

#     def __eq__(self, item):
#         if not isinstance(item, AnnotatedSequence):
#             return False
#         return (    self.annotation == item.annotation
#                 and self.sequence   == item.sequence
#                 and self._seqstart  == item._seqstart)
