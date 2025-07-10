# (C) 2024, Tom Eulenfeld, MIT license
"""
Several core helper classes and functions, such as `~.cane.translate()`, `~.cane.match()`, and `~.cane.find_orfs()`
"""

import collections
import warnings

from sugar.data import gcode
from sugar.core.fts import FeatureList


def _keyfuncs(objs, keys, attr=None):
    from collections.abc import Iterable
    if isinstance(keys, str):
        keys = keys.split()
    if not isinstance(keys, Iterable):
        keys = [keys]
    keyfuncs = [
        (lambda ft, key=key: (getattr(ft, key, None) if attr is None else
                              getattr(getattr(ft, attr), key, None))) if isinstance(key, str) else
        key
        for key in keys
        ]
    return keyfuncs


def _flatten(d, pk=()):
    d2 = {}
    for k, v in d.items():
        k = pk + (k,)
        if isinstance(v, dict):
            d2.update(_flatten(v, k).items())
        else:
            d2[k] = v
    return d2


def _groupby(objs, keys, attr=None, flatten=False):
    """
    Group objects, used by several objects in sugar.core

    :param keys: Tuple of keys or functions to use for grouping.
        May also be a single string or callable.
    :param attr: Attribute where to look for keys
    :param flatten: Return a flattened instead of a nested dict
    :return: Nested dict structure
    """
    if keys is None:
        return {None: objs}
    keyfuncs = _keyfuncs(objs, keys, attr=attr)
    d = {}
    cls = objs.__class__
    for obj in objs:
        d2 = d
        for keyfunc in keyfuncs[:-1]:
            d2 = d2.setdefault(keyfunc(obj), {})
        d2.setdefault(keyfuncs[-1](obj), cls()).append(obj)
    if flatten:
        d = _flatten(d)
    return d


def _sorted(objs, keys=None, reverse=False, attr=None):
    """
    Sort objects, used by several objects in sugar.core

    :param keys: Tuple of keys or functions to use for sorting.
        None can be used as a single value or in the tuple
        to apply the default sorting.
        May also be a single string or callable.
    :param reverse: Use reversed order (default: False)
    :param attr: Attribute where to look for keys

    :return: Sorted objects
    """
    keyfuncs = _keyfuncs(objs, keys, attr=attr)
    for keyfunc in keyfuncs[::-1]:
        objs = sorted(objs, key=keyfunc, reverse=reverse)
    return objs


def _getattr_warn(obj, key, default=None):
    """
    Like getattr with default value, but warns if obj does not have the attribute
    """
    if not hasattr(obj, key):
        warnings.warn(f'Attribute "{key}" not present in {obj.__class__.__name__} object')
    return getattr(obj, key, default)



def _select(objs, attr='meta', **kwargs):
    r"""
    Select objects, used by several objects in sugar.core

    :param attr: Attribute where to look for keys
    :param \*\*kw: All kwargs need to be of the form
        ``key_op=value``, where op is one of
        the operators from the `python:operator` module.
        Additionally, the operators
        ``'in'`` (membership),
        ``'lowerin'`` (membership of lower case word),
        ``'lowereq'`` (equality of lower case word),
        ``'max'`` (alias for le),
        ``'min'`` (alias for ge)
        are supported.
        The different selection criteria are combined with
        the *and* operator.
    :return: Selected objects
    """
    import operator
    ops = {'max': operator.le,
           'min': operator.ge,
           'in': lambda a, b: a in b,
           }
    allowed_funcs = {'len': len}
    for kw, value in kwargs.items():
        if '_' in kw:
            key, kop = kw.rsplit('_', 1)
        else:
            key, kop = kw, 'eq'
        applynot = kop.startswith('not')
        kop = kop.removeprefix('not')
        applylower = kop.startswith('lower')
        kop = kop.removeprefix('lower')
        op = ops.get(kop)
        if op is None:
            op = getattr(operator, kop)
        getv = ((lambda obj: allowed_funcs[key](obj)) if key in allowed_funcs else
                (lambda obj: _getattr_warn(obj, key, None)) if attr is None else
                (lambda obj: _getattr_warn(getattr(obj, attr), key, None)))
        match (applynot, applylower):
            case (False, False): filt = lambda obj: op(getv(obj), value)
            case (True, False): filt = lambda obj: not op(getv(obj), value)
            case (False, True): filt = lambda obj: op(getv(obj).lower(), value)
            case (True, True): filt = lambda obj: not op(getv(obj).lower(), value)
        objs = [obj for obj in objs if filt(obj)]
    return objs


class BioMatch(object):
    """
    The BioMatch object is returned by `~.cane.match()` and the different match methods.

    It is designed to behave like the original `python:re.Match` object.
    See there for available methods.
    It also has the `BioMatch.rf` attribute, which holds the
    reading frame (between -3 and 2, inclusive) of the match.

    .. rubric:: Example:

    >>> match = read()[0].match('AT.')
    >>> match
    <sugar.BioMatch object; seqid=AB047639; rf=2; span=(11, 14); match=ATA>
    >>> print(match.group(), match.rf)
    ATA 2
    """
    def __init__(self, match, rf=None, seqlen=None, seqid=None):
        #: original Match object
        self._match = match
        #: reading frame (-3 to 2)
        self.rf = rf
        self.seqlen = seqlen
        #: sequence id
        self.seqid = seqid
    def __getattr__(self, attr):
        return getattr(self._match, attr)
    def __repr__(self):
        return f'<sugar.BioMatch object; match={self.group()}; span={self.span()}; rf={self.rf}; seqid={self.seqid}>'
    def start(self):
        if self.rf is not None and self.rf < 0:
            return self.seqlen - self._match.end()
        else:
            return self._match.start()
    def end(self):
        if self.rf is not None and self.rf < 0:
            return self.seqlen - self._match.start()
        else:
            return self._match.end()
    def span(self):
        return self.start(), self.end()


class BioMatchList(collections.UserList):
    """
    List of `BioMatch` objects
    """
    def groupby(self, keys='rf'):
        """
        Group matches

        :param keys: Tuple of meta keys or functions to use for grouping.
            Can also be a single string or a callable.
            By default the method groups by rf.
        :return: Nested dict structure
        """
        return _groupby(self, keys)

    def select(self, **kw):
        """
        Select matches

        :param keys: Tuple of meta keys or functions to use for grouping.
            Can also be a single string or a callable.
            By default the method groups by seqid only.
        :return: Nested dict structure
        """
        return _select(self, attr=None, **kw)

    @property
    def d(self):
        """
        Group matches by seqid, alias for ``BioMatchList.groupby('seqid')``
        """
        return self.groupby('seqid')

    def tostr(self):
        lines = [f'{len(self)} BioMatches']
        for m in self:
            lines.append(f'match={m.group()} span={m.span()} rf={m.rf} seqid={m.seqid}')
        return '\n'.join(lines)

    def __str__(self):
        return self.tostr()


def match(seq, sub, *, rf='fwd',
          start=0, gap='-', matchall=False):
    """
    Return `BioMatch` object for first found match of regex sub, None if not found.

    Args:
        sub (str): regex or ``'start'`` or ``'stop'`` to find start/stop codon,
            please specify different codons like
        rf (int): Can be set to an integer between -3 and 2
            inclusive to respect the corresponding reading frame.
            Rfs 0 to 2 are on the forward strand,
            rfs -3 to -1 are on the backward strand,
            You can also specify a set or tuple of reading frames.
            Additionally you can use one of ('fwd', 'bwd', 'all') to select
            all reading frames on the specified strands.
            Defaults to ``'fwd'`` -- all three reading frames on the forward strand.
            You may set rf to ``None`` to ignore reading frames (i.e. for aa seqs)
        start (int): Index of the nucleobase to start matching. Defaults to 0.
        gap (str): Consider gaps of given character, Defaults to '-'. The
            character is inserted between each two letters of the regex.
            Be careful, this approach does not work for arbitrary regexes.
        matchall (bool): False will return first match of type `BioMatch`,
            True will return all matches in a `BioMatchList`.
            Defaults to False.

    Returns:
        match (`BioMatch` or `BioMatchList` of matches or None),
        the list will be sorted by match position,
        matches on the forward strand first,
        then matches on the backward strand.
    """
    from bisect import bisect
    import re
    from sugar.core.seq import BioSeq

    if isinstance(rf, int):
        rf = (rf,)
    elif isinstance(rf, str):
        assert rf in ('fwd', 'bwd', 'all', 'both')
        if rf == 'fwd':
            rf = (0, 1, 2)
        elif rf == 'bwd':
            rf = (-1, -2, -3)
        else:
            rf = (0, 1, 2, -1, -2, -3)
    if isinstance(sub, BioSeq):
        sub = sub.data
    if sub == 'start':
        sub = 'AUG|ATG'
        # if gap is not None:
        #     sub = f'(A[{gap}]*U[{gap}]*G|A[{gap}]*T[{gap}]*G)'
    elif sub == 'stop':
        sub = 'UAG|UAA|UGA|TAG|TAA|TGA'
    if gap is not None:
        gapstr = f'[{gap}]*'
        sub = ''.join(ch + (gapstr if ((ch.isalpha() or ch=='.') and i+1 < len(sub) and
                                       (sub[i+1].isalpha() or sub[i+1] == '.')) else '')
                      for i, ch in enumerate(sub))
    if gap is None or rf is None:
        gaps = None
    else:
        gaps = [i for i, nt in enumerate(str(seq)) if nt in gap if i >= start]
    matches = []
    if rf is not None:
        fwd_rfs = set(rf) & {0, 1, 2}
        bwd_rfs = set(rf) & {-1, -2, -3}
    if rf is None or len(fwd_rfs) > 0:
        for m in re.finditer(sub, str(seq)):
            if (i := m.start()) >= start:
                # bisect(gaps, i) gives number of gaps before index i
                this_rf = (i - start - (bisect(gaps, i) if gaps else 0)) % 3 if rf is not None else None
                if this_rf is None or this_rf in rf:
                    m = BioMatch(m, rf=this_rf, seqlen=len(seq), seqid=seq.id)
                    if matchall:
                        matches.append(m)
                    else:
                        return m
    if rf is not None and len(bwd_rfs) > 0:
        seq = seq.copy().rc()
        for m in re.finditer(sub, str(seq)):
            if (i := m.start()) >= start:
                # bisect(gaps, i) gives number of gaps before index i
                this_rf = (i - start - (bisect(gaps, i) if gaps else 0)) % 3
                if -1*this_rf-1 in rf:
                    m = BioMatch(m, rf=-1*this_rf-1, seqlen=len(seq), seqid=seq.id)
                    if matchall:
                        matches.append(m)
                    else:
                        return m
    if matchall:
        return BioMatchList(matches)


class ORFList(FeatureList):
    """
    List of open reading frames (ORFs)
    """
    pass


def _inds2orf(i1, i2, rf, lensec, ftype='ORF', seqid=None, has_start=None, has_stop=None):
    """
    Create ORF feature from indices
    """
    from sugar import Feature
    if rf is None or rf >= 0:
        strand = '+'
    else:
        i1, i2 = lensec - i2, lensec - i1
        strand = '-'
    assert i1 < i2
    ft = Feature(ftype, start=i1, stop=i2, strand=strand)
    ft.seqid = seqid
    ft.meta.rf = rf
    ft.meta.has_start = has_start
    ft.meta.has_stop = has_stop
    return ft


def _get_from_list(l, i, default=None):
    try:
        return l[i]
    except IndexError:
        return default


def find_orfs(seq, rf='all', start='start', stop='stop', need_start='always', need_stop=True,
              nested='other', gap=None, len_ge=0, ftype='ORF'):
    """
    Find open reading frames (ORFs)

    :param seq: `.BioSeq` sequence
    :param rf: reading frame, possible values: int, string or tuple. See also
        `~match`. Default is ``'all'``.
    :param start: regular expression defining the start codons, defaults to ATG/AUG
    :param stop: regular expression defining the stop codons, defaults to stop codons
        in the default translation table
    :param need_start: One of ``('always', 'once', 'never')``.
        Always: Each ORF starts with a start codon.
        Once: Only the first ORF on the forward and backward strand starts with a start codon.
        Never: ORFs can start at each codon.
    :param need_stop: Whether the last ORF in each RF must end with a stop codon.
    :param nested: Allow nested ORFs fully contained within other ORFs,
        one of ``('no', 'other', 'all')``.
        No: No nested ORFs.
        Other: Nested ORFs allowed in other reading frames (default).
        All: Nested ORFs allowed in all reading frames.
    :param gap: Gap character inserted into the start and stop codon regexes,
        default is None.
    :param len_ge: Return only ORFs with length greater equal, default: 0.
    :param ftype: Feature type for found ORFs, default is ``'ORF'``

    :returns: Returns a `~ORFList` of all found ORFs.
        You can attach these features to sequences using `.BioSeq.add_fts()` or `.BioBasket.add_fts()`.
        Use the `.BioSeq.fts` and `.BioBasket.fts` properties to overwrite features with the found ORFs.

    .. note::
        Python's :func:`python:re.finditer` is used internally to find start and stop codons.
        The limitations of this function apply; for example, matches cannot overlap.
        Care must be taken in special cases.
        For instance, if ORFs do not need to start with a start codon,
        do not use the regular expression ``start='...'``;
        use the ``need_start='never'`` option instead.
    """
    # rf  0, 1, 2, -1, -2, -3, 'fwd', 'bwd', 'all'
    assert need_start in ('never', 'always', 'once')
    assert nested in ('no', 'other', 'all')
    if isinstance(rf, int):
        rf = (rf,)
    if rf == 'fwd':
        rf = (0, 1, 2)
    elif rf == 'bwd':
        rf = (-1, -2, -3)
    elif rf in ('all', 'both'):
        rf = (0, 1, 2, -1, -2, -3)
    rfsort = lambda i: [0, 1, 2, -1, -2, -3].index(i)
    rf = sorted(rf, key=rfsort)
    all_starts = seq.matchall(start, rf=rf, gap=gap)
    all_starts.reverse()
    starts = all_starts.groupby('rf')
    all_stops = seq.matchall(stop, rf=rf, gap=gap)
    all_stops.reverse()
    stops = all_stops.groupby('rf')
    orfs = []
    for frame in rf:
        start_pos = {m._match.start() for m in starts.get(frame, [])}
        framefwd = frame if frame >= 0 else abs(frame) - 1
        i2 = None
        while i2 is None or i2 < len(seq) - 2:
            # find start position of ORF
            if i2 is None and need_start in ('never', 'once'):
                if need_start == 'once':
                    kw = dict(rf=rf, start=start, stop=stop, need_start='always', need_stop=need_stop,
                              nested='other', gap=gap, len_ge=len_ge, ftype='ORF')
                    # TODO this can be enhanced, we only need the first ORF on the fwd and bwd strand.
                    start_orfs = seq.find_orfs(**kw).select(rf_in=((0, 1, 2) if frame >= 0 else (-1, -2, -3)))
                    if len(start_orfs) == 0:
                        break
                    first_start = start_orfs[0]
                    i1 = first_start.loc.start - 1
                    frame0 = first_start.meta.rf
                    if frame0 < 0:
                        i1 = len(seq) - first_start.loc.stop - 1
                        frame0 = abs(frame0) - 1
                else:
                    i1 = -1
                    frame0 = 0
                # find in-frame start position accounting for gaps
                for i in range(i1+1, len(seq)):
                    if gap is None or seq.data[i if frame >= 0 else len(seq) - 1 - i] not in gap:
                        i1 += 1
                        if i1 % 3 == framefwd:
                            break
                else:
                    break
            elif i2 is not None and need_start in ('never', 'once'):
                i1 = i2
            else:
                try:
                    i1 = starts[frame].pop()._match.start()
                except (KeyError, IndexError):  # no new start codon found
                    break
            # find end position of ORF
            if i2 is not None and i1 < i2:  # start codon before last stop codon (nested ORF)
                if nested in ('other', 'no'):
                    continue
            else:
                while stops.get(frame):
                    i2 = stops[frame].pop()._match.end()
                    if i1 < i2:  # stop codons before ORF starts are ignored
                        has_stop = True
                        break
                else:  # did not found any stop codon
                    if need_stop:
                        break
                    has_stop = False
                    if gap is not None and any(g in seq.data for g in gap):
                        warnings.warn('The stop position of the last ORF is set to the end of the sequence for `need_stop=False` in gap mode')
                        i2 = len(seq)
                        # TODO: need to account for gaps
                        # num_chars_from_end = (len(seq) - sum(seq.str.count(g) for g in gap) - framefwd) % 3
                        # while num_chars_from_end > 0:
                        #     i2 -= 1
                        #     if seq.data[i2] not in gap:
                        #         num_chars_from_end -= 1
                    else:
                        i2 = len(seq) - (len(seq) - framefwd) % 3
            orf = _inds2orf(i1, i2, frame, len(seq), seqid=seq.id, ftype=ftype, has_start=i1 in start_pos, has_stop=has_stop)
            if len(orf) >= len_ge:
                orfs.append(orf)
    orfs = ORFList(orfs)
    if nested == 'no':
        orfs.remove_nested()
    return orfs


def translate(seq, *, complete=False, check_start=None, check_stop=False,
              final_stop=None,
              warn=False, astop='X', gap='-', gap_after=2, tt=1):
    """
    Translate a string or `.BioSeq` object into an amino acid string

    :param bool complete: If set to ``True`` ignore stop codons,
        otherwise the translation is stopped at the first stop codon
    :param bool check_start: Check that the first codon is a start codon,
        default is True for ``complete=False`` otherwise False
    :param bool check_stop: Check that the sequence ends with the first stop
        codon, default is False
    :param bool final_stop: Append * for the final stop character,
        defaults to False for ``complete=False`` and True for ``complete=True``
    :param bool warn: Warn if the first codon might not be a start codon,
        warn for ambiguous stop codons,
        warn if the sequence does not end with a stop codon,
        default is False
    :param str astop: Symbol for ambiguous stop codons
    :param str gap: Gap character, default ``'-'``, set to ``None``
       to raise an error for non-nucleotide characters
    :param int gap_after: A single gap in the amino acid string is
        written after the first ``gap_after`` gaps in the
        nucleotide sequence and after every third gap thereafter,
        default is 2
    :param int tt: the number of the translation table, default is 1

    :return: Translated string
    """
    gc = gcode(tt)
    aas = []
    ngap = 0
    check_start = check_start if check_start is not None else not complete
    check_start_warn = warn
    final_stop = final_stop if final_stop is not None else complete
    codon = ''
    for i, nt in enumerate(str(seq).replace('U', 'T')):
        if nt == gap:
            ngap += 1
        else:
            codon = codon + nt
        if gap and gap_after is not None and ngap == gap_after:
            aas.append(gap)
            ngap -= 3
        if len(codon) == 3:
            if check_start_warn or check_start:
                if codon not in gc.starts and codon not in gc.astarts:
                    msg = f'Codon {codon} is not a start codon {gc.starts}'
                    if check_start:
                        raise ValueError(msg)
                    else:
                        warnings.warn(msg)
                        check_start_warn = False
                check_start = False
            if check_start_warn:
                check_start_warn = False
                if codon not in gc.starts:
                    msg = f'Codon {codon} possibly is not a start codon.'
                    warnings.warn(msg)
            try:
                aa = gc.tt[codon]
            except (KeyError):
                aa = 'X'
            if codon in gc.astops:
                aa = astop
                if warn:
                    warnings.warn(f'Codon {codon} might be a stop codon.')
            if codon in gc.stops:
                if (check_stop or warn) and i < len(seq) - 3:
                    msg = 'First stop codon is not at the end of the sequence.'
                    if check_stop:
                        raise ValueError(msg)
                    else:
                        warnings.warn(msg)
                if i >= len(seq) - 3 or not complete:
                    if final_stop:
                        aas.append(aa)
                    break
            aas.append(aa)
            codon = ''
    else:
        if (check_stop or warn) and codon not in gc.astops:
            msg = f'Last codon {codon} is not a stop codon {gc.stops}'
            if check_stop:
                raise ValueError(msg)
            else:
                warnings.warn(msg)
        elif warn and codon in gc.astops:
            msg = f'Last codon {codon} possibly is not a stop codon'
            warnings.warn(msg)
    return ''.join(aas)
