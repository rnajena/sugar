# (C) 2024, Tom Eulenfeld, MIT license
"""
Several core helper classes and functions, like `~.cane.translate()`, `~.cane.match()` and `.find_orfs()`
"""

import collections
import warnings

from sugar.data import gcode
from sugar.core.fts import FeatureList


class BioMatch(object):
    def __init__(self, match, rf=None, lenseq=None, seqid=None):
        self._match = match
        self.rf = rf
        self.lenseq = lenseq
        self.seqid = seqid
    def __getattr__(self, attr):
        return getattr(self._match, attr)
    def __repr__(self):
        return f'<sugar.BioMatch object; seqid={self.seqid}; rf={self.rf}; span={self.span()}; match={self.group()}>'
    def span(self):
        start, stop = self._match.span()
        if self.rf is not None and self.rf < 0:
            start, stop = self.lenseq - stop, self.lenseq - start
        return start, stop


class BioMatchList(collections.UserList):
    def groupby(self, attr='rf'):
        assert attr in ('rf', 'seqid', 'both')
        d = {}
        if attr == 'both':
            for m in self:
                d.setdefault(m.seqid, {}).setdefault(m.rf, BioMatchList()).append(m)
        else:
            for m in self:
                d.setdefault(getattr(m, attr), BioMatchList()).append(m)
        return d

    @property
    def d(self):
        return self.groupby('seqid')


class ORFList(FeatureList):
    def groupby(self, attr='rf'):
        assert attr in ('rf', 'seqid', 'both')
        d = {}
        if attr == 'both':
            for ft in self:
                d.setdefault(ft.seqid, {}).setdefault(ft.meta.rf, ORFList()).append(ft)
        else:
            for ft in self:
                d.setdefault(getattr(ft.meta, attr), ORFList()).append(ft)
        return d

    def filter(self, minlen=0, rfs=None):
        orfs = []
        for orf in self:
            if len(orf) > minlen and (rfs is None or orf.meta.rf in rfs):
                orfs.append(orf)
        self.data = orfs
        return self

    def sort(self, s=('len',), reverse=False):
        for k in s[::-1]:
            self.data = sorted(self, key=len if k == 'len' else (lambda ft: ft.meta.get(k)), reverse=reverse)
        return self

    @property
    def d(self):
        return self.groupby('seqid')




def _inds2orf(i1, i2, rf, lensec, ftype='ORF', seqid=None):
    from sugar import Feature
    if rf is None or rf >= 0:
        strand = '+'
    else:
        i1, i2 = lensec - i2, lensec - i1
        strand = '-'
    assert i1 < i2
    ft = Feature(ftype, start=i1, stop=i2)
    ft.seqid = seqid
    ft.loc.strand = strand
    ft.meta.rf = rf
    return ft


def find_orfs(seq, rf='fwd', start='start', stop='stop', need_start='always', need_stop=True, gap='-', minlen=0, ftype='ORF'):
    """
    Find open reading frames (ORFS)

    :param seq: `.BioSeq` sequence
    :param rf: reading frame, possible values: int, string or tuple. See also
        `~match`. Default is ``'fwd'``.
    :param start: regular expression defining the start codons, defaults to ATG/AUG
    :param stop: regular expression defining the stop codons, defaults to stop codons
        in default translation table
    :param need_start: One of ``('always', 'once', 'never')``.
        Always: Ech ORF starts with a start codon.
        Once: Only the first ORF in each RF start with a start codon.
        Never: ORFs can start at each codon.
    :param need_stop: Weather the last ORF in each RF needs to end with a stop codon.
    :param gap: gap character inserted into the start and top codon regexes,
        default is ``'-'``.
    :param minlin: Minimum length of ORFs
    :param ftype: Feature type for found ORFS, default is ``'ORF'``

    :returns: Returns a `~ORFList` of all found orfs.
        You can attach these features to sequences using `.BioSeq.add_fts()` or `.BioBasket.add_fts()`.
        Use the `.BioSeq.fts` and `.BioBasket.fts` properties to overwrite all features with the found ORFs.
    """
    # rf  0, 1, 2, -1, -2, -3, 'fwd', 'bwd', 'both'
    assert need_start in ('never', 'always', 'once')
    if need_start in ('once', 'always'):
        starts = seq.matchall(start, rf=rf, gap=gap).groupby()
    stops = seq.matchall(stop, rf=rf, gap=gap).groupby()
    if rf == 'fwd':
        rf = (0, 1, 2)
    elif rf == 'bwd':
        rf = (-1, -2, -3)
    elif rf == 'both':
        rf = (0, 1, 2, -1, -2, -3)
    orfs = []
    for frame in rf:
        i2 = None
        while need_start == 'never' or len(starts[frame]) > 0 or (need_start=='once' and i2 is not None):
            i1 = (frame if need_start == 'never' and i2 is None else
                  i2 if need_start in ('never', 'once') and i2 is not None else
                  starts[frame].pop(0).start())
            if i2 is not None and i1 < i2:  # start codon before last stop codon (alread present in another ORF)
                continue
            while len(stops[frame]) > 0:
                i2 = stops[frame].pop(0).end()
                if i2 > i1:  # stop codon before ORF is ignored
                    break
            else:
                i2 = None if need_stop else len(seq)
            if i2 is not None:
                orf = _inds2orf(i1, i2, frame, len(seq), seqid=seq.id, ftype=ftype)
                if len(orf) >= minlen:
                    orfs.append(orf)
            if i2 in (None, len(seq)):
                break
    return ORFList(orfs)


def match(seq, sub, *, rf='fwd',
          start=0, gap='-', matchall=False):
    """
    Return match object for first found occurence of regex sub, None if not found

    Args:
        sub (str): regex or ``'start'`` or ``'stop'`` to find start/stop codon,
            please specify different codons like
        rf (int): May be set to an integer between -3 and 2
            inclusive to respect the corresponding reading frame.
            Rfs 0 to 2 are on the forward starnd,
            rfs -3 to -1 are on the backward strand,
            You may also specify a set or tuple of reading frames.
            Additionally you can use one of ('fwd', 'bwd', 'both') to select
            all reading frames on the specified strands.
            Defaults to ``'fwd'`` -- all three reading frames on the forward strand.
            You may set rf to ``None`` to ignore reading frames (i.e. for aa seqs)
        start (int): Index of nucleobase to start matching. Defaults to 0.
        gap (str): Consider gaps of given character, Defaults to '-'. The
            character is inserted between each two letters of the regex.
            Be careful, this approach does not work for arbitrary regexes.
        matchall (bool): False will return first match of type `~BioMatch`,
            True will return all matches in a `~BioMatchList`.
            Defaults to False.

    Returns:
        match (match or list of matches or None)
    """
    from bisect import bisect
    import re
    from sugar.core.seq import MutableMetaString

    if isinstance(rf, int):
        rf = (rf,)
    elif isinstance(rf, str):
        assert rf in ('fwd', 'bwd', 'both')
        if rf == 'fwd':
            rf = (0, 1, 2)
        elif rf == 'bwd':
            rf = (-1, -2, -3)
        elif rf == 'both':
            rf = (0, 1, 2, -1, -2, -3)
    if isinstance(sub, MutableMetaString):
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
        gaps = [i for i, nt in enumerate(str(seq)) if nt == gap if i >= start]
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
                    m = BioMatch(m, rf=this_rf, lenseq=len(seq), seqid=seq.id)
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
                    m = BioMatch(m, rf=-1*this_rf-1, lenseq=len(seq), seqid=seq.id)
                    if matchall:
                        matches.append(m)
                    else:
                        return m

    if matchall:
        return BioMatchList(matches)


def translate(seq, *, complete=False, check_start=None, check_stop=False,
              final_stop=None,
              warn=False, astop='X', gap='-', gap_after=2, tt=1):
    """
    Translate a string or `.BioSeq` object into an amino acid string

    :param bool complete: If set to ``True`` ignores stop codons,
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
    :param str gap: gap character, default ``'-'``, set to ``None``
       to raise an error for non nucleotide characters
    :param int gap_after: A single gap in the amino acids string is
        written after the first ``gap_after`` gaps in the
        nucleotide sequence and afterwards after each third gap,
        default is 2
    :param int tt: the number of the translation table, default is 1
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
