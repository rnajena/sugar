# (C) 2024, Tom Eulenfeld, MIT license
"""
Generic tabular feature reader that is used by several plugins
"""

from collections import namedtuple
from sugar.core.fts import Feature, FeatureList, Location
from sugar.core.meta import Meta
from sugar._io.util import _add_fmt_doc


_BHV = namedtuple('_BlastHeaderValue', 'name long_name type')
_MHV = namedtuple('_MMseqsHeaderValue', 'name type blast_equivalent')
_IHV = namedtuple('_InfernalHeaderValue', 'name type blast_equivalent')
_HEADER = {
    'blast': [
        _BHV('qseqid', 'query id', str),
        _BHV('qgi', 'query gi', str),
        _BHV('qacc', 'query acc.', str),
        _BHV('qaccver', 'query acc.ver', str),
        _BHV('qlen', 'query length', int),
        _BHV('sseqid', 'subject id', str),
        _BHV('sallseqid', 'subject ids', str),
        _BHV('sgi', 'subject gi', str),
        _BHV('sallgi', 'subject gis', str),
        _BHV('sacc', 'subject acc.', str),
        _BHV('saccver', 'subject acc.ver', str),
        _BHV('sallacc', 'subject accs.', str),
        _BHV('slen', 'subject length', int),
        _BHV('qstart', 'q. start', int),
        _BHV('qend', 'q. end', int),
        _BHV('sstart', 's. start', int),
        _BHV('send', 's. end', int),
        _BHV('qseq', 'query seq', str),
        _BHV('sseq', 'subject seq', str),
        _BHV('evalue', 'evalue', float),
        _BHV('bitscore', 'bit score', float),
        _BHV('score', 'score', float),
        _BHV('length', 'alignment length', int),
        _BHV('pident', '% identity', float),
        _BHV('nident', 'identical', int),
        _BHV('mismatch', 'mismatches', int),
        _BHV('positive', 'positives', int),
        _BHV('gapopen', 'gap opens', int),
        _BHV('gaps', 'gaps', int),
        _BHV('ppos', '% positives', float),
        _BHV('frames', 'query/sbjct frames', str),
        _BHV('qframe', 'query frame', int),
        _BHV('sframe', 'sbjct frame', int),
        _BHV('btop', 'BTOP', str),
        _BHV('staxid', 'subject tax id', str),
        _BHV('ssciname', 'subject sci name', str),
        _BHV('scomname', 'subject com name', str),
        _BHV('sblastname', 'subject blast name', str),
        _BHV('sskingdom', 'subject super kingdom', str),
        _BHV('staxids', 'subject tax ids', str),
        _BHV('sscinames', 'subject sci names', str),
        _BHV('scomnames', 'subject com names', str),
        _BHV('sskingdoms', 'subject super kingdoms', str),
        _BHV('stitle', 'subject title', str),
        _BHV('salltitles', 'subject titles', str),
        _BHV('sstrand', 'subject strand', str),
        _BHV('qcovs', '% query coverage per subject', float),
        _BHV('qcovhsp', '% query coverage per hsp', float),
        _BHV('qcovus', '% query coverage per uniq subject', float),
    ],
    'mmseqs': [
        _MHV('query', str, 'qseqid'),
        _MHV('qlen', int, 'qlen'),
        _MHV('target', str, 'sseqid'),
        _MHV('tlen', int, 'slen'),
        _MHV('qstart', int, 'qstart'),
        _MHV('qend', int, 'qend'),
        _MHV('tstart', int, 'sstart'),
        _MHV('tend', int, 'send'),
        _MHV('qseq', str, 'qseq'),
        _MHV('tseq', str, 'sseq'),
        _MHV('evalue', float, 'evalue'),
        _MHV('bits', float, 'bitscore'),
        _MHV('alnlen', int, 'length'),
        _MHV('pident', float, 'pident'),
        _MHV('nident', int, 'nident'),
        _MHV('mismatch', int, 'mismatch'),
        _MHV('gapopen', int, 'gapopen'),
        _MHV('ppos', float, 'ppos'),
        _MHV('qframe', str, 'qframe'),
        _MHV('tframe', str, 'sframe'),
        _MHV('qcov', float, 'qcovs'),
        _MHV('fident', float, None),
        _MHV('raw', str, None),
        _MHV('cigar', str, None),
        _MHV('qheader', str, None),
        _MHV('theader', str, None),
        _MHV('qaln', str, None),
        _MHV('taln', str, None),
        _MHV('tcov', float, None),
        _MHV('qset', str, None),
        _MHV('qsetid', str, None),
        _MHV('tset', str, None),
        _MHV('tsetid', str, None),
        _MHV('taxid', str, None),
        _MHV('taxname', str, None),
        _MHV('taxlineage', str, None),
        _MHV('qorfstart', int, None),
        _MHV('qorfend', int, None),
        _MHV('torfstart', int, None),
        _MHV('torfend', int, None),
    ],
    'infernal': [
        # 18 fields from tblout fmt 1
        _IHV('target', str, 'sseqid'),
        _IHV('target_acc', str, None),
        _IHV('query', str, 'qseqid'),
        _IHV('query_acc', str, None),
        _IHV('model', str, None),
        _IHV('mstart', int, 'qstart'),
        _IHV('mend', int, 'qend'),
        _IHV('sstart', int, 'sstart'),
        _IHV('send', int, 'send'),
        _IHV('sstrand', str, 'sstrand'),
        _IHV('trunc', str, None),
        _IHV('pass', int, None),
        _IHV('gc', float, None),
        _IHV('bias', float, None),
        _IHV('bitscore', float, 'bitscore'),
        _IHV('evalue', float, 'evalue'),
        _IHV('inc', str, None),
        _IHV('description', str, None),
        # 11 additional field from tblout fmt 2
        _IHV('idx', int, None),
        _IHV('clan', str, None),
        _IHV('overlap', str, None),
        _IHV('anyidx', float, None),
        _IHV('anyfrct1', float, None),
        _IHV('anyfrct2', float, None),
        _IHV('winidx', float, None),
        _IHV('winfrct1', float, None),
        _IHV('winfrct2', float, None),
        _IHV('mlen', int, 'qlen'),
        _IHV('slen', int, 'slen'),
    ],
}
_CONVERTH = {
    'blast': {hv.name: hv.name for hv in _HEADER['blast']},
    'mmseqs': {hv.blast_equivalent: hv.name for hv in _HEADER['mmseqs']},
    'infernal': {hv.blast_equivalent: hv.name for hv in _HEADER['infernal']},
}
_DEFAULT_OUTFMT = {
    'blast': (  # outfmt 6
        'qseqid sseqid pident length mismatch gapopen '
        'qstart qend sstart send evalue bitscore'
    ).split(),
    'mmseqs': (  # fmtmode 0
        'query target fident alnlen mismatch gapopen '
        'qstart qend tstart tend evalue bits'
    ).split(),
    'infernal_1': (  # fmt 1
        'target target_acc query query_acc '
        'model mstart mend sstart send sstrand '
        'trunc pass gc bias bitscore evalue inc '
        'description'
    ).split(),
    'infernal_2': (  # fmt 2
        'idx target target_acc query query_acc clan '
        'model mstart mend sstart send sstrand '
        'trunc pass gc bias bitscore evalue inc '
        'overlap anyidx anyfrct1 anyfrct2 winidx winfrct1 winfrct2 '
        'mlen slen description'
    ).split(),
    # fmt 2, old? version found in the Infernal guide
    # and here https://docs.rfam.org/en/latest/genome-annotation.html
    'infernal_2old': (
        'idx target target_acc query query_acc clan '
        'model mstart mend sstart send sstrand '
        'trunc pass gc bias bitscore evalue inc '
        'overlap anyidx anyfrct1 anyfrct2 winidx winfrct1 winfrct2 '
        'description'
    ).split(),
    'infernal_3': (  # fmt 3
        'target target_acc query query_acc '
        'model mstart mend sstart send sstrand '
        'trunc pass gc bias bitscore evalue inc '
        'mlen slen description'
    ).split(),
}
_MMSEQS_HEADER_NAMES = [hv.name for hv in _HEADER['mmseqs']]

copyattrs = [('bitscore', 'score'),
             ('evalue', 'evalue'),
             ('saccver', 'seqid'),
             ('sacc', 'seqid'),
             ('sseqid', 'seqid'),
             ('qaccver', 'name'),
             ('qacc', 'name'),
             ('qseqid', 'name'),
             ]


def extract_seqs(fts, x='source', change_id=False):
    """
    Extract sequences from BLAST features and return `.BioBasket`

    :param x: One of ``'source'`` or ``'query'``
    """
    # TODO remove change_id parameter
    # TODO move these functions to new helper.py module
    # TODO make them more general for all blast-like tools
    from sugar import BioSeq, BioBasket
    x = x[0]
    seqs = []
    for ft in fts:
        seq = BioSeq(ft.meta._blast[f'{x}seq'],
                     id=ft.meta._blast[f'{x}seqid'],
                     meta={'fts': FeatureList([ft])})
        if change_id:
            evalue = ft.meta.evalue
            from math import log10
            if abs(evalue) > 0:
                evalue = log10(evalue)
            start = min(ft.meta._blast[f'{x}start'], ft.meta._blast[f'{x}end'])
            seq.id = f'{seq.id}/{start}{ft.loc.strand}/e{evalue:.0f}'
        seqs.append(seq)
    return BioBasket(seqs)


def get_query_fts(fts):
    """
    Convert BLAST source features to query features
    """
    fts = fts.copy()
    for ft in fts:
        _b = ft.meta._blast
        if 'qseqid' in _b:
            ft.meta.seqid = _b.qseqid
        if 'sseqid' in _b:
            ft.meta.name = _b.sseqid
        start, stop = _b.qstart, _b.qend
        if start > stop:
            start, stop, strand = stop, start, '-'
        elif start < stop:
            strand = '+'
        else:
            strand = '.'
        loc = Location(start-1, stop, strand)
        ft.locs = [loc]
    return fts


def _headers_from_fmtstrings(fmtstrs, fmt='blast', attr='name'):
    """
    Takes a list of format header names and returns a list of header objects
    """
    headers = []
    for fmtstr in fmtstrs:
        for h in _HEADER[fmt]:
            if getattr(h, attr) == fmtstr:
                headers.append(h)
                break
        else:
            raise ValueError(f'Unknown header string: {fmtstr}')
    return headers


def read_tabular(f, *, sep='\t', outfmt=None, ftype=None, fmt='blast',
                  comments=None):
    """
    Reader function used by BLAST, Infernal and MMseqs2 readers
    """
    assert fmt in ('blast', 'infernal', 'mmseqs')
    fts = []
    headers = None
    maxsplit = -1
    if outfmt is not None:
        headers = _headers_from_fmtstrings(outfmt.split(), fmt=fmt)
    for line in f:
        if comments is not None and line.startswith('#'):
            comments.append(line)
        if fmt == 'blast' and outfmt is None and line.startswith('# Fields:'):
            headerline = line.removeprefix('# Fields:')
            headers = _headers_from_fmtstrings(map(str.strip, headerline.split(',')),
                                               fmt=fmt, attr='long_name')
        elif fmt == 'infernal' and headers is None and '--' in line:
            nheaders = len(line.lstrip('#').split())
            fmtv = {18: '1', 29: '2', 20: '3', 27: '2old'}[nheaders]
            headers = _headers_from_fmtstrings(_DEFAULT_OUTFMT[f'{fmt}_{fmtv}'], fmt=fmt)
            maxsplit = nheaders - 1
        elif line.startswith('#') or line.strip() == '':
            pass
        elif (fmt == 'mmseqs' and
                  len(header_fmtstrings := line.strip().split(sep)) > 1 and
                  set(header_fmtstrings) <= set(_MMSEQS_HEADER_NAMES)):
                if headers is None:
                    headers = _headers_from_fmtstrings(header_fmtstrings, fmt=fmt)
        else:
            if headers is None:
                # outfmt not defined in header
                # (as with BLAST outfmt 7 or MMseqs2 fmtmode 4)
                # and not supplied by user
                # assume default
                headers = _headers_from_fmtstrings(_DEFAULT_OUTFMT[fmt], fmt=fmt)
            splitted_line = line.strip().split(sep, maxsplit=maxsplit)
            if len(splitted_line) != len(headers):
                raise ValueError('Number of header fields does not equal '
                                 'number of data values')
            attrs = {}
            for i, v in enumerate(splitted_line):
                try:
                    v = headers[i].type(v)
                except ValueError:
                    pass
                attrs[headers[i].name] = v
            c = _CONVERTH[fmt]
            start = attrs[c['sstart']]
            stop = attrs[c['send']]
            qstart = attrs[c['qstart']]
            qstop = attrs[c['qend']]
            if attrs.get('sstrand') == 'N/A':  # only applicable for BLAST
                if start > stop:
                    # here old BLAST versions were inconsequenct
                    # sstrand is N/A, but start is smaller than stop indicating a ft on the - strand
                    from warnings import warn
                    warn('Inconsequent values: sstrand N/A, but start larger than stop. Strand is set to -.'
                         'Consider updating BLAST, as this was a problem with old BLAST versions. '
                         'If you are already using the latest BLAST version, '
                         'please submit a bug report to sugar developers.')
                    start, stop = stop, start
                    strand = '-'
                else:
                    # in some cases, setting '+' would be more appropriate
                    # but from the available data it is impossible to say
                    strand = '.'
            elif (start > stop and qstart < qstop) or (start < stop and qstart > qstop):
                if start > stop:
                    start, stop = stop, start
                else:
                    qstart, qstop = qstop, qstart
                strand = '-'
                if attrs.get('sstrand') not in (None, '-', 'minus'):
                    raise ValueError('Expected strand -, got + in sstrand')
            elif (start < stop and qstart < qstop) or (start > stop and qstart > qstop):
                if start > stop:
                    start, stop = stop, start
                    qstart, qstop = qstop, qstart
                strand = '+'
                if attrs.get('sstrand') not in (None, '+', 'plus'):
                    raise ValueError('Expected strand +, got - in sstrand')
            else:
                # single residuu features
                strand = attrs.get('sstrand', '.')
            loc = Location(start-1, stop, strand)
            _fmt = '_' + fmt
            ft = Feature(attrs.get(ftype, ftype), [loc],
                         meta=Meta({_fmt: attrs}))
            if 'pident' in ft.meta[_fmt] and 'fident' not in ft.meta[_fmt]:
                ft.meta[_fmt].fident = ft.meta[_fmt].pident / 100
            if 'fident' in ft.meta[_fmt] and 'pident' not in ft.meta[_fmt]:
                ft.meta[_fmt].pident = ft.meta[_fmt].fident * 100
            # TODO adapt
            for blastattr, metaattr in copyattrs:
                if blastattr in c and c[blastattr] in ft.meta[_fmt]:
                    ft.meta[metaattr] = ft.meta[_fmt][c[blastattr]]
            fts.append(ft)
    return FeatureList(fts)
