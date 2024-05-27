# (C) 2024, Tom Eulenfeld, MIT license
"""
BLAST reader for output generated with option outfmt 7 (preferred), 6, or 11
"""

from sugar.core.fts import Feature, FeatureList, Location
from sugar.core.meta import Meta
from sugar._io.util import _add_fmt_doc


OUTFMT = """
qseqid qgi qacc qaccver qlen
sseqid sallseqid sgi sallgi sacc
saccver sallacc slen qstart qend
sstart send qseq sseq evalue bitscore score
length pident nident mismatch positive
gapopen gaps ppos frames qframe sframe
btop staxid ssciname scomname sblastname
sskingdom staxids sscinames scomnames
sskingdoms stitle salltitles sstrand
qcovs qcovhsp qcovus
""".split()

# Fields:
HEADERFMT = list(map(str.strip, """
query id, query gi, query acc., query acc.ver, query length,
subject id, subject ids, subject gi, subject gis, subject acc.,
subject acc.ver, subject accs., subject length, q. start, q. end,
s. start, s. end, query seq, subject seq, evalue, bit score, score,
alignment length, % identity, identical, mismatches, positives,
gap opens, gaps, % positives, query/sbjct frames, query frame, sbjct frame,
BTOP, subject tax id, subject sci name, subject com names, subject blast name,
subject super kingdom, subject tax ids, subject sci names, subject com names,
subject super kingdoms, subject title, subject titles, subject strand,
% query coverage per subject, % query coverage per hsp, % query coverage per uniq subject
""".split(',')))

TYPES = (
    str, str, str, str, int,
    str, str, str, str, str,
    str, str, int, int, int,
    int, int, str, str, float, float, float,
    int, float, int, int, int,
    int, int, float, str, int, int,
    int, str, str, str, str,
    str, str, str, str,
    str, str, str, str,
    float, float, float)

copyattrs = [('bitscore', 'bitscore'),
             ('evalue', 'evalue'),
             ('sseqid', 'seqid')]

assert len(OUTFMT) == len(HEADERFMT) == len(TYPES)


def extract_seqs(fts, x='source', change_id=False):
    """
    Extract sequences from BLAST features and return `.BioBasket`

    :param x: One of ``'source'`` or ``'query'``
    """
    # TODO remove change_id parameter
    from sugar import BioSeq, BioBasket
    x = x[0]
    seqs = []
    for ft in fts:
        seq = BioSeq(ft.meta._blast[f'{x}seq'],
                     id=ft.meta._blast[f'{x}seqid'],
                     meta={'features': FeatureList([ft])})
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
        ft.meta.seqid = _b.qseqid
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


def is_format_fts(f, *, sep='\t', outfmt=None, **kw):
    content = f.read(1000)
    if content.startswith('#') and 'BLAST' in content:
        return True
    if outfmt is not None:
        line = content.splitlines()[0]
        return len(line.split(sep)) == len(outfmt.split())


def _convert2ind(headers, l=OUTFMT):
    inds = []
    for h in headers:
        inds.append(l.index(h.strip()))
    return inds


@_add_fmt_doc('read_fts')
def read_fts(f, *, sep='\t', outfmt=None, ftype=None):
    """
    BLAST reader for output generated with option outfmt 7 (preferred), 6, or 11

    :param str sep: Separator of fields, use ``','`` for outfmt 11, default ``'\\t'``
    :param str outfmt: The outfmt string passed to BLAST, can be ommited for outfmt 7
    :param str ftype: Parameter used as ftype
    """
    fts = []
    inds = None
    if outfmt is not None:
        inds = _convert2ind(outfmt.split())
    for line in f:
        if outfmt is None and line.startswith('# Fields:'):
            header = line.removeprefix('# Fields:')
            inds = _convert2ind(header.split(','), HEADERFMT)
        elif line.startswith('#') or line.strip() == '':
            pass
        else:
            attrs = {}
            for i, v in enumerate(line.strip().split(sep)):
                attrs[OUTFMT[inds[i]]] = TYPES[inds[i]](v)
            start = attrs['sstart']
            stop = attrs['send']
            if attrs.get('sstrand') == 'N/A':
                strand = '.'
            elif start > stop:
                start, stop, strand = stop, start, '-'
                if attrs.get('sstrand') not in (None, 'minus'):
                    raise ValueError('Expected strand -, got + in sstrand')
            elif start < stop:
                strand = '+'
                if attrs.get('sstrand') not in (None, 'plus'):
                    raise ValueError('Expected strand +, got - in sstrand')
            else:
                strand = attrs.get('sstrand', '.')
            loc = Location(start-1, stop, strand)

            ft = Feature(attrs.get(ftype, ftype), [loc],
                         meta=Meta(_blast=attrs))
            for blastattr, metaattr in copyattrs:
                if blastattr in ft.meta._blast:
                    ft.meta[metaattr] = ft.meta._blast[blastattr]
            fts.append(ft)
    return FeatureList(fts)
