# (C) 2024, Tom Eulenfeld, MIT license
"""
Generic feature format (`GFF`_) IO
"""


from urllib.parse import quote, unquote
from warnings import warn

from sugar.core.fts import FeatureList, Location, Feature
from sugar.core.meta import Attr, Meta
from sugar._io import read as read_seqs, write as write_seqs
from sugar._io.util import _add_fmt_doc


filename_extensions = ['gff']
filename_extensions_fts = ['gff']


def is_gff(f, **kw):
    content = f.read(100)
    if content.strip().startswith('##gff-version 3'):
        return True
    start, stop, _, strand, phase = content.split('\t', maxsplit=9)[3:8]
    int(start)
    int(stop)
    return strand in '+-.?' and phase in '.012'


is_fts_gff = is_gff


copyattrs = [('Name', 'name'), ('ID', 'id'), ('score', 'score'),
             ('evalue', 'evalue'),
             ('seqid', 'seqid'), ('phase', 'phase'),
             ('type', 'type')]


@_add_fmt_doc('read_fts')
def read_fts_gff(f, filt=None, default_ftype=None, comments=None):
    """
    Read a GFF file and return `.FeatureList`

    :param list filt: Return only Features of given ftypes, default: all
    :param str default_ftype: default ftype for entries without type
    :param list comments: comment lines inside the file are stored in
        the comments list (optional)
    """
    fts = []
    lastid = None
    for line in f:
        if line.startswith('##FASTA'):
            break
        elif line.startswith('#') or line.strip() == '':
            if comments is not None:
                comments.append(line)
            continue
        *cols, attrcol = line.strip().split('\t')
        seqid, source, type_, start, stop, score, strand, phase = map(unquote, cols)
        if type_ == '.':
            type_ = default_ftype
        if filt and type_ not in filt:
            continue
        loc = Location(int(start)-1, int(stop), strand=strand)
        attrs = {}
        if attrcol != '.':
            for kv in attrcol.split(';'):
               k, v = kv.strip().split('=')
               attrs[unquote(k.strip())] = (
                   unquote(v.strip()) if ',' not in v else
                   [unquote(vv.strip()) for vv in v.strip().split(',')])
        if seqid != '.':
            attrs['seqid'] = seqid
        if source != '.':
            attrs['source'] = source
        if score != '.':
            attrs['score'] = float(score)
        if phase != '.':
            attrs['phase'] = int(phase)
        if 'ID' in attrs:
            id_ = (attrs['ID'], type_, seqid)
        else:
            id_ = None
        if id_ is not None and id_ == lastid:
            for k, v in attrs.items():
                if fts[-1].meta._gff.get(k) != v:
                    # TODO: add proper meta attribute to Location?
                    if not hasattr(loc.meta, '_gff'):
                        loc.meta._gff = Attr()
                    loc.meta._gff[k] = v
            fts[-1].locs = fts[-1].locs + (loc,)
        else:
            meta = Meta(_gff=attrs)
            fts.append(Feature(type_, locs=[loc], meta=meta))
        lastid = id_
    for ft in fts:
        for gffattr, metaattr in copyattrs:
            if gffattr in ft.meta._gff:
                ft.meta[metaattr] = ft.meta._gff[gffattr]
    return FeatureList(fts)


@_add_fmt_doc('read')
def read_gff(f, **kw):
    """
    Read sequences and their features from GFF file
    """
    fts = read_fts_gff(f, **kw)
    seqs = read_seqs(f, fmt='fasta')
    seqs.fts = fts
    return seqs


@_add_fmt_doc('write_fts')
def write_fts_gff(fts, f, header=None):
    """
    Write features to GFF file

    :param str header: Optionally write a header at the top of file
    """
    f.write('##gff-version 3\n')
    if header:
        f.write(header)
    for ft in fts:
        meta = ft.meta.copy()
        meta.setdefault('_gff', {})
        for gffattr, metaattr in copyattrs:
            if metaattr in meta:
                meta._gff[gffattr] = meta[metaattr]
        if ft.locs[0]._meta and hasattr(ft.locs[0].meta, '_gff'):
            meta._gff = {k: v for k, v in meta._gff.items() + ft.locs[0].meta._gff.items()}
        if len(ft.locs) > 1 and 'ID' not in meta._gff:
            # TODO warn, test
            from random import choices
            from string import ascii_lowercase
            meta._gff.ID = ''.join(choices(ascii_lowercase, k=10))
        seqid = quote(meta._gff.pop('seqid', '.'))
        source = quote(meta._gff.pop('source', '.'))
        score = meta._gff.pop('score', '.')
        phase = meta._gff.pop('phase', '.')
        type_ = meta._gff.pop('type', None) or meta.get('type') or '.'
        for i, loc in enumerate(ft.locs):
            if i == 0:
                gff_meta = meta._gff
                nscore = score
                nphase = phase
            else:
                if loc._meta and hasattr(loc.meta, '_gff'):
                    gff_meta = {k: v for k, v in loc.meta._gff.items() if meta._gff.get(k) != v}
                else:
                    gff_meta = {}
                gff_meta['ID'] = meta._gff.ID
                nscore = gff_meta.pop('score', score)
                nphase = gff_meta.pop('phase', phase)
            if len(gff_meta) == 0:
                attrstr = '.'
            else:
                attrstr = ';'.join(quote(k) + '=' + (','.join(quote(vv) for vv in v) if isinstance(v, (list, tuple)) else
                                                     quote(str(v)))
                                   for k, v in gff_meta.items())
            f.write(f'{seqid}\t{source}\t{type_}\t{loc.start+1}\t{loc.stop}\t{nscore}\t{loc.strand}\t{nphase}\t{attrstr}\n')


@_add_fmt_doc('write')
def write_gff(seqs, f, **kw):
    """
    Write sequences and their features to GFF file
    """
    write_fts_gff(seqs.fts, f, **kw)
    f.write('##FASTA\n')
    write_seqs(seqs, f, fmt='fasta')
