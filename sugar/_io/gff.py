# (C) 2024, Tom Eulenfeld, MIT license
"""
Generic feature format IO
"""


from urllib.parse import quote, unquote
from warnings import warn

from sugar.core.fts import FeatureList, Location, Feature
from sugar.core.meta import Attr, Meta
from sugar._io import read as read_seqs, write as write_seqs
from sugar._io.util import _add_fmt_doc


filename_extensions = ['gff']
filename_extensions_fts = ['gff']


def is_format(f, **kw):
    content = f.read(100)
    if content.strip().startswith('##gff-version 3'):
        return True
    start, stop, _, strand, phase = content.split('\t', maxsplit=9)[3:8]
    int(start)
    int(stop)
    return strand in '+-.?' and phase in '.012'


is_format_fts = is_format


copyattrs = [('Name', 'name'), ('ID', 'id'), ('score', 'score'),
             ('evalue', 'evalue'),
             ('seqid', 'seqid'), ('phase', 'phase')]


@_add_fmt_doc('read_fts')
def read_fts(f, filt=None, default_ftype=None, comments=None):
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
        attrs['seqid'] = seqid
        attrs['source'] = source
        # attrs['type'] = type_
        attrs['score'] = score
        attrs['phase'] = phase

        if 'ID' in attrs:
            id_ = (attrs['ID'], type_, seqid)
        else:
            id_ = None
        if id_ is not None and id_ == lastid:
            for k, v in attrs.items():
                if fts[-1].meta._gff.get(k) != v:
                    if not hasattr(loc, '_gff'):
                        loc._gff = Attr()
                    loc._gff[k] = v
            fts[-1].locs.append(loc)
        else:
            meta = Meta(_gff=attrs)
            fts.append(Feature(type_, locs=[loc], meta=meta))
        lastid = id_

    for ft in fts:
        meta = ft.meta
        if meta._gff.seqid == '.':
            del meta._gff.seqid
        if meta._gff.score != '.':
            meta._gff.score = float(meta._gff.score)
        else:
            del meta._gff.score
        if meta._gff.source == '.':
            del meta._gff.source
        if meta._gff.phase != '.':
            meta._gff.phase = int(meta._gff.phase)
        else:
            del meta._gff.phase
        for gffattr, metaattr in copyattrs:
            if gffattr in meta._gff:
                meta[metaattr] = meta._gff[gffattr]
    return FeatureList(fts)


@_add_fmt_doc('read')
def read(f, **kw):
    """
    Read sequences and their features from GFF file
    """
    fts = read_fts(f, **kw).todict()
    seqs = read_seqs(f, fmt='fasta')
    for seq in seqs:
        if seq.id in fts:
            seq.fts = fts.pop(seq.id, None)
    if len(fts) > 0:
        missing_ids = ', '.join(fts.keys())
        warn(f'Features for seqids {missing_ids} could not be '
             'attached to any sequence')
    return seqs


@_add_fmt_doc('write_fts')
def write_fts(fts, f, header=None):
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
        if hasattr(ft.locs[0], '_gff'):
            meta._gff = {k: v for k, v in meta._gff.items() + ft.locs[0]._gff.items()}
        if len(ft.locs) > 1 and 'ID' not in meta._gff:
            # TODO warn, test
            from random import choices
            from string import ascii_lowercase
            meta._gff.ID = ''.join(choices(ascii_lowercase, k=10))
        seqid = quote(meta._gff.pop('seqid', '.'))
        source = quote(meta._gff.pop('source', '.'))
        score = meta._gff.pop('score', '.')
        phase = meta._gff.pop('phase', '.')
        # meta._gff.pop('type', None)
        type_ = '.' if ft.type is None else ft.type
        for i, loc in enumerate(ft.locs):
            if i == 0:
                gff_meta = meta._gff
                nscore = score
                nphase = phase
            else:
                if hasattr(loc, '_gff'):
                    gff_meta = {k: v for k, v in loc._gff.items() if meta._gff.get(k) != v}
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
def write(seqs, f, **kw):
    """
    Write sequences and their features to GFF file
    """
    for seq in seqs:
        if seq.id is not None and 'features' in seq.meta:
            for ft in seq.fts:
                ft.meta = ft.meta.copy()
                ft.meta.seqid = seq.id
    fts = [ft for seq in seqs for ft in seq.fts]
    write_fts(fts, f, **kw)
    f.write('##FASTA\n')
    write_seqs(seqs, f, fmt='fasta')
