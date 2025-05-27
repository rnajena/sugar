# (C) 2024, Tom Eulenfeld, MIT license
"""
Generic feature format (`GFF`_) and Gene transfer format (`GTF`_) IO

For a review of the different versions of GFF and GTF, see `here`_.

The readers in this module also support older versions of the GFF format.

.. _here: https://agat.readthedocs.io/en/latest/gxf.html
"""

from urllib.parse import unquote
from warnings import warn

from sugar.core.fts import FeatureList, Location, Feature
from sugar.core.meta import Attr, Meta
from sugar._io import read as read_seqs, write as write_seqs
from sugar._io.util import _add_fmt_doc


filename_extensions_gff = ['gff']
filename_extensions_fts_gff = ['gff']
filename_extensions_fts_gtf = ['gtf']


def is_gff(f, **kw):
    content = f.read(100)
    return content.strip().startswith('##gff-version')


def _detect_gff_version(f):
    fpos = f.tell()
    try:
        content = f.read(100).strip()
        if content.startswith('##gff-version'):
            version = content.split(maxsplit=2)[1]
            if version in '123':
                return version
    except Exception:
        pass
    finally:
        f.seek(fpos)


is_fts_gff = is_gff


def is_fts_gtf(f, **kw):
    content = f.read(500)
    if not content.strip().startswith('##gff-version'):
        for line in content.splitlines():
            if line.startswith('#'):
                continue
            start, stop, _, strand, phase = line.split('\t', maxsplit=9)[3:8]
            int(start)
            int(stop)
            return strand in '+-.' and phase in '.012'


copyattrs = {
    'gff' : [
        ('Name', 'name'), ('ID', 'id'), ('E_value', 'evalue'), ('evalue', 'evalue'),
        ('seqid', 'seqid'), ('type', 'type'), ('score', 'score'),
    ],
    'gtf': [
        ('gene_id', 'gene_id'), ('transcript_id', 'transcript_id'), ('E_value', 'evalue'),
        ('seqid', 'seqid'), ('type', 'type'), ('score', 'score'),
    ]}

copyattrs_inv = {
    'gff': {v: k for k, v in copyattrs['gff']},
    'gtf': {v: k for k, v in copyattrs['gtf']},
    }

_QUOTE_MAP = {chr(n): f'%{n:02X}' for n in (37, 38, 44, 59, 61, 127) + tuple(range(32))}


def _quote_gff3(s):
    # Alternatively re.sub could be used for a single iteration over the string
    for char, quoted_char in _QUOTE_MAP.items():
        s = s.replace(char, quoted_char)
    return s


@_add_fmt_doc('read_fts')
def read_fts_gff(f, **kw):
    """
    Read a GFF file and return `.FeatureList`

    :param list filt: Return only Features of given ftypes, default: all
    :param str filt_fast: Read only lines which include this string
    :param str default_ftype: default ftype for entries without type
    :param list comments: comments inside the file are stored in
        the comments list (optional)
    :param str gff_version: The GFF version of the file,
        one of ``'1'``, ``'2'``, ``'3'``, by default it is auto-detected
        from the header
    """
    return _read_fts_gxf(f, flavor='gff', **kw)


@_add_fmt_doc('read_fts')
def read_fts_gtf(f, **kw):
    """
    Read a GTF file and return `.FeatureList`

    :param list filt: Return only Features of given ftypes, default: all
    :param str filt_fast: Read only lines which include this string
    :param str default_ftype: default ftype for entries without type
    :param list comments: comments inside the file are stored in
        the comments list (optional)

    .. note::
        The GTF reader is experimental. Alternatively, convert your file
        to GFF using `AGAT`_ and read that instead.
    """
    return _read_fts_gxf(f, flavor='gtf', **kw)


def _read_fts_gxf(f, filt=None, filt_fast=None, default_ftype=None, comments=None,
                  gff_version=None, flavor='gff'):
    assert flavor in ('gff', 'gtf')
    assert gff_version in (None, '1', '2', '3')
    if flavor == 'gff' and gff_version is None:
        gff_version = _detect_gff_version(f)
        if gff_version is None:
            warn('GFF version could not be auto-detected from header, assume version 3')
            gff_version = '3'
    fts = []
    lastid = None
    for lineno, line in enumerate(f):
        try:
            line2 = line
            line = line.strip()
            if flavor == 'gff' and line.startswith('##FASTA'):
                break
            if filt_fast is not None and filt_fast.lower() not in line.lower():
                continue
            if '#' in line:
                line, comment = line.split('#', 1)
                line = line.strip()
                if comments is not None:
                    comments.append('#' + comment)
            if  line == '':
                continue
            cols = line.split('\t')
            if len(cols) not in (8, 9):
                # be forgiving here and also except any whitespace instead of tab
                cols = line.split(maxsplit=8)
            if len(cols) == 9:
                *cols, attrcol = cols
            else:
                attrcol = ''
            if flavor == 'gff' and gff_version == '3':
                cols = list(map(unquote, cols))
            seqid, source, type_, start, stop, score, strand, phase = cols
            if type_ == '.':
                type_ = default_ftype
            if filt and type_ not in filt:
                continue
            loc = Location(int(start)-1, int(stop), strand=strand)
            attrs = {}
            attrcol = attrcol.strip()
            if attrcol not in ('', '.'):
                if flavor == 'gff' and gff_version == '1':
                    attrs['group'] = attrcol
                elif flavor == 'gff' and gff_version == '2' or flavor == 'gtf':
                    for kv in attrcol.split(';'):
                        if kv.strip() == '':
                            continue
                        elif '"' in kv:
                            k, v = kv.strip().split('"', maxsplit=1)
                            if v[-1] == '"':
                                v = v[:-1]
                        else:
                            k, v = kv.split(maxsplit=1)
                        attrs[k.strip()] = v
                elif flavor == 'gff' and gff_version == '3':
                    for kv in attrcol.split(';'):
                        if kv.strip() == '':
                            continue
                        k, v = kv.strip().split('=')
                        attrs[unquote(k.strip())] = (
                            unquote(v.strip()) if ',' not in v else
                            [unquote(vv.strip()) for vv in v.strip().split(',')])
                else:
                    assert False
                for k in ('evalue', 'E_value'):
                    if k in attrs:
                        try:
                            attrs[k] = float(attrs[k])
                        except ValueError:
                            pass
            if seqid != '.':
                attrs['seqid'] = seqid
            if source != '.':
                attrs['source'] = source
            if score != '.':
                attrs['score'] = float(score)
            if phase != '.':
                if flavor == 'gff':
                    attrs['phase'] = int(phase)
                else:
                    attrs['frame'] = int(phase)
            if 'ID' in attrs:
                id_ = (attrs['ID'], type_, seqid)
            else:
                id_ = None
            mf = '_' + flavor
            if id_ is not None and id_ == lastid:
                for k, v in attrs.items():
                    if fts[-1].meta[mf].get(k) != v:
                        # TODO: add proper meta attribute to Location?
                        if not hasattr(loc.meta, mf):
                            loc.meta[mf] = Attr()
                        loc.meta[mf][k] = v
                fts[-1].locs = fts[-1].locs + (loc,)
            else:
                meta = Meta({mf: attrs})
                fts.append(Feature(type_, locs=[loc], meta=meta))
            lastid = id_
        except Exception as ex:
            msg = (f'Parsing error on line number {lineno+1} of '
                   f"{flavor.upper()}{gff_version if flavor == 'gff' else ''} file.\n"
                   'If the detected file format is not correct, please use the fmt parameter. '
                   'If the file is valid, or if you think this case should be covered, please contact the developers. '
                   f'Problematic line:\n{line}')
            raise ValueError(msg) from ex
    for ft in fts:
        for gxfattr, metaattr in copyattrs[flavor]:
            if gxfattr in ft.meta[mf]:
                ft.meta[metaattr] = ft.meta[mf][gxfattr]
    return FeatureList(fts)


@_add_fmt_doc('read')
def read_gff(f, **kw):
    r"""
    Read sequences and their features from GFF file

    :param \*\*kw: All kwargs are passed to `read_fts_gff()`
    """
    fts = read_fts_gff(f, **kw)
    seqs = read_seqs(f, fmt='fasta')
    seqs.fts = fts
    return seqs


@_add_fmt_doc('write_fts')
def write_fts_gff(fts, f, **kw):
    """
    Write features to GFF file

    :param bool header_sugar: Add a comment to the header with sugar version, default is True
    :param str header: Optionally write additional header at the top of file
    """
    _write_fts_gxf(fts, f, flavor='gff', **kw)


@_add_fmt_doc('write_fts')
def write_fts_gtf(fts, f, **kw):
    """
    Write features to GTF file

    :param bool header_sugar: Add a comment to the header with sugar version, default is True
    :param str header: Optionally write additional header at the top of file

    .. note::
        The GTF writer is experimental. For complicated files use the GFF format and convert with `AGAT`_.
    """
    _write_fts_gxf(fts, f, flavor='gtf', **kw)


def _write_fts_gxf(fts, f, header=None, header_sugar=True, flavor='gff'):
    assert flavor in ('gff', 'gtf')
    if flavor == 'gff':
        f.write('##gff-version 3\n')
    if header_sugar:
        from sugar import __version__
        f.write(f'# {flavor.upper()} written by sugar v{__version__}\n')
    mf = '_' + flavor
    if header:
        f.write(header)
    for ft in fts:
        meta = ft.meta.copy()
        meta.setdefault(mf, {})
        for metaattr, gxfattr in copyattrs_inv[flavor].items():
            if metaattr in meta:
                meta[mf][gxfattr] = meta[metaattr]
        if flavor == 'gff' and ft.locs[0]._meta and hasattr(ft.locs[0].meta, '_gff'):
            meta._gff = {k: v for k, v in meta._gff.items() + ft.locs[0].meta._gff.items()}
        if flavor == 'gff' and len(ft.locs) > 1 and 'ID' not in meta._gff:
            # TODO warn, test
            from random import choices
            from string import ascii_lowercase
            meta._gff.ID = ''.join(choices(ascii_lowercase, k=10))
        if flavor == 'gtf' and len(ft.locs) > 1:
            warn('Only a single location is supported when writing GTF, use the first location')
        seqid = meta[mf].pop('seqid', '.')
        source = meta[mf].pop('source', '.')
        if flavor == 'gff':
            seqid = _quote_gff3(seqid)
            source = _quote_gff3(source)
        score = meta[mf].pop('score', '.')
        phase = meta[mf].pop('phase' if flavor == 'gff' else 'frame', '.')
        type_ = meta[mf].pop('type', None) or meta.get('type') or '.'
        for i, loc in enumerate(ft.locs):
            if i == 0:
                nscore = score
                nphase = phase
                gxf_meta = meta[mf]
            elif flavor == 'gff':
                if loc._meta and hasattr(loc.meta, '_gff'):
                    gxf_meta = {k: v for k, v in loc.meta._gff.items() if meta._gff.get(k) != v}
                else:
                    gxf_meta = {}
                gxf_meta['ID'] = meta._gff.ID
                nscore = gxf_meta.pop('score', score)
                nphase = gxf_meta.pop('phase', phase)
            if len(gxf_meta) == 0:
                attrstr = '.'
            elif flavor == 'gff':
                attrstr = ';'.join(_quote_gff3(k) + '=' + (','.join(_quote_gff3(vv) for vv in v) if isinstance(v, (list, tuple)) else
                                                     _quote_gff3(str(v)))
                                   for k, v in gxf_meta.items())
            else:
                assert flavor == 'gtf'
                attrstr = '; '.join(k + ' "' + (' '.join(str(vv) for vv in v) if isinstance(v, (list, tuple)) else
                                               str(v) + '"')
                                   for k, v in gxf_meta.items())
            f.write(f'{seqid}\t{source}\t{type_}\t{loc.start+1}\t{loc.stop}\t{nscore}\t{loc.strand}\t{nphase}\t{attrstr}\n')


@_add_fmt_doc('write')
def write_gff(seqs, f, **kw):
    r"""
    Write sequences and their features to GFF file using a FASTA directive

    :param \*\*kw: All kwargs are passed to `write_fts_gff()`
    """
    write_fts_gff(seqs.fts, f, **kw)
    f.write('##FASTA\n')
    write_seqs(seqs, f, fmt='fasta')
