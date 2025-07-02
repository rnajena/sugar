# (C) 2023, Tom Eulenfeld, MIT license
"""
`GenBank`_ reader
"""

from sugar.core.fts import Defect, Feature, FeatureList, Location
from sugar.core.seq import Attr, BioSeq, Meta
from sugar._io.util import _add_fmt_doc


def is_genbank(f, **kw):
    content = f.read(5)
    return content.lower() == 'locus'


is_fts_genbank = is_genbank


def _parse_locs(loc: str):
    # See https://www.insdc.org/submitting-standards/feature-table/#3.4
    if loc.startswith('complement'):
        locs = _parse_locs(loc[loc.index('(')+1:loc.rindex(')')])[::-1]
        for locobj in locs:
            locobj.strand = {'-': '+', '+': '-'}.get(locobj.strand, locobj.strand)
    elif loc.startswith(('join', 'order')):
        locs = [_parse_locs(subloc.strip())
                for subloc in
                loc[loc.index('(')+1:loc.rindex(')')].split(',')]
        locs = [l for ll in locs for l in ll]
    else:
        locs = [_parse_single_loc(loc)]
    return locs


def _parse_single_loc(loc: str):
    if ':' in loc:
        from warnings import warn
        warn('Found seqid inside GenBank loc field, '
             'the seqid is saved in Location.meta._genbank.seqid. '
             'Other parts of sugar may ignore this information.')
        seqid, loc = loc.split(':')
        meta = {'_genbank': {'seqid': seqid}}
    else:
        meta = None
    defect = Defect.NONE
    if loc[0] == '<':
        defect |= Defect.BEYOND_LEFT
        loc = loc[1:]
    if '>' in loc:
        defect |= Defect.BEYOND_RIGHT
        loc = loc.replace('>', '')
    if ".." in loc:
        splitter = ".."
    elif "." in loc:
        splitter = "."
        defect |= Defect.UNKNOWN_SINGLE_BETWEEN
    elif "^" in loc:
        splitter = "^"
        defect |= Defect.BETWEEN_CONSECUTIVE
    else:
        # single base
        return Location(int(loc)-1, int(loc), defect=defect, meta=meta)
    # base range
    start, stop = loc.split(splitter)
    return Location(int(start)-1, int(stop), defect=defect, meta=meta)


@_add_fmt_doc('read_fts')
def read_fts_genbank(f, exclude=()):
    """
    Read GenBank feature records from file into `.FeatureList`

    :param tuple exclude: Tuple of feature names to exclude,
        possible options are ``'translation', 'fts'``,
        sequences are excluded anyway.
    """
    fts = FeatureList()
    for seq in iter_genbank(f, exclude=('seq',) + exclude):
        fts.extend(seq.fts)
    return fts


@_add_fmt_doc('read')
def iter_genbank(f, exclude=()):
    """
    Read GenBank records and sequences from file into `.BioBasket`

    :param tuple exclude: Tuple of feature names to exclude,
        possible options are ``'seq', 'translation', 'fts'``.

    """
    # allowed entries in exclude: features, translation, seq
    meta = Meta()
    attrs = Attr()
    fts = []
    misc = []
    fttype = None
    ftmeta = None
    locs = None
    val = None
    key = None
    subkey = None
    parse = 'header'
    seq = ''
    for line in f:
        line = line.rstrip()
        if line.strip() == '':
            continue
        if fttype is not None and (
                line.strip() == '//' or parse == 'fts' and len(line[:20].strip()) > 0):
            # create a new feature
            ft = Feature(type=fttype, locs=_parse_locs(locs),
                            meta=Meta(_genbank=ftmeta))
            fts.append(ft)
            fttype = None
            meta.fts = FeatureList(fts)
        if line.strip() == '//':
            # create a new sequence object
            assert len(misc) == 0
            assert fttype is None
            meta._genbank = attrs
            if 'accession' in meta._genbank:
                meta.id = meta._genbank.accession.split()[0]
                if 'fts' in meta:
                    for ft in meta.fts:
                        ft.meta.seqid = meta.id
            if 'translation' in exclude and 'fts' in meta:
                for feature in meta.fts:
                    try:
                        del feature.meta._genbank.translation
                    except Exception:
                        pass
            yield BioSeq(seq.upper(), meta=meta)
            meta = Meta()
            attrs = Attr()
            fts = []
            misc = []
            fttype = None
            ftmeta = None
            locs = None
            val = None
            key = None
            subkey = None
            parse = 'header'
            seq = ''
            continue
        if parse in ('header', 'reference'):
            if len(line)>0:
                if line[:12].strip() != '':
                    key = line[:12].lower().strip()
                if key == 'features':
                    parse = 'fts'
                    key = None
                    continue
                try:
                    val = line.strip().split(maxsplit=1)[1]
                except Exception:
                    val = ''
                if key == 'locus':
                    val = ', '.join(val.split())
                if key == 'reference':
                    parse = 'reference'
                if line[0] == ' ' and parse == 'reference':
                    key = 'reference'
                if key == 'organism' and line[:12].strip() == '':
                    key = 'taxonomy'
                attrs[key] = val if key not in attrs else attrs[key] + '; ' + val
        elif parse == 'fts':
            if 'fts' in exclude:
                continue
            if len(line[:20].strip()) > 0:
                assert fttype is None
                key = line[:20].strip().split()[0]
                # subkey = None
                if key.lower() == 'origin':
                    parse = 'origin'
                    meta.fts = FeatureList(fts)
                    continue
                fttype = key
                locs = ''
                ftmeta = {}
                key2 = None
                val = line.strip()
                try:
                    val = val.split(maxsplit=1)[1]
                except Exception:
                    pass
                else:
                    locs = val
            elif len(line.strip()) > 0:
                line = line.strip()
                if line.startswith('/'):
                    line = line.removeprefix('/')
                    if '=' in line:
                        key2, val = line.split('=')
                        if not val.startswith('"'):
                            try:
                                val = int(val)
                            except Exception:
                                pass
                        else:
                            val = val.strip('"')
                        ftmeta[key2] = val
                    else:
                        ftmeta.setdefault('misc', []).append(line)
                elif key2 is None:
                    # location spanning multiple lines
                    locs = locs + line
                else:
                    ftmeta[key2] = ftmeta[key2] + line.strip('"')
        elif parse == 'origin':
            if 'seq' in exclude:
                seq = ''
                continue
            if len(line) > 10:
                seq = seq + line[10:].replace(' ', '')
        else:
            assert False
