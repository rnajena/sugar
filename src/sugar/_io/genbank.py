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
    locs = []
    if loc.startswith(('join', 'order', 'complement')):
        from warnings import warn
        warn('Parsing of genbank loc join, order, complement is untested')
        # TODO: add some tests
        locs = [_parse_locs(subloc.strip())
                for subloc in
                loc[loc.index('(')+1:loc.rindex(')')].split(',')]
        locs = [l for ll in locs for l in ll]
        if loc.startswith('complement'):
            for loc in locs:
                loc.strand = {'-': '+', '+': '-'}.get(loc.strand, loc.strand)
    else:
        locs = [_parse_single_loc(loc)]
    return locs


def _parse_single_loc(loc: str):
    if ':' in loc:
        from warnings import warn
        warn('Parsing of seqids inside genbank loc fields is not yet supported, ignore the seqid')
        _, loc = loc.split(':')
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
        return Location(int(loc)-1, int(loc), defect=defect)
    # base range
    start, stop = loc.split(splitter)
    return Location(int(start)-1, int(stop), defect=defect)


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
        elif line.strip() == '//':
            if fttype is not None:
                ft = Feature(type=fttype, locs=_parse_locs(locs),
                                meta=Meta(_genbank=ftmeta))
                fts.append(ft)
                fttype = None
                meta.fts = FeatureList(fts)
            assert len(misc) == 0
            assert fttype is None
            meta._genbank = attrs
            if 'accession' in meta._genbank:
                meta.id = meta._genbank.accession.split()[0]
                if 'fts' in meta:
                    for ft in meta.fts:
                        ft.meta.seqid = meta.id
            try:
                del meta._genbank.reference  # references could be parsed in a list, not implemented
            except Exception:
                pass
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
        if parse == 'header':
            if not line.startswith(' ') and len(line)>0:
                key = line[:12].lower().strip()
                subkey = None
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
                attrs[key] = val
            elif line.startswith(' ' * 12) and len(line.strip()) > 0:
                if subkey:
                    attrs[key][subkey] = attrs[key][subkey] + ' ' + line.strip()
                else:
                    attrs[key] = attrs[key] + ' ' + line.strip()
            elif line.startswith(' ') and len(line.strip()) > 0:
                subkey = line[:12].lower().strip()
                try:
                    val = line.strip().split(maxsplit=1)[1]
                except Exception:
                    val = ''
                attrs[key] = Attr(id=attrs[key])
                attrs[key][subkey] = val
        elif parse == 'fts':
            # if 'fts' in exclude and 'seq' in exclude:
            #     continue
            if 'fts' in exclude:
                continue
            if len(line[:20].strip()) > 0:
                if fttype is not None:
                    ft = Feature(type=fttype, locs=_parse_locs(locs),
                                 meta=Meta(_genbank=ftmeta))
                    fts.append(ft)
                    fttype = None
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
