# (C) 2023, Tom Eulenfeld, MIT license
"""
Genbank reader

The format is defined `here <https://www.insdc.org/submitting-standards/feature-table/>`_.
"""

from sugar.core.seq import Attr, BioSeq, Meta, Feature, FeatureList
from sugar._io._genbank_loc import _parse_locs
from sugar._io.util import _add_fmt_doc


def is_format(f, **kw):
    content = f.read(5)
    return content.lower() == 'locus'


is_format_fts = is_format


@_add_fmt_doc('read_fts')
def read_fts(f, exclude=()):
    """
    Read Genbank feature records from file into `.FeatureList`

    :param tuple exclude: Tuple of feature names to exclude,
        possible options are ``'translation', 'features'``,
        sequences are excluded anyway.
    """
    fts = FeatureList()
    for seq in iter_(f, exclude=('seq',) + exclude):
        fts.extend(seq.fts)
    return fts


@_add_fmt_doc('read')
def iter_(f, exclude=()):
    """
    Read Genbank records and sequences from file into `.BioBasket`

    :param tuple exclude: Tuple of feature names to exclude,
        possible options are ``'seq', 'translation', 'features'``.

    """
    # allowed entries in exclude: features, translation, seq
    meta = Meta()
    attrs = Attr()
    features = []
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
            assert len(misc) == 0
            assert fttype is None
            meta._genbank = attrs
            if 'accession' in meta._genbank:
                meta.id = meta._genbank.accession.split()[0]
                for ft in meta.features:
                    ft.meta.seqid = meta.id
            try:
                del meta._genbank.reference  # TODO: references should be parsed in a list, not yet done
            except Exception:
                pass
            if 'translation' in exclude and 'features' in meta:
                for feature in meta.features:
                    try:
                        del feature.meta._genbank.translation
                    except Exception:
                        pass
            yield BioSeq(seq.upper(), meta=meta)
            meta = Meta()
            attrs = Attr()
            features = []
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
                    parse = 'features'
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
        elif parse == 'features':
            # if 'features' in exclude and 'seq' in exclude:
            #     continue
            if 'features' in exclude:
                continue
            if len(line[:20].strip()) > 0:
                if fttype is not None:
                    ft = Feature(type=fttype, locs=_parse_locs(locs),
                                 meta=Meta(_genbank=ftmeta))
                    features.append(ft)
                    fttype = None
                key = line[:20].strip().split()[0]
                # subkey = None
                if key.lower() == 'origin':
                    parse = 'origin'
                    meta.features = FeatureList(features)
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
