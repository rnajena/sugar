# (C) 2023, Tom Eulenfeld, MIT license
# https://www.insdc.org/submitting-standards/feature-table/

from sugar.core.seq import Attr, BioSeq, Meta, Feature, FeatureList
from sugar.io.genbank_loc import _parse_locs


def is_format(f):
    content = f.read(5)
    return content.lower() == 'locus'

is_format_fts = is_format

# def _parse_feature_location(loc):
#     try:
#         # TODO, this does not cover all stuff
#         # join, complement, <, >, ^, etc
#         start, stop = map(int, loc.split('..'))
#     except Exception:
#         return {}
#     else:
#         assert start <= stop
#         return dict(start=start-1, stop=stop)
        # if start <= stop:
        #     kw = dict(start=start-1, stop=stop)
        # else:
        #     kw = dict(start=start-1, stride=-1,
        #               stop=stop-2 if stop > 1 else None)


def read_fts(f, exclude=()):
    fts = FeatureList()
    for seq in iter_(f, exclude=('seq',) + exclude):
        fts.extend(seq.fts)
    return fts


def iter_(f, exclude=()):
    # allowed entries in exclude: features, translation, seq
    content = f.read()
    for record in content.split('//'):
        if len(record.strip()) == 0:
            continue
        meta = Meta()
        features = []
        misc = []
        fttype = None
        ftmeta=None
        locs = None
        val = None
        key = None
        subkey = None
        parse = 'header'
        seq = ''
        for line in record.splitlines():
            if line.strip() == '':
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
                    meta[key] = val
                elif line.startswith(' ' * 12) and len(line.strip()) > 0:
                    if subkey:
                        meta[key][subkey] = meta[key][subkey] + ' ' + line.strip()
                    else:
                        meta[key] = meta[key] + ' ' + line.strip()
                elif line.startswith(' ') and len(line.strip()) > 0:
                    subkey = line[:12].lower().strip()
                    try:
                        val = line.strip().split(maxsplit=1)[1]
                    except Exception:
                        val = ''
                    meta[key] = Attr(id=meta[key])
                    meta[key][subkey] = val
            elif parse == 'features':
                if 'features' in exclude and 'seq' in exclude:
                    break
                elif 'features' in exclude:
                    continue
                if len(line[:20].strip()) > 0:
                    if fttype is not None:
                        ft = Feature(type=fttype, locs=_parse_locs(locs),
                                          meta=ftmeta)
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
                    break
                if len(line) > 10:
                    seq = seq + line[10:].replace(' ', '')
            else:
                assert False
        assert len(misc) == 0
        assert fttype is None
        if 'accession' in meta:
            meta.id = meta.accession.split()[0]
            del meta.accession
        try:
            del meta.reference  # TODO: references should be parsed in a list, not yet done
        except Exception:
            pass
        if 'translation' in exclude and 'features' in meta:
            for feature in meta.features:
                try:
                    del feature.meta.translation
                except Exception:
                    pass
        yield BioSeq(seq.upper(), meta=meta)
