# (C) 2023, Tom Eulenfeld, MIT license
# https://www.insdc.org/submitting-standards/feature-table/

from sugar.core.seq import Attr, BioSeq, Meta, Feature, FeatureList


def is_format(f):
    content = f.read(5)
    return content.lower() == 'locus'


def _parse_feature_location(loc):
    try:
        # TODO, this does not cover all stuff
        # join, complement, <, >, ^, etc
        start, stop = map(int, loc.split('..'))
    except Exception:
        return {}
    else:
        assert start < stop
        return dict(start=start-1, stop=stop)
        # if start <= stop:
        #     kw = dict(start=start-1, stop=stop)
        # else:
        #     kw = dict(start=start-1, stride=-1,
        #               stop=stop-2 if stop > 1 else None)


def iter_(f, translation=True):
    content = f.read()
    for record in content.split('//'):
        if len(record.strip()) == 0:
            continue
        meta = Meta()
        features = []
        misc = []
        feature = None
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
                if len(line[:20].strip()) > 0:
                    if feature is not None:
                        features.append(feature)
                        feature = None
                    key = line[:20].strip().lower()
                    # subkey = None
                    if key == 'origin':
                        parse = 'origin'
                        meta.features = FeatureList(features)
                        continue
                    val = line.strip()
                    try:
                        val = val.split(maxsplit=1)[1]
                    except Exception:
                        pass
                    else:
                        feature = Feature(type=key, loc=val,
                                          **_parse_feature_location(val))
                elif len(line.strip()) > 0:
                    line = line.strip()
                    if line.startswith('/'):
                        if '=' in line:
                            key, val = line[1:].split('=')
                            if not val.startswith('"'):
                                try:
                                    val = int(val)
                                except Exception:
                                    pass
                            else:
                                val = val.strip('"')
                            feature[key] = val
                        else:
                            feature.setdefault('misc', []).append(line[1:])
                    else:
                        feature[key] = feature[key] + line.strip('"')
            elif parse == 'origin':
                if len(line) > 10:
                    seq = seq + line[10:].replace(' ', '')
            else:
                assert False
        assert len(misc) == 0
        assert feature is None
        if 'accession' in meta:
            meta.id = meta.accession
            del meta.accession
        try:
            del meta.reference  # TODO: references should be parsed in a list, not yet done
        except Exception:
            pass
        if not translation and 'features' in meta:
            for feature in meta.features:
                try:
                    del feature.translation
                except Exception:
                    pass
        yield BioSeq(seq.upper(), meta=meta)
