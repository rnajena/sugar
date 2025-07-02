# (C) 2023, Tom Eulenfeld, MIT license
"""
`EMBL`_ flat file reader for ENA and UniProt

``EMBL`` stands for the European Molecular Biology Laboratory.
``ENA`` stands for the European Nucleotide Archive.
``UniProt`` stands for the Universal Protein Database.
"""

from sugar.core.fts import Defect, Feature, FeatureList, Location
from sugar.core.seq import Attr, BioSeq, Meta
from sugar._io.util import _add_fmt_doc
from sugar._io.genbank import _parse_locs


def is_embl(f, **kw):
    content = f.read(5)
    return content == 'ID   '


is_fts_embl = is_embl


_EMBL2GENBANK = {
    'ID': 'locus',
    'AC': 'accession',
    'DT': 'date',
    'DE': 'definition',
    'GN': 'gene_name',
    'KW': 'keywords',
    'DR': 'dbsource',
    'CC': 'comment',
    'OS': 'organism',
    'OC': 'taxonomy',
    'OG': 'organelle',
    'PA': 'parent_accession',
    'PR': 'project',
    'R': 'reference',
    }


@_add_fmt_doc('read_fts')
def read_fts_embl(f, exclude=()):
    """
    Read EMBL feature records from file into `.FeatureList`

    :param tuple exclude: Tuple of feature names to exclude,
        possible options are ``'translation', 'fts'``,
        sequences are excluded anyway.
    """
    fts = FeatureList()
    for seq in iter_embl(f, exclude=('seq',) + exclude):
        fts.extend(seq.fts)
    return fts


@_add_fmt_doc('read')
def iter_embl(f, exclude=(), genbank=True):
    """
    Read EMBL records and sequences from file into `.BioBasket`

    :param tuple exclude: Tuple of feature names to exclude,
        possible options are ``'seq', 'translation', 'fts'``
        or line keys (e.g. ``'CC'``).
    :param genbank: By default, use genbank like key names in the ``_embl`` meta data,
        if set to ``False`` will use EMBL two character keys, except for references which will be
        saved in the ``_embl.R`` attribute.
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
    parse = 'header'
    seq = ''
    for line in f:
        line = line.rstrip()
        if fttype is not None and (not line.startswith('FT') or len(line[2:20].strip()) > 0):
            # create a new feature
            ft = Feature(type=fttype, locs=_parse_locs(locs),
                            meta=Meta(_embl=ftmeta))
            fts.append(ft)
            fttype = None
        if line.strip() == '' or line[:2] in exclude or line[:2] in ('XX', 'FH'):
            continue
        if line.strip() == '//':
            # create a new sequence object
            assert len(misc) == 0
            assert fttype is None
            if 'AC' in attrs:
                meta.id = attrs['AC'].split(';')[0].strip()
                if 'fts' in meta:
                    for ft in meta.fts:
                        ft.meta.seqid = meta.id
            if genbank:
                meta._embl = {gbkey: attrs[emblkey] for emblkey, gbkey in _EMBL2GENBANK.items() if emblkey in attrs}
            else:
                meta._embl = attrs
            meta.fts = FeatureList(fts)
            if 'translation' in exclude and 'fts' in meta:
                for feature in meta.fts:
                    try:
                        del feature.meta._embl.translation
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
            parse = 'header'
            seq = ''
            continue
        elif line.startswith('FT'):
            parse = 'fts'
        elif line.startswith('SQ'):
            parse = 'origin'
        if parse == 'header':
            key = line[:2] if not line.startswith('R') else 'R'
            line = line[2:].strip()
            if line == '':
                continue
            if key == 'ID':
                line = ' '.join(line.split())
            attrs[key] = line if key not in attrs else attrs[key] + ' ' + line
        elif parse == 'fts':
            if 'fts' in exclude:
                continue
            line = line[2:]
            if len(line[:18].strip()) > 0:
                key = line[:18].strip().split()[0]
                fttype = line[:18].strip().split()[0]
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
        elif parse in ('origin', 'ena', 'uniprot'):
            if 'seq' not in exclude and not line.startswith('SQ') and line.strip() != '':
                if parse == 'origin':
                    if line[-7:-5] == '  ':
                        parse = 'ena'  # char counts at the end of line
                        try:
                            int(line[-1])
                        except Exception:
                            from warnings import warn
                            warn('Inconsistent EMBL file, please check the read content, contact developers')
                    else:
                        parse = 'uniprot'
                if parse == 'ena':
                    line = line[:-10]
                seq = seq + line.replace(' ', '')
        else:
            assert False
