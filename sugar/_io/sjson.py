# (C) 2024, Tom Eulenfeld, MIT license
"""
SJson IO, custom lossless sugar format
"""
import sugar
from sugar import __version__
import json


from sugar.core.fts import Location, Defect, Strand, Feature, FeatureList
from sugar.core.meta import Attr, Meta
from sugar.core.seq import BioBasket, BioSeq
from sugar._io.util import _add_fmt_doc


SUGAR = (Location, Defect, Strand, Feature, FeatureList,
         Attr, Meta,
         BioBasket, BioSeq
         )

filename_extensions = ['sjson', 'json']
COMMENT = f'sugar JSON format written by sugar v{__version__}'


class _SJSONEncoder(json.JSONEncoder):
    def default(self, o):
        if isinstance(o, (Strand, Defect)):
            obj = {'_cls': type(o).__name__,
                   'value': o.value}
            return obj
        elif isinstance(o, SUGAR):
            obj = {k: v for k, v in o.__dict__.items()
                   if not k.startswith("_") or k == '_fmtcomment'}
            obj['_cls'] = type(o).__name__
            return obj
        else:
            # Let the base class default method raise the TypeError
            return json.JSONEncoder.default(self, o)


def _json_hook(d):
    if cls := d.pop('_cls', None):
        d.pop('_fmtcomment', None)
        cls = globals()[cls]
        if isinstance(cls, (Strand, Defect)):
            return cls(d['value'])
        else:
            return cls(**d)
    else:
        return d


def is_format(f, **kw):
    content = f.read(51)
    return COMMENT[:17].lower() in content.lower()


@_add_fmt_doc('read')
def read(f):
    """
    Read SJson file
    """
    return json.load(f, object_hook=_json_hook)


@_add_fmt_doc('write')
def write(seqs, f):
    """
    Write sequences into SJson file
    """
    seqs.__dict__ = dict(_fmtcomment=COMMENT,
                         **seqs.__dict__)
    json.dump(seqs, f, cls=_SJSONEncoder)
    del seqs._fmtcomment
