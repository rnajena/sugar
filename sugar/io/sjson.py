# (C) 2024, Tom Eulenfeld, MIT license
import sugar
from sugar import __version__
from sugar.core.seq import SUGAR
import json


EXT = ['sjson', 'json']
COMMENT = f'sugar JSON format written by sugar v{__version__}'


class SJSONEncoder(json.JSONEncoder):
    def default(self, o):
        if isinstance(o, SUGAR):
            obj = {k: v for k, v in o.__dict__.items()
                   if not k.startswith("_") or k == '_fmtcomment'}
            obj['_cls'] = type(o).__name__
            return obj
        else:
            # Let the base class default method raise the TypeError
            return json.JSONEncoder.default(self, o)


def json_hook(d):
    if cls := d.pop('_cls', None):
        d.pop('_fmtcomment', None)
        return getattr(sugar, cls)(**d)
    else:
        return d


def is_format(f):
    content = f.read(51)
    return COMMENT[:17].lower() in content.lower()


def read(f):
    return json.load(f, object_hook=json_hook)


def write(seqs, f):
    seqs.__dict__ = dict(_fmtcomment=COMMENT,
                         **seqs.__dict__)
    json.dump(seqs, f, cls=SJSONEncoder)
    del seqs._fmtcomment
