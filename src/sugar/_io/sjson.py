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

filename_extensions_sjson = ['sjson', 'json']
filename_extensions_fts_sjson = ['sjson', 'json']
COMMENT = f'sugar JSON format written by sugar v{__version__}'
COMMENT_FTS = f'sugar JSON feature format written by sugar v{__version__}'


def _SJSONEncoder_factory(private=False):
    class _SJSONEncoder(json.JSONEncoder):
        _private = private
        def default(self, o):
            if isinstance(o, (Strand, Defect)):
                obj = {'_cls': type(o).__name__,
                    'value': o.value}
                return obj
            elif isinstance(o, Location):
                obj = {'_cls': type(o).__name__,
                    'start': o.start,
                    'stop': o.stop,
                    'strand': o.strand,
                    'defect': o.defect
                    }
                if len(o.meta) > 0:
                    obj['meta'] = o.meta
                return obj
            elif isinstance(o, SUGAR):
                obj = {}
                if '_fmtcomment' in o.__dict__:
                    obj['_fmtcomment'] = o.__dict__['_fmtcomment']
                obj['_cls'] = type(o).__name__
                if isinstance(o, Feature):
                    obj['locs'] = o.locs
                obj.update({
                    k: v for k, v in o.__dict__.items()
                    if not k.startswith('_') or
                    self._private and not k.startswith('__') and isinstance(o, Meta)
                    })
                if 'meta' in obj and len(obj['meta']) == 0:
                    del obj['meta']
                return obj
            else:
                # Let the base class default method raise the TypeError
                return json.JSONEncoder.default(self, o)
    return _SJSONEncoder


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


def is_sjson(f, **kw):
    content = f.read(51)
    return COMMENT[:17].lower() in content.lower()


def is_fts_sjson(f, **kw):
    content = f.read(51)
    return COMMENT_FTS[:25].lower() in content.lower()


@_add_fmt_doc('read')
def read_sjson(f):
    """
    Read SJson file

    .. note::
        You can use this function directly to load arbitrary
        objects containing sugar objects.
        Just ignore the warning and use a file descriptor ;)
    """
    return json.load(f, object_hook=_json_hook)


@_add_fmt_doc('read_fts')
def read_fts_sjson(f):
    """
    Read features from SJson file
    """
    return json.load(f, object_hook=_json_hook)


@_add_fmt_doc('write')
def write_sjson(seqs, f, *, private=False, indent=None):
    """
    Write sequences into SJson file

    :param bool private: Also write private metadata (mostly format-related, default False)
    :param int indent: Indent in JSON file.

    .. note::
        You can use this function directly to write arbitrary
        objects containing sugar objects.
        Just ignore the warning and use a file descriptor ;)

        The following example writes an object to JSON and reads it again. ::

            from sugar import read
            from sugar._io.sjson import read_sjson, write_sjson
            seqs = read()
            with open('test.json', 'w') as f:
                write_sjson([1, seqs], f)
            with open('test.json') as f:
                obj2 = read_sjson(f)
    """
    try:
        seqs._fmtcomment = COMMENT
    except AttributeError:
        pass
    json.dump(seqs, f, cls=_SJSONEncoder_factory(private=private), indent=indent)
    try:
        del seqs._fmtcomment
    except AttributeError:
        pass


@_add_fmt_doc('write_fts')
def write_fts_sjson(fts, f, *, private=False, indent=None):
    """
    Write features into SJson file

    :param bool private: Also write private metadata (mostly format-related, default False)
    :param int indent: Indent in JSON file.
    """
    fts._fmtcomment = COMMENT_FTS
    json.dump(fts, f, cls=_SJSONEncoder_factory(private=private), indent=indent)
    del fts._fmtcomment
