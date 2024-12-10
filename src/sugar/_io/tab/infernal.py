# (C) 2024, Tom Eulenfeld, MIT license
"""
`Infernal`_ reader for output generated with tblout fmt 1, 2, 3
"""

from sugar._io.tab.core import read_tabular
from sugar._io.util import _add_fmt_doc


_HEADER_KW = set("""idx target name          accession query name
accession clan name mdl mdl from   mdl to seq from   seq to strand trunc pass
gc  bias  score   E-value inc olp anyidx afrct1 afrct2 winidx wfrct1 wfrct2
mdl len seq len description of target""".split())


def is_fts_infernal(f, **kw):
    content = f.read(1000)
    lines = content.splitlines()
    return (set(lines[0].lstrip('#').split()) <= _HEADER_KW and
            len(lines[1].lstrip('#').split()) in (18, 29, 20, 27))


@_add_fmt_doc('read_fts')
def read_fts_infernal(f, ftype=None, comments=None):
    """
    Infernal reader for output generated with tblout fmt 1, 2, 3

    :param str ftype: Parameter used as ftype
    :param list comments: comment lines inside the file are stored in
        the comments list (optional)
    """
    return read_tabular(f, sep=None, ftype=ftype, comments=comments, fmt='infernal')
