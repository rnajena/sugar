# (C) 2024, Tom Eulenfeld, MIT license
"""
`BLAST`_ reader for output generated with outfmt 7 (preferred), 6 or 10
"""

from sugar._io.tab.core import read_tabular
from sugar._io.util import _add_fmt_doc


def is_fts_blast(f, outfmt=None, **kw):
    content = f.read(1000)
    if content.startswith('#') and 'BLAST' in content:
        return True
    # just try to read first line
    line = content.splitlines()[0]
    fts = read_fts_blast([line], outfmt=outfmt, **kw)
    return len(fts) == 1 and (outfmt is not None or 0 <= fts[0].meta._blast.pident <= 100)


@_add_fmt_doc('read_fts')
def read_fts_blast(f, *, sep='\t', outfmt=None, ftype=None, comments=None):
    """
    BLAST reader for output generated with option outfmt 7 (preferred), 6, or 10

    :param str sep: Separator of fields, use ``','`` for outfmt 10, default ``'\\t'``,
        can be set to ``None`` for any whitespace.
    :param str outfmt: The outfmt string passed to BLAST, can be omitted for outfmt 7
        or default outfmt 6 or 10 output.
    :param str ftype: Parameter used as ftype, if parameter is not present use the value
        of ftype itself.
    :param list comments: comment lines inside the file are stored in
        the comments list (optional)
    """
    return read_tabular(f, sep=sep, outfmt=outfmt, ftype=ftype, comments=comments, fmt='blast')
