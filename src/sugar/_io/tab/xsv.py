"""
TSV, CSV and XSV file formats IO

TSV and CSV use the same functions under the hood,
and just use a different default value for the ``sep`` option defining the separator character.
In general, arbitrary characters can be used to separate the values (XSV).
"""
import csv
from sugar._io.util import _add_fmt_doc
from sugar import FeatureList


filename_extensions_fts_csv = ['csv']
filename_extensions_fts_tsv = ['tsv']


def _is_fts_xsv(f, sep, **kw):
    lines = f.read(1000).splitlines()[:-1]
    line0 = lines[0].split(sep)
    return ('start' in line0) + ('stop' in line0) + ('len' in line0) >= 2 and all(len(line.split(sep)) == len(line0) for line in lines)


def is_fts_csv(f, sep=',', **kw):
    return _is_fts_xsv(f, sep)


def is_fts_tsv(f, sep='\t', **kw):
    return _is_fts_xsv(f, sep)


@_add_fmt_doc('read_fts')
def read_fts_tsv(f, sep='\t', **kw):
    r"""
    Read tab separated value (TSV) files with feature information

    :param str sep: Separator of the fields, defaults to ``'\t'``
    :param str ftype: Parameter used as feature type, if parameter is not present use the value
        of ftype itself.
    :param \*\*kw: All other kwargs are passed to :func:`pandas.read_csv()`

    """
    return _read_fts_xsv(f, sep=sep, **kw)


@_add_fmt_doc('read_fts')
def read_fts_csv(f, sep=',', **kw):
    r"""
    Read comma separated value (CSV) files with feature information

    :param str sep: Separator of the fields, defaults to ``','``
    :param str ftype: Parameter used as feature type, if parameter is not present use the value
        of ftype itself.
    :param \*\*kw: All other kwargs are passed to :func:`pandas.read_csv()`
    """
    return _read_fts_xsv(f, sep=sep, **kw)


@_add_fmt_doc('write_fts')
def write_fts_tsv(fts, f, sep='\t', **kw):
    r"""
    Write tab separated value (TSV) files with feature information

    :param str sep: Separator of the fields, defaults to ``'\t'``
    :param vals: Parameters from the metadata or location to write,
        ``'len'`` is also allowed,
        might be a string or tuple, defaults to ``'type start stop strand'``
    :param \*\*kw: All other kwargs are passed to :meth:`pandas.DataFrame.to_csv()`

    """
    return _write_fts_xsv(fts, f, sep=sep, **kw)


@_add_fmt_doc('write_fts')
def write_fts_csv(fts, f, sep=',', **kw):
    r"""
    Write comma separated value (CSV) files with feature information

    :param str sep: Separator of the fields, defaults to ``','``
    :param vals: Parameters from the metadata or location to write,
        ``'len'`` is also allowed,
        might be a string or tuple, defaults to ``'type start stop strand'``
    :param \*\*kw: All other kwargs are passed to :meth:`pandas.DataFrame.to_csv()`
    """
    return _write_fts_xsv(fts, f, sep=sep, **kw)


def _read_fts_xsv(f, ftype=None, **kw):
    try:
        import pandas
    except ImportError:
        raise ImportError('TSV/CSV reader depends on pandas module')
    df = pandas.read_csv(f, **kw)
    return FeatureList.frompandas(df, ftype=ftype)


def _write_fts_xsv(fts: FeatureList, f, keys='type start stop strand', dtype=None, index=False, **kw):
    try:
        df = fts.topandas(keys, dtype=dtype)
    except ImportError:
        raise ImportError('TSV/CSV writer depends on pandas module')
    return df.to_csv(f, index=index, **kw)
