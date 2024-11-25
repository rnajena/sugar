"""
TSV, CSV and similar file formats IO
"""
import csv
from sugar._io.util import _add_fmt_doc
from sugar import FeatureList


filename_extensions_fts_tsv = ['tsv']
filename_extensions_fts_csv = ['csv']


def is_fts_tsv(f, **kw):
    sep = kw.pop('sep', '\t')
    lines = f.read(1000).splitlines()[:-1]
    return 'start' in lines[0] and 'stop' in lines[0] and all(sep in line for line in lines)


def is_fts_csv(f, **kw):
    sep = kw.pop('sep', ',')
    lines = f.read(1000).splitlines()[:-1]
    return 'start' in lines[0] and 'stop' in lines[0] and all(sep in line for line in lines)


@_add_fmt_doc('read_fts')
def read_fts_tsv(f, sep='\t', **kw):
    r"""
    Read tab separated value (TSV) files with feature information

    :param str sep: Separator of the fields, defaults to ``'\t'``
    :param str ftype: Parameter used as feature type, if parameter is not present use the value
        of ftype itself.
    :param \*\*kw: All other kwargs are passed to `pandas.read_csv()`

    """
    return _read_xsv(f, sep=sep, **kw)


@_add_fmt_doc('read_fts')
def read_fts_csv(f, sep=',', **kw):
    r"""
    Read comma separated value (CSV) files with feature information

    :param str sep: Separator of the fields, defaults to ``','``
    :param str ftype: Parameter used as feature type, if parameter is not present use the value
        of ftype itself.
    :param \*\*kw: All other kwargs are passed to `pandas.read_csv()`
    """
    return _read_xsv(f, sep=sep, **kw)


@_add_fmt_doc('write_fts')
def write_fts_tsv(fts, f, sep='\t', **kw):
    r"""
    Write tab separated value (TSV) files with feature information

    :param str sep: Separator of the fields, defaults to ``'\t'``
    :param vals: Parameters from the metadata or location to write,
        might be a string or tuple, defaults to ``'type start stop strand'``
    :param \*\*kw: All other kwargs are passed to `pandas.DataFrame.to_csv()`

    """
    return _write_xsv(fts, f, sep=sep, **kw)


@_add_fmt_doc('write_fts')
def write_fts_csv(fts, f, sep=',', **kw):
    r"""
    Write comma separated value (CSV) files with feature information

    :param str sep: Separator of the fields, defaults to ``','``
    :param vals: Parameters from the metadata or location to write,
        might be a string or tuple, defaults to ``'type start stop strand'``
    :param \*\*kw: All other kwargs are passed to `pandas.DataFrame.to_csv()`
    """
    return _write_xsv(fts, f, sep=sep, **kw)


def _read_xsv(f, ftype=None, **kw):
    try:
        import pandas
    except ImportError:
        raise ImportError('TSV/CSV reader depends on pandas module')
    df = pandas.read_csv(f, **kw)
    return FeatureList.frompandas(df, ftype=ftype)


def _write_xsv(fts: FeatureList, f, vals='type start stop strand', dtype=None, index=False, **kw):
    try:
        df = fts.topandas(vals, dtype=dtype)
    except ImportError:
        raise ImportError('TSV/CSV writer depends on pandas module')
    return df.to_csv(f, index=index, **kw)


