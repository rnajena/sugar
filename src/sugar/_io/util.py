# (C) 2024, Tom Eulenfeld, MIT license
"""
Helper functions for _io module
"""
from importlib.metadata import entry_points
import shutil

# Lists in this dictionary define the detection order when reading
# a file without specifying the format.
# More common and easier to detect formats should be earlier in the list.
FMTS = {'seqs': ['fasta', 'genbank', 'embl', 'clustal', 'stockholm', 'gff', 'sjson'],
        'fts': ['gff', 'gtf', 'genbank', 'embl', 'infernal', 'mmseqs', 'blast', 'tsv', 'csv', 'sjson']}


def _epsname_key(epsname, what='seqs'):
    assert what in ('seqs', 'fts')
    try:
        return FMTS[what].index(epsname)
    except ValueError:
        return len(FMTS[what])

def _epsname_fts_key(epsname):
    return _epsname_key(epsname, what='fts')


# All formats including format plugins defined in other modules
EPS = {'seqs': entry_points(group='sugar.io'),
       'fts': entry_points(group='sugar.io.fts')}
FMTS_ALL = {'seqs': sorted(EPS['seqs'].names, key=_epsname_key),
            'fts': sorted(EPS['fts'].names, key=_epsname_fts_key)}

ARCHIVE_EXTS = [ext.removeprefix('.') for fmt in shutil.get_unpack_formats()
                for ext in fmt[1]]

def _create_format_plugin_table(what, io, ws=0, fancy=False):

    assert what in ('seqs', 'fts')
    assert io in ('in', 'out', 'io')
    rc = 'i' in io  # display read columns
    wc = 'o' in io  # display write columns
    sep = ['==='] * (3 + rc + wc)
    header = ['format', 'module'] + ['read'] * rc + ['write'] * wc + ['description']
    table = [sep, header, sep]
    yes = b'\xf0\x9f\x91\x8d'.decode('utf-8')   # thumbs up
    yes = '\u2705'
    for fmt in FMTS_ALL[what]:
        module = EPS[what][fmt].load()
        row = [fmt, f'`{module.__name__}`' if fancy else module.__name__]
        if what == 'seqs':
            if rc and ((hr:=hasattr(module, f'read_{fmt}')) or hasattr(module, f'iter_{fmt}')):
                if fancy:
                    t1 = f'`{yes}<{module.__name__}.read_{fmt}()>`'
                    t2 = f'`{yes}<{module.__name__}.iter_{fmt}()>`'
                    row.append(t1 if hr else t2)
                else:
                    row.append('yes')
            elif rc and wc:
                row.append('..' if fancy else 'no')
            if wc and ((hw:=hasattr(module, f'write_{fmt}')) or hasattr(module, f'append_{fmt}')):
                if fancy:
                    t1 = f'`{yes}<{module.__name__}.write_{fmt}()>`'
                    t2 = f'`{yes}<{module.__name__}.append_{fmt}()>`'
                    row.append(t1 if hw else t2)
                else:
                    row.append('yes')
            elif rc and wc:
                row.append('..' if fancy else 'no')
        elif what == 'fts':
            if rc and hasattr(module, f'read_fts_{fmt}'):
                if fancy:
                    t1 = f'`{yes}<{module.__name__}.read_fts_{fmt}()>`'
                else:
                    t1 = 'yes'
                row.append(t1)
            elif rc and wc:
                row.append('..' if fancy else 'no')
            if wc and hasattr(module, f'write_fts_{fmt}'):
                if fancy:
                    t1 = f'`{yes}<{module.__name__}.write_fts_{fmt}()>`'
                else:
                    t1 = 'yes'
                row.append(t1)
            elif rc and wc:
                row.append('..' if fancy else 'no')
        if len(row) > 2:
            table.append(row + [module.__doc__.lstrip().splitlines()[0].strip() or '..'])
    if len(table) < 5:  # we need at least 2 rows in table body, add empty row
        table.append(['..'] * (3 + rc + wc))
    table.append(sep)
    # fancy = False
    lens = [max(map(len, col))+fancy for col in zip(*table)]
    text = ('\n' + ws * ' ').join('  '.join(col.strip().ljust(l-(yes in col)*fancy, '=' if col[0] == '=' else ' ')
                                            for l, col in zip(lens, row))
                         for row in table)
    return text


def _insert_format_plugin_table(what, io, ws=None, fancy=False):
    """
    A function to populate the docstring of a function with its plugin table.
    """
    def deco(func):
        try:
            doc = func.__doc__
        except AttributeError:
            return
        try:
            doc = func._original_doc
        except AttributeError:
            pass
        if '{format_table}' in doc:
            func._original_doc = doc
            ws2 = ws
            if ws2 is None:
                ws2 = len(doc) - len(doc.lstrip()) - 1
            func.__doc__ = doc.format(format_table=_create_format_plugin_table(
                what, io, ws=ws2, fancy=fancy))
        return func
    return deco

_ADD_FMT_DOC = """
    .. warning::
        This function should NOT be called directly, it registers via {call}, call this instead.
"""
FULLCALL = {
    'read': '`.read()` and `.iter_()`',
    'read_fts': '`.read_fts()`',
    'write': '`.BioBasket.write()`',
    'write_fts': '`.FeatureList.write()`',
    }

def _add_fmt_doc(call):
    def deco(func):
        try:
            doc = func.__doc__
        except AttributeError:
            doc = None
        try:
            doc = func._original_doc
        except AttributeError:
            pass
        func._original_doc = doc
        func.__doc__ = (doc or '') + _ADD_FMT_DOC.format(call=FULLCALL[call])
        return func
    return deco
