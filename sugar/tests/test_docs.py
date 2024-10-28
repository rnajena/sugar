# (C) 2024, Tom Eulenfeld, MIT license

from contextlib import redirect_stdout
import io
from sugar import read, read_fts


def doctest_module(m):
    from doctest import testmod, ELLIPSIS
    raised = False
    flags = ELLIPSIS
    globs = {'read': read, 'read_fts': read_fts}
    try:
        testmod(m, raise_on_error=True, optionflags=flags, extraglobs=globs)
    except Exception:
        raised = True
    if raised:
        report = io.StringIO()
        with redirect_stdout(report):
            testmod(m, optionflags=flags, report=True, extraglobs=globs)
        assert report.getvalue() == ''


def test_docs_data():
    from sugar import data
    doctest_module(data)

def test_docs_core():
    from sugar import core
    doctest_module(core)

def test_docs_core_seq():
    from sugar.core import seq
    doctest_module(seq)

def test_docs_core_fts():
    from sugar.core import fts
    doctest_module(fts)

def test_docs_core_meta():
    from sugar.core import meta
    doctest_module(meta)

def test_docs_io():
    from sugar import _io
    doctest_module(_io)
    doctest_module(_io.main)
