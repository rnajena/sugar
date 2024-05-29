# (C) 2024, Tom Eulenfeld, MIT license

from contextlib import redirect_stdout
import io


def doctest_module(m):
    from doctest import testmod
    raised = False
    try:
        testmod(m, raise_on_error=True)
    except Exception:
        raised = True
    if raised:
        report = io.StringIO()
        with redirect_stdout(report):
            testmod(m, report=True)
        assert report.getvalue() == ''


def test_docs_data():
    from sugar import data
    doctest_module(data)


def test_docs_core():
    from sugar import core
    doctest_module(core)

def test_docs_io():
    from sugar import _io
    doctest_module(_io)
    doctest_module(_io.main)
