# (C) 2024, Tom Eulenfeld, MIT license
"""
Helper functions for _core module
"""


_ADD_INPLACE_DOC = """
        .. note::
            This function works in place and modifies the data.
            If you want to keep the original data use the `copy()` method first.
"""

def _add_inplace_doc(func):
    try:
        doc = func.__doc__
    except AttributeError:
        doc = None
    try:
        doc = func._original_doc
    except AttributeError:
        pass
    func._original_doc = doc
    func.__doc__ = (doc or '') + _ADD_INPLACE_DOC
    return func

