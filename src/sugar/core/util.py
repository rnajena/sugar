# (C) 2024, Tom Eulenfeld, MIT license
"""
Helper functions for _core module
"""
import functools
import warnings


_ADD_INPLACE_DOC = """
        .. note::
            This function works in place and modifies the data.
            If you want to keep the original data,
            use the `copy()` method first.
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


class SugarDeprecationWarning(UserWarning):
    pass


def deprecated(msg):
    """
    Decorator factory for deprecated functions

    :param msg: The msg will be emitted as warning, when the function is called.
        Additionally, a warning is displayed in the doc string,
        if the term 'deprecated' is not found in the doc string.
    """
    def _deprecated(func):
        @functools.wraps(func)
        def dfunc(*args, **kw):
            warnings.warn(msg, category=SugarDeprecationWarning, stacklevel=2)
            return func(*args, **kw)
        if hasattr(dfunc, '__doc__') and 'deprecated' not in dfunc.__doc__.lower():
            ws = (len(dfunc.__doc__.lstrip('\n')) - len(dfunc.__doc__.lstrip())) * ' '
            dfunc.__doc__ = f'{ws}.. warning::\n{ws}\n{ws}    {msg}\n' + dfunc.__doc__
        return dfunc
    return _deprecated
