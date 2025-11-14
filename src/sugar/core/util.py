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


# def _add_doc_from_other(ofunc):
#     def decorator(func):
#         func.__doc__ = ofunc.__doc__
#         return func
#     return decorator


class SugarDeprecationWarning(UserWarning):
    pass


def deprecated(msg):
    """
    Decorator factory for deprecated functions

    :param msg: The msg will be emitted as warning, when the function is called.
        Additionally, a warning is displayed in the doc string,
        if the term 'deprecated' is not found in the doc string.

    Use it the following way:

    ```
    from sugar.core.util import deprecated
    @deprecated("old_func is deprecated, use new_func instead")
    def old_func(...):
        ...
    ```
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


def _func_kws(func):
    from inspect import signature
    params = signature(func).parameters
    return [pname for pname, p in params.items() if p.default is not p.empty]


def _get_kws_for_func(kw, func):
    return {k: v for k, v in kw.items() if k in _func_kws(func)}


def _pop_kws_for_func(kw, func):
    kw2 = _get_kws_for_func(kw, func)
    for k in kw2:
        kw.pop(k)
    return kw2
