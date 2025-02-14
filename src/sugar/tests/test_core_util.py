# (C) 2024, Tom Eulenfeld, MIT license
from sugar.core.util import _get_kws_for_func, _pop_kws_for_func


def _f1(a, b=5, *, c=6):
    return


def test_kw_for_func():
    kw = {'a': 1, 'b': 2, 'c': 3}
    kw2 = _get_kws_for_func(kw, _f1)
    assert len(kw2) == 2
    assert len(kw) == 3
    kw3 = _pop_kws_for_func(kw, _f1)
    assert len(kw2) == 2
    assert len(kw) == 1
    assert kw2 == kw3
