# (C) 2024, Tom Eulenfeld, MIT license
import pytest

from sugar import read
from  sugar.core.translate import translate


def test_translate():
    s = 'CCC-AT-GAT-NCC--CCCCT---ANT-A--GGGN'
    aas = 'P-MX-PP-X-*G'

    assert translate(s, complete=True, warn=True) == aas
    with pytest.warns(UserWarning, match='might be a stop'):
        assert translate(s[3:-3], warn=True) == aas[1:-2]
    assert translate(s[3:-3], warn=False) == aas[1:-2]
    with pytest.warns(UserWarning, match='not a start|might be a stop'):
        assert translate(s[8:-3], warn=True) == aas[3:-2]
    with pytest.raises(ValueError, match='not a start'):
        translate(s)


def test_translate_method():
    seqs = read()['cds'].translate()
    assert str(seqs[0]) == seqs[0].fts.get('cds').meta._genbank.translation
