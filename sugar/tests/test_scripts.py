# (C) 2023, Tom Eulenfeld, MIT license
from subprocess import check_output

from sugar import read
from sugar.tests.util import tempfilename



def test_print_convert():
    with tempfilename() as fname:
        read().write(fname, 'fasta')
        out = check_output(f'sugar print {fname}'.split())
        assert b'ACCTG' in out
        with tempfilename() as fname2:
            check_output(f'sugar convert {fname} -o {fname2} -fo sjson'.split())
