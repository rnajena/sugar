"""
sugar.data -- Use genetic code and substitution matrices

For reference the IUPAC nucleotide code:

.. code-block:: text

    IUPAC nucleotide code 	Base
    A 	Adenine
    C 	Cytosine
    G 	Guanine
    T (or U) 	Thymine (or Uracil)
    R 	A or G
    Y 	C or T
    S 	G or C
    W 	A or T
    K 	G or T
    M 	A or C
    B 	C or G or T
    D 	A or G or T
    H 	A or C or T
    V 	A or C or G
    N 	any base
    . or - 	gap

And the amino acid codes:

.. code-block:: text

    IUPAC amino acid code 	Three letter code 	Amino acid
    A 	Ala 	Alanine
    C 	Cys 	Cysteine
    D 	Asp 	Aspartic Acid
    E 	Glu 	Glutamic Acid
    F 	Phe 	Phenylalanine
    G 	Gly 	Glycine
    H 	His 	Histidine
    I 	Ile 	Isoleucine
    K 	Lys 	Lysine
    L 	Leu 	Leucine
    M 	Met 	Methionine
    N 	Asn 	Asparagine
    P 	Pro 	Proline
    Q 	Gln 	Glutamine
    R 	Arg 	Arginine
    S 	Ser 	Serine
    T 	Thr 	Threonine
    V 	Val 	Valine
    W 	Trp 	Tryptophan
    Y 	Tyr 	Tyrosine
"""

from functools import lru_cache
from importlib.resources import files
import json
from os.path import exists

# IUPAC nucleotid code
CODES = {'A': 'A', 'C': 'C', 'G': 'G', 'T':'T',
         'R': 'AG', 'Y': 'CT', 'S': 'GC', 'W': 'AT', 'K': 'GT', 'M': 'AC',
         'B': 'CGT', 'D': 'AGT', 'H': 'ACT', 'V': 'ACG', 'N': 'ACGT',
         '.': '.', '-': '-'}


def _submat_files():
    return sorted(f.name for f in files('sugar.data.data_submat').iterdir()
                  if not f.name.startswith('README'))

@lru_cache
def submat(fname):
    """
    Return substitution matrix as a dict of dicts

    >>> from sugar.data import submat
    >>> bl = submat('blosum62')
    >>> bl['A']['A']
    4

    :param fname: One of the following values: ``{}``. Or use your own file.
    """
    if not exists(fname):
        fname2 = str(files('sugar.data.data_submat').joinpath(fname.upper()))
        if not exists(fname2):
            fnames = ', '.join(_submat_files())
            msg = f'No file at {fname} or {fname2}, available matrices: {fnames}'
            raise FileNotFoundError(msg)
        fname = fname2
    with open(fname) as f:
        content = f.read()
    mat = {}
    letters = None
    for line in content.splitlines():
        if line.strip().startswith('#') or line.strip() == '':
            continue
        if letters is None:
            letters = line.split()
        else:
            l1, rest = line.split(maxsplit=1)
            if l1 not in letters:
                from warnings import warn
                warn(f'Letter {l1} not found in table header "{" ".join(letters)}"')
            vals = map(float if '.' in rest else int, rest.split())
            mat[l1] = {l2: v for l2, v in zip(letters, vals)}
    return mat


@lru_cache
def scale_submat(sm, scale=1):
    """
    Return Scaled substition matrix

    The matrix values are divided by the sum of all entries and
    multiplied with the scale factor.

    :param scale: scale factor

    .. warning::
        It is not clear if this function is useful. It might be removed in
        a later version of sugar without further notice.
    """
    s = sum(abs(v) for row in sm.values() for v in row.values())
    for k1 in sm:
        for k2 in sm[k1]:
            sm[k1][k2] = sm[k1][k2] / s * scale
    return sm


@lru_cache
def gcode(tt=1):
    """
    Return a genetic code object

    :param tt: number of the translation table (default: 1)

    >>> from sugar.data import gcode
    >>> gc = gcode()
    >>> gc.tt['TAG']
    '*'
    >>> gc.starts  # doctest: +SKIP
    {'ATG', 'CTG', 'TTG'}
    """
    from sugar import Attr
    fname = files('sugar.data.data_gcode').joinpath('gc.json')
    with open(fname) as f:
        gcs = json.load(f)
    gc = Attr(**gcs[str(tt)])
    gc.ttinv = {k: set(v) for k, v in gc.ttinv.items()}
    gc.starts = set(gc.starts)
    gc.astarts = set(gc.astarts)
    gc.stops = set(gc.stops)
    gc.astops = set(gc.astops)
    return gc


if hasattr(submat, '__doc__'):
    submat.__doc__ = submat.__doc__.format(', '.join(_submat_files()))
