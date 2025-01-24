Plotting alignments
===================

Alignments can be plotted using the `.BioBasket.plot_alignment()` method.
Some code examples follow:

.. These figures are generated manually with the test suite

>>> from sugar import read
>>> seqs = read('https://osf.io/download/j2wyv')
>>> seqs.plot_alignment(show=True, figsize=(10, 4))

.. image:: ../_static/ali1.png
    :width: 60%

>>> seqs[:, 70:120].plot_alignment(fname='ali.pdf', color=None, figsize=(10,8),
...                                symbols=True, aspect=2, alpha=0.5)

.. image:: ../_static/ali2.png
    :width: 40%

The following example shows an alignment of amino acids
using the ``'flower'`` color scheme, a color scheme that
displays amino acids with a high BLOSSUM62 score in similar colors.
Color schemes originate from Biotite (Gecos), Jalview and ClustalX,
see `here`_ for an overview.
The available color schemes are listed in `.get_color_scheme()`.

>>> seqs2 = seqs[:5, :150].copy()
>>> seqs2.translate(complete=True).plot_alignment(
...     show=True, color='flower', figsize=(10,8),  symbols=True,
...     aspect=2, alpha=0.5, edgecolors='w')

.. image:: ../_static/ali3.png
    :width: 40%

The plotting function has many options,
including marking of plotting feature regions with different colors.
Multi-line plots are not supported.
If you need these, consider converting the sequences to a biotite ``Alignment`` object via
`seqs.tobiotite(msa=True) <.BioBasket.tobiotite>`
and using
`Biotite's plotting capabilities <https://www.biotite-python.org/latest/examples/gallery/sequence>`_.

.. _here: https://www.biotite-python.org/latest/examples/gallery/sequence/misc/color_schemes_protein.html