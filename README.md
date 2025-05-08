<img src="https://raw.github.com/rnajena/sugar/logo/sugar_logo.png" alt="logo" width="200">

## A Python framework for bioinformatics
[![build status](https://github.com/rnajena/sugar/workflows/tests/badge.svg)](https://github.com/rnajena/sugar/actions)
[![docs status](https://readthedocs.org/projects/rnajena-sugar/badge/?version=latest)](https://rnajena-sugar.readthedocs.io)
[![codecov](https://codecov.io/gh/rnajena/sugar/branch/master/graph/badge.svg)](https://codecov.io/gh/rnajena/sugar)
[![pypi version](https://img.shields.io/pypi/v/rnajena-sugar.svg)](https://pypi.python.org/pypi/rnajena-sugar)
[![python version](https://img.shields.io/pypi/pyversions/rnajena-sugar.svg)](https://python.org)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.11388074.svg)](https://doi.org/10.5281/zenodo.11388074)

The sugar project aims to provide a lightweight framework to facilitate rapid application development for bioinformatics.

It thus provides classes and functions for working with DNA and RNA sequences, as well as related annotations, and it provides parsers and writers for various file formats using a plugin interface.

### Installation

Use pip, e.g.

```
pip install rnajena-sugar
```

Run tests with the `sugar test` command.

Other options are described in the [documentation](https://rnajena-sugar.readthedocs.io/en/latest/src/tutorial_install.html).

### Getting started

Read about how to get started in the [tutorial](https://rnajena-sugar.readthedocs.io/en/latest/src/tutorial_install.html) section of the documentation.

```python
from sugar import read
seqs = read()  # load example GenBank file
print(seqs)
print(seqs[1].fts)  # show features attached to second sequence
seqs.plot_ftsviewer(show=True)  # plot attached features
aas = seqs['cds'].translate()
print(aas)
aas.write('translated_cds.fasta')
```

### Documentation and Changelog

Documentation can be found at [Read the Docs](https://rnajena-sugar.readthedocs.io).
The detailed changelog is available [here](https://github.com/rnajena/sugar/blob/master/CHANGELOG).

### Contributions

Contributions are welcome! -- e.g. report or fix bugs, discuss or add features, improve the documentation.
See the [CONTRIBUTING.md](https://github.com/rnajena/sugar/blob/master/CONTRIBUTING.md) file for details.
