<img src="https://raw.github.com/rnajena/sugar/logo/sugar_logo.png" alt="logo" width="200">

## A Python framework for bioinformatics
[![build status](https://github.com/rnajena/sugar/workflows/tests/badge.svg)](https://github.com/rnajena/sugar/actions)
[![docs status](https://readthedocs.org/projects/rnajena-sugar/badge/?version=latest)](https://rnajena-sugar.readthedocs.io)
[![codecov](https://codecov.io/gh/rnajena/sugar/branch/master/graph/badge.svg)](https://codecov.io/gh/rnajena/sugar)
[![pypi version](https://img.shields.io/pypi/v/rnajena-sugar.svg)](https://pypi.python.org/pypi/rnajena-sugar)
[![python version](https://img.shields.io/pypi/pyversions/rnajena-sugar.svg)](https://python.org)
[![DOI](https://joss.theoj.org/papers/10.21105/joss.08122/status.svg)](https://doi.org/10.21105/joss.08122)

The sugar project aims to provide a lightweight framework to facilitate rapid application development for bioinformatics.

It thus provides classes and functions for working with DNA and RNA sequences, as well as related annotations, and it provides parsers and writers for various file formats using a plugin interface.

### Installation

Sugar requires a Python version `>=3.11`.
Use pip to install the package, e.g.

```
pip install rnajena-sugar
```

Other installation options are described in the [documentation](https://rnajena-sugar.readthedocs.io/en/latest/src/tutorial_install.html).

Run tests with the `sugar test` command. Note that tests that require an Internet connection are skipped by default, this behavior can be turned off with the `--web` option. Also, depending on the installation and platform, some tests may be skipped or have an expected failure, to see details about the reasons please use the `--verbose` flag.

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

### Citation

If you found this package useful, please consider citing it.

Eulenfeld T (2025),
Sugar: A Python framework for bioinformatics,
*Journal of Open Source Software*, 10(111), 8122,
doi:[10.21105/joss.08122](https://doi.org/10.21105/joss.08122)
