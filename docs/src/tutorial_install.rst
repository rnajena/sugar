Installation
============

.. toctree::
    :hidden:

    self
    Sequences <tutorial_seqs>
    Features <tutorial_fts>
    Plotting <tutorial_imaging>
    More bits <tutorial_misc>

Sugar is hosted on PyPi and can be installed with ``pip``::

    pip install rnajena-sugar

To upgrade an existing installation, use the ``-U`` flag, e.g. ``pip install rnajena-sugar -U``.

Note, that the package name is different from the import name ``sugar``.


.. rubric:: Inside a conda environment

If you want to install sugar in a dedicated conda environment,
you may want to install most of its dependencies using conda first, e.g. ::

    conda install -c conda-forge matplotlib pandas platformdirs pytest requests seaborn
    pip install rnajena-sugar


.. rubric:: Development version

If you want to install the latest bleeding-edge version you can clone the GitHub repository
and install the code with ``pip``, one way to do this, is::

    pip install git+https://github.com/rnajena/sugar.git


.. rubric:: Editable install

If you want to contribute code or want to hack into sugar quickly,
an editable install is generally preferable::

    git clone https://github.com/rnajena/sugar.git
    pip install -e sugar


.. rubric:: Minimal install

For advanced users, it is possible to get a base installation of sugar using pip's ``--no-deps`` flag.
A subset of dependencies can be installed manually, depending on the functionality needed.


.. rubric:: Running the test suite

After installation, the test suite can be run with::

    sugar test

You may want to check out its options with ``sugar test -h`` and ``pytest -h``.
Note that tests that require an Internet connection are skipped by default,
this behavior can be turned off with the ``--web`` option.
Also, depending on the installation and platform, some tests may be skipped or have an expected failure,
to see details about the reasons please use the ``--verbose`` flag.
