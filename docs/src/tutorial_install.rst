Installation
============

.. toctree::
    :hidden:

    self
    Sequences <tutorial_seqs>
    Features <tutorial_fts>
    Plotting <tutorial_imaging>
    tutorial_misc

Sugar is hosted on PyPi and can be installed with ``pip``::

    pip install rnajena-sugar

Note, that the package name is different from the import name ``sugar``.


.. rubric:: Inside a conda environment

If you want to install sugar in a dedicated conda environment,
you might want to install most of its dependencies using conda beforehand, e.g. ::

    conda install -c conda-forge matplotlib pandas platformdirs pytest requests seaborn
    pip install rnajena-sugar


.. rubric:: Development version

If you want to install the latest bleeding-edge version you can clone the GitHub repository
and install the code with ``pip``, one way to do it, is::

    pip install git+https://github.com/rnajena/sugar.git


.. rubric:: Editable install

If you want to contribute code or quickly want to hack into sugar,
an editable install is generally preferable::

    git clone https://github.com/rnajena/sugar.git
    pip install -e sugar


.. rubric:: Running the test suite

After installation the test suite can be run with::

    sugar test

You might want to check out its options with ``sugar test -h`` and ``pytest -h``.
