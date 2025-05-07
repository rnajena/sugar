How to create a new IO plugin
=============================

.. rubric:: Create a sequence plugin in the sugar package

Let's assume you want to create an IO plugin for the sequence format *fancy*.

First, take a look at other plugins in the ``sugar/_io`` folder.

1. Fork the repository to your account and create a new branch
2. Create a new module or package inside the ``sugar/_io`` folder, e.g. ``fancymodule.py``.
3. Use the following template to provide read and/or write functionality for your format::

    """
    My Fancy Seq Plugin

    This is an example sequence file format. The layout is as follows:
    #MyFancySeqFormat
    seq1 AAATTGGGCCC
    seq2 ATGGCT
    """

    # To add read support you must define either an iter_fancy or read_fancy function, or both.
    # To add write support you must define either an append_fancy or write_fancy function, or both.

    from sugar import BioBasket, BioSeq
    from sugar._io.util import _add_fmt_doc

    # Use the following flag to indicate, that your file format is binary rather than text-based,
    # the passed file handlers will be opened in binary mode.
    #binary_fmt = True

    # optional, filename extensions for automatic detection of file format
    # when writing
    filename_extensions = ['fancy']

    def is_fancy(f, **kw):
        """
        Function is optional, used for auto-detection of format when reading

        It should return True if the format is detected,
        otherwise it may raise any exception or return False.
        """
        content = f.read(50)
        return content.strip().lower().startswith('#myfancyseqformat')

    # The function decorators are used to automatically add a warning
    # to the docstring, that this function should be called via the main
    # iter_ or read functions.
    @_add_fmt_doc('read')
    def iter_fancy(f, optional_argument=None):
        """
        The iter_fancy function expects a file handler and has to yield BioSeq objects.

        You can define optional arguments.
        """
        for line in f:
            if line.strip() != '' and not line.startswith('#'):
                seqid, data = line.split()
                yield BioSeq(data, id=seqid)

    @_add_fmt_doc('read')
    def read_fancy(f, **kw):
        """
        The read_fancy function expects a file handler and has to return a BioBasket object
        """
        # We are lazy here and reuse iter_fancy
        return BioBasket(list(iter_fancy(f, **kw)))

    @_add_fmt_doc('write')
    def append_fancy(seq, f, **kw):
        """
        Write a single seq to file handler
        """
        f.write(f'{seq.id} {seq.data}\n')

    @_add_fmt_doc('write')
    def write_fancy(seqs, f, **kw):
        """
        Write a BioBasket object to file handler
        """
        f.write('#MyFancySeqFormat 3.14159\n')
        for seq in seqs:
            # be lazy again
            append_fancy(seq, f, **kw)

4. Add your plugin to the ``FMTS`` list in ``sugar/_io/util.py``
5. Register the plugin in the ``pyproject.toml`` file::

    [project.entry-points."sugar.io"]
    fancy = "sugar._io.fancymodule"

6. Write some tests in a new file ``sugar/tests/test_io_fancy.py``.
7. Re-Install your branch of sugar and check that everything is working::

    from sugar import read
    seqs = read('example.fancy')
    print(seqs)

8. Run your tests.
9. Create a pull request to get your plugin into the main repository.

.. rubric:: Create a sequence plugin that can be used with sugar, but is in an external package

Create your own package by following only step 3 above. Register the plugin in the ``pyproject.toml`` of your own project::

    [project.entry-points."sugar.io"]
    fancy = "myfancypackage.fancymodule"


When your package is installed you can still read seq files using the commands in point 7 above.

.. rubric:: Create a features plugin in the sugar package or in an external package

This is analogous to the sequence plugin. It can even be located in the same file as the sequence plugin.

Instead, use the following variables and function definitions::

    from sugar.core.fts import FeatureList, Location, Feature
    from sugar._io.util import _add_fmt_doc

    #binary_fmt_fts
    filename_extensions_fts
    def is_fts_fancy(f, **kw):
        ...

    @_add_fmt_doc('read_fts')
    def read_fts_fancy(f, **kw):
        ...
        return FeatureList(fts)

    @_add_fmt_doc('write_fts')
    def write_fts_fancy(fts, f, **kw):
        f.write(...)

.. rubric:: Format test files

Test files for IO plugins can be placed in the ``src/sugar/tests/data`` folder.
These files can be read using the ``!data`` magic.
For example, ``read('!data/example.fasta')`` will read the file ``example.fasta`` in the above folder.
