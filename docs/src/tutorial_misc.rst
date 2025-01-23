What else
=========

.. rubric:: Adaptors

Sugar provides adaptors to convert to the corresponding sequence objects
in the Biopython and Biotite libraries, and vice versa.

.. runblock:: pycon

    >>> from sugar import BioBasket, read
    >>> seqs = read()
    >>> print(seqs)
    >>> bios = seqs.tobiopython()
    >>> print(bios[0])
    >>> print(BioBasket.frombiopython(bios))

.. runblock:: pycon

    >>> from sugar import BioBasket, read
    >>> seqs = read()
    >>> tites = seqs.tobiotite()
    >>> print(repr(tites[0])[:80])
    >>> print(BioBasket.frombiotite(tites))  # Biotite does not use ids

.. rubric:: Indexing of FASTA files

sugar provides an indexing tool to quickly retrieve
sequences or subsequences from large FASTA files.
Use it via the `.FastaIndex` class or the ``sugar index`` command.

The following example uses the ``sugar index create`` and ``add`` commands to create
the index and queries the index in python. Both tasks can be performed
either on the command line or with Python code. ::

    sugar index create index.db
    sugar index add *.fasta

::

    >>> from sugar import FastaIndex
    >>> index = FastaIndex('index.db')
    >>> print(index)  # display information about index
    >>> seq = index.get_seq('NC_081844.1')

The Fasta index either uses a binary search file or
a data base via Pythons ``dbm`` module.
Depending on the use case, one or the other option might be preferred,
usually you want the binary search file.

.. rubric:: Downloading sequences from NCBI

The `.Entrez` class can be used to fetch sequences
from the NCBI online database::

    >>> from sugar.web import Entrez
    >>> client = Entrez()
    >>> seqs = client.get_basket(['AF522874', 'NC_077015.1'])
    >>> print(seqs)
    2 seqs in basket
    AF522874    19k  CGGACACACAAAAAGAAAAAAGGTTTTTTAAGACTTTTTGTGTGCGAG...  GC:40.63%
    NC_077015   12k  TGCATAACCCTGATTGTAATTGGCTGGGTTATGCATGTGAGAACGCAA...  GC:43.03%
    customize output with BioBasket.tostr() method

Use the ``ENTREZ_PATH`` environment variable or
the ``path`` option of `.Entrez` or its methods
for a caching of sequence files on the disk.
An Entrez API key can be used via the ``ENTREZ_API_KEY`` environment variable
or the ``api_key`` option of `.Entrez`.


.. rubric:: Miscellaneous

For command line junkies
    The command line interface of sugar can be used for common tasks.
    Call ``sugar -h`` for an overview of available commands.
The sugar.data package
    sugar also provides access to codon translation tables and substitution
    matrices (i.e. BLOSUM62) within the `sugar.data` package.
Find open reading frames
    The `~.BioBasket.find_orfs()` method can be used to find open reading frames,
    see the advanced example in the :doc:`Sequences Tutorial <tutorial_seqs>`.
Attributes behave like dictionaries
    Metadata and attributes in sugar are instances of `.Attr`,
    therefore attributes can be used as dictionary keys,
    e.g. ``seq.meta.id`` and ``seq.meta['id']`` have the same effect.
Shortcuts are convenient
    Using the shortcuts ``seq.fts``, ``seq.id``, ``ft.seqid``, e.t.c., has the advantage,
    that additional checks are performed,
    e.g. assigning ``seq.fts = my_fts`` automatically converts ``my_fts`` to a `.FeatureList`,
    accessing ``seq.id`` returns an emtpy string if no ``id`` attribute is present in the metadata.
    ``ft.loc`` is a shortcut for ``ft.locs[0]``.
Overloading of operators
    Adding two sequences ``seq1 + seq2``, concatenates the data.
    Adding two sequence lists ``seqs1 + seqs2``, concatenates the two lists.
    To concatenate sequences inside a list use the `~.BioBasket.merge()` method.
    Adding two feature lists ``fts1 + fts2`` also works as expected.
