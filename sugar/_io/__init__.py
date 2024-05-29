"""
sugar._io -- Input/Output support for different sequence and feature file formats

IO support is discovered via plugins, you can simply install your own plugins.

..
    TODO update and move to wiki
    i.e. custom_sugar_io.py
    Some of the following features can be used.
    EXT  # list of extensions for autodetection of format upon writing
    is_format(f)  # function for autodetection of format upon reading
    read()  # function for reading th full file
    iter\_()  # function, yielding sequence while reading the file
    write()  # function for writing all sequences
    append()  # function for appending a single sequence to a file
    # You need to define read or iter\_ (or both) for reading support
    # You need to define write or append (or both) for writing support
    # See for example sugar.io.fasta
"""

from sugar._io.main import (detect, detect_ext,
                           iter_, read, write, read_fts, write_fts)
