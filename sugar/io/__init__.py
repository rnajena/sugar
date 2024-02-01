"""
sugar.io

Input/Output support for different sequence file formats

IO support is discovered via plugins, you can simply install your own plugins.

i.e. custom_sugar_io.py
Some of the following features can be used.

EXT  # list of extensions for autodetection of format upon writing
is_format(f)  # function for autodetection of format upon reading
read()  # function for reading th full file
iter_()  # function, yielding sequence while reading the file
write()  # function for writing all sequences
append()  # function for appending a single sequence to a file

# You need to define read or iter_ (or both) for reading support
# You need to define write or append (or both) for writing support
# See for example sugar.io.fasta
"""

from sugar.io.main import detect, detect_ext, iter_, read, write
