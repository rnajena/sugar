"""
sugar._io -- Input/Output support for different sequence and feature file formats

IO support is discovered via plugins, you can simply install your own plugins.
"""

from sugar._io.main import (detect, detect_ext,
                           iter_, read, write, read_fts, write_fts)
