# Configuration file for the Sphinx documentation builder.
# Build doc locally with
# sphinx-build -anE -c . src _build

import os

os.environ['SPHINX_BUILD'] = '1'


def parse_version():
    from pathlib import Path
    import re
    init_path = Path(__file__).parent.parent / 'sugar/__init__.py'
    with open(init_path) as f:
        content = f.read()
    regex = r"""__version__\s*=\s*['"]([^\s]+)['"]"""
    match = re.search(regex, content)
    if match:
        return match.group(1)


# def parse_imports():
#     from pathlib import Path
#     import re
#     root_path = Path(__file__).parent.parent / 'sugar/'
#     regex = r'^\s*import\s+(\w+)(?:.\w+)*|^\s*from\s+(\w+)(?:.\w+)*\s+import'
#     modules = set()
#     for fname in root_path.glob('**/*.py'):
#         with open(fname) as f:
#             content = f.read()
#         for match in re.finditer(regex, content, flags=re.MULTILINE):
#             modules |= set(match.groups())
#     modules -= {None, 'anchorna'}
#     return sorted(modules)


def download_logo():
    from pathlib import Path
    import urllib.request
    logo = Path(__file__).parent / '_static/sugar_logo.png'
    url = 'https://raw.githubusercontent.com/rnajena/sugar/logo/sugar_logo.png'
    urllib.request.urlretrieve(url, logo)
    return logo.name


project = 'sugar'
copyright = '2024, Tom Eulenfeld'
author = 'Tom Eulenfeld'
release = parse_version()

default_role = 'py:obj'
templates_path = ['_templates']
exclude_patterns = ['_build']
rst_epilog = """
.. _BinarySeachFile: https://github.com/trichter/binarysearchfile
.. _BLAST: https://blast.ncbi.nlm.nih.gov
.. _FASTA: https://en.wikipedia.org/wiki/FASTA_format
.. _Genbank: https://www.insdc.org/submitting-standards/feature-table
.. _GFF: https://github.com/The-Sequence-Ontology/Specifications/blob/master/gff3.md
.. _Infernal: http://eddylab.org/infernal/
.. _MMseqs2: https://github.com/soedinglab/MMseqs2
.. _Stockholm: https://en.wikipedia.org/wiki/Stockholm_format
"""

extensions = ['sphinx.ext.autodoc',
              'sphinx.ext.doctest',
              'sphinx.ext.intersphinx',
              'sphinx.ext.autosummary',
              'sphinx.ext.viewcode',
              'sphinx.ext.napoleon',
              ]

# autodoc_mock_imports = parse_imports()
# print(f'set {autodoc_mock_imports=}')
autodoc_default_options = {
    'members': True,
    'undoc-members': True,
    'private-members': False,
    'show-inheritance': True,
}

html_theme = 'furo'
html_static_path = ['_static']
html_title = f'S*R for your RNA <br>v{release} docs'
html_logo = '_static/' + download_logo()
html_show_sphinx  = True
html_theme_options = {
    'footer_icons' : [],
    "source_edit_link": " https://github.com/rnajena/sugar/edit/master/docs/src/{filename}",
    "source_view_link": "_sources/{filename}.txt",
    #'top_of_page_buttons': [],
}


intersphinx_mapping = {
    'python': ('https://docs.python.org/3/', None)
    }
