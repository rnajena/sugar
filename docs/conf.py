# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information

import sys
import os
from pathlib import Path
from sphinx.ext import apidoc

os.environ['SPHINX_BUILD'] = '1'
rootpath = Path(__file__).parent.resolve()


### First run sphinx-apidoc to generate the documentation
# usage: sphinx-apidoc [OPTIONS] -o <OUTPUT_PATH> <MODULE_PATH> [EXCLUDE_PATTERN, ...]
#   module_path           path to module to document
#   exclude_pattern       fnmatch-style file and/or directory patterns to exclude from generation
#   -o DESTDIR            directory to place all output
#   -d MAXDEPTH           maximum depth of submodules to show in the TOC (default: 4)
#   -f, --force           overwrite existing files
#   -e, --separate        put documentation for each module on its own page
#   -P, --private         include "_private" modules
#   -T, --no-toc          don't create a table of contents file
#   -E, --no-headings     don't create headings for the module/package packages (e.g. when the docstrings already contain them)
#   -M, --module-first    put module documentation before submodule documentation

spath = rootpath.parent / 'sugar'
src = rootpath / 'src_api'
exclude = ['tests', 'scripts.py']
call = f"sphinx-apidoc -d 3 -f -M -e -T -P -o {src} {spath} " + ' '.join(str(spath / ex) for ex in exclude)
print(call)
apidoc.main(call.split()[1:])


### Now modify some stuff in the generated apidoc
index_orig = src / 'sugar.rst'
# index = src / 'index.rst'
# print(f'move and modify {index_orig} -> {index}')
# data = index_orig.read_text()
# index.write_text(data)
index_orig.unlink()
import shutil
shutil.copy(rootpath / 'src/index.rst', src)
#shutil.copy(rootpath / 'src/sidebar_toc.rst', src)


### Now configure sphinx-build
sys.path.insert(0, rootpath)
import sugar

project = 'sugar'
copyright = '2024, Tom Eulenfeld'
author = 'Tom Eulenfeld'
release = sugar.__version__

### General configuration
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration
root_doc = 'index'
default_role = 'py:obj'
nitpicky = True  # Show warnings for unreferenced targets
templates_path = ['_templates']
exclude_patterns = ['_build', 'Thumbs.db', '.DS_Store']
rst_epilog = """
.. _BinarySeachFile: https://github.com/trichter/binarysearchfile
"""

# Add any Sphinx extension module names here
extensions = ['sphinx.ext.autodoc',
              'sphinx.ext.doctest',
              'sphinx.ext.intersphinx',
              'sphinx.ext.autosummary',
              'sphinx.ext.viewcode',
              'sphinx.ext.napoleon',
              'sphinx_rtd_theme',
              ]

autodoc_default_options = {
    'members': True,
    'private-members': False,
    'show-inheritance': True,
}

def autodoc_skip_member(app, what, name, obj, skip, options):
    return
    print(name, what, skip)
    if name == '_io':
        return False

def setup(app):
    app.connect('autodoc-skip-member', autodoc_skip_member)

autosummary_generate = True



### Options for HTML output
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output
def download_logo():
    import urllib.request
    logo = rootpath / '_static/sugar_logo.png'
    url = 'https://raw.githubusercontent.com/rnajena/sugar/logo/sugar_logo.png'
    urllib.request.urlretrieve(url, logo)
    return logo.name

html_theme = 'alabaster'
#html_theme = 'sphinx_rtd_theme'
html_static_path = ['_static']
html_show_sphinx  = False
html_show_sourcelink = True
html_sidebars = {
    '**': [
        'about.html',
        'navigation.html',
        'searchbox.html'
    ]
}

html_sidebars = { '**': ['about.html', 'searchbox.html', 'mynavigation.html']}#, 'globaltoc.html', 'localtoc.html'] }
html_theme_options = {
    'logo': download_logo(),
    'description': 'S*R for your RNA',
    'description_font_style': 'text-align: center;',
    'page_width': '1240px',
    #'sidebar_width': '220px',
    'fixed_sidebar': True,
    'sidebar_collapse': True,
    'extra_nav_links': {'Homepage': 'https://github.com/rnajena/sugar'},
    # 'show_related': True,
    # 'show_relbars': True,
    # 'navigation_depth': 4,
    # 'sidebar_maxdepth': -1,
}


intersphinx_mapping = {'python': ('https://docs.python.org/3/', None)
                       }
