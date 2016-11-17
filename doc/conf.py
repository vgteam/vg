# -*- coding: utf-8 -*-
#

# Sphinx configuration
# Cribbed from <https://github.com/ondrejsika/sphinx-autodoc-example> and
# <http://breathe.readthedocs.io/en/latest/readthedocs.html>

import sys
import os
import subprocess

read_the_docs_build = os.environ.get('READTHEDOCS', None) == 'True'

if read_the_docs_build:

    subprocess.call('cd ..; doxygen', shell=True)

sys.path.insert(0, os.path.abspath('../'))

extensions = ['breathe']
breathe_projects = {'vg': 'doxygen/xml' }
breathe_default_project = 'vg'

source_suffix = '.rst'
master_doc = 'index'
project = u'vg'
copyright = u'vgteam'
exclude_patterns = ['_build']
pygments_style = 'sphinx'
html_theme = 'default'
autoclass_content = 'both'
html_theme = 'sphinx_rtd_theme'
html_static_path = ['static']

def setup(app):
    # Fix up included file list for files so they don't all run together
    app.add_stylesheet('custom.css')
