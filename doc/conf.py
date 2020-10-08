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
    # Make sure to build Doxygen XML
    subprocess.call('cd ..; doxygen', shell=True)
else:
    # Make sure to override the default Sphinx theme when not on RTD
    html_theme = 'sphinx_rtd_theme'

sys.path.insert(0, os.path.abspath('../'))

extensions = ['breathe']
breathe_projects = {'vg': 'doxygen/xml' }
breathe_default_project = 'vg'

source_suffix = '.rst'
master_doc = 'index'
project = 'vg'
copyright = 'vgteam'
exclude_patterns = ['_build']
pygments_style = 'sphinx'

html_static_path = ['static']

def setup(app):
    # Fix up included file list for files so they don't all run together
    app.add_stylesheet('custom.css')
