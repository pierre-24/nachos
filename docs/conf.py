# Configuration file for the Sphinx documentation builder.
#
# This file only contains a selection of the most common options. For a full
# list see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Path setup --------------------------------------------------------------

# If extensions (or modules to document with autodoc) are in another directory,
# add these directories to sys.path here. If the directory is relative to the
# documentation root, use os.path.abspath to make it absolute, like shown here.
#
import os
import sys
sys.path.insert(0, os.path.abspath('..'))

import nachos

# -- Project information -----------------------------------------------------

project = nachos.__name__
copyright = '2026, {} (University of Namur)'.format(nachos.__author__)
author = nachos.__author__
# The short X.Y version.
version = nachos.__version__
# The full version, including alpha/beta/rc tags.
release = nachos.__version__

# -- General configuration ---------------------------------------------------

# Add any Sphinx extension module names here, as strings. They can be
# extensions coming with Sphinx (named 'sphinx.ext.*') or your custom
# ones.
extensions = [
    'sphinx.ext.autodoc',
    'sphinx.ext.mathjax',
    'sphinx.ext.napoleon',
    'sphinx.ext.githubpages',
    'sphinxcontrib.autoprogram'
]

# Add any paths that contain templates here, relative to this directory.
templates_path = ['_templates']

# List of patterns, relative to source directory, that match files and
# directories to ignore when looking for source files.
# This pattern also affects html_static_path and html_extra_path.
exclude_patterns = ['_build', 'Thumbs.db', '.DS_Store']


# -- Options for HTML output -------------------------------------------------

# Napoleon settings
napoleon_google_docstring = True
napoleon_include_init_with_doc = True

# The theme to use for HTML and HTML Help pages.  See the documentation for
# a list of builtin themes.
#
html_theme = 'shibuya'

html_theme_options = {
  "github_url": "https://github.com/pierre-24/nachos"
}

html_context = {
    "source_type": "github",
    "source_user": "pierre-24",
    "source_repo": "nachos",
    "source_edit_template": "https://github.com/pierre-24/nachos/blob/dev/docs/{0}",
}

html_sidebars = {
  "**": [
    "sidebars/localtoc.html",
    "sidebars/repo-stats.html",
    "sidebars/edit-this-page.html",
  ]
}