# -*- coding: utf-8 -*-
#
# Configuration file for the Sphinx documentation builder.
#
# This file does only contain a selection of the most common options. For a
# full list see the documentation:
# http://www.sphinx-doc.org/en/master/config

# -- Path setup --------------------------------------------------------------

# If extensions (or modules to document with autodoc) are in another directory,
# add these directories to sys.path here. If the directory is relative to the
# documentation root, use os.path.abspath to make it absolute, like shown here.
#
import os
import sys
from importlib.metadata import PackageNotFoundError
from importlib.metadata import version as _get_version

sys.path.insert(0, os.path.abspath("../.."))

# -- Project information -----------------------------------------------------

project = "ShapePipe"
copyright = "2022, CosmoStat"
author = "CosmoStat"

# The version is read from the installed package metadata, so the docs never
# carry a hand-maintained (and inevitably stale) version string.
try:
    # The full version, including any alpha/beta/rc tags.
    release = _get_version("shapepipe")
except PackageNotFoundError:
    release = "0.0.0"
# The short X.Y version.
version = ".".join(release.split(".")[:2])


# -- General configuration ---------------------------------------------------

# If your documentation needs a minimal Sphinx version, state it here.
#
# needs_sphinx >= '3.3'

# Add any Sphinx extension module names here, as strings. They can be
# extensions coming with Sphinx (named 'sphinx.ext.*') or your custom
# ones.
extensions = [
    "sphinx.ext.autodoc",
    "sphinx.ext.autosummary",
    "sphinx.ext.coverage",
    "sphinx.ext.doctest",
    "sphinx.ext.ifconfig",
    "sphinx.ext.intersphinx",
    "sphinx.ext.mathjax",
    "sphinx.ext.napoleon",
    "sphinx.ext.todo",
    "sphinx.ext.viewcode",
    "sphinxcontrib.bibtex",
    "myst_parser",
    "numpydoc",
]

# Include module names for objects
add_module_names = False

# Include init in class documentation.
autoclass_content = "class"

# Audodoc options
autodoc_default_options = {
    "member-order": "bysource",
    "private-members": True,
    "special-members": "__call__",
    "show-inheritance": True,
    "exclude-members": "urand,randint",
}

# Generate summaries
autosummary_generate = True

# Suppress class members in toctree.
numpydoc_show_class_members = False

# The suffix(es) of source filenames.
# You can specify multiple suffix as a list of string:
#
source_suffix = [".rst", ".md"]

# The master toctree document.
master_doc = "index"

# If true, sectionauthor and moduleauthor directives will be shown in the
# output. They are ignored by default.
show_authors = True

# The name of the Pygments (syntax highlighting) style to use.
pygments_style = "default"

# If true, `todo` and `todoList` produce output, else they produce nothing.
todo_include_todos = True

# -- Options for Markdown files ----------------------------------------------

myst_enable_extensions = ["html_image"]

# Generate slug anchors for h1-h3 so in-page links like
# `[Mask images](#mask-images)` resolve instead of warning.
myst_heading_anchors = 3

# -- Options for HTML output -------------------------------------------------

# The theme to use for HTML and HTML Help pages.  See the documentation for
# a list of builtin themes.
#
html_theme = "sphinx_book_theme"

# Theme options are theme-specific and customize the look and feel of a theme
# further.  For a list of options available for each theme, see the
# documentation.
#
html_theme_options = {
    "repository_url": "https://github.com/CosmoStat/shapepipe",
    "use_repository_button": True,
    "use_issues_button": True,
    "use_download_button": False,
    "use_fullscreen_button": False,
    "use_edit_page_button": True,
    "path_to_docs": "docs/source",
    "extra_navbar": "<p></p>",
    # -- Version switcher ----------------------------------------------------
    # The published docs are versioned: stable (master) lives at the site root,
    # `latest` (develop) and tagged releases live in sub-directories, and this
    # dropdown moves between them. `switcher.json` (deployed to the site root)
    # lists the versions; DOCS_VERSION tells the build which one it is, so the
    # dropdown highlights the right entry. The CI deploy sets it to the slug it
    # publishes ("stable", "latest", or a tag like "v1.1.0"); local builds
    # default to "latest".
    "switcher": {
        "json_url": "https://cosmostat.github.io/shapepipe/switcher.json",
        "version_match": os.environ.get("DOCS_VERSION", "latest"),
    },
    "navbar_end": ["version-switcher", "navbar-icon-links"],
}
html_collapsible_definitions = True
html_awesome_headerlinks = True

# The name for this set of Sphinx documents.  If None, it defaults to
# "<project> v<release> documentation".
html_title = "{0} v{1}".format(project, version)

# If not '', a 'Last updated on:' timestamp is inserted at every page bottom,
# using the given strftime format.
html_last_updated_fmt = "%d %b, %Y"

# If true, SmartyPants will be used to convert quotes and dashes to
# typographically correct entities.
html_use_smartypants = True

# If true, "Created using Sphinx" is shown in the HTML footer. Default is True.
html_show_sphinx = False

# If true, "(C) Copyright ..." is shown in the HTML footer. Default is True.
html_show_copyright = True


# -- Options for HTMLHelp output ---------------------------------------------

# Output file base name for HTML help builder.
htmlhelp_basename = "ShapePipedoc"

# -- Intersphinx Mapping ----------------------------------------------

# Refer to the package libraries for type definitions
intersphinx_mapping = {
    "python": ("http://docs.python.org/3", None),
    "numpy": ("https://numpy.org/doc/stable/", None),
    "scipy": ("https://docs.scipy.org/doc/scipy/reference", None),
    "matplotlib": ("https://matplotlib.org", None),
    "astropy": ("http://docs.astropy.org/en/latest/", None),
}

# -- BibTeX Setting  ----------------------------------------------

bibtex_bibfiles = ["refs.bib"]
bibtex_default_style = "unsrt"
bibtex_reference_style = "author_year"
