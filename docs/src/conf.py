# Configuration file for the Sphinx documentation builder.
#
# For the full list of options see:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Project information -----------------------------------------------------

project = "MolOrient"
author = "MolOrient Developers"
copyright = "2026, MolOrient Developers"

# -- General configuration ---------------------------------------------------

extensions = [
    "myst_parser",
    "sphinxcontrib.bibtex",
    "sphinxcontrib.mermaid",
    "autoapi.extension",
]

source_suffix = {".rst": "restructuredtext", ".md": "markdown"}
exclude_patterns = ["_build"]

# MyST: render ```mermaid fences as diagrams and enable $...$ math.
myst_fence_as_directive = ["mermaid"]
myst_enable_extensions = ["dollarmath"]

# Bibliography
bibtex_bibfiles = ["CITATION.bib"]
bibtex_default_style = "unsrt"
bibtex_reference_style = "author_year"

# API documentation, generated from the source's docstrings
autoapi_dirs = ["../../src"]
autoapi_root = "api"
autoapi_add_toctree_entry = False
autoapi_options = [
    "members",
    "undoc-members",
    "show-inheritance",
    "show-module-summary",
]

# -- Options for HTML output -------------------------------------------------

html_theme = "pydata_sphinx_theme"
html_title = "MolOrient"
html_theme_options = {
    "github_url": "https://github.com/rmrresearch/molorient",
    "navbar_align": "left",
    "show_toc_level": 2,
}
