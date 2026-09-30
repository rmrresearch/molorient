# Building the Documentation

The MolOrient documentation is built with [Sphinx](https://www.sphinx-doc.org)
using the [PyData Sphinx Theme](https://pydata-sphinx-theme.readthedocs.io).
Pages are written in Markdown (via [MyST](https://myst-parser.readthedocs.io)),
and the API reference is generated from the docstrings in `src/molorient`. The
documentation source lives in `docs/src`.

## Building Locally

From the root of the repository:

```bash
# (Optional, but highly recommended)
python3 -m venv .venv
source .venv/bin/activate

# Build the docs
pip install -r docs/requirements.txt
sphinx-build -b html docs/src docs/_build/html
```

Then open `docs/_build/html/index.html` in a browser.

CI builds the documentation with warnings treated as errors. To check that
your changes will pass, build with:

```bash
sphinx-build -W --keep-going -b html docs/src docs/_build/html
```

## Publishing

The documentation is published to GitHub Pages automatically whenever changes
are merged into `main`.
