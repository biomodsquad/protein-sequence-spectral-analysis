Contributing to the documentation
=================================

The site follows the MISTIC documentation pattern: Sphinx, Furo, guided tutorials,
reproducible figures, downloadable examples, and an API reference generated from
public docstrings. Protspect's existing Markdown guides remain the authoritative
technical notes rather than being copied into separate pages.

Build and verify
----------------

From a development checkout:

.. code-block:: bash

   python -m pip install -e '.[docs]'
   python -m docs.generate_tutorial_figures
   python -m sphinx -W --keep-going -b html docs docs/_build/html
   python -m ci.check_docs_links docs/_build/html
   python -m sphinx -W --keep-going -b doctest docs docs/_build/doctest
   python -m ci.check_tutorials --output docs/_build/notebooks
   python -m http.server 8000 --directory docs/_build/html

Open ``http://localhost:8000`` in your browser. Generated HTML and executed
notebooks are ignored. Figure source lives in ``docs/generate_tutorial_figures.py``.
The downloadable source notebooks in ``protspect/examples`` contain no stale
outputs. CI executes them in a writable temporary directory with no network
requirements. Their synthetic sequences teach API usage, not biological accuracy.

``PROTSPECT_DOCS_REF`` controls GitHub source links for a branch preview;
``PROTSPECT_DOCS_BASEURL`` optionally sets the canonical deployed URL. The default
build is a development snapshot and does not claim a new PyPI release.

Review and publishing
---------------------

The documentation workflow builds PR previews and uploads HTML artifacts. It
also supports a manual Pages publication request; publication is off by default.
Publishing needs repository GitHub Pages configured to use GitHub Actions and
the ``github-pages`` environment. No repository Pages settings are changed by
building these sources.

Keep numerical examples runnable and link every new page in a toctree. Run the
strict HTML/doctest builds before opening a PR. Preserve the distinction between
full-coefficient exports and historical Fourier closure/dyadic exports. Never
bundle private peptide datasets or downloaded AAindex files. General contribution
instructions are in the repository's ``CONTRIBUTING.md``.
