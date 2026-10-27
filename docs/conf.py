"""Sphinx/Furo documentation for the development Protspect API."""

import os
from pathlib import Path
import re
import sys

try:
    import tomllib
except ModuleNotFoundError:  # Python 3.10
    import tomli as tomllib

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
project = "Protspect"
author = "Biomolecular Systems Modeling and Design Group"
copyright = "2026, Biomolecular Systems Modeling and Design Group"
release = tomllib.loads((ROOT / "pyproject.toml").read_text())["project"]["version"]
version = release
extensions = ["myst_parser", "sphinx.ext.autodoc", "sphinx.ext.autosummary",
              "sphinx.ext.napoleon", "sphinx.ext.viewcode", "sphinx.ext.mathjax",
              "sphinx.ext.doctest"]
source_suffix = {".rst": "restructuredtext", ".md": "markdown"}
master_doc = "index"
exclude_patterns = ["_build", "_static/**", ".DS_Store"]
myst_heading_anchors = 4
myst_enable_extensions = ["colon_fence", "dollarmath"]
autodoc_member_order = "bysource"
autodoc_typehints = "description"
napoleon_numpy_docstring = True
html_theme = "furo"
html_title = "Protspect · sequence spectral analysis"
html_static_path = ["_static"]
html_css_files = ["custom.css"]
html_logo = "_static/protspect-logo-spectrum.png"
html_favicon = "_static/protspect-logo-spectrum.png"
html_baseurl = os.environ.get("PROTSPECT_DOCS_BASEURL", "")
source_ref = os.environ.get("PROTSPECT_DOCS_REF", "main")
html_theme_options = {
    "source_repository": "https://github.com/biomodsquad/protein-sequence-spectral-analysis/",
    "source_branch": source_ref, "source_directory": "docs/",
    "announcement": "Protspect 0.1.0a2 alpha: APIs may change.",
    "light_css_variables": {"color-brand-primary": "#205b87", "color-brand-content": "#205b87"},
    "dark_css_variables": {"color-brand-primary": "#7ac4ef", "color-brand-content": "#7ac4ef"},
}
html_last_updated_fmt = None


def repository_links(app, docname, source):
    """Keep repository-relative Markdown source links valid in rendered HTML.

    Links to documents within docs/ are handled by MyST. Links to repository
    source/data outside the docs tree go to their source revision on GitHub.
    """
    if not app.env.doc2path(docname).endswith('.md'):
        return
    directory = Path(app.env.doc2path(docname)).parent
    docs_root = ROOT / 'docs'

    def replace(match):
        target = match.group(2)
        if '://' in target or target.startswith(('#', 'mailto:')):
            return match.group(0)
        path, _, anchor = target.partition('#')
        resolved = (directory / path).resolve()
        if resolved.is_relative_to(docs_root) or not resolved.is_relative_to(ROOT):
            return match.group(0)
        ref = os.environ.get('PROTSPECT_DOCS_REF', 'main')
        url = f'https://github.com/biomodsquad/protein-sequence-spectral-analysis/blob/{ref}/'
        url += resolved.relative_to(ROOT).as_posix()
        if anchor:
            url += '#' + anchor
        return f'{match.group(1)}({url})'

    source[0] = re.sub(r'(\[[^\]]*\])\(([^)]+)\)', replace, source[0])


def setup(app):
    """Adapt existing technical Markdown without duplicating its contents."""
    app.connect('source-read', repository_links)
