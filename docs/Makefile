SPHINXBUILD ?= python -m sphinx
.PHONY: html doctest clean
html:
	$(SPHINXBUILD) -W --keep-going -b html . _build/html
doctest:
	$(SPHINXBUILD) -W --keep-going -b doctest . _build/doctest
clean:
	$(SPHINXBUILD) -M clean . _build
