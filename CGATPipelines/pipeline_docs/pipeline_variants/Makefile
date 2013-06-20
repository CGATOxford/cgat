# Makefile for Sphinx documentation
#

# You can set these variables from the command line.
SPHINXOPTS    =
SPHINXCLEAN   = sphinxreport-clean
SPHINXBUILD   = sphinxreport-build --num-jobs=5 sphinx-build
PAPER         =

# Internal variables.
PAPEROPT_a4     = -D latex_paper_size=a4
PAPEROPT_letter = -D latex_paper_size=letter
ALLSPHINXOPTS   = -d _build/doctrees $(PAPEROPT_$(PAPER)) $(SPHINXOPTS) .

VERSION=version6refseq

DIR_PUBLISH=/net/fgu-project/collaborations/web/andreas/mousestrains/$(VERSION)
URL_SUB="s/\/net\/cpp-compute\/backup\/andreas\/projects\/mousestrains\/$(VERSION)\/doc/http:\/\/wwwfgu.anat.ox.ac.uk\/~andreas\/mousestrains\/$(VERSION)/" 
DIR_EXPORT=wwwrun@fgu1147:/backup/live/mousestrains/$(VERSION)

#html port
PORT=8100

.PHONY: help clean html dirhtml pickle json htmlhelp qthelp latex changes linkcheck doctest

help:
	@echo "Please use \`make <target>' where <target> is one of"
	@echo "  setup     to setup directories for the analysis"
	@echo "  html      to make standalone HTML files"
	@echo "  clean     to clean previous run"
	@echo "  check     to check the trackers in the python directory and collect all the data."
	@echo "  gallery   to update the thumbnail gallery - you might need to run 'make html' again afterwards."
	@echo "  dirhtml   to make HTML files named index.html in directories"
	@echo "  pickle    to make pickle files"
	@echo "  json      to make JSON files"
	@echo "  htmlhelp  to make HTML files and a HTML help project"
	@echo "  qthelp    to make HTML files and a qthelp project"
	@echo "  latex     to make LaTeX files, you can set PAPER=a4 or PAPER=letter"
	@echo "  changes   to make an overview of all changed/added/deprecated items"
	@echo "  linkcheck to check all external links for integrity"
	@echo "  doctest   to run all doctests embedded in the documentation (if enabled)"
	@echo "  server-start    to start a web server."
	@echo "  server-stop     to start a web server."


setup:
	@if test ! -d analysis; then mkdir analysis; fi
	@for x in templates/*.tpl_rst; do \
		xx=$${x#templates/}; \
		fn=$${xx%.tpl_rst}; \
		if test ! -e analysis/$${fn}.rst; then \
			echo "creating $${fn}"; \
			cp $${x} analysis/$${fn}.rst; \
		fi; \
	done

html:
	$(SPHINXBUILD) -b html $(ALLSPHINXOPTS) _build/html
	@echo
	@echo "Build finished. The HTML pages are in _build/html."

clean:
	$(SPHINXCLEAN) clean

distclean:
	$(SPHINXCLEAN) distclean

server-start:
	nohup nice -19 ./serve.py $(PORT) 2>&1 >> server.log &

server-stop:
	pkill "serve.py $(PORT)"

check:
	sphinxreport-test
	@echo
	@echo "Check finished. All data has been saved in the cache."

gallery:
	sphinxreport-gallery
	@echo
	@echo "The gallery has been rebuilt. Running 'make html' again."
	$(MAKE) html

dirhtml:
	$(SPHINXBUILD) -b dirhtml $(ALLSPHINXOPTS) _build/dirhtml
	@echo
	@echo "Build finished. The HTML pages are in _build/dirhtml."

pickle:
	$(SPHINXBUILD) -b pickle $(ALLSPHINXOPTS) _build/pickle
	@echo
	@echo "Build finished; now you can process the pickle files."

json:
	$(SPHINXBUILD) -b json $(ALLSPHINXOPTS) _build/json
	@echo
	@echo "Build finished; now you can process the JSON files."

htmlhelp:
	$(SPHINXBUILD) -b htmlhelp $(ALLSPHINXOPTS) _build/htmlhelp
	@echo
	@echo "Build finished; now you can run HTML Help Workshop with the" \
	      ".hhp project file in _build/htmlhelp."

qthelp:
	$(SPHINXBUILD) -b qthelp $(ALLSPHINXOPTS) _build/qthelp
	@echo
	@echo "Build finished; now you can run "qcollectiongenerator" with the" \
	      ".qhcp project file in _build/qthelp, like this:"
	@echo "# qcollectiongenerator _build/qthelp/Test.qhcp"
	@echo "To view the help file:"
	@echo "# assistant -collectionFile _build/qthelp/Test.qhc"

latex:
	$(SPHINXBUILD) -b latex $(ALLSPHINXOPTS) _build/latex
	@echo
	@echo "Build finished; the LaTeX files are in _build/latex."
	@echo "Run \`make all-pdf' or \`make all-ps' in that directory to" \
	      "run these through (pdf)latex."

changes:
	$(SPHINXBUILD) -b changes $(ALLSPHINXOPTS) _build/changes
	@echo
	@echo "The overview file is in _build/changes."

linkcheck:
	$(SPHINXBUILD) -b linkcheck $(ALLSPHINXOPTS) _build/linkcheck
	@echo
	@echo "Link check complete; look for any errors in the above output " \
	      "or in _build/linkcheck/output.txt."

doctest:
	$(SPHINXBUILD) -b doctest $(ALLSPHINXOPTS) _build/doctest
	@echo "Testing of doctests in the sources finished, look at the " \
	      "results in _build/doctest/output.txt."

publish:
	-mkdir -p $(DIR_PUBLISH)
	@echo "copying data"
	rsync --size-only -acvz _build/html/ $(DIR_PUBLISH)
	rsync --size-only --copy-links -acvz export $(DIR_PUBLISH)
	@echo "setting permissions"
	cp htaccess $(DIR_PUBLISH)/.htaccess
	chmod -R a+x $(DIR_PUBLISH)
	chmod -R a+r $(DIR_PUBLISH)
	@echo "fixing URLs"
	find $(DIR_PUBLISH)/ -name "*.html" -exec \
		perl -p -i -e $(URL_SUB) {} \;
