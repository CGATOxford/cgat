## Note that this solution requires runipy to be installed
## see: http://stackoverflow.com/questions/17905350/running-an-ipython-notebook-non-interactively

## Note: To build the notebooks, you need to have jupyter installed:
## conda install jupyter nb_conda nb_conda_kernels nb_anacondacloud
## conda install r-gplots r-rcolorbrewer r-ggplot2
##
## pip install bx-python
##
## mymodule load bio/picard-tools
## Dealing with errors:

## 1.readline can be an issue:
## ImportError: /ifs/home/andreas/miniconda3/envs/cgat-production-devel/bin/../lib/libreadline.so.6: undefined symbol: PC
## To fix:
##    conda uninstall readline --force

## In order to build the recipes, make sure that the CGAT 
## code collection has been installed as the recipes use
## the cgat frontend. Building the recipes in a development
## environment is possible, but there are ugly warnings
## in the notebooks about modules being loaded multiple
## times.

## primary targets:
##
## 
## make build
##    rebuilds all notebooks and convrets them to html. 
##    Use for testing new release.
##
## make export
##    converts all notebooks to html. Use for publishing.
##

NOTEBOOKS_= $(wildcard Recipe*)
NOTEBOOKS := $(filter-out %.build, $(NOTEBOOKS_))

build: $(NOTEBOOKS:%=%.build)
export: $(NOTEBOOKS:%=%.export)


# change into working directory
# link all data files
# run ipython notebook
%.build: %
	cd $*; \
	ln -fs /ifs/public/cgat/recipes/* .; \
	jupyter nbconvert \
		--ExecutePreprocessor.timeout=600 \
		--to html \
		--output ../../doc/recipes/$*.html \
		--template=full \
		--execute \
		*.ipynb >& ../$@

# export all ipython notebooks without recomputing
%.export: %
	cd $*; \
	ipython nbconvert --to html $*.ipynb; \
	cp -f $*.html ../../doc/recipes/$*.html; \
	touch $@

clean:
	rm -f */*.bam */*.bam.bai */*.fastq* */README */Makefile */metrics* */*.log
