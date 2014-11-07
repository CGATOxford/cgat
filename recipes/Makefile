## Note that this solution requires runipy to be installed
## see: http://stackoverflow.com/questions/17905350/running-an-ipython-notebook-non-interactively

## primary targets:

## 
## make build
##    rebuilds all notebooks and convrets them to html. 
##    Use for testing new release.
##
## make export
##    converts all notebooks to html. Use for publishing.
##

NOTEBOOKS= $(wildcard Recipe*)

build: $(NOTEBOOKS:%=%.build)
export: $(NOTEBOOKS:%=%.export)


# change into working directory
# link all data files
# run ipython notebook
%.build: %
	cd $*; \
	ln -s /ifs/public/cgat/recipes/* .; \
	runipy --template=full --html ../../doc/recipes/$*.html *.ipynb >& ../$@

# export all ipython notebooks without recomputing
%.export: %
	cd $*; \
	ipython nbconvert --to html $*.ipynb; \
	cp -f $*.html ../../doc/recipes/$*.html; \
	touch $@

clean:
	rm -f */*.bam */*.bam.bai */*.fastq* */README */Makefile */metrics* */*.log
