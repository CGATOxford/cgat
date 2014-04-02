## Note that this solution requires runipy to be installed

## see: http://stackoverflow.com/questions/17905350/running-an-ipython-notebook-non-interactively

NOTEBOOKS= $(wildcard Recipe*)

HTML=$(NOTEBOOKS:%=%.html)

html: $(HTML)

# change into working directory
# link all data files
# run ipython notebook
%.html: %
	cd $*; \
	ln -s /ifs/public/cgat/recipes/* .; \
	runipy --html ../$* *.ipynb >& rebuild.log

clean:
	rm -f */*.bam */*.bam.bai */*.fastq*
