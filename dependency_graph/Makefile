NAME = cgat
ROOT = ../scripts ../CGAT ../../CGATPipelines/CGATPipelines
PDFS = $(NAME)_raw.pdf 
SVGS = $(NAME)_raw.svg
CHECK = $(NAME).check


FOOD_FLAGS?=--internal --follow

.SUFFIXES: .deps .dot .pdf .svg .clusters

all: $(PDFS) $(SVGS)

check: $(CHECK)

raw.deps: $(ROOT)
	sfood $(FOOD_FLAGS) $(ROOT) > $@

$(NAME).clusters: $(ROOT)
	cd $(ROOT) ; ls -1d * > $(shell pwd)/$@

$(NAME).deps: $(NAME).clusters raw.deps
	cat raw.deps | sfood-cluster -f $< > $@

$(NAME)_raw.deps: raw.deps
	cp $< $@

.deps.pdf:
	cat $< | python graph.py | dot -Tps | ps2pdf - $@

.deps.svg:
	cat $< | python graph.py | dot -Tsvg > $@

$(NAME).check:
	sfood-checker $(ROOT) >& $@

clean:
	rm -f *.clusters *.dot *.pdf *.svg
	ls -1 *.deps | grep -v ^raw.deps | xargs rm -f

realclean: clean
	rm -f raw.deps
