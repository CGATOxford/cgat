################################################################################
#   Gene prediction pipeline 
#
#   $Id: Makefile.plot 15 2005-08-09 15:24:35Z andreas $
#
#   Copyright (C) 2004 Andreas Heger
#
#   This program is free software; you can redistribute it and/or
#   modify it under the terms of the GNU General Public License
#   as published by the Free Software Foundation; either version 2
#   of the License, or (at your option) any later version.
#
#   This program is distributed in the hope that it will be useful,
#   but WITHOUT ANY WARRANTY; without even the implied warranty of
#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#   GNU General Public License for more details.
#
#   You should have received a copy of the GNU General Public License
#   along with this program; if not, write to the Free Software
#   Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.
#################################################################################

## plotting of histograms
%.png: %.hist
	python $(DIR_SCRIPTS_TOOLS)plot_histogram.py --hardcopy=$@ --title=$* $(PLOT_OPTIONS) < $*.hist

%.plot: %.hist
	python $(DIR_SCRIPTS_TOOLS)plot_histogram.py --title=$* $(PLOT_OPTIONS) < $*.hist

## plotting of scatterplots
%.png: %.data
	python $(DIR_SCRIPTS_TOOLS)plot_data.py --hardcopy=$@ --title=$* $(PLOT_OPTIONS) < $*.data

%.plot: %.data
	python $(DIR_SCRIPTS_TOOLS)plot_data.py --title=$* $(PLOT_OPTIONS) < $*.data


combined_%.hist:
	rm -f $@
	python $(DIR_SCRIPTS_TOOLS)combine_histograms.py \
	--headers=`echo *_$*.hist | perl -pe "s/.hist//g; s/ +/,/g" ` *_$*.hist > $@

combined_%.stats:
	rm -f $@
	python $(DIR_SCRIPTS_TOOLS)combine_tables.py \
	--headers=`echo *_$*.stats | perl -pe "s/.stats//g; s/ +/,/g" ` *_$*.stats > $@
