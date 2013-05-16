*******************************
Transcript comparison pipeline
*******************************

Purpose
-------

Map 454 reads onto a genome and assemble overlapping
transcripts into transcript models.

The pipeline currently does not use base quality information
during mapping and does not consider alternative transcripts.

Setting up
----------

To set up the pipeline in the current directory run::

   python setup.py --method=compare_transcripts > setup.log


Link towards the genome from /net/cpp-data/backup/databases/indexed_fasta and
call the files genome.fasta and genome.idx. For example::
       
   ln -s /net/cpp-data/backup/databases/indexed_fasta/hs_ncbi36_softmasked.fasta genome.fasta
   ln -s /net/cpp-data/backup/databases/indexed_fasta/hs_ncbi36_softmasked.idx genome.idx

Input (required): 

%.gtf
   gtf files with (experimental) transcripts. The % denotes the
   track name, for example heart.gtf, kidney.gtf, sample1.gtf, ...

genome.fasta, genome.idx
   an indexed genome :term:`PARAM_GENOME`. See also index_fasta.py.
   
ensembl.gtf
   a gtf file with a reference sequence set: the default is ``ensembl``, 
   but can be changed in :term:`PARAM_MASTER_SET_GENES`

annotations.gff
   a gff file with annotated genomic regions. See 
   :term:`PARAM_GENOME_REGIONS`. Use gtf2gff.py to create this file.

The pipeline includes additional information if it is present:

%.coverage
   table with coverage information for a track. The output 
   is from blat2assembly.py

%.polyA 
   information about polyA tails. The output is from 
   blat2assembly.py

%.readstats
   a table with read alignment statistics after filtering
   (see output from MapTranscripts)

%.readmap
   a table mapping gene_ids to read_ids after filtering
   (see output from MapTranscripts)

%.readinfo
   a table with read information.

%.readgtf
   mapped locations of reads after filtering.

:term:`PARAM_FILE_REPEATS_RATES`
   gff formatted file of ancesctral repeats.
   The score field contains the rate (see Makefile.ancestral_repeats)

:term:`PARAM_FILE_REPEATS`
   gff file with repeats in genome. These are used for masking in coding 
   potential predictions. 

:term:`PARAM_FILE_REPEATS_GC`
   gff formatted file of ancesctral repeats.
   The score field contains the G+C content (see Makefile.ancestral_repeats)

:term:`PARAM_FILE_ALIGNMENTS`
   psl formatted file with genomic alignments
   between this species in query and another at appropriate evolutionary 
   distance in target. 
   

:term:`PARAM_FILENAME_GO` and :term:`PARAM_FILENAME_GOSLIM`
   GO annotations for genes in the reference set. Example format is::

       cell_location   ENSPPYG00000000676      GO:0016020      membrane        NA

:term:`PARAM_FILENAME_TERRITORIES`
   gene territories. GTF formatted file, an example entry would be::

       chr1    protein_coding  exon    3979975 4199559 .       -       .       transcript_id "ENSPPYG00000000050"; gene_id "ENSPPYG00000000050";#
 
:term:`PARAM_CPC_UNIREF`
   uniref database to use for coding potential predictions.

Output from the mapTranscripts454 project can be imported with a single command::

   make PATH_TO_MAPPING_DIR.add-tracks

Configuration
-------------

Edit the :file:`Makefile` to configure the pipeline. See Parameters_ below.

Usage
-----

The pipeline is controlled by running `make`_ targets. The results of the pipeline
computation are stored as tab separated tables in the working directory. Most of these
tables are then imported into an `sqlite`_ database called ``csvdb`` (see :term:`PARAM_DATABASE`).

Annotation
~~~~~~~~~~

Type::

   make all

to do all.

Fine grained control
++++++++++++++++++++

A more complete list of targets:

all
   make all

build
   only build, but do not import.

import
   import


Visualization
~~~~~~~~~~~~~

The following targets aid visualizatiov:

   ucsc-tracks-gtf
      export the segments as compressed gtf files. Can be viewed as
      user tracks in the `ucsc`_ genome browser.

GO analysis
~~~~~~~~~~~

GO analysis will compute the relative enrichment/depletion gene sets.

Requires :term:`PARAM_FILENAME_TERRITORIES`, :term:`PARAM_FILENAME_GO`
and :term:`PARAM_FILENAME_GOSLIM` to be set.

There are two counting methods. The first method (``go``) assigns GO terms associated with 
the reference gene set to TLs and counts these. The second method (``territorygo``) assigns TLs to genes
in the reference set and then does a GO analysis on theses.

.. note::
   The convential GO analysis based on gene list is the ``territorygo`` method. 

Usage
+++++

Usage::

   make <track>:<slice>:<subset>:<background>.<go>.<method>analysis

The fields are:

track
   the data track to be chosen. 

slice
   the slices correspond to flags in the table <track>_annotation. Use ``all``
   to use all segments in a ``track``.

subset
   the subset corresponds to a table that is joined with <track>_annotation to
   restrict segments to a user-specified set. Use ``all`` for no restriction.

background
   the background gene set

go
   either ``go`` or ``goslim``

method
   either ``go`` or ``goterritory``

Results will be in the directory :file:`<track>:<slice>:<subset>:<background>.<go>.<method>analysis.dir`.

For example::

   make thoracic:known:all:thoracic.go.goanalysis

will compute the enrichment of protein coding TL in the track ``thoracic`` using
all ``thoracic`` genes as the background. 

The command::

   make thoracic:known:all:ensembl.goslim.territorygoanalysis

will compute  ``goslim`` term enrichment. The foreground set are genes from the reference set (``ensembl``) overlapping 
protein coding TL in the track ``thoracic``. The background is the complete reference gene
set (``ensembl``).

Annotator analysis
~~~~~~~~~~~~~~~~~~

Annotator computes the statistical significance of enrichment/depletion
of genomic features (called segments) within genomic regions (called annotations).

To run annotator analysis, two files need to be present:

1. A workspace
2. A collection of annotations on the genome

Building workspaces
+++++++++++++++++++

Workspaces are built using makefile targets. For example to build :file:``genome.workspace``, type::
   make genome.workspace

All workspaces exclude contigs called matching ``random``.

genome.workspace
   full genome

intergenic.workspace
   only intergenic regions

intronic.workspace
   only intronic regions

unknown.workspace
   both intergenic and intronic regions

territories.workspace
   workspace of territories

alignable.workspace
   only segments that can be aligned to a reference genome (see :term:`...`).

There is a convenience target::

   make annotator-workspaces

that will build all available workspaces.

Annotations
+++++++++++

Annotations are built using makefile targets.

all.annotations:
   all subsets (all/known/unknown) for each track.

architecture.annotations:
   annotations according to genes (intronic, intergenic, ...).

{all,known,unknown}_sets.annotations
   annotations of known, unknown, all transcripts

allgo_territories.annotations
   territories annotation with GO categories

allgoslim_territories.annotations
   territories annotation with GOSlim categories

intronicgo_territories.annotations
   territories annotation with GO categories

intronicgoslim_territories.annotations
   territories annotation with GOSlim categories

intergenicgo_territories.annotations
   territories annotation with GO categories

intergenicgoslim_territories.annotations
   territories annotation with GOSlim categories

There is a convenience target::

   make annotator-annotations

that will build all available annotations.

Usage
+++++

In order to perform :term:`Annotator` analyses, you run a make target::

   make <track>:<slice>:<subset>:<workspace>:<workspace2>_<annotations>.annotators

The fields determine which segments are used for the enrichment analysis.

track
   the data track to be chosen. 

slice
   the slices correspond to flags in the table <track>_annotation. Use ``all``
   to use all segments in a ``track``.

subset
   the subset corresponds to a table that is joined with <track>_annotation to
   restrict segments to a user-specified set. Use ``all`` for no restriction.

workspace
   the workspace to be used

workspace2
   a second workspace. The actual workspace will be the intersection of both workspaces.

annotations
   annotations to use.

.. note:: 
   Annotations, segments and the workspace need to be chosen carefully for each experiment.
   For example, failing to use territories for goterritory analysis will measure enrichment
   of segments within goterritories in general, and not necessarily relative enrichment 
   between go territories.

The results will be in the file :file:`<track>:<slice>:<subset>:<workspace>:<workspace2>_<annotations>.annotators`.

Examples
++++++++

The command::

   make thoracic:unknown:all:intergenic:all_unknownsets.annotators

will test for enrichment among ``unknown`` transcripts in the track ``thoracic``
with intergenic segments the other sets. The command::

   make thoracic:intronic:all:intronic:territories_intronicgoslimterritories.annotators

will check for enrichment of ``intronic`` transcripts from the track ``merged``
within intronic genomic segments that also have GO assignments (intersection
of workspaces ``intronic`` and ``territories``. It will label GO territories
by GOslim territories.

Association analysis
~~~~~~~~~~~~~~~~~~~~

Association analysis computes the significance of finding segments close
to annotations.

Type::

   make annotator-distance-run

to run all association analyses.

Parameters
----------

The following parameters can be set in the :file:`Makefile`:

.. report:: Trackers.MakefileParameters
   :render: glossary
   :tracks: Makefile.compare_transcripts
   :transpose:

   Overview of pipeline parameters.

.. _make: http://www.gnu.org/software/make

.. _sqlite: http://www.sqlite.org 

.. _ucsc: http://genome.ucsc.edu
