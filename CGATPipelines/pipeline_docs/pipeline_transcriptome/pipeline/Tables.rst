*****************************
Administrative matters
*****************************

This page presents sanity checks on the tables in 
the SQL database. It will present an overview
of which tables are present in the database and
how complete these are.

Completeness of tables
======================

The following table `PipeTableCheckTables`_ lists the
number of transcripts (identified by the column *gene_id*)
in each table.

.. _PipeTableCheckTables:

.. report:: Trackers.TrackerSQLCheckTablesGeneId
   :render: matrix

   Number of transcripts in each table.

The figure PipeFigCheckTables_ shows the same data
as `PipeTableCheckTables`_, but the data is normalized
by row.

.. _PipeFigCheckTables:

.. report:: Trackers.TrackerSQLCheckTablesGeneId
   :render: matrix-plot
   :transform-matrix: normalized-row-max
   :palette: Reds
   :reverse-palette:
   
   Number of transcripts in each table.


Computation of rates
====================

The table `PipeTableCheckTableEvol`_ lists the number of 
rates and other entries in the tables ending on *evol*.
This tables contain a summary of evolutionary rates that
could be computed for each track.

.. _PipeTableCheckTableEvol:

.. report:: Trackers.TrackerSQLCheckTableEvol
   :render: matrix

   Number of rates for each transcript.

The figure PipeFigCheckTableEvol_ shows the same data
as `PipeTableCheckTableEvol`_, but the data is normalized
by row.

.. _PipeFigCheckTableEvol:

.. report:: Trackers.TrackerSQLCheckTableEvol
   :render: matrix-plot
   :transform-matrix: normalized-row-max
   :palette: Reds
   :reverse-palette:
   
   Number of rates for each transcript.

