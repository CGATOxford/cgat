=======
Methods
=======

Tiling strategies
==================

Differentially methylated regions are identified via three tiling
strategies:

1. Variable width tiles - densities are clustered across all samples
and tiles are defined based on these clusters.

2. Non-overlapping fixed width tiles - fixed width tiles of a certain
size (default = 1kb) are tiled over the whole genome. Tiles are
touching, but non-overlapping.

3. Overlapping fixed width tiles - fixed width tiles of a certain size
(default = 1kb) are tiled over the whole genome. Tiles are overlapping
by a certain amount (default = 500bp).

Differential methylation
========================

Differential methylation is assessed per window. Read counts are
collected per window. Tag counting (DESeq and EdgeR) are then employed
to estimate differential methylation.

Merging tiles
=============

After defininig differentially methylated regions, adjacent and
overlapping regions are merged that are called the same (DMR or not
DMR) and show the same direction of fold change.
