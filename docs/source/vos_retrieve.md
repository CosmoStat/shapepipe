## Retrieve files from VOspace

This page describes how ShapePipe output files can be retrieved via the Virtual Observatory Space
on canfar. This system was used for the CFIS v0 and v1 runs, and is now obsolete.

1. Retrieve ShapePipe result files 

  For a local run on the same machine as for post-processing, nothing needs to be done. In some cases, the run was carried out on a remote machine or cluster, and the resulting ShapePipe output files need to be retrieved.

  In the specific case of canfar_avail_results.py, this is done as follows.

    A. Check availability of results

      A canfar job can submit a large number of tiles, whose processing time can vary a lot. We assume that the submitted tile ID list is available locally via the ascii file tile_numbers.txt. To check which tiles have finished running, and whose results have been uploaded, use
