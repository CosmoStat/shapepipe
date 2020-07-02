# Post-processing

This page specifically deals with the post-processing step of results obtained
on [canfar](./canfar.md).

1. Check availability of results

   A `canfar` job can submit a large number of tiles, and many but not all might be processed
   simultaneously. To check which tiles are finished, and whose results have been uploaded, use
   ```bash
   scripts/python/canfar_avail_results.py -i <ID_files> -v
   ```
   E.g  with `-i aux/CFIS/tiles_202007/tiles_W3.txt` as input ID file.

1. Download results
