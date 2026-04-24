Work flow for upcoming v2.0 ShapePipe run:
- Run per tile
  The pipeline should be able to run on a single tile end-to-end, or to any user-specified point.
- Used exposures to be stored in central place, linked from tile, only used once for entire run
  Currently there is a mix of storage in central place (per patch), and tile directory.
- Runs no longer organised in patches.
- Use scratch
  For processing speed, we should be able to use the /scratch directory.
- File structure
  tiles/
   123/
        123.456
        123.666
	...
   124/
   ...
  exp/
   12/
   13/
   ..
 Explanation: The very large number of files (order 20000 tiles, a few 10000 single-HDU exposure files)
 requires smart organisation into subdirs. The above is a suggestion.
- Instead of symbolic links, could have (yaml) metadata of used exposures.
  Currently all required input files are stored in the tile dir or accessed with symlinks.
  It would probabyl be better to have a yaml or otherwise based database of used files per tile.
- Update ngmix
  Can come later.
- Ext catalogue instead of SExtractor detection
  Replace current sun of SExtracor on tiles with external ASCII catalogue. On exposures, we still run
  SExtractor to identify stars.
- Remove SM
  remve spread model.
