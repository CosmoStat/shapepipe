# job_list_help.bash
# Shared job-list description sourced by run_job_sp_canfar_v2.0.bash and
# job_sp_canfar_v2.0.bash so both -h outputs stay in sync.
#
# Defines: JOB_LIST_HELP  (ready to embed in a usage string passed to echo -ne)

JOB_LIST_HELP="\
   \t   1: retrieve tile images and weights (online if method=vos)\n\
   \t   2: uncompress weights (offline)\n\
   \t   4: find exposures (offline)\n\
   \t   8: retrieve exposures (online if method=vos)\n\
   \t  16: split exposures, get WCS header (offline)\n\
   \t  32: mask exposures (online if star_cat_for_mask=onthefly)\n\
   \t  64: process stars on exposures, PSF model (offline)\n\
   \t 128: merge exposure WCS headers into tile-level sqlite log\n\
   \t 256: object selection on tiles (online if UNIONS catalogue or tar_cat_for_mask=onthefly)\n\
   \t 512: create final catalogue\n\
   \t1024: process tiles (PSF interp, vignet, shape measurement)\n"
