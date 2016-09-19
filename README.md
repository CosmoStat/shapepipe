ShapePipe v1.0
This Pipeline is derived from the CEA-EPFL pipeline used for Great3.
This was originally coded by Marc Gentile (marc.gentile@epfl.ch) with
contributions by Florent Sureau (florent.sureau@cea.fr).

This is the main branch.

The pipeline currently includes the following modules:

- I/Os modules:
   * scatalog_package
   * sconfig_package
   * slogger_package

- Data/Processing Layers:
   * mpfg_package
   * mpfx_package
   * mpfg3_package
   * mpfcs82_package

- Software modules:
   * pse_package
   * mkpsf_package
   * mksim_package

- Main shape measurement module:
   * gfit_package/gfit_common_package (including wgfit)
   * isap_package (C++ iSAP code with python interface for
     wavelet transform)
   * multifit_package
   * scdm_package

