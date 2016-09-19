wgfit v1.0.0
===========

Galaxy shape measurement in wavelet space

Author: Florent Sureau (florent.sureau@cea.fr), from gift v1.0.0 written by Marc Gentile (marc.gentile@epfl.ch)

Description
-----------

The gfit method uses a forward model fitting algorithm to measure galaxy shapes. 

An earlier version of gfit, used to participate in the GREAT10 galaxy Challenge
(Kitching et al., 2011, Applied Statistics, Kitching et al, 2012, MNRAS), has been described
in Gentile et al, 2012, astro-ph/1211.4847 (http://arxiv.org/abs/1211.4847).

The following galaxy models are currently supported:
  - scsersic: an implementation of galaxy model base on a pure disk Sérsic model 
              (http://en.wikipedia.org/wiki/Sersic_profile)  
  - gscshsersic:
  - gbdsersic: an implementation of galaxy model based on the combination of an exponential 
               Sérsic profile (sersic index 1) to model the galaxy disk and a de Vaucouleur Sérsic 
               profile (Sérsic index 4) to model the galaxy bulge. Both the disk and the bulge 
               share the same centroid and the same elipticity. The underlying implementation is 
               done using GalSim v1.x, the "modular galaxy image simulation toolkit", see
               https://github.com/GalSim-developers/GalSim.
  - gbdshsersic:
  - gcbdsersic: an implementation of galaxy model identical to that of gbdsersic except that
                the disk Sérsic index can also vary.
  - gccbdsersic: an implementation of galaxy model identical to that of gcbdsersic except that
                 the disk and the bulge Sérsic index can also vary.


Fitting can currently be performed using two minimizers:
  - scileastsq: based on the Levenberg-Marquardt non-linear least-squares implementation, 
                available with the Python SciPy library 
                (http://docs.scipy.org/doc/scipy/reference/generated/scipy.optimize.leastsq.html))
  - scdmin: a wrapper to SCDM (Simple Coordinate Descent Minimizer, M. Gentile, 2011-2013)
            an implementation of the Coordinate Descent algorithm 
            (http://en.wikipedia.org/wiki/Coordinate_descent).

gfit is written in Python (http://www.python.org/) and takes advantage of the language
flexibility and richness of features. It relies heavily on the numpy and scipy 
Python C/C++/fortran libraries. 

The software has been designed in such a way that additional galaxy models and minimizers can be 
plugged-in fairly quickly.
 
The implementation of this version relies mainly on the following modules: 

- mpfg (Multi-Processing Framework, M. Gentile, 2013): a Python module for executing code on 
  both SMP (Symmetric Multiprocessing) or on serveurs supporting the Message Passing Interface 
  (MPI). Refer to the documentation of mpf module for additional details.
- mpfx (M. Gentile, 2013): adaptation layer based on mpfg allowing gfit to access images and 
  catalogs
- multifit (MulTipurpose Fitting framework, M. Gentile, 2013): a python module that offers a 
  common fitting interface to different parametric models and fitting algorithms. Refer to the
  documentation of Multifit for additional details.
- scdm (Simple Coordinate Descent Minimizer, M. Gentile, 2011-2013) an implementation of the 
  Coordinate Descent algorithm (http://en.wikipedia.org/wiki/Coordinate_descent). 
- Numpy v1.7.x., SciPy v0.9.x and pyfits Python modules 

Source code
-----------

The source code is organized as follows:

- Configuration files:
  - config/gfit.cfg (or any other name), containing the top gfit configuration parameters.  

  - config/models: configuration for model: scsersic, gbdsersic, gcbdsersic
    - scsersic.cfg 
  - config/methods: configuration for methods: scdmin and scileastsq
    - scileastsq.cfg
    - scdmin.cfg
  - config/fitting: extra fitting configuration that combine model and method information. 
                    At present, only the scdmin method uses such configuration. 
    - scdmin_fitting.cfg: extra configuration for scdmin: scsersic, gbdsersic, gcbdsersic

- Python code:
  - gfit_MPI.py: main entry point of the MPI executable 
  - gfit_SMP.py: main entry point of the SMP executable
  - gfit_args.py: handle command-line arguments and options
  - gfit_helper.py: helper class
  - gfit_job: job handling, specializations of job handling classes available from mpf 
              and mpfg3
  - gfit_shape: shape measurement
  - gfit_plot: plotting classes 
  - gfit_version: version information

Prerequisites 
-------------

- Mandatory:
   - Python 2.6 or 2.7 on Unix or Mac OS
   - Numpy v1.7.x, Scipy v0.9.x, Pyfits v3.1.x, Matplotlib 1.2.1 or greater
   - mpfg v1.x or greater: base multiprocessing framework
   - mpfx v1.x: adaptation of mpf for the GREAT3 dataset
   - multifit v1.x, MulTipurpose Fitting framework 
   - sconfig v1.x: configuration file parser module
   - slogger v1.x: logging module
   - scatalog v1.x: catalog management for astronomical data
     Note that this module also requires the following:
      - pyfits v3.1.x (http://www.stsci.edu/institute/software_hardware/pyfits)
      - AstroAsciiData v1.1 or greater (http://www.stecf.org/software/PYTHONtools/astroasciidata)

- Optional:
   - GalSim v1.x, if the gbdsersic and gcbdsersic galaxy models are used
   - scdm v0.6.x: if the Simple Coordinate Descent Minimizer (SCDM) is used   
   - dwtwiener v0.5.x: if the dwtwiener denoising module (Nurbaeva et al., 2011, A&A,
     http://www.aanda.org/articles/aa/full_html/2011/07/aa16556-11/aa16556-11.html) is used

Installation:
-------------

Download and untar the file gfit-1.0.0.tar.gz in the directory of your choice. 
The gfit executables are located in directory gfit-1.0.0/gfit. 

Execution
---------

To run gfit_MPI.py, one has to use the mpirun executable as:

mpirun -np N gfit_MPI.py [gfit options]or 
mpirun -np N python [python options] gfit_MPI.py [gfit options]@endcode
 
where N is the number of processors (or nodes) to use. MPI must has been 
installed and configured appropriately. 
On a cluster with a queuing system managed by a job scheduler, mpirun will have to be integrated 
in a shell script, invoked using a command like qsub. 
   
Since the manager process of gfit uses one processor for itself, workers will share N-1 processors. 

For example: mpirun -np 6 gfit_MPI -d mydir -c myconfig.cfg will run gfit
on 6 processors, 1 for the manager and 5 for the workers. The manager and each
of the workers will look for a configuration file: ./mydir/myconfig.cfg.

@note Before running gfit, edit the configuration file to set the BASE_INPUT_DIR value, which
      should point to the directory where the GREAT3 files reside. 
@note In the configuration file, any environment variable prefixed with $ such as $HOME will
      be expanded when read.

Configuration
-------------

Refer to the HTML documentation in directory gfit-1.0.0/doc. The main entry point is index.html.

