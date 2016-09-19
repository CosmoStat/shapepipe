MultiFit v1.0.3
===============

This is the MultiFit (MulTipurpose Fitting module) README file.

Author: Marc Gentile (marc.gentile@epfl.ch)

Description
-----------

The main objective of the MultiFit framework is to provide a common fitting interface to different
parametric models and fitting algorithms

To this end, MultiFit provides a common programming interface (API) to different minimizers and 
parametric models. 
This allows in particular to:
- Test how well different minimizers perform when fitting the same model. This is useful, 
  for instance to evaluate different implementations of the Levenberg-Marquardt non-linear least 
  squares algorithm (LVM) (e.g. levmar, http://users.ics.forth.gr/~lourakis/levmar/ or 
  scipy.leastsq, http://docs.scipy.org/doc/scipy/reference/generated/scipy.optimize.leastsq.html). 
- Evaluate minimization algorithms other than LVM such as Coordinate Descent (CD), 
  http://en.wikipedia.org/wiki/Coordinate_descent.  
- Keeping the same minimizer and the same data, select the model that give the most accurate 
  results. In the context of astronomy, one can for instance compare models based on 
  single-component Sérsics and Bulge+Disk Sérsics.

This version is still in beta. Only the following minimizers and astronomical models are supported:
- Minimizers:
  - scileastsq: based on the scipy Levenberg-Marquardt non-linear least available with the 
                Python scipy library 
                (http://docs.scipy.org/doc/scipy/reference/generated/scipy.optimize.leastsq.html))
  - scdmin: an implementation of the Coordinate Descent algorithm, based on the
               SCDM (Simple Coordinate Descent Minimizer) minimizer (Marc Gentile, 2013).

- Parametric models:
  - scsersic: an implementation of galaxy model base on a pure disk Sérsic model 
              (http://en.wikipedia.org/wiki/Sersic_profile)  
  - gbdsersic: an implementation of galaxy model based on the combination of an exponential 
               Sérsic profile (sersic index 1) to model the galaxy disk and a de Vaucouleur Sérsic 
               profile (Sérsic index 4) to model the galaxy bulge. Both the disk and the bulge share 
               the same centroid and the same elipticity. The underlying implementation is done 
               using GalSim v1.x, the "modular galaxy image simulation toolkit", 
               https://github.com/GalSim-developers/GalSim.
  - gcbdsersic: an implementation of galaxy model identical to that of gbdsersic except that the
                disk Sérsic index can also vary. 
  - profile2D: model of a generic two-dimensional profile, mostly for testing purpose.

The software is designed so that additional models and fitting methods can be addded fairly quickly.

Documentation
-------------

A Doxygen reference documentation is available in the 'doc' directory.


Source code
-----------

The MultiFit code is organized as follows:

- Configuration files:
  - config/models: configuration for model: scsersic, gbdsersic, gcbdsersic, profile2D
    - scsersic.cfg 
    - profile2D.cfg
  - config/methods: configuration for methods: scdmin and scileastsq
    - scileastsq.cfg
    - scdmin.cfg
    - profile2D.cfg
  - config/fitting: extra fitting configuration that combine model and method information. 
                    At present, only the scdmin method uses such configuration. 
    -  scdmin_fitting: extra configuration for scdmin: scsersic, gbdsersic, gcbdsersic

- Module files:
  - modules/models: implementation of scsersic, gbdsersic, gcbdsersic, profile2D
    - scsersic.py
    - gbdsersic.py
    - gcbdsersic.py  
    - profile2D.py 
  - modules/methods: implementation of scileastsq, scdmin
    - scileastsq.py
    - scdmin.py

- Main multifit code:
  - multifit.py


Prerequisites 
-------------

- Mandatory:
   - Python 2.6 or 2.7 on Unix or Mac OS, with the corresponding numpy and scipy modules
   - sconfig v0.5.0 or greater: configuration file parser module
   - slogger v0.5.0 or greater: logging module

- Optional:
   - GalSim v1.0 or greater, if the gbdsersic and gcbdsersic galaxy models are used
   - scdm v0.5.0 or greater: if the SCDM minimizer is used   


Installation:
-------------

- Download and unpack the file multifit-1.0.3.tar.gz in the directory of your
  choice
- mpf is a pure Python module and can then be installed using the standard
  procedure: 

     1. create a console or terminal session
     2. cd to the directory where multifit-1.0.3.tar.gz was detared
     3. enter the command: sudo python setup.py install

  If no root access, use: python setup install --home=<target directory> to 
  install the module to <target directory>. This directory must be included in
  the PYTHONPATH environment variable. 

