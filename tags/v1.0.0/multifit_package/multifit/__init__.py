"""!
@mainpage MultiFit (MulTipurpose Fitting framework)

Author: Marc Gentile (marc.gentile@epfl.ch)

@section Description

The main objective of the MultiFit framework is to provide a common fitting interface to different
parametric models and fitting algorithms.

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
  single-component Sersics and Bulge+Disk Sersics.

This version is still in beta. Only the following minimizers and astronomical models are supported:
- Minimizers:
  - @c scileastsq: based on the scipy Levenberg-Marquardt non-linear least available with the 
                Python scipy library 
                (http://docs.scipy.org/doc/scipy/reference/generated/scipy.optimize.leastsq.html))
  - @c scdmin: an implementation of the Coordinate Descent algorithm, based on 
               SCDM (Simple Coordinate Descent Minimizer, Marc Gentile, 2013).

- Parametric models:
  - @c scsersic: an implementation of galaxy model base on a pure disk Sersic model 
              (http://en.wikipedia.org/wiki/Sersic_profile)  
  - @c gbdsersic: an implementation of galaxy model based on the combination of an exponential Sersic
               profile (sersic index 1) to model the galaxy disk and a de Vaucouleur Sersic profile
               (Sersic index 4) to model the galaxy bulge. Both the disk and the bulge share the
               same centroid and the same elipticity. The underlying implementation is done using
               @c GalSim v1.x, the "modular galaxy image simulation toolkit", 
               https://github.com/GalSim-developers/GalSim.
  - @c gcbdsersic: an implementation of galaxy model identical to that of @c gbdsersic except the
                disk Sersic index can also vary. 
  - @c profile2D: model of a generic two-dimensional profile, mostly for testing purpose.

The software is designed so that additional models and fitting methods can be addded fairly quickly.

@section code Source code

The MultiFit code is organized as follows:

- Configuration files:
  - config/models: configuration for model: scsersic, gbdsersic, gcbdsersic, profile2D
    - scsersic.cfg 
    - profile2D.cfg
  - config/methods: configuration for methods: @c scdmin and @c scileastsq
    - scileastsq.cfg
    - scdmin.cfg
    - profile2D.cfg
  - config/fitting: extra fitting configuration that combine model and method information. 
                    At present, only the scdmin method uses such configuration. 
    -  scdmin_fitting: extra configuration for scdmin: @c scsersic, @c gbdsersic, @c gcbdsersic

- Module files:
  - modules/models: implementation of @c scsersic, @c gbdsersic, @c gcbdsersic, @c profile2D
    - scsersic.py
    - gbdsersic.py
    - gcbdsersic.py  
    - profile2D.py 
  - modules/methods: implementation of @c scileastsq, @c scdmin
    - scileastsq.py
    - scdmin.py

- Main multifit code:
  - multifit.py
  
@section Prerequisites 

- Mandatory:
   - Python 2.6 or 2.7 on Unix or Mac OS, with the corresponding numpy and scipy modules
   - sconfig v0.5.0 or greater: configuration file parser module
   - slogger v0.5.0 or greater: logging module

- Optional:
   - GalSim v1.0 or greater, if the gbdsersic and gcbdsersic galaxy models are used
   - scdm v0.5.0 or greater: if the SCDM minimizer is used   

@section Installation

Download and unpack the file multifit-1.0.2.tar.gz in the directory of your
choice. The multifit module is a pure Python module and can thus be installed using the standard
procedure: 

<ol>
<li> create a console or terminal session </li>
<li> cd to the directory where multifit-1.0.2.tar.gz was detared </li>
<li> enter the command: command: @code sudo python setup.py install @endcode </li>

<li> if no root access to th default installation direcory, 
      use: @code python setup install --home=target_directory 
           @endcode to install the module to target_directory. This directory must be
           included in the @c PYTHONPATH environment variable </li>
</ol> 

"""

from multifit import *
from multifit_version import __version__


