"""!
@mainpage Galaxy Shape Measurement

Author: Marc Gentile (marc.gentile@epfl.ch)

@section Description

The @c gfit method uses a forward model fitting algorithm to measure galaxy shapes. 

An earlier version of gfit, used to participate in the @c GREAT10 galaxy Challenge
(Kitching et al., 2011, Applied Statistics, Kitching et al, 2012, MNRAS), has been described
in Gentile et al, 2012, astro-ph/1211.4847 (http://arxiv.org/abs/1211.4847).

The present version has been customized to participate in the @c GREAT3 Challenge 
(http://great3.projects.phys.ucl.ac.uk/leaderboard/).

The following galaxy models are currently supported:
  - @c scsersic: an implementation of galaxy model base on a pure disk Sérsic model 
              (http://en.wikipedia.org/wiki/Sersic_profile)  
  - @c gbdsersic: an implementation of galaxy model based on the combination of an exponential 
               Sérsic profile (sersic index 1) to model the galaxy disk and a de Vaucouleur Sérsic 
               profile (Sérsic index 4) to model the galaxy bulge. Both the disk and the bulge 
               share the same centroid and the same elipticity. The shear estimator is the 
               distortion |e| = (a^2 - b^2)/(a^2 + b^2). The underlying implementation is 
               done using @c GalSim v1.x, the "modular galaxy image simulation toolkit", see
               https://github.com/GalSim-developers/GalSim.
  - @c gbdshsersic:
  - @c gbdshsersic: an implementation of galaxy model identical to that of @c gbdsersic except that
               the reduced shear |g| = (a - b)/(a + b) is used insead of the distortion |e|. 
  - @c gcbdsersic: an implementation of galaxy model identical to that of @c gbdsersic but with
                   a varying disk Sérsic index.
  - @c gccbdsersic: an implementation of galaxy model identical to that of gcbdsersic except that
                    the disk and the bulge Sérsic index can also vary.

Fitting can currently be performed using two minimizers:
  - @c scileastsq: based on the <tt>Levenberg-Marquardt</tt> non-linear least-squares 
                   implementation, available with the Python SciPy library 
                (http://docs.scipy.org/doc/scipy/reference/generated/scipy.optimize.leastsq.html))
  - @c scdmin: a wrapper to SCDM (Simple Coordinate Descent Minimizer, M. Gentile, 2011-2013)
               an implementation of the Coordinate Descent algorithm 
               (http://en.wikipedia.org/wiki/Coordinate_descent).

@c gfit is written in @c Python (http://www.python.org/) and takes advantage of the language
flexibility and richness of features. It relies heavily on the @c numpy and @c scipy 
@c Python C/C++/fortran libraries. 

The software has been designed in such a way that additional galaxy models and minimizers can be 
plugged-in fairly quickly.
 
The implementation of this version relies mainly on the following modules: 

- ../../mpf_package/mpf (Multi-Processing Framework, M. Gentile, 2013): a Python module for executing code on 
  both SMP (Symmetric Multiprocessing) or on serveurs supporting the Message Passing Interface 
  (MPI). Refer to the documentation of @c mpf module for additional details.
- @c mpfg3 (M. Gentile, 2013): an extension of @c mpf that provides an adaptation layer to access
  the @c GREAT3 datasets. Refer to the documentation of @c mpfg3  module for additional details.
- @c multifit (MulTipurpose Fitting framework, M. Gentile, 2013): a python module that offers a 
  common fitting interface to different parametric models and fitting algorithms. Refer to the
  documentation of @c Multifit for additional details.
- @c scdm (Simple Coordinate Descent Minimizer, M. Gentile, 2011-2013) an implementation of the 
  Coordinate Descent algorithm (http://en.wikipedia.org/wiki/Coordinate_descent). 
- Numpy v1.7.x., SciPy v0.9.x and pyfits Python modules 

@section code Source code

The source code is organized as follows:

- Configuration files:
  - @c config/gfit.cfg (or any other name), containing the top @c gfit configuration parameters.  

  - config/models: configuration for model: @c scsersic, @c gbdsersic, @c gcbdsersic
    - @c scsersic.cfg 
  - config/methods: configuration for methods: @c scdmin and @c scileastsq
    - @c scileastsq.cfg
    - @c scdmin.cfg
  - config/fitting: extra fitting configuration that combine model and method information. 
                    At present, only the @c scdmin method uses such configuration. 
    - @c scdmin_fitting.cfg: extra configuration for @c scdmin: @c scsersic, @c gbdsersic, 
                                                     @c gcbdsersic

- Python code:
  - @c gfit_MPI.py: main entry point of the MPI executable 
  - @c gfit_SMP.py: main entry point of the SMP executable
  - @c gfit_args.py: handle command-line arguments and options
  - @c gfit_helper.py: helper class
  - @c gfit_job: job handling, specializations of job handling classes available from @c mpf 
              and @c mpfg3
  - @c gfit_shape: shape measurement
  - @c gfit_filter: object filtering
  - @c gfit_plot: plotting classes 
  - @c gfit_version: version information
  

@section Prerequisites 

- Mandatory:

   - @c Python 2.6 or 2.7 on Unix or Mac OS
   - @c Numpy v1.7.x, Scipy v0.9.x, Pyfits v3.1.x, Matplotlib 1.2.1 or greater
   - @c mpfg v1.x or greater: base multiprocessing framework
   - @c mpfx v1.x or greater: adaptation of mpf for the GREAT3 dataset
   - @c multifit v1.x, MulTipurpose Fitting framework 
   - @c sconfig v1.x: configuration file parser module
   - @c slogger v1.x: logging module
   - @c scatalog v1.x: catalog management for astronomical data
     Note that this module also requires the following:
      - @c pyfits v3.1.x (http://www.stsci.edu/institute/software_hardware/pyfits)
      - @c AstroAsciiData v1.1 or greater (http://www.stecf.org/software/PYTHONtools/astroasciidata)

- Optional:
   - @c GalSim v1.x, if the @c gbdsersic @c gbdshsersic or @c gcbdsersic galaxy models are used
   - @c scdm v0.6.x: if the Simple Coordinate Descent Minimizer (SCDM) is used   
   - @c dwtwiener v0.5.x: if uisng the dwtwiener denoising module (Nurbaeva et al., 2011, A&A,
     http://www.aanda.org/articles/aa/full_html/2011/07/aa16556-11/aa16556-11.html)

@section Installation

Download and untar the file <tt>gfit-1.0.0.tar.gz</tt> in the directory of your choice. 
The gfit executables are located in directory gfit-1.0.0/gfit. 

@section Execution

The general syntax is: 

- SMP version: @code gfit_SMP.py [gfit options] @endcode 
- MPI version: @code gfit_MPI.py [gfit options] @endcode

These executables are located in the installation directory under <tt>gfit-0.7.5/gfit</tt>.
By default, gfit will look for a configuration file with default name @c gfit.cfg, located in the
<tt>./gfit/config directory</tt>. These settings can be changed using the @c -c and 
@c -d options (see below).

The supported @c gfit options are:

@verbatim
 -h, --help          Display the usage syntax and the list of supported options.

 -d, --config-dir    Directory where to find the configuration file (default: ./config)

 -c, --config-file   Name of the configuration file (default: gfit.cfg) 
@endverbatim

To run @c gfit_MPI.py, one has to use the mpirun executable as:

@code mpirun -np N gfit_MPI.py [gfit options]@endcode or 
@code mpirun -np N python [python options] gfit_MPI.py [gfit options]@endcode
 
where @c N is the number of processors (or nodes) to use. @c MPI must has been 
installed and configured appropriately. 
On a cluster with a queuing system managed by a job scheduler, @c mpirun will have to be integrated 
in a shell script, invoked using a command like @c qsub. 
   
Since the manager process of gfit uses one processor for itself, workers will 
share <tt>N-1</tt> processors. 

For example: @code mpirun -np 6 gfit_MPI -d mydir -c myconfig.cfg @endcode will run @c gfit
on 6 processors, 1 for the manager and 5 for the workers. The manager and each
of the workers will look for a configuration file: <tt>./mydir/myconfig.cfg </tt>.

@note Before running @c gfit, edit the configuration file to set the @c BASE_INPUT_DIR value, which
      should point to the directory where the @c GREAT3 files reside. 
@note In the configuration file, any environment variable prefixed with @c $ such as @c $HOME will
      be expanded when read.

@section gfit_config Configuration

The behaviour of @c gfit can be controlled though various configuration options located in its
configuration file, <tt>gfit.cfg</tt> by default, in the <tt>gfit/config</tt> directory.

The format of the configuration file is similar to that of MS Windows&trade; .ini text files, with
@c sections embedded within brackets, like <tt>[section]</tt>, each section containing a set of
(@c key, @c value) pairs. The <tt>.ini</tt> format has been extended to support subsections and
section inheritance. Configuration files are managed by the @c sconfig python module, wiich should
be separately installed (see section <b>Prerequisites</b> above).   

@subsection dirs Directories

The input, output and corrresponsing sub-directories are specified under the <tt>[DIR]</tt> section.

@verbatim
# --- This section contains the directories from where gfit will read or store data                                  
[DIR]    

[DIR.INPUT]                            # base input directories: currently unused 
BASE_INPUT_DIR = ./input               # location of input objects to process     

[DIR.OUTPUT]                           # base output directories
BASE_OUTPUT_DIR = ./output             # where to store results, etc.
OUTPUT_PLOT_DIR   = plots              # where to save plots
OUTPUT_MAP_DIR    = maps               # where to save maps
OUTPUT_STAT_DIR   = stats              # where to record statistics
OUTPUT_LOG_DIR    = logs               # where to write logs
OUTPUT_RESULT_DIR = results            # where to write results of runs
OUTPUT_ERROR_DIR  = errors             # where to write errors during runs
@endverbatim

Each time @c gfit is run, a directory of the form: 
<tt>run_DD.MM.YY_HH.MM.SS.branch.obs_type.data_type</TT> is created under the output directory
specified by the key @c BASE_OUTPUT_DIR under the section @c DIR.DIR_OUPUT. The @c branch, 
@c obs_type and @ data_type components are named after the first capitalized letter of the branch,
data type and observation types respectively. 

For instance a run in the control-space-variable 
branch will have its data contained in a directory of the form <tt>run_DD.MM.YY_HH.MM.SS.CSV</TT>. 
Or for a run in both the variable_psf-space constant and the full-ground-variable branches,
a directory like <tt>run_DD.MM.YY_HH.MM.SS.VSC_FGV</TT> will be created. The created directory
will in turn contains the sub-directories @c plots, @c maps, @c stats @c logs, @c results and
@c errors. The @c plots, @c maps, @c stats and @c results will have a directory tree-like structure
that mimics the structure of the dataset.

The <tt>DIR.INPUT</tt> section is currently unused.

@subsection datasets Datasets

The @c [DATASET] section specifies where input images, catalogs an other requred files to process 
are located. It also specifies Unix-like file matching patterns used to identify those files. 
At a minimum, a @c dataset has a name, a base directory and a a list of file matching patterns. 

In this version @c gfit uses 3 datasets: one to access the @c GREAT3 images, catalogs and other 
files, @c SExtractor galaxy and PSF catalog files and a @c PSF dataset with PSF images at the 
positions of galaxies. These datasets are described in more detail below.

@subsubsection g3_dataset GREAT3 dataset

The @c GREAT3 dataset is identified by the <tt>DATASET</tt> section:

@verbatim
# --- Main source dataset, in this case GREAT3
[DATASET]                        # note: cannot process both deep and normal images in the same run

NAME = great3                    # name of dataset to use
BASE_DIR = $HOME/GREAT3          # dataset base input directory

# --- gfit needs galaxy images and catalogs
FILE_PATTERNS = ["image*.fits", "galaxy_catalog*"]  
                                                    
# Note: to process deep images, use: FILE_PATTERNS = ["deep_image*", "deep_galaxy_catalog*"]
 
# --- Search criteria, branches: ["control", "real_galaxy", "variable_psf", "multiepoch", "full"]
BRANCHES   = ["control"]               

# --- Observation type: ["ground", "space"]
OBS_TYPES  = ["ground"]             
                           
# --- Data type: ["constant", "variable"]                                           
DATA_TYPES = ["constant"]  

# --- Set of Image numbers (use [] to process all images (IMAGE_LIST has priority over IMAGE_RANGE) 
IMAGE_LIST  = [0]                   
                                    
# --- [Min,Max] image number range (-1 or [] to disable criterion)       
IMAGE_RANGE = []    

# --- Actual Galaxy stamp size to use for shape measurement
GALAXY_EFFECTIVE_STAMP_SIZE = {"control":{"ground":26, "space":54},\
                                "multiepoch":{"ground":44, "space":44},\
                                "variable_psf":{"ground":44, "space":88},\
                                "real_galaxy":{"ground":44, "space":88},\
                                "full":{"ground":44, "space":88}\
                              }          

# --- Actual PSF stamp size to use for shape measurement
PSF_EFFECTIVE_STAMP_SIZE = {"control":{"ground":26, "space":54},\
                             "multiepoch":{"ground":44, "space":44},\
                             "variable_psf":{"ground":44, "space":88},\
                             "real_galaxy":{"ground":44, "space":88},\
                             "full":{"ground":44, "space":88}\
                           }  
@endverbatim

The <tt>BASE_DIR</tt> key specifies where the base @c GREAT3 directories are located, i.e. where
the branch directories (@c control, @c real_galaxy, @c multiepoch, etc.) are to be read from. 

The <tt>FILE_PATTERNS</tt> key contains a list of pattern matching strings, used to identify which
files are collectively needed to perform shape measurement. A set of such files is processed by 
@c gfit as a @c job. 
For instance, in the context of @c GREAT3, a job  will contain all files related to image @c 000 
in the <tt>control-ground-constant</tt> branch tree: <tt>image-000.0.fits</tt> and
<tt>galaxy_catalog-000.fits</tt>. On a multiprocessor machine or cluster, each job will be assigned 
to one processor.

The <tt>BRANCHES</tt>, <tt>OBS_TYPE</tt> and <tt>DATA_TYPE</tt> keys respectively specify the lists
of branches (e.g. "control"), observation types (e.g. "ground", "space") and data types (e.g. 
"constant", "variable") to process.  

The <tt>IMAGE_LIST</tt> key is a list of image numbers to process, like for instance  
<tt>[0, 10, 67]</tt>. 
An empty list <tt>[]</tt> indicates that @e all images in the underlying branch, observation type 
and data type have to be processed.

The <tt>IMAGE_RANGE</tt> key specifies a range of image numbers to process. for instance, [1, 50] 
will have @c gfit access all images from 1 to 50, 1 and 50 included. Note that the 
<tt>IMAGE_LIST</tt> has priority other <tt>IMAGE_RANGE</tt>: the image range will be ignored unless
<tt>IMAGE_RANGE=[]</tt>.

The <tt>GALAXY_EFFECTIVE_STAMP_SIZE</tt> and <tt>PSF_EFFECTIVE_STAMP_SIZE</tt> keys respectively 
control the sizes of the postage stamps to cut-out from galaxy and PSF images. @c gfit will 
estimate the centroid of the postage stamp image and cut-out a stamp with the specified size.
Note that since most centroids are shifted by 1 or more pixels, cutting-out a square stamp of the 
exact size specified in <tt>GALAXY_EFFECTIVE_STAMP_SIZE</tt> and <tt>PSF_EFFECTIVE_STAMP_SIZE</tt> 
may not be possible. In such cases, @c gfit will automatically reduce the size of the stamp to the
maximum possible size. A warning will be issued in the log file for that image.

@subsubsection se_dataset SEXtractor dataset

The @c SExtractor dataset, <tt>[SEXTRACTOR_DATASET]</tt>, points to the location of 
SExtractor catalogs. These catalogs have been generated beforehand for every galaxy and PSF image
using the @c pse "Parallel SExtractor" utility (M. Gentile, 2011-2013). @c gfit makes use of 
these catalogs to provide guess values of object centroids, fluxes and half-light radii to the 
minimizer prior to fitting. 

@verbatim
[SEXTRACTOR_DATASET]

# --- gfit needs to access galaxy catalogs created by SEXtractor
[SEXTRACTOR_DATASET.GALAXY]         
NAME = gal_sextractor                      
BASE_DIR = $HOME/mg-dev/trunk/pse_package/pse/completed_runs/run_gal_19.11.2013/results
FILE_PATTERNS = ["image*.cat"]   

# Note: to process deep images, use: FILE_PATTERNS = ["deep_image*.cat"]

# --- gfit needs to access star catalogs created by SEXtractor
[SEXTRACTOR_DATASET.PSF]
NAME = psf_sextractor                      
BASE_DIR = $HOME/mg-dev/trunk/pse_package/pse/completed_runs/run_psf_21.11.2013/results
FILE_PATTERNS = ["starfield*.cat"]

# Note: to process deep images, use: FILE_PATTERNS = ["deep_starfield*.cat"]   
@endverbatim

@subsubsection psf_dataset PSF dataset

The @c PSF dataset specifies the location the PSF images interpolated at the positions of 
galaxies. In the case of @c GREAT3 constant PSF branches, it contains only the images of the 
PSF centered at coordinates <tt>(0, 0)</tt> in the @c starfields_image-*-*.fits images (the PSF
that is centered). For the "variable_PSF" branch, it contains images similar to starfield images 
(@c 10000 postage stamps on a grid).  

@verbatim
# --- PSF dataset (PSFs at galaxy positions)
[PSF_DATASET]
NAME = psf 
BASE_DIR = $HOME/mg-dev/trunk/mkpsf_package/mkp/completed_runs/run_21.11.2013/results
FILE_PATTERNS = ["starfield*.fits"]  # or  ["deep_starfield*.cat"] 
#FILE_PATTERNS = ["deep_starfTrueield*.fits"]
@endverbatim

@subsection img_properties Image properties

The <tt>[IMAGE_PROPERTIES]</tt> section details a few important image properties required by 
@c gfit.

@verbatim
[IMAGE_PROPERTIES]

[IMAGE_PROPERTIES.GALAXY]
PIXEL_FLOAT_SIZE = float32             # size of float for pixel values 
IMAGE_HDU_NO = 0                       # HDU no corresponding to the image data
CATALOG_HDU_NO = 1                     # HDU no corresponding to the catalog data
CATALOG_EXTENSION = .fits              # ".txt" or ".fits"

[IMAGE_PROPERTIES.PSF]
PIXEL_FLOAT_SIZE = float32             # size of float for pixel values 
IMAGE_HDU_NO = 0                       # HDU no corresponding to the image data
CATALOG_HDU_NO = 1                     # HDU no corresponding to the catalog data
CATALOG_EXTENSION = .fits              # ".txt" or ".fits"
@endverbatim

The <tt>[CATALOG_EXTENSION]</tt> is required because the @c GREAT3 dataset contains star and galaxy
catalogs in both .txt and .fits format.

@subsection fitting Galaxy Fitting

The main purpose of <tt>[GALAXY_FITTING]</tt> section is to indicate which galaxy model and 
fitting method to use.
@c gfit currently supports the @c scsersic, @c gbdsersic, and @c gcbdsersic galaxy models already 
mentioned.

Two minimizers are available in ths version as mentioned before: @c scileastsq and @c scdmin. 

@verbatim
# --- Galaxy fitting                          
[GALAXY_FITTING]

GALAXY_MODEL_NAME     = gbdsersic      # among: [scsersic, gbdsersic, gbdshsersic, gcbdsersic]
GALAXY_FITTING_METHOD = scdmin         # among: [scileastsq, scdmin]
@endverbatim

@subsection multifit Multifit

The @c Multifit Python module provides @c gfit with a unified interface to galaxy models and 
fitting methods. This makes the @c gfit code independent of a particular galaxy model or minimizer. 

The <tt>[MULTIFIT]</tt> section specifies where the models and fitting methods configuration files
are located. It also indicates from where the corresponding python modules can be dynamically 
loaded. 

@verbatim
# --- Multifit module
[MULTIFIT]
MODEL_CONFIG_PATH   = ./config/models             # path to multifit configuration files for models
METHOD_CONFIG_PATH  = ./config/methods            # path to multifit configuration files for methods 
FITTING_CONFIG_PATH = ./config/methods            # path to multifit configuration files for methods 

MODEL_MODULE_PATH  = $MFIT_DIR/modules/models     # path to multifit modules of fitting models
METHOD_MODULE_PATH = $MFIT_DIR/modules/methods    # path to multifit modules of fitting methods
@endverbatim

@subsection shape Shape measurement

The <tt>[SHAPE_MEASUREMENT]</tt> controls the shape measurement operation.

@verbatim
[SHAPE_MEASUREMENT] 
NORMALIZE_PSF = True                   # if True, normalize the PSF image before convolution
SUBSTRACT_GALAXY_SKY = False           # if True estimate and remove the sky noise in galaxy images
SUBSTRACT_PSF_SKY = False              # if True estimate and remove the sky noise in PSF images
FAILED_ELLIPTICITY_VALUE = 10.0        # if -1,.0 take the value of the fit, even though it failed
GALAXY__SUBSCALE = 1                   # resolution of constructed galaxy models
@endverbatim

The key <tt>FAILED_ELLIPTICITY_VALUE</tt> sets the value of ellipticity used to flag wrong fits in
the results submission files produced for each processed image. The @c GREAT3 presubmission script 
will ignore entries with ellipticity values greater than or equal to @c 10. If set to @c -1, 
@c gfit will not flag entries with bad fits. This is not recommended.

Note that @c gfit also always produce, for each image, a result submission file that 
includes @e all fits either good or bad. So the preferred approach is to choose
<tt>FAILED_ELLIPTICITY_VALUE = 10.0</tt>.

The results files for submission are stored in the directory specified in the 
<tt>OUTPUT_RESULT_DIR</tt> key in the <tt>[DIR.OUTPUT]</tt> section. Their names are of the form: 
<tt>"results_branch_obs_type_data_type-img_no.txt"</tt>, for instance 
<tt>"results_control_ground_constant-000.txt"</tt>. The "full" results file containing all
the fits is prefixed by <tt>full</tt>, like 
<tt>"full_results_control_ground_constant-000.txt"</tt>.

Lastly, @c gfit also produces for each processed image a files containing only the bad fits. These
files are stored in the directory pointed by the <tt>OUTPUT_ERROR_DIR</tt> key in the 
<tt>[DIR.OUTPUT]</tt> section.

@subsection denoising Denoising

The <tt>[DENOISING]</tt> section is intended to include any denoising scheme that may be applied 
to each postage stamp. 

In this version, @c gfit only supports the DWT Wiener wavelet-based denoising algorithm from
Nurbaeva et al., 2011, whose settings can be set in the <tt>[DENOISING.DWTWIENER]</tt> section.

@verbatim
[DENOISING]

[DENOISING.DWTWIENER]                  # Denoising with dwtwiener (Nurbaeva et al., A&A, 2011)

[DENOISING.DWTWIENER.GALAXY]
APPLY_DWTWIENER = False                # if True apply dwtwiener to denoise the images
DECOMPOSITION_LEVEL = 2                # decomposition level (from l=1, default l=2)
CYCLE_SPINNING = False                 # if True, use cycle spinning
SIGMA_NOISE =-1                        # noise standard deviation (estimated if set to -1)
DENOISING_STRENGTH = 0.6               # denoising strength fron 0 to 1
@endverbatim
 
@subsection collect Data collection

The <tt>[DATA_COLLECTION]</tt> section relates to extra data that can e.g. be collected to help 
diagnose problems, improve gfit's performance, perform training for bias correction, etc.  

@verbatim
[DATA_COLLECTION]

COLLECT_DATA = True                    # if True collect shape measurement data, besides parameters
DUMP_COLLECTED_DATA = True             # if True, create a catalog with collected data per image
PLOT_COLLECTED_DATA = True             # it True, plot collected data per image
@endverbatim

@c gfit will always record fitted parameters values and will allow to dump, plot or map these.

Additionally, if <tt>COLLECT_ERROR_DATA</tt> is set to @c True, @c gfit will collect additional data 
besides fitted parameters, such as @c SExtractor data (from SExtrator catalogs) and statistics on 
residuals.

If <tt>PLOT_COLLECTED_DATA</tt> is set to @c True, fitted parameter data and possibly extra data 
collected if <tt>COLLECT_DATA=True</tt> will be plotted and stored in the directory specified by
key <tt>OUTPUT_PLOT_DIR</tt> in section <tt>[DIR.OUTPUT]</tt>. Otherwise, no plot will be made.

If <tt>DUMP_COLLECTED_DATA</tt> is set to true, fitted parameter data and extra data collected
if <tt>COLLECT_DATA=True</tt> will dumped to a file with name of the form
"data_branch_obs_type_data_type-img_no.txt", like for example 
"data_control_ground_constant-000.txt". The file will be stored in the directory specified in the 
<tt>OUTPUT_RESULT_DIR</tt> key in the <tt>[DIR.OUTPUT]</tt> section. Otherwise, no file will be
created.
  
@subsection filtering Object Filtering

This section provides a way to filter-out of "clean" unwanted objects, i.e.
galaxies or stars. A filtering criterion takes the form of a range of allowed values, like 
<tt>[0, 600]</tt>. The value @c None denotes that there is no lower or upper bound. For instance,
<tt>[0, None]</tt> allows for values greater than zero, with no upper limit.

Filtering takes place <em>before</em> fitting ("pre-filtering") or <em>after</em> 
("post-filtering"). 
- In this version, pre-filtering is applied on selected SExtractor catalog keywords. The
name of the criterion is the name of the SExtractor keyword prefixed by @c "SE_", like 
@c SE_FLUX_RADIUS.
- Post-filtering can currently be done on the sum of the residuals between observed and built model, 
normalized by the total number of pixels in the postage stamp used for fitting.

The keyword @c FILTERING_POLICY specifies what should be done with the "filtered" objects. In this 
version, objects are effectively "removed" from the final object catalog 
(<tt>FILTERING_POLICY = remove</tt>). A future version may allow to ascribe a weighting factor
to filtered objects instead of removing them.

Filtering can be enabled or disabled by specifying <tt>ENABLE_FILTERING=True</tt> or 
<tt>ENABLE_FILTERING=False</tt>. If filtering is enabled, a catalog with un-filtered objects will 
still be produced in order to better investigate the effect of filtering.

@verbatim
[OBJECT_FILTERING]

# --- Galaxy ---
[OBJECT_FILTERING.GALAXY]

ENABLE_FILTERING = True    # set to True to enable galaxy filtering
FILTERING_POLICY = remove  # how to handle filtered data, among: [remove, weight]

FILTER = {"SE_HLR" : [0, None],\     # allowed SExtractor FLUX_RADIUS range
          "SE_FWHM_IMAGE"  : [0, None],\     # allowed SExtractor FWHM_IMAGE range
          "SE_FLUX"        : [0, None],\     # allowed SExtractor FLUX_BEST range
          "SE_MAG"         : [None, None],\  # allowed SExtractor MAG_BEST range
          "SE_SNR"         : [0, None],\     # allowed SExtractor SNR (FLUX_BEST/FLUXERR_BEST)
          "NORM_RESIDUALS" : [0, None] \     # allowed range of normalized sum of residuals
         }         

# --- PSF ---
[OBJECT_FILTERING.PSF]

ENABLE_FILTERING = False   # set to True to enable PSF filtering
FILTERING_POLICY = remove  # how to handle filtered data, among: [remove, weight]

FILTER = {"SE_HLR"         : [0, None],\     # allowed SExtractor FLUX_RADIUS range
          "SE_FWHM_IMAGE"  : [0, None],\     # allowed SExtractor FWHM_IMAGE range
          "SE_FLUX"        : [0, None],\     # allowed SExtractor FLUX_BEST range
          "SE_MAG"         : [None, None],\  # allowed SExtractor MAG_BEST range
          "SE_SNR"         : [0, None],\     # allowed SExtractor SNR (FLUX_BEST/FLUXERR_BEST)
          "NORM_RESIDUALS" : [0, None] \     # allowed range of normalized sum of residuals
         } 
@endverbatim

@subsection logging Logging

The <tt>[LOGGING]</tt> section tells if logging is enabled, and which logging mode is used. 

@verbatim
[LOGGING]
MASTER_LOGGING_ENABLED = True          # quiet mode if set to False
CREATE_WORKER_LOGGERS = per_job        # among: {none, per_job, per_worker}
@endverbatim

If the @c MASTER_LOGGING_ENABLED key value is set to False, no display output nor file logging will 
be done.

In the case <tt>MASTER_LOGGING_ENABLED=True</tt>, if the @c CREATE_WORKER_LOGGERS key value is set
to @c none, then no log file will be created for workers. 

If not @c none, the @c CREATE_WORKER_LOGGERS key specifies two logging modes: @c per_job and 
@c per_worker. In the @c per_job mode, one log file will be created for each job processed (i.e.
each @c GREAT3 image). In the @c per_worker mode, one log file will be created for each worker, not
for each image, so the worker's log will trace the shear measurement output of all the images 
processed by a particular worker. 

@subsection debugging Debugging

The <tt>[DEBUGGING]</tt> section controls the amount and type of debugging information.

@verbatim
[DEBUGGING]

MAX_NB_STAMPS = 5                       # maximum number of stamps to process (-1: all) 

EXTRACT_GALAXY_STAMPS = True            # if True, write individual galaxy stamps as .FITS
EXTRACT_CONVOLVED_GALAXY_STAMPS = False # if true, write individual convolved galaxy stamps as .FITS
EXTRACT_GALAXY_RESIDUALS = False        # if True, extract galaxy residuals after fitting
EXTRACT_PSF_STAMPS = False              # if True, write individual PSF stamps as .FITS files
CREATE_CONVOLVED_GALAXY_MOSAIC = True   # if True, create a convolved stamp mosaic
CREATE_RESIDUALS_GALAXY_MOSAIC = True   # if True, create a convolved stamp mosaic
MOSAIC_STAMP_SIZE = 48                  # stamp size in generated mosaic

PLOT_GALAXY_STAMPS = True               # if True, plot individual galaxy stamps
PLOT_CONVOLVED_GALAXY_STAMPS = False    # if True, plot individual convolved galaxy stamps 
PLOT_GALAXY_RESIDUALS = False           # if True, plot galaxy residuals after fitting
PLOT_PSF_STAMPS = False                 # True to plot individual PSF stamps
PLOT_ENABLE_3D = True                   # if True, also plot in 3D

RECORD_FAILED_FITS = True               # Keep a copy of the galaxy stamp in the errors directory

PROGRESS_COUNT = 5                      # display progress every <PROGRESS_COUNT> galaxies
@endverbatim

The number of postage stamps to fit within an image is set by the <tt>MAX_NB_STAMPS</tt>. A value
of @c -1 means that @e all postage stamps will be processed. otherwise, @c gfit will process the
@e first number of stamps specified in <tt>MAX_NB_STAMPS</tt>. The sorting order is y then x, so 
@c gfit will proceed row after row in an image, processing each column of a given row.  

The @c EXTRACT_GALAXY_STAMPS, @c EXTRACT_CONVOLVED_GALAXY_STAMPS, @c EXTRACT_GALAXY_RESIDUALS, 
@c EXTRACT_PSF_STAMPS boolean flags tell whether postage stamps should be extracted as .FITS
files for observed galaxies, convolved galaxies, galaxy residuals and observed PSFs, in this order.
The target directory is that specified as <tt>OUTPUT_RESULT_DIR</tt> in section 
<tt>[DIR.OUTPUT]</tt>.

If @c CREATE_CONVOLVED_GALAXY_MOSAIC is @c True, @c gfit will create a mosaic containing the 
best-fitted convolved galaxy models inside postage stamps of size set by @c MOSAIC_STAMP_SIZE. The
coordinate of galaxies within the mosaic will be the same as in the source image and catalog. Failed
fits are indicated by a postage stamp with zero pixel values at the corresponding coordinates.<br>
Similarly, if @c CREATE_RESIDUALS_GALAXY_MOSAIC is set to @c True, a mosaic of the residuals 
postage stamps will also be created, where residual postage stamps are the difference between 
observed galaxies and model-convolved galaxies.<br>
Note that, in practice, most postage stamps will have an effective size smaller than the specified mosaic
stamp size (e.g. 44 instead of <tt>MOSAIC_STAMP_SIZE = 48) </tt>. Such stamps will be aligned on 
the bottom-left corner of the mosaic stamp.

The @c PLOT_GALAXY_STAMPS, @c PLOT_CONVOLVED_GALAXY_STAMPS, @c PLOT_GALAXY_RESIDUALS and
@c PLOT_PSF_STAMPS specify whether postage stamps should be plotted for observed galaxies, 
convolved galaxies, galaxy residuals and observed PSFs, in this order.
The target directory is that specified as <tt>OUTPUT_PLOT_DIR</tt> in section 
<tt>[DIR.OUTPUT]</tt>. If @c PLOT_ENABLE_3D is set to @c True, three-dimensional plots are also 
drawn.

The @c RECORD_FAILED_FITS flag indicates if postrage stamps for which the fit failed should be
stored for later investigation in the directory pointed by <tt>OUTPUT_ERROR_DIR</tt> in section 
<tt>[DIR.OUTPUT]</tt>.

@c PROGRESS_COUNT specifies how often @c gfit should display or log output information. For 
instance, a value of @c 10 will output every time @c 10 postage stamp images have been fitted.

@subsection multiproc SMP Multiprocessing

The <tt>[MULTIPROC_SMP]</tt> only concerns the @c SMP version of @c gfit and tells how many processors
should be allocated through the @c WORKER_PROCESSES_COUNT key. A value of @c -1 let the system
decide the optimal number of processors to use. 

@verbatim
[MULTIPROC_SMP]
WORKER_PROCESSES_COUNT = 0   # nb. worker processes for SMP (0 automatic detection, > 0 otherwise)
@endverbatim

@section multifit_config Multifit configuration

@c gfit interacts with the Multifit module through a set of configuration files. There are 
configuration files for fitting methods and for parametric models.

@subsection model_config Galaxy Model configuration

A Multifit configuration file is present for each supported galaxy model and bear the same name as 
that of the model. Thus, there are @c 3 model configuration files in this version of @c gfit:
- <tt>scsersic.cfg</tt> for the single Sérsic component model
- <tt>gbdsersic.cfg</tt> for the bulge+disk model
- <tt>gcbdsersic.cfg</tt> for the "complex" bulge+disk model

For each galaxy model, the list of model parameters is given, specifing the name of each parameter,
its guess value and its bounds (min and max alowed values). 

@subsection method_config Fitting Method configuration

A Multifit configuration file is also present for each supported fitting method and bear the same
name as that of the method. Since @c 2 minimizers are supported in this version, there are also 
@c 2 configuration files:
- <tt>scdmin.cfg</tt> for the @c SCDM minimizer
- <tt>scileastsq.cfg</tt> for the <tt>Levenberg-Marquardt</tt> non-linear least-squares 
  implementation available with @c SciPy.  

Sometimes, it is mecessary to provide fitting information that is common to both the model and
the method. This is the case for the @c SCDM minimizer, which uses such information for optimizing
the fitting process. For instance one can specify the initial step value to fit a given parameter.
Such extra fitting information can be kept in a file prefixed by the minimizer name and suffixed by
@c "_fitting". In the case of SCDM, known as @c scdmin by @c gfit, the fitting configuration file 
is  <tt>scdmin_fitting.cfg</tt>. It contains extra fitting information for each galaxy model.

"""

from gfit import *
from gfit.gfit_version import __version__

