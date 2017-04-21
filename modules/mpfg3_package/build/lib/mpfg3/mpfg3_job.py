"""! 
   @package mpfg3.mpfg3_job Job Management
   @author Marc Gentile
   @file mpfg3_job.py
   Job Management
""" 

# -- Python imports
import os, sys
import time

# -- External imports
from mpfx.mpfx_job import *         # base job processing

# --- Module-specific imports
from mpfg3_data import *            # dataset management

# --------------------------------------------------------------------------------------------------
class Mpfg3JobProcessor(MpfxJobProcessor):
   
   """! 
      Job processor: submit jobs and process associated job results. Based on mpfx.MpfxJobProcessor.   
   """

   def __init__(self, master):
      """ JobProcessor default constuctor """

      MpfxJobProcessor.__init__(self, master)   
      
   # ~~~~~~~~~~~~~~~
   # Public methods 
   # ~~~~~~~~~~~~~~~

   # -----------------------------------------------------------------------------------------------
   def preprocess_job(self, job, worker):
      """! Invoked by worker to perform some optional pre-processing on the job. It may be for 
         instance some format conversion. or other preparation step before the actual
         processing takes place in process_job().

         @param job an object of class Job to pre-process
         @param worker instance of the worker process

         @note this method is called <b>before</b> process_job_result()
         @note this method is to be overridden by subclasses to perform
               any processing required <b>before</b> process_job() is called
         @see process_job, postprocess_job, Job
      """
      pass

   # -----------------------------------------------------------------------------------------------
   def process_job(self, job, worker):
      """!
         Process a job and return the corresponding results to the Master (here sample data). 

         @param job an object of class Job to process
         @param worker instance of the worker process
         @return an object of class JobResult containing the data of the processed job
         @note this method should be overridden in a subclass of mpf 
         @see preprocess_job, postprocess_job, Job, JobResult
      """

      if worker.logging_enabled():
         worker.logger.log_info_p("{0} - job {1} - File types included: {2}".format(
                                                              worker, job, job.get_file_types()))

      data_dico = {"result": "result from job: {0}".format(job.name)}    # sample result data

      time.sleep(1)     # simulate processing duration

      return Mpfg3JobResult(worker, job, data_dico)  # store the job results

   # -----------------------------------------------------------------------------------------------
   def postprocess_job(self, job_result, worker):
      """! 
         Invoked by the Worker to perform some optional post-processing on the job.

         @param job_result object of class JobResult with processed data
         @param worker instance of the worker process

         @note this method is called <b>after</b> process_job_result()
         @note this method is to be overridden by subclasses to perform
               any processing required <b>after</b> process_job() is called
         @see process_job, postprocess_job, JobResult
      """
      pass

   # -----------------------------------------------------------------------------------------------
   def process_job_result(self, job_result, master):
      """! 
         Process the result associated with a processed job, which may consist, for instance, in 
         storing files to the output result directory tree, draw plots, compute statistics... 

         @param job_result object of class JobResult with processed data
         @param master Master object instance
         @note This method should be overridden in a subclass of mpf  
         @see Job, JobResult 
      """

      if master.logging_enabled():
         master.logger.log_info_p("{0} - Processing result from job {1}...".format(
                                                                 master, job_result.job))

   # -----------------------------------------------------------------------------------------------
   def create_jobs(self, master):
      """! 
         Locate all objects to process and create the corresponding jobs. 
         @param master Master object instance
         @return the list of created jobs 
      """

      return MpfxJobProcessor.create_jobs(self, master)

   # -----------------------------------------------------------------------------------------------
   def create_job(self, master, dataset, *args):
      """! 
         Factory method for creating a Job object.
         Each subclass can define a job class derived from mpf.Job and instanciate it here.
         @param master the master process
         @param dataset data source used to create the job
         @param args a list of arguments chosen by the caller and whose nature depends on the 
                data necessary to process the job
         @return a new Mpfg3Job job object
         @see Job, Mpfg3Job
      """

      return Mpfg3Job(master, dataset, *args)
      
   # -----------------------------------------------------------------------------------------------
   def create_dataset(self, master, dataset_name, dataset_type, 
                                    dataset_base_dir, dataset_dir_list, dataset_recurse):
      """! 
         Create a Dataset object that represent the data source for images and catalogs.
         @param master master process instance
         @param dataset_name dataset name
         @param dataset_type prefix of a Dataset class, assumed to be of the form 
                @code <dataset_type>Dataset @endcode, like @c MpfxDataset
         @param dataset_base_dir dataset base directory
         @param dataset_dir_list [optional] a list of specific dirs to search under base directory
         @param dataset_recurse  [optional] tell whether to walk down directories (default @c True)
         @return the Dataset instance 
      """

      return Mpfg3Dataset(master, dataset_name, 
                                  dataset_base_dir, dataset_dir_list, dataset_recurse)

   # -----------------------------------------------------------------------------------------------
   def record_job_result(self, job_result, master):
      """! 
         Add a job result to the list of job results. This method is invoked by the Master. 
         @param master instance of the Master
         @param job_result object of class JobResult containing the data of the processed job  
         @see Job, JobResult       
      """

      MpfxJobProcessor.record_job_result(self, job_result, master)

   # -----------------------------------------------------------------------------------------------
   def all_jobs_processed(self, master):
      """! 
         This method is called by the Master once all the jobs have been processed. 
         @param master instance of the Master
         @note Job results can be selectively obtained by calling the get_job_results() methods,
               specifying the relevant query criterion: branch, observation type, data type, image 
         @note the entire list of Job results can aldo be obtained with by calling 
               JobProcessor.job_result_list()
         @see Job, JobResult, get_job_result()
      """

      MpfxJobProcessor.all_jobs_processed(self, master)

   # -----------------------------------------------------------------------------------------------
   def get_job_results(self, **kwargs):
      """! 
         Return the list of Mpfg3JobResult objects for a specified branch, observation type, 
         data_type and image number. These parameters must be specified using the 
         "branch=", "img_no=" and "epoch=" keywords
         @return list of Mpfg3JobResult objects
         @code
            Examples
            - self.get_job_results(branch="control/space") to retrieve job results for all images 
              and epochs in branch "control/space"
            - self.get_job_results(branch="control/space/variable", img_no=1) to retrieve 
              job results of branch "control/space/variable", image number 1, all epochs 
            - self.get_job_results(branch="control/space/variable", img_no=1, epoch=2) to retrieve 
              job results of branch "control/space/variable", image number 1, epoch index 2 
         @endcode
         @note Mpfg3JobResult objects can also be directly obtained from the 
               mpfg3.mpfg3_job.Mpfg3JobProcessor.job_result_dico dictionary
      """

      return MpfxJobProcessor.get_job_results(self, **kwargs)


   # ~~~~~~~~~~~~~~~
   # Private methods 
   # ~~~~~~~~~~~~~~~



# -------------------------------------------------------------------------------------------------
class Mpfg3Job(MpfxJob):      

   """! 
       Represents a Mpfg3 job 
   """

   def __init__(self, master, dataset, img_path_dico, *args):
      """! 
         Construct a new job for identifying a set of files in the dataset to process together.

         To obtain the path associated with a given type, just use:

         @code self.img_path_dico[type]@endcode with one of the aforementioned types.
         @param master instance of the Master
         @param dataset data source used to create the job
         @param img_path_dico file absolute path dictionary  for every file types referenced in 
                the job
         @param args a list with the branch tree, image number and epoch number

         @return an initialized Mpfg3Job object 
      """

      MpfxJob.__init__(self, master, dataset, img_path_dico, *args)
      
   # ~~~~~~~~~~~~~~
   # Public methods 
   # ~~~~~~~~~~~~~~

   # -----------------------------------------------------------------------------------------------
   def get_catalog_from_image(self, img_name, cat_extension=".txt"):
      """! 
         Convenience method to find the catalog name matching an image name, for instance
         "starfield_image" -> "star_catalog", or  "image.fits" -> "galaxy_catalog"  
         @param img_name name or image file type of an image 
         @param cat_extension expected file extension for the catalog (default: .txt) 
         @return a catalog name among:
         - "deep_star_catalog"
         - "star_catalog"
         - "deep_galaxy_catalog"
         - "galaxy_catalog"
         @see see get_file_types. __init__
         @note: the same catalog filename may have more than one extension (like .txt and .fits)
      """
      cat_name = None   

      if "deep_starfield_image" in img_name:
         cat_name = "deep_star_catalog"         # a deep star catalog
      elif "starfield_image" in img_name:
         cat_name = "star_catalog"              # a star catalog
      elif "deep_image" in img_name:
         cat_name = "deep_galaxy_catalog"       # a deep galaxy catalog
      elif "image" in img_name:                 
         cat_name = "galaxy_catalog"            # a galaxy catalog

      if not cat_name is None:
         cat_name += cat_extension

      return cat_name

   # -----------------------------------------------------------------------------------------------
   def get_image_type_from_catalog(self, cat_name):
      """!
         Convenience method to find the image matching a given catalog file type, 
         for instance
         "galaxy_catalog.txt" -> "image.fits" or "galaxy_catalog.fits" -> "image.fits"
         @param cat_name name of an catalog type
         @return an image file type among:
         - "deep_starfield_image.fits"
         - "starfield_image.fits"
         - "deep_image.fits"
         - "image.fits"
         @see see get_file_types. __init__
      """  

      img_name = None

      if "deep_star_catalog" in cat_name:
         img_name = "deep_starfield_image.fits"         # a deep star image
      elif "star_catalog" in cat_name:
         img_name = "starfield_image.fits"              # a star image
      elif "deep_galaxy_catalog" in cat_name:
         img_name = "deep_image.fits"                   # a deep galaxy image
      elif "galaxy_catalog" in cat_name:                 
         img_name = "image.fits"                        # a galaxy image

      return img_name

   # -----------------------------------------------------------------------------------------------
   def get_stamp_size(self):
      """! 
         Return the image postage stamp size based on the branch
         @return postage stamp size for the job's branch  
      """    
      #print "branch:", self.branch, "ground in branch:",  "ground" in self.branch 
      
      if self.branch != "multiepoch" and self.branch != "full":
         if "ground" in self.branch:
            return 48
         elif "space" in self.branch:
            return 96
         else:
            return 48

      elif self.branch == "multiepoch":
         # Assumes G3 space images are coadded with size 96
         if "ground" in self.branch:
            return 48
         elif "space" in self.branch:
            return 96
         else:
            return 48

      elif self.branch == "full":
         if "ground" in self.branch:
            return 48
         elif "space" in self.branch:
            return 96
         else:
            return 48

# -------------------------------------------------------------------------------------------------
class Mpfg3JobResult(MpfxJobResult):      

   """! 
       Represents the result of a GREAT3 job.
   """

   def __init__(self, worker, job, result):   
      """! Construct a Mpfg3JobResult job result object """

      MpfxJobResult.__init__(self, worker, job, result)




# --------------------------------------------------------------------------------------------------
class Mpfg3MultiEpochJobProcessor(MpfxMultiEpochJobProcessor):
   
   """! 
      Job processor: submit jobs and process associated job results, taking care
      of multiepoch (multiple exposure) images and related catalogs. 
      Based on mpfx.MpfxJobProcessor.   
   """

   def __init__(self, master):
      """ JobProcessor default constuctor """

      MpfxMultiEpochJobProcessor.__init__(self, master)   
      
   # ~~~~~~~~~~~~~~~
   # Public methods 
   # ~~~~~~~~~~~~~~~

   # -----------------------------------------------------------------------------------------------
   def create_jobs(self, master):
      """! 
         Locate all objects to process and create the corresponding jobs. 
         @param master Master object instance
         @return the list of created jobs 
      """
   
      return MpfxMultiEpochJobProcessor.create_jobs(self, master)

   # -----------------------------------------------------------------------------------------------
   def create_job(self, master, dataset, *args):
      """! 
         Factory method for creating a Job object.
         Each subclass can define a job class derived from mpf.Job and instanciate it here.
         @param master the master process
         @param dataset data source used to create the job
         @param args a list of arguments chosen by the caller and whose nature depends on the 
                data necessary to process the job
         @return a new Mpfg3Job job object
         @see Job, Mpfg3Job
      """

      return Mpfg3MultiEpochJob(master, dataset, *args)  # create a Mpfg3Job instance, based on parent class mpf.Job
      
   # -----------------------------------------------------------------------------------------------
   def create_dataset(self, master, dataset_name, dataset_type, 
                                    dataset_base_dir, dataset_dir_list, dataset_recurse):
      """! 
         Create a Dataset object that represent the data source for images and catalogs.
         @param master master process instance
         @param dataset_name dataset name
         @param dataset_type prefix of a Dataset class, assumed to be of the form 
                @code <dataset_type>Dataset @endcode, like @c MpfxDataset
         @param dataset_base_dir dataset base directory
         @param dataset_dir_list [optional] a list of specific dirs to search under base directory
         @param dataset_recurse  [optional] tell whether to walk down directories (default @c True)
         @return the Dataset instance 
      """

      return Mpfg3MultiEpochDataset(master, dataset_name, 
                                    dataset_base_dir, dataset_dir_list, dataset_recurse)

   # -----------------------------------------------------------------------------------------------
   def process_job(self, job, worker):
      """!
         Process a job and return the corresponding results to the Master (here sample data). 

         @param job an object of class Job to process
         @param worker instance of the worker process
         @return an object of class JobResult containing the data of the processed job
         @note this method should be overridden in a subclass of mpf 
         @see preprocess_job, postprocess_job, Job, JobResult
      """

      if worker.logging_enabled():
         worker.logger.log_info_p("{0} - job {1} - File types included: {2}".format(
                                                              worker, job, job.get_file_types()))

      data_dico = {"result": "result from Mpfg3MultiEpochJob: {0}".format(job.name)}    # sample result data

      # simulate processing duration
      time.sleep(1)

      return Mpfg3MultiEpochJobResult(worker, job, data_dico)  # store the job results

   # -----------------------------------------------------------------------------------------------
   def process_job_result(self, job_result, master):
      """! 
         Process the result associated with a processed job, which may consist, for instance, in 
         storing files to the output result directory tree, draw plots, compute statistics... 

         @param job_result object of class JobResult with processed data
         @param master Master object instance
         @note This method should be overridden in a subclass of mpf  
         @see Job, JobResult 
      """

      if master.logging_enabled():
         master.logger.log_info_p("{0} - Processing result from job {1}...".format(
                                                                 master, job_result.job))

   # -----------------------------------------------------------------------------------------------
   def record_job_result(self, job_result, master):
      """! 
         Add a job result to the list of job results. This method is invoked by the Master. 
         @param master instance of the Master
         @param job_result object of class JobResult containing the data of the processed job  
         @see Job, JobResult       
      """

      MpfxMultiEpochJobProcessor.record_job_result(self, job_result, master)

   # -----------------------------------------------------------------------------------------------
   def all_jobs_processed(self, master):
      """! 
         This method is called by the Master once all the jobs have been processed. 
         @param master instance of the Master
         @note Job results can be selectively obtained by calling the get_job_results() methods,
               specifying the relevant query criterion: branch, observation type, data type, image 
         @note the entire list of Job results can aldo be obtained with by calling 
               JobProcessor.job_result_list()
         @see Job, JobResult, get_job_result()
      """

      MpfxMultiEpochJobProcessor.all_jobs_processed(self, master)

   # -----------------------------------------------------------------------------------------------
   def get_job_results(self, **kwargs):
      """! 
         Return the list of Mpfg3JobResult objects for a specified branch, observation type, 
         data_type and image number. These parameters must be specified using the 
         "branch=", "img_no=" and "epoch=" keywords
         @return list of Mpfg3JobResult objects
         @code
            Examples
            - self.get_job_results(branch="control/space") to retrieve job results for all images 
              and epochs in branch "control/space"
            - self.get_job_results(branch="control/space/variable", img_no=1) to retrieve 
              job results of branch "control/space/variable", image number 1, all epochs 
            - self.get_job_results(branch="control/space/variable", img_no=1, epoch=2) to retrieve 
              job results of branch "control/space/variable", image number 1, epoch index 2 
         @endcode
         @note Mpfg3JobResult objects can also be directly obtained from the 
               mpfg3.mpfg3_job.Mpfg3JobProcessor.job_result_dico dictionary
      """

      return MpfxMultiEpochJobProcessor.get_job_results(self, **kwargs)



# -------------------------------------------------------------------------------------------------
class Mpfg3MultiEpochJob(MpfxMultiEpochJob):      

   """! 
       Represents a Mpfg3MultiEpochJob job, an extension of a Mpg3Job for holding a group of
       images with multiple exposures, with related catalogs. 
   """

   def __init__(self, master, dataset, img_path_dico, *args):
      """! 
         Construct a new job for identifiying a set of files in the dataset to process together.

         To obtain the path associated with a given type, just use:

         @code self.img_path_dico[type]@endcode with one of the aforementioned types.
         @param master instance of the Master
         @param dataset data source used to create the job
         @param img_path_dico file absolute path dictionary  for every file types referenced in 
                the job
         @param args a list with the branch tree, image number and epoch number

         @return an initialized Mpfg3MultiEpochJob object 
      """
 
      # Parent class initialization      
      MpfxMultiEpochJob.__init__(self, master, dataset, img_path_dico, *args)

# -------------------------------------------------------------------------------------------------
class Mpfg3MultiEpochJobResult(Mpfg3JobResult):      

   """! 
       Represents the result of a job with multiple epochs (exposures) images and related catalogs.
   """

   def __init__(self, worker, job, result):   
      """! Construct a MpfxJobResult job result object """
      Mpfg3JobResult.__init__(self, worker, job, result)


# -- EOF mpfg3_job.py
