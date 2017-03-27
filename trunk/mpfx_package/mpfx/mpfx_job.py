"""!
   @package mpfx.mpfx_job Job Management
   @author Marc Gentile
   @file mpfx_job.py
   Job Management
"""

# -- Python imports
import os, sys
import re
import time
from operator import itemgetter, attrgetter
from itertools import chain

# -- External imports
from mpfg.mp_job import *         # base job processing
from mpfg.mp_data import *        # data access

# --- Module-specific imports
from mpfx_data import *           # dataset management

# --------------------------------------------------------------------------------------------------
class MpfxJobProcessor(JobProcessor):

   """!
      Job processor: submit jobs and process associated job results. Based on mpf.JobProcessor.
   """

   def __init__(self, master):
      """ MpfxJobProcessor default constuctor """
      JobProcessor.__init__(self, master)

      self._job_result_dico = {}          # dictionary holding all job results
      self._job_branches = None           # all branches to process (initialized at jobs' creation)

   # ~~~~~~~~~~
   # Properties
   # ~~~~~~~~~~

   @property
   def job_result_dico(self):
      """! @return the results of the processed jobs as a dictionary """
      return self._job_result_dico

   @property
   def job_branches(self):
      """! @return a list of all the branches (path trees) to process """
      return self._job_branches

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

      data_dico = {"result": "result from MpfxJob: {0}".format(job.name)}    # sample result data

      # simulate processing duration
      time.sleep(1)

      return MpfxJobResult(worker, job, data_dico)  # store the job results

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
         @param master master process instance
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
         @param master Master process instance
         @return the list of created jobs
      """

      # --- Create primary source Dataset
      dataset_name = master.config.get_as_string("NAME", "PRIMARY_DATASET")
      if master.config.has_key("TYPE", "PRIMARY_DATASET"):
         dataset_type = master.config.get_as_string("TYPE", "PRIMARY_DATASET")
      else:
         dataset_type = "Mpfx"      # default Dataset class prefix
      dataset_base_dir = master.config.get_as_string("BASE_DIR", "PRIMARY_DATASET")
      if master.config.has_key("DIR_INPUT_LIST", "PRIMARY_DATASET"):
         dataset_dir_list = master.config.get_as_list("DIR_INPUT_LIST", "PRIMARY_DATASET")
      else:
         dataset_dir_list = []
      dataset_image_list = master.config.get_as_list("IMAGE_LIST", "PRIMARY_DATASET")
      dataset_image_range = master.config.get_as_list("IMAGE_RANGE", "PRIMARY_DATASET")
      dataset_file_patterns = master.config.get_as_list("FILE_PATTERNS", "PRIMARY_DATASET")
      if master.config.has_key("DIR_RECURSE", "PRIMARY_DATASET"):
         dataset_recurse_dirs = master.config.get_as_boolean("DIR_RECURSE", "PRIMARY_DATASET")
      else:
         dataset_recurse_dirs = True

      dataset = self.create_dataset(master, dataset_name, dataset_type,
                                            dataset_base_dir, dataset_dir_list,
                                            dataset_recurse_dirs)
      if dataset is None:
         if master.logging_enabled():
            master.logger.log_error_p("{0} - Dataset {1} could not be instantiated...".format(
                                                                master, dataset_name))

         # --- Dataset could not be instantiated, return no jobs
         self._job_list = []
         self._job_branches = []
         return self._job_list

      # --- Search for matching files to process

      if master.logging_enabled():
         master.logger.log_info_p("{0} - Looking for files to process in directory: {1}...".format(
                                                             master, dataset_base_dir))
         if len(dataset_dir_list) > 0:
            master.logger.log_info_p("{0} - Input dir(s): {1}".format(master, dataset_dir_list))
         else:
            master.logger.log_info_p("{0} - Input dir(s): {1}".format(master, dataset_base_dir))

         master.logger.log_info_p("{0} - Input files: {1}".format(master, dataset_file_patterns))

      filepath_dico = dataset.query(master, dataset_file_patterns,
                                    dir_list=dataset_dir_list,
                                    image_list=dataset_image_list,
                                    image_range=dataset_image_range,
                                    sort=True, recurse=dataset_recurse_dirs)

      # --- For the multiepoch branch, make sure the catalog paths are also included
      #     for images with non-zero epoch indice
      for branch_key in filepath_dico:

         if len(filepath_dico[branch_key].keys()) < 2:

            # --- Extract catalog file_types and paths from image with 0 epoch index
            #     and insert them in filepath_dico[branch_key] where they are missing

            new_branch_key = branch_key[:-1] + tuple([0])

            branch_items = filepath_dico[new_branch_key].items()
            for (file_type, file_path) in branch_items:

               if dataset.is_catalog(master, file_type):
                  filepath_dico[branch_key][file_type] = file_path

      # --- Create a Job object for each image <img_no> in the relevant branch structure,
      #     with the filepaths for each type of file (images, catalogs...)

      self._job_list = [self.create_job(master, dataset, filepath_dico, *key) \
                        for key in sorted(filepath_dico.keys())]

      # --- Record the list of branches to which jobs belongs
      self._job_branches = list(set([j.branch for j in self._job_list]))

      # --- Create output directories for job outputs:
      output_dir_dico = master.config.get_section_data("DIR.OUTPUT")
      for (dir_key, dir_value) in output_dir_dico.items():
         if dir_key.startswith("OUTPUT_"):
            if dir_key.find("LOG") == -1:  # no tree for log dirs
               base_dir = os.path.join(master.run_output_dir, dir_value)
               #self.helper.make_dir(base_dir)
               self._create_dir_tree(base_dir, self._job_branches)

      return sorted(self._job_list, key=attrgetter('img_no'))   # sort by img_no

   # -----------------------------------------------------------------------------------------------
   def create_job(self, master, dataset, *args):

      """!
         Factory method for creating a Job object.
         Each subclass can define a job class derived from mpf.Job and instanciate it here.
         @param master master process instance
         @param dataset data source used to create the job
         @param args a list of arguments chosen by the caller and whose nature depends on the
                data necessary to process the job
         @return a new MpfxJob job object
         @see Job, MpfxJob
      """

      return MpfxJob(master, dataset, *args)  # create a MpfxJob instance, based on parent class mpf.Job

   # -----------------------------------------------------------------------------------------------
   def create_dataset(self, master, dataset_name, dataset_type, dataset_base_dir,
                                    dataset_dir_list, dataset_recurse_dirs):
      """!
         Create a Dataset object that represent the data source for images and catalogs.
         @param master master process instance
         @param dataset_name dataset name
         @param dataset_type prefix of a Dataset class, assumed to be of the form
                @code <dataset_type>Dataset @endcode, like @c MpfxDataset
         @param dataset_base_dir dataset base directory
         @param dataset_dir_list dataset list of directories to search for input files
         @param dataset_recurse_dirs true for a recursive input file search, false otherwise
         @return the Dataset instance
      """

      return MpfxDataset(master, dataset_name, dataset_base_dir,
                                 dataset_dir_list, dataset_recurse_dirs)

   # -----------------------------------------------------------------------------------------------
   def get_dataset_module(self, master, dataset_name, dataset_type, dataset_base_dir):
      """!
         Load and return a specific Dataset class, ready to be instantiated to create a dataset.
         Such a class may be derived fronm class @c MpxDataset to support a specific survey Dataset.
         @param master master process instance
         @param dataset_name dataset name
         @param dataset_type prefix of a Dataset class, assumed to be of the form
                @code <dataset_type>Dataset @endcode, like @c MpfxDataset
         @param dataset_base_dir dataset base directory
         @return the Dataset object
      """

      if dataset_type == "Mpfx":

         dataset_module = sys.modules["mpfx.mpfx_data"]
         class_name = "{0}Dataset".format(dataset_type)
         dataset_ctor = eval("dataset_module." + class_name)
         if master.logging_enabled():
            master.logger.log_info_p(
              "{0} - Using Dataset type {1} [{2}]".format(
                     master.name, dataset_ctor, os.path.dirname(dataset_module.__file__)))
         return dataset_ctor

      else:

         try:
            prefix = "{0}.{0}_data".format(dataset_type.lower())
            class_name = "{0}_data.{1}Dataset".format(dataset_type.lower(), dataset_type)
            dataset_module = __import__(prefix, globals(), locals(), [], -1)
            dataset_ctor = eval("dataset_module." + class_name)

            if master.logging_enabled():
               master.logger.log_info_p(
                 "{0} - Loading Dataset type {1} [{2}]".format(
                        master.name, dataset_ctor, os.path.dirname(dataset_module.__file__)))
            return dataset_ctor

         except:
            if master.logging_enabled():
               master.logger.log_error_p(
                 "{0} - the Dataset {1}Dataset could not be loaded. "\
                  "Check that it has been installed "\
                  "and that it is referenced in PYTHONPATH.".format(master.name, dataset_type))
               master.logger.log_warning_p(
                 "{0} - Using {1} as default Dataset type "\
                  "instead of {2}Dataset ...".format(master.name, "MpfxDataset", dataset_type))
            return MpfxDataset

   # -----------------------------------------------------------------------------------------------
   def record_job_result(self, job_result, master):
      """!
         Add a job result to the list of job results. This method is invoked by the Master.
         @param master instance of the Master
         @param job_result object of class JobResult containing the data of the processed job
         @see Job, JobResult
      """

      img_no = job_result.job.img_no
      epoch = job_result.job.epoch

      result_dico = self.job_result_dico
      for branch_dir in job_result.job.branch_dirs:
         if not branch_dir in result_dico:
            result_dico[branch_dir] = {}
         result_dico = result_dico[branch_dir]

      if not img_no in result_dico:
         result_dico[img_no] = {}

      if not epoch is None and not epoch in result_dico[img_no]:
         result_dico[img_no][epoch] = {}

      # --- Record job result with the relevant branch tree structure
      if epoch is None:
         result_dico[img_no][None] = job_result
      else:
         result_dico[img_no][epoch] = job_result

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
      # --- Dump all results for checking
      if master.logging_enabled():
         for branch in self.job_branches:
            job_results = self.get_job_results(branch=branch)
            master.logger.log_info_p("{0} - Results from branch {1}: {2}".format(
                                     master, branch, [r.result for r in job_results]))
      else:
         pass

   # -----------------------------------------------------------------------------------------------
   def get_job_results(self, **kwargs):
      """!
         Return the list of MpfxJobResult objects for a specified branch, observation type,
         data_type and image number. These parameters must be specified using the
         "branch=", "img_no=" and "epoch=" keywords
         @return list of MpfxJobResult objects
         @code
            Examples
            - self.get_job_results(branch="control/space") to retrieve job results for all images
              and epochs in branch "control/space"
            - self.get_job_results(branch="control/space/variable", img_no=1) to retrieve
              job results of branch "control/space/variable", image number 1, all epochs
            - self.get_job_results(branch="control/space/variable", img_no=1, epoch=2) to retrieve
              job results of branch "control/space/variable", image number 1, epoch index 2
         @endcode
         @note MpfxJobResult objects can also be directly obtained from the
               mpfx.mpfx_job.MpfxJobProcessor.job_result_dico dictionary
      """
      job_result_list = []

      branch = kwargs.get("branch", None)    # branch tree path, e.g. "control/ground/constant"
      img_no = kwargs.get("img_no", None)
      epoch  = kwargs.get("epoch", None)

      if branch is not None:

         branch.strip("/")                # strip first "/" just in case...
         branch_dirs = branch.split("/")  # extract branch dir components

         # Point to the dictionary level of the leaf branch directory
         job_result_dico = self._job_result_dico
         for branch_dir in branch_dirs:
            job_result_dico = job_result_dico[branch_dir]

         if img_no is not None:
            if img_no in job_result_dico:
               if epoch is not None:
                  # image and epoch
                  if epoch in job_result_dico[img_no]:
                     job_result_list.append(job_result_dico[img_no][epoch])
               else:
                  # all epochs of the image
                  job_result_list = list(self._flatten_dict(job_result_dico[img_no]))
         else:
            if epoch is not None:
               self.helper.print_error(
                       "get_job_results(): an image number must be specified")
            else:
               # all image and epochs in the branch
               job_result_list = list(self._flatten_dict(job_result_dico))

      else:
         # Unspecified branch
         if img_no is None:
            job_result_list = list(self._flatten_dict(self._job_result_dico))
         else:
            self.helper.print_error(
                             "get_job_results(): a branch must be specified")

      return job_result_list


   # -----------------------------------------------------------------------------------------------
   def import_module(self, master, top_module, module_name):
      """! Find and load a module @c module_name from a directory @c module_dir """

      try:

         file_obj, filename, data = imp.find_module(top_module)

         return imp.load_module(module_name, file_obj, filename, data)

      except:
         if master.logging_enabled():
            master.logger.log_error_p(
              "{0} - The module {1} is not available: "\
              "check it has been installed and that it is referenced in PYTHONPATH)".format(
                                                                        master.name, module_name))


   # ~~~~~~~~~~~~~~~
   # Private methods
   # ~~~~~~~~~~~~~~~

   # -----------------------------------------------------------------------------------------------
   def _create_dir_tree(self, base_dir, branches):
      """!
         Create the directory tree for storing the results.
         @param base_dir base directory
         @param branches list of ordered directory names, making a path tree
         @note example of directory layout: control/ground/constant
         @see configuration file
      """

      for branch in branches:
         dir_names = branch.split("/")
         branch_path = base_dir
         for dir_name in dir_names:
            branch_path = os.path.join(branch_path, dir_name)
            self.helper.make_dir(branch_path)

   # -----------------------------------------------------------------------------------------------
   def _flatten_dict(self, d):
      """!
         Flatten a nested dictionary in order to return all values irrespective of nesting levels
      """
      for v in d.values():
         if isinstance(v, dict):
            for val in self._flatten_dict(v):
               yield val
         else:
            yield v


# -------------------------------------------------------------------------------------------------
class MpfxJob(Job):

   """!
       Represents a Mpfx job
   """

   def __init__(self, master, dataset, img_path_dico, *args):
      """!
         Construct a new job for identifiying a set of files in the dataset to process together.

         To obtain the path associated with a given type, just use:
         @code self.img_path_dico[type]@endcode with one of the aforementioned types.

         @param master master process instance
         @param dataset data source used to create the job
         @param img_path_dico file absolute path dictionary for every file types referenced in
                the job
         @param args a list with the branch tree, image number and epoch number

         @return an initialized MpfxJob object
      """

      img_no = args[-2]
      epoch_no = args[-1]
      branch_dirs = args[:-2]

      branch_key = list(branch_dirs)
      branch_key.extend([img_no, epoch_no])

      if epoch_no is None:
         # A file without epoch index, e.g. a catalog like star_catalog-000.txt
         job_name = "img_{0:03d}-*".format(img_no)
      else:
         job_name = "img_{0:03d}-{1:1d}".format(img_no, epoch_no)

      # Parent class initialization
      Job.__init__(self, master, job_name, dataset)

      # Set job attributes
      self._branch = "/".join(branch_dirs)   # branch (path) to the files of the job
      self._branch_dirs = branch_dirs        # list of directory names in the branch
      self._branch_key = tuple(branch_key)   # tuple (branch_dirs, im_no, epoch)
      self._img_no = img_no                  # image number associated with the job
      self._epoch = epoch_no                 # epoch index: only one value in this case
      self._img_path_dico = img_path_dico    # dictionary of file paths per file type
      self._file_names = self.get_file_names()

   # ~~~~~~~~~~~
   # Properties
   # ~~~~~~~~~~~

   # --- Getters

   @property
   def branch(self):
      """! @return the branch (path from base directory) containing the image file of the job. """
      return self._branch

   @property
   def img_no(self):
      """! @return the image number associated with the job. """
      return self._img_no

   @property
   def epoch(self):
      """! @return the epoch index associated with the job if any, None otherwise. """
      return self._epoch

   @property
   def epoch_nos(self):
      """! @return the epoch (multiple exposure) indice associated with the job
                   (only one in this class). """
      return [self.epoch]

   @property
   def branch_dirs(self):
      """! @return the list of directory names in the branch associated with the job. """
      return self._branch_dirs

   @property
   def branch_key(self):
      """! @return a tuple (branch_dirs, im_no, epoch). """
      return self._branch_key

   @property
   def img_path_dico(self):
      """! @return the file path dictionary associated with the job. """
      return self._img_path_dico

   @property
   def file_names(self):
      """! @return the list of file names associated with the job. """
      return self._file_names

   # --- Setters

   @branch.setter
   def branch(self, branch):
      """!
         Specify the branch (path from base directory) containing the image file of the job.
         @param branch the branch containing the image file of the job
      """
      self._branch = branch

   @img_no.setter
   def img_no(self, img_no):
      """!
         Specify the image number associated with the job.
         @param img_no the image number associated with the job
      """
      self._img_no = img_no

   @epoch.setter
   def epoch(self, epoch):
      """!
         Specify the epoch (exposure) index associated with the job.
         @param epoch the epoch (exposure) index associated with the job
      """
      self._epoch = epoch

   @branch_dirs.setter
   def branch_dirs(self, branch_dirs):
      """!
         Specify the list of directory names in the branch associated with the job.
         @param branch_dirs the list of directory names in the branch associated with the job
      """
      self._branch_dirs = branch_dirs

   @branch_key.setter
   def branch_key(self, branch_key):
      """!
         Specify a key made of a tuple (branch_dirs, im_no, epoch).
         @param branch_key a tuple (branch_dirs, im_no, epoch)
      """
      self._branch_key = branch_key

   @img_path_dico.setter
   def img_path_dico(self, img_path_dico):
      """!
         Specify the file path dictionary with the job.
         @param img_path_dico the file path dictionary with the job
      """
      self._img_path_dico = img_path_dico

   @file_names.setter
   def file_names(self, file_names):
      """!
         Specify the list of file names associated with the job.
         @param file_names list of file names associated with the job.
      """
      self._file_names = file_names

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
         @return the corresponding catalog name
         @see see get_file_types. __init__
      """
      return img_name.replace(".fits", cat_extension)

   # -----------------------------------------------------------------------------------------------
   def get_image_type_from_catalog(self, cat_name):
      """!
         Convenience method to find the image matching a given catalog file type,
         for instance
         "galaxy_catalog.txt" -> "image.fits" or "galaxy_catalog.fits" -> "image.fits"
         @param cat_name name of an catalog type
         @return the image name that corresponds to the catalog

         @see see get_file_types. __init__
      """

      return cat_name.replace(".txt", ".fits")

   # -----------------------------------------------------------------------------------------------
   def get_file_types(self):
      """!
         Return the list of object types managed by the job:
         @return the list of file types managed by the job
         @see list in comments for __init__
         @see get_file_types()
      """
      path_dico = self.img_path_dico.get(self.branch_key, None)
      if not path_dico is None:
         return path_dico.keys()
      else:
         return []

   # -----------------------------------------------------------------------------------------------
   def get_image_file_types(self, process):
      """!
         Return the list of image file types managed by the job:
         @param process master or worker process
         @return the list of image file types (stars or galaxies) managed by the job
         @see list in comments for __init__
         @see get_file_types(), get_image_file_types
      """

      file_types = []
      path_dico = self.img_path_dico.get(self.branch_key, None)
      if not path_dico is None:
         file_types = [p for p in path_dico.keys() if self.dataset.is_image(process, p)]
      return file_types

   # -----------------------------------------------------------------------------------------------
   def get_star_image_file_types(self, process):
      """!
         Return the list of star image file types managed by the job:
         @param process master or worker process
         @return the list of star image file types managed by the job
         @see list in comments for __init__
         @see get_file_types(), get_image_file_types(), get_galaxy_image_file_types()
      """

      file_types = []
      path_dico = self.img_path_dico.get(self.branch_key, None)
      if not path_dico is None:
         file_types = [p for p in path_dico.keys() if self.dataset.is_star_image(process, p)]
      return file_types

   # -----------------------------------------------------------------------------------------------
   def get_galaxy_image_file_types(self, process):
      """!
         Return the list of galaxy image file types managed by the job:
         @param process master or worker process
         @return the list of galaxy image file types managed by the job
         @see list in comments for __init__
         @see get_file_types(), get_image_file_types(), get_star_image_file_types()
      """

      file_types = []
      path_dico = self.img_path_dico.get(self.branch_key, None)
      if not path_dico is None:
         file_types = [p for p in path_dico.keys() if self.dataset.is_galaxy_image(process, p)]
      return file_types

   # -----------------------------------------------------------------------------------------------
   def get_catalog_file_types(self, process):
      """!
         Return the list of catalog file types managed by the job:
         @param process master or worker process
         @return the list of catalog file types (stars or galaxies) managed by the job
         @see list in comments for __init__
         @see get_file_types(), get_image_file_types()
      """

      file_types = []
      path_dico = self.img_path_dico.get(self.branch_key, None)
      if not path_dico is None:
         file_types = [p for p in path_dico.keys() if self.dataset.is_catalog(process, p)]
      return file_types

   # -----------------------------------------------------------------------------------------------
   def get_star_catalog_file_types(self, process):
      """!
         Return the list of star catalog file types managed by the job:
         @param process master or worker process
         @return the list of star catalog file types managed by the job
         @see list in comments for __init__
         @see get_file_types(), get_image_file_types(), get_galaxy_catalog_file_types()
      """

      file_types = []
      path_dico = self.img_path_dico.get(self.branch_key, None)
      if not path_dico is None:
         file_types = [p for p in path_dico.keys() if self.dataset.is_star_catalog(process, p)]
      return file_types

   # -----------------------------------------------------------------------------------------------
   def get_galaxy_catalog_file_types(self, process):
      """!
         Return the list of galaxy catalog file types managed by the job:
         @param process master or worker process
         @return the list of galaxy catalog file types managed by the job
         @see list in comments for __init__
         @see get_file_types(), get_image_file_types(), get_star_catalog_file_types()
      """

      file_types = []
      path_dico = self.img_path_dico.get(self.branch_key, None)
      if not path_dico is None:
         file_types = [p for p in path_dico.keys() if self.dataset.is_galaxy_catalog(process, p)]
      return file_types

   # -----------------------------------------------------------------------------------------------
   def get_file_names(self, file_type=None, epoch=None):
      """!
         return the list of file names managed by the job, possible for a given file type
         @param file_type optional type of the file (like "star_catalog.txt")
         @param epoch optional epoch index: if None, the epoch index is ignored
         @return the list of files managed by the job
      """

      file_names = []
      path_dico = self .img_path_dico.get(self.branch_key, None)
      if file_type is None:
         if epoch is None:
            file_names = [os.path.split(f)[1] for (t, f) in path_dico.items()]
         else:
            file_names = [os.path.split(f)[1] for (t, f) in path_dico.items() \
                                              if self.get_file_epoch_from_path(f) is not None and \
                                              self.get_file_epoch_from_path(f) == epoch]
      else:
         if epoch is None:
            file_names = [os.path.split(f)[1] for (t, f) in path_dico.items() if t == file_type]
         else:
            file_names = [os.path.split(f)[1] for (t, f) in path_dico.items() \
                                              if t == file_type and \
                                                self.get_file_epoch_from_path(f) is not None and \
                                                self.get_file_epoch_from_path(f) == epoch]

      return sorted(file_names)


   # -----------------------------------------------------------------------------------------------
   def get_file_paths(self, file_type=None, epoch=None):
      """!
         Return the list of object paths associated with the job
         @param file_type optional type of the file (like "star_catalog.txt")
         @param epoch optional epoch index: if None, the epoch index is ignored
         @return list of file paths associated with the job, possibly matching an epoch index
         @see get_file_types()
      """

      file_paths = []
      path_dico = self.img_path_dico.get(self.branch_key, None)
      if not path_dico is None:
         if file_type is None:
               if epoch is None:
                  file_paths = path_dico.values()
               else:
                  file_paths =  [ p for p in path_dico.values()
                                    if self.get_file_epoch_from_path(p) is not None and \
                                       self.get_file_epoch_from_path(p) == epoch ]
         else:
            if epoch is None:
               file_paths =  [ p for (t, p) in path_dico.items() if t == file_type ]
            else:
               file_paths =  [ p for (t, p) in path_dico.items() \
                                 if t == file_type and \
                                    self.get_file_epoch_from_path(p) is not None and \
                                    self.get_file_epoch_from_path(p) == epoch]
      return sorted(file_paths)

   # -----------------------------------------------------------------------------------------------
   def get_file_path(self, file_type=None, epoch=None):
      """!
         Return the path of a given file type. For instance, get_file_path("star_catalog.txt")
         @param file_type optional type type of the file (like "star_catalog.txt")
         @param epoch optional epoch index or None for any epoch
         @return path of the file identified by the file_type, None otherwise
         @see get_file_types()
      """

      path_list = [ p for p in self.get_file_paths(epoch=epoch)\
                            if self.get_file_type_from_path(p) == file_type]

      if len(path_list) > 0:
         return path_list[0]
      else:
         return None

   # -----------------------------------------------------------------------------------------------
   def get_catalog_file_paths(self, process):
      """!
         Return the list of catalog file paths managed by the job:
         @param process master or worker process
         @return the list of catalog file paths (stars or galaxies) managed by the job
         @see get_file_paths(), get_image_file_paths()
      """
      return [p for p in self.get_file_paths() if self.dataset.is_catalog(process, p)]

   # -----------------------------------------------------------------------------------------------
   def get_galaxy_catalog_file_paths(self, process):
      """!
         Return the list of galaxy catalog file paths managed by the job:
         @param process master or worker process
         @return the list of galaxy catalog file paths managed by the job
         @see get_catalog_file_paths(), get_star_catalog_file_paths()
      """
      return [p for p in self.get_file_paths() if self.dataset.is_galaxy_catalog(process, p)]

   # -----------------------------------------------------------------------------------------------
   def get_star_catalog_file_paths(self, process):
      """!
         Return the list of star catalog file paths managed by the job:
         @param process master or worker process
         @return the list of star catalog file paths managed by the job
         @see get_catalog_file_paths(), get_galaxy_catalog_file_paths()
      """
      return [p for p in self.get_file_paths() if self.dataset.is_star_catalog(process, p)]

   # -----------------------------------------------------------------------------------------------
   def get_image_file_paths(self, process):
      """!
         Return the list of image file paths managed by the job:
         @param process master or worker process
         @return the list of image file paths (stars or galaxies) managed by the job
         @see get_file_paths(), get_catalog_file_paths()
      """
      return [p for p in self.get_file_paths() if self.dataset.is_image(process, p)]

   # -----------------------------------------------------------------------------------------------
   def get_galaxy_image_file_paths(self, process):
      """!
         Return the list of galaxy image file paths managed by the job:
         @param process master or worker process
         @return the list of galaxy image file paths managed by the job
         @see get_image_file_paths(), get_star_image_file_paths()
      """
      return [p for p in self.get_file_paths() if self.dataset.is_galaxy_image(process, p)]

   # -----------------------------------------------------------------------------------------------
   def get_star_image_file_paths(self, process):
      """!
         Return the list of star image file paths managed by the job:
         @param process master or worker process
         @return the list of star image file paths managed by the job
         @see get_image_file_paths(), get_galaxy_image_file_paths()
      """
      return [p for p in self.get_file_paths() if self.dataset.is_star_image(process, p)]

   # -----------------------------------------------------------------------------------------------
   def get_file_components_from_path(self, file_path):
      """!
         Return the list of file components of a file absolute path: file type, image number
         and possibly epoch index.
         @code
         - If the file is of the form /<path>/<file_type>-<image_no>-<epoch>.<extension> a list
           [file_type, img_no, epoch] is returned
         - If the file is of the form /<path>/<file_type>-<image_no>.<extension>, a list
           [file_type, img_no, None] is returned
         @endcode
         @param file_path absolute file path of a file (image, catalog, etc.)
         @code
           Examples:
           - get_file_compoments_from_path("starfield-000.1.fits") => [starfield_image.fits, 0, 1]
           - get_file_compoments_from_path("star_catalog-000.txt") => [star_catalog.txt, 0, None]
         @endcode
         @return
         @see get_file_type_from_path(), get_file_epoch_from_path()
      """

      _, filename = os.path.split(file_path)
      _, fileext = os.path.splitext(filename)

      regex_dash = re.compile("(.+)\-(.*[0-9])\-(.*[0-9])\.")
      result = regex_dash.search(filename)
      if not result is None:
         return [result.group(1)+fileext, int(result.group(2)), int(result.group(3))]
      else:
         regex_dot = re.compile("(.+)\-(.*[0-9])\.")
         result = regex_dot.search(filename)
         return [result.group(1)+fileext, int(result.group(2)), None]


   # -----------------------------------------------------------------------------------------------
   def get_file_type_from_path(self, file_path):
      """!
         Return the file type of a file with a given absolute file path
         @param file_path absolute file path of a file (image, catalog, etc.)
         @return the file type of a file (like star_catalog.fits)
         @note: None is returned if no file type could be found
         @see get_file_components_from_path()
      """
      components = self.get_file_components_from_path(file_path)

      if len(components) > 0:
         return components[0]
      else:
         return None

   # -----------------------------------------------------------------------------------------------
   def get_file_image_no_from_path(self, file_path):
      """!
         Return the image number of a file with a given absolute file path
         @param file_path absolute file path of a file (image, catalog, etc.)
         @return the file image number
         @note: None is returned if no file type could be found
         @see get_file_components_from_path()
      """
      components = self.get_file_components_from_path(file_path)
      if len(components) > 1:
         return components[1]
      else:
         return None

   # -----------------------------------------------------------------------------------------------
   def get_file_epoch_from_path(self, file_path):
      """!
         Return the epoch index of a file with a given absolute file path
         @param file_path absolute file path of a file (image, catalog, etc.)
         @return the file epoch index
         @note: None is returned if no file type could be found
         @see get_file_components_from_path()
      """
      components = self.get_file_components_from_path(file_path)
      if len(components) > 2:
         return components[2]
      else:
         return None

   # -----------------------------------------------------------------------------------------------
   def get_branch_tree(self):
      """!
         Convenience method to return the jobs's branch tree structure, for example
         "control/space/constant".
         @return branch directory tree for this job
      """
      return self.branch



# -------------------------------------------------------------------------------------------------
class MpfxJobResult(JobResult):

   """!
       Represents the result of a job.
   """

   def __init__(self, worker, job, result):
      """! Construct a MpfxJobResult job result object """
      JobResult.__init__(self, worker, job, result)




# --------------------------------------------------------------------------------------------------
class MpfxMultiEpochJobProcessor(MpfxJobProcessor):

   """!
      Extension of MpfxJobProcessor: submit jobs and process associated job results, taking care
      of multiepoch (multiple exposure) images and related catalogs.
      Based on mpf.MpfxJobProcessor.
   """

   def __init__(self, master):
      """ MpfxMultiEpochJobProcessorProcessor default constuctor """

      MpfxJobProcessor.__init__(self, master)


   # ~~~~~~~~~~~~~~~
   # Public methods
   # ~~~~~~~~~~~~~~~

   # -----------------------------------------------------------------------------------------------
   def create_dataset(self, master, dataset_name, dataset_type, dataset_base_dir,
                                    dataset_dir_list, dataset_recurse_dir):
      """!
         Create a Dataset object that represent the data source for images and catalogs.
         @param master master process instance
         @param dataset_name dataset name
         @param dataset_type prefix of a Dataset class, assumed to be of the form
                @code <dataset_type>Dataset @endcode, like @c MpfxDataset
         @param dataset_base_dir dataset base directory
         @param dataset_dir_list dataset list of directories to search for input files
         @param dataset_recurse_dir true for a recursive input file search, false otherwise
         @return the Dataset instance, None in case of error
      """

      return MpfxMultiEpochDataset(master, dataset_name, dataset_base_dir,
                                           dataset_dir_list, dataset_recurse_dir)

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
         @param master master process instance
         @param dataset data source used to create the job
         @param args a list of arguments chosen by the caller and whose nature depends on the
                data necessary to process the job
         @return a new MpfxJob job object
         @see Job, MpfxJob
      """

      # --- Create a MpfxMultiEpochJob instance, based on parent class mpfx.MpfxJob
      return MpfxMultiEpochJob(master, dataset, *args)

   # -----------------------------------------------------------------------------------------------
   def record_job_result(self, job_result, master):
      """!
         Add a job result to the list of job results. This method is invoked by the Master.
         @param master instance of the Master
         @param job_result object of class JobResult containing the data of the processed job
         @see Job, JobResult
         @note individual epochs are not recorded, unlike for mono-epoch jobs
      """

      img_no = job_result.job.img_no

      result_dico = self.job_result_dico
      for branch_dir in job_result.job.branch_dirs:
         if not branch_dir in result_dico:
            result_dico[branch_dir] = {}
         result_dico = result_dico[branch_dir]

      if not img_no in result_dico:
         result_dico[img_no] = {}

      # --- Record job result with the relevant branch tree structure
      result_dico[img_no] = job_result

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

      data_dico = {"result": "result from MpfxMultiEpochJob: {0}".format(job.name)}    # sample result data

      return MpfxMultiEpochJobResult(worker, job, data_dico)  # store the job results

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

#       # --- Dump all results for checking
#       if master.logging_enabled():
#          for branch in self.job_branches:
#             job_results = self.get_job_results(branch=branch, img_no=0)
#             master.logger.log_info_p("{0} - Results from branch {1}: {2}".format(
#                                      master, branch, [r.result for r in job_results]))
#       else:
#          pass
      pass

   # -----------------------------------------------------------------------------------------------
   def get_job_results(self, **kwargs):
      """!
         Return the list of MpfxMultiEpochJobResult objects for a specified branch, observation type,
         data_type and image number. These parameters must be specified using the
         "branch=", "img_no=" keywords
         @return list of MpfxMultiEpochJobResult objects
         @code
            Examples
            - self.get_job_results(branch="control/space") to retrieve job results for all images
              and epochs in branch "control/space"
            - self.get_job_results(branch="control/space/variable", img_no=1) to retrieve
              job results of branch "control/space/variable", image number 1, all epochs
         @endcode
         @note MpfxMultiEpochJobResult objects can also be directly obtained from the
               mpfx.mpfx_job.MpfxMultiEpochJobProcessor.job_result_dico dictionary
      """

      job_result_list = []

      branch = kwargs.get("branch", None)    # branch tree path, e.g. "control/ground/constant"
      img_no = kwargs.get("img_no", None)

      if branch is not None:

         branch.strip("/")                # strip first "/" just in case...
         branch_dirs = branch.split("/")  # extract branch dir components

         # Point to the dictionary level of the leaf branch directory
         job_result_dico = self._job_result_dico
         for branch_dir in branch_dirs:
            job_result_dico = job_result_dico[branch_dir]

         if img_no is not None:
            if img_no in job_result_dico:
               # all epochs of the image
               job_result_list.append(job_result_dico[img_no])
         else:
            # all image and epochs in the branch
            job_result_list = list(self._flatten_dict(job_result_dico))

      else:
         # Unspecified branch
         if img_no is None:
            job_result_list = list(self._flatten_dict(self._job_result_dico))
         else:
            self.helper.print_error(
                             "get_job_results(): a branch must be specified")

      return job_result_list



# -------------------------------------------------------------------------------------------------
class MpfxMultiEpochJob(MpfxJob):

   """!
       Represents a MpfxMultiEpochJob job, an extension of a MpfxJob for holding a group of
       images with multiple exposures, with related catalogs.
   """

   def __init__(self, master, dataset, img_path_dico, *args):
      """!
         Construct a new job for identifying a set of files in the dataset to process together.

         To obtain the path associated with a given type, just use:
         @code self.img_path_dico[type]@endcode with one of the aforementioned types.

         @param master master process instance
         @param dataset data source used to create the job
         @param img_path_dico file absolute path dictionary  for every file types referenced in
                the job
         @param args a list with the branch tree, image number and epoch number

         @return an initialized MpfxMultiEpochJob object
      """

      # ~~~~~~~~~~~
      # Properties
      # ~~~~~~~~~~~

      img_no = args[-1]
      branch_dirs = args[:-1]

      # Set job attributes
      self.branch = "/".join(branch_dirs)    # branch (path) to the files of the job
      self.branch_dirs = branch_dirs         # list of directory names in the branch
      self.img_no = img_no                   # image number associated with the job
      self.img_path_dico = img_path_dico     # dictionary of file paths per file type

      branch_key = list(branch_dirs)
      branch_key.extend([img_no])
      self.branch_key = tuple(branch_key)
      self._file_names = []

      # --- Build epoch list
      self._epoch_nos = []

      reg_epoch = re.compile("\-(.*[0-9]\-(.*[0-9])\.")

      path_dico = img_path_dico.get(self.branch_key, None)
      if path_dico is not None:
         file_lists = [l for l in path_dico.values() if dataset.is_image(master, l[0])]
         if len(file_lists) > 0:
            self._epoch_nos = sorted([ int(reg_epoch.search(os.path.split(t)[1]).group(2)) \
                                       for t in  file_lists[0] ])
         else:
            if master.logging_enabled():
               master.helper.print_error("/{0}/image-{1:03d} - "\
                                         "Job creation failure: empty image file list".format(
                                                                             self.branch, img_no))

         # --- Set job name
         if len(dataset.dir_list) > 0:
            self.job_name = "{0}_img_{1:03d}-[{2}]".format("_".join(branch_dirs), img_no,
                                                     "-".join([str(e) for e in self._epoch_nos]) )
         else:
            self.job_name = "img_{0:03d}-[{1}]".format(img_no,
                                                     "-".join([str(e) for e in self._epoch_nos]) )

         # Base class initialization (Note: this is NOT the parent class)
         Job.__init__(self, master, self.job_name, dataset)
         self._file_names = self.get_file_names()

      else:
         # Job creation failure
         if master.logging_enabled():
            master.helper.print_error("/{0}/image-{1:03d} - "\
                                      "Job creation failure-: unknown directory path {2}".format(
                                                            self.branch, img_no, self.branch_key))


   # ~~~~~~~~~~~
   # Properties
   # ~~~~~~~~~~~

   # --- Getters

   @property
   def epoch_nos(self):
      """! @return the epoch (multiple exposure) indice associated with the job. """
      return self._epoch_nos

   @property
   def epoch(self):
      """! @return the first epoch index associated with the job. """
      return self._epoch_nos[0]

   @property
   def file_names(self):
      """! @return the list of file names associated with the job. """
      return self._file_names

   # --- Setters

   @epoch_nos.setter
   def epoch_nos(self, epoch_nos):
      """!
         Specify the epoch (multiple exposure) indice associated with the job.
         @param epoch_nos the epoch indice associated with the job.
      """
      self._epoch_nos = epoch_nos

   # -----------------------------------------------------------------------------------------------
   def get_file_names(self, file_type=None, epoch=None):
      """!
         return the list of file names managed by the job, possible for a given file type
         @param file_type optional type of the file (like "star_catalog.txt")
         @param epoch optional epoch index: if None, the epoch index is ignored
         @return the list of files managed by the job
      """

      path_dico = self.img_path_dico.get(self.branch_key, None)
      if file_type is None:
         file_list = [l for (t, l) in path_dico.items()]
      else:
         file_list = [l for (t, l) in path_dico.items() if t == file_type]

      file_paths = []
      [file_paths.extend(l) for l in file_list]
      if epoch is None:
         return sorted([os.path.split(f)[1] for f in file_paths])
      else:
         return sorted([os.path.split(f)[1] for f in file_paths \
                                            if self.get_file_epoch_from_path(f) is not None and \
                                               self.get_file_epoch_from_path(f) == epoch])

   # -----------------------------------------------------------------------------------------------
   def get_file_paths(self, file_type=None, epoch=None):
      """!
         Return the list of object paths (like associated with the job
         @param file_type optional type of the file (like "star_catalog.txt")
         @param epoch optional epoch index: if None, the epoch index is ignored
         @return list of file paths associated with the job, possibly matching an epoch index
         @see get_file_types()
      """

      path_dico = self.img_path_dico.get(self.branch_key, None)
      if file_type is None:
         file_list = [l for (t, l) in path_dico.items()]
      else:
         file_list = [l for (t, l) in path_dico.items() if t == file_type]

      file_paths = []
      [file_paths.extend(l) for l in file_list]
      if epoch is None:
         return sorted([f for f in file_paths])
      else:
         return sorted([f for f in file_paths \
                                            if self.get_file_epoch_from_path(f) is not None and \
                                               self.get_file_epoch_from_path(f) == epoch])


# -------------------------------------------------------------------------------------------------
class MpfxMultiEpochJobResult(JobResult):

   """!
       Represents the result of a job with multiple epochs (exposures) images and related catalogs.
   """

   def __init__(self, worker, job, result):
      """! Construct a MpfxJobResult job result object """
      JobResult.__init__(self, worker, job, result)


# -- EOF mpfx_job.py
