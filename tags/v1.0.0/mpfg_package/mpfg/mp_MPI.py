#!/usr/bin/env python 

"""! 
   @package mpfg.mp_MPI MPI main module
   @author Marc Gentile
   @file mp_MPI.py

   mp_MPI.py: MPI main module

   Template for multiprocessing using MPI (Message Passing Interface). 
   This code implements a simple producer-consumer scheme for processing a 
   large number of images (or objects in general) using the MPI multiprocessing architecture.
   A producer process (the "master", or "Manager") is responsible for providing the images to
   a pool of consumer processes ("salves" or "workers") and for collecting their output. 
   A worker process is responsible for processing the input data provided by
   the master as well as returning the corresponding results.
"""

# -- Python imports
import os, sys
from mpi4py import MPI  # MPI interface

# -- Module-specific imports
import mp_args as arg         # command-line arguments and options
import mp_calc_MPI as mpi     # calculation module (MPI)

# ~~~~~~~~~~~~~~~~~
# Main entry point
# ~~~~~~~~~~~~~~~~~


if __name__ == '__main__':
   """! Main entry point, either for the master or for a particular worker """

   comm = MPI.COMM_WORLD
   rank = comm.Get_rank()  # rank: 0 for the main process, > 0 for workers

   args = arg.Args()    # command-line arguments and options

   if rank == 0:

      # --- Check if Help requested
      if args.options["-h"] or args.options["--help"]:
         args.print_usage(args.helper.get_main_basename())
         sys.exit(0)

      # --- Master process
      master = mpi.MasterMPI(args, comm)
      if master is None:
         sys.exit(2)
      else:
         try:
            master.run()
         except Exception as detail:
            print_error("The master process ended unexpectedly: {0}", detail)
         finally:
            master.shutdown()
            sys.exit(0)

   else:

      # --- Check if Help requested
      if args.options["-h"] or args.options["--help"]:
         sys.exit(0)

      # --- Worker process
      worker = mpi.WorkerMPI(args, comm, rank)  
      if worker is None:
         sys.exit(2)
      else:
         try:
            worker.run()
         except Exception as detail:
            print_error("The worker process {0} ended unexpectedly {1}".format(
                                                                             worker.name, detail))
         finally:
            sys.exit(0)

# -- EOF mp_MPI.py


