#!/usr/bin/env python

"""!
   @package mpfg.mpSMP SMP main module
   @author Marc Gentile
   @file mp_MP.py

   mp_smp.py: SMP main module

   Template for multiprocessing using SMP (Symmetric Multiprocessing). 
   This code implements a simple producer-consumer scheme for processing a 
   large number of images (or objects in general) uisng a SMP multiprocessing architecture.
   A producer process (the "master", or "Manager") is responsible for providing the images to
   a pool of consumer processes ("salves" or "workers") and for collecting their output. 
   A worker process is responsible for processing the input data provided by
   the master as well as returning the corresponding results.
"""

# -- Python imports
import sys

# -- Module-specific imports
import mp_args as arg            # command-line arguments handling
import mp_calc_SMP as smp        # calculation module (SMP)
import mp_helper                 # helper functions

# ~~~~~~~~~~~~~~~~~
# Main entry point
# ~~~~~~~~~~~~~~~~~

if __name__ == '__main__':
   """! Main entry point, either for the master or for a particular worker """
 
   # --- Command line arguments and options
   args = arg.Args() 

   if args.options["-h"] or args.options["--help"]:
      args.print_usage(args.helper.get_main_basename())
      sys.exit(0)

   # --- Invoke the master calculator in module: mp_calc_SMP.py
   master = smp.MasterSMP(args)
   if master is None:
      sys.exit(2)
   else:
      try:
         master.run()
      except Exception as detail:
         mp_helper.Helper.print_error("The master process ended unexpectedly: {0}", detail)
      finally:
         master.shutdown()
         sys.exit(0)


# -- EOF mp_SMP.py
