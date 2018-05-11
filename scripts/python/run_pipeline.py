#!/usr/bin/env python

"""Script run_pipeline.py

Run ShapePipe pipeline

:Authors: Martin Kilbinger

:Date: 17/04/2018
"""


# Compability with python2.x for x>6
from __future__ import print_function


import sys
import os
import re
import copy
import glob

import numpy as np

from optparse import OptionParser, IndentedHelpFormatter, OptionGroup

import stuff


"""
Call each module
	- check whether required executables are present
	- check whether input files are present
	  if not prepare module
	- check config files, pointing to right dirs? 
	- run
	- ask to adopt/not adopt/postpone run

Adpot last run
	- create link output/<module>/adopt

Prepare module
	- in data create links to previous module's adopted run results

Requires on input (config)
	- List of modules (for each pipeline sequence mode)

"""

# Global variable, solve differenly
modules_glob = {}
modules_glob['std'] = ['mask', 'SExtractor', 'SETools', 'PSFExRun', 'PSFExInterpolation']

path_sp     = '{}/ShapePipe'.format(os.environ['HOME'])
path_spmod  = '{}/modules'.format(path_sp)
path_output = 'output' 
path_data   = '{}/data'.format(os.environ['HOME'])

name_adopted = 'adopted'
name_results = 'results'


def params_default():
    """Set default parameter values.

    Parameters
    ----------
    None

    Returns
    -------
    p_def: class tuff.param
        parameter values
    """

    p_def = stuff.param(
        scheme  = 'std',
        band  = 'r',
        image_type  = 'tile',
    )

    return p_def



def parse_options(p_def):
    """Parse command line options.

    Parameters
    ----------
    p_def: class tuff.param
        parameter values

    Returns
    -------
    options: tuple
        Command line options
    args: string
        Command line string
    """

    usage  = "%prog [OPTIONS]"
    parser = OptionParser(usage=usage, formatter=stuff.IndentedHelpFormatterWithNL())

    #parser.add_option('-s', '--scheme', dest='scheme', type='string', default=p_def.scheme,
     #    help='Scheme, default={}'.format(p_def.scheme))
    parser.add_option('-M', '--module', dest='module', type='string', default=None,
         help='Pipeline module, see \'-l\' for list of all modules')

    # Run mode
    parser.add_option('-m', '--mode', dest='mode', default=None,
         help='run mode, one of [l|r|a]:\n'
          ' l: list modules of all scheme and exit\n'
          ' r: run module given by \'-M\'\n'
          ' a: adopt (last) run for module given by \'-M\'\n'
          ' s: set input file links in <data> to results from adopted run of module given by \'-M\'\n')

    # Monitoring
    parser.add_option('-v', '--verbose', dest='verbose', action='store_true', help='verbose output')

    options, args = parser.parse_args()

    return options, args



def check_options(options):
    """Check command line options.

    Parameters
    ----------
    options: tuple
        Command line options

    Returns
    -------
    erg: bool
        Result of option check. False if invalid option value.
    """

    #if not options.scheme in modules_glob:
        #print('Invalid scheme {}, not found in '.format(options.scheme), modules_glob.keys())
        #return False

    if options.mode is None:
        print('No run mode given (option \'-r\')')

    if len(options.mode) > 1:
        print('Invalid run mode \'{}\', only one allowed'.format(options.mode))
        return False

    if options.mode == 'r' or options.mode == 'a':
        if not options.module:
            print('Module needs to be given, via \'-M\'')
            return False

    return True



def update_param(p_def, options):
    """Return default parameter, updated and complemented according to options.
    
    Parameters
    ----------
    p_def:  class stuff.param
        parameter values
    optiosn: tuple
        command line options
    
    Returns
    -------
    param: class stuff.param
        updated paramter values
    """

    param = copy.copy(p_def)

    # Update keys in param according to options values
    for key in vars(param):
        if key in vars(options):
            setattr(param, key, getattr(options, key))

    # Add remaining keys from options to param
    for key in vars(options):
        if not key in vars(param):
            setattr(param, key, getattr(options, key))

    return param



def list_modules():
    """List modules of all schemes in order.

    Parameters
    ----------
    None

    Returns
    -------
    None
    """

    for s in modules_glob:
        print('Scheme {}'.format(s))
        for i, m in enumerate(modules_glob[s]):
            print(i, m)
        print()



def run_module(module, verbose=False):
    """Run module.
    
    Parameters
    ----------
    module: sring
        module name
    verbose: bool, optional, default=False

    Returns
    -------
    None
    """

    package     = '{}_package'.format(module)
    launch_path = '{}/{}/config/launch.cmd'.format(path_spmod, package)

    stuff.run_cmd(launch_path, run=True, verbose=verbose, devnull=False)
    


def adopt_run(module, verbose=False):
    """Adopt (last) run for given module.

    Parameters
    ----------
    module: string
        module name
    verbose: bool, optional, default=False
        verbose output if True

    Returns
    -------
    None
    """

    path_runs = '{}/{}'.format(path_output, module)
    f = glob.glob('{}/run_*'.format(path_runs))
    f = sorted(f)
    last_run  = f[-1]
    source          = os.path.basename(last_run)
    source_to_check = last_run
    link_name = '{}/{}'.format(path_runs, name_adopted)

    stuff.ln_s(source, link_name, orig_to_check=source_to_check, verbose=verbose, force=True)



def set_results(module, verbose=False):
    """Set file links in <data> to adopted run for given module.

    Parameters
    ----------
    module: string
        module name
    verbose: bool, optional, default=False
        verbose output if True

    Returns
    -------
    None
    """

    path_results = '{}/{}/{}/{}/{}'.format(os.getcwd(), path_output, module, name_adopted, name_results)

    files = glob.glob('{}/*'.format(path_results))
    for f in files:
        source    = f
        link_name = '{}/{}'.format(path_data, os.path.basename(f))
        stuff.ln_s(source, link_name, orig_to_check=source, verbose=verbose, force=True)



def main(argv=None):
    """Main program.
    """

    # Set default parameters
    p_def = params_default()

    # Command line options
    options, args = parse_options(p_def)
    # Without option parsing, this would be: args = argv[1:]

    if check_options(options) is False:
        return 1

    param = update_param(p_def, options)


    # Save calling command
    stuff.log_command(argv)
    if param.verbose:
        stuff.log_command(argv, name='sys.stderr')

    if param.verbose is True:
        print('Start of program {}'.format(os.path.basename(argv[0])))


    ### Start main program ###


    if param.mode == 'l':
        list_modules()

    elif param.mode == 'r':
        run_module(param.module, verbose=param.verbose)

    elif param.mode == 'a':
        adopt_run(param.module, verbose=param.verbose)

    elif param.mode == 's':
        set_results(param.module, verbose=param.verbose)

    ### End main program

    if param.verbose is True:
        print('End of program {}'.format(os.path.basename(argv[0])))


    return 0



if __name__ == "__main__":
    sys.exit(main(sys.argv))


