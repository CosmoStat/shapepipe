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
modules_glob['std'] = ['select', 'mask', 'SExtractor', 'SETools', 'PSFExRun', 'PSFExInterpolation']

path_sp     = '{}/ShapePipe'.format(os.environ['HOME'])
path_spmod  = '{}/modules'.format(path_sp)
path_sppy   = '{}/scripts/python'.format(path_sp)
path_output = 'output' 
path_data   = '{}/data'.format(os.environ['HOME'])

path_CFIS_data = '{}/astro/data/CFIS'.format(os.environ['HOME'])

name_adopted = 'adopted'
name_results = 'results'


class modules_local:
    """Has modules as subroutines, priority over pipeline packages with same name.
    """

    def init(self):
        pass

    def select(self, param):
        """Area selection.
        Parameters
        ----------
        param: class param
            parameter values

        Returns
        -------
        None
        """

        # Check whether data dir does not exist
        if os.path.isfile(path_data):
            stuff.error('Path {} exists, please remove for field selection step'.format(path_data))
        os.mkdir(path_data)

        # Find fields in given area
        fixed_options = '-i {}/tiles -t tile -m a --plot -v'.format(path_CFIS_data)
        launch_path = '{}/cfis_field_select.py {} {}'.format(path_sppy, fixed_options, param.options)
        stuff.run_cmd(launch_path, run=not param.dry_run, verbose=param.verbose, devnull=False) 

        # Output file base name
        m = re.search('-o\s+(\S+)', param.options)
        if m:
            name = m.groups()[0]
        else:
            stuff.error('No output file basename in options \'{}\' found'.format(param.options))

        # Create links to original files
        launch_path = '{}/create_image_links.py -i {}.txt -v -o {}'.format(path_sppy, name, path_data)
        stuff.run_cmd(launch_path, run=not param.dry_run, verbose=param.verbose, devnull=False)


def params_default():
    """Set default parameter values.

    Parameters
    ----------
    None

    Returns
    -------
    p_def: class param
        parameter values
    """

    p_def = stuff.param(
        scheme  = 'std',
        band  = 'r',
        image_type  = 'tile',
        job   = 'manual',
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

    # Pipeline and modules
    parser.add_option('-s', '--scheme', dest='scheme', type='string', default=p_def.scheme,
            help='scheme, default=\'{}\''.format(p_def.scheme))
    parser.add_option('-M', '--module', dest='module', type='string', default=None,
            help='pipeline module, see \'-l\' for list of all modules')
    parser.add_option('-O', '--options', dest='options', type='string', default=None,
            help='options for (local) modules')

    # Job execution
    parser.add_option('-j', '--job', dest='job', type='string', default=p_def.job,
            help='job exceution, one of [\'manual\'|\'qsub\'], default=\'{}\''.format(p_def.job))
    parser.add_option('-n', '--dry-run', dest='dry_run', action='store_true', default=False,
            help='dry run, only print commands')
 
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

    if not options.scheme in modules_glob:
        print('Invalid scheme {}, not found in '.format(options.scheme), modules_glob.keys())
        return False

    if options.mode is None:
        print('No run mode given (option \'-r\')')
        return False

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
    param: class param
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

    if param.module is not None:
        try:
            param.module = int(param.module)
            if param.module < 0:
                raise IndexError
            param.smodule = modules_glob[param.scheme][param.module]

        except ValueError:
            if not param.module in modules_glob[param.scheme]:
                stuff.error('Module \'{}\' not found, list available options with \'-m l\''.format(param.module))
            param.smodule = param.module

        except IndexError:
            stuff.error('Invalid module number {}, list available options with \'-m l\''.format(param.module))

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



def create_qsub_script(module):
    """Create bash script to be submitted with qsub.

    Parameters
    ----------
    module: string
        module name

    Returns:
    qsub_script_name: string
        job script name
    """


    qsub_script_name = 'job_{}.sh'.format(module)

    f = open(qsub_script_name, 'w')
    print('#!/usr/bin/bash\n', file=f)
    print('#PBS -S /usr/bin/bash\n', file=f)
    print('#PBS -N {}\n'.format(module), file=f)
    print('#PBS -o job.out', file=f)
    print('#PBS -j oe\n', file=f)
    print('#PBS -l nodes=1:ppn=1,walltime=10:00:00\n', file=f)

    cwd = os.getcwd()
    print('#PBS -d {}\n'.format(cwd), file=f)

    print('module load intelpython/2', file=f)
    print('export PATH="$PATH:/softs/astromatic/bin/:$HOME/.local/bin"', file=f)
    print('export LD_LIBRARY_PATH="$HOME/.local/lib/"', file=f)
    print('export PYTHONPATH="$HOME/.local/lib/python2.7/site-packages:$HOME/.local/bin:/opt/intel/intelpython2/lib/python2.7/site-packages"\n', file=f)

    print('echo -n "pwd = "', file=f)
    print('pwd\n', file=f)

    print('$HOME/ShapePipe/modules/{}_package/config/launch.cmd'.format(module), file=f)
    print('ex=$?', file=f)
    print('exit $ex\n', file=f)

    f.close()
    return qsub_script_name


def run_module(param):
    """Run module param.module as subroutine from modules_local,
       if exists, or it not, as pipeline package.
    
    Parameters
    ----------
    param: class param
        paramter values

    Returns
    -------
    None
    """

    # Module name as string
    module = param.smodule

    # Class containing modules as methods
    ml     = modules_local()

    if hasattr(ml, module):
        method = getattr(ml, module, None)
        if callable(method):
            # Run module as method of local class
            if param.verbose:
                print('Running local module \'{}\''.format(module))
            getattr(ml, module)(param)

    else:
        # Run pipeline module
        if param.verbose:
            print('Running pipeline module \'{}\''.format(module))

        # Create output dir if necessary
        stuff.mkdir_p(path_output, verbose=param.verbose)

        package     = '{}_package'.format(module)
        launch_path = '{}/{}/config/launch.cmd'.format(path_spmod, package)

        # Run module
        if param.job == 'manual':
            ex = stuff.run_cmd(launch_path, run=not param.dry_run, verbose=param.verbose, devnull=False)
            if param.verbose:
                print('Sum[exit codes] = {}'.format(ex))

        elif param.job == 'qsub':

            qsub_script_name = create_qsub_script(module)
            stuff.run_cmd('qsub {}'.format(qsub_script_name), run=not param.dry_run, verbose=param.verbose, devnull=False)

        else:
            stuff.error('Invalid job execution mode \'{}\''.format(param.job))

    

def adopt_run(module, verbose=False):
    """Adopt (last) run for given moduleo (set symbolic link).

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


def discard_run(module, verbose=False):
    """Discard adopted run for given module (remove symbolic link)

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

    path = '{}/{}/adopted'.format(path_output, module)
    if os.path.islink(path):
        if verbose:
            print('Discarding adopted run of module \'{}\''.format(module))
        os.remove(path)
    else:
        stuff.error('No adopted run found for module \'{}\''.format(module))



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
        run_module(param)

    elif param.mode == 'a':
        adopt_run(param.smodule, verbose=param.verbose)

    elif param.mode == 'd':
        discard_run(param.smodule, verbose=param.verbose)

    elif param.mode == 's':
        set_results(param.smodule, verbose=param.verbose)

    ### End main program

    if param.verbose is True:
        print('End of program {}'.format(os.path.basename(argv[0])))


    return 0



if __name__ == "__main__":
    sys.exit(main(sys.argv))


