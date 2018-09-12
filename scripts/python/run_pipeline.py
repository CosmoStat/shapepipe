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
import shutil

import numpy as np

from optparse import OptionParser, IndentedHelpFormatter, OptionGroup

import stuff


# Global variables, TODO: do differenly
modules_glob = {}
types_glob = {}
modules_glob['std'] = ['select', 'mask', 'SExtractor', 'SETools', 'PSFExRun', 'PSFExInterpolation']
types_glob['std'] = ['tile', 'tile', 'tile', 'tile', 'tile', 'tile']

modules_glob['tiles_exp'] = ['select', 'mask', 'SExtractor', 'find_exp', 'mask']
types_glob['tiles_exp'] = ['tile', 'tile', 'tile', 'exposure', 'exposure']

# Basic paths for pipeline codes
path_sp     = '{}/ShapePipe'.format(os.environ['HOME'])
path_spmod  = '{}/modules'.format(path_sp)
path_sppy   = '{}/scripts/python'.format(path_sp)
config_name = 'package_config_smp.cfg'

# Run-time create paths
path_config = {}
path_config['tile']     = 'input_tile'
path_config['exposure'] = 'input_exp'
path_output = {}
path_output['tile']     = 'output_tile' 
path_output['exposure'] = 'output_exp'

# Data path(s)
path_data = {}
path_data['tile']     = '{}/data'.format(os.environ['HOME'])
path_data['exposure'] = '{}/hdu'.format(path_data['tile'])

path_CFIS_data = '{}/astro/data/CFIS'.format(os.environ['HOME'])
data_file_base = {}
data_file_base['tile'] = 'CFIS'
data_file_base['exposure'] = 'cfisexp'.format(data_file_base['tile'])

name_adopted = 'adopted'
name_results = 'results'



class modules_local:
    """Has modules as subroutines, priority over pipeline packages with same name.
    """

    def init(self):
        pass

    def select(self, param):
        """Area selection. Create links to tiles.

        Parameters
        ----------
        param: class param
            parameter values

        Returns
        -------
        None
        """

        # Create tile links directory
        if os.path.isfile(path_data['tile']):
            stuff.error('Path {} exists, please remove for field selection step (creating links)'.format(path_data['tile']))
        stuff.mkdir_p(path_data['tile'])

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
        launch_path = '{}/create_image_links.py -i {}.txt -v -o {}'.format(path_sppy, name, path_data['tile'])
        stuff.run_cmd(launch_path, run=not param.dry_run, verbose=param.verbose, devnull=False)


    def find_exp(self, param):
        """Find exposures used in selected tiles and create links.

        Parameters
        ----------
        param: class param
            parameter values

        Returns
        -------
        None
        """

        # Create exposure links directory
        if os.path.isdir(path_data['exposure']):
            stuff.error('Path {} exists, please remove before setting links to exposures'.format(path_data['exposure']))
        os.mkdir(path_data['exposure'])

        if param.verbose:
            verbose_flag = ' -v'
        else:
            verbose_flag = ''
        cmd = '{}/cfis_create_exposures.py -i {} -o {} -p \'CFIS-\' --exp_base_new=\'{}\' -O hdu{}'.\
            format(path_sppy, path_data['tile'], path_data['exposure'], data_file_base['exposure'], verbose_flag)
        stuff.run_cmd(cmd, run=not param.dry_run, verbose=param.verbose, devnull=False) 


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
            help='pipeline module, see \'-m l\' for list of all modules')
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
            param.smodule    = modules_glob[param.scheme][param.module]
            param.image_type = types_glob[param.scheme][param.module]

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
            print('{} {} [{}]'.format(i, m, types_glob[s][i]))
        print()



def create_qsub_script(module):
    """Create bash script to be submitted with qsub.

    Parameters
    ----------
    module: string
        module name

    Returns:
    qsub_script_base: string
        job script base name
    """


    qsub_script_base = 'job_{}'.format(module)

    import platform
    hostname = platform.node()

    if 'candid' in hostname:
        create_qsub_script_candide(qsub_script_base, module)

    else:
        print('Warning: Using candide qsub job script for different machine')
        create_qsub_script_candide(qsub_script_base, module)

    return qsub_script_base



def create_qsub_script_candide(qsub_script_base, module):
    """Create a bash script to be submitted with qsub, working on iap:candide

    Parameters
    ----------
    qsub_script_base: string
        script name
    module: string
        module name (only used for qsub job name)

    Returns:
    None
    """

    qsub_script_name = '{}.sh'.format(qsub_script_base)

    f = open('{}'.format(qsub_script_name), 'w')
    print('#!/usr/bin/bash\n', file=f)
    print('#PBS -S /usr/bin/bash\n', file=f)
    print('#PBS -N {}\n'.format(module), file=f)
    print('#PBS -o {}.out'.format(qsub_script_base), file=f)
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

        package = '{}_package'.format(module)

        # Creat input dir (at the moment used for config files)
        path = '{}/{}/{}'.format(os.getcwd(), path_config[param.image_type], package)
        stuff.mkdir_p(path, verbose=param.verbose)

        # Copy config files: copy entire config tree. Need to remove existing one first.
        path_src  = '{}/{}/config'.format(path_spmod, package)
        path_dest = '{}/config'.format(path)
        if os.path.isdir(path_dest):
            shutil.rmtree(path_dest)
        shutil.copytree(path_src, path_dest)

        # Perform substitutions in package config file
        do_substitutions(path_dest, param.image_type)

        # Create output dir if necessary
        stuff.mkdir_p(path_output[param.image_type], verbose=param.verbose)

        launch_path = '{}/{}/config/launch.cmd'.format(path_spmod, package)

        # Run module
        if param.job == 'manual':
            my_env = os.environ.copy()
            my_env['CONFIG_DIR'] = path_dest
            ex, out, err = stuff.run_cmd(launch_path, run=not param.dry_run, verbose=param.verbose, devnull=False, env=my_env)
            if param.verbose:
                print('Sum[exit codes] = {}'.format(ex))
                if ex != 0:
                    print('output, error = ', out, err)

        elif param.job == 'qsub':

            qsub_script_base = create_qsub_script(module)
            stuff.run_cmd('qsub {}.sh'.format(qsub_script_base), run=not param.dry_run, verbose=param.verbose, devnull=False)

        else:
            stuff.error('Invalid job execution mode \'{}\''.format(param.job))

    

def adopt_run(module, image_type, verbose=False):
    """Adopt (last) run for given moduleo (set symbolic link).

    Parameters
    ----------
    module: string
        module name
    image_type: string
        image type, one in 'tile', 'exposure'
    verbose: bool, optional, default=False
        verbose output if True

    Returns
    -------
    None
    """

    path_runs = '{}/{}'.format(path_output[image_type], module)
    f = glob.glob('{}/run_*'.format(path_runs))
    f = sorted(f)
    last_run  = f[-1]
    source          = os.path.basename(last_run)
    source_to_check = last_run
    link_name = '{}/{}'.format(path_runs, name_adopted)

    stuff.ln_s(source, link_name, orig_to_check=source_to_check, verbose=verbose, force=True)


def discard_run(module, image_type, verbose=False):
    """Discard adopted run for given module (remove symbolic link)

    Parameters
    ----------
    module: string
        module name
    image_type: string
        image type, one in 'tile', 'exposure'
    verbose: bool, optional, default=False
        verbose output if True

    Returns
    -------
    None
    """

    path = '{}/{}/adopted'.format(path_output[image_type], module)
    if os.path.islink(path):
        if verbose:
            print('Discarding adopted run of module \'{}\''.format(module))
        os.remove(path)
    else:
        stuff.error('No adopted run found for module \'{}\''.format(module))



def set_results(module, image_type, verbose=False):
    """Set file links in <data> to adopted run for given module.

    Parameters
    ----------
    module: string
        module name
    image_type: string
        image type, one in 'tile', 'exposure'
    verbose: bool, optional, default=False
        verbose output if True

    Returns
    -------
    None
    """

    path_results = '{}/{}/{}/{}/{}'.format(os.getcwd(), path_output[image_type], module, name_adopted, name_results)

    files = glob.glob('{}/*'.format(path_results))
    for f in files:
        source    = f
        link_name = '{}/{}'.format(path_data[image_type], os.path.basename(f))
        stuff.ln_s(source, link_name, orig_to_check=source, verbose=verbose, force=True)


def do_substitutions(path_dest, image_type):
    """Perform basic substitutions of values for certain keys in package config file

    Parameters
    ----------
    path_dest: string
        path of config file
    image_type: string
        type of image, 'exposure' or 'image'

    Returns
    -------
    None
    """

    in_name  = '{}/{}'.format(path_dest, config_name)
    tmp_name = 'temp.tmp'
    fin  = open(in_name)
    fout = open(tmp_name, 'w')

    dat = fin.read()
    fin.close()

    # Set keys specifiy to image type (tile, exposure)
    dat = stuff.substitute(dat, 'BASE_DIR', '\$HOME/data', path_data[image_type])
    dat = stuff.substitute(dat, 'BASE_OUTPUT_DIR', 'output', path_output[image_type])

    if image_type == 'exposure':
        dat = stuff.substitute_arr(dat, 'INPUT_FILENAME_FORMATS', data_file_base['tile'], data_file_base['exposure'])
        dat = stuff.add_to_arr(dat, 'INPUT_FILENAME_FORMATS', '\'{}.flag.fits\''.format(data_file_base['exposure']))
    elif image_type == 'tile':
        pass


    fout.write(dat)
    fout.close()

    os.rename(tmp_name, in_name)



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
        adopt_run(param.smodule, param.image_type, verbose=param.verbose)

    elif param.mode == 'd':
        discard_run(param.smodule, param.image_type, verbose=param.verbose)

    elif param.mode == 's':
        set_results(param.smodule, param.image_type, verbose=param.verbose)

    ### End main program ###


    if param.verbose is True:
        print('End of program {}'.format(os.path.basename(argv[0])))


    return 0



if __name__ == "__main__":
    sys.exit(main(sys.argv))


