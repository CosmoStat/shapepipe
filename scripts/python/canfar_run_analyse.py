#!/usr/bin/env python3

"""Script canfar_run_analyse.py

Analyse a run (directory) on canfar using the job log files.

:Author: Martin Kilbinger

:Date: 05/2020
"""

import re
import os
import sys
import copy
from collections import Counter

from optparse import OptionParser


class param:
    """General class to store (default) variables
    """

    def __init__(self, **kwds):
        self.__dict__.update(kwds)

    def print(self, **kwds):
        print(self.__dict__)

    def var_list(self, **kwds):
        return vars(self)


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

    p_def = param(
        input_job = 'job_tile.sh',
    )

    return p_def


def parse_options(p_def):
    """Parse command line options.

    Parameters
    ----------
    p_def: class param
        parameter values

    Returns
    -------
    options: tuple
        Command line options
    args: string
        Command line string
    """

    usage  = "%prog [OPTIONS]"
    parser = OptionParser(usage=usage)

    # I/O
    parser.add_option('-i', '--input_job', dest='input_job', type='string', default=p_def.input_job,
         help='input job file, default=\'{}\''.format(p_def.input_job))
    parser.add_option('-o', '--output_fail', dest='output_fail', type='string',
         help='output file for failed jobs, none if not given')

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

    return True


def update_param(p_def, options):
    """Return default parameter, updated and complemented according to options.

    Parameters
    ----------
    p_def:  class param
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

    # Do extra stuff if necessary

    return param


# Results
res_ok = 0

res_unk = -1
res_wait = -2
res_subm = -3
res_abort = -4
res_err_mod = -5

res_noout = 1

# Failures
fail_unknown = -10
fail_vos_not_found = -11
fail_corrupt_fits = -12
fail_time_out = -13
fail_res_err_mod = -14


def get_status(tile_num):

    base_name = 'log_canfar_sp_'
    #base_name = 'log_sp_tile_'

    log_name = '{}{}.log'.format(base_name, tile_num)
    out_name = '{}{}.out'.format(base_name, tile_num)
    err_name = '{}{}.err'.format(base_name, tile_num)

    status = res_unk, 'unknown status'

    if not os.path.exists(log_name):
        status = res_wait,  'waiting for submission'
    else:
        if os.path.exists(out_name):
            final_cat_found = False
            with open(out_name) as out_file:
                for line in out_file:
                    m = re.match('Upload final_cat to.*(\d+) files', line)
                    if m:
                        final_cat_found = True
                        if m[1] == '0':
                            status = res_noout, '0 final_cat output files'
                        else:
                            status = res_ok, 'success, {} final_cat output files'.format(m[1])
                        break

            if final_cat_found == False:
                status = res_unk, 'Failed before final_cat'

                # Look for known errors in error log file
                with open(err_name) as err_file:
                    for line_err in err_file:
                        mm = re.search('NodeNotFound', line_err)
                        if mm:
                            status = status + (fail_vos_not_found, 'vos file not found')
                            break
                        mm = re.search('Empty or corrupt FITS file', line_err)
                        if mm:
                            status = status + (fail_corrupt_fits, 'corrupt input FITS file')
                            break
                        mm = re.search('ERROR:: HTTPSConnectionPool.*Max retries', line_err)
                        if mm:
                            status = status + (fail_time_out, 'vos time out')
                            break

                if len(status) == 2:

                    # Look for some errors in output file (pipeline error message)
                    module_last = None
                    with open(out_name) as out_file:
                        for line in out_file:
                            mmm = re.search('\- Module: (.*)', line)
                            if mmm:
                                module_last = mmm[1]
                            mmm = re.search('A total of (\d+) errors were recorded.', line)
                            if mmm and int(mmm[1]) != 0:
                                status = status + (fail_res_err_mod, '{} errors recorded, last module was {}'.format(mmm[1], module_last))
                                break

                if len(status) == 2:
                    status = status + (fail_unknown, 'unknown error')


        else:
            log_file = open(log_name)
            lines = log_file.readlines()
            log_file.close()
            for line in lines:
                if re.match('.*aborted', line):
                    status = res_abort, 'job aborted by user'
                    return status
            for line in lines:
                if re.match('.*executing', line):
                    status = res_wait, 'job running'
                    return status
            status = res_subm, 'job submitted, not started yet'

    return status


def output(status):

    print('## Issues')
    n_issue = 0
    for tile_num in sorted(status.keys()):
        if status[tile_num][0] == res_noout or status[tile_num][0] == res_unk:
            print('   ', tile_num, status[tile_num])
            n_issue = n_issue + 1
    if n_issue == 0:
        print('   none')

    hist = Counter(status.values())
    print('## Summary')
    print('# Nb: status (code)') 
    for s in hist:
        #print('{:2d}, {}: {}'.format(int(s[0]), s[1], hist[s]))
        print('{:6d}: {} ({})'.format(hist[s], s[1], int(s[0])), end='')
        if len(s) == 4:
            print('; {} ({})'.format(s[3], int(s[2])), end='')
        print()
        #print(hist[s], s)


def output_failed(output_fail, status):

    if output_fail:

        with open(output_fail, 'w') as f_out:
            for tile_num in status.keys():
                if status[tile_num][0] == res_noout or status[tile_num][0] == res_unk:
                    print(tile_num, file=f_out)
            

def main(argv=None):

    # Set default parameters
    p_def = params_default()

    # Command line options
    options, args = parse_options(p_def)
    # Without option parsing, this would be: args = argv[1:]

    if check_options(options) is False:
        return 1

    param = update_param(p_def, options)


    ### Start main program ###

    job_tile_path = param.input_job

    status = {}
    with open(job_tile_path) as job_file:
        for line in job_file:
            m = re.match('(\d{3}\.\d{3})', line)
            if m:
                tile_num = m[0]

                if len(argv) == 2 and argv[1] != tile_num:
                    continue
                status[tile_num] = get_status(tile_num)

    output(status)

    output_failed(param.output_fail, status)


if __name__ == "__main__":
    sys.exit(main(sys.argv))

