# -*- coding: utf-8 -*-

"""PSFEXINTERP RUNNER

This file is the pipeline runner for the PSFExInterpolation package.

:Author: Axel Guinot

"""

from shapepipe.modules.module_decorator import module_runner
from shapepipe.modules.PSFExInterpolation_package import interpolation_script as interp


@module_runner(input_module='setools_runner', version='1.0',
               file_pattern=['star_selection', 'galaxy_selection'],
               file_ext=['.psf', '.fits'], depends=['numpy', 'astropy', 'galsim'])
def psfexinterp_runner(input_file_list, output_dir, file_number_string,
                       config, w_log):

    psfcat_path, galcat_path = input_file_list
    
    pos_params = config.getlist('PSFEXINTERP_RUNNER', 'POSITION_PARAMS')
    get_shapes = config.getboolean('PSFEXINTERP_RUNNER', 'GET_SHAPES')
    

    inst = interp.PSFExInterpolator(psfcat_path, galcat_path, output_dir, pos_params, get_shapes)
    inst.write_output()

    return None, None