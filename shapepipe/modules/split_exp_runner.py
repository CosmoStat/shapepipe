# -*- coding: utf-8 -*-

"""SPLIT EXP RUNNER

This module splits the different CCDs (= hdus in FITS files) of a
single exposure into separate files.

:Author: Axel Guinot

:Date: 2019, 2020

:Package: ShapePipe

"""


from shapepipe.modules.module_decorator import module_runner
from shapepipe.modules.SplitExp_package import SplitExp_script as split


@module_runner(
    version='1.0',
    input_module='get_images_runner',
    file_pattern=['image', 'weight', 'flag'],
    file_ext=['.fz', '.fz', '.fz'],
    depends=['numpy', 'astropy', 'sip_tpv'],
    run_method='parallel'
)
def split_exp_runner(
    input_file_list,
    run_dirs,
    file_number_string,
    config,
    w_log
):

    output_suffix = config.getlist("SPLIT_EXP_RUNNER", "OUTPUT_SUFFIX")
    n_hdu = config.getint("SPLIT_EXP_RUNNER", "N_HDU")

    inst = split.SplitExposures(
        input_file_list,
        run_dirs['output'],
        file_number_string,
        output_suffix,
        n_hdu
    )

    inst.process()

    return None, None
