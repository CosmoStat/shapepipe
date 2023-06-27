"""MERGE HEADERS.

This module merges the output *header* files of the ``split_exp``
module or external scamp headers. It creates a binary file that 
contains the WCS of each CCD for each exposure.

:Author: Axel Guinot, Lucie Baumont

"""

import os
import re

import numpy as np
from sqlitedict import SqliteDict
from astropy.io.fits import Header
from astropy.wcs import WCS

def merge_headers(input_file_list, output_dir, pattern, ext_header, n_hdu):
    """Merge Headers.

    This function opens the files in the input file list and merges them into
    a SqliteDict file provided they match the appropriate pattern.

    Parameters
    ----------
    input_file_list : list
        List of input files
    output_dir : str
        Output path
    pattern : str
        File pattern
    ext_header : bool
        Use external scamp header if ``True``
    n_hdu : int
        number of ccds
    Raises
    ------
    TypeError
        For invalid ``output_dir`` type

    """
    if not isinstance(output_dir, str):
        raise TypeError(
            'Output directory for merge headers must be a string '
            + f'not {type(output_dir)}.'
        )

    # Open SqliteDict file
    final_file = SqliteDict(f'{output_dir}/log_exp_headers.sqlite')
    # define file properties for bookkeeping
    header_dir = os.path.split(input_file_list[0][0])[0]
    ext = os.path.splitext(input_file_list[0][0])[1]
    keys = unique_exposure_list(input_file_list, header_dir, pattern, ext, ext_header)
       
    for key in keys:
        # Load header, numpy binary or external header
        if ext_header:
            final_file[key] = create_joint_header(n_hdu, key, header_dir, pattern, ext)
        else:
            file_path = datadir+pattern+"key"+ext
            final_file[key] = np.load(file_path, allow_pickle=True)
    # Close file
    final_file.commit()
    final_file.close()

def unique_exposure_list(input_file_list, header_dir, file_pattern, ext, ext_header):
    """unique_exposure_list.

    Extract unique exposure ids from file list.

    Parameters
    ----------
    input_file_list : list
        list of input files
    header_dir : str
        directory of headers
    file_pattern : str
        text in filename 
    ext : str
        file extension, eg .npy, .head
    ext_header : bool
        Uses external scamp header if ``True`` 

    Returns
    -------
    list
        List of unique exposure numbers
    """
    file_list = []
    if ext_header:
        ext = '\-\d{1,2}' + ext 
    # extract keys
    for file in input_file_list:
        file_path_scalar = file[0]
        file_name = os.path.split(file_path_scalar)[1]
        root = re.split(ext, file_name)[0]
        pattern_split = re.split(file_pattern, root)
        if len(pattern_split) < 2:
            raise IndexError(
                f'Regex "{pattern}" not found in base name "{file_base_name}".'
            )
        key = pattern_split[1]
        
        file_list.append(key)
    # only keep unique keys
    unique_keys = list(set(file_list))

    return unique_keys

def create_joint_header(n_hdu, key, header_dir, pattern, ext):
    """create_joint_header.

    Packages scamp headers for ccds into a numpy array per exposure.
     
    Parameters
    ----------
    n_hdu : int
        total number of ccds
    key : str
        exposure number
    header_dir : str
        directory where headers are stored
    pattern : str
        file pattern, e.g. headers-
    ext : str
        header extension 
    Returns
    -------
    list
        compilation of headers corresponding to exposure number
    """

    header_file = np.zeros(n_hdu, dtype='O')
    for idx in range(1, n_hdu + 1):
        filepath= header_dir+'/'+pattern+key+'-'+str(idx)+ext
        print(filepath)
        # check for file                                                
        if os.path.exists(filepath):
            h=Header.fromfile(filepath,sep='\n',endcard=False,padding=False)
            try:
                w = WCS(h)
            except Exception:
                print(f'WCS error for file {exp_path}')
                raise
            header_file[idx - 1] = {'WCS': w, 'header': h.tostring()}
        else:
            header_file[idx - 1] = 'empty'

    return header_file
