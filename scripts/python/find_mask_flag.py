# -*- coding: utf-8 -*-

"""FIND FLAG

This script look for mask flag in image headers.

:Authors: Axel Guinot

:Date: 07/03/2018

"""

import scatalog as sc
import numpy as np

import re
import argparse
import sys
import os


def get_files(path, file_name):
    """Get files

    Make a list of all the files to check.

    Parameters
    ----------
    path : str
        Path to the directory containing the file
    file_name : str
        Generic name of file to check (optional)

    Returns
    -------
    list
        List of the file to check.

    """

    l = []

    for i in os.listdir(path):
        file_ext = os.path.splitext(i)
        if file_ext[1] != '.fits':
            continue

        if file_name is not None:
            s = re.split('-', i)
            if s[0]+file_ext[1] == file_name:
                l.append(i)
            else:
                continue
        else:
            l.append(i)

    return l


def check_files(path, file_list):
    """Check files

    Read the header of the files in file_list and build a dictionary.

    Parameters
    ----------
    path : str
        Path to the directory containing the files
    file_list : list
        List of the files to check

    Returns
    -------
    dict
        Dictionary containing names, ratios and flags for each file.

    """

    d = {'file': [], 'ratio': [], 'flag': []}

    for i in file_list:
        f=sc.FITSCatalog(os.path.join(path, i), hdu_no= 0)
        f.open()
        h=f.get_header()
        f.close()

        d['file'].append(i)
        d['ratio'].append(h['MRATIO'])
        d['flag'].append(h['MFLAG'])

    return d


def save_to_file(dict, path, file_name):
    """Save to file

    Save the dictionary to a ASCII file.

    Parameters
    ----------
    dict : dict
        Dictionary to save
    path : str
        Path to the directory containing the files
    file_list : list
        List of the files to check

    Notes
    -----

    The file also contain the average ratio of files check.

    """

    if file_name == None:
        file_name = '*'
    f = open(path+'/check_flags.txt', 'w+')
    f.write('# Dir : {0}\n'.format(path))
    f.write('# File_name : {0}\n'.format(file_name))
    f.write('# Mean ratio : {0}\n'.format(np.mean(dict['ratio'])))
    f.write('# name   ratio   flag\n')
    for i in range(len(dict['file'])):
        f.write('{0}   {1}   {2}\n'.format(dict['file'][i], dict['ratio'][i], dict['flag'][i]))
    f.close()

def main(argv):
    """Main

    Main function

    Parameters
    ----------
    argv : sys.argv
        Arguments provide as input

    """

    parser = argparse.ArgumentParser(description= 'Get mask flag from image header')
    parser.add_argument('--path', help= 'Path to the directory containing the files to check')
    parser.add_argument('--file_name', help= "Generic name of files to check. Example : 'CFIS_flag-000-0.fits, CFIS_flag-001-0.fits, ...' -> 'CFIS_flag.fits' ")

    args = parser.parse_args()

    if args.path is not None:
        path = os.path.realpath(args.path)
    else:
        path = os.path.split(os.path.realpath(__file__))[0]

    if args.file_name is not None:
        file_name = args.file_name
    else:
        file_name = None

    file_list = get_files(path, file_name)

    d = check_files(path, file_list)

    save_to_file(d, path, file_name)


if __name__ == '__main__':
    main(argv = sys.argv)
