# -*- coding: utf-8 -*-
"""SHAPELENS SCRIPT

This script contain a class to handle shapelens software.

:Authors: Axel Guinot

:Date: 06/02/2017

"""

import numpy as np
import os
import subprocess
import re

import scatalog as sc


class shapelens(object):
    """Shapelens class

    Class to handle Shapelens software.

    Parameters
    ----------
    exec_path : str
        Path to the executable
    exec_option : str
        Option for the execution
    gal_cat_path : str
        Path to SExtractor catalog of galaxies
    psf_cat_path : str
        Path to PSF catalog
    output_dir : str
        Path to the output directory
    temp_dir : str
        Path to the temporary directory

    """

    def __init__(self, exec_path, exec_option, gal_cat_path, psf_cat_path, output_dir, temp_dir):

        self._gal_cat_path = gal_cat_path
        self._psf_cat_path = psf_cat_path

        self._exec_path = exec_path
        self._exec_mode = os.path.split(self._exec_path)[1]
        self._exec_option = exec_option

        self._output_dir = output_dir

        self._temp_dir = temp_dir
        if not os.path.isdir(self._temp_dir):
            os.system('mkdir {}'.format(self._temp_dir))

        s=re.split("\-([0-9]*)\-([0-9]+)\.",self._gal_cat_path)
        self._img_number='-{0}-{1}'.format(s[1],s[2])


    def process(self):
        """Process data

        Main function to process the data.

        """

        if self._exec_mode == 'get_shapes':
            self._set_data_format()
            self.get_shapes(self._gal_image_path, self._psf_image_path, self._stamp_size, self._nx, self._squared)
            self._save(self._shape_cat)
        else:
            raise ValueError("Exec mode has to be in ['get_shapes']")


    def get_shapes(self, gal_image_path, psf_image_path, stamp_size, n_stamp, remove_last):
        """Get shapes

        This function call the "get_shapes" function from Shapelens.

        Parameters
        ----------
        gal_image_path : str
            Path to the galaxies image with the vignet mapped
        psf_image_path : str
            Path to the PSFs image with the vignet mapped
        stamp_size : int
            Size of stamps (assume square stamp)
        n_stamp : int
            Number of vignet each dimension (assume square image)
        remove_last : int
            Number of stamp add to have a square image

        """

        self._shape_cat = {}
        for i in ['id', 'x', 'y', 'gamma1', 'gamma2', 'scale', 'SNR']:
            self._shape_cat[i] = []

        output = subprocess.check_output('{0} {1} -p {2} -s {3} -g {4} {5}'.format(self._exec_path, gal_image_path, psf_image_path, stamp_size, n_stamp, self._exec_option), shell=True)
        os.system('rm {}'.format(gal_image_path))
        os.system('rm {}'.format(psf_image_path))

        output = re.split('\n', output)
        output.remove('')
        if remove_last:
            output = output[:-remove_last]
        for i in output:
            items = re.split('\t', i)
            for j,key in zip(items, ['id', 'x', 'y', 'gamma1', 'gamma2', 'scale', 'SNR']):
                self._shape_cat[key].append(j)


    def _set_data_format(self):
        """Set format of data

        Map vignet on one single image.

        """

        gal_cat = sc.FITSCatalog(self._gal_cat_path, SEx_catalog=True)
        gal_cat.open()
        try:
            gal_vign = gal_cat.get_data()['VIGNET']
        except:
            raise ValueError('No vignet find in : {}'.format(self._gal_cat_path))
        gal_cat.close()

        sf_cat = sc.FITSCatalog(self._psf_cat_path, SEx_catalog=True)
        psf_cat.open()
        try:
            psf_vign = psf_cat.get_data()['VIGNET']
        except:
            raise ValueError('No psf find in : {}'.format(self._psf_cat_path))
        psf_cat.close()

        gal_image = self._map_vignet(gal_vign)
        self._gal_image_path = self._temp_dir + '/gal_temp{0}.fits'.format(self._img_number)
        psf_image = self._map_vignet(psf_vign)
        self._psf_image_path = self._temp_dir + '/psf_temp{0}.fits'.format(self._img_number)

        f_gal = sc.FITSCatalog(self._gal_image_path, open_mode = sc.BaseCatalog.OpenMode.ReadWrite)
        f_gal.save_as_fits(gal_image, image=True)
        f_psf = sc.FITSCatalog(self._psf_image_path, open_mode = sc.BaseCatalog.OpenMode.ReadWrite)
        f_psf.save_as_fits(psf_image, image=True)


    def _map_vignet(self, I):
        """Map vignet

        Map vignet on one single image.

        Parameters
        ----------
        I : numpy.ndarray
            Array of vinget to map

        """

        n_obj = I.shape[0]
        xs = I[0].shape[0]
        ys = I[0].shape[1]

        self._stamp_size = xs

        nx = int(np.sqrt(n_obj))
        if nx*nx != n_obj:
            nx += 1
            self._squared = nx*nx - n_obj
        else:
            self._squared = 0
        ny = nx
        self._nx = nx

        img_map=np.ones((xs*nx,ys*ny))

        ii=0
        jj=0
        for i in range(n_obj):
            if jj>nx-1:
                jj=0
                ii+=1
            img_map[ii*xs:(ii+1)*xs,jj*ys:(jj+1)*ys]=I[i]
            jj+=1

        return img_map


    def _save(self, cat):
        """Save

        Save output shear catalog

        """

        output_path = self._output_dir + '/' + self._exec_mode + self._img_number + '.fits'

        f = sc.FITSCatalog(output_path, open_mode= sc.BaseCatalog.OpenMode.ReadWrite, SEx_catalog=True)
        f.save_as_fits(cat, sex_cat_path= self._gal_cat_path)


    def _is_null(self, x):
        """Null-ness check

        Chek if the input is 0.

        Parameters
        ----------
        x : int or float

        Returns
        -------
        int (bool)
            If the input is 0 return 0 and 1 otherwise.
        """

        if x == 0:
            return 0
        else:
            return 1
