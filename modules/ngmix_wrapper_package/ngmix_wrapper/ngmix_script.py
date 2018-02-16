# -*- coding: utf-8 -*-
"""NGMIX SCRIPT

This script contain a class to run ngmix for metacalibration.

:Authors: Axel Guinot

:Date: 16/02/2017

"""

import numpy as np
import ngmix
import scatalog as sc
from astropy import stats

import os

import re

class ngmix_wrapper(object):
    """
    """

    def __init__(self, gal_cat_path, psf_cat_path, output_dir, option_dict):

        self._gal_cat_path = gal_cat_path
        self._psf_cat_path = psf_cat_path
        self._output_dir = output_dir
        self._option_dict = option_dict

        s=re.split("\-([0-9]{3})\-([0-9]+)\.",self._gal_cat_path)
        self._img_number='-{0}-{1}'.format(s[1],s[2])

        self._gal_vign = self._load_data(self._gal_cat_path)
        self._psf_vign = self._load_data(self._psf_cat_path)


    def _load_data(self, path):
        """
        """

        f=sc.FITSCatalog(path, SEx_catalog=True)
        f.open()
        try:
            return f.get_data()['VIGNET']
        except:
            raise TypeError("Data in the wrong format. Assume SExtractor like catalog with the vignet in 'VIGNET'")

    def _psf_fitter(self, psf_vign):
        """
        """

        psf_obs=ngmix.Observation(psf_vign)
        pfitter=ngmix.fitting.LMSimple(psf_obs,'gauss')

        shape = psf_vign.shape
        psf_pars = np.array([shape[0]/2., shape[1]/2., 0., 0., 1., 1.])
        pfitter.go(psf_pars)

        psf_gmix_fit=pfitter.get_gmix()
        psf_obs.set_gmix(psf_gmix_fit)

        return psf_obs

    def _get_weight(self, gal_vign):
        """
        """

        shape = gal_vign.shape

        std = stats.mad_std(gal_vign)

        w = 1./(np.random.normal(0.,std,(shape[0], shape[1])))**2.

        return w

    def process(self):
        """
        """

        n_obj = len(self._gal_vign)

        output_dict = {'psf': []}
        for i in self._option_dict['TYPES']:
            output_dict[i] = []

        for i in range(n_obj):
            temp_dict = self.make_metacal(self._gal_vign[i], self._psf_vign, self._option_dict)
            for j in self._option_dict['TYPES']:
                output_dict[j].append(temp_dict[j].image)
            try:
                output_dict['psf'].append(temp_dict['noshear'].get_psf().image)
            except:
                output_dict['psf'].append(temp_dict[j].get_psf().image)

        self._save(output_dict)


    def make_metacal(self, gal_vign, psf_vign, option_dict):
        """
        """

        psf_obs = self._psf_fitter(psf_vign)

        weight = self._get_weight(gal_vign)

        obs = ngmix.Observation(gal_vign, psf=psf_obs, weight=weight)

        obs_out = ngmix.metacal.get_all_metacal(obs,
                                                types= option_dict['TYPES'],
                                                fixnoise= option_dict['FIXNOISE'],
                                                cheatnoise= option_dict['CHEATNOISE'],
                                                symmetrize_psf= option_dict['SYMMETRIZE_PSF'],
                                                step= option_dict['STEP'])

        return obs_out

    def _save(self, output_dict):
        """
        """

        output_dir_path = {}
        for i in output_dict.keys():
            output_dir_path[i] = self._output_dir + '/' + i
            if not os.path.isdir(output_dir_path[i]):
                os.system('mkdir {}'.format(output_dir_path[i]))

        for i in output_dict.keys():
            output_file_path = output_dir_path[i] + '/' + 'ngmix_' + i + self._img_number + '.fits'
            f = sc.FITSCatalog(output_file_path, SEx_catalog= True, open_mode= sc.BaseCatalog.OpenMode.ReadWrite)
            f.save_as_fits(output_dict[i], ['VIGNET'], sex_cat_path= self._gal_cat_path)
