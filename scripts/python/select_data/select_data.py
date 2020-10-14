# -*- coding: utf-8 -*-

"""SELECT DATA

This module selection the data to process.

:Author: Axel Guinot

"""

from shapepipe.pipeline.config import CustomParser
from shapepipe.utilities.file_system import mkdir

from astropy import coordinates as c
from astropy import units as u

import os
import sys
import re
from glob import iglob

import numpy as np

from scipy.spatial import cKDTree

import matplotlib.pyplot as plt
from matplotlib import patches

class Selection(object):
    """ Selection

    This class contain method to select the single exposures in a specified area.
    It can also create the config files required for swarp in order to create coadd.

    Parameters
    ----------
    config_path : str
        Path to the config file.

    """

    def __init__(self, config_path):

        self.get_config(config_path)

        ra_single, dec_single, names_single = self.select_single()

        if self.config['COADD']['MAKE']:
            self.do_tiles(ra_single, dec_single, names_single)

    def get_config(self, config_path):
        """ Get config

        Read the config file and get parameters to run the script.

        Parameters
        ----------
        config_path : str
            Path to the config file.

        """
        config = CustomParser()
        config.read(config_path)

        self.config = {'SINGLE': {}, 'COADD': {}}

        self.config['SINGLE']['LOG_PATH'] = config.getexpanded('SINGLE_EPOCH_PARAMETERS', 'PATH_LOG')

        ra = config.getlist('SINGLE_EPOCH_PARAMETERS', 'RA')
        dec = config.getlist('SINGLE_EPOCH_PARAMETERS', 'DEC')
        self.config['SINGLE']['FIELD'] = np.array([[float(ra[0]), float(ra[1])],        # ra_min , ra_max      Deg
                                                   [float(dec[0]), float(dec[1])]])     # dec_min, dec_max     Deg

        self.config['SINGLE']['IMAGES'] = config.getboolean('SINGLE_EPOCH_PARAMETERS', 'GET_IMAGE')
        self.config['SINGLE']['WEIGHTS'] = config.getboolean('SINGLE_EPOCH_PARAMETERS', 'GET_WEIGHT')
        self.config['SINGLE']['FLAGS'] = config.getboolean('SINGLE_EPOCH_PARAMETERS', 'GET_FLAG')

        self.config['SINGLE']['LINK'] = config.getboolean('SINGLE_EPOCH_PARAMETERS', 'GET_LINK')
        self.config['SINGLE']['PLOT'] = config.getboolean('SINGLE_EPOCH_PARAMETERS', 'MAKE_PLOT')

        self.config['SINGLE']['OUT_DIR'] = config.getexpanded('SINGLE_EPOCH_PARAMETERS', 'OUTPUT_DIR')

        self.config['COADD']['MAKE'] = config.getboolean('COADD_PARAMETERS', 'MAKE_TILE')
        if self.config['COADD']['MAKE']:
            self.config['COADD']['SINGLE_DIR'] = config.getexpanded('COADD_PARAMETERS', 'PATH_SINGLE_DIR')

            ra = config.getlist('COADD_PARAMETERS', 'RA')
            dec = config.getlist('COADD_PARAMETERS', 'DEC')
            self.config['COADD']['FIELD'] = np.array([[float(ra[0]), float(ra[1])],        # ra_min , ra_max      Deg
                                                      [float(dec[0]), float(dec[1])]])     # dec_min, dec_max     Deg
            self.config['COADD']['SPACE'] = config.getfloat('COADD_PARAMETERS', 'SPACE')
            self.config['COADD']['SIZE'] = config.getint('COADD_PARAMETERS', 'SIZE')
            self.config['COADD']['PIXEL_SCALE'] = config.getfloat('COADD_PARAMETERS', 'PIXEL_SCALE')

            self.config['COADD']['PLOT'] = config.getboolean('COADD_PARAMETERS', 'MAKE_PLOT')

            self.config['COADD']['OUT_DIR'] = config.getexpanded('COADD_PARAMETERS', 'OUTPUT_DIR')

    def select_single(self):
        """ Select single

        Select the single exposures required.

        Returns
        -------
        Ra, Dec, names : list, list, list
            List of Ra, Dec center positions and names of the single exposures.

        """
        log_file = open(self.config['SINGLE']['LOG_PATH'])
        lines = log_file.readlines()
        log_file.close()

        ra_min, ra_max = self.config['SINGLE']['FIELD'][0]
        dec_min, dec_max = self.config['SINGLE']['FIELD'][1]

        # Select all exposures
        ra = []
        dec = []
        names = []
        seeing = []
        for l in lines:
            it = re.split('\|', l)

            name = it[0]
            name = name.replace(' ', '')

            tmp = re.split('\s+', it[3])
            ra_tmp = tmp[1]
            dec_tmp = tmp[2]
            coord = c.SkyCoord(ra_tmp, dec_tmp, unit=(u.hourangle, u.deg))

            tmp = re.split('\s+', it[7])[1]
            s = float(tmp)

            tmp = it[10]
            quality = re.split(' ', tmp)

            if (coord.ra.value > ra_max) | (coord.ra.value < ra_min) | (coord.dec.value > dec_max) | (coord.dec.value < dec_min):
                continue

            # Close to S. Gwyn criteria DR1 (version 2)
            # if (quality[0] not in ['P']) | (int(quality[1]) > 2) | (quality[2] != 'V') | (quality[3] != 'Q') | (quality[4] != 'R'):
            #     continue

            # J.C. Cuillandre criteria
            if (quality[2] != 'V'):
                continue

            names.append(name)
            ra.append(coord.ra.value)
            dec.append(coord.dec.value)
            seeing.append(s)

        # Remove repetition
        for ra_tmp, dec_tmp in zip(ra, dec):
            r = list(np.where((np.array(ra_tmp) == ra) & (np.array(dec_tmp) == dec))[0])
            if len(r) > 1:
                keep_ind = np.argmin(np.array(seeing)[np.array(r)])
                r.remove(r[keep_ind])
                for i in r:
                    del ra[i]
                    del dec[i]
                    del names[i]
                    del seeing[i]

        print('N exposures : {}'.format(len(ra)))

        # Write exposures to file
        if self.config['SINGLE']['IMAGES'] | self.config['SINGLE']['LINK']:
            self._write_single('images', names)
        if self.config['SINGLE']['WEIGHTS'] | self.config['SINGLE']['LINK']:
            self._write_single('weights', names)
        if self.config['SINGLE']['FLAGS'] | self.config['SINGLE']['LINK']:
            self._write_single('flags', names)

        # Plot single exposures
        if self.config['SINGLE']['PLOT']:
            self._plot_single(ra, dec)

        return ra, dec, names

    def do_tiles(self, ra_single, dec_single, names_single):
        """ Do tiles

        Create the config files for swarp to create the coadded images.

        Parameters
        ----------
        ra_single : list
            List of Ra center positions of single exposures.
        dec_single : list
            List of Ra center positions of single exposures.
        names_single : list
            List of single exposures names.

        """
        if not os.path.isdir(self.config['COADD']['OUT_DIR']):
            mkdir(self.config['COADD']['OUT_DIR'])

        output_dir = self.config['COADD']['OUT_DIR'] + '/tiles_config'
        if not os.path.isdir(output_dir):
            mkdir(output_dir)
        else:
            f_list = iglob(output_dir + '/tile_*')
            for f in f_list:
                os.remove(f)

        ra_min, ra_max = self.config['COADD']['FIELD'][0]
        dec_min, dec_max = self.config['COADD']['FIELD'][1]

        # Make tiles
        space = self.config['COADD']['SPACE']    # Deg
        size = self.config['COADD']['SIZE']     # Pix

        single_dir = self.config['COADD']['SINGLE_DIR']

        ra_corr = np.cos(np.mean(self.config['COADD']['FIELD'][1])*np.pi/180.)
        tree = cKDTree(np.array([np.array(ra_single) * ra_corr, dec_single]).T)
        count = 0
        n_exp = []
        for dec in np.arange(dec_min + space, dec_max, space):
            # ra_corr2 = np.cos(dec*np.pi/180.)
            for ra in np.arange(ra_min + space/ra_corr, ra_max, space/ra_corr):
                res = tree.query(np.array([ra * ra_corr, dec]), k=30, distance_upper_bound=1.2, n_jobs=-1)
                exp_names = np.array(names_single)[res[1][np.where(res[0] != np.inf)]]
                f_tile = open(output_dir + '/tile_{:.1f}_{:.1f}-{}.txt'.format(ra, dec, count), 'w')
                for n in exp_names:
                    f_tile.write(single_dir + '/image-{}.fits\n'.format(n))
                f_tile.close()
                n_exp.append(len(exp_names))
                count += 1

        print('N tiles : {}'.format(count))
        print('Average exp per tile : {}'.format(int(round(np.mean(n_exp)))))

        if self.config['COADD']['PLOT']:
            self._plot_coadd()

    def _write_single(self, img_type, names):
        """ Write single

        Write single exposures information to files.

        Parameters
        ----------
        img_type : str
            In ['images', 'weights', 'flags'].
        names : list
            List of exposure names.

        """
        if img_type not in ['images', 'weights', 'flags']:
            raise ValueError("img_type has to be in ['images', 'weights', 'flags']")

        if not os.path.isdir(self.config['SINGLE']['OUT_DIR']):
            mkdir(self.config['SINGLE']['OUT_DIR'])

        if self.config['SINGLE'][img_type.upper()]:
            f_select = open(self.config['SINGLE']['OUT_DIR'] + '/selected_{}.txt'.format(img_type), 'w')
        if self.config['SINGLE']['LINK']:
            f_link = open(self.config['SINGLE']['OUT_DIR'] + '/link_{}.txt'.format(img_type), 'w')

        if img_type == 'images':
            img_dir = 'pitcairn'
            suffix = ''
        else:
            img_dir = img_type
            suffix = '.' + img_type[:-1]

        for i in range(len(names)):
            if self.config['SINGLE'][img_type.upper()]:
                f_select.write(names[i] + 'p{}.fits.fz\n'.format(suffix))
            if self.config['SINGLE']['LINK']:
                f_link.write('https://www.cadc-ccda.hia-iha.nrc-cnrc.gc.ca/files/vault/cfis/{}/'.format(img_dir) + names[i] + 'p{}.fits.fz\n'.format(suffix))
        if self.config['SINGLE'][img_type.upper()]:
            f_select.close()
        if self.config['SINGLE']['LINK']:
            f_link.close()

    def _plot_single(self, ra, dec):
        """ Plot single exposures

        Plot the single exposures in the desired field.

        Parameters
        ----------
        ra : list
            List of Ra center positions.
        dec : list
            List of Ra center positions.

        """
        if not os.path.isdir(self.config['SINGLE']['OUT_DIR']):
            mkdir(self.config['SINGLE']['OUT_DIR'])

        ra_min, ra_max = self.config['SINGLE']['FIELD'][0]
        dec_min, dec_max = self.config['SINGLE']['FIELD'][1]

        ra_size = ra_max - ra_min
        dec_size = dec_max - dec_min
        ratio = dec_size/ra_size

        fig, ax = plt.subplots(figsize=(15,int(round(15*ratio))))
        plt.plot(ra, dec, '+', c='k')
        for xx, yy in zip(ra, dec):
            r = patches.Rectangle((xx-0.5, yy-0.5), 1, 1, fill=None)
            ax.add_patch(r)
        plt.xlim(ra_min-1, ra_max+1)
        plt.ylim(dec_min-1, dec_max+1)
        plt.xlabel('Ra (deg)')
        plt.ylabel('Dec (deg)')
        plt.gca().invert_xaxis()
        plt.savefig(self.config['SINGLE']['OUT_DIR'] + '/single_exp_field.png')
        plt.close(fig)

    def _plot_coadd(self):
        """ Plot coadd

        Plot the coadded images in the desired field.

        """
        if not os.path.isdir(self.config['COADD']['OUT_DIR']):
            mkdir(self.config['COADD']['OUT_DIR'])

        ra_min, ra_max = self.config['COADD']['FIELD'][0]
        dec_min, dec_max = self.config['COADD']['FIELD'][1]

        space = self.config['COADD']['SPACE']    # Deg
        size = self.config['COADD']['SIZE']      # Pix
        pixel_scale = self.config['COADD']['PIXEL_SCALE']   # arcsec/pix

        fig, ax= plt.subplots(figsize=(15,15))
        for dec in np.arange(dec_min + space, dec_max, space):
            ra_corr2 = np.cos(dec*np.pi/180.)
            for ra in np.arange(ra_min + space/ra_corr2, ra_max, space/ra_corr2):
                plt.plot(ra, dec, '+', c='k')
                rect = patches.Rectangle((ra-(space/2./ra_corr2), dec-space/2.),
                                      size*pixel_scale/3600/ra_corr2,
                                      size*pixel_scale/3600,
                                      fill=None, color="blue")
                ax.add_patch(rect)
        plt.xlabel('Ra (deg)')
        plt.ylabel('Dec (deg)')
        plt.gca().invert_xaxis()
        plt.savefig(self.config['COADD']['OUT_DIR'] + '/tiles_field.png')
        plt.close(fig)


if __name__ == '__main__':

    argv = sys.argv

    try:
        config_path = argv[1]
    except:
        raise ValueError('No config file provided')

    Selection(config_path)
