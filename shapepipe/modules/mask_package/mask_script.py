# -*- coding: utf-8 -*-

"""MASK SCRIPT

This module contain a class to create star mask for an image.

:Authors: Axel Guinot, Martin Kilbinger

:Date: 20/12/2017

:Version: 1.1

"""

import shapepipe.pipeline.file_io as sc
from shapepipe.pipeline.config import CustomParser
from shapepipe.pipeline.execute import execute
from shapepipe.utilities.file_system import mkdir

import numpy as np
from astropy.coordinates import SkyCoord
from astropy import wcs
from astropy import units as u
import re
import os


class mask(object):
    """Mask class

    Class to create mask based on a star catalog.

    Parameters
    ----------
    image_path : str
        Path to image (fits format)
    weight_path : str
        Path to the weight image (fits format)
    config_filepath : str
        Path to the *.mask config file
    output_dir : str
        Path to the output directory
    hdu : int, optional, default = 0
        HDU number

    """

    def __init__(self, image_path, weight_path, image_suffix, image_num, config_filepath,
                 output_dir, path_external_flag=None, outname_base='flag',
                 star_cat_path=None, hdu=0):

        self._image_fullpath = image_path                                       # Path to the image to mask
        self._weight_fullpath = weight_path                                     # Path to the weight associated to the image
        self._config_filepath = config_filepath
        self._output_dir = output_dir                                           # Path to the output directory
        self._path_external_flag = path_external_flag                           # Path to an external flag file
        self._outname_base = outname_base                                       # Output file base name

        self._img_number = image_num
        if (image_suffix.lower() != 'none') & (image_suffix != ''):
            self._img_suffix = image_suffix + '_'
        else:
            self._img_suffix = ''

        self._star_cat_path = None
        if star_cat_path is not None:
            self._star_cat_path = star_cat_path

        self._hdu = hdu

        self._get_config()

        self._set_parameters()                                       # Set parameters needed for the stars detection

        self._err = False

    def _get_config(self):
        """Get config value

        Read the config file and set parameters.

        Raises
        ------
        ValueError : not config file name
        IOError : no config file

        """

        if self._config_filepath is None:
            raise ValueError('No path to config file')

        if not os.path.exists(self._config_filepath):
            raise IOError('Config file \'{}\' not found'
                          ''.format(self._config_filepath))

        conf = CustomParser()
        conf.read(self._config_filepath)

        self._config = {'PATH': {},
                        'BORDER': {},
                        'HALO': {},
                        'SPIKE': {},
                        'MESSIER': {},
                        'MD': {},
                        'EXTERNAL_FLAG': {}}

        self._config['PATH']['WW'] = conf.getexpanded('PROGRAM_PATH', 'WW_PATH')
        self._config['PATH']['WW_configfile'] = conf.getexpanded('PROGRAM_PATH', 'WW_CONFIG_FILE')
        # self._config['PATH']['CDSclient'] = conf.getexpanded('PROGRAM_PATH', 'CDSCLIENT_PATH')
        if conf.has_option('PROGRAM_PATH', 'CDSCLIENT_PATH'):
            self._config['PATH']['CDSclient'] = conf.getexpanded('PROGRAM_PATH', 'CDSCLIENT_PATH')
        elif self._star_cat_path is not None:
            self._config['PATH']['star_cat'] = self._star_cat_path
        else:
            raise ValueError('Either CDSCLIENT_PATH or STAR_CAT needs to be given in the [PROGRAM_PATH] section of the mask config file')

        self._config['PATH']['temp_dir'] = self._get_temp_dir_path(conf.getexpanded('OTHER', 'TEMP_DIRECTORY'))
        self._config['BORDER']['make'] = conf.getboolean('BORDER_PARAMETERS', 'BORDER_MAKE')
        if self._config['BORDER']['make']:
            self._config['BORDER']['width'] = conf.getint('BORDER_PARAMETERS', 'BORDER_WIDTH')
            self._config['BORDER']['flag'] = conf.get('BORDER_PARAMETERS', 'BORDER_FLAG_VALUE')
        for mask_shape in ['HALO', 'SPIKE']:
            self._config[mask_shape]['make'] = conf.getboolean(mask_shape + '_PARAMETERS', mask_shape + '_MAKE')
            self._config[mask_shape]['individual'] = conf.getboolean('OTHER', 'KEEP_INDIVIDUAL_MASK')
            if self._config[mask_shape]['make']:
                self._config[mask_shape]['maskmodel_path'] = conf.getexpanded(mask_shape+'_PARAMETERS', mask_shape+'_MASKMODEL_PATH')
                self._config[mask_shape]['mag_lim'] = conf.getfloat(mask_shape+'_PARAMETERS', mask_shape+'_MAG_LIM')
                self._config[mask_shape]['scale_factor'] = conf.getfloat(mask_shape+'_PARAMETERS', mask_shape+'_SCALE_FACTOR')
                self._config[mask_shape]['mag_pivot'] = conf.getfloat(mask_shape+'_PARAMETERS', mask_shape+'_MAG_PIVOT')
                self._config[mask_shape]['flag'] = conf.getint(mask_shape+'_PARAMETERS', mask_shape+'_FLAG_VALUE')
                if conf.getboolean('OTHER', 'KEEP_REG_FILE'):
                    self._config[mask_shape]['reg_file'] = self._config['PATH']['temp_dir'] + '/{0}{1}.reg'.format(re.split(".reg", conf.getexpanded(mask_shape+'_PARAMETERS', mask_shape+'_REG_FILE'))[0], self._img_number)
                else:
                    self._config[mask_shape]['reg_file'] = None
        self._config['MESSIER']['make'] = conf.getboolean('MESSIER_PARAMETERS', 'MESSIER_MAKE')
        if self._config['MESSIER']['make']:
            self._config['MESSIER']['cat_path'] = conf.getexpanded('MESSIER_PARAMETERS', 'MESSIER_CAT_PATH')
            self._config['MESSIER']['pixel_scale'] = conf.getfloat('MESSIER_PARAMETERS', 'MESSIER_PIXEL_SCALE')
            self._config['MESSIER']['size_plus'] = conf.getfloat('MESSIER_PARAMETERS', 'MESSIER_SIZE_PLUS')
            self._config['MESSIER']['flag'] = conf.getint('MESSIER_PARAMETERS', 'MESSIER_FLAG_VALUE')
        self._config['MD']['make'] = conf.getboolean('MD_PARAMETERS', 'MD_MAKE')
        if self._config['MD']['make']:
            self._config['MD']['thresh_flag'] = conf.getfloat('MD_PARAMETERS', 'MD_THRESH_FLAG')
            self._config['MD']['thresh_remove'] = conf.getfloat('MD_PARAMETERS', 'MD_THRESH_REMOVE')
            self._config['MD']['remove'] = conf.getboolean('MD_PARAMETERS', 'MD_REMOVE')
        self._config['EXTERNAL_FLAG']['make'] = conf.getboolean('EXTERNAL_FLAG', 'EF_MAKE')
        if self._config['EXTERNAL_FLAG']['make']:
            if self._path_external_flag is None:
                raise ValueError('External flag file has to be provided')
            self._config['EXTERNAL_FLAG']['path'] = self._path_external_flag

    def _set_parameters(self):
        """Set parameters

        Set the parameters for the stars detection.

        """

        img = sc.FITSCatalog(self._image_fullpath, hdu_no=0)
        img.open()
        self._header = img.get_header()
        img_shape = img.get_data().shape
        img.close()
        del(img)

        self._wcs = wcs.WCS(self._header)

        # Note: get_data().shape corresponds to (n_y, n_x)
        pix_center = [img_shape[1]/2., img_shape[0]/2.]
        wcs_center = self._wcs.all_pix2world([pix_center], 1)[0]

        self._fieldcenter = {}
        self._fieldcenter['pix'] = np.array(pix_center)
        self._fieldcenter['wcs'] = SkyCoord(ra=wcs_center[0], dec=wcs_center[1], unit='deg')

        self._img_radius = self._get_image_radius()

    def make_mask(self):
        """Make mask

        Main function to create the mask.

        """

        if self._config['MD']['make']:
            self.missing_data()

        if self._config['HALO']['make'] | self._config['SPIKE']['make']:
            stars = self.find_stars(np.array([self._fieldcenter['wcs'].ra.value, self._fieldcenter['wcs'].dec.value]), radius=self._img_radius)

        if not self._err:
            for i in ['HALO', 'SPIKE']:
                if self._config[i]['make']:
                    self._create_mask(stars=stars, types=i, mag_limit=self._config[i]['mag_lim'], scale_factor=self._config[i]['scale_factor'], mag_pivot=self._config[i]['mag_pivot'])

        if not self._err:
            if self._config['BORDER']['make']:
                border_mask = self.mask_border(width=self._config['BORDER']['width'])
            else:
                border_mask = None

        if not self._err:
            mask_name = []
            if self._config['HALO']['make'] & self._config['SPIKE']['make']:
                self._exec_WW(types='ALL')
                mask_name.append(self._config['PATH']['temp_dir'] + 'halo_spike_flag' + self._img_number + '.fits')
                mask_name.append(None)
                # mask_name.append(self._config['PATH']['temp_dir'] + 'spike_flag' + self._img_number + '.fits')
            else:
                for i in ['HALO', 'SPIKE']:
                    if self._config[i]['make']:
                        self._exec_WW(types=i)
                        mask_name.append(self._config['PATH']['temp_dir'] + i.lower()+'_flag' + self._img_number + '.fits')
                    else:
                        mask_name.append(None)

        if not self._err:
            if self._config['MESSIER']['make']:
                messier_mask = self.mask_messier(self._config['MESSIER']['cat_path'], size_plus=self._config['MESSIER']['size_plus'], flag_value=self._config['MESSIER']['flag'])
            else:
                messier_mask = None

        if not self._err:
            try:
                im_pass = self._config['MD']['im_remove']
            except:
                im_pass = True

        if not self._err:
            if self._config['EXTERNAL_FLAG']['make']:
                path_external_flag = self._config['EXTERNAL_FLAG']['path']
            else:
                path_external_flag = None

        if not self._err:
            if im_pass:
                final_mask = self._build_final_mask(path_mask1=mask_name[0], path_mask2=mask_name[1],
                                                    border=border_mask, messier=messier_mask,
                                                    path_external_flag=path_external_flag)

                if not self._config['HALO']['individual']:
                    if mask_name[0] is not None:
                        # os.system('rm {0}'.format(mask_name[0]))
                        self._rm_fits1_stdout, self._rm_fits1_stderr = execute('rm {0}'.format(mask_name[0]))
                    if mask_name[1] is not None:
                        # os.system('rm {1}'.format(mask_name[1]))
                        self._rm_fits2_stdout, self._rm_fits2_stderr = execute('rm {0}'.format(mask_name[1]))

                self._mask_to_file(input_mask=final_mask, output_fullpath=self._output_dir+'/'+self._img_suffix+self._outname_base+self._img_number+'.fits')

        # Handle stdout / stderr
        general_stdout = '\nCDSClient\n' + self._CDS_stdout
        general_stderr = ''
        if self._CDS_stderr != '':
            general_stderr += '\nCDSClient\n' + self._CDS_stderr
        if hasattr(self, '_WW_stdout') | hasattr(self, '_WW_stdout'):
            general_stdout += '\n\nWeightWatcher\n' + self._WW_stdout
            if self._WW_stderr != '':
                general_stderr += '\n\nWeightWatcher\n' + self._WW_stderr
        if hasattr(self, '_rm_reg_stderr') | hasattr(self, '_rm_reg_stdout'):
            general_stdout += '\n\nrm reg file\n' + self._rm_reg_stdout
            if self._rm_reg_stderr != '':
                general_stderr += '\n\nrm reg file\n' + self._rm_reg_stderr
        if hasattr(self, '_rm_fits1_stderr') | hasattr(self, '_rm_fits1_stdout'):
            general_stdout += '\n\nrm fits1 file\n' + self._rm_fits1_stdout
            if self._rm_fits1_stderr != '':
                general_stderr += '\n\nrm fits1 file\n' + self._rm_fits1_stderr
        if hasattr(self, '_rm_fits2_stderr') | hasattr(self, '_rm_fits2_stdout'):
            general_stdout += '\n\nrm fits2 file\n' + self._rm_fits2_stdout
            if self._rm_fits2_stderr != '':
                general_stderr += '\n\nrm fits2 file\n' + self._rm_fits2_stderr

        return general_stdout, general_stderr

    def find_stars(self, position, radius=None):
        """Find stars

        Return GSC (Guide Star Catalog) objects for a field with center (ra,dec) and radius r.

        Parameters
        ----------
        position : numpy.ndarray
          Position of the center of the field
        radius : float
          Radius in which the query is done (in arcmin)

        Returns
        -------
        dict
          Stars dicotionnary for GSC objects in the field.

        """

        if 'CDSclient' in self._config['PATH']:

          ra = position[0]
          dec = position[1]
          if dec > 0.:
              sign = '+'
          else:
              sign = ''
          cmd_line = '{0} {1} {2}{3} -r {4} -n 1000000'.format(self._config['PATH']['CDSclient'], ra, sign, dec, radius)
          # self._w_log.info('Calling command \'{}\''.format(cmd_line))
          self._CDS_stdout, self._CDS_stderr = execute(cmd_line)

        elif 'star_cat' in self._config['PATH']:

          # self._w_log.info('Reading star catalogue file \'{}\''.format(self._config['PATH']['star_cat']))
          f = open(self._config['PATH']['star_cat'], 'r')
          self._CDS_stdout = f.read()
          self._CDS_stderr = ''
          f.close()

        else:

          raise ValueError('Either CDSCLIENT_PATH or STAR_CAT needs to be given in the [PROGRAM_PATH] section of the mask config file')

        if self._CDS_stderr != '':
          self._err = True
          return None

        return self._make_star_cat(self._CDS_stdout)

    def mask_border(self, width=100, flag_value=4):
        """Create mask border

        Mask 'width' pixels around the image.

        Parameters
        ----------
        width : int
            Width of the mask mask border
        flag_value : int
            Value of the flag for the border (power of 2)

        Returns
        -------
        numpy.ndarray
            Array containing the mask.

        """

        if width is None:
            raise ValueError('Width not provided')

        # Note that python image array is [y, x]
        flag = np.zeros((int(self._fieldcenter['pix'][1]*2), int(self._fieldcenter['pix'][0]*2)), dtype='uint16')

        flag[0:width, :] = flag_value
        flag[-width:, :] = flag_value
        flag[:, 0:width] = flag_value
        flag[:, -width:] = flag_value

        return flag

    def mask_messier(self, cat_path, size_plus=0.1, flag_value=8):
        """Create mask for Messier objects

        Create a circular patch for Messier objects.

        Parameters
        ----------
        cat_path : str
            Path to the Messier catalog
        size_plus : float
            Increase the size of the mask by this factor (Example : 0.1 means 10%)
        flag_value : float
            Value of the flag, some power of 2

        Returns
        -------
        numpy.ndarray/None
            If no Messier objectds find in the field return None and the flag map otherwise.

        """

        if size_plus < 0:
            raise ValueError('size_plus has to be larger than 0')

        if cat_path is None:
            raise ValueError('cat_path has to be provided')

        m_cat = np.load(cat_path, allow_pickle=True)
        m_sc = SkyCoord(ra=m_cat['ra'] * u.degree, dec=m_cat['dec'] * u.degree)

        nx = self._fieldcenter['pix'][0] * 2
        ny = self._fieldcenter['pix'][1] * 2

        # Get the four corners of the image
        corners = self._wcs.calc_footprint()
        corners_sc = SkyCoord(ra=corners[:, 0] * u.degree, dec=corners[:, 1] * u.degree)

        # Loop through all Messier objects and check whether any corner is closer
        # than the object's radius
        ind = []
        for i, m_obj in enumerate(m_cat):
            r = max(m_obj['size']) * u.arcmin
            r_deg = r.to(u.degree)
            if np.any(corners_sc.separation(m_sc[i]) < r_deg):
                ind.append(i)

        if len(ind) == 0:
            # No closeby Messier object found
            return None

        # Note: python image array is [y, x]
        flag = np.zeros((int(self._fieldcenter['pix'][1]*2), int(self._fieldcenter['pix'][0]*2)), dtype='uint16')

        for i in ind:
            m_center = np.hstack(self._wcs.all_world2pix(m_cat['ra'][i], m_cat['dec'][i], 0))
            r_pix = max(m_cat['size'][i])/60. * (1 + size_plus) / np.abs(self._wcs.pixel_scale_matrix[0][0])

            # The following accounts for Messier centers outside of image, without
            # creating masks for coordinates out of range
            y_c, x_c = np.ogrid[0:ny, 0:nx]
            mask_tmp = (x_c - m_center[0])**2 + (y_c - m_center[1])**2 <= r_pix**2

            flag[mask_tmp] = flag_value

        return flag

    def missing_data(self):
        """Find missing data

        Look for 0 value in the image and flag it depending of the configuration.

        """

        img = sc.FITSCatalog(self._image_fullpath, hdu_no=0)
        img.open()

        im_shape = img.get_data().shape
        tot = float(im_shape[0] * im_shape[1])

        missing = float(len(np.where(img.get_data() == 0.)[0]))

        self._ratio = missing/tot

        if self._ratio >= self._config['MD']['thresh_flag']:
            self._config['MD']['im_flagged'] = True
        else:
            self._config['MD']['im_flagged'] = False
        if self._config['MD']['remove']:
            if self._ratio >= self._config['MD']['thresh_remove']:
                self._config['MD']['im_remove'] = True
                for i in ['HALO', 'SPIKE', 'MESSIER', 'BORDER']:
                    self._config[i]['make'] = False
            else:
                self._config['MD']['im_remove'] = False

        img.close()

    def SphereDist(self, position1, position2):
        """Compute spherical distance

        Compute spherical distance between 2 points.

        Parameters
        ----------
        position1 : numpy.ndarray
            [x,y] first point (in pixel)
        position2 : numpy.ndarray
            [x,y] second point (in pixel)

        Returns
        -------
        float
            The distance in degree.

        """

        if (type(position1) is not np.ndarray) & (type(position2) is not np.ndarray):
            raise ValueError('Positions need to be a numpy.ndarray')

        p1 = (np.pi/180.)*np.hstack(self._wcs.all_pix2world(position1[0], position1[1], 1))
        p2 = (np.pi/180.)*np.hstack(self._wcs.all_pix2world(position2[0], position2[1], 1))

        dTheta = p1 - p2
        dLong = dTheta[0]
        dLat = dTheta[1]

        dist = 2*np.arcsin(np.sqrt(np.sin(dLat/2.)**2. + np.cos(p1[1])*np.cos(p2[1])*np.sin(dLong/2.)**2.))

        return dist*(180./np.pi)*3600.

    def _get_image_radius(self, center=None):
        """Get image radius

        Compute the diagonal distance of the image in arcmin.

        Parameters
        ----------
        center : numpy.ndarray
            Coordinates of the center of the image (in pixel)

        Returns
        -------
        float
            The diagonal distance of the image in arcmin.

        """

        if center is None:
            return self.SphereDist(self._fieldcenter['pix'], np.zeros(2))/60.
        else:
            if type(center) is np.ndarray:
                return self.SphereDist(center, np.zeros(2))/60.
            else:
                raise TypeError('center has to be a numpy.ndarray')

    def _make_star_cat(self, CDSclient_output):
        """Make star catalog

        Make a dicotionnary from 'findgsc2.2' output.

        Parameters
        ----------
        CDSclient_output : str
            Output of 'findgsc2.2'

        Returns
        -------
        dict
            Stars dicotionnary containing all informations.

        """

        h = []
        stars = {}
        # get header
        for i in CDSclient_output.splitlines()[3].split(' '):
            if (i != '') & (i != ';'):
                # cleaning output
                i = i.replace(' ', '')
                for v in re.split(',|#|;', i):
                    if v != '':
                        i = v
                h.append(i)
                stars[i] = []

        # get data
        for i in range(4, len(CDSclient_output.splitlines()) - 5):
            k = 0
            for j in CDSclient_output.splitlines()[i].split(' '):
                if (j != '') & (j != ';'):
                    # cleaning output
                    j = j.replace(' ', '')
                    for v in re.split(',|#|;', j):
                        if v != '':
                            j = v
                    # handle missing data
                    try:
                        j = float(j)
                        stars[h[k]].append(j)
                    except:
                        if j == '---':
                            stars[h[k]].append(None)
                        else:
                            stars[h[k]].append(j)
                    k += 1

        return stars

    def _create_mask(self, stars, types=None, mag_limit=18., mag_pivot=13.8, scale_factor=0.3):
        """Create mask

        Apply mask from model to stars and save into DS9 region file.

        Parameters
        ----------
        stars : dict
            Stars dictionary (output of find_stars)
        types : str
            Type of mask in ['HALO', 'SPIKE']
        mag_limit : float
            Higher magnitude to apply the mask
        mag_pivot : float
            Pivot magnitude for the model
        scale_factor : float
            Scaling for the model

        Returns
        -------
        None
        """

        if stars is None:
            raise ValueError('No stars catalog provided')

        if self._config[types]['reg_file'] is None:
            reg = self._config['PATH']['temp_dir'] + types.lower() + self._img_number + '.reg'
        else:
            reg = self._config[types]['reg_file']

        if types == 'HALO':
            mask_model = np.loadtxt(self._config['HALO']['maskmodel_path']).transpose()
            mask_reg = open(reg, 'w')
        elif types == 'SPIKE':
            mask_model = np.loadtxt(self._config['SPIKE']['maskmodel_path']).transpose()
            mask_reg = open(reg, 'w')
        else:
            ValueError("types need to be in ['HALO', 'SPIKE']")

        stars_used = [[], [], []]
        for ra, dec, Fmag, Jmag, Vmag, Nmag, clas in zip(stars['RA(J2000)'], stars['Dec(J2000)'], stars['Fmag'], stars['Jmag'], stars['Vmag'], stars['Nmag'], stars['Clas']):
            mag = 0.
            i = 0.
            if Fmag is not None:
                mag += Fmag
                i += 1.
            if Jmag is not None:
                mag += Jmag
                i += 1.
            if Vmag is not None:
                mag += Vmag
                i += 1.
            if Nmag is not None:
                mag += Nmag
                i += 1.
            if i == 0.:
                mag = None
            else:
                mag /= i

            if (ra is not None) & (dec is not None) & (mag is not None) & (clas is not None):
                if (mag < mag_limit) & (clas == 0):
                    scaling = 1. - scale_factor * (mag - mag_pivot)
                    pos = self._wcs.all_world2pix(ra, dec, 0)
                    stars_used[0].append(pos[0])
                    stars_used[1].append(pos[1])
                    stars_used[2].append(scaling)

        for i in range(len(stars_used[0])):
            poly = 'polygon('
            for x, y in zip(mask_model[0], mask_model[1]):
                angle = np.arctan2(y, x)
                ll = stars_used[2][i] * np.sqrt(x ** 2. + y ** 2.)
                xnew = ll * np.cos(angle)
                ynew = ll * np.sin(angle)

                poly = poly + str(stars_used[0][i] + xnew + 0.5) + ' ' + str(stars_used[1][i] + ynew + 0.5) + ' '
            poly = poly + ')\n'
            mask_reg.write(poly)

        mask_reg.close()

    def _exec_WW(self, types=None):
        """Execute WeightWatcher

        Execute WeightWatcher to transform '.reg' to '.fits' flag map.

        Parameters
        ----------
        types : str
            Type of mask to make in ['HALO','SPIKE']

        """

        if types in ['HALO', 'SPIKE']:
            default_reg = self._config['PATH']['temp_dir'] + types.lower() + self._img_number + '.reg'
            defaul_out = self._config['PATH']['temp_dir'] + types.lower() + '_flag' + self._img_number + '.fits'
            if self._config[types]['reg_file'] is None:
                reg = default_reg
                if not sc.BaseCatalog(reg)._file_exists(reg):
                    raise sc.BaseCatalog.CatalogFileNotFound(reg)
                cmd = '{0} -c {1} -WEIGHT_NAMES {2} -POLY_NAMES {3} -POLY_OUTFLAGS {4} -FLAG_NAMES "" -OUTFLAG_NAME {5} -OUTWEIGHT_NAME ""'.\
                      format(self._config['PATH']['WW'], self._config['PATH']['WW_configfile'], self._weight_fullpath, reg, self._config[types]['flag'], defaul_out)
                # os.system(cmd)
                self._WW_stdout, self._WW_stderr = execute(cmd)
                # os.system('rm {0}'.format(reg))
                self._rm_reg_stdout, self._rm_reg_stderr = execute('rm {0}'.format(reg))
            else:
                reg = self._config[types]['reg_file']
                if not sc.BaseCatalog(reg)._file_exists(reg):
                    raise sc.BaseCatalog.CatalogFileNotFound(reg)
                cmd = '{0} -c {1} -WEIGHT_NAMES {2} -POLY_NAMES {3} -POLY_OUTFLAGS {4} -FLAG_NAMES "" -OUTFLAG_NAME {5} -OUTWEIGHT_NAME ""'.format(self._config['PATH']['WW'], self._config['PATH']['WW_configfile'], self._weight_fullpath, reg, self._config[types]['flag'], defaul_out)
                # os.system(cmd)
                self._WW_stdout, self._WW_stderr = execute(cmd)

        elif types == 'ALL':
            default_reg = [self._config['PATH']['temp_dir'] + 'halo' + self._img_number + '.reg', self._config['PATH']['temp_dir'] + 'spike' + self._img_number + '.reg']
            defaul_out = self._config['PATH']['temp_dir'] + 'halo_spike_flag' + self._img_number + '.fits'
            if self._config['HALO']['reg_file'] is None:
                reg = default_reg
                for i in range(2):
                    if not sc.BaseCatalog(reg[i])._file_exists(reg[i]):
                        raise sc.BaseCatalog.CatalogFileNotFound(reg[i])
                cmd = '{0} -c {1} -WEIGHT_NAMES {2} -POLY_NAMES {3},{4} -POLY_OUTFLAGS {5},{6} -FLAG_NAMES "" -OUTFLAG_NAME {7} -OUTWEIGHT_NAME ""'.format(self._config['PATH']['WW'], self._config['PATH']['WW_configfile'], self._weight_fullpath, reg[0], reg[1], self._config['HALO']['flag'], self._config['SPIKE']['flag'], defaul_out)
                # os.system(cmd)
                self._WW_stdout, self._WW_stderr = execute(cmd)
                # os.system('rm {0} {1}'.format(reg[0],reg[1]))
                self._rm_reg_stdout, self._rm_reg_stderr = execute('rm {0} {1}'.format(reg[0], reg[1]))
            else:
                reg = [self._config['HALO']['reg_file'], self._config['SPIKE']['reg_file']]
                for i in range(2):
                    if not sc.BaseCatalog(reg[i])._file_exists(reg[i]):
                        raise sc.BaseCatalog.CatalogFileNotFound(reg[i])
                cmd = '{0} -c {1} -WEIGHT_NAMES {2} -POLY_NAMES {3},{4} -POLY_OUTFLAGS {5},{6} -FLAG_NAMES "" -OUTFLAG_NAME {7} -OUTWEIGHT_NAME ""'.format(self._config['PATH']['WW'], self._config['PATH']['WW_configfile'], self._weight_fullpath, reg[0], reg[1], self._config['HALO']['flag'], self._config['SPIKE']['flag'], defaul_out)
                # os.system(cmd)
                self._WW_stdout, self._WW_stderr = execute(cmd)
        else:
            ValueError("types must be in ['HALO','SPIKE','ALL']")

        if (self._WW_stderr != '') | (self._rm_reg_stderr != ''):
            self._err = True

    def _build_final_mask(self, path_mask1, path_mask2=None, border=None, messier=None, path_external_flag=None):
        """Create final mask

        Create the final mask by combination of individual mask.

        Parameters
        ----------
        path_mask1 : str
            Path to a mask (fits format)
        path_mask2 : str
            Path to a mask (fits format)
        border : numpy.ndarray
            Array containing the border mask
        messier : numpy.ndarray
            Array containing the messier mask
        path_external_flag : str
            Path to an external flag file

        Returns
        -------
        numpy.ndarray
            Array containing the final mask.

        """

        final_mask = None

        if (path_mask1 is None) & (path_mask2 is None) & (border is None) & (messier is None):
            raise ValueError('No path to a mask, border or messier provided')

        if path_mask1 is not None:
            mask1 = sc.FITSCatalog(path_mask1, hdu_no=self._hdu)
            mask1.open()
            dat = mask1.get_data()
            final_mask = dat[:, :]

        if path_mask2 is not None:
            mask2 = sc.FITSCatalog(path_mask2, hdu_no=self._hdu)
            mask2.open()
            if final_mask is not None:
                final_mask += mask2.get_data()[:, :]
            else:
                final_mask = mask2.get_data()[:, :]

        if border is not None:
            if type(border) is np.ndarray:
                if final_mask is not None:
                    final_mask += border
                else:
                    final_mask = border
            else:
                raise TypeError('border has to be a numpy.ndarray')

        if messier is not None:
            if type(messier) is np.ndarray:
                if final_mask is not None:
                    final_mask += messier
                else:
                    final_mask = messier
            else:
                raise TypeError('messier has to be a numpy.ndarray')

        if path_external_flag is not None:
            external_flag = sc.FITSCatalog(path_external_flag, hdu_no=self._hdu)
            external_flag.open()
            if final_mask is not None:
                final_mask += external_flag.get_data()[:, :]
            else:
                final_mask = external_flag.get_data()[:, :]
            external_flag.close()

        return final_mask.astype(np.int16, copy=False)

    def _mask_to_file(self, input_mask, output_fullpath):
        """Mask to file

        Save the mask to a fits file.

        Parameters
        ----------
        input_mask : numpy.ndarray
            Mask to save
        output_fullpath : str
            Path of the output file

        """

        if input_mask is None:
            raise ValueError('input_mask not provided')
        if output_fullpath is None:
            raise ValueError('fullpath not provided')

        out = sc.FITSCatalog(output_fullpath, open_mode=sc.BaseCatalog.OpenMode.ReadWrite, hdu_no=0)
        out.save_as_fits(data=input_mask, image=True)

        if self._config['MD']['make']:
            out.open()
            out.add_header_card('MRATIO', self._ratio, 'ratio missing_pixels/all_pixels')
            out.add_header_card('MFLAG', self._config['MD']['im_flagged'], 'threshold value {:.3}'.format(self._config['MD']['thresh_flag']))

            # Write WCS information to header
            if self._wcs:
                header_wcs = self._wcs.to_header()
                for card in header_wcs:
                    out.add_header_card(card, header_wcs[card], header_wcs.comments[card])
            out.close()

    def _get_temp_dir_path(self, temp_dir_path):
        """Get temporary directory path

        Create the path and the directory for temporary file.

        Parameters
        -----------
        temp_dir_path : str
            Path to the temporary directory

        Returns
        -------
        str
            Path to the temporary directory

        Notes
        -----
        It is possible to have it in the output directory of the run, just
        enter 'OUTPUT'.

        """

        if temp_dir_path is None:
            raise ValueError('temp directory path not parovided')

        path = temp_dir_path.replace(' ', '')

        if path == 'OUTPUT':
            path = self._output_dir + '/temp'

        path += '/'
        if not os.path.isdir(path):
            mkdir(path)

        return path
