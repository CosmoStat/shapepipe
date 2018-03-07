# -*- coding: utf-8 -*-

"""MASK SCRIPT

This module contain a class to create star mask for an image.

:Authors: Axel Guinot

:Date: 20/12/2017

"""

import scatalog as sc
import sconfig

import numpy as np
import subprocess
import astropy.coordinates as coord
from astropy import wcs
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

    """

    def __init__(self, image_path, weight_path, config_filepath, output_dir):

        self._image_fullpath = image_path                                       # Path to the image to mask
        self._weight_fullpath = weight_path                                     # Path to the weight associated to the image
        self._config_filepath = config_filepath
        self._output_dir = output_dir                                           # Path to the output directory

        s=re.split("\-([0-9]{3})\-([0-9]+)\.",self._image_fullpath)
        self._img_number='-{0}-{1}'.format(s[1],s[2])                           # Needed for temporary file
        self._img_name = os.path.split(s[0])[1]

        self._get_config(self._config_filepath)                                 # Get parameters from config file

        self._set_parameters()                                       # Set parameters needed for the stars detection


    def _get_config(self, config_filepath):
        """Get config value

        Read the config file and set parameters.

        Parameters
        ----------
        config_filepath : str
            Path to the *.mask config file

        """

        if config_filepath is None:
            raise ValueError('No path to config file')

        conf=sconfig.SConfig(config_filepath)

        self._config={'PATH': {}, 'BORDER': {}, 'HALO': {}, 'SPIKE': {}, 'MESSIER': {}, 'MD': {}}

        self._config['PATH']['WW'] = conf.get_as_string('WW_PATH','PROGRAM_PATH')
        self._config['PATH']['WW_configfile'] = conf.get_as_string('WW_CONFIG_FILE','PROGRAM_PATH')
        self._config['PATH']['CDSclient'] = conf.get_as_string('CDSCLIENT_PATH','PROGRAM_PATH')
        self._config['PATH']['temp_dir'] = self._get_temp_dir_path(conf.get_as_string('TEMP_DIRECTORY','OTHER'))
        self._config['BORDER']['make'] = conf.get_as_boolean('BORDER_MAKE','BORDER_PARAMETERS')
        if self._config['BORDER']['make']:
            self._config['BORDER']['width'] = conf.get_as_int('BORDER_WIDTH','BORDER_PARAMETERS')
            self._config['BORDER']['flag'] = conf.get_as_string('BORDER_FLAG_VALUE','BORDER_PARAMETERS')
        for i in ['HALO','SPIKE']:
            self._config[i]['make'] = conf.get_as_boolean(i+'_MAKE',i+'_PARAMETERS')
            self._config[i]['individual'] = conf.get_as_boolean('KEEP_INDIVIDUAL_MASK','OTHER')
            if self._config[i]['make']:
                self._config[i]['maskmodel_path'] = conf.get_as_string(i+'_MASKMODEL_PATH',i+'_PARAMETERS')
                self._config[i]['mag_lim'] = conf.get_as_float(i+'_MAG_LIM',i+'_PARAMETERS')
                self._config[i]['scale_factor'] = conf.get_as_float(i+'_SCALE_FACTOR',i+'_PARAMETERS')
                self._config[i]['mag_pivot'] = conf.get_as_float(i+'_MAG_PIVOT',i+'_PARAMETERS')
                self._config[i]['flag'] = conf.get_as_int(i+'_FLAG_VALUE',i+'_PARAMETERS')
                if conf.get_as_boolean('KEEP_REG_FILE','OTHER'):
                    self._config[i]['reg_file'] = self._config['PATH']['temp_dir'] + '/{0}{1}.reg'.format(re.split(".reg",conf.get_as_string(i+'_REG_FILE',i+'_PARAMETERS'))[0],self._img_number)
                else:
                    self._config[i]['reg_file'] = None
        self._config['MESSIER']['make'] = conf.get_as_boolean('MESSIER_MAKE','MESSIER_PARAMETERS')
        if self._config['MESSIER']['make']:
            self._config['MESSIER']['cat_path'] = conf.get_as_string('MESSIER_CAT_PATH', 'MESSIER_PARAMETERS')
            self._config['MESSIER']['pixel_scale'] = conf.get_as_float('MESSIER_PIXEL_SCALE', 'MESSIER_PARAMETERS')
            self._config['MESSIER']['size_plus'] = conf.get_as_float('MESSIER_SIZE_PLUS', 'MESSIER_PARAMETERS')
            self._config['MESSIER']['flag'] = conf.get_as_int('MESSIER_FLAG_VALUE', 'MESSIER_PARAMETERS')
        self._config['MD']['make'] = conf.get_as_boolean('MD_MAKE', 'MD_PARAMETERS')
        if self._config['MD']['make']:
            self._config['MD']['thresh_flag'] = conf.get_as_float('MD_THRESH_FLAG', 'MD_PARAMETERS')
            self._config['MD']['thresh_remove'] = conf.get_as_float('MD_THRESH_REMOVE', 'MD_PARAMETERS')
            self._config['MD']['remove'] = conf.get_as_boolean('MD_REMOVE', 'MD_PARAMETERS')


    def _set_parameters(self):
        """Set parameters

        Set the parameters for the stars detection.

        """

        img = sc.FITSCatalog(self._image_fullpath, hdu_no=0)
        img.open()
        self._header = img.get_header()
        img.close()
        del(img)

        self._wcs = wcs.WCS(self._header)

        self._fieldcenter={}
        self._fieldcenter['pix']=np.array([self._header['CRPIX1'],self._header['CRPIX2']])
        self._fieldcenter['wcs']=coord.SkyCoord(ra=self._header['CRVAL1'], dec=self._header['CRVAL2'], unit='deg')

        self._img_radius=self._get_image_radius()


    def make_mask(self):
        """Make mask

        Main function to create the mask.

        """

        if self._config['MD']['make']:
            self.missing_data()

        if self._config['HALO']['make'] | self._config['SPIKE']['make']:
            stars=self.find_stars(np.array([self._fieldcenter['wcs'].ra.value,self._fieldcenter['wcs'].dec.value]), radius=self._img_radius)

        for i in ['HALO', 'SPIKE']:
            if self._config[i]['make']:
                self._create_mask(stars=stars, types=i, mag_limit=self._config[i]['mag_lim'], scale_factor=self._config[i]['scale_factor'], mag_pivot=self._config[i]['mag_pivot'])

        if self._config['BORDER']['make']:
            border_mask=self.mask_border(width=self._config['BORDER']['width'])
        else:
            border_mask=None

        mask_name=[]
        if self._config['HALO']['make'] & self._config['SPIKE']['make']:
            self._exec_WW(types='ALL')
            mask_name.append(self._config['PATH']['temp_dir'] + 'halo_spike_flag' + self._img_number + '.fits')
            mask_name.append(None)
            # mask_name.append(self._config['PATH']['temp_dir'] + 'spike_flag' + self._img_number + '.fits')
        else:
            for i in ['HALO','SPIKE']:
                if self._config[i]['make']:
                    self._exec_WW(types=i)
                    mask_name.append(self._config['PATH']['temp_dir'] + i.lower()+'_flag' + self._img_number + '.fits')
                else:
                    mask_name.append(None)

        if self._config['MESSIER']['make']:
            messier_mask = self.mask_messier(self._config['MESSIER']['cat_path'], size_plus= self._config['MESSIER']['size_plus'], flag_value= self._config['MESSIER']['flag'])
        else:
            messier_mask = None

        try:
            im_pass = self._config['MD']['im_remove']
        except:
            im_pass = True

        if im_pass:
            final_mask=self._build_final_mask(path_mask1=mask_name[0],path_mask2=mask_name[1], border=border_mask, messier=messier_mask)

            if not self._config['HALO']['individual']:
                if mask_name[0] is not None:
                    os.system('rm {0}'.format(mask_name[0]))
                if mask_name[1] is not None:
                    os.system('rm {1}'.format(mask_name[1]))

            self._mask_to_file(input_mask=final_mask, output_fullpath=self._output_dir+'/'+self._img_name+'_flag'+self._img_number+'.fits')


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

        ra=position[0]
        dec=position[1]

        #check ra dec types

        if dec>0. :
            sign='+'
        else:
            sign=''

        output=subprocess.check_output('{0} {1} {2}{3} -r {4} -n 1000000'.format(self._config['PATH']['CDSclient'], ra, sign, dec, radius), shell=True)

        return self._make_star_cat(output)


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

        flag = np.zeros((int(self._fieldcenter['pix'][0]*2),int(self._fieldcenter['pix'][1]*2)),dtype='uint8')

        flag[0:width,:]=flag_value
        flag[-width:-1,:]=flag_value
        flag[:,0:width]=flag_value
        flag[:,-width:-1]=flag_value

        return flag


    def mask_messier(self, cat_path, size_plus= 0.1, flag_value= 8):
        """Create mask Messier

        Create a circular patch for Messier objects.

        Parameters
        ----------
        cat_path : str
            Path to the Messier catalog
        size_plus : float
            Increase the size of the mask by this factor (Example : 0.1 means 10%)
        flag_value : intMessier objects (power of 2)
            Value of the flag for

        Returns
        -------
        numpy.ndarray/None
            If no Messier objectds find in the field return None and the flag map otherwise.

        """

        if cat_path == None:
            raise ValueError('cat_path has to be provided')

        if (size_plus < 0):
            raise ValueError('size_plus has to be in [0, 1]')

        nx = self._fieldcenter['pix'][0]*2
        ny = self._fieldcenter['pix'][1]*2

        m_cat = np.load(cat_path)

        ra_max = np.hstack(self._wcs.all_pix2world(0, self._fieldcenter['pix'][1], 1))[0]
        ra_min = np.hstack(self._wcs.all_pix2world(nx, self._fieldcenter['pix'][1], 1))[0]
        dec_min = np.hstack(self._wcs.all_pix2world(self._fieldcenter['pix'][0], 0, 1))[1]
        dec_max = np.hstack(self._wcs.all_pix2world(self._fieldcenter['pix'][0], ny, 1))[1]

        ind = np.where((m_cat['ra'] > ra_min) & (m_cat['ra'] < ra_max) & (m_cat['dec'] > dec_min) & (m_cat['dec'] < dec_max))[0]

        if len(ind) == 0:
            return None

        flag = np.zeros((int(self._fieldcenter['pix'][0]*2),int(self._fieldcenter['pix'][1]*2)),dtype='uint8')

        for i in ind:
            m_center = np.hstack(self._wcs.all_world2pix(m_cat['ra'][i], m_cat['dec'][i], 0))
            r_pix = max(m_cat['size'][i])/60. * (1 + size_plus) / np.abs(self._wcs.pixel_scale_matrix[0][0])
            y_c, x_c = np.ogrid[-int(m_center[1]):ny-int(m_center[1]), -int(m_center[0]):nx-int(m_center[0])]
            mask_tmp = x_c*x_c + y_c*y_c <= r_pix*r_pix
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

        if not (type(position1) is np.ndarray) & (type(position2) is np.ndarray):
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

        h=[]
        stars={}
        #get header
        for i in CDSclient_output.splitlines()[3].split(' '):
            if (i != '') & (i != ';'):
                #cleaning output
                i=i.replace(' ','')
                for v in re.split(',|#|;',i):
                    if v != '':
                        i=v
                h.append(i)
                stars[i]=[]

        #get data
        for i in range(4,len(CDSclient_output.splitlines())-5):
            k=0
            for j in CDSclient_output.splitlines()[i].split(' '):
                if (j != '') & (j != ';'):
                    #cleaning output
                    j=j.replace(' ','')
                    for v in re.split(',|#|;',j):
                        if v != '':
                            j=v
                    #handle missing data
                    try:
                        j=float(j)
                        stars[h[k]].append(j)
                    except:
                        if j == '---':
                            stars[h[k]].append(None)
                        else:
                            stars[h[k]].append(j)
                    k+=1

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

        """

        if stars is None:
            raise ValueError('No stars catalog provided')

        if self._config[types]['reg_file'] is None:
            reg = self._config['PATH']['temp_dir'] + types.lower()+ self._img_number + '.reg'
        else:
            reg = self._config[types]['reg_file']

        if types == 'HALO':
            mask_model = np.loadtxt(self._config['HALO']['maskmodel_path']).transpose()
            mask_reg = open(reg,'w')
        elif types == 'SPIKE':
            mask_model = np.loadtxt(self._config['SPIKE']['maskmodel_path']).transpose()
            mask_reg = open(reg,'w')
        else:
            ValueError("types need to be in ['HALO', 'SPIKE']")

        stars_used=[[],[],[]]
        for ra, dec, Fmag, Jmag, Vmag, Nmag, clas in zip(stars['RA(J2000)'], stars['Dec(J2000)'], stars['Fmag'], stars['Jmag'], stars['Vmag'], stars['Nmag'], stars['Clas']):
            mag=0.
            i=0.
            if Fmag!=None:
                mag+=Fmag
                i+=1.
            if Jmag!=None:
                mag+=Jmag
                i+=1.
            if Vmag!=None:
                mag+=Vmag
                i+=1.
            if Nmag!=None:
                mag+=Nmag
                i+=1.
            if i==0.:
                mag=None
            else:
                mag/=i

            if (ra!=None) & (dec!=None) & (mag!=None) & (clas!=None):
                if (mag<mag_limit) & (clas==0):
                    scaling = 1. - scale_factor * (mag - mag_pivot)
                    pos = self._wcs.all_world2pix(ra,dec,0)
                    stars_used[0].append(pos[0])
                    stars_used[1].append(pos[1])
                    stars_used[2].append(scaling)

        for i in range(len(stars_used[0])):
            poly = 'polygon('
            for x,y in zip(mask_model[0],mask_model[1]):
                angle = np.arctan2(y,x)
                l = stars_used[2][i] * np.sqrt(x**2. + y**2.)
                xnew = l * np.cos(angle)
                ynew = l * np.sin(angle)

                poly = poly + str(stars_used[0][i] + xnew + 0.5) + ' ' + str(stars_used[1][i] + ynew + 0.5) + ' '
            poly = poly + ')\n'
            mask_reg.write(poly)

        mask_reg.close()


    def _exec_WW(self,types=None):
        """Execute WeightWatcher

        Execute WeightWatcher to transform '.reg' to '.fits' flag map.

        Parameters
        ----------
        types : str
            Type of mask to make in ['HALO','SPIKE']

        """

        if types in ['HALO','SPIKE']:
            default_reg = self._config['PATH']['temp_dir'] +  types.lower()+ self._img_number + '.reg'
            defaul_out = self._config['PATH']['temp_dir'] +  types.lower() + '_flag' + self._img_number + '.fits'
            if self._config[types]['reg_file'] is None:
                reg=default_reg
                if not sc.BaseCatalog(reg)._helper.file_exists(reg):
                    raise sc.BaseCatalog.CatalogNotFound(reg)
                os.system('{0} -c {1} -WEIGHT_NAMES {2} -POLY_NAMES {3} -POLY_OUTFLAGS {4} -FLAG_NAMES "" -OUTFLAG_NAME {5} -OUTWEIGHT_NAME ""'.format(self._config['PATH']['WW'],self._config['PATH']['WW_configfile'],self._weight_fullpath,reg,self._config[types]['flag'],defaul_out))
                os.system('rm {0}'.format(reg))
            else:
                reg=self._config[types]['reg_file']
                if not sc.BaseCatalog(reg)._helper.file_exists(reg):
                    raise sc.BaseCatalog.CatalogNotFound(reg)
                os.system('{0} -c {1} -WEIGHT_NAMES {2} -POLY_NAMES {3} -POLY_OUTFLAGS {4} -FLAG_NAMES "" -OUTFLAG_NAME {5} -OUTWEIGHT_NAME ""'.format(self._config['PATH']['WW'],self._config['PATH']['WW_configfile'],self._weight_fullpath,reg,self._config[types]['flag'],defaul_out))

        elif types == 'ALL':
            default_reg = [self._config['PATH']['temp_dir'] + 'halo' + self._img_number + '.reg', self._config['PATH']['temp_dir'] + 'spike' + self._img_number + '.reg']
            defaul_out = self._config['PATH']['temp_dir'] + 'halo_spike_flag' + self._img_number + '.fits'
            if self._config['HALO']['reg_file'] is None:
                reg=default_reg
                for i in range(2):
                    if not sc.BaseCatalog(reg[i])._helper.file_exists(reg[i]):
                        raise sc.BaseCatalog.CatalogNotFound(reg[i])
                os.system('{0} -c {1} -WEIGHT_NAMES {2} -POLY_NAMES {3},{4} -POLY_OUTFLAGS {5},{6} -FLAG_NAMES "" -OUTFLAG_NAME {7} -OUTWEIGHT_NAME ""'.format(self._config['PATH']['WW'],self._config['PATH']['WW_configfile'],self._weight_fullpath,reg[0],reg[1],self._config['HALO']['flag'],self._config['SPIKE']['flag'],defaul_out))
                os.system('rm {0} {1}'.format(reg[0],reg[1]))
            else:
                reg=[self._config['HALO']['reg_file'], self._config['SPIKE']['reg_file']]
                for i in range(2):
                    if not sc.BaseCatalog(reg[i])._helper.file_exists(reg[i]):
                        raise sc.BaseCatalog.CatalogNotFound(reg[i])
                os.system('{0} -c {1} -WEIGHT_NAMES {2} -POLY_NAMES {3},{4} -POLY_OUTFLAGS {5},{6} -FLAG_NAMES "" -OUTFLAG_NAME {7} -OUTWEIGHT_NAME ""'.format(self._config['PATH']['WW'],self._config['PATH']['WW_configfile'],self._weight_fullpath,reg[0],reg[1],self._config['HALO']['flag'],self._config['SPIKE']['flag'],defaul_out))
        else:
                ValueError("types must be in ['HALO','SPIKE','ALL']")


    def _build_final_mask(self, path_mask1, path_mask2=None, border=None, messier=None):
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

        Returns
        -------
        numpy.ndarray
            Array containing the final mask.

        """

        final_mask=None

        if (path_mask1 is None) & (path_mask2 is None) & (border is None) & (messier is None):
            raise ValueError('No path to a mask, border or messier provided')

        if path_mask1 is not None:
            mask1=sc.FITSCatalog(path_mask1, hdu_no=0)
            mask1.open()
            final_mask=mask1.get_data()[:,:]

        if path_mask2 is not None:
            mask2=sc.FITSCatalog(path_mask2, hdu_no=0)
            mask2.open()
            if final_mask is not None:
                final_mask+=mask2.get_data()[:,:]
            else:
                final_mask=mask2.get_data()[:,:]

        if border is not None:
            if type(border) is np.ndarray:
                if final_mask is not None:
                    final_mask+=border
                else:
                    final_mask=border
            else:
                raise TypeError('border has to be a numpy.ndarray')

        if messier is not None:
            if type(messier) is np.ndarray:
                if final_mask is not None:
                    final_mask+=messier
                else:
                    final_mask=messier
            else:
                raise TypeError('messier has to be a numpy.ndarray')

        return final_mask.astype(np.int16,copy=False)


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

        out=sc.FITSCatalog(output_fullpath, open_mode= sc.BaseCatalog.OpenMode.ReadWrite, hdu_no= 0)
        out.save_as_fits(data=input_mask,image=True)

        if self._config['MD']['make']:
            out.open()
            out.add_header_card('MRATIO', self._ratio, 'ratio missing_pixels/all_pixels')
            out.add_header_card('MFLAG', self._config['MD']['im_flagged'], 'threshold value {:.3}'.format(self._config['MD']['thresh_flag']))
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

        path = temp_dir_path.replace(' ','')

        if path == 'OUTPUT':
            path = self._output_dir + '/temp'

        path += '/'
        if not os.path.isdir(path):
            os.system('mkdir ' + path)

        return path
