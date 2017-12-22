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
    """!
        Mask creation module
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

        self._set_parameters()                                                  # Set parameters needed for the stars detection


    def _get_config(self, config_filepath=None):
        """!
            Read config file and set parameters
            @param config_filepath path to 'config.mask'
        """

        if config_filepath is None:
            raise ValueError('No path to config file')

        conf=sconfig.SConfig(config_filepath)

        self._config={'PATH': {}, 'BORDER': {}, 'HALO': {}, 'SPIKE': {}}

        self._config['PATH']['WW'] = conf.get_as_string('WW_PATH','PROGRAM_PATH')
        self._config['PATH']['WW_configfile'] = conf.get_as_string('WW_CONFIG_FILE','PROGRAM_PATH')
        self._config['PATH']['CDSclient'] = conf.get_as_string('CDSCLIENT_PATH','PROGRAM_PATH')
        self._config['PATH']['temp_dir'] = conf.get_as_string('TEMP_DIRECTORY','OTHER')
        self._config['BORDER']['make'] = conf.get_as_boolean('BORDER_MAKE','BORDER_PARAMETERS')
        if self._config['BORDER']['make'] is True:
            self._config['BORDER']['width'] = conf.get_as_int('BORDER_WIDTH','BORDER_PARAMETERS')
            self._config['BORDER']['flag'] = conf.get_as_string('BORDER_FLAG_VALUE','BORDER_PARAMETERS')
        for i in ['HALO','SPIKE']:
            self._config[i]['make'] = conf.get_as_boolean(i+'_MAKE',i+'_PARAMETERS')
            self._config[i]['individual'] = conf.get_as_boolean('KEEP_INDIVIDUAL_MASK','OTHER')
            if self._config[i]['make'] == True:
                self._config[i]['maskmodel_path'] = conf.get_as_string(i+'_MASKMODEL_PATH',i+'_PARAMETERS')
                self._config[i]['mag_lim'] = conf.get_as_float(i+'_MAG_LIM',i+'_PARAMETERS')
                self._config[i]['scale_factor'] = conf.get_as_float(i+'_SCALE_FACTOR',i+'_PARAMETERS')
                self._config[i]['flag'] = conf.get_as_int(i+'_FLAG_VALUE',i+'_PARAMETERS')
                if conf.get_as_boolean('KEEP_REG_FILE','OTHER') == True:
                    self._config[i]['reg_file'] = self._config['PATH']['temp_dir'] + '{0}{1}.reg'.format(re.split(".reg",conf.get_as_string(i+'_REG_FILE',i+'_PARAMETERS'))[0],self._img_number)
                else:
                    self._config[i]['reg_file'] = None


    def _set_parameters(self):
        """!
            Set the parameters for the stars detection
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
        """!
            Main function to create the mask
        """

        stars=self.find_stars(np.array([self._fieldcenter['wcs'].ra.value,self._fieldcenter['wcs'].dec.value]), radius=self._img_radius)

        for i in ['HALO', 'SPIKE']:
            if self._config[i]['make']:
                self._create_mask(stars=stars, types=i, mag_limit=self._config[i]['mag_lim'], scale_factor=self._config[i]['scale_factor'])

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

        final_mask=self._build_final_mask(path_mask1=mask_name[0],path_mask2=mask_name[1], border=border_mask)

        if not self._config['HALO']['individual']:
            if mask_name[0] is not None:
                os.system('rm {0}'.format(mask_name[0]))
            if mask_name[1] is not None:
                os.system('rm {1}'.format(mask_name[1]))

        self._mask_to_file(input_mask=final_mask, output_fullpath=self._output_dir+'/'+self._img_name+'_flag'+self._img_number+'.fits')


    def find_stars(self, position, radius=None):
        """!
            Return GSC (Guide Star Catalog) objects for a field with center (ra,dec) and radius r
            @param ra right ascention astropy.wcs oject
            @param dec declinaison astropy.wcs oject
            @param r radius in arcmin
            @return stars dicotionnary for GSC objects in the field
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
        """!
            Mask 'width' pixels around the image
            @param width width of the mask mask border
            @return array containing the mask
        """

        if width is None:
            raise ValueError('Width not provided')

        flag = np.zeros((int(self._fieldcenter['pix'][0]*2),int(self._fieldcenter['pix'][1]*2)),dtype='uint8')

        flag[0:width,:]=flag_value
        flag[-width:-1,:]=flag_value
        flag[:,0:width]=flag_value
        flag[:,-width:-1]=flag_value

        return flag


    def SphereDist(self, position1, position2):
        """!
            Compute spheric distance between 2 points
            @param p1 array [x,y] first point (in pixel)
            @param p2 array [x,y] second point (in pixel)
            @param wcs astropy.wcs object containing wcs image informations
            @return the distance in degree
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
        """!
            Compute the diagonal distance of the image in arcmin.
            @param center coordinates of the center of the image in pixel
            @return the diagonal distance in arcmin
        """

        if center is None:
            return self.SphereDist(self._fieldcenter['pix'], np.zeros(2))/60.
        else:
            if type(center) is np.ndarray:
                return self.SphereDist(center, np.zeros(2))/60.
            else:
                raise TypeError('center has to be a numpy.ndarray')


    def _make_star_cat(self, CDSclient_output):
        """!
            Make a dicotionnary from 'findgsc2.2' output
            @param CDSclient_output output 'findgsc2.2'
            @return stars dicotionnary containing all informations
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


    def _create_mask(self, stars=None, types=None, mag_limit=18., mag_pivot=13.8, scale_factor=0.3):
        """!
            Apply mask from model to stars and save into DS9 region file
            @param stars stars dico (output of find_stars)
            @param types type of mask in ['halo', 'spike']
            @param mag_limit higher magnitude to apply the mask
            @param mag_pivot pivot magnitude for the model
            @param scale_factor scaling for the model
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
        """!
            Execute WW to transform '.reg' to '.fits' flag map
            @param types the type of mask to make in ['HALO','SPIKE']
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


    def _build_final_mask(self, path_mask1=None, path_mask2=None, border=None):
        """!
            Create the final mask by combination of individual mask
            @param path to a mask (fits format)
            @param path to a mask (fits format)
            @param border array containing the border mask
            @return final_mask array containing the final mask
        """

        final_mask=None

        if (path_mask1 is None) & (path_mask2 is None) & (border is None):
            raise ValueError('No path to a mask or border provided')

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

        return final_mask.astype(np.int16,copy=False)


    def _mask_to_file(self, input_mask=None, output_fullpath=None):
        """!
            Save the mask to a fits file
            @param input_mask mask to save
            @param output_fullpath path of the output file
        """

        if input_mask is None:
            raise ValueError('input_mask not provided')
        if output_fullpath is None:
            raise ValueError('fullpath not provided')

        out=sc.FITSCatalog(output_fullpath, open_mode= sc.BaseCatalog.OpenMode.ReadWrite)
        out.save_as_fits(data=input_mask,image=True)
