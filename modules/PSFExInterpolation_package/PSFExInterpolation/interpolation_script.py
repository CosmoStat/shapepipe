import numpy as np
import psfex
import scatalog as sc
import re
from astropy.io import fits

class PSFExInterpolator(object):
    def __init__(self, dotpsf_path, galcat_path, output_path, pos_params=''):
        self._dotpsf_path = dotpsf_path # Path to PSFEx output file
        self._galcat_path = galcat_path # Path to catalog containing galaxy positions
        self._output_path = output_path+'galaxy_psf'   # Path to output file to be written
        if pos_params:
            if not len(pos_params)==2:
                raise ValueError('{} position parameters were passed on; there should be exactly two.'.format(len(pos_params)))
            self._pos_params = pos_params
        else:
            self._pos_params = None
        self.gal_pos = None
        self.interp_PSFs = None
        
        # get number naming convention for this particular run
        s=re.split("\-([0-9]{3})\-([0-9]+)\.",self._galcat_path)
        self._img_number='-{0}-{1}'.format(s[1],s[2])      
        
    def _get_position_parameters(self):
        dotpsf = sc.FITSCatalog(self._dotpsf_path)
        dotpsf.open()
        self._pos_params = [dotpsf.get_header()['POLNAME1'], dotpsf.get_header()['POLNAME2']] 
        dotpsf.close()
        
    def _get_galaxy_positions(self):
        if self._pos_params is None:
            self._get_position_parameters()
        
        galcat = sc.FITSCatalog(self._galcat_path, SEx_catalog=True)
        galcat.open()
        try:
            self.gal_pos = np.array([[x,y] for x,y in zip(galcat.get_data()[self._pos_params[0]],
                                     galcat.get_data()[self._pos_params[1]])])
        except KeyError as detail:
            # extract erroneous position parameter from original exception
            err_pos_param = detail.args[0][4:-15]
            pos_param_err = 'Required position parameter '+err_pos_param+\
            'was not found in galaxy catalog. Leave pos_params (or EXTRA_CODE_OPTION) blank to read them from .psf file.'
            raise KeyError(pos_param_err)
        galcat.close()
    
    def _interpolate(self):
        if self.gal_pos is None:
            self._get_galaxy_positions()
        
        pex = psfex.PSFEx(self._dotpsf_path)
        self.interp_PSFs = np.array([pex.get_rec(x,y) for x,y in zip(self.gal_pos[:,0],
                                     self.gal_pos[:,1])])
        
    def write_output(self):
        if self.interp_PSFs is None:
            self._interpolate()
        output = fits.ImageHDU(self.interp_PSFs)
        output.writeto(self._output_path+self._img_number+'.fits', overwrite=True)
