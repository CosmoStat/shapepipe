import numpy as np
import psfex
import scatalog as sc

class PSFExInterpolator(object):
    def __init__(self, dotpsf_path, galcat_path, output_path):
        self._dotpsf_path = dotpsf_path # Path to PSFEx output file
        self._galcat_path = galcat_path # Path to catalog containing galaxy positions
        self._output_path = output_path   # Path to output file to be written
        self._pos_params = None
        self.gal_pos = None
        self.interp_PSFs = None
        
    def _get_position_parameters(self):
        self._pos_params = ['XWIN_IMAGE', 'YWIN_IMAGE'] #TEMP: hardcoded pos params
        
    def _get_galaxy_positions(self):
        if self._pos_params is None:
            self._get_position_parameters()
        
        galcat = sc.FITSCatalog(self._galcat_path, SEx_catalog=True)
        galcat.open()
        gal_data = galcat.get_data()
        self.gal_pos = np.array([[x,y] for x,y in zip(gal_data[self._pos_params[0]],
                                 gal_data[self._pos_params[1]])])
        galcat.close()
    
    def _interpolate(self):
        if self.gal_pos is None:
            self._get_galaxy_positions()
        
        pex = psfex.PSFEx(self._dotpsf_path)
        self.interp_PSFs = np.array([pex.get_rec(x,y) for x,y in zip(self.gal_pos[:,0],
                                     self.gal_pos[:,1])])
        
    def _write_output(self):
        if self.interp_PSFs is None:
            self._interpolate()
        
        output = sc.FITSCatalog(self._output_path, 
                                open_mode=sc.BaseCatalog.OpenMode.ReadWrite)
        output.save_as_fits(self.interp_PSFs, self.gal_pos)