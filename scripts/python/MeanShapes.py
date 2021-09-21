#!/usr/bin/env python

"""Script MeanShapes.py

Unzips and removes first (empty) HDU of CFIS tile weight
such that they can be read by SExtractor,

:Authors: Axel Guinot, Morgan Schmitz, Martin Kilbinger

:Date: 2019, 2020

:Package: ShapePipe
"""

# Compability with python2.x for x>6
from __future__ import print_function

import os
import sys
import copy

import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
from optparse import OptionParser

import galsim


class param:
    """General class to store (default) variables
    """

    def __init__(self, **kwds):
        self.__dict__.update(kwds)

    def print(self, **kwds):
        print(self.__dict__)

    def var_list(self, **kwds):
        return vars(self)


def params_default():
    """Set default parameter values.

    Parameters
    ----------
    None

    Returns
    -------
    p_def: class tuff.param
        parameter values
    """

    p_def = param(
        nx = 20,
        input_path = './psf_cat_full.fits',
        output_dir = './psf_validation',
        hdu = 2,
    )

    return p_def


def parse_options(p_def):
    """Parse command line options.

    Parameters
    ----------
    p_def: class param
        parameter values

    Returns
    -------
    options: tuple
        Command line options
    args: string
        Command line string
    """

    usage  = "%prog [OPTIONS]"
    parser = OptionParser(usage=usage)

    parser.add_option('-i', '--input', dest='input_path', type='string', default=p_def.input_path,
         help='input file name, default=\'{}\''.format(p_def.input_path))
    parser.add_option('-o', '--output_dir', dest='output_dir', type='string', default=p_def.output_dir,
         help='output_directory, default=\'{}\''.format(p_def.output_dir))

    parser.add_option('-x', '--npix_x', dest='nx', type='int', default=p_def.nx,
         help='number of pixels per CCD in the x-direction, default=\'{}\''.format(p_def.nx))
    parser.add_option('-y', '--npix_y', dest='ny', type='int',
         help='number of pixels per CCD in the y-direction, default to provide square pixels given nx')
    parser.add_option('', '--max_e', dest='max_e', type='float',
         help='max value for ellipticity plots (model, star)')
    parser.add_option('', '--max_d', dest='max_d', type='float',
         help='max value for ellipticity residuals plots (model, star)')

    parser.add_option('', '--hdu', dest='hdu', type='int', default=p_def.hdu,
         help='HDU number on input, default={}'.format(p_def.hdu))

    parser.add_option('-v', '--verbose', dest='verbose', action='store_true', help='verbose output')

    options, args = parser.parse_args()

    return options, args


def check_options(options):
    """Check command line options.

    Parameters
    ----------
    options: tuple
        Command line options

    Returns
    -------
    erg: bool
        Result of option check. False if invalid option value.
    """

    if not os.path.isdir(options.output_dir):
        print('Output directory \'{}\' does not exist'.format(options.output_dir))
        return False

    return True


def update_param(p_def, options):
    """Return default parameter, updated and complemented according to options.

    Parameters
    ----------
    p_def:  class param
        parameter values
    optiosn: tuple
        command line options

    Returns
    -------
    param: class param
        updated paramter values
    """

    param = copy.copy(p_def)

    # Update keys in param according to options values
    for key in vars(param):
        if key in vars(options):
            setattr(param, key, getattr(options, key))

    # Add remaining keys from options to param
    for key in vars(options):
        if not key in vars(param):
            setattr(param, key, getattr(options, key))

    if not param.ny:
        npix_x, npix_y = MegaCamNpix()
        param.ny = int(param.nx * npix_y / npix_x)

    return param


def MegaCamNpix():
    """Return number of pixels of a MegaCam CCD.

    Parameters
    ----------
    None

    Returns
    -------
    nx, ny: int
        number of pixels along x, y
    """

    return 2048, 4612


# MegaCam -> plt.subplot correspondance, as given by:
'''        'COMMENT Unique detector IDs for MegaCam (North on top, East to the left)',
           'COMMENT    --------------------------',
           'COMMENT    ba ba ba ba ba ba ba ba ba',
           'COMMENT    00 01 02 03 04 05 06 07 08',
           'COMMENT --------------------------------',
           'COMMENT ba ba ba ba ba ba ba ba ba ba ba',
           'COMMENT 36 09 10 11 12 13 14 15 16 17 37',
           'COMMENT --------------------------------',
           'COMMENT 38 18 19 20 21 22 23 24 25 26 39',
           'COMMENT ab ab ab ab ab ab ab ab ab ab ab',
           'COMMENT --------------------------------',
           'COMMENT    27 28 29 30 31 32 33 34 35',
           'COMMENT    ab ab ab ab ab ab ab ab ab',
           'COMMENT    __________________________'
'''
def MegaCamPos(j):
    if j < 9:
        # first row - shift by one
        return j+1
    elif j < 18:
        # second row, non-ears
        return j+3
    elif j < 27:
        # third row non-ears
        return j+5
    elif j < 36:
        # fourth row
        return j+7
    else:
        MegaCamCoords = {36:11, 37:21, 38:22, 39:32}
        return MegaCamCoords[j]

def MegaCamFlip(xbins, ybins, ccd_nb, nb_pixel):

    if ccd_nb < 18 or ccd_nb in [36,37]:
        # swap x axis so origin is on top-right
        xbins = nb_pixel[0] - xbins + 1
    else:
        # swap y axis so origin is on bottom-left
        ybins = nb_pixel[1] - ybins + 1

    return xbins, ybins


def MeanShapesPlot(ccd_maps, filename, title='', colorbar_ampl=1., wind=None, cmap='bwr'):

    # colorbar amplitude
    if wind is None:
        vmax = max(np.nanmax(ccd_maps), np.abs(np.nanmin(ccd_maps))) * colorbar_ampl
        vmin = -vmax * colorbar_ampl
    else:
        vmin, vmax = wind[0]*colorbar_ampl, wind[1]*colorbar_ampl

    # create full plot
    fig, axes = plt.subplots(nrows=4,ncols=11)
    # remove corner axes (above and below ears)
    for j in [0, 10, -1, -11]:
        axes.flat[j].axis('off')
    for ccd_nb,ccd_map in enumerate(ccd_maps):
        ax = axes.flat[MegaCamPos(ccd_nb)]
        im = ax.imshow(ccd_map.T,cmap=cmap,interpolation='Nearest',vmin=vmin,vmax=vmax)
        ax.set_xticks([])
        ax.set_yticks([])
    plt.suptitle(title) #TODO: fix title
    fig.subplots_adjust(right=0.8)
    cbar_ax = fig.add_axes([0.85, 0.15, 0.05, 0.7])
    fig.colorbar(im, cax=cbar_ax)

    plt.savefig('{}.png'.format(filename))
    plt.close()


def MeanWhiskerPlot(ccd_maps_e1, ccd_maps_e2, filename, title=''):
    """
    """

    nx, ny = ccd_maps_e1[0].shape

    # create full plot
    fig, axes = plt.subplots(nrows=4,ncols=11)
    # remove corner axes (above and below ears)
    for j in [0, 10, -1, -11]:
        axes.flat[j].axis('off')

    n_ccd = len(ccd_maps_e1)
    for ccd_nb in range(n_ccd):
        ax = axes.flat[MegaCamPos(ccd_nb)]
        for x in range(nx):
            for y in range(ny):
                true_x = x + 0.5
                true_y = y + 0.5

                trushap = galsim.Shear(g1=ccd_maps_e1[ccd_nb][x,y],
                                       g2=ccd_maps_e2[ccd_nb][x,y])
                U = trushap.g * np.cos(trushap.beta)
                V = trushap.g * np.sin(trushap.beta)

                q = ax.quiver(true_x, true_y, U, V, headwidth=0, headlength=0, headaxislength=0, alpha=1.,
                           angles='xy',
                           scale_units='inches', scale=.07)
        # im = ax.imshow(ccd_map.T,cmap=cmap,interpolation='Nearest',vmin=vmin,vmax=vmax)
        ax.set_xticks([])
        ax.set_yticks([])
        # leg = ax.quiverkey(q, X=.94, Y=1.02, U=0.05,
                           # label=r'$e$=0.05', labelpos='N')
        ax.set_xlim(0, nx)
        ax.set_ylim(0, ny)
        # ax.set_xlabel(r'$x$-coordinate position (deg)')
        # ax.set_ylabel(r'$y$-coordinate position (deg)')
    plt.suptitle(title) #TODO: fix title
    fig.subplots_adjust(right=0.8)

    plt.savefig('{}.png'.format(filename))
    plt.close()


def log_command(argv, name=None, close_no_return=True):
    """Write command with arguments to a file or stdout.
       Choose name = 'sys.stdout' or 'sys.stderr' for output on sceen.

    Parameters
    ----------
    argv: array of strings
        Command line arguments
    name: string
        Output file name (default: 'log_<command>')
    close_no_return: bool
        If True (default), close log file. If False, keep log file open
        and return file handler

    Returns
    -------
    log: filehandler
        log file handler (if close_no_return is False)
    """

    if name is None:
        name = 'log_' + os.path.basename(argv[0])

    if name == 'sys.stdout':
        f = sys.stdout
    elif name == 'sys.stderr':
        f = sys.stderr
    else:
        f = open(name, 'w')

    for a in argv:

        # Quote argument if special characters
        if ']' in a or ']' in a:
            a = '\"{}\"'.format(a)

        print(a, end='', file=f)
        print(' ', end='', file=f)

    print('', file=f)

    if close_no_return == False:
        return f

    if name != 'sys.stdout' and name != 'sys.stderr':
        f.close()


def main(argv=None):
    """ Compute and plot average shapes (e1, e2, R^2) on both stars and PSF model, and
    plot them on the MegaCam mosaic. Syntax:

    > python MeanShapes.py gridsize_x gridsize_y path/to/starcat
    Where gridsize_x and gridsize_y determine the number of bins per CCD in the horizontal
    and vertical directions, and path/to/starcat is the path to the full validation star
    catalog.
    """

    # Set default parameters
    p_def = params_default()

    # Command line options
    options, args = parse_options(p_def)

    if check_options(options) is False:
        return 1

    param = update_param(p_def, options)

    # Save calling command
    log_command(argv)
    if param.verbose:
        log_command(argv, name='sys.stderr')


    nb_pixel = param.nx, param.ny
    starcat_path = param.input_path


    # MegaCam: each CCD is 2048x4612
    npix_x, npix_y = MegaCamNpix()
    grid = np.linspace(0, npix_x, nb_pixel[0]+1), np.linspace(0, npix_y, nb_pixel[1]+1)

    # Read full star catalogue
    starcat = fits.open(starcat_path)[param.hdu].data

    # Flag mask
    star_flags = starcat['FLAG_STAR_HSM']
    psf_flags = starcat['FLAG_PSF_HSM']
    flagmask = np.abs(star_flags-1) * np.abs(psf_flags-1)

    # convert sigma to R^2's
    all_star_shapes = np.array([starcat['E1_STAR_HSM'],starcat['E2_STAR_HSM'],
                            2.*starcat['SIGMA_STAR_HSM']**2])
    all_psf_shapes = np.array([starcat['E1_PSF_HSM'],starcat['E2_PSF_HSM'],
                            2.*starcat['SIGMA_PSF_HSM']**2])

    ccd_maps = np.ones((40, 2, 3)+nb_pixel) * np.nan # CCDs x star/model x (e1,e2,R2) x xpos x ypos
    for ccd_nb,ccd_map in enumerate(ccd_maps):
        ccd_mask = ((starcat['CCD_NB']==str(ccd_nb)) * flagmask).astype(bool)

        star_shapes = all_star_shapes[:,ccd_mask]
        psf_shapes = all_psf_shapes[:,ccd_mask]

        xs, ys = starcat['X'][ccd_mask], starcat['Y'][ccd_mask]
        xbins = np.digitize(xs, grid[0])
        ybins = np.digitize(ys, grid[1])
        # swap axes to match CCD orientation and origin convention
        xbins, ybins = MegaCamFlip(xbins, ybins, ccd_nb, nb_pixel)
        for xb in range(nb_pixel[0]):
            for yb in range(nb_pixel[1]):
                bin_star_shapes = star_shapes[:,(xbins==xb+1) * (ybins==yb+1)]
                bin_psf_shapes = psf_shapes[:,(xbins==xb+1) * (ybins==yb+1)]
                ccd_map[0,:,xb,yb] = np.mean(bin_star_shapes,axis=1)
                ccd_map[1,:,xb,yb] = np.mean(bin_psf_shapes,axis=1)

    # Create plots
    colorbar_ampl = 1

    # e_1
    if not param.max_e:
        vmax = max(np.nanmax(ccd_maps[:,:,0]), np.abs(np.nanmin(ccd_maps[:,:,0])))
    else:
        vmax = param.max_e
    vmin = -vmax
    wind_e = [vmin, vmax]
    MeanShapesPlot(ccd_maps[:,0,0], '{}/e1s'.format(param.output_dir), r'$e_1$ (stars)', wind=wind_e)
    MeanShapesPlot(ccd_maps[:,1,0], '{}/e1m'.format(param.output_dir), r'$e_1$ (model)', wind=wind_e)

    if param.max_d:
        vmax = param.max_d
        vmin = -vmax
        wind_d = [vmin, vmax]
    else:
        wind_d = wind_e
    MeanShapesPlot(ccd_maps[:,0,0]-ccd_maps[:,1,0], '{}/e1res'.format(param.output_dir), r'$e_1$ residual', wind=wind_d,
            colorbar_ampl=colorbar_ampl)

    # e_2
    if not param.max_e:
        vmax = max(np.nanmax(ccd_maps[:,:,1]), np.abs(np.nanmin(ccd_maps[:,:,0])))
    else:
        vmax = param.max_e
    vmin = -vmax
    wind_e = [vmin, vmax]
    MeanShapesPlot(ccd_maps[:,0,1], '{}/e2s'.format(param.output_dir), r'$e_2$ (stars)', wind=wind_e)
    MeanShapesPlot(ccd_maps[:,1,1], '{}/e2m'.format(param.output_dir), r'$e_2$ (model)', wind=wind_e)

    if param.max_d:
        vmax = param.max_d
        vmin = -vmax
        wind_d = [vmin, vmax]
    else:
        wind_d = wind_e
    MeanShapesPlot(ccd_maps[:,0,1]-ccd_maps[:,1,1], '{}/e2res'.format(param.output_dir), r'$e_2$ residual', wind=wind_d,
            colorbar_ampl=colorbar_ampl)

    # Whisker
    MeanWhiskerPlot(ccd_maps[:,0,0], ccd_maps[:,0,1], '{}/whisker_star'.format(param.output_dir), 'Whisker plot star')
    MeanWhiskerPlot(ccd_maps[:,1,0], ccd_maps[:,1,1], '{}/whisker_model'.format(param.output_dir), 'Whisker plot model')
    MeanWhiskerPlot(ccd_maps[:,0,0]-ccd_maps[:,1,0], ccd_maps[:,0,1]-ccd_maps[:,1,1],
                    '{}/whisker_resi'.format(param.output_dir), 'Whisker plot resi')

    # R^2
    wind = [0,np.nanmax(ccd_maps[:,:,2])]

    MeanShapesPlot(ccd_maps[:,0,2], '{}/R2s'.format(param.output_dir), r'$R^2$ (stars)', wind=wind, cmap='Reds')
    MeanShapesPlot(ccd_maps[:,1,2], '{}/R2m'.format(param.output_dir), r'$R^2$ (model)', wind=wind, cmap='Reds')
    wind = [-np.nanmax(ccd_maps[:,:,2]), np.nanmax(ccd_maps[:,:,2])]
    MeanShapesPlot(ccd_maps[:,0,2]-ccd_maps[:,1,2], '{}/R2res'.format(param.output_dir), r'$R^2$ residual', wind=wind,
            colorbar_ampl=colorbar_ampl)
    
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))
