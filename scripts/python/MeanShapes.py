import numpy as np
import matplotlib.pyplot as plt
import os
from astropy.io import fits
import sys

import galsim

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

def MeanWhiskerPlot(ccd_maps_e1, ccd_maps_e2, filename, title='', wind=None):
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


def main():
    """ Compute and plot average shapes (e1, e2, R^2) on both stars and PSF model, and
    plot them on the MegaCam mosaic. Syntax:

    > python MeanShapes.py gridsize_x gridsize_y path/to/starcat
    Where gridsize_x and gridsize_y determine the number of bins per CCD in the horizontal
    and vertical directions, and path/to/starcat is the path to the full validation star
    catalog.
    """
    nb_pixel = int(sys.argv[1]),int(sys.argv[2])
    starcat_path = sys.argv[3]
    # MegaCam: each CCD is 2048x4612
    grid = np.linspace(0, 2048, nb_pixel[0]+1), np.linspace(0, 4612, nb_pixel[1]+1)

    # READ FULL STARCAT
    starcat = fits.open(starcat_path)[2].data

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

    # e_1
    vmax = max(np.nanmax(ccd_maps[:,:,0]), np.abs(np.nanmin(ccd_maps[:,:,0])))
    vmin = -vmax
    wind = [vmin, vmax]
    colorbar_ampl = 1
    MeanShapesPlot(ccd_maps[:,0,0], 'e1s', r'$e_1$ (stars)', wind=wind)
    MeanShapesPlot(ccd_maps[:,1,0], 'e1m', r'$e_1$ (model)', wind=wind)
    MeanShapesPlot(ccd_maps[:,0,0]-ccd_maps[:,1,0], 'e1res', r'$e_1$ residual', wind=wind,
            colorbar_ampl=colorbar_ampl)

    # e_2
    vmax = max(np.nanmax(ccd_maps[:,:,1]), np.abs(np.nanmin(ccd_maps[:,:,0])))
    colorbar_ampl = 1
    vmin = -vmax
    wind = [vmin, vmax]
    MeanShapesPlot(ccd_maps[:,0,1], 'e2s', r'$e_2$ (stars)', wind=wind)
    MeanShapesPlot(ccd_maps[:,1,1], 'e2m', r'$e_2$ (model)', wind=wind)
    MeanShapesPlot(ccd_maps[:,0,1]-ccd_maps[:,1,1], 'e2res', r'$e_2$ residual', wind=wind,
            colorbar_ampl=colorbar_ampl)

    # Whisker
    MeanWhiskerPlot(ccd_maps[:,0,0], ccd_maps[:,0,1], 'whisker_star', 'Whisker plot star', wind=wind)
    MeanWhiskerPlot(ccd_maps[:,1,0], ccd_maps[:,1,1], 'whisker_model', 'Whisker plot model', wind=wind)
    MeanWhiskerPlot(ccd_maps[:,0,0]-ccd_maps[:,1,0], ccd_maps[:,0,1]-ccd_maps[:,1,1], 'whisker_resi', 'Whisker plot resi', wind=wind)

    # R^2
    wind = [0,np.nanmax(ccd_maps[:,:,2])]
    colorbar_ampl = 1
    MeanShapesPlot(ccd_maps[:,0,2], 'R2s', r'$R^2$ (stars)', wind=wind, cmap='Reds')
    MeanShapesPlot(ccd_maps[:,1,2], 'R2m', r'$R^2$ (model)', wind=wind, cmap='Reds')
    wind = [-np.nanmax(ccd_maps[:,:,2]), np.nanmax(ccd_maps[:,:,2])]
    MeanShapesPlot(ccd_maps[:,0,2]-ccd_maps[:,1,2], 'R2res', r'$R^2$ residual', wind=wind,
            colorbar_ampl=colorbar_ampl)


if __name__ == "__main__":
    main()
