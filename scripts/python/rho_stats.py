import numpy as np
import stile
from astropy.io import fits
import time
import matplotlib.pyplot as plt
import DES_rhostats as Drs
import sys

def NegDash(x_in, y_in, yerr_in, plot_name='', vertical_lines=True, 
            xlabel='', ylabel='',
            ylim=None, semilogx=False, semilogy=False,
            **kwargs):
    """ This function is for making plots with vertical errorbars, where negative values
    are shown in absolute value as dashed lines. The resulting plot can either be saved
    by specifying a file name as `plot_name', or be kept as a pyplot instance (for instance
    to combine several NegDashes).
    """
    x = np.copy(x_in)
    y = np.copy(y_in)
    yerr = np.copy(yerr_in)
    # catch and separate errorbar-specific keywords from Lines2D ones
    safekwargs = dict(kwargs)
    errbkwargs = dict()
    if 'linestyle' in kwargs.keys():
        print("""Warning: linestyle was provided but that would kind of defeat the purpose, 
         so I'll just ignore it. Sorry.""")
        del safekwargs['linestyle']
    for errorbar_kword in ['fmt', 'ecolor', 'elinewidth', 'capsize', 'barsabove', 'errorevery']:
        if errorbar_kword in kwargs.keys():
            #posfmt = '-'+kwargs['fmt']
            #negfmt = '--'+kwargs['fmt']
            errbkwargs[errorbar_kword] = kwargs[errorbar_kword]
            del safekwargs[errorbar_kword]
    errbkwargs = dict(errbkwargs, **safekwargs)
    """else:
        posfmt = '-'
        negfmt = '--'
        safekwargs = kwargs"""
    # plot up to next change of sign
    current_sign = np.sign(y[0])
    first_change = np.argmax(current_sign*y < 0)
    while first_change:
        if current_sign>0:
            plt.errorbar(x[:first_change], y[:first_change], yerr=yerr[:first_change], 
                         linestyle='-', **errbkwargs)
            if vertical_lines:
                plt.vlines(x[first_change-1], 0, y[first_change-1], linestyle='-', **safekwargs)
                plt.vlines(x[first_change], 0, np.abs(y[first_change]), linestyle='--', **safekwargs)
        else:
            plt.errorbar(x[:first_change], np.abs(y[:first_change]), yerr=yerr[:first_change], 
                         linestyle='--', **errbkwargs)
            if vertical_lines:
                plt.vlines(x[first_change-1], 0, np.abs(y[first_change-1]), linestyle='--', 
                            **safekwargs)
                plt.vlines(x[first_change], 0, y[first_change], linestyle='-', **safekwargs)
        x = x[first_change:]
        y = y[first_change:]
        yerr = yerr[first_change:]
        current_sign *= -1
        first_change = np.argmax(current_sign*y < 0)
    # one last time when `first_change'==0 ie no more changes:
    if current_sign>0:
        plt.errorbar(x, y, yerr=yerr, linestyle='-', **errbkwargs)
    else:
        plt.errorbar(x, np.abs(y), yerr=yerr, linestyle='--', **errbkwargs)
    if semilogx:
        plt.xscale('log')
    if semilogy:
        plt.yscale('log')
    if ylim is not None:
        plt.ylim(ylim)
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    if plot_name:
        plt.savefig(plot_name)
        plt.close()


def main():
    """ Compute and plot the five rho statistics. Syntax:
    
    > python rho_stats.py path/to/starcat rhofile_name [rho_def]
    Where path/to/starcat is the path to the full validation star
    catalog. rhofile_name is the name of the file where all rho_stats 
    will be written.
    rho_def is optional and can be either 'HSC' or 'DES'
    depending on the desired definition to use for tho statistics.
    """
    starcat_path = sys.argv[1] 
    rhofile_name = sys.argv[2]
    if len(sys.argv)>3:
        rho_def = sys.argv[3] 
    else:
        rho_def = 'HSC'

    # Read starcat
    starcat = fits.open(starcat_path)

    # Convert HSM flags to 0/1 weights
    star_flags = starcat[2].data['FLAG_STAR_HSM']
    psf_flags = starcat[2].data['FLAG_PSF_HSM']
    w = np.abs(star_flags-1) * np.abs(psf_flags-1)

    # Convert to Stile-compatible and change sigmas to R^2 (up to constant)
    stilecat = np.rec.fromarrays([w, starcat[2].data['RA'], starcat[2].data['DEC'],
                                  starcat[2].data['E1_STAR_HSM'], starcat[2].data['E2_STAR_HSM'], 
                                  starcat[2].data['SIGMA_STAR_HSM']**2,
                                  starcat[2].data['E1_PSF_HSM'], starcat[2].data['E2_PSF_HSM'], 
                                  starcat[2].data['SIGMA_PSF_HSM']**2],
                                  names=['w','ra','dec','g1','g2','sigma','psf_g1','psf_g2','psf_sigma'])

    # TreeCorr config:
    TreeCorrConfig = {'ra_units': 'degrees', 'dec_units': 'degrees', 
                      'max_sep': '3e2', 'min_sep': 5e-1, 'sep_units': 'arcminutes',
                      'nbins': 32}

    # Ininitialize all 5 rho stats
    if rho_def == 'HSC':
        rho_stats = [stile.CorrelationFunctionSysTest('Rho{}'.format(j)) for j in range(1,6)]
    elif rho_def == 'DES':
        rho_stats = [stile.CorrelationFunctionSysTest('Rho1'), Drs.DESRho2SysTest(), 
                      Drs.DESRho3SysTest(), Drs.DESRho4SysTest(), Drs.DESRho5SysTest()]
    for rho in rho_stats: 
        print rho.required_quantities

    # Compute them!
    print ' > Computing rho statistics...'
    start = time.time()
    rho_results = [rho_stat(stilecat, config=TreeCorrConfig) for rho_stat in rho_stats]
    print ' > Done in {}s.'.format(time.time()-start)
    np.save(rhofile_name+'.npy', np.array(rho_results))

    # in brackets are DES Y1's ylims
    ylim_l = None#(5e-9,5e-6)
    ylim_r = None#(5e-8,1e-5)

    # Plots
    ylims = [ylim_l, ylim_r, ylim_l, ylim_l, ylim_r]
    colors = ['blue', 'blue', 'green', 'orange', 'green']
    markers =['o', 'o', 's', '^', 's']

    for j,rhores in enumerate(rho_results):
        NegDash(rhores['meanR [arcmin]'], rhores['xip'], rhores['sigma_xi'], 
                'rho_{}.png'.format(j+1), 
                semilogx=True, semilogy=True,
                color=colors[j], capsize=3, fmt=markers[j], alpha=.7, ylim=ylims[j],
                xlabel=r'$\theta$ (arcmin)', ylabel=r'$\rho_i(\theta)$')


    for j,rhores in enumerate(rho_results):
        if j in [0,2,3]:
            NegDash(rhores['meanR [arcmin]'], rhores['xip'], rhores['sigma_xi'], 
                    semilogx=True, semilogy=True,
                    color=colors[j], capsize=3, fmt=markers[j], alpha=.7, ylim=ylims[j],
                    xlabel=r'$\theta$ (arcmin)', ylabel=r'$\rho_i(\theta)$')
    plt.savefig('lefthand_rhos.pdf')
    plt.close()

    for j,rhores in enumerate(rho_results):
        if j in [1,4]:
            NegDash(rhores['meanR [arcmin]'], rhores['xip'], rhores['sigma_xi'], 
                    semilogx=True, semilogy=True,
                    color=colors[j], capsize=3, fmt=markers[j], alpha=.7, ylim=ylims[j],
                    xlabel=r'$\theta$ (arcmin)', ylabel=r'$\rho_i(\theta)$')
    plt.savefig('righthand_rhos.pdf')
    plt.close()


if __name__ == "__main__":
    main()
