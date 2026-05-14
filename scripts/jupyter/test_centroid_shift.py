#!/usr/bin/env python
# coding: utf-8

# In[14]:


import numpy as np
import galsim

import ngmix

from shapepipe.modules.ngmix_package.ngmix import (
    get_noise,
)

#print("Use shapepipe version of metacal")
#from shapepipe.modules.ngmix_package.ngmix import do_ngmix_metacal

print("Use this version of metacal")

from modopt.math.stats import sigma_mad
from ngmix import Observation, ObsList
from numpy.random import uniform as urand


# ### Simu

# In[15]:


def make_data(rng, shear, noise=1e-5, n_epochs=1, share_shift=False):
    """
    simulate an exponential object with moffat psf

    Parameters
    ----------
    rng: np.random.RandomState
        The random number generator
    noise: float
        Noise for the image
    shear: (g1, g2)
        The shear in each component

    Returns
    -------
    Lists
    """

    psf_noise = 1.0e-6

    img_size = 201
    scale = 0.1857
    wcs = galsim.PixelScale(scale)

    psf_fwhm = 0.55
    gal_hlr = 0.3

    if share_shift:
        dy, dx = rng.uniform(low=-scale/2, high=scale/2, size=2)


    gals = []
    psfs = []
    psfs_sigmas = []
    weights = []
    flags = []
    jacob_lists = []
    for epoch in range(n_epochs):
        if not share_shift:
            dy, dx = rng.uniform(low=-scale/2, high=scale/2, size=2)
        psf = galsim.Moffat(
            beta=2.5, fwhm=psf_fwhm,
        )

        obj0 = galsim.Exponential(
            half_light_radius=gal_hlr,
            flux=1000,
        ).shear(
            g1=shear[0],
            g2=shear[1],
        )

        obj = galsim.Convolve(psf, obj0)
        obj = obj.shift(dx, dy)

        psf_im_ = psf.drawImage(nx=img_size, ny=img_size, wcs=wcs)
        s = galsim.hsm.FindAdaptiveMom(psf_im_)
        psfs_sigmas.append(s.moments_sigma)
        psf_im = psf_im_.array.astype(np.float64)
        im = obj.drawImage(nx=img_size, ny=img_size, wcs=wcs).array.astype(np.float64)

        psf_im += rng.normal(scale=psf_noise, size=psf_im.shape)
        im += rng.normal(scale=noise, size=im.shape)


        wt = im*0 + 1.0/noise**2
        flag = im*0

        weights.append(wt)
        psfs.append(psf_im)
        flags.append(flag)
        gals.append(im)
        jacob_lists.append(wcs.jacobian())

    return gals, psfs, psfs_sigmas, weights, flags, jacob_lists


# ### From shapepipe

# In[16]:


def get_prior(pixel_scale, rng=None):
    """Get Prior.

    Return prior for the different parameters.

    Returns
    -------
    ngmix.priors
        Priors for the different parameters

    """
    if rng is None:
        rng = np.random.default_rng()

    # Prior on ellipticity. Details do not matter, as long
    # as it regularizes the fit. From Bernstein & Armstrong 2014
    g_sigma = 0.4
    g_prior = ngmix.priors.GPriorBA(g_sigma, rng=rng)

    # 2-d Gaussian prior on the center row and column center
    # (relative to the center of the jacobian, which
    # would be zero) and the sigma of the Gaussians.
    # Units same as jacobian, probably arcsec
    row, col = 0.0, 0.0
    row_sigma, col_sigma = pixel_scale, pixel_scale
    cen_prior = ngmix.priors.CenPrior(row, col, row_sigma, col_sigma, rng=rng)

    # Size prior. Instead of flat, two-sided error function (TwoSidedErf)
    # could be used
    Tminval = -10.0  # arcsec squared
    Tmaxval = 1.0e6
    T_prior = ngmix.priors.FlatPrior(Tminval, Tmaxval, rng=rng)

    # Flux prior. Bounds need to make sense for
    # images in question
    Fminval = -1.0e4
    Fmaxval = 1.0e9
    F_prior = ngmix.priors.FlatPrior(Fminval, Fmaxval, rng=rng)

    # Joint prior, combine all individual priors
    prior = ngmix.joint_prior.PriorSimpleSep(
        cen_prior, g_prior, T_prior, F_prior
    )

    return prior


# In[17]:


def make_galsimfit(obs, model, guess0, prior=None, ntry=5):
    """Make GalSim Fit.

    Fit image using simple GalSim model.

    Parameters
    ----------
    obs : ngmix.observation.Observation
        Image to fit
    model : str
        Model for fit
    guess0 : numpy.ndarray
        Parameters of first model guess
    prior : ngmix.prior, optional
        Prior for fit paraemeters
    ntry : int, optional
        Number of tries for fit, the default is ``5``

    Returns
    -------
    dict
        Results

    Raises
    ------
    ngmix.BootGalFailure
        Failure to bootstrap galaxy

    """
    limit = 0.1

    guess = np.copy(guess0)
    fres = {}
    for it in range(ntry):
        guess[0:5] += urand(low=-limit, high=limit)
        guess[5:] *= 1 + urand(low=-limit, high=limit)
        fres["flags"] = 1
        try:
            fitter = ngmix.fitting.GalsimFitter(model, prior=prior)
            fres = fitter.go(obs, guess)
            if fres["flags"] == 0 and "T" not in fres:
                fres["T"] = fres["pars"][4]
                fres["T_err"] = np.sqrt(fres["pars_cov"][4, 4])
        except:
            continue

        if fres["flags"] == 0:
            break

    if fres["flags"] != 0:
        raise ngmix.gexceptions.BootGalFailure(
            "Failed to fit galaxy image with galsimfit"
        )

    fres["ntry"] = it + 1

    return fres


# In[18]:


def do_ngmix_metacal(
    gals, psfs, psfs_sigma, weights, flags, jacob_list, prior, pixel_scale, sig_noise=None,
):
    """Do Ngmix Metacal.

    Perform the metacalibration on a multi-epoch object and return the joint
    shape measurement with NGMIX.

    Parameters
    ----------
    gals : list
        List of the galaxy vignets
    psfs : list
        List of the PSF vignets
    psfs_sigma : list
        List of the sigma PSFs
    weights : list
        List of the weight vignets
    flags : list
        List of the flag vignets
    jacob_list : list
        List of the Jacobians
    prior : ngmix.priors
        Priors for the fitting parameters
    pixel_scale : float
        pixel scale in arcsec

    Returns
    -------
    dict
        Dictionary containing the results of NGMIX metacal

    """
    n_epoch = len(gals)

    if n_epoch == 0:
        raise ValueError("0 epoch to process")

    # Make observation
    gal_obs_list = ObsList()
    T_guess_psf = []
    psf_res_gT = {
        "g_PSFo": np.array([0.0, 0.0]),
        "g_err_PSFo": np.array([0.0, 0.0]),
        "T_PSFo": 0.0,
        "T_err_PSFo": 0.0,
    }
    gal_guess = []
    gal_guess_flag = True
    wsum = 0
    for n_e in range(n_epoch):

        psf_jacob = ngmix.Jacobian(
            row=(psfs[0].shape[0] - 1) / 2,
            col=(psfs[0].shape[1] - 1) / 2,
            wcs=jacob_list[n_e],
        )

        psf_obs = Observation(psfs[n_e], jacobian=psf_jacob)

        psf_T = psfs_sigma[n_e] * 1.17741 * pixel_scale

        weight_map = np.copy(weights[n_e])
        weight_map[np.where(flags[n_e] != 0)] = 0.0
        weight_map[weight_map != 0] = 1

        psf_guess = np.array([0.0, 0.0, 0.0, 0.0, psf_T, 1.0])
        try:
            psf_res = make_galsimfit(psf_obs, "gauss", psf_guess)
        except Exception:
            continue

        # Gal guess via adaptive moments
        try:
            _gim = galsim.Image(gals[n_e], scale=pixel_scale)
            _hsm = galsim.hsm.FindAdaptiveMom(_gim, strict=False)
            if _hsm.error_message != "":
                raise galsim.hsm.GalSimHSMError(_hsm.error_message)
            _cen = _hsm.moments_centroid - _gim.center   # pixels (centroid_unit="img")
            _size = 2 * (_hsm.moments_sigma * pixel_scale) ** 2  # T in arcsec^2
            gal_guess_tmp = np.array([_cen.x, _cen.y, 0.0, 0.0, _size, _hsm.moments_amp])
        except Exception:
            gal_guess_flag = False
            gal_guess_tmp = np.array([0.0, 0.0, 0.0, 0.0, 1, 100])

        # Recenter jacobian if necessary
        gal_jacob = ngmix.Jacobian(
            row=(gals[n_e].shape[0] - 1) / 2 + gal_guess_tmp[1],
            col=(gals[n_e].shape[1] - 1) / 2 + gal_guess_tmp[0],
            wcs=jacob_list[n_e],
        )

        # Noise handling
        if sig_noise is None:
            if gal_guess_flag:
                sig_noise = get_noise(
                    gals[n_e],
                    weight_map,
                    gal_guess_tmp,
                    pixel_scale,
                )
            else:
                sig_noise = sigma_mad(gals[n_e])

        noise_img = np.random.randn(*gals[n_e].shape) * sig_noise
        noise_img_gal = np.random.randn(*gals[n_e].shape) * sig_noise

        gal_masked = np.copy(gals[n_e])
        if len(np.where(weight_map == 0)[0]) != 0:
            gal_masked[weight_map == 0] = noise_img_gal[weight_map == 0]

        weight_map *= 1 / sig_noise**2

        # Original PSF fit
        w_tmp = np.sum(weight_map)
        psf_res_gT["g_PSFo"] += psf_res["g"] * w_tmp
        psf_res_gT["g_err_PSFo"] += (
            np.array([psf_res["pars_err"][2], psf_res["pars_err"][3]]) * w_tmp
        )
        psf_res_gT["T_PSFo"] += psf_res["T"] * w_tmp
        psf_res_gT["T_err_PSFo"] += psf_res["T_err"] * w_tmp
        wsum += w_tmp

        gal_obs = Observation(
            gal_masked,
            weight=weight_map,
            jacobian=gal_jacob,
            psf=psf_obs,
            noise=noise_img,
        )

        if gal_guess_flag:
            # gal_guess_tmp[:2] = 0
            final_gal_guess = np.copy(gal_guess_tmp)
            final_gal_guess[:2] = 0
            gal_guess.append(final_gal_guess)

        gal_obs_list.append(gal_obs)
        T_guess_psf.append(psf_T)
        gal_guess_flag = True

    if wsum == 0:
        raise ZeroDivisionError("Sum of weights = 0, division by zero")

    # Normalize PSF fit output
    for key in psf_res_gT.keys():
        psf_res_gT[key] /= wsum

    # Gal guess handling
    fail_get_guess = False
    if len(gal_guess) == 0:
        fail_get_guess = True
        gal_pars = [0.0, 0.0, 0.0, 0.0, 1, 100]
    else:
        gal_pars = np.mean(gal_guess, 0)

    psf_model = "gauss"
    gal_model = "gauss"

    # metacal specific parameters
    metacal_pars = {
        "types": ["noshear", "1p", "1m"],
        "step": 0.01,
        "psf": "gauss",
        "fixnoise": True,
        "use_noise_image": True,
    }

    Tguess = np.mean(T_guess_psf)

    # retry the fit twice
    obs_dict_mcal = ngmix.metacal.get_all_metacal(gal_obs_list, **metacal_pars)
    res = {"mcal_flags": 0}

    ntry = 5

    for key in sorted(obs_dict_mcal):

        fres = make_galsimfit(
            obs_dict_mcal[key], gal_model, gal_pars, prior=prior
        )
        # print(fres["pars"][0:2])
        # fres = make_fit(
        #     obs_dict_mcal[key], gal_model, gal_pars, prior=prior
        # )

        res["mcal_flags"] |= fres["flags"]
        tres = {}

        for name in fres.keys():
            tres[name] = fres[name]
        tres["flags"] = fres["flags"]

        wsum = 0
        Tpsf_sum = 0
        gpsf_sum = np.zeros(2)
        npsf = 0
        for obs in obs_dict_mcal[key]:

            if hasattr(obs, "psf_nopix"):
                try:
                    psf_res = make_galsimfit(
                        obs.psf_nopix,
                        psf_model,
                        np.array([0.0, 0.0, 0.0, 0.0, Tguess, 1.0]),
                        ntry=ntry,
                    )
                except Exception:
                    continue
                g1, g2 = psf_res["g"]
                T = psf_res["T"]
            else:
                try:
                    psf_res = make_galsimfit(
                        obs.psf,
                        psf_model,
                        np.array([0.0, 0.0, 0.0, 0.0, Tguess, 1.0]),
                    )
                except Exception:
                    continue
                g1, g2 = psf_res["g"]
                T = psf_res["T"]

            # TODO we sometimes use other weights
            twsum = obs.weight.sum()

            wsum += twsum
            gpsf_sum[0] += g1 * twsum
            gpsf_sum[1] += g2 * twsum
            Tpsf_sum += T * twsum
            npsf += 1

        tres["gpsf"] = gpsf_sum / wsum
        tres["Tpsf"] = Tpsf_sum / wsum

        res[key] = tres

    # result dictionary, keyed by the types in metacal_pars above
    metacal_res = res

    metacal_res.update(psf_res_gT)
    metacal_res["moments_fail"] = fail_get_guess

    return metacal_res


# ### Extra func modified from: [metacal example](https://github.com/esheldon/ngmix/blob/master/examples/metacal/metacal.py)

# In[19]:


def progress(total, miniters=1):
    last_print_n = 0
    last_printed_len = 0
    sl = str(len(str(total)))
    mf = '%'+sl+'d/%'+sl+'d %3d%%'
    for i in range(total):
        yield i

        num = i+1
        if i == 0 or num == total or num - last_print_n >= miniters:
            meter = mf % (num, total, 100*float(num) / total)
            nspace = max(last_printed_len-len(meter), 0)

            print('\r'+meter+' '*nspace, flush=True, end='')
            last_printed_len = len(meter)
            if i > 0:
                last_print_n = num

    print(flush=True)


# In[20]:


def make_struct(res, shear_type):
    """
    make the data structure

    Parameters
    ----------
    res: dict
        With keys 's2n', 'e', and 'T'
    obs: ngmix.Observation
        The observation for this shear type
    shear_type: str
        The shear type

    Returns
    -------
    1-element array with fields
    """
    dt = [
        ('flags', 'i4'),
        ('shear_type', 'U7'),
        ('s2n', 'f8'),
        ('g', 'f8', 2),
        ('T', 'f8'),
        ('Tpsf', 'f8'),
    ]
    data = np.zeros(1, dtype=dt)
    data['shear_type'] = shear_type
    data['flags'] = res['flags']
    if res['flags'] == 0:
        s2n = res["flux"]/res["flux_err"]
        data['s2n'] = s2n
        # for moments we are actually measureing e, the elliptity
        data['g'] = res['g']
        data['T'] = res['T']
    else:
        data['s2n'] = np.nan
        data['g'] = np.nan
        data['T'] = np.nan
        data['Tpsf'] = np.nan

        # we only have one epoch and band, so we can get the psf T from the
        # observation rather than averaging over epochs/bands
        data['Tpsf'] = res["Tpsf"]

    return data


# In[21]:


def select(data, shear_type):
    """
    select the data by shear type and size

    Parameters
    ----------
    data: array
        The array with fields shear_type and T
    shear_type: str
        e.g. 'noshear', '1p', etc.

    Returns
    -------
    array of indices
    """

    w, = np.where(
        (data['flags'] == 0) & (data['shear_type'] == shear_type)
    )
    return w


# In[22]:


gals, psfs, psfs_sigmas, weights, flags, jacob_lists = make_data(rng=np.random.RandomState(10), noise=1.e-10, shear=(0.02, 0.0), n_epochs=3, share_shift=False)


# In[23]:


pixel_scale = 0.1857
prior = get_prior(pixel_scale, rng=np.random.default_rng(42))


# In[24]:


res = do_ngmix_metacal(gals, psfs, psfs_sigmas, weights, flags, jacob_lists, prior=prior, pixel_scale=pixel_scale, sig_noise=1e-10)


# In[25]:


sig_noise = 1e-10
ntrial = 50
seed = 42
shear_types = ["noshear", "1p", "1m"]
rng = np.random.RandomState(seed)
seeds = rng.randint(0, 2**30, size=ntrial)

dlist_p_ = []
dlist_m_ = []
for i in progress(ntrial, miniters=10):

    d_ = []
    for shear_true in [(0.02, 0.00), (-0.02, 0.00)]:

        gals, psfs, psfs_sigmas, weights, flags, jacob_lists = make_data(rng=np.random.RandomState(seeds[i]), noise=sig_noise, shear=shear_true, n_epochs=3, share_shift=False)
        try:
            resdict = do_ngmix_metacal(gals, psfs, psfs_sigmas, weights, flags, jacob_lists, prior=prior, pixel_scale=pixel_scale, sig_noise=sig_noise)
            stt = []
            for stype in shear_types:
                st = make_struct(res=resdict[stype], shear_type=stype)
                stt.append(st)
        except:
            continue
        d_.append(np.hstack(stt))
    if len(d_) != 2:
        continue
    dlist_p_.extend(d_[0])
    dlist_m_.extend(d_[1])

print()

data_p = np.hstack(dlist_p_)
data_m = np.hstack(dlist_m_)

# P shear
w = select(data=data_p, shear_type='noshear')
w_1p = select(data=data_p, shear_type='1p')
w_1m = select(data=data_p, shear_type='1m')

g_p = data_p['g'][w]
g1_1p_p = data_p['g'][w_1p, 0]
g1_1m_p = data_p['g'][w_1m, 0]
R11_p = np.atleast_2d((g1_1p_p - g1_1m_p)/0.02).T

# M shear
w = select(data=data_m, shear_type='noshear')
w_1p = select(data=data_m, shear_type='1p')
w_1m = select(data=data_m, shear_type='1m')

g_m = data_m['g'][w]
g1_1p_m = data_m['g'][w_1p, 0]
g1_1m_m = data_m['g'][w_1m, 0]
R11_m = np.atleast_2d((g1_1p_m - g1_1m_m)/0.02).T

shear_ = (g_p - g_m) / (R11_p + R11_m)
shear = np.mean(shear_, axis=0)
shear_err = np.std(shear_, axis=0)/np.sqrt(len(shear_))

m = shear[0]/0.02-1
merr = shear_err[0]/0.02

s2n = data_p['s2n'][w].mean()

print('S/N: %g' % s2n)
print('R11: %g %g' % (np.mean(R11_p), np.mean(R11_m)))
print('m[1e-3, 3sigmas]: %g +/- %g (99.7%% conf)' % (m/1e-3, merr*3/1e-3))
print('c[1e-5, 3sigmas]: %g +/- %g (99.7%% conf)' % (shear[1]/1e-5, shear_err[1]*3/1e-5))


# # Below is from Axel

# ## Before the fix
# 
# S/N: 342.343
# 
# R11: 0.926435 0.926423
# 
# m[1e-3, 3sigmas]: -27.8983 +/- 6.50568 (99.7% conf)
# 
# c[1e-5, 3sigmas]: 0.0068871 +/- 0.18067 (99.7% conf)

# ## After the fix
# 
# 50/50 100%
# 
# S/N: 1006.71
# 
# R11: 0.922675 0.922739
# 
# m[1e-3, 3sigmas]: 0.308435 +/- 0.597766 (99.7% conf)
# 
# c[1e-5, 3sigmas]: 0.0981643 +/- 0.435279 (99.7% conf)
# 

# 
