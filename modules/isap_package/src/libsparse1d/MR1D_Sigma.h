/******************************************************************************
**                   Copyright (C) 1996 by CEA
*******************************************************************************
**
**    UNIT
**
**    Version: 0.1
**
**    Author: Rene Gastaud & Jean-Luc Starck
**
**    Date:  22-OCT-1996
**    
**    File:  MR1D_NoiseModel.h
**
*******************************************************************************
**
**    DESCRIPTION  Noise determination routine
**    ----------- 
**                 
******************************************************************************/

#ifndef _SIGMA1DNOISE_H_
#define _SIGMA1DNOISE_H_

#include "MR1D_Obj.h"  
 
#define NBR_SIGMA_METHOD 4
enum type_sigma_method_1d {SIGMA_MEDIAN_1D, SIGMA_BSPLINE_1D, SIGMA_CLIPIMA_1D,
                        SIGMA_SUPPORT_1D, SIGM_UNDEFINED_1D=-1};
#define DEF_SIGMA_METHOD SIGMA_MEDIAN



#define DEF_NSIG 3.   // Default Nsigma detection
#define DEF_NSIG0 4. // Default Nsigma detection at the first scale


void  mr1d_noise_compute (int Nbr_Plan, type_trans_1d Transform, 
                          type_border Border=DEFAULT_BORDER, 
                          int MedianWinSize=DEFAULT_MEDIAN_1D_WINSIZE,
                          sb_type_norm  Norm = NORM_L1);

// set the internal table (Tab1DSignifLevel) following the input parameter

float mr1d_level_noise (MR_1D &MR_Data, int s, 
                        float N_Sigma=DEF_NSIG, float N_Sigma0=DEF_NSIG0);
// return for a given scale, the detection level

float mr1d_tab_noise(int s);
// return the noise standard deviation of the noise at a given scale
// considering a white noise of sigma=1 in the data

float detect1d_noise_from_med (fltarray &Signal, int NiterSigmaClip=3);
// estimate the noise standard deviation in the data by sigma clipping
// on the data set = Signal - median(Signal)

float mr1d_detect_noise(fltarray &Sig, int NbrScale=5, float N_Sigma=3, 
                        int Niter=10, int NiterClip=3);
// estimate the noise standard deviation from the multiresolution support

#endif
