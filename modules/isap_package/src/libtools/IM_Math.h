/******************************************************************************
**                   Copyright (C) 1994 by CEA
*******************************************************************************
**
**    UNIT
**
**    Version: 3.2
**
**    Author: Jean-Luc Starck
**
**    Date:  96/05/07 
**    
**    File:  IM_Math.h
**
*******************************************************************************
**
**    DESCRIPTION  
**    ----------- 
**                 
**    PARAMETRES    
**    ----------    
** 
**    RESULTS      
**    -------  
**
**
******************************************************************************/

#ifndef _IM_MATH_H_
#define _IM_MATH_H_

#include "GlobalInc.h"
#include "IM_Obj.h"

inline double im_b3_spline (double x) {return b3_spline (x);}

const Ifloat im_gaussian (int Nl, int Nc, float Fwhm, int Indi=-1, int Indj=-1);
void make_gaussian(Ifloat & Result, float Sigma, int Indi=-1, int Indj=-1);
void make_gaussian_fwhm(Ifloat & Result, float Fwhm, int Indi=-1, int Indj=-1);

//const fltarray randomn(int N, float Sigma=1., unsigned int InitRnd=100);
void  randomn(fltarray& tab, int N, float Sigma=1., unsigned int InitRnd=100);
void  im_noise_poisson (Ifloat & NoiseIma, float Gain=1, unsigned int InitRnd=100);
const Ifloat im_noise (float Sigma, int Nl, int Nc, unsigned int InitRnd=100);
void  im_noise_rayleigh(Ifloat & NoiseIma, int  NbrImage=1, unsigned int InitRnd=100);
void  im_noise_laplace(Ifloat & NoiseIma, int  NbrImage=1, unsigned  int InitRnd=100);
void  im_noise_gaussian (Ifloat & NoiseIma, float Sigma=1, unsigned  int InitRnd=100);
void im_moment4(Ifloat & Ima, double &Mean, double &Sigma, 
                double &Skew, double & Curt, float & Min, 
	        float & Max, int BorderSize=0);
void im_moment4(fltarray & Tab, double &Mean, double &Sigma, 
                double &Skew, double & Curt, float & Min, 
	        float & Max, int BorderSize=0);
void im_hc_test(Ifloat & Ima, float &HC1, float &HC2, int BorderSize=0, Bool ZeroMean=False);
void im_hc_test(fltarray & Tab,  float &HC1, float &HC2, int BorderSize=0, Bool ZeroMean=False);
void im_gaussianity_test(fltarray &Band, float &T1, float &T2, int BorderSize=0);
void im_gaussianity_test(Ifloat  &Band, float &T1, float &T2, int BorderSize=0);

#endif
