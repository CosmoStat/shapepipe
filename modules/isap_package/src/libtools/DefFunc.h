/******************************************************************************
**                   Copyright (C) 2003 by CEA
*******************************************************************************
**
**    UNIT
**
**    Version: 1.0
**
**    Author: Jean-Luc Starck
**
**    Date:  30/03/03 
**    
**    File:  DefFunc.h
**
******************************************************************************/

#ifndef _DEF_FUNC_H_
#define _DEF_FUNC_H_

#include "GlobalInc.h"

void simu_gaussian_noise(fltarray & NoiseData, float Sigma=1., unsigned int InitRnd=0);
// Gaussian noise simulation. NoiseData can be of dimension 1,2 or 3 and
// must be allocated
void simu_gaussian_noise3d (fltarray & NoiseData, float Sigma=1., unsigned int InitRnd=0);

void make_gaussian1d(fltarray & Result, float Sigma, int Indi=-1);
void make_gaussian2d(fltarray & Result, float Sigma, int Indi=-1, int Indj=-1);
void make_gaussian3d(fltarray & Result, float Sigma, int Indi=-1, int Indj=-1, int Indk=-1);
void make_gaussian(fltarray & Result, float Sigma);
// Create a Gaussian function. By default, it is centered at the center
// of the array. Result  must be allocated and and be of dimension 1,2 or 3.

void msvst_transform (fltarray &data, fltarray &vstdata, int scale, Bool coupled=True);
void msvst_inv_transform(fltarray &vstdata, fltarray &data, int scale, Bool Biais=True);
void msvst_transform (Ifloat &data, Ifloat &vstdata, int scale, Bool coupled=True);
void msvst_inv_transform(Ifloat &vstdata,  Ifloat &data, int scale, Bool Biais=True);
double msvst_var (int dim, int scale);
void survival (fltarray & Data, fltarray &Survival, int NpixSurv=1001, float NSigma=10., Bool Norm=True, float SigmaNorm=0., float MeanNorm=-1.);
void get_stat(double *x, int N, dblarray &TabStat);

// Compute the survival function of the Data.
// By default, the data are centered:  X = ( Data - Mean(Data) ) / Sigma(Data) 
// The mean and the standard deviation can by be given.
#endif
