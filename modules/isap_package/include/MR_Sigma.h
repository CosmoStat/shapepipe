/******************************************************************************
**                   Copyright (C) 1995 by CEA
*******************************************************************************
**
**    UNIT
**
**    Version: 3.1
**
**    Author: Jean-Luc Starck
**
**    Date:  96/05/02 
**    
**    File:  MR_Sigma.h
**
*******************************************************************************
**
**    DESCRIPTION  
**    ----------- 
**                 
**    PARAMETRES    
**    ----------    
** 
** 
**    RESULTS      
**    -------  
**
**
******************************************************************************/


#ifndef __SIGMA__
#define __SIGMA__

#include "IM_Obj.h"
#include "IM_Sigma.h"
#include "MR_Obj.h"
#include "MR_Support.h"
#include "IM_Noise.h"
#include "MR_NoiseModel.h"
#include "MR_Noise.h"


float detect_noise_from_support(Ifloat &Imag, int Niter=DEFAULT_MAX_ITER_METHOD,
                                int Nbr_Plan=DEFAULT_NBR_SCALE);
float detect_noise_from_support(Ifloat &Imag, MultiResol & MR_Data, 
                                MRNoiseModel & ModelData, 
                                int Niter=DEFAULT_MAX_ITER_METHOD);

float get_noise (Ifloat &Imag, type_method Method=DEFAULT_SIGMA_METHOD, 
                 int Niter=DEFAULT_MAX_ITER_METHOD,
                 int Nbr_Plan=DEFAULT_NBR_SCALE);

float get_gain(Ifloat &Imag, float Readout, 
               type_method Method=DEFAULT_SIGMA_METHOD,
               float GainMin=0.01, float GainMax=100);

float get_readout_noise(Ifloat &Imag, float Gain, 
               type_method Method=DEFAULT_SIGMA_METHOD, 
               int Nbr_Plan=DEFAULT_NBR_SCALE);

#endif




