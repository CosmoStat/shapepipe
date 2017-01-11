/******************************************************************************
**                   Copyright (C) 1994 by CEA
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
**    File:  MR_Noise.h
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

#ifndef __NOISE__
#define __NOISE__

#include "IM_Obj.h"
#include "MR_Obj.h"
#include "IM_Noise.h"
#include "MR_Support.h"

extern float TabSignificantLevel[MAX_BAND];         

void noise_compute (MultiResol &MR_Data);
void noise_compute (int Nbr_Plan=DEFAULT_NBR_SCALE, 
                    type_transform Transform=DEFAULT_TRANSFORM,
                    int Nl=DEFAULT_NL, int Nc=DEFAULT_NC,
                    int WindowMedianSize=DEFAULT_MEDIAN_WINDOW_SIZE);

float mr_noise_estimation (MultiResol &MR_Data);
float mr_level_noise (MultiResol &MR_Data, int s, 
                      float N_Sigma=DEFAULT_N_SIGMA);
void mr_threshold (MultiResol &MR_Data, float &Noise_Ima, 
                  float N_Sigma = DEFAULT_N_SIGMA);
void mr_support_threshold (MultiResol &MR_Data, float Noise_Ima, 
                           float N_Sigma = DEFAULT_N_SIGMA, 
                           Bool WriteSupport=DEFAULT_WRITE_SUPPORT);
Bool mallat_significant (MultiResol &MR_Data, int s, int i, int j,
                         float Level, details Det=D_NULL);
#endif
