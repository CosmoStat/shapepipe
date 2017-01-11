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
**    File:  MR_Support.h
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

#ifndef __SUPPORT__
#define __SUPPORT__

#include "MR_Obj.h"
#include "IM_Noise.h"

#define DEFAULT_WRITE_SUPPORT False

/* Multiresolution support procedures */
void mr_create_support (MultiResol &MR_Data, float &Noise_Ima, 
                        float N_Sigma=DEFAULT_N_SIGMA);
void mr_dilate_support (Bool SmoothSupport=False);
void mr_smooth_support ();
void mr_hierar_support ();
void mr_add_map_support (Ifloat &MapPos);
void mr_support (Ifloat &Imag, float &Noise_Ima, type_transform Transform, 
                 int Nbr_Plan, float N_Sigma=DEFAULT_N_SIGMA, 
                 type_noise Stat_Noise=DEFAULT_STAT_NOISE);
#endif
