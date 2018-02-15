/******************************************************************************
**                   Copyright (C) 2003 by CEA
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
**    File:  IM_Morpho.h
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

#ifndef _MORPHO_H_
#define _MORPHO_H_

#include "GlobalInc.h"
#include "Border.h"
#include <cmath>


void morphoi_erosion (int *Imag1, int *Imag2, int Nl, int Nc, int Window_Size = 3);
void morphoi_dilation (int *Imag1, int *Imag2, int Nl, int Nc, int Window_Size = 3);
void morphoi_cercle_erosion (int *Imag1, int *Imag2, int Nl, int Nc, int Window_Size=5);
void morphoi_cercle_dilation (int *Imag1, int *Imag2, int Nl, int Nc, int Window_Size=5);
 
//void morphoi_dilation_class (Iint &Imag1, Iint &Imag2, int Window_Size=3);
//void cerclei_erosion_class (Iint &Imag1, Iint &Imag2, int Window_Size=5);
//void cerclei_dilation_class (Iint &Imag1, Iint &Imag2, int Window_Size=5);

void morpho_erosion (Ifloat &Imag1, Ifloat& Imag2, int Window_Size = 3);
void morpho_dilation (Ifloat &Imag1, Ifloat &Imag2, int Window_Size = 3);
void morpho_cercle_erosion (Ifloat &Imag1, Ifloat& Imag2, int Window_Size=5);
void morpho_cercle_dilation (Ifloat &Imag1, Ifloat &Imag2, int Window_Size=5);

void im_detect(Ifloat &Imag_in, Ifloat &Imag_Detect, 
               float **x_ef, float **y_ef,  float **f_ef, int *nmax);
void im_segment (Ifloat &Imag, Ifloat &Segment, int &NLabel, 
                float Level=0.5, Bool CleanBord = False, int FirstLabel=1);
void im_thin(Ifloat &Imag,float Level);

void morpho4_dilation (Ifloat &Imag1, Ifloat &Imag2, int Step_trou=0);
void morpho4_erosion (Ifloat &Imag1, Ifloat& Imag2, int Step_trou=0);

#endif

