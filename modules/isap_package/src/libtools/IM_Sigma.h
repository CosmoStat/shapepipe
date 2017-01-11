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


#ifndef __IM_SIGMA__
#define __IM_SIGMA__

#include"IM_Noise.h"

#define NBR_METHODS 6

enum type_method {METHOD_SIGMA_CLIPPING,METHOD_MEDIANE, METHOD_BSPLINE, METHOD_MR_SUPPORT, METHOD_BLOCK, METHOD_MAD, METHOD_SCATTER};

inline const char * StringMethod (type_method type)
{
    switch (type)
    {
	case METHOD_SIGMA_CLIPPING:
		return("3_sigma_clipping");break;
        case METHOD_MEDIANE: 
              return ("Mediane + 3_sigma_clipping ");break;
        case METHOD_BSPLINE: 
              return ("Bspline + 3_sigma_clipping");break;
        case METHOD_MR_SUPPORT: 
              return ("Multiresolution support");break;
        case METHOD_BLOCK: 
              return ("Block method");break;
        case METHOD_MAD: 
              return ("MAD method (Median of Absolute Deviation)");break;
        case METHOD_SCATTER: 
              return ("Scatter method");break;
	    }
        return ("Error: bad method fpr noise estimation");
  }

inline void sigma_method_usage(type_method Method)
{
  fprintf(OUTMAN, "         [-m type_of_methods]\n");
  for (int i = 0; i < NBR_METHODS; i++)
    fprintf(OUTMAN, "              %d: %s \n",i+1,StringMethod((type_method)i));
  fprintf(OUTMAN, "              default is %s.\n", StringMethod((type_method) Method));
  
}

#define DEFAULT_SIGMA_METHOD METHOD_MEDIANE
#define DEFAULT_MAX_ITER_METHOD 10
#define DEFAULT_SIZE_BLOCK 7

float detect_noise_sigma ( Ifloat &, Bool Average_Non_Null=True, int Nit=3);
float detect_noise_sigma (const Ifloat &, Bool Average_Non_Null=True, int Nit=3);
float detect_noise_from_med (const Ifloat &Image);
float detect_noise_from_bspline (const Ifloat &Image);
float detect_noise_from_mad (Ifloat &Image, Bool CenterOnly=False);
float detect_noise_sigma (const Iint &, Bool Average_Non_Null=True, int Nit=3);
float detect_noise_from_block (Ifloat &Ima, int SizeBlock=DEFAULT_SIZE_BLOCK);
void sigma_clip(const Ifloat &Image, float &Mean, float &Sigma, int Nit=3);
void sort_bulle (float *Data, int Np);

void block_transform(Ifloat &Ima,float * & PtrMean, float * & PtrSigma, 
                     int &Nb, int SizeBlock=DEFAULT_SIZE_BLOCK);
float sigma_clip_int (int *Image, int Nl, int Nc, int Opt=1, int Niter=3);
void  sigma_clip_int (int *Image, float & Moyenne, float & Sigma, 
                     int Nl, int Nc, int Option=1, int Nit=3);
#endif




