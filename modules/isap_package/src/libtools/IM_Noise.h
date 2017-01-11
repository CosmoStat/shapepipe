/******************************************************************************
**                   Copyright (C) 1995 by CEA
*******************************************************************************
**
**    UNIT
**
**    Version: 1.0
**
**    Author: Jean-Luc Starck
**
**    Date:  96/06/13 
**    
**    File:  IM_Noise.h
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

#ifndef __IM_NOISE__
#define __IM_NOISE__

#include "IM_Obj.h"


#define NBR_NOISE 10
#define NBR_GAUSS_TRANSF_NOISE 8

enum type_noise {NOISE_GAUSSIAN, NOISE_POISSON, NOISE_GAUSS_POISSON,
                 NOISE_MULTI, NOISE_NON_UNI_ADD, NOISE_NON_UNI_MULT, 
                 NOISE_UNI_UNDEFINED, NOISE_UNDEFINED, NOISE_CORREL,
		 NOISE_EVENT_POISSON, NOISE_SPECKLE};

const char * StringNoise (type_noise type);


inline void noise_usage()
{
    fprintf(OUTMAN, "         [-m type_of_noise]\n");
    for (int i = 0; i < NBR_NOISE; i++)
    fprintf(OUTMAN, "              %d: %s \n",i+1,
                                            StringNoise((type_noise )i));
    fprintf(OUTMAN, "             default is Gaussian noise\n");
}

inline void rms_noise_usage()
{
    fprintf(OUTMAN, "\n         [-R RMS_Map_File_Name]\n");
    fprintf(OUTMAN, "              RMS Map (only used with -m 5 and -m 9 options). \n");
}

#define DEFAULT_STAT_NOISE NOISE_GAUSSIAN
#define DEFAULT_SIZE_BLOCK_SIG 7

/* procedures for noise treatment */
void noise_poisson_transform (fltarray &Resi, fltarray &Sn);
void noise_poisson_transform (int *Resi, int *Sn, int Nl, int Nc);

void noise_inverse_poisson_transform (fltarray &Sn, fltarray &Resi);
void noise_inverse_poisson_transform (int *Sn, int *Resi, int Nl, int Nc);

void noise_poisson_transform (Ifloat &Imag_in, Ifloat &Imag_out);
void noise_inverse_poisson_transform (Ifloat &Imag_in, Ifloat &Imag_out);
void noise_residu_estimation (Ifloat &In, Ifloat &Sn, Ifloat &Rn);
void noise_threshold_imag (Ifloat &Imag, float Level);
void noise_log_transform (Ifloat &Imag, Ifloat &Result);
void noise_inv_log_transform (Ifloat &Imag, Ifloat &Result);
void noise_support_threshold_imag (Ifloat &Imag, float Level, 
                                   Ifloat &Support, 
                                   Bool WriteSupport=False);
float ran1(int *idum);
float gammln(float xx);
float poidev(float xm,int *idum);
void im_poisson_noise(Ifloat &Data);
void im_sigma(Ifloat &Data, Ifloat &Ima_Sigma, 
              int SizeBlock=DEFAULT_SIZE_BLOCK_SIG, int Nit=1);
void im_sigma_mad(Ifloat &Data, Ifloat &Ima_Sigma, int SizeBlock);
void im_sigma_block_mad(Ifloat &Data, Ifloat &Ima_Sigma, int SizeBlock, int Resol);
void im_sigma_block(Ifloat &Data, Ifloat &Ima_Sigma, int SizeBlock, int Nit, int Resol);
#endif
