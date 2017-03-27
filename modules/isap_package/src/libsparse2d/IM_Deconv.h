/******************************************************************************
**                   Copyright (C) 1995 by CEA
*******************************************************************************
**
**    UNIT
**
**    Version: 3.2
**
**    Author: Jean-Luc Starck
**
**    Date:  96/05/22 
**    
**    File:  IM_Deconv.h
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

#ifndef __IM_DECONV__
#define __IM_DECONV__

#include "GlobalInc.h"

#define NBR_DECONV 16     // number of implemented methods
#define NBR_STD_DECONV 11  // number of standard deconvolution methods
#define NBR_STD_WAVE 4    // number of standard deconvolution methods
                          // regularized by the wavelets
#define NBR_MEM_WAVE 2    // number of MEM wavelet methods

enum type_deconv {DEC_CITTERT, DEC_GRADIENT, DEC_INVERSE,
                  DEC_LUCY, DEC_CLEAN, DEC_MEM, DEC_MEM_MODEL,
		  DEC_TIKHONOV, DEC_MAP, DEC_MARKOV, DEC_MARKOV_LUCY,
		  DEC_MR_CITTERT, DEC_MR_GRADIENT,
                  DEC_MR_LUCY, DEC_MR_MAP, DEC_MR_MEM, DEC_MR_MEM_NOISE,
 		  DEC_MR_CLEAN, DEC_MR_VAGUELET, DEC_MR_INVERSE};

#define DEFAULT_STD_DECONV  DEC_GRADIENT
          
inline const char * StringDeconv (type_deconv type)
{
    switch (type)
    {
     case DEC_MARKOV_LUCY:
      return ("Deconvolution by Lucy's algorithm + Markov Random Field Regularization");break;
     case DEC_MARKOV:
      return ("Deconvolution by the Gradient algorithm + Markov Random Field Regularization");break;
     case DEC_MR_MEM:
      return ("Deconvolution by the Multiscale entropy method (N1-MSE)");break;
     case DEC_MR_MEM_NOISE:
      return ("Deconvolution by the Multiscale entropy method  (N2-MSE)");break; 
     case DEC_MEM: 
      return ("Deconvolution by the MEM method (Maximize -sum O log O)");break;
     case DEC_MEM_MODEL: 
      return ("Deconvolution by the MEM method (Gull entropy)");break; 
     case DEC_TIKHONOV: 
      return ("Deconvolution using Tikhonov Regularization");break;
     case DEC_MAP: 
      return ("Deconvolution using MAP method");break;
    case DEC_CITTERT: 
      return ("Deconvolution by Van Cittert's algorithm");break;
    case DEC_GRADIENT: 
      return ("Deconvolution by gradient algorithm");break;
    case DEC_INVERSE: 
      return ("Deconvolution by division in Fourier space");break;
    case DEC_LUCY: 
      return ("Deconvolution by Lucy's algorithm");break;
    case DEC_CLEAN:
      return ("Deconvolution by CLEAN algorithm");break;
    case DEC_MR_CITTERT: 
      return ("Deconvolution by multiresolution Van Cittert's algorithm");break;
    case DEC_MR_GRADIENT:
      return ("Deconvolution by multiresolution gradient's algorithm");break;
    case DEC_MR_LUCY:
      return ("Deconvolution by multiresolution Lucy algorithm");break;
    case DEC_MR_MAP:
      return ("Deconvolution by multiresolution MAP algorithm");break;
    case DEC_MR_CLEAN:
      return ("Deconvolution by multiresolution CLEAN algorithm");break;
    case DEC_MR_VAGUELET:
      return ("Deconvolution by the division in Fourier space + Wavelet filtering");break;
    case DEC_MR_INVERSE:
      return ("Deconvolution by multiresolution division in Fourier space");break;
    }
    return ("Error: bad method of deconvolution");
}

// return ("Deconvolution by optimal gradient method and multiresolution regularization");break;

inline void deconv_usage(type_deconv Deconv)
{
    fprintf(OUTMAN, "         [-d type_of_deconvolution]\n");
    for (int i = 0; i < NBR_DECONV; i++)
    fprintf(OUTMAN, "              %d: %s \n",i+1,StringDeconv((type_deconv)i));
    fprintf(OUTMAN, "              default is %s\n", StringDeconv((type_deconv) Deconv));
}

#define DEFAULT_MAX_ITER_DECONV 500
#define DEFAULT_EPSILON_DECONV 0.001
#define DEFAULT_DECONV DEC_MR_LUCY 

#define DEFAULT_BLIND_ITER_PER_STEP 10
#define DEFAULT_BLIND_MAX_ITER 50
#define DEFAULT_BLIND_FWHM 3.

#define DEFAULT_CLEAN_NITER 100000
#define DEFAULT_CLEAN_KAVER 1
#define DEFAULT_CLEAN_ADDRESI False
#define DEFAULT_CLEAN_GAMMA 0.01
#define DEFAULT_CLEAN_FWHM 0.1

void dec_pos_max (Ifloat &Tab, int &Ind_i, int &ind_j , float &Val_Max, 
                  Bool SearchPositiv=False);
void dec_line_column (int N, int &N_Out);
void dec_center_psf (Ifloat &D_Beam, Ifloat &Psf);
void dec_convol_conj (Ifloat &Resi, Icomplex_f &Psf_cf);
void dec_convol (Ifloat &Resi, Icomplex_f &Psf_cf, Ifloat &Result);
void dec_inverse (Ifloat &Resi, Icomplex_f &Psf_cf, Ifloat &Result, float Eps=0.001);

void psf_get(Ifloat &InPsf, Icomplex_f &Psf_cf, int ImaNl, int ImaNc, 
             Bool PsfMaxShift);
void psf_convol (Ifloat &Imag, Icomplex_f &Psf_cf, Ifloat &Imag_out, Bool FluxNorm=True);
void psf_convol (Ifloat &Imag, Icomplex_f &Psf_cf, Bool FluxNorm=True);
void psf_convol (Ifloat &Imag, Ifloat &Psf, Ifloat &Imag_out, Bool FluxNorm=True,Bool PsfMaxShift=True);
void psf_convol (Ifloat &Imag, Ifloat &Psf, Bool FluxNorm=True);
void psf_convol_conj (Ifloat &Imag, Icomplex_f &Psf_cf, Ifloat &Imag_out);
void psf_convol_conj (Ifloat &Imag, Icomplex_f &Psf_cf);

void dec_im_standart (Ifloat &Imag, Ifloat &Obj, Icomplex_f &Psf_cf,
                     float Eps_cv = DEFAULT_EPSILON_DECONV,
                     int Nbr_Iter = DEFAULT_MAX_ITER_DECONV, 
                     type_deconv Deconv = DEFAULT_DECONV, 
                     float TotalFlux=0.,
                     float Noise_Ima=0.);

void dec_im_support (Ifloat &Imag, Ifloat &Obj, Icomplex_f &Psf_cf,
               Ifloat &Support, float Eps_cv = DEFAULT_EPSILON_DECONV, 
               int Nbr_Iter = DEFAULT_MAX_ITER_DECONV, 
               type_deconv Deconv = DEFAULT_DECONV);


void dec_clean (Ifloat & Imag, Ifloat & Psf, Ifloat & Imag_Out, 
                float Fwhm, float Noise, float Gamma=DEFAULT_CLEAN_GAMMA,
                int Niter=DEFAULT_CLEAN_NITER, int K_Aver=DEFAULT_CLEAN_KAVER,
                Bool AddResi=True, Bool SearchPositiv=True,
		Bool UseEnergy=False, Bool Verbose=False);
                
#endif
