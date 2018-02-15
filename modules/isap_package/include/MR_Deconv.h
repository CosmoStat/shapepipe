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
**    Date:  96/06/13 
**    
**    File:  MR_Deconv.h
**
*******************************************************************************
**
**    DESCRIPTION  Definition for multiresolution deconvolution
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

#ifndef __MR_DECONV__
#define __MR_DECONV__

#include "GlobalInc.h"
#include "IM_Deconv.h"
#include "MR_Obj.h"
#include "MR_Obj.h"
#include "MR_Support.h"
#include "MR_Noise.h"
#include "MR_Filter.h"
#include "IM_Regul.h"
#include "IM_IO.h"
#include "IM_Math.h"
#include "IM_Noise.h"
#include "IM_Deconv.h"
#include "MR_Sigma.h"
#include "MR_NoiseModel.h"
#include "MR_Filter.h"
#include "NR.h"
#include "FFTN_2D.h"
#include "MR_SoftRegul.h"


enum type_init_deconv {DEC_INIT_ZERO, DEC_INIT_FLAT, DEC_INIT_IMA, DEC_INIT_GUESS};


class MRDeconv {
  protected:
   FFTN_2D FFT2D;

  // compute the residual
   virtual void compute_resi();
      
   // float find_optim(Ifloat &Gradient, Ifloat & Buff);
   void im_grad(Ifloat &gradient);
   		    
   void init_param();

   int Nl,Nc;         // working image size
   Bool UseICF;       // if True, the ICF is given by the user
   Ifloat Ima_ICF;    // ICF image
   Bool KeepImagn;    // Keep the variable Imag_n
   Bool WaveFilterResi; // if true the residual is filtered in the wavelet
                        // space at each iteration using the multiresolution
			// support

   Ifloat Mult;         // Multiplyer for MEM methods
   float  FluxImag;     // Image flux
   type_init_deconv TypeInit; // Initialisation type
   Bool AdjointRec;     // Use the adjoint for the reconstruction
   Icomplex_f *TabCB_cf; // CLEAN BEAM WT in Fourier space (for MRC only)
   
   // optimization routine
   float find_optim_tikhonov(Ifloat &Gradient);
   float find_optim_xi2(Ifloat &Gradient);  
   float find_optim_poisson(Ifloat &Gradient); 
   float im_find_optim(Ifloat &Gradient);
   
   // compute fonctional value
   float fonctional();
   
  // iterative deconvolution (init_deconv must have been called before)
   void im_iter_deconv();
  
  void init_mrc()      ; // initialization for MultiResolution CLEAN
  void mrc(Ifloat *ICF); // MultiResolution CLEAN deconvolution
  void  mrc_rec_obj(MultiResol &MR_SOL, Icomplex_f *TabCB_cf);
                         // solution reconstruction in Obj
			 // using the dirac list in MR_SOL, and CLEAN BEAM
			 // table TabCB_cf
  void  vaguelet();
  void  va_inverse_data(Ifloat &NoiseSimu);
  
  public:
 
      
   Ifloat Obj;      // solution
   Ifloat Resi;     // residual
   Ifloat Imag;     // input data
   Ifloat Imag_n;  // only used by Lucy method
   Ifloat Psf;      // Point spread function
   // Ifloat ICF;      // Intrinsic Correlation Funciton
   Icomplex_f Psf_cf; // PSF Fourier transform
   
   Bool PositivConstraint; // if true, the solution must positive
   Bool NormFlux;   // if True, the flux is renormalized at each iteration
   float IterCvg;   // congergence paramter applied to the gradient
                    // at each iteration
   Bool OptimParam; // Optimize the convergence parameter at each iteration
                    // if OptimParam == True then IterCvg is mot used
   float EpsCvg;     // congergence parameter for the final convergence
   type_deconv DecMethod; // type of deconvolution method
  
   float RegulParam;  // User regularization parameter
   int MaxIter;       // maximum number of iterations
   int MaxMRCleanIter;// maximum number of iterations for MR CLEAN method
   Bool UseMRCEnergy; // Use standard clean in MRC if UseMRCEnergy = False
                      // otherwise use the ENERGY map in MRC
   Bool CleanLastScale; // if True, CLEAN is also applied on the last scale
                        // otherwise apply a simple Van Citter method
   int CleanFirstScale;// First Scale to CLEAN
   Bool SupportDilate; // Dilate the multiresolution support
                       // (used by MRC method)
		       
   Bool Verbose;      // verbose mode
   Bool PsfMaxShift;  // by default the PSF is shifted to place its max 
                      // at the center.
		      // if PsfMaxShift == False, this is not done
   Bool GaussConv;    // If true, the solution has a limited resolution
                      // (i.e. it must a convolution product between 
		      // a gaussian (fixed by Fwhm) and a hidden solution)
   float Fwhm;        // Full Width at Half Maximum of the Gaussian ICF
   type_noise StatNoise; // type of noise
   float Noise_Ima;   // Gaussian noise standard deviation
   float N_Sigma;     // CLEAN parameter: stop at N_Sigma above the noise
   Ifloat MemModel;   // Model image for mem restoration using Gulling entropy
   Bool UseModel;     // True when a model is used with MEM MODEL method
   type_border Border; // Border type used for multiscale transform calculation
   
   // Parameter for wavelet methods
   Bool KillLastScale;     // Do not restore the last scale
   MultiResol MR_Data;
   MRNoiseModel *ModelData; // noise model class
    
   void convol_gauss(Ifloat &Data, float FWHM);
   // convolve (using the FFT an image by a Gaussian of size FWHM
   
   MRDeconv()  {init_param(); }

   // initialize the deconvolution   
   void init_deconv(Ifloat *FirstGuess=NULL,  Ifloat *ICF=NULL);
       
   // deconvolution 
   void im_deconv(Ifloat *FirstGuess, Ifloat *ICF);
    
   RegulIma RIm; // Regularization Class
   MR_Regul MRIm;

   virtual ~MRDeconv(){};
   };

                 
#endif
