/******************************************************************************
**                   Copyright (C) 1998 by CEA
*******************************************************************************
**
**    UNIT
**
**    Version: 1.0
**
**    Author: Jean-Luc Starck
**
**    Date:  23/12/98 
**    
**    File:  MW_Deconv.h
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


#ifndef  __MWDECONV__
#define  __MWDECONV__

const int DEF_DEC_MAX_ITER = 500;
const int DEF_NB_ITER_SUP_CAL = 5;

// different type of data weighting
const int DEC_NBR_WEIGHT = 3;
enum  type_mem_weight {DEC_NO_WEIGHT,  DEC_WEIGHT_PROBA, DEC_WEIGHT_SUPPORT};
const type_mem_weight DEF_DEC_WEIGHT = DEC_WEIGHT_SUPPORT;

// different type of regularization
const int DEC_NBR_REGUL = 5;
enum type_mem_regul {DEC_NO_REGUL, DEC_REGUL_NO_PROTECT, DEC_REGUL_PROBA, 
                     DEC_REGUL_SUPPORT,DEC_REGUL_PROBA_SUPPORT}; 

const type_mem_regul DEF_DEC_REGUL = DEC_REGUL_SUPPORT;

const type_deconv DEF_DEC_ENTROP = DEC_MR_MEM;
const type_deconv DEF_DEC_ENTROP_GAUSS = DEC_MR_MEM_NOISE;

class MWDeconv: public MRDeconv 
{
   float sigma_obj(float SigmaNoise);
   // create a noisy map of standard deviation SigmaNoise
   // convolve it by the PSF
   // return the standard deviation of convolved map
      		    
   void init_mem_param()
   {
      RegulParam = 0;
      MaxIter = DEF_DEC_MAX_ITER;
      DecMethod = DEF_DEC_ENTROP;
      Nscale = 0;
      NbrBand = 0;
      Transform = DEFAULT_TRANSFORM;
      TypeWeight = DEC_NO_WEIGHT;
      TypeRegul = DEC_REGUL_SUPPORT;
      WaveFilterResi=True;
      GetMR_FirstGuess=True;
   }

   // Apply first another regularized iterative method
   // in order to improve the convergence.
   void get_mr_support_first_guess(type_deconv Deconv=DEC_MR_CITTERT, int Niter=10);
   
   CMemWave MemWObj;                      // Gradient class 
   void compute_mem_resi();               // compute the residual
   void mw_grad(Ifloat & Gradient);       // compute the gradient
   float mw_find_optim(Ifloat &Gradient); // optimize the convergence paramter
   void find_support_obj();               // find the multiresolution support
                                          // of the object
					  
   MultiResol MR_Obj; // Multiresolution buffer used in some calculation
   MultiResol MR_RegulObj; // Multiresolution buffer used in some calculation
   int NbIter_SupportCalc; // Number of iterations between two reestimations
                           // of the multiresolution support (def is 5).
   Ifloat Buff;       // Buffer image used during some calcualtion
   Ifloat ObjN;
   Ifloat ObjS;
   float Noise_Obj;   // Noise in the object space
   fltarray TabSigmaObj; // object noise estimation per scale
   fltarray TabAlpha;    // Regularization parameter per scale

    int Nscale; 
    int NbrBand;
    type_transform Transform;   // type of transform
    float SigmaResi         ;   // residual standard deviation
    
   public:
    type_mem_regul TypeRegul;     // type of regularization weighting
    type_mem_weight TypeWeight;   // type of weighting applied to the data
    Bool GetMR_FirstGuess;        // regularized Van Cittert method  is 
                                  // first apply to improve the convergence
    float NSigmaObj;              // Nsigma for object multiresolution 
                                  // support detection

    MWDeconv(): MRDeconv() { init_mem_param();}

    void mw_deconv(Ifloat *FirstGuess, Ifloat *ICF);
   // deconvolution by the multiscale entropy method
   
    void mw_vaguelet(Bool DataSNR, Ifloat *FirstGuess, Ifloat *ICF);
    
   ~MWDeconv(){};
   };

#endif




