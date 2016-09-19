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
**    File:  MR_Filter.h
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


#ifndef __FILTER__
#define __FILTER__

#include "IM_Obj.h"
#include "MR_Obj.h"
#include "MR_Support.h"
#include "IM_Noise.h"
#include "MR_Noise.h"
#include "MR_Rayleigh.h"
#include "MR_NoiseModel.h"
 
// Only the 18 first are available
// The last three did not produce very good results
#define NBR_FILTERING 10

enum type_filter {FILTER_THRESHOLD, 
                  FILTER_SOFT_THRESHOLD,
                  FILTER_ITER_THRESHOLD,
                  FILTER_ITER_ADJOINT_COEFF,
                  FILTER_BIVARIATE_SHRINKAGE,
                  FILTER_MULTI_RES_WIENER,
 		  FILTER_TV_CONSTRAINT,
		  FILTER_WAVELET_CONSTRAINT,
		  FILTER_MULTI_HARD_MAD,
                  FILTER_MULTI_SOFT_MAD,
                  FILTER_HIERARCHICAL_WIENER,
                  FILTER_HIERARCHICAL_TRESHOLD,  
		  FILTER_BAYES_BESSEL};

inline const char * StringFilter (type_filter type)
{
    switch (type)
    {
        case FILTER_BIVARIATE_SHRINKAGE:
	      return ("Bivariate Shrinkage");break;
        case FILTER_BAYES_BESSEL:
              return ("Bayesian Bessel Method");break;
        case FILTER_WAVELET_CONSTRAINT:
              return ("Wavelet Constraint Iterative Methods");break;
        case FILTER_TV_CONSTRAINT:
              return ("Total Variation + Wavelet Constraint");break;
        case FILTER_MULTI_HARD_MAD: 
              return ("Median Absolute Deviation (MAD) Hard Thesholding");break;
        case FILTER_MULTI_SOFT_MAD: 
              return ("Median Absolute Deviation (MAD) Soft Thesholding");break;
        case FILTER_THRESHOLD:
              return ("Multiresolution Hard K-Sigma Thresholding");break;
        case FILTER_SOFT_THRESHOLD: 
              return ("Multiresolution Soft K-Sigma Thresholding");break;
        case FILTER_ITER_THRESHOLD: 
              return ("Iterative Multiresolution Thresholding");break;
        case FILTER_ITER_ADJOINT_COEFF: 
              return ("Adjoint operator applied to the multiresolution support ");break;
        case FILTER_HIERARCHICAL_TRESHOLD: 
              return ("Hierarchical Hard Thresholding");break;
        case FILTER_HIERARCHICAL_WIENER: 
              return ("Hierarchical Wiener Filtering");break;
        case FILTER_MULTI_RES_WIENER:
              return ("Multiresolution Wiener Filtering");break;
    }
   return ("Error: bad type of filtering");
}

inline void filter_usage(type_filter Filter)
{
    fprintf(OUTMAN, "         [-f type_of_filtering]\n");
    for (int i = 0; i < NBR_FILTERING; i++)
    fprintf(OUTMAN, "              %d: %s \n",i+1,StringFilter((type_filter)i));
    fprintf(OUTMAN, "              default is %s.\n", StringFilter(Filter));
}

inline void filterp_usage(type_filter Filter)
{
    fprintf(OUTMAN, "         [-f type_of_filtering]\n");
    fprintf(OUTMAN, "              %d: %s \n",1,StringFilter(FILTER_ITER_THRESHOLD));
    fprintf(OUTMAN, "              %d: %s \n",2,StringFilter(FILTER_ITER_ADJOINT_COEFF));
    fprintf(OUTMAN, "              %d: %s \n",3,StringFilter(FILTER_THRESHOLD));

    fprintf(OUTMAN, "              default is %s.\n", StringFilter((type_filter) Filter));

}


#define DEFAULT_MR_TRANS TO_PAVE_BSPLINE
#define DEFAULT_MAX_ITER_FILTER 10
#define DEFAULT_EPSILON_FILTERING 1e-3
 
#define DEFAULT_FILTERING FILTER_ITER_THRESHOLD
#define DEFAULT_CVG_COEF_MEM 0.05      /* default convergence parameter
                                        for mem filtering */
#define DEFAULT_MODEL_COEFF 1./100.  /* default model coefficient
                                         for mem filtering */

#define DEF_EPS_EVENT_FILTERING 1e-5
#define DEF_WINDOW_SIZE_FILTERING 5
#define DEF_POSITIV_CONSTRAINT True
#define DEF_MAX_CONSTRAINT False
#define DEX_MAX_IMA 255

void mr_nmfilter(Ifloat &Imag, Ifloat &Result, MRNoiseModel &NoiseModel, 
                 int Max_Iter=DEFAULT_MAX_ITER_FILTER, 
                 float Epsilon=DEFAULT_EPSILON_FILTERING, 
                 Bool Sup_Set=False, 
                 type_border Border=DEFAULT_BORDER,
                 Bool PositivIma=True, 
                 Bool KillLastScale=False);

void mr_grad_adj_filter(Ifloat &Imag, Ifloat &Result, MRNoiseModel &NoiseModel,
                     int Max_Iter, type_border Border, Bool KillLastScale=False);

void mr_sfilter(Ifloat &Imag, Ifloat &Result, StatRayleigh & NoiseModel, 
                fltarray & Tab_Epsilon_Threshold,
		float Epsilon_Ratio, int MaxIter, float Sigma_Converge);

inline float linear_weight(float Coef, float Sigma, float Nsigma)
{
   float w = (ABS(Coef) / Sigma -  Nsigma/2.) / (float) Nsigma;
   if (w < 0) w = 0;
   else if (w > 1) w = 1.;
   return w;
}
		
class MRFiltering {  
  void init_param(type_noise TypeNoise, type_filter FilterMethod)
   {
      StatNoise = TypeNoise;
      Max_Iter = DEFAULT_MAX_ITER_FILTER;
      Filter = FilterMethod;
      Epsilon = DEFAULT_EPSILON_FILTERING;
      Border = DEFAULT_BORDER;
      PositivIma = DEF_POSITIV_CONSTRAINT;
      MaxIma = DEF_MAX_CONSTRAINT;
      KillLastScale = False;
      Sup_Set = False;
      Cvg_Coef = 1.;
      IsDataModelled = False;
      UseFlatField = False;
      BackgroundImage = False;
      RegulParam = 0.;
      NoiseModel = NULL;
      Verbose = False;
      UseAdjointRec = True;
      if (Filter == FILTER_ITER_THRESHOLD)  Sup_Set = True;
       
      if (StatNoise == NOISE_EVENT_POISSON) 
      {
         Sup_Set = True; 
         Border = I_MIRROR;
         PositivIma = True;
	 MaxIma = False;
         Epsilon = DEF_EPS_EVENT_FILTERING;
      }
      
      if (StatNoise == NOISE_POISSON) PositivIma = True;
   }      
   MRNoiseModel *NoiseModel;
   
   // filtering by the multiresolution support and optimal step computation
   void grad_adj_filter(MultiResol &MR_Data, Ifloat &Result);
   
   // iterative filtering, minimizing the total variation under constraint
   // on the wavelet coefficient.
   void mr_tv_filter (Ifloat &Imag, MultiResol & MR_Data, Ifloat &Result);
   
   // iterative filtering, minimizing the l_1 norm of the wavelet coefficients
   // under constraint on the wavelet coefficient.
   void mr_wc_filter (Ifloat &Imag, MultiResol & MR_Data, Ifloat &Result);


   // filtering by hierarchical filtering (MR wiener and thresholding).                                
   void hierar_filter(MultiResol & MR_Data, Ifloat &Result);
   
   // filtering by MRwiener, and donoho filtering
   void mr_direct_filter(MultiResol & MR_Data, Ifloat &Result);
  
   public:
   // iterative filtering by the multiresolution support
   void mr_support_filter(Ifloat &Imag, MultiResol & MR_Data, 
                                 Ifloat &Result);   

   // suppress the last scale
   void kill_last_scale(MultiResol & MR_Data);
   
   Bool iter_method()
   {
      Bool Val = False;
      switch(Filter)
      {
         case FILTER_ITER_THRESHOLD:
         case FILTER_ITER_ADJOINT_COEFF:
	 case FILTER_WAVELET_CONSTRAINT:
         case FILTER_TV_CONSTRAINT:
                   Val = True;break;
         default: Val = False;break;
      }
      return Val;
   }
   
   int Max_Iter;        // Maximum number of iterations, for iterative filtering
   type_filter Filter;  // Filtering method
   Bool KillLastScale;  // Kill last scale during the filtering
   type_border Border;  // type of border
   float Epsilon;       // Convergence parameter
   Bool Sup_Set;        // reinitialize the support at each iteration
   Bool PositivIma;     // Apply the positivity constraint
   Bool MaxIma;         // Apply the max constraint
   float Cvg_Coef;      // iterative convergence parameter
   Bool IsDataModelled; // True is the class noise model is already initialized,
                        // then the command ModelData.model(Imag) is not 
                        // performed.
   float RegulParam;    // Regularizarion parameter; 
   type_noise StatNoise; // noise statistic  
   Bool Verbose;         // if set, print informations on stdout
   Bool UseAdjointRec;   // When true (which is the default value)
                         // the adjoint reconstruction is used for the
                         // reconstruction instead of a sinple reconstruction
                         // (for A trou algo and pyramidal methods)
   Bool BackgroundImage; // If true, a background image is given
   Ifloat BackgroundData; // Background image
   Ifloat MaskIma;
   Bool MissingData;     // It true, the input image is masked, and MaskIma = 1/0, 1 for valid pixels.
   
   MRFiltering(MRNoiseModel &ModelData, type_filter FilterMethod)
   {       
      init_param(ModelData.which_noise(), FilterMethod);
      NoiseModel  = &ModelData;
   }
   MRFiltering(type_noise TypeNoise, type_filter FilterMethod)
                            { init_param(TypeNoise, FilterMethod); }
   MRFiltering(type_filter FilterMethod)
                            { init_param(NOISE_GAUSSIAN, FilterMethod); }                          
   void filter(Ifloat &Imag, Ifloat &Result,Bool ModifiedATWT=False);
   
   Bool UseFlatField;  // If true, the input image must be flatfielded. 
   Ifloat FlatField;   // Flat field image
   
   ~MRFiltering(){};
};
   
#endif
