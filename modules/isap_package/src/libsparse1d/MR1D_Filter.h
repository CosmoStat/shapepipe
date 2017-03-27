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
**    Date:  96/05/02
**    
**    File:  MR1D_Filter.h
**
**    Modification history:
**           25-OCT-1996 R Gastaud add mr1d_tab_noise
**
*******************************************************************************
**
**    DESCRIPTION 
**    -----------    
**
**    this module contains the routines for 1D filtering
**
**************************************************************************/ 

#ifndef __MR1D_FIL_1D__
#define __MR1D_FIL_1D__

#include "MR1D_Sigma.h"
#include "MR1D_NoiseModel.h"


// **************************************************************************

#define NBR_1D_FILTERING 4

enum type_1d_filter {FIL_1D_THRESHOLD, 
                     FIL_1D_SOFT_THRESHOLD,
		     FIL_1D_ITER_THRESHOLD,
                     FIL_1D_TV_CONSTRAINT};

inline const char * String1DFilter (type_1d_filter type)
{
    switch (type)
    {
        case FIL_1D_TV_CONSTRAINT:
              return ("Total Variation + Wavelet Constraint");break;
        case FIL_1D_THRESHOLD:
              return ("Multiscale Hard Thresholding");break;
        case FIL_1D_SOFT_THRESHOLD: 
              return ("Multiscale Soft Thresholding");break;
        case FIL_1D_ITER_THRESHOLD: 
              return ("Iterative Multiscale Thresholding");break;
        default: return ("Error: bad type of filtering"); break;
    }
   return ("Error: bad type of filtering");
}

inline void filter_1d_usage(type_1d_filter Filter)
{
    fprintf(OUTMAN, "         [-f type_of_filtering]\n");
    for (int i = 0; i < NBR_1D_FILTERING; i++)
    fprintf(OUTMAN, "              %d: %s \n",i+1,String1DFilter((type_1d_filter)i));
    fprintf(OUTMAN, "              default is %s.\n", String1DFilter(Filter));
}

inline void filterp_usage(type_1d_filter Filter)
{
    fprintf(OUTMAN, "         [-f type_of_filtering]\n");
    fprintf(OUTMAN, "              %d: %s \n",1,String1DFilter(FIL_1D_ITER_THRESHOLD));
    fprintf(OUTMAN, "              %d: %s \n",2,String1DFilter(FIL_1D_THRESHOLD));

    fprintf(OUTMAN, "              default is %s.\n", String1DFilter((type_1d_filter) Filter));
}

#define DEF_EPS 1e-6
#define DEF_ITER 10

#define DEF_EPS_EVENT_FIL_1D 1e-5
#define DEF_WINDOW_SIZE_FIL_1D 5
#define DEF_1D_POSITIV_CONSTRAINT True
#define DEF_1D_MAX_CONSTRAINT False
#define DEF_1D_MAX_SIGNAL 255
#define DEF_EPS_EVENT_1D_FILTERING 1e-3

class MR1DFiltering {  
  void init_param(type_noise TypeNoise, type_1d_filter FilterMethod)
   {
      SmoothWindowSize = DEF_WINDOW_SIZE_FIL_1D;
      StatNoise = TypeNoise;
      Max_Iter = DEF_ITER;
      Filter = FilterMethod;
      Epsilon = DEF_EPS_EVENT_FIL_1D;
      Border = DEFAULT_BORDER;
      PositivIma = DEF_1D_POSITIV_CONSTRAINT;
      MaxIma = DEF_1D_MAX_CONSTRAINT;
      KillLastScale = False;
      Sup_Set = False;
      Cvg_Coef = 1.;
      IsDataModelled = False;
      RegulParam = 0.;
      NoiseModel = NULL;
      Verbose = False;
      ShrinkParam=1;
      UseAdjointRec = True;
      WriteFilterMR = False;
      if (Filter == FIL_1D_ITER_THRESHOLD)  Sup_Set = True;
      NameSupport = NULL;

      if (StatNoise == NOISE_EVENT_POISSON) 
      {
         Sup_Set = True; 
         Border = I_MIRROR;
         PositivIma = True;
	 MaxIma = False;
         Epsilon = DEF_EPS_EVENT_1D_FILTERING;
      }
      
      if (StatNoise == NOISE_POISSON) PositivIma = True;
   }      
   MR1DNoiseModel *NoiseModel;
   
   // filtering by MRwiener, and donoho filtering
   void mr_direct_filter(MR_1D & MR_Data, fltarray &Result);
   
   // simple filtering (mediane, bspline, ...)
   void im_filter(fltarray &Imag, fltarray &Result);
   
   // multiresolution filtering
   void mr_filter(fltarray &Imag, fltarray &Result);
  
   // estimate the threshold estimate in a given band of the wavelet transform
   float sure_estimation(MR_1D &MR_Data);
   float multi_sure_estimation(MR_1D & MR_Data, int b);
  
   void threshold(fltarray &Data, float T=0.);
   void MaxThreshold(fltarray &Data, float T=255.);

   public:
   // iterative filtering by the multiresolution support
   void mr_support_filter(fltarray &Imag, MR_1D & MR_Data, 
                                 fltarray &Result);   

   // suppress the last scale
   void kill_last_scale(MR_1D & MR_Data);
   
   Bool iter_method()
   {
      Bool Val = False;
      switch(Filter)
      {
         case FIL_1D_ITER_THRESHOLD:
                  Val = True;break;
         default: Val = False;break;
      }
      return Val;
   }
   
   int Max_Iter;        // Maximum number of iterations, for iterative filtering
   type_1d_filter Filter;  // Filtering method
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
   int SmoothWindowSize; // window size for smoothing
   Bool Verbose;         // if set, print informations on stdout
   int ShrinkParam;      // Donoho shrinkage parameter
   Bool UseAdjointRec;   // When true (which is the default value)
                         // the adjoint reconstruction is used for the
                         // reconstruction instead of a sinple reconstruction
                         // (for A trou algo and pyramidal methods)
   Bool WriteFilterMR;   // Write the multiresolution support to the disk
   char *NameSupport;    // Name of the multiresolution support

   MR1DFiltering(MR1DNoiseModel &ModelData, type_1d_filter FilterMethod)
   {       
      init_param(ModelData.which_noise(), FilterMethod);
      NoiseModel  = &ModelData;
   }
   MR1DFiltering(type_noise TypeNoise, type_1d_filter FilterMethod)
                            { init_param(TypeNoise, FilterMethod); }
   MR1DFiltering(type_1d_filter FilterMethod)
                            { init_param(NOISE_GAUSSIAN, FilterMethod); }                          
   void filter(fltarray &Imag, fltarray &Result);
   ~MR1DFiltering(){};
};


// **************************************************************************

void mr1d_create_support (MR_1D &MR_Data, float Noise_Ima, 
                             float N_Sigma=DEF_NSIG, float N_Sigma0=DEF_NSIG0);
void mr1d_support_threshold (MR_1D &MR_Data, float Noise_Ima, 
                             float N_Sigma=DEF_NSIG, float N_Sigma0=DEF_NSIG0);
void mr1d_support_set(MR_1D &MR_Data, float Noise_Ima, 
                             float N_Sigma=DEF_NSIG, float N_Sigma0=DEF_NSIG0);

void mr1d_iter_threshold (fltarray &Imag, fltarray &Result, float &Noise_Ima,
                   type_trans_1d Transform, int Nbr_Plan, 
                   float N_Sigma, float Epsilon, int Max_Iter);
             
void mr1d_filter (fltarray &Imag, fltarray &Result, float &Sigma_Noise,
                     type_trans_1d Transform, int Nbr_Plan,
                     type_noise Stat_Noise, float N_Sigma, 
                     float Epsilon, int Max_Iter, 
                     Bool WriteFilterMr=False);

#endif
