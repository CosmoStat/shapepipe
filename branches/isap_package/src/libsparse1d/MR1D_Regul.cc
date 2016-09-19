/*****************************************************************
**
**    UNIT
**
**    Version: 3.1
**
**    Author: Jean-Luc Starck
**
**    Date:  96/05/02 
**    
**    File:  MR_Haar.cc
**
*****************************************************************/ 



#include <stdio.h>
#include <math.h>
#include <string.h>

#include "GlobalInc.h"
#include "MR1D_Obj.h"
#include "MR1D_Regul.h"

  
/****************************************************************************/

static int max_scale_number (int N)
{
    int ScaleMax;
    ScaleMax  = iround((float)log((float) (N / 4. * 3.) / log(2.)));
    return ScaleMax;
}

void RegulSig::mr1d_soft_threshold(fltarray &Sig, fltarray &Rec, float Lambda, float NoiseLevel)
{
   int Nx = Sig.nx();
   int b,i; 
   if (NbrScale < 2) NbrScale =  max_scale_number(Nx);
   FilterAnaSynt *SelectFilter; // Selected filter bank
   SubBandFilter *SB1D;  // Pointer to the fiter bank decomposition class
   fltarray WT_Trans;           // Wavelet bases decomposition
   SelectFilter = new FilterAnaSynt(TypeFilter);
   SB1D = new SubBandFilter(*SelectFilter, NORM_L2);
   PAVE_1D_WT  WT(*SB1D); // Pointer to the wavelet transform class
   
    WT.transform(Sig, WT_Trans,NbrScale);
    float Tau = Lambda;
   for (b=0; b < NbrScale-1; b++)
   {
      int s = b;
      if (ExpDecreasingLambda == True)
      {
         if (s != 0) Tau = Lambda / pow(sqrt(2.), s); 
         // cout << "  Band " << b+1 << " TV threshold = " << Tau << endl;
      }
      for (i=0; i <  Nx; i++)
      {
        float Coef = WT_Trans(i,b);
	float SoftLevel = Tau;
        if (NoiseLevel != 0)
	{
	    float SNR = ABS(Coef) / NoiseLevel -1;
	    if (SNR > 0) SoftLevel *= erfc(SNR);
	}
        WT_Trans(i,b) = soft_threshold(Coef, SoftLevel);
      }
   }
    WT.recons(WT_Trans, Rec, NbrScale);  
    if (SB1D != NULL) delete SB1D;
   if (SB1D != NULL) delete SelectFilter;  
 }

/****************************************************************************/
 
