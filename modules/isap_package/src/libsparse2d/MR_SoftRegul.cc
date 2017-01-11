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
#include "MR_Obj.h"
#include "MR_SoftRegul.h"

/****************************************************************************/

void regul_tv_haarwt(Ifloat &Ima, Ifloat &Rec, float LambdaTV, int NbrScaleTV)
{
   int Nl = Ima.nl();
   int Nc = Ima.nc();
   int b,i,j;
   FilterAnaSynt *TV_SelectFilter; // Selected filter bank
   SubBandFilter *TV_SB1D;  // Pointer to the fiter bank decomposition class
   HALF_DECIMATED_2D_WT *TV_WT; // Pointer to the wavelet transform class
   Ifloat *TV_WT_Trans;           // Wavelet bases decomposition
   int NbrUndecimatedScale=-1;
   TV_SelectFilter = new FilterAnaSynt(F_HAAR);
   TV_SB1D = new SubBandFilter(*TV_SelectFilter, NORM_L2);
   TV_WT = new HALF_DECIMATED_2D_WT(*TV_SB1D);
   int TV_NbrBand = TV_WT->alloc(TV_WT_Trans, Nl, Nc, NbrScaleTV, NbrUndecimatedScale);

   TV_WT->transform(Ima, TV_WT_Trans, NbrScaleTV, NbrUndecimatedScale);
   
   float Tau = LambdaTV;
   for (b=0; b < TV_NbrBand-1; b++)
   {
      int s = b / 3;
      if ((s != 0) && ((b % 3) == 0)) Tau = LambdaTV / pow(sqrt(2.), s); 
      // cout << "  Band " << b+1 << " TV threshold = " << Tau << endl;
      for (i=0; i <  (TV_WT_Trans[b]).nl(); i++)
      for (j=0; j <  (TV_WT_Trans[b]).nc(); j++)
      {
        (TV_WT_Trans[b])(i,j) = soft_threshold((TV_WT_Trans[b])(i,j), Tau);
      }
   }
   TV_WT->recons(TV_WT_Trans, Rec, NbrScaleTV, NbrUndecimatedScale);   

   TV_WT->free(TV_WT_Trans, NbrScaleTV);
   if (TV_SB1D != NULL) delete TV_SB1D;
   if (TV_SB1D != NULL) delete TV_SelectFilter;   
}

/****************************************************************************/

static int max_scale_number (int N)
{
    int ScaleMax;
    ScaleMax  = iround((float)log((float) (N / 4. * 3.) / log(2.)));
    return ScaleMax;
}

void MR_Regul::im_soft_non_ortho_threshold(Ifloat &Ima, Ifloat &Rec, float Lambda, float NoiseLevel)
{
   int Nl = Ima.nl();
   int Nc = Ima.nc();
   int b,i,j; 
   if (NbrScale < 2) NbrScale =  max_scale_number(MIN(Nl,Nc));

   MultiResol MR_Data;
   FilterAnaSynt FAS;
   FilterAnaSynt *PtrFAS = NULL;
   MR_Data.alloc (Nl, Nc, NbrScale, TO_UNDECIMATED_NON_ORTHO, PtrFAS, NORM_L1, -1, U_B3SPLINE_2);
   MR_Data.ModifiedATWT = True;
   MR_Data.transform(Ima);
   float Tau = Lambda;
   for (b=0; b < MR_Data.nbr_band()-1; b++)
   {
      int s = b / 3;
      if (ExpDecreasingLambda == True)
      {
         if ((s != 0) && ((b % 3) == 0)) Tau = Lambda / pow(sqrt(2.), s); 
         // cout << "  Band " << b+1 << " TV threshold = " << Tau << endl;
      }
      for (i=0; i < Nl; i++)
      for (j=0; j < Nc; j++)
      {
        float Coef = MR_Data(b,i,j);
	float SoftLevel = Tau*MR_Data.band_norm(b);
        if (NoiseLevel != 0)
	{
	    float SNR = ABS(Coef) / (NoiseLevel*MR_Data.band_norm(b)) -1;
	    if (SNR > 0) SoftLevel *= erfc(SNR);
	}
        MR_Data(b,i,j) = soft_threshold(Coef, SoftLevel);
      }
   }
   MR_Data.recons(Rec);   
}

void MR_Regul::im_soft_iwt_threshold(Ifloat &Ima, Ifloat &Rec, float Lambda, float NoiseLevel)
{
   int Nl = Ima.nl();
   int Nc = Ima.nc();
   int b,i,j; 
   if (NbrScale < 2) NbrScale =  max_scale_number(MIN(Nl,Nc));

   MultiResol MR_Data;
   FilterAnaSynt FAS;
   FilterAnaSynt *PtrFAS = NULL;
   MR_Data.alloc (Nl, Nc, NbrScale, TO_PAVE_BSPLINE, PtrFAS, NORM_L1, -1, U_B3SPLINE_2);
   MR_Data.ModifiedATWT = True;
   MR_Data.transform(Ima);
   float Tau = Lambda;
   for (b=0; b < MR_Data.nbr_band()-1; b++)
   {
      int s = b / 3;
      if (ExpDecreasingLambda == True)
      {
         if ((s != 0) && ((b % 3) == 0)) Tau = Lambda / pow(sqrt(2.), s); 
         // cout << "  Band " << b+1 << " TV threshold = " << Tau << endl;
      }
      for (i=0; i < Nl; i++)
      for (j=0; j < Nc; j++)
      {
        float Coef = MR_Data(b,i,j);
	float SoftLevel = Tau*MR_Data.band_norm(b);
        if (NoiseLevel != 0)
	{
	    float SNR = ABS(Coef) / (NoiseLevel*MR_Data.band_norm(b)) -1;
	    if (SNR > 0) SoftLevel *= erfc(SNR);
	}
        MR_Data(b,i,j) = soft_threshold(Coef, SoftLevel);
      }
   }
   MR_Data.recons(Rec);   
}


/****************************************************************************/


void MR_Regul::im_soft_threshold(Ifloat &Ima, Ifloat &Rec, float Lambda, float NoiseLevel)
{
   int Nl = Ima.nl();
   int Nc = Ima.nc();
   int b,i,j; 
   if (NbrScale < 2) NbrScale =  max_scale_number(MIN(Nl,Nc));

   FilterAnaSynt *TV_SelectFilter; // Selected filter bank
   SubBandFilter *TV_SB1D;  // Pointer to the fiter bank decomposition class
   HALF_DECIMATED_2D_WT *TV_WT; // Pointer to the wavelet transform class
   Ifloat *TV_WT_Trans;           // Wavelet bases decomposition
   TV_SelectFilter = new FilterAnaSynt(TypeFilter);
   TV_SB1D = new SubBandFilter(*TV_SelectFilter, NORM_L2);
   TV_WT = new HALF_DECIMATED_2D_WT(*TV_SB1D);
   int TV_NbrBand = TV_WT->alloc(TV_WT_Trans, Nl, Nc, NbrScale, NbrUndecimatedScale);

   TV_WT->transform(Ima, TV_WT_Trans, NbrScale, NbrUndecimatedScale);
   
   float Tau = Lambda;
   for (b=0; b < TV_NbrBand-1; b++)
   {
      int s = b / 3;
      if (ExpDecreasingLambda == True)
      {
         if ((s != 0) && ((b % 3) == 0)) Tau = Lambda / pow(sqrt(2.), s); 
         // cout << "  Band " << b+1 << " TV threshold = " << Tau << endl;
      }
      for (i=0; i <  (TV_WT_Trans[b]).nl(); i++)
      for (j=0; j <  (TV_WT_Trans[b]).nc(); j++)
      {
        float Coef = (TV_WT_Trans[b])(i,j);
	float SoftLevel = Tau;
        if (NoiseLevel != 0)
	{
	    float SNR = ABS(Coef) / NoiseLevel -1;
	    if (SNR > 0) SoftLevel *= erfc(SNR);
	}
        (TV_WT_Trans[b])(i,j) = soft_threshold(Coef, SoftLevel);
      }
   }
   TV_WT->recons(TV_WT_Trans, Rec, NbrScale, NbrUndecimatedScale);   
   TV_WT->free(TV_WT_Trans, NbrScale);
   if (TV_SB1D != NULL) delete TV_SB1D;
   if (TV_SB1D != NULL) delete TV_SelectFilter;      
}

/****************************************************************************/

void MR_Regul::ima_regul(Ifloat &Obj, Ifloat &Grad, float Lambda)
{
   int i,j;
   int Nl = Obj.nl();
   int Nc = Obj.nc();
   im_soft_threshold(Obj, Grad, Lambda);
   for (i=0;i < Nl; i++)
   for (j=0;j < Nc; j++)
         Grad(i,j) = (Obj(i,j) - Grad(i,j));  
}

/****************************************************************************/

void MR_Regul::obj_regul(Ifloat &Obj, Ifloat &Grad, float  Lambda)
{
      int i,j;
      int Nl = Obj.nl();
      int Nc = Obj.nc();
      Ifloat ImaAux(Nl,Nc,"aux");
      im_soft_threshold(Obj, ImaAux, Lambda);
      for (i=0;i < Nl; i++)
      for (j=0;j < Nc; j++)
          Grad(i,j) += (ImaAux(i,j) - Obj(i,j));
}
   
/****************************************************************************/
