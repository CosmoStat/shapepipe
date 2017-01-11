/*******************************************************************************
**                   Copyright (C) 2002 by CEA
*******************************************************************************
**
**    UNIT
**
**    Version: 1.0
**
**    Author: Jean-Luc Starck
**
**    Date:  18/01/02 
**    
**    File:  MR_Contrast.cc
**
*******************************************************************************
**
**    DESCRIPTION  multiresolution constrast enhancement
**    -----------  
**                 
**
******************************************************************************/

#include "IM_Obj.h"
#include "MR_Obj.h"
#include "MR_Contrast.h"

/****************************************************************************/

const char * StringContrast ( type_contrast  type)
{ 
    switch (type)
    {
        case  CT_HISTO:
	      return ("Histogram Equalization"); 
        case  CT_RET: 
              return ("Retinex Method"); 
        case  CT_MRET: 
              return ("Multiscale Retinex Method"); 
        case  CT_ATROU_RET:
	      return ("Multiscale a Trous Retinex Method");
        case  CT_WT:
	      return ("Wavelet Coefficients Enhancement");
        case  CT_CUR:
	      return ("Curvelet Coefficients Enhancement");
        case  CT_ATROU_LOG:
	      return ("Wavelet-Log function: f(w) = log(|w|+L)");
	case  CT_ATROU_LOG_SGN:
	      return ("Wavelet-Log function: f(w) = sgn(w).log(|w|+L)");
	case  CT_CLIPPING:
	      return ("K-Sigma clipping.");
        case  CT_LAPLACIAN:
	      return ("Add the Laplacian to the data.");
	default: 
	      return ("Undefined contrast enhancement method");
     } 
}

/****************************************************************************/
  
static float get_laplacian(Ifloat &Obj, int i, int j)
{
  type_border type=I_CONT;
  float Val =  Obj(i,j) - 0.25*( Obj(i+1,j,type)
                    +  Obj(i-1,j,type) +  Obj(i,j+1,type) +  Obj(i,j-1,type));
  return Val;
}   

/****************************************************************************/

void Contrast::laplacian(Ifloat &Data)
{ 
   int i,j,Nc = Data.nc();
   int Nl = Data.nl(); 
   Ifloat  Buff(Nl,Nc, "buff");

   Buff = Data;
   for (i = 0; i < Nl; i ++)
   for (j = 0; j < Nc; j ++)
         Data(i,j) +=  LaplacianParam*get_laplacian(Buff,i,j);
}

/****************************************************************************/

void Contrast::retinex(Ifloat &Data)
{ 
   int i,j;
   int Nc = Data.nc();
   int Nl = Data.nl();
   float Sigma = 80.;
   Ifloat Gaussian(Nl,Nc,"gauss");
   FFTN_2D CFFT;
   Ifloat  Buff(Nl,Nc, "buff");
 
   Gaussian.init();
   make_gaussian(Gaussian, (float) (Sigma/sqrt(2.)));
   CFFT.convolve(Data, Gaussian, Buff);
   for (i = 0; i < Nl; i++)
   for (j = 0; j < Nc; j++) 
   {
      if ((Data(i,j) > 0) && (Buff(i,j) > 0))
	        Data(i,j) = log(Data(i,j)) - log( Buff(i,j));
   }
}

/*************************************************************************/
  
void Contrast::wt_enhance(Ifloat &Data, int Nbr_Plan)
{
   int Nl = Data.nl();
   int Nc = Data.nc();
   int k,i,j;
   // float MeanU, MeanV, MeanTU, MeanTV;
   float TabMin, TabMax;
   float TabTMin, TabTMax;    
   type_transform Transform = TO_DIADIC_MALLAT;
   MultiResol MR_Data(Nl, Nc, Nbr_Plan, Transform, "diadic");
   double Q100 = pow((double) 100., Contrast_Q_Param);

    // Non linear mapping of the luminance
    if (Contrast_Q_Param != 0)
    for (i=0; i < Nl; i++)
    for (j=0; j < Nc; j++)
 	 Data(i,j) = (float) (pow((double) Data(i,j), 1.-Contrast_Q_Param)*Q100);
 
    // Wavelet Transform
    float Gamma;
    //TabMin = min(Data);
    //TabMax = max(Data);
    if (Verbose == True) cout <<"Wavelet Transform" <<endl;
    MR_Data.transform(Data);
 
    if (Verbose == True) cout <<"Contrast function calculation " <<endl;
    
    // Wavelet coefficients correction
    for (k=0; k < MR_Data.nbr_scale()-1; k++) 
    {   
        float Norm = MR_Data.band_norm(2*k);
	float Noise = Noise_Ima; // Noise_Ima*Norm;
	Contrast_C_Param = Noise;
 	if (Verbose == True) cout <<"Scale " << k+1 << endl;
	for (i=0; i < Nl; i++)
        for (j=0; j < Nc; j++) 
	{
	   // We need to divive the wavelet coefficients because of 
	   // the normalization of the g filter. Indeed, in Velde's paper
	   // g = (-1/2,1/2) and in our program g = (2, -2).
	   // G-normarmalization has no effect as long as we take into
	   // account the normalization value in the contrast function.
 	   Gamma = sqrt(sqr(MR_Data(2*k,i,j)/4.)+ sqr(MR_Data(2*k+1,i,j)/4.));
  	   Gamma = contrast_function(Gamma);
 	   MR_Data(2*k,i,j) *= Gamma;
	   MR_Data(2*k+1,i,j) *= Gamma;
        }
     }
 
     MR_Data.recons(Data);
     
//     TabTMin = min(Data);
//     TabTMax = max(Data);
    
//     cout << "Min " << TabMin << " Max = " << TabMax << endl;
//     cout << "MinR " << TabTMin << " MaxR = " << TabTMax << endl;
//     INFO(Data,"REC");
//     for (i=0; i < Nl; i++)
//     for (j=0; j < Nc; j++)
//     {
//        float Scale = (TabMax - TabMin) / (TabTMax - TabTMin);
//        Data(i,j) =  (Data(i,j) - TabTMin) * Scale + TabMin;
//     }
//     INFO(Data,"REC1");
    // inverse mapping of the Luminance
    if (Contrast_Q_Param != 0)
    {
       if (Verbose == True) 
           cout << "inverse mapping of the Luminance " <<  endl;

       double  Q1 = 1. / (1-Contrast_Q_Param);
       Q100 = pow((double) 100., -Contrast_Q_Param*Q1);
       for (i=0; i < Nl; i++)
       for (j=0; j < Nc; j++)
 	  Data(i,j) = (float) (pow((double) Data(i,j), Q1)*Q100);
    }
}

/***************************************************************************/

void Contrast::multiscale_retinex(Ifloat &Data)
{ 
   int i,j,b;
   int Nc = Data.nc();
   int Nl = Data.nl();
   float TabSigma[3] = {15.,80.,250.};
   float G = 192;
   float B =  -30.;
   float Alpha = 125.;
   float Beta = 46.;
   Ifloat Gaussian(Nl,Nc,"gauss");
   Ifloat Ret(Nl,Nc,"ret");
   Ifloat Ci(Nl,Nc,"ci");
   FFTN_2D CFFT;

   // Ci calculation
   for (i = 0; i < Nl; i++)
   for (j = 0; j < Nc; j++) 
      if (Data(i,j) > 0) Ci(i,j) = Beta*(log(Alpha*Data(i,j)));
   Ifloat Buff(Nl,Nc, "buff");
 
   int N = 3;
   for (b = 0; b < N; b++)
   {
      Gaussian.init();
      make_gaussian(Gaussian, (float) (TabSigma[b]/sqrt(2.)));
      CFFT.convolve(Data, Gaussian, Buff);
      for (i = 0; i < Nl; i++)
      for (j = 0; j < Nc; j++) 
      {
            if ((Data(i,j) > 0) && (Buff(i,j) > 0))
 	     Ret(i,j) += 
 	        G*(Ci(i,j)*(log(Data(i,j)) - log( Buff(i,j))) + B) / (float) N;
      }
   }
   Data = Ret;
}
   
/***************************************************************************/

void Contrast::atrou_retinex(Ifloat &Data, int Nbr_Plan)
{
   int i,j,b;
   int Nc = Data.nc();
   int Nl = Data.nl();
   type_transform Transform = TO_PAVE_BSPLINE;
   MultiResol MR_Data(Nl, Nc, Nbr_Plan, Transform, "diadic");
   Ifloat Ret(Nl,Nc,"ret");
   Ifloat Buff(Nl,Nc, "buff");

   MR_Data.transform(Data);
   Buff = MR_Data.band(Nbr_Plan-1);
   for (b = Nbr_Plan-2; b >= 0; b--)
   {
      for (i = 0; i < Nl; i++)
      for (j = 0; j < Nc; j++) 
       {
            if ((Data(i,j) > 0) && (Buff(i,j) > 0))
              Ret(i,j) += (log(Data(i,j)) - log( Buff(i,j))) / (float) (Nbr_Plan-1);
         }
	 Buff += MR_Data.band(b);
   }
   Data = Ret;    
}

/*************************************************************************/

void Contrast::atrou_log(Ifloat &Data, int Nbr_Plan, Bool Sign)
{
   int i,j,b;
   int Nc = Data.nc();
   int Nl = Data.nl();
   type_transform Transform = TO_PAVE_BSPLINE;
   MultiResol MR_Data(Nl, Nc, Nbr_Plan, Transform, "diadic");
   Ifloat Buff(Nl,Nc, "buff");

   MR_Data.transform(Data);
   Data.init();
   for (b = Nbr_Plan-1; b >= 0; b--)
   {
      float Noise = Noise_Ima*MR_Data.band_norm(b);
      if (Noise < FLOAT_EPSILON) Noise = 1.;
      for (i = 0; i < Nl; i++)
      for (j = 0; j < Nc; j++)  
      {
         float Coef = MR_Data(b,i,j)/Noise;
         if (Sign == False) Data(i,j) += log(LogParam+ABS(Coef));
	 else
	 {
	   if (MR_Data(b,i,j) > 0) Data(i,j) += log(LogParam+Coef);
	   else Data(i,j) -= log(LogParam-Coef);
	 }
      }
   }
}

/*************************************************************************/

void Contrast::sature_clipping(Ifloat &Data)
{
   int i,j;
   int Nc = Data.nc();
   int Nl = Data.nl();
   float MeanClip, SigmaClip;
   
   sigma_clip (Data, MeanClip, SigmaClip,  3);
   float MinVal = MeanClip - ClipVal*SigmaClip;
   float MaxVal = MeanClip + ClipVal*SigmaClip;
   for (i = 0; i < Nl; i++)
   for (j = 0; j < Nc; j++) 
   {
      if (Data(i,j) > MaxVal) Data(i,j) = MaxVal;
      else if (Data(i,j) < MinVal ) Data(i,j) = MinVal;
   }
}

/*************************************************************************/

void Contrast::histo_equal(Ifloat &Data, int NbrBin)
{
   int i,j,k,Nl = Data.nl();
   int Nc = Data.nc();
   float Min, Max, Scale;
   Min = Max = Data(0,0);
   
   for (i=0; i< Nl; i++)
   for (j=0; j< Nc; j++)
   {
      if (Data(i,j) < Min) Min = Data(i,j);
      else if (Data(i,j) > Max) Max = Data(i,j);
   }
   if (Max>Min) Scale = (NbrBin-1)/(Max-Min);
   else  Scale = 0;
   intarray Map(NbrBin);
   Map.init();
   for (i=0; i< Nl; i++)
   for (j=0; j< Nc; j++)
   {
      int Ind = (int)((Data(i,j)-Min)*Scale);
      Map(Ind) ++;
   }
   for (k=1; k<NbrBin; k++) Map(k) += Map(k-1);

   for (i=0; i<NbrBin; i++)
         Map(i) = (int)(Min + (Max-Min)*Map(i)/(float)(Nl*Nc));
	 
   for (i=0; i< Nl; i++)
   for (j=0; j< Nc; j++)
   {
      int Ind = (int)((Data(i,j)-Min)*Scale);
      Data(i,j) = Map(Ind);
   }
}

/****************************************************************************/

