/******************************************************************************
**                   Copyright (C) 1996 by CEA
*******************************************************************************
**
**    UNIT
**
**    Version: 3.3
**
**    Author: Jean-Luc Starck
**
**    Date:  96/07/08
**    
**    File:  MR_Sigma.cc
**
*******************************************************************************
**
**    DESCRIPTION  
**    -----------  
**    
*******************************************************************************
**
** float detect_noise_from_support(Ifloat &Imag, int Nbr_Plan)
**
** Estimates the noise from the multiresolution support
**
*******************************************************************************
**
** float detect_noise_from_bspline(const Ifloat &Image)
** 
** returns an estimation of the standard deviation of the noise in Image.
** This estimation is computed from the difference between Image and 
** Image smoothed by the bspline filter
**
*******************************************************************************
**
** float get_readout_noise(Ifloat &Imag, float Gain, type_method Method)
** 
** Find the readout noise in a image which contains poisson noise + gaussian
** noise (readout noise)
**
**
*******************************************************************************
**
** float get_noise (Ifloat &Imag, type_method Method, int Nbr_Plan)
** 
** Find the standard deviation of the gaussian noise in a image
** by the method specified in Method
**
*******************************************************************************
**
** float get_gain(Ifloat &Imag, float Readout, type_method Method,
**                float GainMin, float GainMax)
**
** Find the gain in a image which contains poisson noise + gaussian
** noise (readout noise)
**
******************************************************************************/

// static char sccsid[] = "@(#)MR_Sigma.cc 3.3 96/07/08 CEA 1996 @(#)";

#include "MR_Obj.h"
#include "MR_Noise.h"
#include "IM_Sigma.h"
#include "MR_NoiseModel.h"
#include "MR_Sigma.h"

/***************************************************************************/

float detect_noise_from_support(Ifloat &Imag, MultiResol &MR_Data, MRNoiseModel &ModelData, int Niter)
{
  int i,j,s,k=0;
  int nb;
  int Nl = ModelData.nl();
  int Nc = ModelData.nc();
  float  ValIma;
  double std, Noise_Ima=0., OldEst, Cvg;
  // int MaxIter=15;
  int Nbr_Plan = MR_Data.nbr_scale();
  unsigned char *Mask;
  float NSigma[MAX_SCALE]; 

  Mask = new unsigned char [Nl*Nc];
  for (s=0; s < MAX_SCALE; s++)
  {
      NSigma[s] = ModelData.NSigma[s];
      ModelData.NSigma[s] = 3.;
  }
 do
 {
     /* creation of the support */
     
     OldEst = Noise_Ima;
     ModelData.SigmaNoise = Noise_Ima;
     ModelData.set_sigma(Imag, MR_Data);
     ModelData.set_support(MR_Data);

     s=0;
     for (i=0;i<Nl;i++)
     for (j=0;j<Nc;j++) 
     {
        s = 0;
        while ( (s < Nbr_Plan-1) && (ModelData(s,i,j) == False) ) s++;
        if (s == Nbr_Plan-1) Mask[i*Nc+j] = 0;
        else Mask[i*Nc+j] = 1;
     }
 
     /* calculation of noise standard deviation with the mask*/
     std = 0.0;
     nb=0;
     for (i=0;i<Nl;i++)
     for (j=0;j<Nc;j++) 
       if ( Mask[i*Nc+j] == 0)
       {
	   ValIma = Imag(i,j) - MR_Data(Nbr_Plan-1,i,j);
	   std = std + ValIma*ValIma;
	   nb ++;
       }

     if (nb == 0)  std = 0.0;
     else std = sqrt(std/(float) nb) / 0.973703; 
     Noise_Ima = std;
     k++;
     Cvg = ABS(OldEst-Noise_Ima)/Noise_Ima;
     cout << "Iter " << k << " Noise = " << Noise_Ima << "  cvg = " << Cvg << endl;
 } while ( (k <= Niter) && (Cvg > 1.e-5));


  for (s=0; s < MAX_SCALE; s++)
  {
      ModelData.NSigma[s] = NSigma[s];
  }

  delete [] Mask;
  return (float) std;
}

/***************************************************************************/

float detect_noise_from_support (Ifloat &Imag, int Niter, int Nbr_Plan)
{
  int Nl = Imag.nl();
  int Nc = Imag.nc();
  float SigmaRet;
  type_noise Stat_Noise = NOISE_GAUSSIAN;
  type_transform Transform = TO_PAVE_BSPLINE;

  MRNoiseModel ModelData(Stat_Noise, Nl, Nc, Nbr_Plan, Transform);
  MultiResol MR_Data(Nl,Nc,Nbr_Plan, Transform,"mr_sigma");
  ModelData.model(Imag, MR_Data);

  SigmaRet = detect_noise_from_support(Imag, MR_Data, ModelData, Niter);

  return SigmaRet;
}


/***************************************************************************/

float get_noise (Ifloat &Imag, type_method Method, int Niter, int Nbr_Plan)
{
  float SigmaRet=0.;

  switch(Method)
  {
     case METHOD_SIGMA_CLIPPING:
		SigmaRet = detect_noise_sigma(Imag);
		break;	
     case METHOD_MEDIANE:
		SigmaRet = detect_noise_from_med(Imag);
		break; 
     case METHOD_BSPLINE:
		SigmaRet = detect_noise_from_bspline (Imag);
		break;
     case METHOD_MR_SUPPORT:
		SigmaRet = detect_noise_from_support(Imag, Niter, Nbr_Plan);
		break;
     case METHOD_BLOCK:
		SigmaRet = detect_noise_from_block(Imag);
		break;
     case METHOD_MAD:
		SigmaRet = detect_noise_from_mad(Imag);
		break;
     default:
		cerr << "Error: Unknown sigma detection method ... " << endl;
		exit(-1);
		break;
  }
  return SigmaRet;
}

/***************************************************************************/

float get_gain(Ifloat &Imag, float Readout, type_method Method, float GainMin, float GainMax)
{
	int i,j,Nc,Nl;
	float Val, Gain;
	float Gain_min=GainMin;
	float Gain_max=GainMax;
	float EstNoise;
	float Ro2 = Readout*Readout;

	Nc = Imag.nc();
	Nl = Imag.nl();
	Ifloat Transfo(Nc,Nl,"mr_sigma");

	while (Gain_max - Gain_min > 1e-1) 
	{
	   Gain = 0.5*(Gain_max + Gain_min);

	   for (i=0;i<Nc;i++)
	   for (j=0;j<Nl;j++) 
	   {
	     Val = Gain*Imag(i,j) + 3*Gain*Gain/8 + Ro2;
	     if ( Val< 0.0 ) Transfo(i,j) = 0.0;
	     else Transfo(i,j) = 2*sqrt(Val)/Gain;
	   }
	   EstNoise = get_noise(Transfo, Method);

	   if (EstNoise > 1.0) Gain_min = Gain;
	   else Gain_max = Gain;
	   cout << "Gain = " << Gain << ", Gain_min = " << Gain_min << ", Gain_max = " << Gain_max << endl; 
	}
	return 0.5*(Gain_max + Gain_min);
}

/***************************************************************************/

float get_readout_noise(Ifloat &Imag, float Gain, type_method Method, int NbrPlan)
{
   int i,j,Nc,Nl;
   float mean, transfo_std, std, std_min=0.0, std_max;
   float val;

   Nc = Imag.nc();
   Nl = Imag.nl();
   Ifloat Transfo(Nc,Nl,"mr_sigma");

   mean = (float) average(Imag);
   std_max = (float) sigma(Imag);
   cout << "Mean = " << mean << "  Sigma = " << std_max << endl;

   while (std_max - std_min > 1e-1)
   {
      std = 0.5*(std_max + std_min);
      for (i=0;i<Nc;i++)
      for (j=0;j<Nl;j++) 
      {
         val = Gain*Imag(i,j) + 3*Gain*Gain/8 + std*std;
         Transfo(i,j) = (val < 0.) ? 0.: 2*sqrt(val)/Gain;
      }
      transfo_std = get_noise(Transfo,Method,DEFAULT_MAX_ITER_METHOD,NbrPlan);

      if ( transfo_std > 1.0) std_min = std;
      else  std_max = std;
      cout << "transfo_std = " << transfo_std << ", std_min = " << std_min << ", std_max = " << std_max << endl;
   }
	
   return 0.5*(std_max + std_min);
}

/*****************************************************************************/
