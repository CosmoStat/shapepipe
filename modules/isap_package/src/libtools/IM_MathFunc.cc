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
**    File:  IM_MathFunc.cc
**
*******************************************************************************
**
** Ifloat & im_noise (float Sigma, int Nl, int Nc)
**
** create a gaussian noise with a standart deviation equal to Sigma
** and return image containing this noise
** Sigma = standard deviation of the noise
** Nl, Nc = image size
**
*******************************************************************************
**
** const Ifloat & im_gaussian (int Nl, int Nc, float Fwhm, int Indi, int Indj)
**
** return an image which contain a gaussian at the position (Indi, Indj) 
** Nl, Nc = image size
** Fwhm = full width at hall maximum
** Indi, Indj = position of the maximum of the gaussian
**
*******************************************************************************
**
** const fltarray randomn(float Sigma, int N)
** 
** generates in a noise in a array
**
*****************************************************************************/
 
// static char sccsid[] = "@(#)IM_MathFunc.cc 3.1 96/05/02 CEA 1994 @(#)";

#include "IM_Obj.h"
#include "Array.h"
#include "IM_Noise.h"
  
/***************************************************************************/

/*const fltarray randomn(int N, float Sigma, unsigned int InitRnd)
{
    int i;
    double x;
    fltarray *Result = new fltarray(N);
  
    if (InitRnd > 0) init_random (InitRnd);
    for (i=0; i < N; i++)
    {
         x =   get_random();
         (*Result)(i) = xerf(x) * Sigma;
    }
    return (*Result);
}*/
void randomn(fltarray& Tab, int N, float Sigma, unsigned int InitRnd)
{
    int i;
    double x;
  
    if (InitRnd > 0) init_random (InitRnd);
    for (i=0; i < N; i++)
    {
         x =   get_random();
         (Tab)(i) = xerf(x) * Sigma;
    }
}


/***************************************************************************/

const Ifloat im_gaussian (int Nl, int Nc, float Fwhm, int Indi, int Indj)
{
    int Dist,Delta_i;
    float sigma,sigma2,Energ = 0.;
    register int i,j;
    Ifloat *Result = new Ifloat(Nl,Nc, "gaussienne");
    int Ci=Indi, Cj=Indj;

    if (Ci < 0) Ci = Nl / 2;
    if (Cj < 0) Cj = Nc / 2;

    sigma = 0.5 * Fwhm / sqrt (2. * log ((double) 2.));
    sigma2 = -2. * sigma * sigma;

    for (i=0; i < Nl; i++)
    {
        Delta_i = (i - Ci) * (i - Ci);
        for (j = 0; j < Nc; j++)
	{
	    Dist = Delta_i + (j - Cj) * (j - Cj);
	    (*Result)(i,j) = exp((double)((float) Dist / sigma2));
	    Energ +=  (*Result)(i,j);
	 }
    } 

    return (*Result);
}


/***************************************************************************/

void make_gaussian(Ifloat & Result, float Sigma, int Indi, int Indj)
{
   int Nl = Result.nl();
   int Nc = Result.nc();
   int Dist,Delta_i;
   float sigma2,Energ = 0.;
   register int i,j;
   int Ci=Indi, Cj=Indj;
    if (Ci < 0) Ci = Nl / 2;
    if (Cj < 0) Cj = Nc / 2;
    sigma2 = -2. * Sigma * Sigma;
    for (i=0; i < Nl; i++)
    {
        Delta_i = (i - Ci) * (i - Ci);
        for (j = 0; j < Nc; j++)
	{
	    Dist = Delta_i + (j - Cj) * (j - Cj);
	    Result(i,j) = exp((double)((float) Dist / sigma2));
	    Energ +=  Result(i,j);
	 }
    }
    for (i=0; i < Nl; i++)
    for (j=0; j < Nc; j++)  Result(i,j) /= Energ ;
   
}

/***************************************************************************/

void make_gaussian_fwhm(Ifloat & Result, float Fwhm, int Indi, int Indj)
{
   float sigma;
   sigma = 0.5 * Fwhm / sqrt (2. * log ((double) 2.));
   make_gaussian(Result, sigma, Indi,Indj);
}

/***************************************************************************/

void im_noise_poisson (Ifloat & NoiseIma, float  Gain, unsigned int InitRnd)
{
  int i;
  int number= (int) InitRnd;

  for (i=0;i< NoiseIma.n_elem();i++) 
        NoiseIma(i) = poidev(NoiseIma(i),&number) * Gain;
}

/***************************************************************************/

void im_noise_gaussian (Ifloat & NoiseIma, float Sigma, unsigned int InitRnd)
{
    int i,j;
    double x;
    int Nl = NoiseIma.nl();
    int Nc = NoiseIma.nc();
     
    if (InitRnd > 0) init_random (InitRnd);
    for (i=0; i < Nl; i++)
    for (j=0; j < Nc; j++)
    {
       x = get_random();
       NoiseIma(i,j) = xerf(x) * Sigma;
    }
}

/***************************************************************************/

const Ifloat im_noise (float Sigma, int Nl, int Nc, unsigned int InitRnd)
/* create a gaussian noise */
{
    Ifloat *Result = new Ifloat(Nl,Nc, "im_noise");
  
    im_noise_gaussian(*Result, Sigma, InitRnd);
    return (*Result);
}

/***************************************************************************/

void im_noise_rayleigh(Ifloat & NoiseIma, int  NbrImage, unsigned int InitRnd)
/* create  noise model
  
      ref: "Digital Image Processing" W.K. PRATT  P.318  
        Rayleigh Noise = image_sum ( 1/(sigma^2) exp(- X^2/(sigma^2)))/NbrImag
        with sigma = 1/alpha
*/
{
   int Nl = NoiseIma.nl();
   int Nc = NoiseIma.nc();
   int i,j,k;
   float sigma2;
   float alpha=1.;
   
   if (InitRnd > 0) init_random (InitRnd);
   sigma2 = 2./(alpha*alpha);
   
   for(i=0; i < Nl;i++)
   for(j=0; j < Nc; j++)
   {
     float V1 = FLOAT_EPSILON;
     float V2 = 1. - FLOAT_EPSILON;
     float Val;
     NoiseIma(i,j)=0.;
     for(k=0; k < NbrImage; k++)
     {     
        Val = 1.- get_random(V1,V2);
        if (Val < V1) Val = V1;
	else if (Val > V2) Val = V2;
	NoiseIma(i,j) +=  sqrt(sigma2*log(1./Val));
     }
     NoiseIma(i,j) /=  (float)NbrImage;
  }
}

/***************************************************************************/

void im_noise_laplace(Ifloat & NoiseIma, int  NbrImage, unsigned int InitRnd)
/* create  noise model
  
      ref: "Digital Image Processing" W.K. PRATT  P.318  
        Laplace Noise  = image_sum ( 1/alpha exp(-X/alpha)  )/NbrImag      
        with alpha = 1 for the model
*/
{
   int Nl = NoiseIma.nl();
   int Nc = NoiseIma.nc();
   int i,j,k;
   float alpha=1.;

   if (InitRnd > 0) init_random (InitRnd);
   for(i=0; i < Nl;i++)
   for(j=0; j < Nc; j++)
   {
     float V1 = FLOAT_EPSILON;
     float V2 = 1. - FLOAT_EPSILON;
     float Val;
     NoiseIma(i,j)=0.;
     for(k=0; k < NbrImage; k++)
     {     
        Val = 1.- get_random(V1,V2);
        if (Val < V1) Val = V1;
	else if (Val > V2) Val = V2;
	NoiseIma(i,j) -=  alpha*log(Val);
        if (NoiseIma(i,j) <= 0)  NoiseIma(i,j) = FLOAT_EPSILON;
      }
      NoiseIma(i,j) /=  (float) NbrImage;
  }
}

/***************************************************************************/

void im_moment4(Ifloat & Ima, double &Mean, double &Sigma, 
                double &Skew, double & Curt, float & Min, 
	        float & Max, int BorderSize)
{
   int i,j,N;
   float *Dat;
   int Nl = Ima.nl();
   int Nc = Ima.nc();
      
   if (BorderSize <= 0)
   {
      Dat = Ima.buffer();
      N = Ima.nl()*Ima.nc();
      moment4(Dat, N, Mean, Sigma,Skew,Curt,Min,Max);
   }
   else
   {
      // cout << "BorderSize = " << BorderSize << endl;
      int Size = BorderSize * (2*Nc+2*Nl - 4*BorderSize);
      int Nx = Ima.nl()*Ima.nc();
      int p = 0;
      fltarray Buff(Nx);
      for (i=BorderSize; i < Ima.nl()-BorderSize; i++)
      for (j=BorderSize; j < Ima.nc()-BorderSize; j++)
           Buff(p++) = Ima(i,j);
      if (p != Nx - Size)
      {
         cout << "Error: SizeBorder = " << BorderSize << endl;
	 cout << "       p = " << p << endl;
	 cout << "       N = " << Nx - Size << endl;
	 cout << "       Size = " << Size << endl;
	 exit(-1);
      }
      Dat = Buff.buffer();
      moment4(Dat, p, Mean, Sigma,  Skew, Curt, Min, Max);
   }
}

/****************************************************************************/

void im_moment4(fltarray & Tab, double &Mean, double &Sigma, 
                double &Skew, double & Curt, float & Min, 
	        float & Max, int BorderSize)
{
   Ifloat Ima;
   Ima.alloc(Tab.buffer(), Tab.ny(), Tab.nx()); 
   im_moment4(Ima, Mean, Sigma, Skew, Curt, Min, Max, BorderSize);
}

/****************************************************************************/
 
void im_hc_test(Ifloat & Ima, float &HC1, float &HC2, int BorderSize, Bool ZeroMean)
{
   int i,j,N;
   float *Dat;
   int Nl = Ima.nl();
   int Nc = Ima.nc();
      
   if (BorderSize <= 0)
   {
      Dat = Ima.buffer();
      N = Ima.nl()*Ima.nc();
      if (ZeroMean == True) hc_test(Dat,N,  HC1,  HC2, 0.);
      else hc_test(Dat,N,  HC1,  HC2);
   }
   else
   {
      // cout << "BorderSize = " << BorderSize << endl;
      int Size = BorderSize * (2*Nc+2*Nl - 4*BorderSize);
      int Nx = Ima.nl()*Ima.nc();
      int p = 0;
      fltarray Buff(Nx);
      for (i=BorderSize; i < Ima.nl()-BorderSize; i++)
      for (j=BorderSize; j < Ima.nc()-BorderSize; j++)
           Buff(p++) = Ima(i,j);
      if (p != Nx - Size)
      {
         cout << "Error: SizeBorder = " << BorderSize << endl;
	 cout << "       p = " << p << endl;
	 cout << "       N = " << Nx - Size << endl;
	 cout << "       Size = " << Size << endl;
	 exit(-1);
      }
      Dat = Buff.buffer();
      if (ZeroMean == True) hc_test(Dat, p,  HC1,  HC2,  0.);
      else hc_test(Dat, p,  HC1,  HC2);
   }
}

/****************************************************************************/

void  im_hc_test(fltarray & Tab,  float &HC1, float &HC2, int BorderSize, Bool ZeroMean)
{
   Ifloat Ima;
   Ima.alloc(Tab.buffer(), Tab.ny(), Tab.nx()); 
   im_hc_test(Ima, HC1, HC2, BorderSize, ZeroMean);
}

/****************************************************************************/

void im_gaussianity_test(Ifloat & Ima, float &T1, float &T2, int BorderSize)
{
   int i,j,N;
   float *Dat;
   int Nl = Ima.nl();
   int Nc = Ima.nc();
   int Size = BorderSize * (2*Nc+2*Nl - 4*BorderSize);
   int Nx = Ima.nl()*Ima.nc();
   int p = 0;
   fltarray Buff(Nx);
   for (i=BorderSize; i < Ima.nl()-BorderSize; i++)
   for (j=BorderSize; j < Ima.nc()-BorderSize; j++) Buff(p++) = Ima(i,j);
   if (p != Nx - Size)
   {
      cout << "Error: SizeBorder = " << BorderSize << endl;
      cout << "       p = " << p << endl;
      cout << "       N = " << Nx - Size << endl;
      cout << "       Size = " << Size << endl;
      exit(-1);
   }
   Dat = Buff.buffer();
   gausstest(Dat, p,  T1,  T2);
}

/****************************************************************************/
 
void im_gaussianity_test(fltarray & Tab,  float &T1, float &T2, int BorderSize)
{
   Ifloat Ima;
   Ima.alloc(Tab.buffer(), Tab.ny(), Tab.nx()); 
   im_gaussianity_test(Ima, T1, T2, BorderSize);
}

/****************************************************************************/
