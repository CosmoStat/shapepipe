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
** float detect_noise_from_med (const Ifloat &Image)
**
** returns an estimation of the standard deviation of the noise in Image.
** This estimation is computed from the difference between Image and 
** Image smoothed by the mediane filter
**
*******************************************************************************
**
** float detect_noise_sigma (Ifloat &Image, Bool Average_Non_Null,int Nit)
**
** returns an estimation of the standard deviation of the noise in an Image
** by doing a 3sigma clipping
**
*******************************************************************************
**
** void sigma_clip(const Ifloat &Image, float &Mean, float &Sigma, int Nit)
**
** Sigma clip of an image: return the mean value of the background and the
** standard deviation
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
** void block_transform(Ifloat &Ima,float * & PtrMean, 
**                      float * & PtrSigma, int &Nb, int SizeBlock)
**
** Compute the sigma and the mean in each block of an image
**
*******************************************************************************
**
** void sort_bulle (float *Data, int Np)
**
**  sort the data by the bulle method
**
*******************************************************************************
**
** float detect_noise_from_block (Ifloat &Ima, int SizeBlock)
**
** detect noise from the block transform
**
*******************************************************************************
*
** float sigma_clip_int (int *Image, int Nl, int Nc, int Option, int Nit)
** 
** sigma clipping on integer 2d array.
** if opt=0, assume that the data have a zero meam.
**
******************************************************************************/

// static char sccsid[] = "@(#)MR_Sigma.cc 3.3 96/07/08 CEA 1996 @(#)";


#include "IM_Obj.h"
#include "IM_Sigma.h"


extern float BadPixalVal; 
extern Bool BadPixel;

/******************************************************************************/

void sigma_clip(const Ifloat &Image, float &Moy, float &Sig, int Nit)
{
    int It, i, j;
    double S0,S1,S2,Sm=0,x,Mean,Sigma=0.;
    int Nl = Image.nl();
    int Nc = Image.nc();

    Mean = 0.;
    for (It = 0; It < Nit; It++)
    {
       S0 = S1 = S2 = 0.;
       for (i = 0; i < Nl; i++)
       for (j = 0; j < Nc; j++)
       {
           x = Image(i,j);
           if ((BadPixel == False) || (ABS(x-BadPixalVal) != 0))
           {
	   if ((It == 0) || (ABS(x - Mean) < Sm))
	   { 
	       S0 ++;
	       S1 += x;
	       S2 += x*x;
	   }}
       }
       if (S0 == 0) S0=1;
       Mean = S1 / S0;
       Sigma = S2/S0-  Mean * Mean;
       if ( Sigma >= 0.) Sigma = sqrt(Sigma);
       else Sigma = 0.;
       Sm = 3. * Sigma;
    }
    Moy= (float) Mean;
    Sig =(float) Sigma;
}
 
/***************************************************************************/
 
void sigma_clip_int (int *Image, float & Moy, float & Sig, 
                     int Nl, int Nc, int Option, int Nit)
{
    int It,Il;
    double S0,S1,S2,Sm=0,x;
    double Val, Moyenne, Sigma;
 
    Sigma=0; 
    Moyenne=0;
    for (It = 0; It < Nit; It++)
    {
       S0 = S1 = S2 = 0.;
       for (Il = 0; Il < Nl*Nc; Il++)
       {
           x = Image[Il];
	   if ((It == 0) || (fabs(x - Moyenne) < Sm))
	   { 
	       S0 ++;
	       S1 += x;
	       S2 += x*x;
	   }
       }       
       if (S0 == 0) S0=1;
       if (Option == 1)
       {
       	   Moyenne = S1 / S0;
       	   Val = S2/S0- Moyenne * Moyenne;
	   if (Val >= 0.) Sigma = sqrt(Val);
	   else Sigma = 0.;
       }
       else
       {
       	    Moyenne = 0.;
	    Sigma = sqrt(S2/S0);
       }
       Sm = 3. * Sigma;
    }
    
    Moy= (float) Moyenne;
    Sig =(float) Sigma;
}

/****************************************************************************/
  
float sigma_clip_int (int *Image, int Nl, int Nc, int Option, int Nit)
{
    float Sigma, Moyenne;
    sigma_clip_int (Image, Moyenne, Sigma, Nl, Nc, Option, Nit);
    return (Sigma);
}

/***************************************************************************/

float detect_noise_sigma (Ifloat &Image, Bool Average_Non_Null, int Nit)
{
    int Nl = Image.nl();
    int Nc = Image.nc();
    double Sigma;
    int N=Nl*Nc;
    Sigma=get_sigma_clip (Image.buffer(), N, Nit, 
                          Average_Non_Null, BadPixel, BadPixalVal);
    return ((float) Sigma);
}
float detect_noise_sigma (const Ifloat &Image, Bool Average_Non_Null, int Nit)
{
    int Nl = Image.nl();
    int Nc = Image.nc();
    double Sigma;
    int N=Nl*Nc;
    Sigma=get_sigma_clip (Image.buffer(), N, Nit, 
                          Average_Non_Null, BadPixel, BadPixalVal);
    return ((float) Sigma);
}
/***************************************************************************/

float detect_noise_from_med (const Ifloat &Image)
{
    double Noise;
    int i,j;
    int Nl = Image.nl();
    int Nc = Image.nc();
    int N=Nl*Nc;
    Bool Average_Non_Null = True;
    int Nit = 3;
    int NpixOk = 0;
    
    Ifloat Buff(Image.nl(), Image.nc(), "Buff noise estimation");
    fltarray Buf2(Image.nl()*Image.nc());
    smooth_mediane (Image, Buff, DEFAULT_BORDER, 0, 3);
    for (i=0; i < Image.nl(); i++)
    for (j=0; j < Image.nc(); j++) 
    {
        Buff(i,j) = Image(i,j) - Buff(i,j);
        if ((BadPixel == True) && ((Image(i,j) == BadPixalVal) || (Image(i-1,j,I_CONT) == BadPixalVal) || (Image(i+1,j,I_CONT) == BadPixalVal)
            || (Image(i,j-1,I_CONT) == BadPixalVal) || (Image(i,j+1,I_CONT)  == BadPixalVal)))
        {
           Buff(i,j) = BadPixalVal;
        }
    }
//    INFO(Buff, "diff");
    
    Noise=get_sigma_clip (Buff.buffer(), N, Nit, 
                          Average_Non_Null, BadPixel, BadPixalVal);
    Noise /= 0.972463; 
//     printf("Calc sigma(%d) = %f\n", N, (float) Noise);
//     if (Noise < FLOAT_EPSILON) Noise = detect_noise_from_mad(Buff);
//     INFO(Buff, "buff");
    return ((float) Noise);
}

/***************************************************************************/

float detect_noise_from_bspline (const Ifloat &Image)
{
    int i,j;
    float Noise;

    Ifloat Buff(Image.nl(), Image.nc(), "Buff noise estimation");
    smooth_bspline (Image, Buff);
    for (i=0; i < Image.nl(); i++)
    for (j=0; j < Image.nc(); j++) Buff(i,j) = Image(i,j) - Buff(i,j);
    Noise = detect_noise_sigma(Buff, False, 3);
    return (Noise/0.889434);
}

/***************************************************************************/

float detect_noise_from_mad (Ifloat &Image, Bool CenterOnly)
{
    float Noise;
    
    if (CenterOnly == False) Noise  = get_sigma_mad(Image.buffer(), Image.n_elem());
    else
    {   
       int Nl = Image.nx();
       int Nc = Image.nc();
       int Nl1 = Nl / 4;
       int Nc1 = Nc / 4;
       int Depi = Nl/2-Nl1/2 ;
       int Depj = Nc/2-Nc1/2 ;
       Ifloat I(Nl1,Nc1);
       for (int i=0; i < Nl1; i++)
       for (int j=0; j < Nc1; j++) I(i,j) = Image(Depi+i,Depj+j);  
       Noise  = get_sigma_mad(I.buffer(), I.n_elem());
    }
    return (Noise);
}

/***************************************************************************/

void block_transform(Ifloat &Ima,float * & PtrMean, float * & PtrSigma, int &Nb, int SizeBlock)
{
    int i,j,Nc,Nl, Nlb,Ncb;
    int ind=0;
    float Moy = 0., Ecart = 0.;
    int Size = SizeBlock*SizeBlock;

    Nc = Ima.nc();
    Nl = Ima.nl();
    Ncb = Nc / SizeBlock;
    Nlb = Nl / SizeBlock;
    Nb = Ncb*Nlb;
    PtrMean = new float [Nb];
    PtrSigma = new float [Nb];

    for (i = 0; i < Nlb; i++)
    for (j = 0; j < Ncb; j++)
    {
        Moy = 0.;
        Ecart = 0.;
        for (int k = 0; k < SizeBlock; k++)
        for (int l = 0; l < SizeBlock; l++)
        {
           float Val;
           Val = Ima(i*SizeBlock+k, j*SizeBlock+l);
           Moy += Val;
           Ecart += Val*Val;
        }
        Moy /= (float) Size;
        Ecart = Ecart / (float) Size - Moy*Moy;
        if (Ecart < 0.) Ecart = 0.;
        Ecart = sqrt(Ecart);

        PtrMean[ind] = Moy;
        PtrSigma[ind] = Ecart;
        ind ++;
    }
}

/***************************************************************************/

void sort_bulle (float *Data, int Np)
{
    int i,j;
    float Aux;

    for (i = 1; i <= Np/2+1; i++)
    for (j = 0; j < Np - i; j++)
    {
       if (Data[j] > Data[j+1])
       {
          Aux = Data[j+1];
          Data[j+1] = Data[j];
          Data[j] = Aux;
        }
    }
}

/***************************************************************************/

float detect_noise_from_block (Ifloat &Ima, int SizeBlock)
{
    int i,Nb;
    float ValReturn=0.;
    float *PtrMean, *PtrSigma;

    block_transform(Ima, PtrMean, PtrSigma, Nb, SizeBlock);
    sort_bulle(PtrMean, Nb);
    sort_bulle(PtrSigma, Nb);
    Nb = (int) (Nb * 0.1);
    
    for (i=0; i < Nb; i++) ValReturn += PtrSigma[i];
    ValReturn /= (float) Nb;

   return(ValReturn);
}

/***************************************************************************/
