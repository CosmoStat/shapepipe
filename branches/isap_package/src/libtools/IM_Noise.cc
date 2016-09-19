/*******************************************************************************
**
**    UNIT
**
**    Version: 3.2
**
**    Author: Jean-Luc Starck
**
**    Date:  96/06/13 
**    
**    File:  IM_Noise.cc
**
*******************************************************************************
**
**    DESCRIPTION routines for noise treatment 
**    -----------  
**                 
*******************************************************************************
**
** void noise_poisson_transform (Ifloat &Imag_in, Ifloat &Imag_out)
**
** compute the Anscombe's transform. The poisson noise contained in Resi
** is transform in a gaussian noise (with a standard deviation equal
** to 1) in image Sn.
**
*******************************************************************************
**
** void noise_inverse_poisson_transform (Ifloat &Imag_in, Ifloat &Imag_out)
**
** compute the inverse transform of Anscombe's transform
**
*******************************************************************************
**
** void noise_residu_estimation (Ifloat &In, Ifloat &Sn, Ifloat &Rn)
**
** compute the signicant residual Rn from the filtered Sn and the the image In
** In: in Ifloat = approximation of an image I
** Sn: in Ifloat = Filtered of T(I) - T(In)
** Rn: out Ifloat = Filtered of I - In 
** This routine extract the significant values in the residual, and is used
** when the image contained a poisson noise, and in the iterative processes
**
**  I = In + Rn
**  I is the image, In is an approximation of I without noise, and Rn 
**  is the error associated to In: Rn = I - In
**   
**  Sn corresponds the image T(I) - T(In) were T is the Anscombe transform
**  Rn filtered is given by:
**  Rn = Sn * ( Sn/4 + T(In) )
**
*******************************************************************************
**
** void noise_threshold_imag (Ifloat &Imag, Ifloat Level)
**
** Threshold all values which have an absolute value less than Level
**
*******************************************************************************
**
** void noise_support_threshold_imag (Ifloat &Imag, float Level, 
**                                    Ifloat &Support, Bool WriteSupport)
**
** Threshold all values which have an absolute value less than Level and 
** which is not in the support
** if WriteSupport is True, then if Imag(i,j) > Level Support(i,j) = 1
**
*******************************************************************************
**
** void noise_log_transform (Ifloat &Imag, Ifloat &Result)
** 
** apply the following transform to the input image
** if (Imag > 0) Result = log(Imag+1) else Result=0
**
*******************************************************************************
**
** void im_sigma(Ifloat &Data, Ifloat &Ima_Sigma, int SizeBlock, int Nit)
**
** computes the standard deviation in a block around each pixel.
** if Nit > 1, then compute the standard deviation by a k-sigma clipping,
** and give an estimation of the local gaussian noise.
**
******************************************************************************/

// static char sccsid[] = "@(#)IM_Noise.cc 3.2 96/06/13 CEA 1995 @(#)";

#include "IM_Noise.h"


/* definition for the poissonian noise, when a gaussian component is added */
float PasCodeur=1.;
float SigmaGauss=0.;
float MeanGauss=0.;
 
/***************************************************************************/

const char * StringNoise (type_noise type)
{
    switch (type)
    {
    case NOISE_GAUSSIAN: 
      return ("Gaussian noise");break;
    case NOISE_POISSON: 
      return ("Poisson noise");break;
    case NOISE_GAUSS_POISSON: 
      return ("Poisson noise + Gaussian noise");break;
    case NOISE_EVENT_POISSON: 
      return ("Poisson noise with few events");break;
    case NOISE_MULTI:
      return ("Multiplicative noise");break;
    case NOISE_NON_UNI_ADD: 
      return ("Non-stationary additive noise");break;
    case NOISE_NON_UNI_MULT:
      return ("Non-stationary multiplicative noise");break;
    case NOISE_UNI_UNDEFINED:
      return ("Undefined stationary noise");break;
    case NOISE_UNDEFINED:
      return ("Undefined noise");break;
    case NOISE_SPECKLE:
      return ("Speckle noise");break;
    case NOISE_CORREL:
      return ("Stationary correlated noise");break;
    }
    return ("Error: bad type of noise");
}

/***************************************************************************/

void noise_poisson_transform (int *Resi, int *Sn, int Nl, int Nc)
{
    int N = Nl*Nc;
    int i,Cpt=0;
    float C1, C2;
    float Val;

    C1 = 2. / PasCodeur;
    C2 = 3. / 8. * PasCodeur*PasCodeur + 
                  SigmaGauss*SigmaGauss - PasCodeur*MeanGauss;

    for (i = 0; i < N; i++)
    {
       Val = Resi[i]*PasCodeur + C2;
       if (Val < 0.) {Val=0;Sn[i]=0;Cpt++;}
       else  Sn[i] = (int) (C1*sqrt(Val) + 0.5);
    }
    if (Cpt > 0) cout << "WARNING: Nbr < 0 = " << Cpt << endl;
}

/***************************************************************************/

void noise_poisson_transform (fltarray &Resi, fltarray &Sn)
{
    int N = Resi.n_elem();
    int i,Cpt=0;
    float C1, C2;
    float Val;
    static int Pas = 0;

    C1 = 2. / PasCodeur;
    C2 = 3. / 8. * PasCodeur*PasCodeur + 
                  SigmaGauss*SigmaGauss - PasCodeur*MeanGauss;
    if (Pas == 0)
    {
        // printf ("\n T(x) = %5.2f sqrt(%5.2f x + %5.2f)\n", C1, PasCodeur, C2);
        Pas = 1;
    }

    for (i = 0; i < N; i++)
    {
       Val = Resi(i)*PasCodeur + C2;
       if (Val < 0.) {Val=0;Sn(i)=0.;Cpt++;}
       else  Sn(i) = C1*sqrt(Val);
    }
    if (Cpt > 0) cout << "WARNING: Nbr < 0 = " << Cpt << endl;
}


/****************************************************************************/

void noise_inverse_poisson_transform (int *Sn, int *Resi, int Nl, int Nc)
{
    int N = Nl*Nc;
    int i;
    float C1, C2;
 
    C1 = 2. / PasCodeur;
    C1 *= C1;
    C2 = 3. / 8. * PasCodeur*PasCodeur + 
                  SigmaGauss*SigmaGauss - PasCodeur*MeanGauss;

    for (i = 0; i < N; i++) Resi[i] = (int) (Sn[i]*Sn[i]/(C1*PasCodeur) - C2/PasCodeur+0.5);
}


/****************************************************************************/

void noise_inverse_poisson_transform (fltarray &Sn, fltarray &Resi)
{
    int N = Sn.n_elem();
    int i;
    float C1, C2;
 
    C1 = 2. / PasCodeur;
    C1 *= C1;
    C2 = 3. / 8. * PasCodeur*PasCodeur + 
                  SigmaGauss*SigmaGauss - PasCodeur*MeanGauss;

    for (i = 0; i < N; i++) Resi(i) = Sn(i)*Sn(i)/(C1*PasCodeur) - C2/PasCodeur;
}

/***************************************************************************/

void noise_poisson_transform (Ifloat &Resi, Ifloat &Sn)
{
    int Nl = Resi.nl();
    int Nc = Resi.nc();
    int i,j,Cpt=0;
    float C1, C2;
    float Val;
    static int Pas = 0;

    C1 = 2. / PasCodeur;
    C2 = 3. / 8. * PasCodeur*PasCodeur + 
                  SigmaGauss*SigmaGauss - PasCodeur*MeanGauss;
    if (Pas == 0)
    {
        // printf ("\n T(x) = %5.2f sqrt(%5.2f x + %5.2f)\n", C1, PasCodeur, C2);
        Pas = 1;
    }

    for (i = 0; i < Nl; i++)
    for (j = 0; j < Nc; j++)
    {
       Val = Resi(i,j)*PasCodeur + C2;
       if (Val < 0.) {Val=0;Sn(i,j)=0.;Cpt++;}
       else  Sn(i,j) = C1*sqrt(Val);
    }
    if (Cpt > 0) cout << "WARNING: Nbr < 0 = " << Cpt << endl;
}

/***************************************************************************/

void noise_inverse_poisson_transform (Ifloat &Sn, Ifloat &Resi)
{
    int Nl = Resi.nl();
    int Nc = Resi.nc();
    int i,j;
    float C1, C2;
 
    C1 = 2. / PasCodeur;
    C1 *= C1;
    C2 = 3. / 8. * PasCodeur*PasCodeur + 
                  SigmaGauss*SigmaGauss - PasCodeur*MeanGauss;

    for (i = 0; i < Nl; i++)
    for (j = 0; j < Nc; j++) 
        Resi(i,j) = Sn(i,j)*Sn(i,j)/(C1*PasCodeur) - C2/PasCodeur;
}

/***************************************************************************/

void noise_log_transform (Ifloat &Imag, Ifloat &Result)
{
    int Nl = Imag.nl();
    int Nc = Imag.nc();
    int i,j;
    for (i = 0; i < Nl; i++)
    for (j = 0; j < Nc; j++) 
    {
       if (Imag(i,j) > 0) Result(i,j) = log(Imag(i,j));
       else Result(i,j) = 0.;
    }
}

/***************************************************************************/

void noise_inv_log_transform (Ifloat &Imag, Ifloat &Result)
{
    int Nl = Imag.nl();
    int Nc = Imag.nc();
    int i,j;

    for (i = 0; i < Nl; i++)
    for (j = 0; j < Nc; j++) 
       Result(i,j) = exp(Imag(i,j));
}

/***************************************************************************/

void noise_residu_estimation (Ifloat &In, Ifloat &Sn, Ifloat &Rn)
{
    int Nl = In.nl();
    int Nc = In.nc();
    int i,j;
    float C1, C2;
 
    C1 = 2. / PasCodeur;
    C2 = 3. / 8. * PasCodeur*PasCodeur + 
                  SigmaGauss*SigmaGauss - PasCodeur*MeanGauss;

    for (i = 0; i < Nl; i++)
    for (j = 0; j < Nc; j++)
        Rn(i,j) = Sn(i,j)*(Sn(i,j)*PasCodeur/4. + sqrt(PasCodeur*In(i,j) + C2));
}

/****************************************************************************/

void noise_threshold_imag (Ifloat &Imag, float Level)
{
    int i, j;

    for (i = 0; i < Imag.nl(); i++)
    for (j = 0; j < Imag.nc(); j++)
       if (fabs (Imag(i,j)) < Level) Imag(i,j) = 0.;
}

/****************************************************************************/

void noise_support_threshold_imag (Ifloat &Imag, float Level, 
                                   Ifloat &Support, Bool WriteSupport)
{
    int i, j;
//  int  cpt=0;

    for (i = 0; i < Imag.nl(); i++)
    for (j = 0; j < Imag.nc(); j++)
    {
       if ((WriteSupport) && (fabs (Imag(i,j)) >= Level))  Support(i,j) = 1.;
       Imag(i,j) *= Support(i,j);
// if (fabs (Imag(i,j)) < FLOAT_EPSILON) cpt++;
    }

// cout << "Thresholding: Level = " << Level << ", Pc = " <<
//            ((float)cpt/(float)(Imag.nl()*Imag.nc())*100.) << "%" << endl;
}

/****************************************************************************/

void im_sigma(Ifloat &Data, Ifloat &Ima_Sigma, int SizeBlock, int Nit)
{
  int i,j,k,l,Nc,Nl,It;
  Nl = Data.nl();
  Nc = Data.nc();
  int Dsize = SizeBlock / 2;
  double x,S0,S1,S2,Mean=0.,Sigma=0.,Sm=0.,Val;
  int Kb,Ke,Lb,Le;

  for (i=0; i< Nl-1; i++)
  for (j=0; j< Nc-1; j++)
  {
      for (It = 0; It < Nit; It++)
      {
         S0 = 0.;
         S1 = 0.;
         S2 = 0.;
         Kb = (i-Dsize >= 0) ? i-Dsize : 0;
         Ke = (i+Dsize <= Nl-1) ? i+Dsize : Nl-1;
         Lb = (j-Dsize >= 0) ? j-Dsize : 0;
         Le = (j+Dsize <= Nc-1) ? j+Dsize : Nc-1;
         for (k=Kb; k<= Ke; k++)
         for (l=Lb; l<= Le; l++)
         {
            x = Data(k,l);
	    if ((It == 0) || (fabs(x - Mean) < Sm))
	    {
               S0 ++;
               S1 += x;
               S2 += x*x;
	    }
         }
         Mean = S1 / S0;
         Val = S2/S0 - Mean * Mean;
         if (Val > 1e-7) Sigma = sqrt(Val);
         Sm = 3. * Sigma;
      }
      Ima_Sigma(i,j) = (float) Sigma;
  }
  // io_write_ima_float("xx_sigma",  Ima_Sigma);

}

/********************************************************************/

void im_sigma_mad(Ifloat &Data, Ifloat &Ima_Sigma, int SizeBlock)
{
  int i,j,k,l,Nc,Nl;
  Nl = Data.nl();
  Nc = Data.nc();
  int Dsize = SizeBlock / 2;
  fltarray Buf(SizeBlock*SizeBlock);
  int Kb,Ke,Lb,Le, Ind;

  for (i=0; i< Nl; i++)
  for (j=0; j< Nc; j++)
  {
     Kb = (i-Dsize >= 0) ? i-Dsize : 0;
     Ke = (i+Dsize <= Nl-1) ? i+Dsize : Nl-1;
     Lb = (j-Dsize >= 0) ? j-Dsize : 0;
     Le = (j+Dsize <= Nc-1) ? j+Dsize : Nc-1;
     Ind=0;
     for (k=Kb; k<= Ke; k++)
     for (l=Lb; l<= Le; l++) Buf(Ind++) = Data(k,l);
     Ima_Sigma(i,j) = get_sigma_mad(Buf.buffer(),Ind); 
  }
}

/********************************************************************/

void im_sigma_block_mad(Ifloat &Data, Ifloat &Ima_Sigma, int SizeBlock, int ResolMad)
{
  int i,j,k,l,Nc,Nl;
  Nl = Data.nl();
  Nc = Data.nc();
  int Dsize = SizeBlock / 2;
  fltarray Buf(SizeBlock*SizeBlock);
  int Kb,Ke,Lb,Le, Ind;
  int Nlb = Nl / ResolMad;
  int Ncb = Nc / ResolMad;
  Ifloat Sig(Nlb,Ncb);
  
  for (i=0; i< Nlb; i++)
  for (j=0; j< Ncb; j++)
  {
     Kb = (i*ResolMad-Dsize >= 0) ? i*ResolMad-Dsize : 0;
     Ke = (i*ResolMad+Dsize <= Nl-1) ? i*ResolMad+Dsize : Nl-1;
     Lb = (j*ResolMad-Dsize >= 0) ? j*ResolMad-Dsize : 0;
     Le = (j*ResolMad+Dsize <= Nc-1) ? j*ResolMad+Dsize : Nc-1;
     Ind=0;
     for (k=Kb; k<= Ke; k++)
     for (l=Lb; l<= Le; l++) Buf(Ind++) = Data(k,l);
     Sig(i,j) = get_sigma_mad(Buf.buffer(),Ind); 
  }
  
  for (i=0; i< Nl; i++)
  for (j=0; j< Nc; j++) Ima_Sigma(i,j) = Sig(i/ResolMad,j/ResolMad);
  
}

/********************************************************************/

void im_sigma_block(Ifloat &Data, Ifloat &Ima_Sigma, int SizeBlock, int Nit, int ResolMad)
{
  int i,j,k,l,Nc,Nl,It;
  Nl = Data.nl();
  Nc = Data.nc();
  int Dsize = SizeBlock / 2;
  double x,S0,S1,S2,Mean=0.,Sigma=0.,Sm=0.,Val;
  int Kb,Ke,Lb,Le;
  int Nlb = Nl / ResolMad;
  int Ncb = Nc / ResolMad;
  Ifloat Sig(Nlb,Ncb);
  
  for (i=0; i< Nlb; i++)
  for (j=0; j< Ncb; j++)
  {
      for (It = 0; It < Nit; It++)
      {
         S0 = 0.;
         S1 = 0.;
         S2 = 0.;
         Kb = (i*ResolMad-Dsize >= 0) ? i*ResolMad-Dsize : 0;
         Ke = (i*ResolMad+Dsize <= Nl-1) ? i*ResolMad+Dsize : Nl-1;
         Lb = (j*ResolMad-Dsize >= 0) ? j*ResolMad-Dsize : 0;
         Le = (j*ResolMad+Dsize <= Nc-1) ? j*ResolMad+Dsize : Nc-1;
         for (k=Kb; k<= Ke; k++)
         for (l=Lb; l<= Le; l++)
         {
            x = Data(k,l);
	    if ((It == 0) || (fabs(x - Mean) < Sm))
	    {
               S0 ++;
               S1 += x;
               S2 += x*x;
	    }
         }
         Mean = S1 / S0;
         Val = S2/S0 - Mean * Mean;
         if (Val > 1e-7) Sigma = sqrt(Val);
         Sm = 3. * Sigma;
      }
      Sig(i,j) = (float) Sigma;
  }
  
  for (i=0; i< Nl; i++)
  for (j=0; j< Nc; j++) Ima_Sigma(i,j) = Sig(i/ResolMad,j/ResolMad);
  
  // io_write_ima_float("xx_sigma",  Ima_Sigma);

}
