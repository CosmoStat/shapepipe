/******************************************************************************
**                   Copyright (C) 1995 by CEA
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
**    File:  IM_OperFloat.cc
**
*******************************************************************************
**
**    DESCRIPTION  
**    -----------  
**    
*******************************************************************************
**
** double sigma(Ifloat & Image)
**
** returns the standard deviation
**
*******************************************************************************
** 
** double average(Ifloat & Image)
** 
** returns the average
**
*******************************************************************************
**
** float min(Ifloat & Image)
** 
** returns the min
**
*******************************************************************************
**
** float max(Ifloat & Image)
** 
** returns the max
**
********************************************************************************
**
** float flux(Ifloat & Image, int N)
** 
** returns the flux
**
********************************************************************************
**
** float energy(Ifloat & Image, int N)
** 
** returns the energy
**
*******************************************************************************
**
** const Ifloat  operator +(const Ifloat &Im1, const Ifloat &Im2)
**
** returns the addition of two images
**
*******************************************************************************
**
** const Ifloat  operator -(const Ifloat &Im1, const Ifloat &Im2)
**
** returns the soustraction of two images
**
*******************************************************************************
**
** const Ifloat operator *(const Ifloat &Im1, const Ifloat &Im2)
**
** returns the multiplication of two images
**
*******************************************************************************
**
** const Ifloat  operator /(const Ifloat &Im1, const Ifloat &Im2)
**
** returns the division of two images
** 
*******************************************************************************
**
** const Ifloat conv_direct(const Ifloat &Im1, const Ifloat &Im2, 
**                     type_border Type_Border)
**
** returns the convolution product  Im1 * Im2
**
*******************************************************************************
**
** float lib_mat_correl (const Ifloat & Tab_Im1, const Ifloat & Tab_Im2)
** 
** returns the corelation between two images
**
*******************************************************************************
**
** void threshold (Ifloat &Image, float T)
**
** threshold an image
** all values less than T are set to zero
**
*******************************************************************************
**
** float total_in_block(Ifloat & Image, int i, int j, int Neighbour)
**
** return the total flux in a box of size Neighbour defined around i,j
**
******************************************************************************/


// static char sccsid[] = "@(#)IM_OperFloat.cc 3.1 96/05/02 CEA 1995 @(#)";


#include "IM_Obj.h"

float BadPixalVal = 0.; 
Bool BadPixel=False;

/******************************************************************************/

void ifloat (Ifloat& iflt, Iint& iint) {

   if (    iflt.n_elem() != iint.n_elem()
        || iflt.naxis()   != iint.naxis()) exit (-1); 
   for (int i=0;i<iflt.n_elem();i++)
      iflt(i) = (float)iint(i);
}


/******************************************************************************/

void threshold (Ifloat &Image, float T)
{
    int i,j;
    for (i = 0; i < Image.nl(); i++)
    for (j = 0; j < Image.nc(); j++) if (Image(i,j) < T) Image(i,j) = 0.;
}

/******************************************************************************/

float total_in_block(Ifloat & Image, int i, int j, int Neighbour)
{
    int k,l;
    float Total=0.;

    for (k =i-Neighbour; k <= i+Neighbour; k++)
    for (l =j-Neighbour; l <= j+Neighbour; l++) Total += Image(k, l, I_ZERO);

    return(Total);
}

/************************************************************************/

double sigma(const Ifloat &Image)
{
    int i,j,Np=0;
    double Sigma, Moy;
    int Nl = Image.nl();
    int Nc = Image.nc();
    
    Moy = Sigma = 0.;
    for (i = 0; i < Nl; i++)
    for (j = 0; j < Nc; j++)
    {
        if ((BadPixel == False) || (ABS(Image(i,j)-BadPixalVal)> FLOAT_EPSILON))
        {
            Moy += Image(i,j);
            Sigma += Image(i,j)*Image(i,j);
            Np++;
        }
    }
    Moy /= (float) Np;
    Sigma /= (float) Np;
    Sigma = sqrt(Sigma - Moy * Moy);
    return Sigma;
};

/******************************************************************************/

double average(const Ifloat &Image)
{
    int i,j;
    double Moy;
    int Nl = Image.nl();
    int Nc = Image.nc();

    Moy = 0.;
    for (i = 0; i < Nl; i++)
    for (j = 0; j < Nc; j++)
    {
        Moy += Image(i,j);
    }
    Moy /= (double) (Nl*Nc);
    return Moy;
}

/******************************************************************************/

float min (const Ifloat &Image) 
{
    int i,j;
    float Val;
    int Nl = Image.nl();
    int Nc = Image.nc();

    Val = Image(0,0);
    for (i = 0; i < Nl; i++)
    for (j = 0; j < Nc; j++)
       if (Val > Image(i,j)) Val = Image(i,j);
    return Val;
}

/******************************************************************************/

float max (const Ifloat &Image)
{
    int i,j;
    float Val;
    int Nl = Image.nl();
    int Nc = Image.nc();

    Val = Image(0,0);
    for (i = 0; i < Nl; i++)
    for (j = 0; j < Nc; j++)
       if (Val < Image(i,j)) Val = Image(i,j);
    return Val;
}

/******************************************************************************/

float flux (const Ifloat &Image)
{
    int i,j;
    float Val;
    int Nl = Image.nl();
    int Nc = Image.nc();

    Val = 0.;
    for (i = 0; i < Nl; i++)
    for (j = 0; j < Nc; j++)
              Val += Image(i,j);
    return Val;
}

/******************************************************************************/

double total (const Ifloat &Image)
{
    int i,j;
    double Val;
    int Nl = Image.nl();
    int Nc = Image.nc();

    Val = 0.;
    for (i = 0; i < Nl; i++)
    for (j = 0; j < Nc; j++) Val += Image(i,j);
    return Val;
}

/******************************************************************************/

void norm_flux (Ifloat &Image, float Val_Norm)
{
    int i,j;
    float F = flux(Image) / Val_Norm;
    for (i = 0; i < Image.nl(); i++)
    for (j = 0; j < Image.nc(); j++) Image(i,j) /= F;
}

/******************************************************************************/

double energy (const Ifloat &Image)
{
    int i,j;
    double Val;
    int Nl = Image.nl();
    int Nc = Image.nc();

    Val = 0.;
    for (i = 0; i < Nl; i++)
    for (j = 0; j < Nc; j++)
            Val += Image(i,j)*Image(i,j);
    return Val;
}
             
/****************************************************************************/


const Ifloat  operator +(const Ifloat &Im1, const Ifloat &Im2)
{
     int Nl = Im1.nl();
     int Nc = Im1.nc();
     Ifloat *Result = new Ifloat(Nl,Nc, "op+");
     int i,j;

     for (i=0; i < Nl; i++) {
         for (j=0; j < Nc; j++) {
             (*Result)(i,j) = Im1(i,j) + Im2(i,j);
         }
     }
     return (*Result);
}

/****************************************************************************/

const Ifloat   operator -(const Ifloat &Im1, const Ifloat &Im2)
{
     int Nl = Im1.nl();
     int Nc = Im1.nc();
     Ifloat *Result = new Ifloat(Nl,Nc, "op-");
     int i,j;

     for (i=0; i < Nl; i++) {
         for (j=0; j < Nc; j++) {
             (*Result)(i,j) = Im1(i,j) - Im2(i,j);
         }
     }
     return (*Result);
}

/****************************************************************************/

const Ifloat operator *(const Ifloat &Im1, const Ifloat &Im2)
{
     int Nl = Im1.nl();
     int Nc = Im1.nc();
     Ifloat *Result = new Ifloat(Nl,Nc, "op*");
     int i,j;

     for (i=0; i < Nl; i++) {
         for (j=0; j < Nc; j++) {
             (*Result)(i,j) = Im1(i,j) * Im2(i,j);
         }
     }
     return (*Result);
}

/****************************************************************************/

const Ifloat operator *(const Ifloat &Im1,  float coef)
{
     int Nl = Im1.nl();
     int Nc = Im1.nc();
     Ifloat *Result = new Ifloat(Nl,Nc, "op*");
     int i,j;

     for (i=0; i < Nl; i++) {
         for (j=0; j < Nc; j++) {
             (*Result)(i,j) = Im1(i,j) * coef;
         }
     }
     return (*Result);
}

/****************************************************************************/

const Ifloat operator *(float coef, const Ifloat &Im1)
{
     int Nl = Im1.nl();
     int Nc = Im1.nc();
     Ifloat *Result = new Ifloat(Nl,Nc, "op*");
     int i,j;

     for (i=0; i < Nl; i++) {
         for (j=0; j < Nc; j++) {
             (*Result)(i,j) = Im1(i,j) * coef;
         }
     }
     return (*Result);
}

/****************************************************************************/

const Ifloat  operator /(const Ifloat &Im1, const Ifloat &Im2)
{
     int Nl = Im1.nl();
     int Nc = Im1.nc();
     Ifloat *Result = new Ifloat(Nl,Nc, "op/");
     int i,j;

     for (i=0; i < Nl; i++) {
         for (j=0; j < Nc; j++) {
             if (ABS(Im2(i,j)) > FLOAT_EPSILON) 
                  (*Result)(i,j) = Im1(i,j) / Im2(i,j);
             else (*Result)(i,j) = 0;
         }
     }
     return (*Result);
}

/****************************************************************************/

const Ifloat conv_direct(const Ifloat &Im1, const Ifloat &Im2, type_border Type_Border)
{
    int Nl1 = Im1.nl();
    int Nc1 = Im1.nc();
    int Nl2 = Im2.nl();
    int Nc2 = Im2.nc();
    int Nl_2,Nc_2;
    int i,j,k,l;
    int ii,jj;

    Ifloat *Result = new Ifloat(Nl1,Nc1, "conv_direct");

    Nl_2 = Nl2 / 2;
    Nc_2 = Nc2 / 2;
    for (i = 0; i < Nl1; i++)
    {
        ii = i + Nl_2;
        for (j = 0; j < Nc1; j++)
        {
            (*Result)(i,j) = 0;
            jj = j + Nc_2;
            for (k = 0; k < Nl2; k++)
            for (l = 0; l < Nc2; l++)
                (*Result)(i,j) += Im1(ii - k, jj - l, Type_Border) * Im2(k, l);
       }
    }
    return (*Result);
}

/****************************************************************************/

void convol(const Ifloat &Im1, const Ifloat &Im2, Ifloat &Result,
            type_border Type_Border)
{
    int Nl1 = Im1.nl();
    int Nc1 = Im1.nc();
    int Nl2 = Im2.nl();
    int Nc2 = Im2.nc();
    int Nl_2,Nc_2;
    int i,j,k,l;
    int ii,jj;

    Nl_2 = Nl2 / 2;
    Nc_2 = Nc2 / 2;
    for (i = 0; i < Nl1; i++)
    {
        ii = i - Nl_2;
        for (j = 0; j < Nc1; j++)
        {
            Result(i,j) = 0;
            jj = j - Nc_2;
            for (k = 0; k < Nl2; k++)
            for (l = 0; l < Nc2; l++)
                Result(i,j) += Im1(ii + k, jj + l, Type_Border) * Im2(k, l);
       }
    }
 }
 
/****************************************************************************/


float correlation (const Ifloat & Im1, const Ifloat & Im2)
 {
     double Sum_X2,Sum_Y2,Sum_XY;
     int i,j;
     int Nl = Im1.nl();
     int Nc = Im1.nc();
     double Coef;
     float ValRet;
 
    // test image size
    if ((Im1.nl() != Im2.nl()) || ( Im1.nc() != Im2.nc()))
    {
       cerr << "Error in correlation routine: images have different sizes ..." << endl ; 
       cerr << "   image 1: " <<  Im1.nl() << "X"  <<  Im1.nc() << endl ;
       cerr << "   image 2: " <<  Im2.nl() << "X"  <<  Im2.nc() << endl ;
       exit(-1);
    }
    
     Sum_X2 = Sum_Y2 = Sum_XY = 0.;
     for (i = 0; i < Nl; i++)
     for (j = 0; j < Nc; j++)
     {
         Sum_X2 += Im1 (i,j) * Im1 (i,j);
         Sum_Y2 += Im2 (i,j) * Im2 (i,j);
         Sum_XY += Im1 (i,j) * Im2 (i,j);
     }
     Coef = sqrt (Sum_X2 * Sum_Y2);
     if (Coef > 0.) ValRet = (float)(Sum_XY / Coef);
     else ValRet = 0;
     return (ValRet);
 }
 
/****************************************************************************/

float correlation (const Ifloat & Im1, const Ifloat & Im2, float ValMin)
{
    double Sum_X2,Sum_Y2,Sum_XY;
    int i,j;
    int Nl = Im1.nl();
    int Nc = Im1.nc();
    double Coef;
    float ValRet;

   // test image size
   if ((Im1.nl() != Im2.nl()) || ( Im1.nc() != Im2.nc()))
   {
      cerr << "Error in correllation routine: images have different sizes ..." << endl ; 
      cerr << "   image 1: " <<  Im1.nl() << "X"  <<  Im1.nc() << endl ;
      cerr << "   image 2: " <<  Im2.nl() << "X"  <<  Im2.nc() << endl ;
      exit(-1);
   }
       
    Sum_X2 = Sum_Y2 = Sum_XY = 0.;
    for (i = 0; i < Nl; i++)
    for (j = 0; j < Nc; j++)
    {
        if (ABS(Im1 (i,j)) > ValMin)
        {
           Sum_X2 += Im1 (i,j) * Im1 (i,j);
           Sum_Y2 += Im2 (i,j) * Im2 (i,j);
           Sum_XY += Im1 (i,j) * Im2 (i,j);
        }
    }
    Coef = sqrt (Sum_X2 * Sum_Y2);
    if (Coef > 0.) ValRet = (float)(Sum_XY / Coef);
    else ValRet = 0;
    return (ValRet);
}

 /****************************************************************************/

float regression (const Ifloat &Im1, const Ifloat & Im2, float ValMin)
// Search the regression coefficient such that:
//   Im1 = Im2 * CoefRegression
{
   int i,j ;
   double Sum_XY=0., Sum_X2=0.;
   float ValRet;
 
   // test image size
   if ((Im1.nl() != Im2.nl()) || ( Im1.nc() != Im2.nc()))
   {
      cerr << "Error in regression routine: images have different sizes ..." << endl ; 
      cerr << "   image 1: " <<  Im1.nl() << "X"  <<  Im1.nc() << endl ;
      cerr << "   image 2: " <<  Im2.nl() << "X"  <<  Im2.nc() << endl ;
      exit(-1);
   }
 	 
   for (i = 0; i <  Im1.nl(); i++)
   for (j = 0; j <  Im1.nc(); j++)
   {
      if ((ABS(Im1(i,j)) > ValMin) && (ABS(Im2(i,j))>  ValMin))
      {
	  Sum_XY+=Im1(i,j)*Im2(i,j);
	  Sum_X2+=Im2(i,j)*Im2(i,j);
      }
   }

   if (Sum_X2 != 0.)  ValRet = (float) (Sum_XY / Sum_X2);
   else  ValRet = 0.;
   
   return ValRet;
}


/***************************************************************************/
 
float regression (const Ifloat &Im1, const Ifloat & Im2)
// Search the regression coefficient such that:
//   Im1 = Im2 * CoefRegression
{
   int i,j ;
   double Sum_XY=0., Sum_X2=0.;
   float ValRet;
 
   // test image size
   if ((Im1.nl() != Im2.nl()) || ( Im1.nc() != Im2.nc()))
   {
      cerr << "Error in regression routine: images have different sizes ..." << endl ; 
      cerr << "   image 1: " <<  Im1.nl() << "X"  <<  Im1.nc() << endl ;
      cerr << "   image 2: " <<  Im2.nl() << "X"  <<  Im2.nc() << endl ;
      exit(-1);
   }
 	 
   for (i = 0; i <  Im1.nl(); i++)
   for (j = 0; j <  Im1.nc(); j++)
   {
      Sum_XY+=Im1(i,j)*Im2(i,j);
      Sum_X2+=Im2(i,j)*Im2(i,j);
   }
   if (Sum_X2 != 0.)  ValRet = (float) (Sum_XY / Sum_X2);
   else  ValRet = 0.;
   
   return ValRet;
}

/***************************************************************************/
