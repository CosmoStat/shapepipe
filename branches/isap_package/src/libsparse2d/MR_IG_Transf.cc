/*******************************************************************************
**
**    UNIT
**
**    Version: 1.1
**
**    Author: Jean-Luc Starck
**
**    Date:  97/09/02 
**    
**    File:  MR_IG_Transf.cc
**
*******************************************************************************
**    DESCRIPTION  routines used for the G-transform algorithm 
**    -----------  
**
*******************************************************************************
**
** void transform_g (Iint &Imag, Iint &Imag_Out, int Nbr_Iter)
** 
** computes the G tranform of an image with Nbr_Iter oterations
**
*******************************************************************************
**
** void inverse_transform_g (Iint &Imag, Iint &Result, int Nbr_Iter)
**
** reconstruct an image from its G transform
**
***************************************************************************/ 


#include <stdio.h>
#include <math.h>
#include <string.h>

#include "IM_Obj.h"

 /***************************************************************************/

void transfo_g_one_scale (int *Imag, int *Imag_Out, int Nl, int Nc)
{
   int i,j;
   int *Ima_Aux, *Ima_Aux1;
   int Nl1,Nc1;
   Ima_Aux = new int [Nl*Nc];
   Ima_Aux1 = new int [Nl*Nc];
  
   Nl1 = (Nl % 2 == 0) ? Nl : Nl - 1;  
   Nc1 = (Nc % 2 == 0) ? Nc : Nc - 1;
   for (i = 0; i < Nl; i++)
   {
       for (j = 0; j < Nc1-1; j+=2)
       {
           if (Imag[i*Nc+j] <= Imag[i*Nc+j+1]) Ima_Aux[i*Nc+j] = Imag[i*Nc+j];
           else Ima_Aux[i*Nc+j] = Imag[i*Nc+j+1];
           Ima_Aux[i*Nc+j+1] = Imag[i*Nc+j] - Imag[i*Nc+j+1];
       }
   }
   // test non square image
   if (Nc1 != Nc)
   {
      j = Nc-1;
      for (i = 0; i < Nl; i++) Ima_Aux[i*Nc+j] = Imag[i*Nc+j];
   }


   for (i = 0; i < Nl1; i+=2)
   {
       for (j = 0; j < Nc; j++)
       {
           if (Ima_Aux[(i+1)*Nc+j] > Ima_Aux[i*Nc+j]) 
                Ima_Aux1[i*Nc+j] = Ima_Aux[i*Nc+j];
           else Ima_Aux1[i*Nc+j] = Ima_Aux[(i+1)*Nc+j];
           Ima_Aux1[(i+1)*Nc+j] = Ima_Aux[i*Nc+j] - Ima_Aux[(i+1)*Nc+j];
       }
   }
   // test non square image
   if (Nl1 != Nl)
   {
      i = Nl-1;
      for (j = 0; j < Nc; j++) Ima_Aux1[i*Nc+j] =  Ima_Aux[i*Nc+j];
   }
   
   
   int Nls_L = (Nl+1)/2;
   int Ncs_L = (Nc+1)/2;
   int Nls_H =  Nls_L;
   int Ncs_H =  Nc /2;
   int Nls_V =  Nl /2;
   int Ncs_V =  Ncs_L;
   int Nls_D =  Nl /2;
   int Ncs_D =  Nc /2;

   // low resolution band
   for (i = 0; i < Nls_L; i++)
   for (j = 0; j < Ncs_L; j++) Imag_Out[i*Nc+j] = Ima_Aux1[2*i*Nc+2*j];
   
   // horizontal band
   for (i = 0; i < Nls_H; i++)
   for (j = 0; j < Ncs_H; j++) Imag_Out[i*Nc+j+Ncs_L] = Ima_Aux1[2*i*Nc+2*j+1];
   
   // vertical band
   for (i = 0; i < Nls_V; i++)
   for (j = 0; j < Ncs_V; j++) Imag_Out[(i+Nls_L)*Nc+j] = Ima_Aux1[(2*i+1)*Nc+2*j];
   
   // diagonal band
   for (i = 0; i < Nls_D; i++)
   for (j = 0; j < Ncs_D; j++) 
              Imag_Out[(i+Nls_L)*Nc+j+Ncs_L] = Ima_Aux1[(2*i+1)*Nc+2*j+1];
         
   delete [] Ima_Aux;
   delete [] Ima_Aux1;
}

/***************************************************************************/

void transform_g (Iint &Imag, Iint &Imag_Out, int Nbr_Iter)
{
   int i,j,k, Nl1, Nc1, Nl2, Nc2;
   int *Ima_Aux, *Ima_Aux1;
   int Nl = Imag.nl();
   int Nc = Imag.nc();

   Ima_Aux = new int [Nl*Nc];
   Ima_Aux1 = new int [Nl*Nc];

   Nl1 = Nl;
   Nc1 = Nc;
   // copy the image in Ima_Aux
   for (i = 0; i < Nl1; i++)
   for (j = 0; j < Nc1; j++)
                      Ima_Aux[i*Nc1+j] = Imag(i,j);

   for (k = 0; k < Nbr_Iter; k++)
   {
       // transform Ima_Aux
       transfo_g_one_scale (Ima_Aux, Ima_Aux1, Nl1, Nc1);

       // copy the transform part in Imag_Out
       for (i = 0; i < Nl1; i++)
       for (j = 0; j < Nc1; j++)
                Imag_Out(i,j) = Ima_Aux1[i*Nc1+j];

       // next step on a reduced image
       Nl2 = (Nl1+1)/2;
       Nc2 = (Nc1+1)/2;

       // copy the smoothed array to Ima_Aux
       for (i = 0; i < Nl2; i++)
           for (j = 0; j < Nc2; j++)
                Ima_Aux[i*Nc2+j] = Ima_Aux1[i*Nc1+j];
      Nl1 = Nl2;
      Nc1 = Nc2;
   }
   delete [] Ima_Aux;
   delete [] Ima_Aux1;
}

/***************************************************************************/

void transfo_g_inv_one_scale (int *Imag, int *Imag_Out, int Nl, int Nc)
{
   int i,j;
   int *Ima_Aux;
   int *Ima_Aux1;

   Ima_Aux = new int [Nl*Nc];
   Ima_Aux1 = new int [Nl*Nc];
   int Nls_L = (Nl+1)/2;
   int Ncs_L = (Nc+1)/2;
   int Nls_H =  Nls_L;
   int Ncs_H =  Nc /2;
   int Nls_V =  Nl /2;
   int Ncs_V =  Ncs_L;
   int Nls_D =  Nl /2;
   int Ncs_D =  Nc /2;
   int Nl1 = (Nl % 2 == 0) ? Nl : Nl - 1;  
   int Nc1 = (Nc % 2 == 0) ? Nc : Nc - 1;
   
   // low resolution band
   for (i = 0; i < Nls_L; i++)
   for (j = 0; j < Ncs_L; j++) Ima_Aux1[2*i*Nc+2*j] = Imag[i*Nc+j];
   
   // horizontal band
   for (i = 0; i < Nls_H; i++)
   for (j = 0; j < Ncs_H; j++) Ima_Aux1[2*i*Nc+2*j+1] = Imag[i*Nc+j+Ncs_L];
   
   // vertical band
   for (i = 0; i < Nls_V; i++)
   for (j = 0; j < Ncs_V; j++) Ima_Aux1[(2*i+1)*Nc+2*j] = Imag[(i+Nls_L)*Nc+j];
   
   // diagonal band
   for (i = 0; i < Nls_D; i++)
   for (j = 0; j < Ncs_D; j++) 
              Ima_Aux1[(2*i+1)*Nc+2*j+1] = Imag[(i+Nls_L)*Nc+j+Ncs_L];
 
   for (i = 0; i < Nl1; i+=2)
   for (j = 0; j < Nc; j++)
   {
        if (Ima_Aux1[(i+1)*Nc+j] > 0.) 
        {
            Ima_Aux[i*Nc+j] = Ima_Aux1[i*Nc+j] + Ima_Aux1[(i+1)*Nc+j];
            Ima_Aux[(i+1)*Nc+j] = Ima_Aux1[i*Nc+j];
        }
        else 
        {
            Ima_Aux[i*Nc+j] = Ima_Aux1[i*Nc+j];
            Ima_Aux[(i+1)*Nc+j] = Ima_Aux1[i*Nc+j] - Ima_Aux1[(i+1)*Nc+j];
        }
   }
   // test non square image
   if (Nl1 != Nl)
   {
      i = Nl-1;
      for (j = 0; j < Nc; j++) Ima_Aux[i*Nc+j] =  Ima_Aux1[i*Nc+j];
   }
   
   
   for (i = 0; i < Nl; i++)
   {
       for (j = 0; j < Nc1-1; j+=2)
       {
           if (Ima_Aux[i*Nc+j+1] > 0.) 
           {
               Imag_Out[i*Nc+j] = Ima_Aux[i*Nc+j] + Ima_Aux[i*Nc+j+1];
               Imag_Out[i*Nc+j+1] = Ima_Aux[i*Nc+j];
           }
           else
           {
               Imag_Out[i*Nc+j] = Ima_Aux[i*Nc+j];
               Imag_Out[i*Nc+j+1] = Ima_Aux[i*Nc+j] - Ima_Aux[i*Nc+j+1];
           }
       }
   }
   // test non square image
   if (Nc1 != Nc)
   {
      j = Nc-1;
      for (i = 0; i < Nl; i++)  Imag_Out[i*Nc+j] = Ima_Aux[i*Nc+j];
   }
   delete [] Ima_Aux;
   delete [] Ima_Aux1;
}

/***************************************************************************/

void inverse_transform_g (Iint &Imag, Iint &Result, int Nbr_Iter)
{
   int i,j,k, Nl1, Nc1;
   int *Ima_Aux, *Imag_Out;
   int Nl = Imag.nl();
   int Nc = Imag.nc();
   int *TabNl = new int [Nbr_Iter+1];
   int *TabNc = new int [Nbr_Iter+1];
   
   Ima_Aux = new int [Nl*Nc];
   Imag_Out = new int [Nl*Nc];

   // set the image size of each size
   Nl1 = Nl;
   Nc1 = Nc;
   TabNl[0] = Nl;
   TabNc[0] = Nc;
   for (k = 0; k < Nbr_Iter; k++)
   {
       TabNl[k+1] = (TabNl[k]+1)/2;
       TabNc[k+1] = (TabNc[k]+1)/2;
   }

   
   Nl1 = TabNl[Nbr_Iter];
   Nc1 = TabNc[Nbr_Iter];
   //cout << "Nl1 = " << Nl1 << " Nc1 = " << Nc1 << endl;
   
   // put the last smoot array in Imag_Out
   for (i = 0; i < Nl1; i++)
   for (j = 0; j < Nc1; j++)
                Imag_Out[i*Nc1+j] = Imag(i,j);

   for (k = Nbr_Iter-1 ; k >= 0; k--)
   {
       Nl1 = TabNl[k];
       Nc1 = TabNc[k];
       // put the multiresolution coefficients in Ima_Aux
       for (i = 0; i < Nl1; i++)
       for (j = 0; j < Nc1; j++)  Ima_Aux[i*Nc1+j] = Imag(i,j);

       // copy the reconstructued last smooth array in Ima_Aux
       for (i = 0; i < TabNl[k+1];i++)
       for (j = 0; j < TabNc[k+1];j++) 
                   Ima_Aux[i*Nc1+j] = Imag_Out[i*TabNc[k+1]+j];

       transfo_g_inv_one_scale (Ima_Aux, Imag_Out, Nl1, Nc1);
 // {
//     for (i = 0; i < Nl1; i++)
//    for (j = 0; j < Nc1; j++)
//                 Result(i,j) = Imag_Out[i*Nc1+j];
//     char ch[80];
//     sprintf(ch,"xx_iter%d", k+1);
//     io_write_ima_int(ch,Result );
//  }
   }

   for (i = 0; i < Nl; i++)
   for (j = 0; j < Nc; j++)
                Result(i,j) = Imag_Out[i*Nc+j];

   delete [] Ima_Aux;
   delete [] Imag_Out;
   delete [] TabNl;
   delete [] TabNc;
}

/***************************************************************************/

    
