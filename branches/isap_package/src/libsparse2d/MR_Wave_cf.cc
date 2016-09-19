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
**    Date:  96/06/13 
**    
**    File:  MR_Wave_cf.cc
**
*******************************************************************************
**
**    DESCRIPTION   routines concerting filters used by the wavelet 
**    -----------   transform which need the FFT   
**
******************************************************************************
**
** float pyr_2d_cf_scaling_function (u, v, Freq_Coup, Nl, Nc)
** float u,v;
** float Freq_Coup;
** int Nl, Nc;
**
**    Computes the theorical transfert function of an instrument
**    which has a  cut off frequency Freq_Coup (0 <Freq_Coup <= 0.5)
**    The assumption is that the transfert function is a b3-spline
**
**    u,v = INPUT frequency point
**    Freq_Coup = INPUT cut off frequency  
**    Nl,Nc = INPUT image size
**
**    Return the value in u,v
**    u and v must verify : - Nl /2 < u and v < Nl / 2
**
******************************************************************************* 
**
** float pyr_2d_cf_filter_h (Freq_u, Freq_v, Freq_Coup, Nl, Nc)
** float Freq_u, Freq_v;
** float Freq_Coup;
** int Nl, Nc;
**
**   computes the coefficients which allow to go from a resolution to
**   the next one.
**
**   Freq_u, Freq_v = INPUT frequency point
**   Freq_Coup = INPUT cut-off frequency
**   Nl,Nc = INPUT image size
**
**   Return the value of the filter h in Freq_u, Freq_v 
**   Freq_u and Freq_v must verify : - Nl /2 < Freq_u and Freq_v < Nl / 2
**
******************************************************************************* 
**
** float pyr_2d_cf_filter_h_tilde (Freq_u, Freq_v, Freq_Coup, 
**                                 Nl, Nc, Type_Wavelet)
** float Freq_u, Freq_v;
** float Freq_Coup;
** int Nl, Nc;
** type_transform Type_Wavelet;
**
**  computes the coefficients which allow to go from a resolution to
**   the previous  one.
**   Freq_u, Freq_v = INPUT frequency point
**   Freq_Coup = INPUT cut-off frequency
**   Nl,Nc = INPUT image size
**
**   Return the value of the filter h_tilde in Freq_u, Freq_v 
**   Freq_u and Freq_v must verify : - Nl /2 < Freq_u and Freq_v < Nl / 2
**
******************************************************************************* 
**
** float pyr_2d_cf_filter_g (Freq_u, Freq_v, Freq_Coup, Nl, Nc, Type_Wavelet)
** float Freq_u, Freq_v;
** float Freq_Coup;
** int Nl, Nc;
** type_transform Type_Wavelet;
**
**   computes the coefficients which allow to go from a resolution to
**   the next one.
**
**   Freq_u, Freq_v = INPUT frequency point
**   Freq_Coup = INPUT cut-off frequency
**   Nl,Nc = INPUT image size
**
**   Return the value of the filter g in Freq_u, Freq_v 
**   Freq_u and Freq_v must verify : - Nl /2 < Freq_u and Freq_v < Nl / 2
**  Type_Wavelet =   TO_PYR_FFT_DIFF_RESOL
**             or    TO_PYR_FFT_DIFF_SQUARE_RESOL
**             or    TO_PAVE_FFT
**
******************************************************************************* 
**
** float pyr_2d_cf_filter_g_tilde (Freq_u, Freq_v, Freq_Coup, 
**                                 Nl, Nc, Type_Wavelet)
** float Freq_u, Freq_v;
** float Freq_Coup;
** int Nl, Nc
** type_transform Type_Wavelet;
**
**   computes the coefficients which allow to go from a resolution to
**   the prevuois one.
**
**   Freq_u, Freq_v = INPUT frequency point
**   Freq_Coup = INPUT cut-off frequency
**   Nl,Nc = INPUT image size
**
**   Return the value of the filter g_tilde in Freq_u, Freq_v 
**   Freq_u and Freq_v must verify : - Nl /2 < Freq_u and Freq_v < Nl / 2
**  Type_Wavelet =   TO_PYR_FFT_DIFF_RESOL
**             or    TO_PYR_FFT_DIFF_SQUARE_RESOL
**             or    TO_PAVE_FFT
**
******************************************************************************* 
**
** float pyr_2d_cf_filter_wavelet (Freq_u, Freq_v, Freq_Coup, 
**                                 Nl, Nc, Type_Wavelet)
** float Freq_u, Freq_v;
** float Freq_Coup;
** int Nl, Nc;
** type_transform Type_Wavelet;
**
**   computes the wavelet in Freq_u, Freq_v
**
**   Freq_u, Freq_v = INPUT frequency point
**   Freq_Coup = INPUT cut-off frequency
**   Nl,Nc = INPUT image size
**
**   Return the value of the wavelet in Freq_u, Freq_v 
**   Freq_u and Freq_v must verify : - Nl /2 < Freq_u and Freq_v < Nl / 2
**  Type_Wavelet =   TO_PYR_FFT_DIFF_RESOL
**             or    TO_PYR_FFT_DIFF_SQUARE_RESOL
**             or    TO_PAVE_FFT
**
******************************************************************************* 
**
**  float pyr_2d_cf_filter(Which_Filter,Freq_u,Freq_v,Fc, Nl, Nc, Type_Wavelet)
**  int Which_Filter;
**  float Freq_u, Freq_v;
**  float Fc;
**  int Nl, Nc;
**  type_transform Type_Wavelet;
**
**   computes the value in Freq_u, Freq_v of the filter defined by Which_Filter  
**   Which_Filter  = INPUT type filter 
**                Which_Filter =  SCALING_FUNCTION
**                             or FILTER_H
**                             or FILTER_H_TILDE
**                             or FILTER_G
**                             or FILTER_G_TILDE
**                             or WAVELET
**   Freq_u, Freq_v = INPUT frequency point
**   Freq_Coup = INPUT cut-off frequency
**   Nl,Nc = INPUT image size
**
**   Return the value of the wavelet in Freq_u, Freq_v 
**   Freq_u and Freq_v must verify : - Nl /2 < Freq_u and Freq_v < Nl / 2
**
**  Type_Wavelet =   TO_PYR_FFT_DIFF_RESOL
**             or    TO_PYR_FFT_DIFF_SQUARE_RESOL
**             or    TO_PAVE_FFT
**
******************************************************************************* 
**
**  void pyr_2d_cf_create_filter (Ifloat &Filter, float Fc, 
**                                int Which_Filter, type_transform Type_Wavelet)
**
**  computes the  filter defined by Which_Filter  
**
**   Which_Filter  = INPUT type filter 
**                Which_Filter =  SCALING_FUNCTION
**                             or FILTER_H
**                             or FILTER_H_TILDE
**                             or FILTER_G
**                             or FILTER_G_TILDE
**                             or WAVELET
**   Freq_Coup = INPUT cut-off frequency
**   Nl,Nc = INPUT image size
**
**   Filter = OUTPUT filter
**
**  Type_Wavelet =   TO_PYR_FFT_DIFF_RESOL
**             or    TO_PYR_FFT_DIFF_SQUARE_RESOL
**             or    TO_PAVE_FFT
**
******************************************************************************* 
**
** void wave_cf_transform (Ifloat &Image, MultiResol &MR_Transf)
**
** computes the wavelet transform in the Fourier space of an image
**
******************************************************************************* 
**
** void wave_cf_recons (MultiResol &MR_Transf, Ifloat &Image)
**
** reconstruct the image from its wavelet transform
**
******************************************************************************* 
**
** void mr_cf_transform (Ifloat &Image, MultiResol &MR_Transf)
** 
** computes the multiresolution transform of an image
** (the scales contains smoothed images, and no details)
**
******************************************************************************/ 

// static char sccsid[] = "@(#)MR_Wave_cf.cc 3.2 96/06/13 CEA 1994 @(#)";

#include <math.h>
#include <stdio.h>
#include <string.h> 

#include "IM_Obj.h"
#include "MR_Obj.h"
#include "IM_IO.h"

static Icomplex_f Ima_h_cf;
static Icomplex_f Ima_b_cf;
static Icomplex_f Buff;

#define DEBUG 1

float pyr_2d_cf_scaling_function (float, float, float, int);
float pyr_2d_cf_filter_h(float, float, float, int, int);
float pyr_2d_cf_filter_h_tilde(float, float, float,int,int, type_transform);
float pyr_2d_cf_filter_g(float, float, float, int, int, type_transform);
float pyr_2d_cf_filter_g_tilde(float, float, float,int,int, type_transform);
float pyr_2d_cf_filter_wavelet (float, float,float,int, type_transform);
float pyr_2d_cf_filter (int Which_Filter, float Freq_u, float Freq_v, 
                        float Fc, int Nl, int Nc, type_transform Type_Wavelet);
void pyr_2d_cf_create_filter (Ifloat &, float, int, type_transform);

/***************************************************************************/

void print_info_cf (Icomplex_f &Image, int s)
{
    int Nl = Image.nl();
    int Nc = Image.nc();
    char Buff_Name[100];
    sprintf (Buff_Name, "%s (%d, %d)", "FFT Ima_Buff", Nl, Nc);
    Ifloat Ima_Buff(Nl, Nc, Buff_Name);

    //Ima_Buff = real(Image);
    real(Ima_Buff, Image);
    cout << "       Sigma reel = " << Ima_Buff.sigma () << endl;
    cout << "       Min reel = " << Ima_Buff.min () << endl;
    cout << "       Max reel = " << Ima_Buff.max () << endl;
    cout << "       Moy reel = " << Ima_Buff.mean() << endl;
}

/****************************************************************************/

float pyr_2d_cf_scaling_function (float u, float v, float Freq_Coup, int Nl)
{
    float F, Fc, Scaling_Function;
 
    /* we fixe the cut-off frequency to N / 2 */
    Fc = Freq_Coup * (float) Nl;

    F = sqrt ((float)(u * u + v * v));
    /*  . we multiply by 3/2 in order that: T(0) = 1
        . we multiply F by 2./Fc (b3_spline(2) = 0) 
          in order that: T(Fc) = 0
    */
    Scaling_Function = 3. / 2. * b3_spline ((double)(2. * F / Fc));

    return (Scaling_Function);
}

/***************************************************************************/

float pyr_2d_cf_filter_h (float Freq_u, float Freq_v, float Freq_Coup,
                        int Nl, int Nc)
{
    float Filter_H;
    int Nl_2, Nc_2;
    float u2,v2;
    float Freq1,Freq2;

    Nl_2 = Nl / 2;
    Nc_2 = Nc / 2;
    u2 = 2. * Freq_u;
    v2 = 2. * Freq_v;

    if ((u2 < - Nl_2 ) || (u2 >= Nl_2) || (v2 < - Nc_2) || (v2 >= Nc_2))
    {
         Filter_H = 0.;
    }
    else
    {
         Freq1 = pyr_2d_cf_scaling_function (Freq_u, Freq_v, Freq_Coup, Nl);
         Freq2 = pyr_2d_cf_scaling_function (u2, v2, Freq_Coup, Nl);  

         if (fabs(Freq1) < FLOAT_EPSILON)
         {
             Filter_H = 0.;
         }
         else 
         {
             Filter_H = Freq2 / Freq1;
         }
    }
    return (Filter_H);
}

/***************************************************************************/   

float pyr_2d_cf_filter_g (float Freq_u, float Freq_v, float Freq_Coup, 
                          int Nl, int Nc, type_transform Type_Wavelet)
{
    float Filter_G=0.;
    float Filter_H;

     Filter_H = pyr_2d_cf_filter_h (Freq_u, Freq_v, Freq_Coup, Nl, Nc);
     switch (Type_Wavelet)
     {
         case TO_PYR_FFT_DIFF_RESOL:
         case TO_PAVE_FFT: 
                      Filter_G = 1 - Filter_H;
                      break;
         case TO_PYR_FFT_DIFF_SQUARE:
                      Filter_G  = sqrt (1. - Filter_H * Filter_H);
                      break;
        default:
            fprintf (stderr, "Error: bad wave in pyr_2d_cf_filter_g\n");
            exit (-1);
            break;
     }

    return (Filter_G);
}

/***************************************************************************/

float pyr_2d_cf_filter_h_tilde (float Freq_u, float Freq_v, float Freq_Coup, 
                                int Nl, int Nc, type_transform Type_Wavelet)
{
    float Filter_H_Tilde=0, Filter_H, Filter_G, Den;


    switch (Type_Wavelet)
    {
        case TO_PYR_FFT_DIFF_RESOL: 
        case TO_PAVE_FFT:
                Filter_H = pyr_2d_cf_filter_h (Freq_u, Freq_v, Freq_Coup, 
                                              Nl, Nc);
                Filter_G = pyr_2d_cf_filter_g (Freq_u, Freq_v, Freq_Coup, 
                                                         Nl, Nc,Type_Wavelet);
                Den = Filter_H * Filter_H + Filter_G * Filter_G;
                if (Den < FLOAT_EPSILON) Filter_H_Tilde = 0.;
                else Filter_H_Tilde = Filter_H / Den;
            break;
        case TO_PYR_FFT_DIFF_SQUARE:
            Filter_H_Tilde = pyr_2d_cf_filter_h(Freq_u, Freq_v, Freq_Coup, 
                                                Nl, Nc);
            break;
        default:
            fprintf (stderr, "Error: bad wave in pyr_2d_cf_filter_h_tilde\n");
            exit (-1);
            break;
    }
    return (Filter_H_Tilde);
}

/***************************************************************************/

float pyr_2d_cf_filter_g_tilde (float Freq_u, float Freq_v, float Freq_Coup, 
                                int Nl, int Nc, type_transform Type_Wavelet)
{
    float Filter_G_Tilde=0, Filter_H, Filter_G, Den;


    switch (Type_Wavelet)
    {
        case TO_PYR_FFT_DIFF_RESOL:
        case TO_PAVE_FFT:
                Filter_H = pyr_2d_cf_filter_h(Freq_u, Freq_v, Freq_Coup, 
                                              Nl, Nc);
                Filter_G = pyr_2d_cf_filter_g(Freq_u, Freq_v, 
                                           Freq_Coup, Nl, Nc,Type_Wavelet);
                Den = Filter_H * Filter_H + Filter_G * Filter_G;
                if (Den < FLOAT_EPSILON) Filter_G_Tilde = 0.;
                else Filter_G_Tilde = Filter_G / Den;
                break;
        case TO_PYR_FFT_DIFF_SQUARE:
                Filter_G_Tilde = pyr_2d_cf_filter_g(Freq_u, Freq_v, Freq_Coup, 
                                              Nl, Nc, Type_Wavelet);
            break;
        default:
            fprintf (stderr, "Error: bad wave in pyr_2d_cf_filter_g_tilde\n");
            exit (-1);
            break;
    }
    return (Filter_G_Tilde);
}
   
/***************************************************************************/   


float pyr_2d_cf_filter_wavelet (float Freq_u, float Freq_v, float Freq_Coup, 
                                int Nl, type_transform Type_Wavelet)
{
    float Filter_Psi=0;
    float u2,v2;
    float Freq1,Freq2;


    Freq1 = pyr_2d_cf_scaling_function (Freq_u, Freq_v, Freq_Coup, Nl);
    u2 = Freq_u / 2.;
    v2 = Freq_v / 2.;
    Freq2 = pyr_2d_cf_scaling_function (u2, v2, Freq_Coup, Nl);

    switch (Type_Wavelet)
    {
        case TO_PYR_FFT_DIFF_RESOL: 
        case TO_PAVE_FFT:
            Filter_Psi = Freq2 - Freq1;
            break;
        case TO_PYR_FFT_DIFF_SQUARE:
            Filter_Psi = Freq2*Freq2 - Freq1*Freq1;
            break;
        default:
            fprintf (stderr, "Error: bad wave in pyr_2d_cf_filter_wavelet\n");
            exit (-1);
            break;
    }
    return (Filter_Psi);
}

/***************************************************************************/

float pyr_2d_cf_filter (int Which_Filter, float Freq_u, float Freq_v, 
                        float Fc, int Nl, int Nc, type_transform Type_Wavelet)
{
    float Val_Return=0.;

    switch (Which_Filter)
    {
        case SCALING_FUNCTION :
              Val_Return = pyr_2d_cf_scaling_function (Freq_u, Freq_v, Fc, Nl);
              break;
        case FILTER_H :
              Val_Return = pyr_2d_cf_filter_h (Freq_u, Freq_v, Fc, Nl, Nc);
              break;
        case FILTER_H_TILDE :
              Val_Return =  pyr_2d_cf_filter_h_tilde (Freq_u, Freq_v, Fc, 
                                                   Nl, Nc, Type_Wavelet);
              break;
        case FILTER_G :
              Val_Return = pyr_2d_cf_filter_g (Freq_u, Freq_v, Fc, 
                                                   Nl, Nc, Type_Wavelet);
              break;
        case FILTER_G_TILDE :
              Val_Return = pyr_2d_cf_filter_g_tilde (Freq_u, Freq_v, Fc, 
                                                   Nl, Nc, Type_Wavelet);
              break;
        case WAVELET :
              Val_Return = pyr_2d_cf_filter_wavelet (Freq_u, Freq_v, Fc, 
                                                   Nl, Type_Wavelet);
              break;
        default:
              fprintf (stderr, "Error: bad filter in pyr_2d_cf_filter\n");
              exit (0);
              break;
    }
    assert ((Val_Return >= -FLOAT_EPSILON) && (Val_Return < 2));
    return (Val_Return);
}

/***************************************************************************/

void pyr_2d_cf_create_filter (Ifloat &Filter, float Fc, 
                              int Which, type_transform Type_Wavelet)
{
    int i,j;
    float u,v;
    int Nl = Filter.nl();
    int Nc = Filter.nc();

    for (i = 0; i < Nl; i++)
    {
        u = (float) i - (float) Nl / 2.;
        for (j = 0; j < Nc; j++)
        {
            v = (float) j - (float) Nc / 2.;
            Filter (i,j) = pyr_2d_cf_filter (Which, u, v, Fc, 
                                                 Nl, Nc, Type_Wavelet);
        }
    }
}

/***************************************************************************/

static void wave_cf_mult (Icomplex_f &Ima_b, Icomplex_f &Ima_h, int s, int N,
                          float Fc, type_transform Type_Transform, Bool Down)
{
    int Nl = N;
    int Nc = N;
    int Nl_s = Ima_b.nl();
    int Nc_s = Ima_b.nc();
    int i,j;
    int u,v;
    int Dep;

    Dep = (int) (pow((double)2., (double) s) + 0.5);

    for (i = 0; i < Nl_s; i++)
    for (j = 0; j < Nc_s; j++)
    {
        u = Dep * (i - Nl_s/2);
        v = Dep * (j - Nc_s/2);
        if ((u+Nl/2 < 0) || (u+Nl/2 >= Nl) || (v+Nc/2 < 0) || (v+Nc/2 >= Nc))
        {
           Ima_b(i,j) = complex_f ((float)0.,(float)0.);
        }
        else
        {
            if (Down)
            {
                 Ima_b (i,j) *= pyr_2d_cf_filter (FILTER_H, (float) u, (float)v,
                                              Fc, Nl, Nc, Type_Transform);
                 Ima_h (i,j) *= pyr_2d_cf_filter (FILTER_G, (float) u, (float)v,
                                                 Fc, Nl, Nc, Type_Transform);
            }
            else
            {
                 Ima_b (i,j) *= pyr_2d_cf_filter (FILTER_H_TILDE, (float) u,
                                                  (float) v, Fc, 
                                                  Nl, Nc, Type_Transform);
                 Ima_h (i,j) *= pyr_2d_cf_filter (FILTER_G_TILDE, (float) u,
                                                  (float) v, Fc, 
                                                  Nl, Nc, Type_Transform);
            }

        }
    }
}

/****************************************************************************/

static void under_sampling_2 (Icomplex_f &Imag)
{
    int Nl = Imag.nl();
    int Nc = Imag.nc();
    int Nl_2 = Nl/2;
    int Nc_2 = Nc/2;
    int i1,j1,i0,j0, Depi, Depj;

    Buff.resize (Nl_2, Nc_2);

    Depi = Nl/4;
    Depj = Nl/4;

    for (i0 = 0; i0 < Nl_2; i0++)
    for (j0 = 0; j0 < Nc_2; j0++)
    {
       i1 = i0 + Depi;
       j1 = j0 + Depj;
       Buff (i0,j0) = Imag (i1,j1);
    }
    Imag.resize(Nl_2, Nc_2);
    Imag = Buff;
}

/****************************************************************************/

static void over_sampling_2 (Icomplex_f &Imag)
{
    int Nl = Imag.nl();
    int Nc = Imag.nc();
    int Nl_2 = 2*Nl;
    int Nc_2 = 2*Nc;
    int i1,j1,i0,j0, Depi, Depj;

    Buff.resize (Nl_2, Nc_2);

    Depi = Nl_2/4;
    Depj = Nc_2/4;

    for (i0 = 0; i0 < Nl_2; i0++)
    for (j0 = 0; j0 < Nc_2; j0++) Buff (i0,j0) = 0;

    for (i1 = 0; i1 < Nl; i1++)
    for (j1 = 0; j1 < Nc; j1++) 
    {
       i0 = i1 + Depi;
       j0 = j1 + Depj;
       Buff (i0,j0) = Imag (i1,j1);
    }
    Imag.resize(Nl_2, Nc_2);
    Imag = Buff;
}

/****************************************************************************/

static void normalize_scale (Ifloat &Scale, int N, int Dir)
{
    float Coef;

    /* Computes the normalisation term */
    if (Dir == -1) Coef = (float)(Scale.nl() * Scale.nc()) / (float)(N * N);
    else Coef = (float)(N * N) / (float)(Scale.nl() * Scale.nc());

    for (int i = 0; i < Scale.nl(); i++)
    for (int j = 0; j < Scale.nc(); j++) Scale (i,j) *= Coef;
}

/****************************************************************************/

void wave_cf_transform (Ifloat &Image, MultiResol &MR_Transf)
{
    int Nl = Image.nl();
    int Nc = Image.nc();
    int Nl_s, Nc_s, s;
    int Nbr_Plan = MR_Transf.nbr_scale();
    float Fc = MR_Transf.Fc;

    Ima_h_cf.resize (Nl, Nc);
    Ima_b_cf.resize (Nl, Nc);
    Buff.resize (Nl, Nc);

    switch (MR_Transf.Type_Transform)
    {
        case TO_PAVE_FFT:
            fft2d (Image, Ima_b_cf, 1);
            for (s = 0; s < Nbr_Plan - 1; s ++)
            {
                Nl_s = MR_Transf.band(s).nl();
                Nc_s = MR_Transf.band(s).nc();
                Ima_h_cf = Ima_b_cf;
                wave_cf_mult (Ima_b_cf, Ima_h_cf, s, Nl, Fc, 
                              MR_Transf.Type_Transform, True);
                fft2d (Ima_h_cf, -1);
		//MR_Transf.scale(s) = real( invfft (Ima_h_cf));
		real( MR_Transf.band(s), Ima_h_cf);
            }
            //MR_Transf.scale(Nbr_Plan - 1) = real( invfft (Ima_b_cf));
	    fft2d (Ima_b_cf, -1);
	    real(MR_Transf.scale(Nbr_Plan - 1), Ima_b_cf);
            break;
        case TO_PYR_FFT_DIFF_RESOL:
        case TO_PYR_FFT_DIFF_SQUARE:
            fft2d (Image, Ima_b_cf, 1);
            for (s = 0; s < Nbr_Plan - 1; s ++)
            {
                Nl_s = MR_Transf.band(s).nl();
                Nc_s = MR_Transf.band(s).nc();
                Ima_h_cf.resize (Nl_s, Nc_s);
                Ima_h_cf = Ima_b_cf;
                wave_cf_mult (Ima_b_cf, Ima_h_cf, s, Nl, Fc, 
                              MR_Transf.Type_Transform, True);
                
                //MR_Transf.scale(s) = real( invfft (Ima_h_cf));
		fft2d (Ima_h_cf, -1);
		real(MR_Transf.band(s), Ima_h_cf);
                normalize_scale (MR_Transf.band(s), Nl, -1);
                under_sampling_2 (Ima_b_cf);
            }
            //MR_Transf.scale(Nbr_Plan - 1) = real( invfft (Ima_b_cf));
            fft2d (Ima_b_cf, -1);
	    real(MR_Transf.scale(Nbr_Plan - 1), Ima_b_cf);
            normalize_scale (MR_Transf.scale(Nbr_Plan - 1), Nl, -1);
            break;
        default:
            fprintf (stderr,"Error: proc. wave_cf_transform.\n");
            fprintf (stderr,"This transform is not computed this procedure\n");
            break;

    }
}

/****************************************************************************/

void wave_cf_recons (MultiResol &MR_Transf, Ifloat &Image)
{
    int Nl = Image.nl();
    int Nc = Image.nc();
    int Nl_s, Nc_s, s;
    int Nbr_Plan = MR_Transf.nbr_scale();
    float Fc = MR_Transf.Fc;

    Ima_h_cf.resize (Nl, Nc);
    Ima_b_cf.resize (Nl, Nc);
    Buff.resize (Nl, Nc);

    switch (MR_Transf.Type_Transform)
    {
        case TO_PAVE_FFT:
            Image = MR_Transf.scale(Nbr_Plan-1);
            for (s = Nbr_Plan -2; s >= 0; s--) 
                          Image += MR_Transf.scale(s);
            break;
        case TO_PYR_FFT_DIFF_RESOL:
        case TO_PYR_FFT_DIFF_SQUARE:
            s = Nbr_Plan - 1;
            Nl_s = MR_Transf.scale(s).nl();
            Nc_s = MR_Transf.scale(s).nc();
            Ima_b_cf.resize (Nl_s, Nc_s);
            normalize_scale (MR_Transf.scale(s), Nl, 1);
            fft2d(MR_Transf.scale(s), Ima_b_cf);
            normalize_scale (MR_Transf.scale(s), Nl, -1);
            for (s = Nbr_Plan - 2; s >= 0; s --)
            {
                Nl_s = MR_Transf.scale(s).nl();
                Nc_s = MR_Transf.scale(s).nc();
                Ima_h_cf.resize (Nl_s, Nc_s);
                normalize_scale (MR_Transf.scale(s), Nl, 1);
                fft2d(MR_Transf.scale(s), Ima_h_cf);

                over_sampling_2 (Ima_b_cf);
                wave_cf_mult (Ima_b_cf, Ima_h_cf, s, Nl, Fc, 
                              MR_Transf.Type_Transform, False);
                Ima_b_cf += Ima_h_cf;
                normalize_scale (MR_Transf.scale(s), Nl, -1);
            }
            //Image = real( invfft (Ima_b_cf));
	    fft2d(Ima_b_cf, -1);
	    real(Image, Ima_b_cf);
            break;
        default:
            fprintf (stderr,"Error: proc. wave_cf_transform.\n");
            fprintf (stderr,"This transform is not computed this procedure\n");
            break;

    }
}

/****************************************************************************/

void mr_cf_transform (Ifloat &Image, MultiResol &MR_Transf)
{
    int Nl = Image.nl();
    int Nc = Image.nc();
    int Nl_s, Nc_s, s;
    int Nbr_Plan = MR_Transf.nbr_scale();
    float Fc = MR_Transf.Fc;

    Ima_h_cf.resize (Nl, Nc);
    Ima_b_cf.resize (Nl, Nc);
    Buff.resize (Nl, Nc);

    switch (MR_Transf.Type_Transform)
    {
        case TO_PAVE_FFT:
            fft2d (Image, Ima_b_cf, 1);
            for (s = 0; s < Nbr_Plan - 1; s ++)
            {
                Nl_s = MR_Transf.band(s).nl();
                Nc_s = MR_Transf.band(s).nc();
                Ima_h_cf = Ima_b_cf;
                //MR_Transf.scale(s) = real( invfft (Ima_h_cf));
		fft2d (Ima_h_cf, -1);
		real(MR_Transf.band(s),Ima_h_cf);
                wave_cf_mult (Ima_b_cf, Ima_h_cf, s, Nl, Fc, 
                              MR_Transf.Type_Transform, True);
            }
            //MR_Transf.scale(Nbr_Plan - 1) = real( invfft (Ima_b_cf));
	    fft2d (Ima_b_cf, -1);
	    real( MR_Transf.band(Nbr_Plan - 1), Ima_b_cf);
            break;
        case TO_PYR_FFT_DIFF_RESOL:
        case TO_PYR_FFT_DIFF_SQUARE:
            fft2d (Image, Ima_b_cf, 1);
            for (s = 0; s < Nbr_Plan - 1; s ++)
            {
                Nl_s = MR_Transf.band(s).nl();
                Nc_s = MR_Transf.band(s).nc();
                Ima_h_cf.resize (Nl_s, Nc_s);
                Ima_h_cf = Ima_b_cf;
                //MR_Transf.scale(s) = real( invfft (Ima_h_cf));
		fft2d (Ima_h_cf, -1);
		real(MR_Transf.scale(s), Ima_h_cf);
                normalize_scale (MR_Transf.band(s), Nl, -1);
                wave_cf_mult (Ima_b_cf, Ima_h_cf, s, Nl, Fc, 
                              MR_Transf.Type_Transform, True);
                under_sampling_2 (Ima_b_cf);
            }
            //MR_Transf.scale(Nbr_Plan - 1) = real( invfft (Ima_b_cf));
	    fft2d (Ima_b_cf, -1);
	    real(MR_Transf.band(Nbr_Plan - 1), Ima_b_cf);
            normalize_scale (MR_Transf.scale(Nbr_Plan - 1), Nl, -1);

            break;
        default:
            fprintf (stderr,"Error: proc. mr_cf_transform.\n");
            fprintf (stderr,"This transform is not computed this procedure\n");
            break;
    }
}

/****************************************************************************/

