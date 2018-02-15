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
**    File:  MR1D_Transf.cc
**
*****************************************************************************
**
**    DESCRIPTION 1D wavelet transform and reconstruction
**    -----------    
**   
*****************************************************************************
** 
** wave_1d_mex (fltarray & Signal, fltarray & W_1D, int N, 
**             int Nbr_Voie, int Nbr_Plan, float Scale_0)
**
** wavelet transform by the mexican hat
** Signal = IN: data
** W_1D = OUT: wavelet transform
** N = IN: number of data in the Signal
** Nbr_Voie = IN: number of channels per octave (in general, we take 12)
** Nbr_Plan = OUT: number of scales
** Scale_0 = OUT: first scale
**
** the wavelet transform W_1D[0..Nbr_Plan-1][0..N] is computed
** W_1D[i] defines the scale i 
**   
****************************************************************************
**
** wave_1d_mex_rec (fltarray &W_1D, fltarray &Signal, int N, 
**                     int Nbr_Voie, int Nbr_Plan, float Scale_0)
**
** wavelet reconstruction from a mexican hat wavelet transform
**
** W_1D = IN: wavelet transform
** Signal = OUT: recosntructed signal
** N = IN: number of data in the Signal
** Nbr_Voie = IN: number of channels per octave (in general, we take 12)
** Nbr_Plan = IN: number of scales
**   
**************************************************************************** 
**   
** wave_1d_french (fltarray & Signal, fltarray & W_1D, int N, 
**                int Nbr_Voie, int Nbr_Plan, float Scale_0)
**
** wavelet transform by the french hat
**
** Signal = IN: data
** W_1D = OUT: wavelet transform
** N = IN: number of data in the Signal
** Nbr_Voie = IN: number of channels per octave (in general, we take 12)
** Nbr_Plan = OUT: number of scales
** Scale_0 = OUT: first scale
**
** the wavelet transform W_1D[0..Nbr_Plan-1][0..N] is computed
** W_1D[i] defines the scale i 
**
**************************************************************************** 
**   
** wave_1d_french_rec (fltarray & W_1D, fltarray & Signal, int N, 
**                     int Nbr_Voie, int Nbr_Plan, float Scale_0)
**
** reconstruction from the wavelet transform with the french hat 
**
** W_1D = IN: wavelet transform
** Signal = OUT: recosntructed signal
** N = IN: number of data in the Signal
** Nbr_Voie = IN: number of channels per octave (in general, we take 12)
** Nbr_Plan = IN: number of scales
**
***************************************************************************
**
** wave_1d_morlet (fltarray & Signal, fltarray & W_1D, int N, 
**                int Nbr_Voie, int Nbr_Plan, float Nu_0, float Scale0)
**
**  Morlet transform (Nu_0 > 0.8)
**
** Signal = IN: data
** W_1D = OUT: real part of wavelet transform from 0 to N
**             imaginary part of wavelet transform from N to 2*N
** N = IN: number of data in the Signal
** Nbr_Voie = IN: number of channels per octave (in general, we take 12)
** Nbr_Plan = OUT: number of scales
** Nu_0 = IN: Morlet's parameter (Nu_0 must be > 0.8)  
** Scale_0 = OUT: first scale
**
** the wavelet transform W_1D(0..Nbr_Plan-1,0..N) is computed
** W_1D(i) defines the scale i 
** 
****************************************************************************
** 
** wave1d_transform (fltarray & Signal, int Np, type_trans_1d Type_Transform, 
**               int Nbr_Voie,  fltarray & W_1D, int Nbr_Plan, 
**               float Scale_0, float Nu_0)
**
** 1D wavelet transform
**
** Signal = IN: data
** W_1D = OUT: wavelet transform
** Nc = IN: number of data in the Signal
** Nbr_Voie = IN: number of channels per octave (in general, we take 12)
** Nbr_Plan = OUT: number of scales
** Nu_0 = IN: Morlet's parameter (Nu_0 must be > 0.8)  
** Scale_0 = OUT: first scale
** Type_Trans = type of wavelet transform
**                    TO1_PAVE_LINEAR 
**                    TO1_PAVE_B1SPLINE
**                    TO1_PAVE_B3SPLINE
**                    TO1_PAVE_MORLET
**                    TO1_PAVE_MEX
**                    TO1_PAVE_FRENCH
**                    TM1_PAVE_MEDIAN
**                    TO1_PYR_LINEAR
**                    TO1_PYR_B3SPLINE
**                    TM1_PYR_MEDIAN
** 
** the wavelet transform W_1D(0..Nbr_Plan-1,0..N) is computed
** W_1D(i) defines the scale i 
** for the Morlet's transform, W_1D[i] define the real part and
** W_1D(i+Nbr_Plan) defines the imaginarty par of the transform
**
****************************************************************************
** 
** wave1d_recons (fltarray &W_1D, int Np, int Nbr_Plan, type_trans_1d T_Transf,
**                 int Nbr_Voie, fltarray & Signal, float Scale_0)
**
** W_1D = IN: wavelet transform
** Signal = OUT: reconstructed signal
** Nc = IN: number of data in the Signal
** Nbr_Voie = IN: number of channels per octave (in general, we take 12)
** Nbr_Plan = IN: number of scales
** Scale_0 = IN: first scale
** T_Transf = IN: type of  transform used
** Nu_0 =  IN:Morlet's parameter 
** 
****************************************************************************/

// static char sccsid[] = "@(#)wave1d.c 6.1.1.1 93/07/16 CEA @(#)";

#include "MR1D_Obj.h"
// #include "NR.h"
 

/***************************************************************************/

void wave_1d_mex_rec (fltarray &W_1D, fltarray &Signal, 
                     int N, type_border Border, 
                     int Nbr_Voie, int Nbr_Plan, float Scale_0)
{
    int i,j,k,l,jj;
    float D_Scale, Scale, x, Step;
    float C_psi_mex = PI;
    float tmp, Wave;

    D_Scale = pow ((double) 2., (double)(1. / (double) Nbr_Voie));
    Step = log(D_Scale);
    Scale = Scale_0;

    for (j = 0; j < N; j++) Signal(j) = 0.;

    for (i = 0; i < Nbr_Plan; i++)
    {
        k = (int)(4. * Scale);
        for (j = 0; j < N; j++)
        {
            tmp = 0.;
            for (l = j-k; l < j+k; l++)
            {
                jj = border_ind_test (l, N, Border);
                x = (float)(j - l)/ Scale;
                x *= x;
                Wave = (1. - x) * exp(-x * .5);
                tmp += Wave * W_1D(jj,i);
            }
            Signal(j) += tmp / (Scale * C_psi_mex) * Step;
        }
        Scale *= D_Scale;
    }
}

/***************************************************************************/

void wave_1d_mex (fltarray & Signal, fltarray & W_1D, int N,type_border Border, 
             int Nbr_Voie, int Nbr_Plan, float Scale_0)
/* wavelet transform by the mexican hat */
{
    int i,j,k,l,jj;
    float D_Scale, Scale, x;

//    printf ("Nbr_Plan = %d\n", Nbr_Plan);
    Scale = Scale_0;
    D_Scale = pow ((double)2., (double)(1. / (double) Nbr_Voie));

    for (i = 0; i < Nbr_Plan; i++)
    {
        k = (int) (4. * Scale);
        for (j = 0; j < N; j++)
        {
            W_1D(j,i) = 0.;
            for (l = j-k; l < j+k; l++)
            {
                jj = border_ind_test (l, N, Border);
                x = (float)(j - l)/ Scale;
                x *= x;
                W_1D(j,i) += (1. - x) * exp(-x * .5) * Signal(jj);
            }
            W_1D(j,i) /= Scale;
        }
        Scale *= D_Scale;
    }
}

/***************************************************************************/

void wave_1d_der_gauss (fltarray & Signal, fltarray & W_1D, int N,type_border Border, 
             int Nbr_Voie, int Nbr_Plan, float Scale_0)
/* wavelet transform by the mexican hat */
{
    int i,j,k,l,jj;
    float D_Scale, Scale, x;

//    printf ("Nbr_Plan = %d\n", Nbr_Plan);
    Scale = Scale_0;
    D_Scale = pow ((double)2., (double)(1. / (double) Nbr_Voie));

    for (i = 0; i < Nbr_Plan; i++)
    {
        k = (int) (4. * Scale);
        for (j = 0; j < N; j++)
        {
            W_1D(j,i) = 0.;
            for (l = j-k; l < j+k; l++)
            {
                jj = border_ind_test (l, N, Border);
                x = (float)(j - l)/ Scale;
                W_1D(j,i) -=  2. * x * exp(-x*x * .5) * Signal(jj);
            }
            W_1D(j,i) /= Scale;
        }
        Scale *= D_Scale;
    }
}

/***************************************************************************/

void wave_1d_french (fltarray & Signal, fltarray & W_1D, int N,
                 type_border Border,
                int Nbr_Voie, int Nbr_Plan, float Scale_0)
/* wavelet transform by the french hat */
{
    int i,j,k,l,k3;
    float D_Scale, Scale;

//    printf ("Nbr_Plan = %d\n", Nbr_Plan);

    D_Scale = pow ((double)2., (double)(1. / (double) Nbr_Voie));
    Scale = Scale_0;

    for (i = 0; i < Nbr_Plan; i++)
    {
        k = (int) Scale;
        k3 = (int) (3. * Scale);
        for (j = 0; j < N; j++)
        {
            W_1D(j,i) = 0.;
            for (l = j-k3; l < j-k; l++)
                W_1D(j,i) -= Signal( border_ind_test (l, N, Border) );
            for (l =  j-k; l <= j+k; l++) 
                W_1D(j,i) += 2*Signal( border_ind_test (l, N, Border) );
            for (l =  j+k+1; l <= j+k3; l++) 
                W_1D(j,i) -= Signal( border_ind_test (l, N, Border) );
            W_1D(j,i) /= Scale;
        }
        Scale *= D_Scale;
    }
}

/***************************************************************************/
 
void wave_1d_french_rec (fltarray & W_1D, fltarray & Signal, int N, 
                         type_border Border,
                         int Nbr_Voie, int Nbr_Plan, float Scale_0)
/* reconstruction */
{
    int i,j,k,l,k3;
    float D_Scale, Scale;
    float C_psi_french = 27.;
    float tmp, Step;

    Scale = Scale_0;
    D_Scale = pow ((double)2., (double)(1. / (double) Nbr_Voie));
    Step = log(D_Scale);

    for (j = 0; j < N; j++) Signal(j) = 0.;

    for (i = 0; i < Nbr_Plan; i++)
    {
        k = (int) Scale;
        k3 = (int) (3. * Scale);

        for (j = 0; j < N; j++)
        {
            tmp = 0.;
            for (l = j-k3; l < j-k; l++)
                tmp -=  W_1D(border_ind_test (l, N, Border), i);
            for (l =  j-k; l <=  j+k; l++)
                tmp += 2. * W_1D(border_ind_test (l, N, Border),i);
            for (l =  j+k+1; l <=  j+k3; l++)
                  tmp -= W_1D( border_ind_test (l, N, Border),i);

            Signal(j) += tmp / (Scale * C_psi_french) * Step;
        }
        Scale *= D_Scale;
    }
}
 
/***************************************************************************/


void wave_1d_morlet (fltarray & Signal, fltarray & W_1D, 
                     int N, type_border Border, 
                int Nbr_Voie, int Nbr_Plan, float Nu_0, float Scale0)
/* Morlet transform (Nu_0 > 0.8) */
{
    int i,j,k,l,jj;
    float D_Scale, x;
    float Norm, Omega, Val, Coef;
    float Scale = Scale0;

    D_Scale = pow ((double)2., (double)(1. / (double) Nbr_Voie));
//    printf ("Nbr_Plan = %d, Scale = %f\n", Nbr_Plan, Scale);
    
    Norm = 1. / sqrt(2.*PI);
    Omega = 2. * PI * Nu_0;

    for (i = 0; i < Nbr_Plan; i++)
    {
        k = (int) (6. * Scale);
        for (j = 0; j < N; j++)
        {
            W_1D(j,i) = 0.;
            W_1D(j, Nbr_Plan+i) = 0.;
            for (l =  j-k; l <  j+k; l++)
            {
                jj = border_ind_test (l, N, Border);
                x = (float)(j - l)/ Scale;
                Coef = Norm * exp(-x*x/2);
                Val = Omega * x;
                /* real part */
                W_1D(j,i) += Coef * cos(Val) * Signal(jj);
                /* imaginary part */
                W_1D(j,Nbr_Plan+i) -= Coef * sin(Val) * Signal(jj);
            }
            W_1D(j,i) /= Scale;
            W_1D(j,Nbr_Plan+i) /= Scale;
        }
        Scale *= D_Scale;
    }
}

/***************************************************************************/

void morlet_mod (fltarray &Wave, fltarray &Wave_Mod)
{
    int i,j;
    int N = Wave.axis(2);
    int Nbr_Plan= Wave.axis(1)/2;

    for (i = 0; i < Nbr_Plan; i++)
    for (j = 0; j < N; j++)
            Wave_Mod(j,i) = sqrt(Wave(j,i)*Wave(j,i) 
                                +Wave(j,Nbr_Plan+i)*Wave(j,Nbr_Plan+i));
}

/***************************************************************************/

void morlet_phase (fltarray &Wave, fltarray &Wave_Phase)
{
    int i,j;
    int N = Wave.axis(2);
    int Nbr_Plan= Wave.axis(1)/2;

    for (i = 0; i < Nbr_Plan; i++)
    {
        for (j = 0; j < N; j++)
        {
            ARG(Wave(j,i),Wave(j,Nbr_Plan+i), Wave_Phase(j,Nbr_Plan+i));
        }
    }
}

/***************************************************************************/

void morlet_re (fltarray &Wave_Mod, fltarray &Wave_Ph, fltarray &Wave)
{
    int i,j;
    int N = Wave_Mod.axis(2);
    int Nbr_Plan= Wave.axis(1)/2;

    for (i = 0; i < Nbr_Plan; i++)
    for (j = 0; j < N; j++)
            Wave(j,i) = Wave_Mod(j,i) * cos(Wave_Ph(j,i));
}

/***************************************************************************/

void morlet_im (fltarray &Wave_Mod, fltarray &Wave_Ph, fltarray &Wave)
{
    int i,j;
    int N = Wave_Mod.axis(2);
    int Nbr_Plan= Wave.axis(1)/2;

    for (i = 0; i < Nbr_Plan; i++)
    for (j = 0; j < N; j++)
            Wave(j,Nbr_Plan+i) = Wave_Mod(j,i) * sin(Wave_Ph(j,i));
}

/***************************************************************************/
// obsolete routine: keep it for compatibility
// see MR1D.transform
void mr1d_transform (fltarray & Signal, int Np, type_trans_1d Type_Transform, 
               type_border Border, int MedianWinSize,
               int Nbr_Voie,  fltarray & W_1D, int Nbr_Plan, 
               float Scale_0, float Nu_0, Bool Interp)
{
    int Iter,i,j;
    fltarray Resi (Np);
    fltarray MR_Iter (Np, Nbr_Plan);

    switch (Type_Transform)
    {
       case TO1_PAVE_LINEAR:
         wave_1d_linear (Signal, W_1D, Np, Nbr_Plan, Border); 
         break;
       case TO1_PAVE_B1SPLINE:
         wave_1d_spline1 (Signal, W_1D, Np, Nbr_Plan, Border);
         // mr1d_medianspline3_bis(Signal,W_1D,Nbr_Plan,MedianWinSize, Border);
         break;
       case TO1_PAVE_B3SPLINE:
         wave_1d_spline3 (Signal, W_1D, Np, Nbr_Plan,Border);
         break;
       case TO1_PAVE_B3_DERIV:
         wave_1d_B3deriv_atrou (Signal, W_1D, Nbr_Plan, Border);
         break;
       case TO1_PAVE_MORLET:
         wave_1d_morlet (Signal, W_1D, Np, Border,
                         Nbr_Voie, Nbr_Plan, Nu_0, Scale_0);
         break;
       case TO1_PAVE_MEX:
         wave_1d_mex (Signal, W_1D, Np, Border, Nbr_Voie, Nbr_Plan, Scale_0);
         break;
	case TO1_PAVE_DERIV_GAUSS:
	 wave_1d_der_gauss (Signal, W_1D, Np, Border,
	                    Nbr_Voie, Nbr_Plan, Scale_0);
	 break;
       case TO1_PAVE_FRENCH:
         wave_1d_french (Signal, W_1D, Np, Border,Nbr_Voie, Nbr_Plan, Scale_0);
         break;
       case TM1_PAVE_MEDIAN:
         mr1d_median (Signal, W_1D, Np, Nbr_Plan, MedianWinSize, Border);
         break;
       case TO1_PYR_LINEAR:
         pyr_1d_linear (Signal, W_1D, Np, Nbr_Plan, Border);
         break;
       case TO1_PYR_B3SPLINE:
         pyr_1d_spline3 (Signal, W_1D, Np, Nbr_Plan, Border);
         break;
       case TM1_PYR_MEDIAN:
         mr1d_pyr_median (Signal, W_1D, Np, Nbr_Plan, MedianWinSize, Border);
         break;
       default:
         cerr << "Error: not implemented here ..." << endl;
         exit(-1);
         break;
   }
   
   if ((Interp==True) && (which_set_is_trans1d(Type_Transform) == TRANS1_PYR))
   {
        for (Iter = 0; Iter < Nbr_Plan; Iter++)
        {  
            for (i = 0; i < Nbr_Plan; i++)
            for (j = 0; j < Np; j++)  MR_Iter(j,i) = 0.;

            for (i = 0; i < Np; i++)  MR_Iter(i,Iter) = W_1D(i,Iter);
                mr1d_pyr_rec (Resi, MR_Iter, Np, Nbr_Plan);
            for (i = 0; i < Np; i++)  W_1D(i,Iter) = Resi(i);
        }
  } 
} 

/*************************************************************************/
// obsolete routine: keep it for compatibility
// see MR1D.transform
void mr1d_recons (fltarray &W_1D, int Np, int Nbr_Plan, type_trans_1d T_Transf,
                 type_border Border, 
                 int Nbr_Voie, fltarray & Signal, float Scale_0)
{
    switch (T_Transf)
    {
       case TO1_PAVE_FRENCH: 
         wave_1d_french_rec (W_1D, Signal, Np, Border,
                             Nbr_Voie, Nbr_Plan, Scale_0); 
         break;
       case TO1_PAVE_MEX: 
         wave_1d_mex_rec (W_1D, Signal, Np, Border,
                             Nbr_Voie, Nbr_Plan, Scale_0);
         break;
       case TO1_PAVE_LINEAR:
       case TO1_PAVE_B1SPLINE :
       case TO1_PAVE_B3SPLINE :
       case TM1_PAVE_MEDIAN:
         wave_1d_algo_trou_rec (W_1D, Signal, Np, Nbr_Plan);
         break;
       case TO1_PAVE_B3_DERIV:
         wave_1d_rec_B3deriv_atrou (W_1D, Signal, Nbr_Plan, Border);
	 break;
       case TO1_PAVE_MORLET:
       case TO1_PAVE_DERIV_GAUSS:
         cerr << "Error: ";
         cerr << "This reconstruction is not implemented" << endl;
         exit (-1);
         break;
       case TO1_PYR_LINEAR:
       case TO1_PYR_B3SPLINE:
       case TM1_PYR_MEDIAN:
           mr1d_pyr_rec (Signal, W_1D, Np, Nbr_Plan);
         break;
       default:
         cerr << "Error: not implemented here ..." << endl;
         exit(-1);
         break;
   }
}

/*************************************************************************/

