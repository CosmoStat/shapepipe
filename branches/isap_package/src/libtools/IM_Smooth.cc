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
**    File:  IM_Smooth.cc
**
*******************************************************************************
**
**    DESCRIPTION  smoothing algorithms
**    -----------  
**                 
**
*******************************************************************************
**
**  void smooth_average (const Iint & Im_in, Iint &Im_out, 
**                     type_border Type, int Step_trou, int Size_Window)
**
**  smoothes the images by averaging a window of size Size_Window
**  Step_trou ==> the distance between two pixels is 2^Step_trou
**  by default, Step_trou is equal to zero
**
*******************************************************************************
**
**  void smooth_average (const Ifloat & Im_in, Ifloat &Im_out, 
**                     type_border Type, int Step_trou, int Size_Window)
**
**  smoothes the images by averaging a window
**  Step_trou ==> the distance between two pixels is 2^Step_trou
**  by default, Step_trou is equal to zero
**
*******************************************************************************
**
**  void smooth_linear (const Iint & Im_in, Iint &Im_out,
**                                   type_border Type, int Step_trou)
**  smoothes the images by the linear interpolation
**  Step_trou ==> the distance between two pixels is 2^Step_trou
**  by default, Step_trou is equal to zero
**
*******************************************************************************
**
**  void smooth_linear (const Ifloat & Im_in, Ifloat &Im_out,
**                                   type_border Type, int Step_trou)
**  smoothes the images by the linear interpolation
**  Step_trou ==> the distance between two pixels is 2^Step_trou
**  by default, Step_trou is equal to zero
**
*******************************************************************************
**
** void smooth_bspline (const Iint & Im_in, Iint &Im_out,
**                                   type_border Type, int Step_trou)
**  smoothes the images by the bspline interpolation
**  Step_trou ==> the distance between two pixels is 2^Step_trou
**  by default, Step_trou is equal to zero
**
*******************************************************************************
**
** void smooth_bspline (const Ifloat & Im_in, Ifloat &Im_out,
**                                   type_border Type, int Step_trou)
**  smoothes the images by the bspline interpolation
**  Step_trou ==> the distance between two pixels is 2^Step_trou
**  by default, Step_trou is equal to zero
**
*******************************************************************************
**
** void smooth_mediane (const Iint &Imag1, Iint &Imag2,  
**                      type_border Type, int Step_trou, int Window_Size)
**
**  smoothes the images by the mediane filtering
**  Step_trou ==> the distance between two pixels is 2^Step_trou
**  by default, Step_trou is equal to zero
**
*******************************************************************************
**
** void smooth_mediane (const Ifloat &Imag1, Ifloat &Imag2, 
**                     type_border Type, int Step_trou, int Window_Size)
**
**  smoothes the images by the mediane filtering
**  Step_trou ==> the distance between two pixels is 2^Step_trou
**  by default, Step_trou is equal to zero
**
******************************************************************************/

// static char sccsid[] = "@(#)IM_Smooth.cc 3.2 96/06/13 CEA 1994 @(#)";


#include "IM_Obj.h"
#include "OptMedian.h"

/****************************************************************************/

void smooth_average (const Iint & Im_in, Iint &Im_out, 
                     type_border Type, int Step_trou, int Size_Window)
{
    int Nl = Im_in.nl();
    int Nc = Im_in.nc();
    int i,j,k,l,Step;
    int Size_2 = Size_Window / 2;
    int Total = Size_Window*Size_Window;

    Step = (int)(pow((double) 2., (double) Step_trou) + 0.5);
    Size_2 *= Step;

    for (i = 0; i < Nl; i ++)
    for (j = 0; j < Nc; j ++)
    {
        Im_out(i,j) = 0;
        for (k = i - Size_2; k <= i+Size_2; k += Step)
        for (l = j - Size_2; l <= j+Size_2; l += Step)
                       Im_out(i,j) += Im_in(k,l,Type);
        Im_out(i,j) /= Total;
    }
}

/****************************************************************************/

void smooth_average (const Ifloat & Im_in, Ifloat &Im_out, 
                     type_border Type, int Step_trou, int Size_Window)
{
    int Nl = Im_in.nl();
    int Nc = Im_in.nc();
    int i,j,k,l,Step;
    int Size_2 = Size_Window / 2;
    float Total = Size_Window*Size_Window;

    Step = (int)(pow((double)2., (double) Step_trou) + 0.5);
    Size_2 *= Step;

    for (i = 0; i < Nl; i ++)
    for (j = 0; j < Nc; j ++)
    {
        Im_out(i,j) = 0.;
        for (k = i - Size_2; k <= i+Size_2; k += Step)
        for (l = j - Size_2; l <= j+Size_2; l += Step)
                       Im_out(i,j) += Im_in(k,l,Type);
        Im_out(i,j) /= Total;
    }
}

/*********************************************************************/

void smooth_haar (const Ifloat & Im_in, Ifloat &Im_out,
                                   type_border Type, int Step_trou, float Norm)
{
    int Nl = Im_in.nl();
    int Nc = Im_in.nc();
    int i,j,Step;
    Ifloat Buff(Nl,Nc,"Buff smooth_linear");

    Step = (int)(pow((double)2., (double) Step_trou) + 0.5);

    for (i = 0; i < Nl; i ++)
    for (j = 0; j < Nc; j ++)
       Buff(i,j) =    (Im_in(i,j) + Im_in (i, j-Step, Type));

    for (i = 0; i < Nl; i ++)
    for (j = 0; j < Nc; j ++)
       Im_out(i,j) =   (Buff(i,j) + Buff (i-Step, j, Type)) / Norm;
}

/****************************************************************************/

void smooth_linear (const Iint & Im_in, Iint &Im_out,
                                   type_border Type, int Step_trou)
{
    int Nl = Im_in.nl();
    int Nc = Im_in.nc();
    int i,j,Step;
    Ifloat Buff(Nl,Nc,"Buff smooth_linear");

    Step = (int)(pow((double)2., (double) Step_trou) + 0.5);

    for (i = 0; i < Nl; i ++)
    for (j = 0; j < Nc; j ++)
       Buff(i,j) = 0.5 *  (float)  Im_in(i,j) +
                   0.25 * (float) ( Im_in (i, j-Step, Type) 
                                  + Im_in (i, j+Step, Type)) + 0.5;

    for (i = 0; i < Nl; i ++)
    for (j = 0; j < Nc; j ++)
       Im_out(i,j) = (int)(0.5 * Buff(i,j) +
                    0.25 *  ( Buff (i-Step, j, Type) 
                            + Buff (i+Step, j, Type)) + 0.5);
}

/****************************************************************************/

void smooth_linear (const Ifloat & Im_in, Ifloat &Im_out,
                                   type_border Type, int Step_trou)
{
    int Nl = Im_in.nl();
    int Nc = Im_in.nc();
    int i,j,Step;
    Ifloat Buff(Nl,Nc,"Buff smooth_linear");

    Step = (int)(pow((double)2., (double) Step_trou) + 0.5);

    for (i = 0; i < Nl; i ++)
    for (j = 0; j < Nc; j ++)
       Buff(i,j) = 0.5 *  Im_in(i,j) +
                  0.25 * ( Im_in (i, j-Step, Type)  + Im_in (i, j+Step, Type));

    for (i = 0; i < Nl; i ++)
    for (j = 0; j < Nc; j ++)
       Im_out(i,j) = 0.5 * Buff(i,j) +
                    0.25 *  ( Buff (i-Step, j, Type)  + Buff (i+Step, j, Type));
}

/***************************************************************************/

void smooth_bspline (const Iint & Im_in, Iint &Im_out,
                                   type_border Type, int Step_trou)
{
    int Nl = Im_in.nl();
    int Nc = Im_in.nc();
    int i,j,Step;
    float Coeff_h0 = 3. / 8.;
    float Coeff_h1 = 1. / 4.;
    float Coeff_h2 = 1. / 16.;
    Ifloat Buff(Nl,Nc,"Buff smooth_bspline");

    Step = (int)(pow((double)2., (double) Step_trou) + 0.5);

    for (i = 0; i < Nl; i ++)
    for (j = 0; j < Nc; j ++)
       Buff(i,j) = Coeff_h0 * (float)  Im_in(i,j)
                 + Coeff_h1 * (float) (  Im_in (i, j-Step, Type) 
                                       + Im_in (i, j+Step, Type)) 
                 + Coeff_h2 * (float) (  Im_in (i, j-2*Step, Type) 
                                       + Im_in (i, j+2*Step, Type));

    for (i = 0; i < Nl; i ++)
    for (j = 0; j < Nc; j ++)
       Im_out(i,j) = (int)(Coeff_h0 * Buff(i,j) +
                 + Coeff_h1 * (float) (  Buff (i-Step, j, Type) 
                                       + Buff (i+Step, j, Type)) 
                 + Coeff_h2 * (float) (  Buff (i-2*Step, j, Type) 
                                       + Buff (i+2*Step, j, Type)) + 0.5);

}

/****************************************************************************/

void smooth_bspline (const Ifloat & Im_in, Ifloat &Im_out,
                                   type_border Type, int Step_trou)
{
    int Nl = Im_in.nl();
    int Nc = Im_in.nc();
    int i,j,Step;
    float Coeff_h0 = 3. / 8.;
    float Coeff_h1 = 1. / 4.;
    float Coeff_h2 = 1. / 16.;
    Ifloat Buff(Nl,Nc,"Buff smooth_bspline");

    Step = (int)(pow((double)2., (double) Step_trou) + 0.5);

    for (i = 0; i < Nl; i ++)
    for (j = 0; j < Nc; j ++)
       Buff(i,j) = Coeff_h0 * Im_in(i,j)
                 + Coeff_h1 * (  Im_in (i, j-Step, Type) 
                               + Im_in (i, j+Step, Type)) 
                 + Coeff_h2 * (  Im_in (i, j-2*Step, Type) 
                               + Im_in (i, j+2*Step, Type));

    for (i = 0; i < Nl; i ++)
    for (j = 0; j < Nc; j ++)
       Im_out(i,j) = Coeff_h0 * Buff(i,j)
                 + Coeff_h1 * (  Buff (i-Step, j, Type) 
                               + Buff (i+Step, j, Type)) 
                 + Coeff_h2 * (  Buff (i-2*Step, j, Type) 
                               + Buff (i+2*Step, j, Type));

}

/****************************************************************************/
/*
static void tri_i (int *fenetre, int Window_Size)
{
    int i,j, Dim_2;
    int Aux;

    Dim_2 = Window_Size * Window_Size;

    for (i = 1; i <= Dim_2/2+1; i++)
    {
        for (j = 0; j < Dim_2 - i; j++)
        {
            if (fenetre[j] > fenetre[j+1])
            {
                Aux = fenetre[j+1];
                fenetre[j+1] = fenetre[j];
                fenetre[j] = Aux;
            }
        }
    }
}
*/
/****************************************************************************/

void smooth_mediane (const Iint &Imag1, Iint &Imag2, 
                     type_border Type, int Step_trou, int Window_Size)
{
    int Nl = Imag1.nl();
    int Nc = Imag1.nc();
    int k,l,i,j,ind_fen;
    int Step;
    int *fenetre;
    int Window2 = (Window_Size - 1) / 2;
    // int Ind_Med = (Window_Size*Window_Size - 1) / 2;
    int Size = Window_Size*Window_Size;
    Step = (int)(pow((double)2., (double) Step_trou) + 0.5);

    fenetre = new int [Size];
    for (i = 0; i < Nl; i++) 
    {
        for (j = 0; j < Nc; j++) 
        {
            ind_fen = 0;
            for (k = i - Window2*Step; k <= i+Window2*Step; k += Step)
            for (l = j - Window2*Step; l <= j+Window2*Step; l += Step)
                                  fenetre[ind_fen++] = Imag1(k , l, Type);
            // tri_i (fenetre, Window_Size);
	    // Imag2(i,j) = fenetre[Ind_Med];
            if (Size == 9) Imag2(i,j) = opt_med9(fenetre);  
	    else Imag2(i,j) = get_median(fenetre, Size);		  
        }
    }
    delete[] fenetre;
}

/***************************************************************************/
 
void smooth_mediane (const Ifloat &Imag1, Ifloat &Imag2,  
                     type_border Type, int Step_trou, int Window_Size)
{
    int Nl = Imag1.nl();
    int Nc = Imag1.nc();
    int k,l,i,j,ind_fen;
    int Step;
    float *fenetre;
    int  Window2 = (Window_Size - 1) / 2;
    int  Size = Window_Size*Window_Size;

    if (Step_trou > 0) Step = (int)(pow((double)2., (double) Step_trou) + 0.5);
    else Step = 1;

    fenetre = new float [Size];
    fenetre[0]=0.;
    for (i = 0; i < Nl; i++) 
    for (j = 0; j < Nc; j++) 
    {
       ind_fen = 0;
       for (k = i - Window2*Step; k <= i+Window2*Step; k += Step)
       for (l = j - Window2*Step; l <= j+Window2*Step; l += Step)
                                  fenetre[ind_fen++] = Imag1(k , l, Type);
       // Imag2(i,j) = hmedian(fenetre, ind_fen); 
       if (Size == 9) Imag2(i,j) = opt_med9(fenetre);  
       else Imag2(i,j) = get_median(fenetre, Size);	       
    }
    delete[] fenetre;
}

/****************************************************************************/
