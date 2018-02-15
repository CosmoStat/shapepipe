/*******************************************************************************
**
**    UNIT
**
**    Version: 3.1
**
**    Author: Jean-Luc Starck
**
**    Date:  96/05/02 
**    
**    File:  IM_SubImag.cc
**
*******************************************************************************
**
**    DESCRIPTION  operation between images with different sizes
**    -----------  
**                 
**
*******************************************************************************
**
** void im_extract (const Ifloat &Imag, Ifloat &Imag_Out)
**
** extracts an image from Imag
** Imag_Out is the center of the image Imag
** we must have:
**            Imag.Nl  > Imag_Out.Nl
**            Imag.Nc  > Imag_Out.Nc
**
*****************************************************************************
**
** void im_extend (const Ifloat &Imag, Ifloat &Imag_Out)
**
** extends by miror effect an image Imag at the center of a  bigger one. 
** all pixels of Imag_Out are modified. A miror is used
** at the board
**
*             Imag.Nl < Imag_Out.Nl
**            Imag.Nc < Imag_Out.Nc
**
**
*****************************************************************************
**
** void im_reduce_size_2 (const Ifloat &Imag, Ifloat &Imag_Out)
**
** Imag_Out = Imag by taken only on pixel over two
** Imag_Out.Nl = Imag.Nl / 2
** Imag_Out.Nc = Imag.Nc / 2
**
** 
*****************************************************************************
**
** void im_reduce_size_2 (const Iint &Imag, Iint &Imag_Out)
**
** Imag_Out = Imag by taken only on pixel over two
** Imag_Out.Nl = Imag.Nl / 2
** Imag_Out.Nc = Imag.Nc / 2
**
*****************************************************************************
**
** void im_increase_size_2 (const Ifloat &Pict_in, Ifloat &Pict_out, 
**                          type_border Border)
**
** interpolates Pict_in by a factor 2 by a Bspline interpolation
**
*****************************************************************************
**
** void im_increase_size_2 (const Iint &Pict_in, Iint &Pict_out
**                          type_border Border)
**
** interpolates Pict_in by a factor 2
** 
*****************************************************************************
**
** void im_bilinear_interp (Ifloat &INimage, Ifloat &OUTimage)
** 
** bilinear interpolation: the size of OUTimage must be superior than
** the size of INimage
**
*****************************************************************************
** 
** void im_block_extend(const Ifloat &Imag, Ifloat &Imag_Out)
**
** block interpolation: the size of OUTimage must be superior than
** the size of INimage
**
**
*****************************************************************************
** 
** void im_cf_interp (Ifloat &Image, int Nls, int Ncs)
**
** image extension by Shannon's interpolation
** the FFT of Image is computed. The Fourier domain is extended by
** zero around the data, then the inverse transform is computed
** the output Image has the size Nls, Ncs
**
******************************************************************************/

// static char sccsid[] = "@(#)IM_SubImag.cc 3.1 96/05/02 CEA 1994 @(#)";

#include "IM_Obj.h"

#define SIZE_TAB_COEF_INTER 10  /* Taille de la table des coefficients */

#define C1 .00245
#define C2 -.00915
#define C3 .03414
#define C4 -0.12720
#define C5 0.60048

#define CS ((C1+C2+C3+C4+C5)*2.)
#define C1n (C1/CS)
#define C2n (C2/CS)
#define C3n (C3/CS)
#define C4n (C4/CS)
#define C5n (C5/CS)

static float Tab_Coef_Inter[SIZE_TAB_COEF_INTER] =
                              {C1n,C2n,C3n,C4n,C5n,C5n,C4n,C3n,C2n,C1n};

/****************************************************************************/

void im_block_extend(const Ifloat &Imag, Ifloat &Imag_Out)
{
    int Nl0 = Imag.nl();
    int Nc0 = Imag.nc();
    int Nl1 = Imag_Out.nl();
    int Nc1 = Imag_Out.nc();
    int i,j;
    float Zoom_i = (float) Nl1 / (float) Nl0;
    float Zoom_j = (float) Nc1 / (float) Nc0;

    for (i = 0; i < Nl1; i++)
    for (j = 0; j < Nc1; j++)
    {
       Imag_Out (i,j) = Imag ((int)(i/Zoom_i),(int)(j/Zoom_j));
    }
}

/****************************************************************************/

void im_extract (const Ifloat &Imag, Ifloat &Imag_Out)
{
    int Nl0 = Imag_Out.nl();
    int Nc0 = Imag_Out.nc();
    int Nl1 = Imag.nl();
    int Nc1 = Imag.nc();
    int i1,j1,i0,j0, Depi, Depj;

    Depi = (Nl1 - Nl0) / 2;
    Depj = (Nc1 - Nc0) / 2;

    for (i0 = 0; i0 < Nl0; i0++)
    for (j0 = 0; j0 < Nc0; j0++)
    {
       i1 = i0 + Depi;
       j1 = j0 + Depj;
       Imag_Out (i0,j0) = Imag (i1,j1);
    }
}

/****************************************************************************/

void im_extract (const Icomplex_f &Imag, Icomplex_f &Imag_Out)
{
    int Nl0 = Imag_Out.nl();
    int Nc0 = Imag_Out.nc();
    int Nl1 = Imag.nl();
    int Nc1 = Imag.nc();
    int i1,j1,i0,j0, Depi, Depj;

    Depi = (Nl1 - Nl0) / 2;
    Depj = (Nc1 - Nc0) / 2;

    for (i0 = 0; i0 < Nl0; i0++)
    for (j0 = 0; j0 < Nc0; j0++)
    {
       i1 = i0 + Depi;
       j1 = j0 + Depj;
       Imag_Out (i0,j0) = Imag (i1,j1);
    }
}

/**************************************************************************/

void im_extract (const Icomplex_f &Imag, Ifloat &Imag_Out, float Norm)
{
    int Nl0 = Imag_Out.nl();
    int Nc0 = Imag_Out.nc();
    int Nl1 = Imag.nl();
    int Nc1 = Imag.nc();
    int i1,j1,i0,j0, Depi, Depj;

    Depi = (Nl1 - Nl0) / 2;
    Depj = (Nc1 - Nc0) / 2;

    for (i0 = 0; i0 < Nl0; i0++)
    for (j0 = 0; j0 < Nc0; j0++)
    {
       i1 = i0 + Depi;
       j1 = j0 + Depj;
       Imag_Out (i0,j0) = Imag (i1,j1).real() / Norm;
    }
}

/**************************************************************************/
static int test_ind (int ind, int N)
{
    int Val;

    if (ind < 0) Val = - ind;
    else
    {
        if (ind >= N) Val = 2 * (N - 1) - ind;
        else Val = ind;
    }
    if ((Val >= N) || (Val < 0)) Val = -1;
    return (Val);
}

/**************************************************************************/

void im_zero_padding (const Ifloat &Imag, Ifloat &Imag_Out)
{
    int Nl0 = Imag.nl();
    int Nc0 = Imag.nc();
    int Nl1 = Imag_Out.nl();
    int Nc1 = Imag_Out.nc();
    int i1,j1,i0,j0, Depi, Depj;

    Depi = (Nl1 - Nl0) / 2;
    Depj = (Nc1 - Nc0) / 2;

    for (i1 = 0; i1 < Nl1; i1++)
    for (j1 = 0; j1 < Nc1; j1++)
    {
       i0 = i1 - Depi;
       j0 = j1 - Depj;
       if ((i0 < 0) || (j0 < 0) || (i0 >= Nl0) || (j0  >= Nc0)) 
            Imag_Out (i1,j1) = 0.;
       else Imag_Out (i1,j1) = Imag(i0,j0);
    }
}

/**************************************************************************/

void im_zero_padding (const Ifloat &Imag,  Icomplex_f &Imag_Out)
{
    int Nl0 = Imag.nl();
    int Nc0 = Imag.nc();
    int Nl1 = Imag_Out.nl();
    int Nc1 = Imag_Out.nc();
    int i1,j1,i0,j0, Depi, Depj;
    complex_f Zero = complex_f(0.,0.);

    Depi = (Nl1 - Nl0) / 2;
    Depj = (Nc1 - Nc0) / 2;

    for (i1 = 0; i1 < Nl1; i1++)
    for (j1 = 0; j1 < Nc1; j1++)
    {
       i0 = i1 - Depi;
       j0 = j1 - Depj;
       if ((i0 < 0) || (j0 < 0) || (i0 >= Nl0) || (j0  >= Nc0)) 
            Imag_Out(i1,j1) = Zero;
       else Imag_Out (i1,j1) = complex_f(Imag(i0,j0), 0.);
    }
}

/**************************************************************************/

void im_extend (const Ifloat &Imag, Ifloat &Imag_Out)
{
    int Nl0 = Imag.nl();
    int Nc0 = Imag.nc();
    int Nl1 = Imag_Out.nl();
    int Nc1 = Imag_Out.nc();
    int i1,j1,i0,j0, Depi, Depj;

    Depi = (Nl1 - Nl0) / 2;
    Depj = (Nc1 - Nc0) / 2;

    for (i1 = 0; i1 < Nl1; i1++)
    for (j1 = 0; j1 < Nc1; j1++)
    {
       i0 = test_ind (i1 - Depi, Nl0);
       j0 = test_ind (j1 - Depj, Nc0);
       if ((i0 < 0) || (j0 < 0)) Imag_Out (i1,j1) = 0.;
       else Imag_Out (i1,j1) = Imag (i0,j0);
    }
}

/**************************************************************************/

void im_extend (const Ifloat &Imag, Icomplex_f &Imag_Out)
{
    int Nl0 = Imag.nl();
    int Nc0 = Imag.nc();
    int Nl1 = Imag_Out.nl();
    int Nc1 = Imag_Out.nc();
    int i1,j1,i0,j0, Depi, Depj;
    complex_f Zero = complex_f(0.,0.);
    
    Depi = (Nl1 - Nl0) / 2;
    Depj = (Nc1 - Nc0) / 2;

    for (i1 = 0; i1 < Nl1; i1++)
    for (j1 = 0; j1 < Nc1; j1++)
    {
       i0 = test_ind (i1 - Depi, Nl0);
       j0 = test_ind (j1 - Depj, Nc0);
       if ((i0 < 0) || (j0 < 0)) Imag_Out(i1,j1) = Zero;
       else Imag_Out (i1,j1) = complex_f(Imag(i0,j0), 0.);
    }
}

/**************************************************************************/

void im_extend (const Icomplex_f &Imag, Icomplex_f &Imag_Out)
{
    int Nl0 = Imag.nl();
    int Nc0 = Imag.nc();
    int Nl1 = Imag_Out.nl();
    int Nc1 = Imag_Out.nc();
    int i1,j1,i0,j0, Depi, Depj;

    Depi = (Nl1 - Nl0) / 2;
    Depj = (Nc1 - Nc0) / 2;

    for (i1 = 0; i1 < Nl1; i1++)
    for (j1 = 0; j1 < Nc1; j1++)
    {
       i0 = test_ind (i1 - Depi, Nl0);
       j0 = test_ind (j1 - Depj, Nc0);
       if ((i0 < 0) || (j0 < 0)) Imag_Out (i1,j1) = complex_f(0., 0.);
       else Imag_Out (i1,j1) = Imag (i0,j0);
    }
}

/**************************************************************************/

void im_extend (const Iint &Imag, Iint &Imag_Out)
{
    int Nl0 = Imag.nl();
    int Nc0 = Imag.nc();
    int Nl1 = Imag_Out.nl();
    int Nc1 = Imag_Out.nc();
    int i1,j1,i0,j0, Depi, Depj;

    Depi = (Nl1 - Nl0) / 2;
    Depj = (Nc1 - Nc0) / 2;

    for (i1 = 0; i1 < Nl1; i1++)
    for (j1 = 0; j1 < Nc1; j1++)
    {
       i0 = test_ind (i1 - Depi, Nl0);
       j0 = test_ind (j1 - Depj, Nc0);
       if ((i0 < 0) || (j0 < 0)) Imag_Out (i1,j1) = 0;
       else Imag_Out (i1,j1) = Imag (i0,j0);
    }
}

/**************************************************************************/

void im_reduce_size_2 (const Ifloat &Imag, Ifloat &Imag_Out)
{
    int Nl = Imag_Out.nl();
    int Nc = Imag_Out.nc();
    int i,j;

    for (i = 0; i < Nl; i++)
    for (j = 0; j < Nc; j++)
         Imag_Out(i,j) = Imag(2*i,2*j, I_CONT);
}

/**************************************************************************/

void im_min_reduce_size_2 (const Ifloat &Imag, Ifloat &Imag_Out)
{
    int Nl = Imag.nl();
    int Nc = Imag.nc();
    int i,j;
    float Min;

    for (i = 0; i < Nl; i+=2)
    for (j = 0; j < Nc; j+=2)
    {
        Min = Imag(i,j);
        if (Min > Imag(i-1, j, I_CONT)) Min = Imag(i-1, j, I_CONT);
        if (Min > Imag(i+1, j, I_CONT)) Min = Imag(i+1, j, I_CONT);
        if (Min > Imag(i, j-1, I_CONT)) Min = Imag(i, j-1, I_CONT);
        if (Min > Imag(i, j+1, I_CONT)) Min = Imag(i, j+1, I_CONT);
/*
        if (Min > Imag(i-1, j-1, I_CONT)) Min = Imag(i-1, j-1, I_CONT);
        if (Min > Imag(i+1, j-1, I_CONT)) Min = Imag(i+1, j-1, I_CONT);
        if (Min > Imag(i-1, j-1, I_CONT)) Min = Imag(i-1, j-1, I_CONT);
        if (Min > Imag(i-1, j+1, I_CONT)) Min = Imag(i-1, j+1, I_CONT);
*/
        Imag_Out(i/2,j/2) = Min;
    }
}

/**************************************************************************/

void im_reduce_size_2 (const Iint &Imag, Iint &Imag_Out)
{
    int Nl = Imag_Out.nl();
    int Nc = Imag_Out.nc();
    int i,j;

    for (i = 0; i < Nl; i++)
    for (j = 0; j < Nc; j++)
         Imag_Out(i,j) = Imag(2*i,2*j, I_CONT);
}

/**************************************************************************/

void im_increase_size_2 (const Ifloat &Pict_in, Ifloat &Pict_out, type_border Border)
{
    int i,j,k,indj;
    int Nl2 = Pict_in.nl();
    int Nc2 = Pict_in.nc();
    int Nl = Pict_out.nl();
    int Nc = Pict_out.nc();
    int pi,pj,lastnc;

    /* Pixels at node interpolation */
    for (i = 0; i < Nl2; i ++)
    {
        pi = 2*i;

        if (pi < Nl)
        {
            for (j = 0; j < Nc2; j ++)
            {
                pj = 2*j;
                if (pj < Nc) Pict_out(pi,pj)  = Pict_in (i,j);
            }
        }
    }

    /* pixel node in column */
    for (i = 0; i < Nl2; i ++)
    {
        pi = 2*i;
        if (pi < Nl)
        {
            for (j = 0; j < Nc2; j ++)
            {
                pj = 2*j+1;
                if (pj < Nc)
                {
                    Pict_out(pi,pj) = 0.;
                    for (k = 0; k < SIZE_TAB_COEF_INTER; k ++)
                    Pict_out(pi,pj) += 
                        Tab_Coef_Inter [k] * Pict_in(i,j+k-4, Border);
                }
            }
        }
    }

    /* pixel node in line */
    for (i = 0; i < Nl2; i ++)
    {
        pi = 2*i+1;
        if (pi < Nl)
        {
            for (j = 0; j < Nc2; j ++)
            {
                pj = 2*j;
                if (pj < Nc)
                {
		    Pict_out(pi,pj) = 0.;
                    for (k = 0; k < SIZE_TAB_COEF_INTER; k ++)
                    Pict_out(pi,pj) += 
                         Tab_Coef_Inter [k] * Pict_in(i+k-4,j, Border);
                }
            }
        }
    }

    /* pixel node in line and column */
    for (i = 0; i < Nl2; i ++)
    {
        pi = 2*i+1;
        if (pi < Nl)
        {
            if (Nc % 2 == 0) lastnc = Nc - 2;
            else lastnc = Nc - 1;
        
            for (j = 0; j < Nc2; j ++)
            {
                pj = 2*j+1;
                if (pj < Nc)
                {
                    Pict_out(pi,pj) = 0.;
                    for (k = 0; k < SIZE_TAB_COEF_INTER; k ++)
                    {
                       indj = 2*j+2*(k-4);
                       if ((Border == I_CONT) && (indj >= Nc)) indj = lastnc;
                       Pict_out(pi,pj) += Tab_Coef_Inter [k] *
                                               Pict_out(pi,indj, Border);
                    }
                }
            }
        }
    }
}

/***************************************************************************/

void im_increase_size_2 (const Iint &Pict_in, Iint &Pict_out,type_border Border)
{
    int i,j,k,indj;
    int Nl2 = Pict_in.nl();
    int Nc2 = Pict_in.nc();
    int Nl = Pict_out.nl();
    int Nc = Pict_out.nc();
    Ifloat Buff(Nl, Nc, "Buff interpol");
    int pi,pj,lastnc=0;

    /* Pixels at node interpolation */
    for (i = 0; i < Nl2; i ++)
    {
        pi = 2*i;

        if (pi < Nl)
        {
            for (j = 0; j < Nc2; j ++)
            {
                pj = 2*j;
                if (pj < Nc) Buff(pi,pj)  = Pict_in (i,j);
            }
        }
    }

    /* pixel node in column */
    for (i = 0; i < Nl2; i ++)
    {
        pi = 2*i;

        if (pi < Nl)
        {
            for (j = 0; j < Nc2; j ++)
            {
                pj = 2*j+1;
                if (pj < Nc)
                {
                    Buff(pi,pj) = 0.;
                    for (k = 0; k < SIZE_TAB_COEF_INTER; k ++)
                    Buff(pi,pj) += 
                        Tab_Coef_Inter [k] * Pict_in(i,j+k-4, Border);
                }
            }
        }
    }

    /* pixel node in line */
    for (i = 0; i < Nl2; i ++)
    {
        pi = 2*i+1;
        if (pi < Nl)
        {
            for (j = 0; j < Nc2; j ++)
            {
                pj = 2*j;
                if (pj < Nc)
                {
		    Buff(pi,pj) = 0.;
                    for (k = 0; k < SIZE_TAB_COEF_INTER; k ++)
                    Buff(pi,pj) += 
                         Tab_Coef_Inter [k] * Pict_in(i+k-4,j, Border);
                }
            }
        
            if (Nc % 2 == 0) lastnc = Nc - 2;
            else lastnc = Nc - 1;
        }
    }


    /* pixel node in line and column */
    for (i = 0; i < Nl2; i ++)
    {
        pi = 2*i+1;
        if (pi < Nl)
        {
            for (j = 0; j < Nc2; j ++)
            {
                pj = 2*j+1;
                if (pj < Nc)
                {
                    Buff(pi,pj) = 0.;
                    for (k = 0; k < SIZE_TAB_COEF_INTER; k ++)
                    {
                       indj = 2*j+2*(k-4);
                       if ((Border == I_CONT) && (indj >= Nc)) indj = lastnc;
                       Buff(pi,pj) += Tab_Coef_Inter [k] *
                                               Buff(pi,indj, Border);
                    }
                }
            }
        }
    }
    for (i = 0; i < Nl; i ++)
    for (j = 0; j < Nc; j ++) Pict_out(i,j) = (int) Buff(i,j);
}

/***************************************************************************/

void im_bilinear_interp (Ifloat &INimage, Ifloat &OUTimage)
{		
    int INlin = INimage.nl();
    int INcol = INimage.nc();
    int OUTlin = OUTimage.nl();
    int OUTcol = OUTimage.nc();
    int i, j, ii, jj, i1, j1, i2, j2;
    int *indice_lin, *indice_col;
    float delta, delta_x, delta_y, delta_yx;

    indice_lin = new int [INlin];
    indice_col = new int [INcol];
    for (i = 0; i < INlin; i++)  indice_lin [i] = (i * (OUTlin-1)) / (INlin-1);
    for (i = 0; i < INcol; i++)  indice_col [i] = (i  *(OUTcol-1)) / (INcol-1);

    for (i=0; i<INlin; i++)
    for (j=0; j<INcol; j++)
            OUTimage (indice_lin[i],indice_col [j]) = INimage (i,j);

    for (i=0; i<INlin; i++)
    {
        ii = indice_lin[i];
        for (j=0; j<(INcol-1); j++)
        {
            j1 = indice_col [j];
            j2 = indice_col [j+1];
            delta_x = (float) j2-j1;
            delta_y = OUTimage (ii, j2) - OUTimage (ii, j1);
            delta_yx = delta_y/delta_x;
            for (jj=j1; jj<j2; jj++)
            {
                delta = (float) jj-j1;
                OUTimage (ii, jj) = OUTimage (ii, j1) + (delta * delta_yx);
            }
        }
     }

     for (jj = 0; jj < OUTcol; jj++)
     {
         for (i = 0; i < (INlin-1); i++)
         {
             i1 = indice_lin[i];
             i2 = indice_lin[i+1];
             delta_x = (float) i2-i1;
             delta_y = OUTimage (i2, jj) - OUTimage (i1, jj);
             delta_yx = delta_y / delta_x;
             for (ii = i1; ii < i2; ii++)
             {
                 delta = (float) ii-i1;
                 OUTimage (ii, jj) = OUTimage (i1, jj) + (delta * delta_yx);
             }
         }
    }
    delete [] indice_lin;
    delete [] indice_col;
}

/*****************************************************************************/
/*
void im_cf_interp (Ifloat &Image, int Nls, int Ncs)
{
    int Nl = Image.nl();
    int Nc = Image.nc();
    Icomplex_f I_cf (Nls, Ncs, "I cf interp ");
    Icomplex_f Buff (Nl, Nc, "Buff cf interp ");
    int i1,j1,i0,j0, Depi, Depj;
    float Coef = (float)(Nls*Ncs) / (float)(Nl*Nc);

#if DEBUG
printf ("Interpolation: %d,%d -> %d,%d\n", Nl, Nc,Nls, Ncs);
printf ("Min = %f\n", min(Image));
printf ("Max = %f\n", max(Image));
#endif
    Depi = (Nls - Nl) / 2;
    Depj = (Ncs - Nc) / 2;

    fft2d (Image, Buff, 1);

    for (i1 = 0; i1 < Nl; i1++)
    for (j1 = 0; j1 < Nc; j1++) 
    {
       i0 = i1 + Depi;
       j0 = j1 + Depj;
       I_cf (i0,j0) = Buff(i1,j1) * Coef;
    }

    Image.resize(Nls, Ncs);
    Image = real (invfft (I_cf));
#if DEBUG
printf ("Min Interp = %f\n", min(Image));
printf ("Max Interp = %f\n", max(Image));
#endif
}
*/
/****************************************************************************/

