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
**    File:  MR_Feauveau.cc
**
*******************************************************************************
**    DESCRIPTION  routines used for the G-transform algorithm 
**    -----------  
**
*******************************************************************************
**
** void tr_feauveau_one_scale (Ifloat &Image, Ifloat &Im_Out, int Nl, int Nc)
**
** computes the wavelet transform with feauveau's algorithm.
** only two scales are computed. 
** Nl, Nc is the size of the subimage contained in Image from which 
** we compute the wavelet transform
*******************************************************************************
**
** void transform_feauveau (Ifloat &Image, Ifloat &Im_Out, int Nbr_Plan)
**
** computes the wavelet transform with feauveau's algorithm.
** A scale contains the wavelet coefficient of this scale and of the half
** resolution and a smooth image
** the number of resolutions is equal to 2*Nbr_Plan+1
**
*******************************************************************************
**
** void rec_feauveau_one_scale (Ifloat &Image, Ifloat &Im_Out, int Nl, int Nc)
** 
** reconstruction of one scale
**
*******************************************************************************
**
** void recons_feauveau (Ifloat &Image, Ifloat &Im_Out, int Nbr_Iter)
**
** Reconstructs an image from its wavelet transform
** Image contains the wavelet transform by feauveau's algorithm.
** and Nbr_Iter = number of scale -1 
**
***************************************************************************/ 

// static char sccsid[] = "@(#)MR_Feauveau.cc 3.1 96/05/02 CEA 1994 @(#)";


#include <stdio.h>
#include <math.h>
#include <string.h>

#include "IM_Math.h"
#include "IM_Obj.h"
#include "MR_Obj.h"

/***************************************************************************/

#define DIM_FILTRE           7
#define DIM_FILTRE_DEMI      3 

#define DIM_FILTRE_DIAG     13
#define DIM_FILTRE_DIAG_DEMI 6

/* Valeurs de coefficients des filtres. Les lettres les representant
 * correspondent a celles donnees par Feauveau dans son article.
 */
#define DIM_FILTER 7
#define DIM_FILTER_DIAG 13

static float A = 0.6433908;
static float B = 0.1433488;
static float C = -0.0289743;
static float D = -0.0125361;
static float E = -0.0091506;
static float F = 0.0027737;
static float G = -0.0021526;
static float H = 0.0016453;
static float L = 0.0010524;
static float M = -0.0004016;
static float O = 0.0;

static float Filtre_U [DIM_FILTER*DIM_FILTER] = {
             M, L, H, G, H, L, M,
             L, F, E, D, E, F, L,
             H, E, C, B, C, E, H,
             G, D, B, A, B, D, G,
             H, E, C, B, C, E, H,
             L, F, E, D, E, F, L,
             M, L, H, G, H, L, M};

static float Filtre_V[DIM_FILTER*DIM_FILTER] = {
              M, -L,  H, -G,  H, -L,  M,
             -L,  F, -E,  D, -E,  F, -L,
              H, -E,  C, -B,  C, -E,  H,
             -G,  D, -B,  A, -B,  D, -G,
              H, -E,  C, -B,  C, -E,  H,
             -L,  F, -E,  D, -E,  F, -L,
              M, -L,  H, -G,  H, -L,  M};

static float Filtre_Diag_U[DIM_FILTER_DIAG*DIM_FILTER_DIAG] = {
  O, O, O, O, O, O, M, O, O, O, O, O, O,
  O, O, O, O, O, L, O, L, O, O, O, O, O,
  O, O, O, O, H, O, F, O, H, O, O, O, O,
  O, O, O, G, O, E, O, E, O, G, O, O, O,
  O, O, H, O, D, O, C, O, D, O, H, O, O,
  O, L, O, E, O, B, O, B, O, E, O, L, O,
  M, O, F, O, C, O, A, O, C, O, F, O, M,
  O, L, O, E, O, B, O, B, O, E, O, L, O,
  O, O, H, O, D, O, C, O, D, O, H, O, O,
  O, O, O, G, O, E, O, E, O, G, O, O, O,
  O, O, O, O, H, O, F, O, H, O, O, O, O,
  O, O, O, O, O, L, O, L, O, O, O, O, O,
  O, O, O, O, O, O, M, O, O, O, O, O, O};

static float Filtre_Diag_V[DIM_FILTER_DIAG*DIM_FILTER_DIAG] = {
  O, O, O, O, O, O, M, O, O, O, O, O, O,
  O, O, O, O, O,-L, O,-L, O, O, O, O, O,
  O, O, O, O, H, O, F, O, H, O, O, O, O,
  O, O, O,-G, O,-E, O,-E, O,-G, O, O, O,
  O, O, H, O, D, O, C, O, D, O, H, O, O,
  O,-L, O,-E, O,-B, O,-B, O,-E, O,-L, O,
  M, O, F, O, C, O, A, O, C, O, F, O, M,
  O,-L, O,-E, O,-B, O,-B, O,-E, O,-L, O,
  O, O, H, O, D, O, C, O, D, O, H, O, O,
  O, O, O,-G, O,-E, O,-E, O,-G, O, O, O,
  O, O, O, O, H, O, F, O, H, O, O, O, O,
  O, O, O, O, O,-L, O,-L, O, O, O, O, O,
  O, O, O, O, O, O, M, O, O, O, O, O, O};

/*************************************************************************/

inline float convolve_filter (Ifloat &Image, int Indi, int Indj, 
                        float *Filter, int Size_Filter, int Nl, int Nc)
{ 
    int k,l;
    float Val_Return = 0.;
    int Half_Filter = Size_Filter/2;
    float *Ptr_Filter = Filter;

    for (k = Indi-Half_Filter; k <= Indi+Half_Filter; k++)
    for (l = Indj-Half_Filter; l <= Indj+Half_Filter; l++)
    {
        if (*Ptr_Filter != 0.)
            Val_Return += Image( test_index_mirror(k,Nl), 
                                 test_index_mirror(l,Nc)) * (*Ptr_Filter);
        Ptr_Filter ++;
    } 
    return Val_Return;
}


/*************************************************************************/

void tr_feauveau_one_scale (Ifloat &Image, Ifloat &Im_Out, int Nl, int Nc)
{
    int i,j;
    char Buff[100];
    sprintf (Buff, "tr feauveau one scale: %d, %d", Nl, Nc);
    Ifloat Im_Aux(Nl, Nc, Buff);

    /* half resolution transform */
     
    for (i = 0; i < Nl; i++)
    for (j = 0; j < Nc; j++)
    {
        if (   ((i%2 == 0) && (j%2 == 0)) 
            || ((i%2 != 0) && (j%2 != 0)))
        {
            Im_Aux(i,j) = convolve_filter(Image,i,j,Filtre_U,DIM_FILTER,Nl,Nc);
         }
        else
        {
           if ((Nc %2 == 0) || (j < Nc-1))
           Im_Out(i, (Nc+1)/2 + j/2) = 
                 convolve_filter(Image, i,j, Filtre_V, DIM_FILTER,Nl,Nc);
        }
     }
       
       
    /* next resolution */
    for (i = 0; i < Nl; i++)
    for (j = 0; j < Nc; j++)
    {
        if ((i%2) && (j%2)) 
        {
            Im_Out(i/2,j/2) =  
             convolve_filter (Im_Aux, i,j,Filtre_Diag_U, DIM_FILTER_DIAG,Nl,Nc);
        }
        if (!(i%2) && !(j%2))
        {
           if ((Nl%2 == 0) || (i < Nl-1))
            Im_Out((Nl+1)/2 + i/2, j/2) =   
               convolve_filter(Im_Aux,i,j,Filtre_Diag_V, DIM_FILTER_DIAG,Nl,Nc);
        }
    }
}


/*************************************************************************/

void transform_feauveau (Ifloat &Image, Ifloat &Im_Out, int Nbr_Plan)
{
    int s;
    int Nl = Image.nl();
    int Nc = Image.nc();
    Ifloat Im_Feauv(Nl, Nc, "tr Feauveau");

    Im_Out = Image;
    for (s = 0; s < Nbr_Plan; s++)
    {
       Im_Feauv = Im_Out;
       tr_feauveau_one_scale (Im_Feauv, Im_Out, Nl , Nc);
/*
Max = 0.;
for (i = 0; i < Nl; i++)
for (j = 0; j < Nc; j++)
    if ((j >= Nc / 2) || (i >= Nl/2))
                  Max =  (Max > Im_Out(i,j)) ? Max:Im_Out(i,j);
for (i = 0; i < Nl; i++)
for (j = 0; j < Nc; j++)
                  Im_Out(i,j) /= Max;

*/
       Nl = (Nl+1)/2;
       Nc = (Nc+1)/2;
    }
/*
Max = 0.;
for (i = 0; i < Nl; i++)
for (j = 0; j < Nc; j++)
                Max =  (Max > Im_Out(i,j)) ? Max:Im_Out(i,j);
for (i = 0; i < Nl; i++)
for (j = 0; j < Nc; j++)
                  Im_Out(i,j) /= Max;
*/

}

/*************************************************************************/

void rec_feauveau_one_scale (Ifloat &Image, Ifloat &Im_Out, int Nl, int Nc)
{ 
    int i,j;
    char Buff[100];
    sprintf (Buff, "Rec feauveau one scale: %d, %d", Nl, Nc);
    Ifloat Im_Res(Nl, Nc, Buff);
    Ifloat Im_Resi(Nl, Nc, "Buffer Smooth");
    Ifloat Im_Det(Nl, Nc, "Buffer Details");

    for (i = 0; i < Nl; i++)
    for (j = 0; j < Nc; j++)
    {
        if ((i%2) && (j%2)) Im_Resi(i,j) = Image(i/2,j/2);
        else Im_Resi(i,j) = 0.;
        if (!(i%2) && !(j%2)) 
        {
         if ((Nl%2 == 0) || (i < Nl-1)) 
                Im_Det (i,j) = Image ((Nl+1)/2+i/2, j/2);
        }
        else Im_Det(i,j) = 0.;
    }

    for (i = 0; i < Nl; i++)
    for (j = 0; j < Nc; j++)
    {
        if (   ((i%2) && (j%2)) 
            || (!(i%2) && !(j%2)))
        {
            Im_Res(i,j) = 2.*(
                          convolve_filter (Im_Resi, i,j, Filtre_Diag_U,
                                           DIM_FILTER_DIAG, Nl, Nc)
                         + convolve_filter (Im_Det, i,j, Filtre_Diag_V,
                                           DIM_FILTER_DIAG, Nl, Nc));
        }
        else Im_Res(i,j) = 0.;
    }

    for (i = 0; i < Nl; i++)
    for (j = 0; j < Nc; j++)
    {
        if (   ((i%2) && !(j%2)) 
            || (!(i%2) && (j%2))) 
        {
            if ((Nc%2 == 0) || (j < Nc-1)) 
                     Im_Det(i,j) = Image (i, (Nc+1)/2+j/2);
        }
        else Im_Det(i,j) = 0.;
    }
    for (i = 0; i < Nl; i++)
    for (j = 0; j < Nc; j++)
    {
        Im_Out(i,j) = 2.*(
                          convolve_filter (Im_Res, i,j, Filtre_U,
                                           DIM_FILTER, Nl, Nc)
                         + convolve_filter (Im_Det, i,j, Filtre_V,
                                           DIM_FILTER, Nl, Nc));
    }
}             


/*************************************************************************/

void recons_feauveau (Ifloat &Image, Ifloat &Im_Out, int Nbr_Iter)
{
    int s;
    int Nl = Image.nl();
    int Nc = Image.nc();
    int Nls,Ncs;
    Ifloat Im_Feauv(Nl, Nc, "Rec Feauv");

    Im_Out = Image;
    for (s = Nbr_Iter-1; s >= 0; s--)
    {
       Nls = size_ima_resol(s, Nl);
       Ncs = size_ima_resol(s, Nc);
  
       Im_Feauv = Im_Out;
       rec_feauveau_one_scale (Im_Feauv, Im_Out, Nls , Ncs);
    }
}

/*************************************************************************/







