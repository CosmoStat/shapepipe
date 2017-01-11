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
**    File:  MR_Pav_Feau.cc
**
*******************************************************************************
**    DESCRIPTION  routines used for the G-transform algorithm 
**    -----------  
**
*******************************************************************************
**
** void transform_feauveau_trou (Ifloat &Image, MultiResol &MR_Transf)
**
** computes the wavelet transform by feauveau's algorithm without
** reducing the sampling
**
***************************************************************************/ 

// static char sccsid[] = "@(#)MR_Pav_Feau.cc 3.2 96/06/13 CEA 1994 @(#)";


#include <stdio.h>
#include <math.h>
#include <string.h>

#include "MR_Obj.h"
#include "IM_IO.h"

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

static void convolve_filter (Ifloat &Image, Ifloat &Image_Out, float *Filter, 
                             int Size_Filter, int Step_trou)
{ 
    int i,j,k,l;
    int Half_Filter = Size_Filter/2;
    int Nl = Image.nl();
    int Nc = Image.nc();
    int Step;

    Step = (int)(pow((double)2., (double) Step_trou) + 0.5);

    for (i = 0; i < Nl; i++) 
    for (j = 0; j < Nc; j++) 
    {
        Image_Out(i,j) = 0.;
        float *Ptr_Filter = Filter;
        for (k = -Half_Filter; k <= Half_Filter; k++)
        for (l = -Half_Filter; l <= Half_Filter; l++)
        {
          if (*Ptr_Filter != 0.)
              Image_Out(i,j) += Image(i+k*Step,j+l*Step,I_MIRROR)*(*Ptr_Filter);
          Ptr_Filter ++;
        }
    } 
}


/*************************************************************************/

void transform_feauveau_trou (Ifloat &Image, MultiResol &MR_Transf)
{
    int Nbr_Plan = MR_Transf.nbr_scale();
    int s;
 
    MR_Transf.scale(0) = Image;
    for (s = 0; s < Nbr_Plan-1; s++)
    {
        if (s%2)  convolve_filter (MR_Transf.scale(s), MR_Transf.scale(s+1), 
                               Filtre_U, DIM_FILTER, s/2);
        else      convolve_filter (MR_Transf.scale(s), MR_Transf.scale(s+1), 
                               Filtre_Diag_U, DIM_FILTER_DIAG, s/2);
    }

    for (s = 0; s < Nbr_Plan-1; s++) MR_Transf.scale(s) -= MR_Transf.scale(s+1);
}

/*************************************************************************/








