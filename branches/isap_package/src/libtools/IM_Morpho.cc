/******************************************************************************
**                   Copyright (C) 1994 CEA
*******************************************************************************
**
**    UNIT
**
**    Version: 3.2
**
**    Author: Jean-Luc Starck and Simona Mei
**
**    Date:  96/06/13
**    
**    File:  IM_Morpho.cc
**
*******************************************************************************
**
**    DESCRIPTION  morphological tools for image processing
**    -----------  
**
******************************************************************************* 
**
** morpho_erosion (Ifloat &Imag1, Ifloat& Imag2)
**
** Imag2 = erosion(Imag1) with a square filter 3x3
** erosion in gray level
**
******************************************************************************* 
**
** morpho_dilation (Ifloat &Imag1, Ifloat &Imag2)
**
** Imag2 = dilation(Imag1) with a square filter 3x3
** dilatation in gray level
** 
*******************************************************************************
** 
** morphoi_erosion (int *Imag1, int *Imag2, int Nl, int Nc, int WindowSize)
**
** Imag2 = erosion(Imag1) with a square filter WSxWS
** erosion in gray level
**
******************************************************************************* 
**
** morphoi_dilation (int *Imag1, int *Imag2, int Nl, int Nc, int WindowSize)
**
** Imag2 = dilation(Imag1) with a square filter WSizexWindowSize
** dilatation in gray level
**
******************************************************************************* 
** Structurant element is a cercle
**
** morpho_cerclei_erosion (int *Imag1, int *Imag2, int Nl, int Nc,int WindowSize)
**
** Imag2 = erosion(Imag1) with a square filter WSxWS
** erosion in gray level
**
******************************************************************************* 
**
** morpho_cerclei_dilation (int *Imag1, int *Imag2, int Nc, Nl,int WindowSize)
**
** Imag2 = dilation(Imag1) with a square filter WSizexWindowSize
** dilatation in gray level
**
******************************************************************************* 
**
** morpho_cercle_erosion (Ifloat &Imag1, Ifloat& Imag2)
**
** Imag2 = erosion(Imag1) with a square filter 3x3
** erosion in gray level
**
******************************************************************************* 
**
** morpho_cercle_dilation (Ifloat &Imag1, Ifloat &Imag2)
**
** Imag2 = dilation(Imag1) with a square filter 3x3
** dilatation in gray level
**
******************************************************************************/ 
  
 
// static char sccsid[] = "@(#)IM_Morpho.cc 3.2 96/06/13 CEA @(#)";


#include "IM_Obj.h"


/***************************************************************************/

void morpho_erosion (Ifloat &Imag1, Ifloat& Imag2, int Window_Size)
{
    int i,j,k,l;
    float Min;
    int Nl = Imag1.nl();
    int Nc = Imag1.nc();
    int Window2 = (Window_Size - 1) / 2;

    for (i = 0; i < Nl; i++) 
    for (j = 0; j < Nc; j++) 
    {
        Min = Imag1(i,j);
        for (k = i - Window2; k <= i + Window2; k++)
        for (l = j - Window2; l <= j + Window2; l++)
               if (Imag1(k,l, I_CONT) < Min) Min = Imag1(k,l, I_CONT);
        Imag2(i,j) = Min;
    }
}

/***************************************************************************/
 
void morpho_dilation (Ifloat &Imag1, Ifloat &Imag2, int Window_Size)
{
    int i,j,k,l;
    float Max;
    int Nl = Imag1.nl();
    int Nc = Imag1.nc();
    int Window2 = (Window_Size - 1) / 2;

    for (i = 0; i < Nl; i++) 
    for (j = 0; j < Nc; j++) 
    { 
        Max = Imag1(i,j);
        for (k = i - Window2; k <= i + Window2; k++)
        for (l = j - Window2; l <= j + Window2; l++)
            if (Imag1(k,l, I_CONT) > Max) Max = Imag1(k,l, I_CONT);
        Imag2(i,j) = Max;
    }
}

/***************************************************************************/

void morpho4_erosion (Ifloat &Imag1, Ifloat& Imag2, int Step_trou)
{
    int i,j;
    float Min;
    int Nl = Imag1.nl();
    int Nc = Imag1.nc();
    int Step;

    Step = (int)(pow((double)2., (double) Step_trou) + 0.5);

    for (i = 0; i < Nl; i++) 
    for (j = 0; j < Nc; j++) 
    {
        Min = Imag1(i,j);
        if (Min > Imag1(i+Step,j, I_CONT) ) Min = Imag1(i+Step, j, I_CONT);
        if (Min > Imag1(i-Step,j, I_CONT) ) Min = Imag1(i-Step,j, I_CONT);
        if (Min > Imag1(i,j+Step, I_CONT) ) Min = Imag1(i,j+Step, I_CONT);
        if (Min > Imag1(i,j-Step, I_CONT) ) Min = Imag1(i,j-Step, I_CONT);
        Imag2(i,j) = Min;
    }
}

/***************************************************************************/
 
void morpho4_dilation (Ifloat &Imag1, Ifloat &Imag2, int Step_trou)
{
    int i,j;
    float Max;
    int Nl = Imag1.nl();
    int Nc = Imag1.nc();
    int Step;

    Step = (int)(pow((double)2., (double) Step_trou) + 0.5);

    for (i = 0; i < Nl-Step_trou; i++) 
    for (j = 0; j < Nc-Step_trou; j++) 
    {
        Max = Imag1(i,j);
        if (Max < Imag1(i+Step,j, I_CONT) ) Max = Imag1(i+Step, j, I_CONT);
        if (Max < Imag1(i-Step,j, I_CONT) ) Max = Imag1(i-Step,j, I_CONT);
        if (Max < Imag1(i,j+Step, I_CONT) ) Max = Imag1(i, j+Step, I_CONT);
        if (Max < Imag1(i,j-Step, I_CONT) ) Max = Imag1(i, j-Step, I_CONT);
        Imag2(i,j) = Max;
    }
}

/***************************************************************************/

static void i_clean_border(int *Imag, int Nl, int Nc, int Bord)
{
   int i,j;
   for (i=0; i < Bord; i++) 
   for (j=0; j < Nc; j++)
   {
      Imag[i*Nc+j] = 0;
      Imag[(Nl-i-1)*Nc+j] = 0;
   }
   for (i=0; i < Nl; i++) 
   for (j=0; j < Bord; j++)
   {
      Imag[i*Nc+j] = 0;
      Imag[i*Nc+Nc-j-1] = 0;
   }
}

/********************************************************************/

void morphoi_erosion (int *Imag1, int *Imag2, int Nl, int Nc, int ElStr_Size)
{
    int i,j,k,l;
    int Min;
    int Window2 = (ElStr_Size - 1) / 2;
         
    i_clean_border(Imag2,Nl,Nc,Window2+1);

    for (i = Window2+1; i < Nl-(Window2+1); i++) 
    for (j = Window2+1; j < Nc-(Window2+1); j++) 
    {
        Min = Imag1[i*Nc+j];

        for (k = i - Window2+1; k <= i + Window2-1; k++)
        for (l = j - Window2+1; l <= j + Window2-1; l++)
               if (Imag1[k*Nc+l] < Min) Min = Imag1[k*Nc+l];

       Imag2[i*Nc+j] = Min;
    }
 }



/***************************************************************************/

void morphoi_dilation (int *Imag1, int *Imag2, int Nl, int Nc, int ElStr_Size)
{
    int i,j,k,l;
    int Max;
    int Window2 = (ElStr_Size - 1) / 2;
         
    i_clean_border(Imag2,Nl,Nc,Window2+1);

    for (i = Window2+1; i < Nl-(Window2+1); i++) 
    for (j = Window2+1; j < Nc-(Window2+1); j++) 
    {
        Max = Imag1[i*Nc+j];

        for (k = i - Window2+1; k <= i + Window2-1; k++)
        for (l = j - Window2+1; l <= j + Window2-1; l++)
            if (Imag1[k*Nc+l] > Max) Max = Imag1[k*Nc+l];

       Imag2[i*Nc+j] = Max;
    }
 }

/***************************************************************************/

void morpho_cercle_erosion (Ifloat &Imag1, Ifloat& Imag2, int ElStr_Size)
{
    int i,j,k,l;
    float Min;
    int Nl = Imag1.nl();
    int Nc = Imag1.nc();
    int Window2 = (ElStr_Size - 1) / 2;

    for (i = 0; i < Nl; i++) 
    for (j = 0; j < Nc; j++) 
    {
        Min = Imag1(i,j);

//internal part of the cercle

        for (k = i - Window2+1; k <= i + Window2-1; k++)
        for (l = j - Window2+1; l <= j + Window2-1; l++)

            if (Imag1(k,l, I_CONT) < Min) Min = Imag1(k,l, I_CONT);

//bords
         
        for (l = j - Window2+1; l <= j + Window2-1; l++)
      if (Imag1(i-Window2,l, I_CONT) < Min) Min = Imag1(i-Window2,l, I_CONT);

        for (l = j - Window2+1; l <= j + Window2-1; l++)
      if (Imag1(i+Window2,l, I_CONT) < Min) Min = Imag1(i+Window2,l, I_CONT);

        for (k = i - Window2+1; k <= i + Window2-1; k++)
      if (Imag1(k,j-Window2, I_CONT) < Min) Min = Imag1(k,j-Window2, I_CONT);

        for (k = i - Window2+1; k <= i + Window2-1; k++)
      if (Imag1(k,j+Window2, I_CONT) < Min) Min = Imag1(k,j+Window2, I_CONT);

// fin bords

        Imag2(i,j) = Min;
               
    }
}

/***************************************************************************/
 
void morpho_cercle_dilation(Ifloat &Imag1, Ifloat &Imag2, int ElStr_Size)
{
    int i,j,k,l;
    float Max;
    int Nl = Imag1.nl();
    int Nc = Imag1.nc();
    int Window2 = (ElStr_Size - 1) / 2;

    for (i = 0; i < Nl; i++) 
    for (j = 0; j < Nc; j++) 
    { 
        Max = Imag1(i,j);

//internal part of the cercle

        for (k = i - Window2+1; k <= i + Window2-1; k++)
        for (l = j - Window2+1; l <= j + Window2-1; l++)

            if (Imag1(k,l, I_CONT) > Max) Max = Imag1(k,l, I_CONT);

//bords
         
        for (l = j - Window2+1; l <= j + Window2-1; l++)
      if (Imag1(i-Window2,l, I_CONT) > Max) Max = Imag1(i-Window2,l, I_CONT);

        for (l = j - Window2+1; l <= j + Window2-1; l++)
      if (Imag1(i+Window2,l, I_CONT) > Max) Max = Imag1(i+Window2,l, I_CONT);

        for (k = i - Window2+1; k <= i + Window2-1; k++)
      if (Imag1(k,j-Window2, I_CONT) > Max) Max = Imag1(k,j-Window2, I_CONT);

        for (k = i - Window2+1; k <= i + Window2-1; k++)
      if (Imag1(k,j+Window2, I_CONT) > Max) Max = Imag1(k,j+Window2, I_CONT);

// fin bords

        Imag2(i,j) = Max;
    }
}

/********************************************************************/


void morphoi_cercle_erosion (int *Imag1, int *Imag2, int Nl, int Nc, int ElStr_Size)
{
    int i,j,k,l;
    int Min;
    int Window2 = (ElStr_Size - 1) / 2;
     
    i_clean_border(Imag2,Nl,Nc,Window2+1);
     
    for (i = Window2+1; i < Nl-(Window2+1); i++) 
    for (j = Window2+1; j < Nc-(Window2+1); j++) 
    {
        Min = Imag1[i*Nc+j];

       //internal part of the cercle

        for (k = i - Window2+1; k <= i + Window2-1; k++)
        for (l = j - Window2+1; l <= j + Window2-1; l++)
                if (Imag1[k*Nc+l] < Min) Min = Imag1[k*Nc+l];
 
       //bords
         
        for (l = j - Window2+1; l <= j + Window2-1; l++)
      if (Imag1[(i-Window2)*Nc+l] < Min) Min = Imag1[(i-Window2)*Nc+l];

        for (l = j - Window2+1; l <= j + Window2-1; l++)
      if (Imag1[(i+Window2)*Nc+l] < Min) Min = Imag1[(i+Window2)*Nc+l];

        for (k = i - Window2+1; k <= i + Window2-1; k++)
      if (Imag1[k*Nc+(j-Window2)] < Min) Min = Imag1[k*Nc+(j-Window2)];

        for (k = i - Window2+1; k <= i + Window2-1; k++)
      if (Imag1[k*Nc+j+Window2] < Min) Min = Imag1[k*Nc+j+Window2];

      // end bords circle

        Imag2[i*Nc+j] = Min;
    }
    
    
 }

/***************************************************************************/

void morphoi_cercle_dilation (int *Imag1, int *Imag2, int Nl, int Nc, int ElStr_Size)
{
    int i,j,k,l;
    int Max;
    int Window2 = (ElStr_Size - 1) / 2;
    
    i_clean_border(Imag2,Nl,Nc,Window2+1);
       
    for (i = Window2+1; i < Nl-Window2-1; i++) 
    for (j = Window2+1; j < Nc-Window2-1; j++) 
    {
        Max = Imag1[i*Nc+j];


        //internal part of the cercle

        for (k = i - Window2+1; k <= i + Window2-1; k++)
        for (l = j - Window2+1; l <= j + Window2-1; l++)
             if (Imag1[k*Nc+l] > Max) Max = Imag1[k*Nc+l];

       //bords
       for (l = j - Window2+1; l <= j + Window2-1; l++)
      if (Imag1[(i-Window2)*Nc+l] > Max) Max = Imag1[(i-Window2)*Nc+l];

        for (l = j - Window2+1; l <= j + Window2-1; l++)
      if (Imag1[(i+Window2)*Nc+l] > Max) Max = Imag1[(i+Window2)*Nc+l];

        for (k = i - Window2+1; k <= i + Window2-1; k++)
      if (Imag1[k*Nc+(j-Window2)] > Max) Max = Imag1[k*Nc+(j-Window2)];

        for (k = i - Window2+1; k <= i + Window2-1; k++)
      if (Imag1[k*Nc+j+Window2] > Max) Max = Imag1[k*Nc+j+Window2];


     Imag2[i*Nc+j] = Max;
   }
}

/***************************************************************************/

