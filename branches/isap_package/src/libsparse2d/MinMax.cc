/*******************************************************************************
**
**    UNIT
**
**    Version: 1.0
**
**    Author: Jean-Luc Starck
**
**    Date:  5/03/2000 
**    
**    File:  MinMax.cc
**
**    Modification history:
**
*******************************************************************************
**
**    DESCRIPTION  1D min-max Sub-Band decomposition
**    -----------  
**                 
******************************************************************************/
 

#include "SB_Filter1D.h"


/*************************************************************/
// G MIN MAX
/*************************************************************/


void G_Minmax::transform(int N, float *High, float *Low, float *Det)
{
   int i;
   for (i = 0; i < N; i += 2)
   {
      int Index = test_index(i+1, N);
      float MinVal = MIN(High[i], High[Index]);
      Det[i/2] =  High[i] - High[Index];
      Low[i/2] =  MinVal;
   }
}

/************************************************************************/

void G_Minmax::recons(int N, float *Low, float *Det, float *High)
{
   int i;    
   for (i = 1; i < N; i += 2)
   {
      if (Det[i/2] >= 0) 
      {
         High[i-1] = Low[i/2] + Det[i/2];
	 High[i] = Low[i/2];
      }
      else
      {
         High[i] = Low[i/2] - Det[i/2];
	 High[i-1] = Low[i/2];
      }
   }
   if (N % 2 == 1) High[N-1] = Low[N/2];
}

/************************************************************************/

void G_Minmax::transform(int N, float *High, float *Low, float *Det, int Step)
{
   cout << "Error: Pave Minmax not implemented ... " << endl;
   exit(-1);
}

/************************************************************************/

void G_Minmax::recons(int N, float *Low, float *Det, float *High, int Step)
{
   cout << "Error: Pave Minmax not implemented ... " << endl;
   exit(-1);
}

/************************************************************************/
