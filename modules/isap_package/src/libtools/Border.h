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
**    Date:  96/05/07 
**    
**    File:  Border.h
**
*******************************************************************************
**
**    DESCRIPTION  
**    ----------- 
**                 
**    PARAMETRES    
**    ----------    
** 
**    RESULTS      
**    -------  
**
**
******************************************************************************/

#ifndef _BORDER_H_
#define _BORDER_H_

#include<stdio.h>
#include<stdlib.h>

#define NBR_BORD 4

enum type_border{I_CONT, I_MIRROR, I_PERIOD, I_ZERO};

#define DEFAULT_BORDER I_CONT

inline int test_index_cont(int i, int N)
{
    int indi = i;
    if (i < 0) indi = 0;
    else if (i >= N) indi = N - 1;
    return indi;
}
inline int test_index_mirror(int i, int N)
{
    int indi = i;
    if (i < 0)
    {
        indi = - i;
	if (indi >= N) indi = N-1;
    }
    else
     if (i >= N)
     {
         indi = 2 * (N - 1) - i;
	 if (indi < 0) indi = 0;
     }
    return indi;
}

inline int test_index_period(int i, int N)
{
    int indi = i;
    if (i < 0) while (indi < 0) indi += N;
    else if (i >= N) while (indi >= N) indi -= N;
    return indi;
}


inline int get_index(int i, int N, type_border TB)
{
    int indi = i;
    switch (TB)
    {
      case I_CONT: indi = test_index_cont(i,N); break;
      case I_MIRROR: indi = test_index_mirror(i,N); break;
      case I_PERIOD: indi = test_index_period(i,N); break;
      case I_ZERO:  
      default:
         printf("Error: bad parameter bord in  get_index");
         break;
    } // end case
    return indi;
}


#endif

