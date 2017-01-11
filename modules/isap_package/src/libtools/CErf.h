/******************************************************************************
**                   Copyright (C) 1998 by CEA
*******************************************************************************
**
**    UNIT
**
**    Version: 1.0
**
**    Author: Jean-Luc Starck
**
**    Date:  98/01/12 
**    
**    File:  Cerf.h
**
*******************************************************************************
**
**    DESCRIPTION  : class for erf and erfc pre-calculation 
**    ----------- 
**                
**
**
******************************************************************************/

#ifndef _CERF_H_
#define _CERF_H_

#include "NR.h"


/* =============== ERF and ERFC ===============================*/

const float DEF_ERF_STEP = 0.01; // default Step  between two values
const float DEF_ERF_MAX = 3.5;   // default maximum value pre-calculated

class CERF {
    double Step;      // Step  between two values
    double *TabErrF;  // table for pre-calculated value

    int SizeTab;      // size of the table
    double MaxVal;    // maximum value pre-calculated in the table
    
   public:
    CERF()
    {
       double Val=0.;
       int i=0;
       Step = DEF_ERF_STEP;
       MaxVal = DEF_ERF_MAX;
       SizeTab = (int) (MaxVal / Step + 1.5);
       TabErrF= new double[SizeTab];
       while (Val < MaxVal)
       {
          TabErrF[i++] = erf(Val);
          Val += Step;
       }
    }
    
    // erf function
    inline double get_erf(double Val)
    {
        int Ind;
        double ValRet = 0.;
        
        if (Val >= MaxVal) ValRet = 1.;
        else if (Val > 0)
        {
            Ind = (int) (Val / Step + 0.5);
            ValRet = TabErrF[Ind];
        }
        return ValRet;
    }
    
    // erfc function
    inline double get_erfc(double Val) { return (1 - erf(Val)); }
    
    ~CERF()
    {
       if (TabErrF != NULL) delete [] TabErrF;
       SizeTab = 0;
    }
  };
  
#endif
