/******************************************************************************
**                   Copyright (C) 2004 by CEA
*******************************************************************************
**
**    UNIT
**
**    Version: 1.0
**
**    Author: Jean-Luc Starck
**
**    Date:  29/11/04 
**    
**    File:  Pyr1D.cc
**
*******************************************************************************
**
**    DESCRIPTION   routines concerting the pyramidal WT in direct space 
**    -----------       
**
**
******************************************************************************/ 

 
#ifndef _WT1D_PYR_H_
#define _WT1D_PYR_H_

#include <math.h>
#include <stdio.h>
#include <string.h> 

#include "IM_Obj.h"
#include "IM_IO.h"
 

enum type_pyr_filter {WT_PYR_BSPLINE,WT_PYR_LINEAR,WT_PYR_MEDIAN};


/****************************************************************************/

class PYR1D_WT {
    
   int TypeWT;  // Type of wavelet transform 
 
   fltarray Ima_h; // buffer for high frequency component
   fltarray Ima_b; // buffer for low frequency component
   fltarray Buff;   // buffer
   fltarray Data_Sub;
   int pyr_tab_pos(int NbrScale, int N);
                       // initialize TabSize and TabPos for a WT transform 
                       // with NbrScale scales and for a data set with  N points
                       // return the number of pixels of the transform

   intarray TabPos;    // position table: TabPos(s) = starting position of
                       // band s in the transform signal

   intarray TabSize;  // size table: :  TabSize(s) = size of
                      // band s in the transform signal

   int PyrPixNbr;     // Number of pixels of the transform
   void extand_sample (fltarray &Data, fltarray &Data_Ext);
   void reduc_sample (fltarray &Data, fltarray &Data_Sub);
   void smooth_b3spline(fltarray &Data, fltarray &Res);
   void smooth_b3spline(fltarray &Data);
   int Nbr_Plan;
   int Np;

   public:
    Bool ModifiedPWT;  
    type_border Border;
    int WinSize;
    int pyr_np () { return  PyrPixNbr;} 
            // return the number of pixels of the transform

    int pyr_pos (int s) { return TabPos(s);}
            // return the position of the band s

    int pyr_size (int s) { return TabSize (s);}
            // return the size of the band s

     // Constructor
    PYR1D_WT() {TypeWT=WT_PYR_BSPLINE; PyrPixNbr=0;ModifiedPWT=False;Nbr_Plan=Np=0;}
    PYR1D_WT(int TWT) {TypeWT=TWT; PyrPixNbr=0;ModifiedPWT=False;Nbr_Plan=Np=0;}
    PYR1D_WT(int TWT, int NbrScale, int N)  { alloc(TWT, NbrScale,N);}
    void alloc(int TWT, int NbrScale, int N);

    void transform (fltarray &Signal, fltarray &Trans);
    // WT transform of a signal with Nbr_Plan scales

    void recons (fltarray  &Trans, fltarray &Signal);
    // WT reconstruction of a signal from its WT
};

/****************************************************************************/

#endif
