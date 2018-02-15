/******************************************************************************
**                   Copyright (C) 2000 by CEA
*******************************************************************************
**
**    UNIT
**
**    Version: 1.0
**
**    Author: Jean-Luc Starck
**
**    Date:  5/03/2000 
**    
**    File:  SB_Filter1D.h
**
**    Modification history :
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

#ifndef _SB_ATROU1D_H_
#define _SB_ATROU1D_H_

#include "Border.h"
#include "Filter.h"
 
 
/***********************************************************************
*  A Trous 1D wavelet transform
************************************************************************/

class ATROUS_1D_WT 
{
   private:
      void b3spline_filtering (fltarray& Sig_in, fltarray& Sig_out, 
                               int Step_trou);
      float update (float CoefSol, float Threshold, float SoftLevel);
      bool is_on_support (float CoefSol, float Threshold, float SoftLevel);

   public:
      Bool ModifiedAWT; // If true, 
      Bool AdjointRec;  // Use the adjoint for the reconstruction
      type_border Bord; // Type of Bord
       
      // Use for MGA1D (cb_mca1d) for the detect
      int NbScale;      // Number of scales
      int NbPts;        // Number of points
      Bool OnlyPositivDetect; // detect only positive coef.
      Bool SuppressIsolatedPixel; // suppress isolated pixels
      int FirstDetectScale;  // First detection scale
      Bool Verbose;          // Verbose mode
      Bool Write;            // Write some results
      float SigmaNoise;      // Noise standard deviation
      Bool UseNormL1;        // Use the L1 norm
      float NSigmaSoft;      // Nsgima used 
      
      ATROUS_1D_WT () {AdjointRec=False;Bord=I_CONT;ModifiedAWT=False;}
      ~ATROUS_1D_WT(){}
      
      void alloc (fltarray * & TabTrans, int Nx, int Nbr_Plan);
      // Allocate TabTrans with Nbr_Plan scales and Nx data points
      
      void transform (fltarray& Image, fltarray *TabBand);
      // Applied the wavelet transform
      
      void recons (fltarray* TabTrans, fltarray& Signal);
      // Reconstruc a signal from its WT
      
      void free (fltarray* TabTrans, int Nbr_Plan);
      // deallocate TabTrans
      
      void KillScaleNotUsed (fltarray* AT_Trans, int FirstDetectScale);
      // Set to zero the FirstDetectScale scales
      
      void KillLastScale (fltarray* AT_Trans);
      // Set to zero  the last scale
      
      float getAbsMaxTransf (fltarray* AT_Trans);
      // Return the coeff which has the maximum absolute value
      
      void Threshold (fltarray* AT_Trans, float NSigma, int IterNumber);
      // Threshold the wavelet coefficient

      float norm_band (int s);
      // Return the normalization value relative the band s
};


#endif
