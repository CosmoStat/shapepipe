/******************************************************************************
**                   Copyright (C) 2003 by CEA
*******************************************************************************
**
**    UNIT
**
**    Version: 1.0
**
**    Author: Jean-Luc Starck
**
**    Date:  15/06/03 
**    
**    File:  MR_SoftRegul.h
**
*******************************************************************************
**
**    DESCRIPTION  multiresolution soft regularization
**    ----------- 
**
******************************************************************************/

#ifndef _MRSREGUL_H_
#define _MRSREGUL_H_

#include "MR_Obj.h"
#include "IM_Regul.h"
 

/***************************************/	     		     

void regul_tv_haarwt(Ifloat &Ima, Ifloat &Rec, float LambdaTV, int NbrScaleTV);
// TV regularization using the soft threshold of undecimated Haar transform.
// LambdaTV = soft threshold
// NbrScaleTV = number of scales in the regularization.
// Standard TV ==> 2 scales.

class MR_Regul
{
  public:
   int NbrScale; // Number of scales used in the soft thresholding
   int NbrUndecimatedScale;
                 // Number of undecimated scale. By default all scales are
		 // undecimated
   type_sb_filter TypeFilter;
                 // Type of filter. By default Haar filter is chosen, which
		 // means that SoftThresholding is equivalent to a TV constaint.
   Bool ExpDecreasingLambda;
                 // If true, the threshold descrease with the scale.
		 
   MR_Regul() {NbrScale=4; NbrUndecimatedScale=10; TypeFilter= F_HAAR; // F_HAAR,F_MALLAT_7_9; 
               ExpDecreasingLambda=False;}
   
   void im_soft_threshold(Ifloat &Ima, Ifloat &Rec, float Lambda, float NoiseLevel=0.);
   // Apply the undecimated bi-orthogonal WT to the image,
   // soft threshold the coefficients with the threshold level Lambda,
   // Reconstruct the filtered image in Rec
   void im_soft_non_ortho_threshold(Ifloat &Ima, Ifloat &Rec, float Lambda, float NoiseLevel=0.);
   void im_soft_iwt_threshold(Ifloat &Ima, Ifloat &Rec, float Lambda, float NoiseLevel=0.);

   void obj_regul(Ifloat &Obj, Ifloat &Grad, float  Lambda);
   // Modify Grad using a soft-thresholding on Obj.
   //  Aux = soft_threshold(Obj, Lambda)
   // Grad += Aux - Obj

   void ima_regul(Ifloat &Obj, Ifloat &Grad, float Lambda);
   // Calculate the gradient relative to the derivative of the TV.
   // Grad = Obj - soft_threshold(Obj, Lambda)
};
                       
#endif
