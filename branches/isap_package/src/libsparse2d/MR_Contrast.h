/******************************************************************************
**                   Copyright (C) 2002 by CEA
*******************************************************************************
**
**    UNIT
**
**    Version: 1.0
**
**    Author: Jean-Luc Starck
**
**    Date:  18/01/02 
**    
**    File:  MR_Contrast.h
**
*******************************************************************************
**
**    DESCRIPTION  multiresolution constrast enhancement
**    ----------- 
**
******************************************************************************/

#ifndef _CMRCONTRAST_H_
#define _CMRCONTRAST_H_

#include "IM_Obj.h"
#include "IM_Math.h"
#include "MR_Obj.h"
#include "FFTN_2D.h"
#include "IM_Sigma.h"

 /***************************************/

#define NBR_CONTRAST 6
enum type_contrast {CT_HISTO, CT_WT, CT_ATROU_LOG, CT_ATROU_LOG_SGN, 
		    CT_CLIPPING, CT_LAPLACIAN, CT_RET, CT_MRET, CT_ATROU_RET, CT_CUR};
#define DEF_CONTRAST CT_WT

const char * StringContrast ( type_contrast  type);

/***************************************/	     		     

class Contrast
{
  public:
   // Contrast function, following K. Vande Velde (IEEE, 1999).
   // y(x) = (m/c)^p if |x| < c
   // y(x) = (m/|x|)^p if c <= |x| < m
   // y(x) = 1 if x >= m
   inline float contrast_function(float x)
   {
       float ValRet=0.;
       if (x > FLOAT_EPSILON)
       {
          if (ABS(x) < Contrast_C_Param) 
	    ValRet = pow( (double) (Contrast_M_Param/Contrast_C_Param), Contrast_P_Param);
          else if (ABS(x) < Contrast_M_Param) 
	    ValRet = pow( (double) (Contrast_M_Param/ABS(x)), Contrast_P_Param);
          else ValRet = 1.;
       }
       return ValRet;
   }
   
  Contrast() {Contrast_C_Param=0;Contrast_M_Param=100; 
                  Contrast_P_Param=0.3; Contrast_Q_Param = 0;  
 		  LogParam=0.1;Verbose=False;ClipVal=3.;LaplacianParam=1.;}
  double LogParam;

  // Parameter for the contrast function		       
  double Contrast_P_Param; // determine the degree of non-linearity
                           // p must be in ]0,1[
  double Contrast_Q_Param; // q must be in [-0.5,0.5] 
                           // q > 0, enhance less the darker part than
			   //        the lighter part
			   // q < 0 ==> enhance more the dark than lighter part
  double Contrast_M_Param; // Transform coefficients larger than m are not
                           // modified
  double Contrast_C_Param; // correspond to the noise level

  // Parameter for the sigma clipping
  float LaplacianParam;
  float ClipVal;
  float Noise_Ima;
  Bool Verbose;
  void retinex(Ifloat &Data);
  void multiscale_retinex(Ifloat &Data);
  void atrou_retinex(Ifloat &Data, int Nbr_Plan);
  void histo_equal(Ifloat &Data, int NbrBin=1024);
  void laplacian(Ifloat &Data);
  void wt_enhance(Ifloat &Data, int Nbr_Plan);
  void atrou_log(Ifloat &Data, int Nbr_Plan, Bool Sign);
  void sature_clipping(Ifloat &Data);
};
                       
#endif
