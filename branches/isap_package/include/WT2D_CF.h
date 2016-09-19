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
**    Date:  25/09/2003 
**    
**    File:  WT2D_CF.h
**
**    Modification history :
**
*******************************************************************************
**
**    DESCRIPTION:  COMPLEX WAVELET TRANSFORM 
**    ----------- 
******************************************************************************/

#ifndef _WTCF2D_H
#define _WTCF2D_H

#include "Border.h"
#include "OptMedian.h"
#include "Filter.h"

/***********************************************************************/

class WTCF2D {
      SubBandFilter *SBF_ODD;
      SubBandFilter *SBF_EVEN;
      Ortho_2D_WT *TreeA;
      Ortho_2D_WT *TreeB;
      Ortho_2D_WT *TreeC;
      Ortho_2D_WT *TreeD;
      HALF_DECIMATED_2D_WT *WT;
      void init();
      public:
       Bool Verbose;
       WTCF2D(){init();}
        void get_z1_re_im(float a, float b, float c, float d, float & zr,  float & zi)
              {zr = a-d; zi = c+b;}
        void get_z2_re_im(float a, float b, float c, float d, float & zr,  float & zi)
              {zr = a+d; zi = b-c;}
        float get_z1(float a, float b, float c, float d)
              {return sqrt( (a-d)*(a-d) + (c+b)*(c+b)); }
        float get_z2(float a, float b, float c, float d)
              {return sqrt( (a+d)*(a+d) + (b-c)*(b-c)); }
	void abcd_to_z(float a, float b, float c, float d, float & z1r,  float & z1i, float & z2r,  float & z2i)
              {get_z1_re_im(a,b,c,d,z1r,z1i);get_z2_re_im(a,b,c,d,z2r,z2i);}      
        void z_to_abcd(float z1r, float z1i, float z2r, float z2i, float & a,  float & b, float & c,  float & d)
	      { a =0.5*(z1r+z2r); b=0.5*(z1i+z2i); c=0.5*(z1i-z2i); d=0.5*(z2r-z1r);}
	~WTCF2D(){delete SBF_ODD; delete SBF_EVEN;delete TreeA;
                delete  TreeB; delete TreeC; delete TreeD;}
       void transform(Ifloat & Data,  Ifloat * & TabTrans, int NbrScale);
       void coeff_to_z(Ifloat *TabTrans, Ifloat *TabZre, Ifloat *TabZim, int Nbr_Plan, int NzBand);
       void alloc_tabz(int Nl, int Nc, Ifloat * & TabZre, Ifloat * & TabZim, int Nbr_Plan, int & NzBand);
       void threshold(Ifloat * & TabZre, Ifloat * & TabZim, int Nz, float SigmaNoise);
       void z_to_coeff(Ifloat * & TabZre, Ifloat * & TabZim, Ifloat *TabTrans, int Nbr_Plan);
       void recons(Ifloat * & TabTrans, Ifloat & Data, int NbrScale);
};

/***********************************************************************/

#endif
