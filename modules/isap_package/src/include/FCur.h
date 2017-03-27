/******************************************************************************
**                   Copyright (C) 2005 by CEA
*******************************************************************************
**
**    UNIT
**
**    Version: 1.0
**
**    Author: Jean-Luc Starck
**
**    Date:  21/01/2005 
**    
**    File:  FCur.h
**
**    Modification history :
**
*******************************************************************************
**
**    DESCRIPTION:  Fast Curvelet Transform 
**    ----------- 
******************************************************************************/

#ifndef _FASTCUR_H
#define _FASTCUR_H

#include "GlobalInc.h"
#include "IM_IO.h"
#include "FFTN.h"
#include "FFTN_2D.h"
#include "MeyerWT.h"
#include "IM_Prob.h"
#include "IM_Math.h"
#include <fstream>
/***********************************************************************/

int fct_real_get_band_size(int Nbr_Scale, int Nl, int Nc, int NDir, intarray &TabSizeNl, intarray &TabSizeNc);


class FCUR: public MEYER_WT
{
protected:
    int DataNl;  // Input image size must be odd. If it is not, we etend it by one line or 
    int DataNc;  // one column, and we save here the original image size
    Bool ModifSize;
    int NewNl;
    int NewNc;
    
	bool isset_tabsigma;	// true if TabCurSigma is filled
    void get_wedges(Icomplex_f * &TabWT);
    void put_wedges(Icomplex_f * &TabWT);
    intarray *TabSizeNl;  // Band size 
    intarray *TabSizeNc;  // Band size 
    
    void get_size();

public:
	bool Undecimated;
   fltarray TabCurSigma;   // Normalization array
   Icomplex_f **TabCF_Band;  // TabCF_Band contains the curvelet transform of the Data
                             // TabCF_Band[s][b].real()  ==> real band at scale s and dir b
			     // TabCF_Band[s][b].imag()  ==> imaginary band at scale s and dir b
			     //    real part and imaginary part can be analyzed independently
   int NbrBandCF;            // Total Number of complex bands
   int NbrBand;              // total Number of Band 
   intarray TabNbrBandPerScale;   // Number of   bands per scale
   intarray TabNbrAnglePerScale;  // Number of directions per scale
   Bool RealBand;                 // If true, the redundancy inside the curvelet transform is used to remove
                                  // the imaginary part.Then the imaginary part equals to zero.
  // public:
    int nbr_dir(int s) {return TabNbrAnglePerScale(s);}
    
    FCUR():MEYER_WT() {TabCF_Band=NULL;NbrBand=0;NbrBandCF=0;RealBand=False;isset_tabsigma=false;Undecimated=false;}
    ~FCUR();
	
    void alloc_with_tab(int Nbr_Scale, int Nl, int Nc, intarray & TabDir, Bool ExtendWT=True, Bool IsotropWT=False, Bool RealCur=False);
    // Allocation of the class
    // Nbr_Scale = number of scales used in the wavelet transform
    // Nl,Nc = input image size
    // TabDir = array of number of directions per scale (minimum 8)
    // ExtendWT = The image is extended for aliasing removal
    // IsotropWT = an isotropic wavelet transform is used instead of the meyer wavelet
    //             in this case ExtendWT is set to false
        
    void alloc_from_coarse(int Nbr_Scale, int Nl, int Nc, int NbrDir, Bool ExtendWT=True, Bool IsotropWT=False, Bool RealCur=False);
    // NbrDir = number of directions at the coarsest resolution (minimum 8)

    
    void alloc_from_fine(int Nbr_Scale, int Nl, int Nc, int NbrDir, Bool ExtendWT=True, Bool IsotropWT=False, Bool RealCur=False);
        // NbrDir = number of directions at the finest resolution (minimum 8)
 
    void cur_write(char *Name);
    void cur_trans(Ifloat &Data);
    void cur_recons(Ifloat &Data);
	
    inline int real() { return RealBand;}
    inline int nbr_band(int s) { return  TabNbrBandPerScale(s);}
    inline int nbr_tot_band() { return  (int) TabNbrBandPerScale.total();}
    // return the  number of bands at scale s
    
    //inline int size_band_nl(int s, int b) { return  (real() == False) ? TabCF_Band[s][b/2].nl() : TabCF_Band[s][b].nl();}
    // inline int size_band_nc(int s, int b) { return  (real() == False) ? TabCF_Band[s][b/2].nc() : TabCF_Band[s][b].nc();}
    inline int size_band_nl(int s, int b) { return  (real() == False) ? TabSizeNl[s](b/2) : TabSizeNl[s](b);}
    inline int size_band_nc(int s, int b) { return  (real() == False) ? TabSizeNc[s](b/2) : TabSizeNc[s](b);}
      
    void get_band(int s, int b, Ifloat &Band);
     // return a given band
    
    void put_band(int s, int b, Ifloat &Band);
    // return a given band
    
	// return a coefficient. If complex transform, returns the real part
    inline float &  operator() (int s, int b, int i, int j) 
    {
        int bb = (real() == False) ? b/2: b;
        complex_f *PtrCF = (TabCF_Band[s][bb]).buffer() + i * (TabCF_Band[s][bb]).nc() + j;
        float *Ptr = (float *) PtrCF;
		if ((real() == False) && (b % 2 != 0)) Ptr++;
        return *Ptr;
    }
     
    void get_stat(fltarray &TabStat);
    void get_norm_coeff(Ifloat &ImaNoise, float N_Sigma);
    void set_noise_level(float N_Sigma);
    void get_norm_coeff(float N_Sigma);
	void import_norm_coef(fltarray * TabSigma);
	void export_norm_coef(fltarray * &TabSigma);
    void threshold(float SigmaNoise, float N_Sigma);
    void wiener_filter(float SigmaNoise, int WienerBlockSize=7);
    float norm_band(int s, int b) { return TabCurSigma(s,b);}
	bool isset_norm() {return isset_tabsigma;}
    void read(char *Name);
    void mr_io_fill_header(fitsfile *fptr);
    void write(char *Name);
};

 
/***********************************************************************/

#endif
