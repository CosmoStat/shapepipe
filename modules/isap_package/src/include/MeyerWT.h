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
**    File:  MeyerWT.h
**
**    Modification history :
**
*******************************************************************************
**
**    DESCRIPTION:  MEYER WAVELET TRANSFORM 
**    ----------- 
******************************************************************************/

#ifndef _MEYER_H
#define _MEYER_H

#include "GlobalInc.h"
#include "FFTN.h"
#include "FFTN_2D.h"
#include "IM_IO.h"

/***********************************************************************/

class MEYER_WT
{
  // Function used to compute the low pass filtering
  Ifloat H;  // H filter
  inline double lowpass_window(double x)
  {
     double l,r;
     l = exp( 1-1 / (1-exp(1-1/(1-x))));
     r = exp(1-1/(1-exp(1-1/x)));
     return (l /= sqrt(l*l+r*r));
  }

   int Nl_Imag;     // Input image number of lines
   int Nc_Imag;     // Input image number of col.
   int NbrScale;    // Nbr scales in the WT

   intarray TabNl;  // TabNl[s] = number of lines at scale s
   intarray TabNc;  // TabNc[s] = number of columns at scale s
   
   Bool NeedOddSize; // IF NeedOddSize==True, all scales will have odd sizes
                     // In this case, the input image must have an odd size

   Bool Extend;     // The image can be extended in order to avoid aliasing
                    // If Extend == True, extend the image NlxNc to an
		    // image 4/3Nl x 4/3 Nc, and take the WT of the extended
		    // image.
   int ExtNl;       // Image size of the image to be transformed
   int ExtNc;         
   double D_ExtNl;  // Image size of the image to be transformed
                    // D_ExtNl = Nl if Extend==False and 4./3. Nl otherwise
   double D_ExtNc;  // D_ExtNc = Nc if Extend==False and 4./3. Nc otherwise 

    Bool Isotrop; // If Isotrop == True, the low pass filter is isotropic. 
                  //    and it is not anymore a Meyer Wavelets.
		  //    The anti-aliasing does not work anymore and the parameter
		  //    extended must be equals to False.
		  
   void get_hfilter(Ifloat &H, double DNl=0., double DNc=0.);
         // Calculate the low pass filtering in the Fourier domain to be 
	 // applied at each step of the wavelet decomposition.

    void get_extFourier(Ifloat &Data, Icomplex_f & TF_ExtData);
         //  Take the FFT of the input image: Data
	 //  Extend it by periodization if Extend == True
	 //  TF_ExtData contains the FFT of the [extended] input image.

    void get_IFFT_Ima(Icomplex_f & TF_ExtData, Ifloat &Data);
         // Take the inv. FFT of the input  TF_ExtData,
	 // with extraction of part of interest if Extend == True

    void ifft_tabima(Icomplex_f * & TabCF_Ima, Ifloat * & Tab_Ima, Bool Alloc=True);
         // Take the inverse FFT of each wavelet scale
	 // If Alloc == True, Tab_Ima is allocated
	 //                 otherwise it must have been allocated before  the call
         // IN: TabCF_Ima[s] = Fourier transform of the scale s
	 // OUT: Tab_Ima[s] = scale s in direct space
	 
    void fft_tabima(Ifloat * & Tab_Ima, Icomplex_f * & TabCF_Ima); 
         // Take the Forward  FFT of each wavelet scale
         // IN: Tab_Ima[s] = scale s in direct space
	 // OUT: TabCF_Ima[s] = Fourier transform of the scale s

    void transform_cf(Icomplex_f & TF_ExtData, Icomplex_f * & TabWT);
         // Decomposition in Fourier space of the input TF_ExtData
	 // TabWT: output: = TabWT[s] = Fourier transform of the scale s
    
    void recons_cf(Icomplex_f * &TabWT,  Icomplex_f & TF_ExtData);
         // Reconstruction the Fourier transform of an image from
	 // the Fourier transforms of its scales.
	 // TabW: input: = TabWT[s] = Fourier transform of the scale s
	 // TF_ExtData: output reconstructed image in Fourier space
           
  public:
    inline int nlima() {return Nl_Imag;}
    inline int ncima() {return Nc_Imag;}
    inline double nld() {return D_ExtNl;}
    inline double ncd() {return D_ExtNc;}
    inline int nls(int s) {return TabNl(s);}
    inline int ncs(int s) {return TabNc(s);}
    inline Bool isotrop() {return Isotrop;}
    inline Bool extend() {return Extend;}
	float get_norm(int s);		// gives the normalizing coefficient at scale s
    
    FFTN_2D FFT2D;   // FTT2D class

    Icomplex_f TF_ExtData;     // Fourier transform of the extended image
                               // If Extend = False then TF_ExtData = FFT(input_data)
			       
    Icomplex_f *Tabcf_WT_Band; // Wavelet transform of the image in Fourier domain
    Ifloat *TabBand;           // Wavelet transform of the image in direct space
    
    Bool Verbose; // Verbose Mode
 		   
    int nbr_scale() {return NbrScale;} // Number of scales of the WT
    
    MEYER_WT() {FFT2D.CenterZeroFreq = True; Verbose = False; Extend=False; 
                 Isotrop=False; NbrScale=0;Tabcf_WT_Band=NULL;}
    
    void init(int Nbr_Scale, int Nl, int Nc, Bool ExtendWT=False, Bool IsotropWT=False, Bool WTNeedOddSize=False);
    // Initialize the WT class for a given number of scales and a give image size
    
    void transform(Ifloat &Data)
    {
       get_extFourier(Data, TF_ExtData);
       transform_cf(TF_ExtData, Tabcf_WT_Band);
    }
    // Computes the wavelet transform in Fourier domain of an image
    // Data: in = input image in direct space
    // TabWT: out = TabWT[s] is the Fourier transform of the wavelet scale s
    
    void transform(Ifloat &Data, Ifloat * & Tab_Ima, Bool Alloc=True)
    {
        transform(Data);
        ifft_tabima(Tabcf_WT_Band, Tab_Ima, Alloc);
    }
    // Computes the wavelet transform  of an image
    // Data: in = input image
    // TabWT: out = TabWT[s] is the scale s of the WT
    // If Alloc == True then Tab_Ima is allocated
    
    void recons(Ifloat * & Tab_Ima, Ifloat &Data, Bool Alloc=True)
    {
       if (Alloc == True) Data.resize(Nl_Imag, Nc_Imag);
       fft_tabima(Tab_Ima, Tabcf_WT_Band);
       recons(Data);
    }
    // Reconstruction of  an image from its WT.
    // Tab_Ima: in = input wavelet scales
    // Data: out = output reconstructed  image
    // If Alloc == True then Data is allocated
    
     void recons(Ifloat &Data)
     {
        recons_cf(Tabcf_WT_Band, TF_ExtData);
        get_IFFT_Ima(TF_ExtData, Data);
     }
    // Reconstruction of an image from the Fourier transforms of its wavelets scales
    // TabWT: in = input wavelet scales in Fourier space
    // Data: out = output reconstructed  image
    
    ~MEYER_WT() {if (Tabcf_WT_Band !=NULL) delete [] Tabcf_WT_Band;}
};
 
/***********************************************************************/

#endif
