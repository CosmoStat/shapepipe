/******************************************************************************
**                   Copyright (C) 2000 by CEA
*******************************************************************************
**
**    UNIT
**
**    Version: 3.2
**
**    Author: Jean-Luc Starck
**
**    Date:  28/08/00
**    
**    File:  FFTN_2D.h
**
*******************************************************************************
**
**    DESCRIPTION  2D FFT Class routine for any size of images
**    -----------  fow square and POW2 image, call the fft_any_power_of_2 
**                 routine which is faster + better precision
**                 for other image sizes, use FFTN class
**                 The drawback is that fft_any_power_of_2 produce 
**                 a Fourier transfrom in which the zero frequency is
**                 the middle of the image, while FFTN zero frequency is
**                 at the bottom left. So we introduced the CenterZeroFreq
**                 parameter in order to specify which arrangement is prefered,
**                 and a swapping of the frequencies componenent is performed
**                 if needed (which adds some computation time).
**
**                 Using the optim option set the CenterZeroFreq in order
**                 not to swap the Fourier pixels for a given image size.
**                 
******************************************************************************/

#ifndef _CFFTN_2D_H_
#define _CFFTN_2D_H_

#include "GlobalInc.h"
#include "IM_Obj.h"
#include "FFTN.h"

class FFTN_2D: public FFTN {
   inline void pix_swap(complex_f &a, complex_f &b) { complex_f temp=a;a=b;b=temp;}
   inline void pix_swap(complex_d &a, complex_d &b) { complex_d temp=a;a=b;b=temp;}

   void center(Ifloat &Image);
   void center(Icomplex_d &Image);
   void uncenter(Ifloat &Image);
   void center(Icomplex_f &Image);
   void uncenter(Icomplex_f &Image);
   void uncenter(Icomplex_d &Image);

   inline Bool is_square_pow2(int Nl, int Nc)
   {
     if ((is_power_of_2(Nl) == True) && (is_power_of_2(Nc) == True) && (Nl==Nc)) return True;
     else return False;
   }

   public:
     Bool CenterZeroFreq; // if True, the zero frequency is in the middle
                          // of the image, else at the bottom left

     FFTN_2D(){CenterZeroFreq=True;} // Constructor

     void fftn2d(Ifloat &Image, Icomplex_f &Buff, Bool Reverse=False, bool normalize=false);
     void fftn2d(Ifloat &Image, Icomplex_d &Im_out, Bool Reverse=False);
     // 2D FFT transform of a real image

     inline void ifftn2d(Ifloat &Image, Icomplex_f &Buff, bool normalize=false)
                                          { fftn2d(Image, Buff, True, normalize);}
     
     void fftn2d (Icomplex_f &Buff, Bool Reverse=False, bool normalize=false);
     void fftn2d(Icomplex_d &Im_out, Bool Reverse=False);

     // 2D FFT transform of a complex image
     inline void ifftn2d (Icomplex_f &Buff, bool normalize=false)
			 								{ fftn2d(Buff, True, normalize); }
     
     void convolve(Ifloat &Ima1, Ifloat &Ima2, Ifloat &Result);
     // Result = convolution of Ima1 and Ima2
     // Result = ITF( TF(Ima1) * TF(Ima2))
     
     inline void convolve(Ifloat &Ima1, Ifloat &Ima2) {convolve(Ima1,Ima2,Ima1);}
	   
     void convolve(Ifloat &Ima1, Icomplex_f &Ima2, Ifloat &Result);
     // Result = convolution of Ima1 and Ima2
     // Result = ITF( TF(Ima1) * Ima2)
     inline void convolve(Ifloat &Ima1, Icomplex_f &Ima2) {convolve(Ima1,Ima2,Ima1);}

     void convolve_conj(Ifloat &Ima1, Icomplex_f &Ima2, Ifloat &Result);
     // Result = convolution of Ima1 and Ima2
     // Result = ITF( TF(Ima1) * conjugate(Ima2))
     inline void convolve_conj(Ifloat &Ima1, Icomplex_f &Ima2) {convolve_conj(Ima1,Ima2,Ima1);}
	
     void swap_buff(Icomplex_f &Buff, Bool Reverse=False);
     void swap_buff(Icomplex_d &Ima, Bool Reverse=False);

     // Swap the complex image (i.e. change the position of the zero 
     // frequency)

     void optim(int Nl, int Nc) {CenterZeroFreq = is_square_pow2(Nl,Nc);}
     // optimize the FFT parameters for a given image size

     ~FFTN_2D() {} // deallocation
};


void get_isotropic_spectrum(Ifloat & Data, fltarray & Spectrum, float Resol=1., float Dens=1.);
// compute the FFT of Data: FDATA
// compute the power spectrum of FDATA: PSDATA
// Average the power spectrum along circle with radius r and centered 
// at the zero frequency.
// r varies between 0 and Data.nl()*sqrt(2), which corresponds to 
// fequencies  between 0 and 0.5*sqrt(2). 
// The number of used radius is 
//                           Nr = (int) (Data.nl()*0.5*sqrt(2) / Resol + 0.5)
// The number of estimated points along a circle of radius r is 
//                           Nc = (int) (Dens*2*PI*r+0.5)
// Data = IN: input image of size Nl x Nc
// Spectrum = OUT: Output power spectrum = PS(0:Nr-1,3)
//                 Spectrum(*,0) = r coordinates
//                 Spectrum(*,1) = Averaged values along the circle of radius r
//                                 using Nc points = (int) (Dens*2*PI*r+0.5)
//                 Spectrum(*,2) = Spectrum(*,1)*Spectrum(*,0)*(Spectrum(*,0)+1)
#endif

