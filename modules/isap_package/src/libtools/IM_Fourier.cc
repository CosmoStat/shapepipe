/*******************************************************************************
**
**    UNIT
**
**    Version: 3.1
**
**    Author: Jean-Luc Starck
**
**    Date:  96/05/02 
**    
**    File:  IM_Fourier.cc
**
*******************************************************************************
**
**    DESCRIPTION  
**    -----------  
**
*******************************************************************************
**
** Ifloat real(const Icomplex_f &Image)
**
** returns the real part of a complexe image
**
*******************************************************************************
**
** Ifloat imag(const Icomplex_f &Image)
**
** returns the imaginary part of a complexe image
**
*******************************************************************************
**
** Ifloat spectrum(const Icomplex_f &Image)
**
** returns the power spectrum  of a complexe image
**
*******************************************************************************
**
** const Icomplex_f conjugate(const Icomplex_f &Image)
**
** returns the conjugate of a complexe image
**
*******************************************************************************
**
** const Icomplex_f fft(const Icomplex_f &Image)
**
** returns the fft of a complexe image
**
*******************************************************************************
**
** const Icomplex_f fft(const Iint &Image)
**
** returns the fft of a int image
**
*******************************************************************************
**
** const Icomplex_f fft(const Ifloat &Image)
**
** returns the fft of a float image
**
*******************************************************************************
**
** const Icomplex_f invfft(const Icomplex_f &Image)
**
** returns the inverse fft of a complexe image
**
*******************************************************************************
**
** const Icomplex_f invfft(const Iint &Image)
**
** returns the inverse fft of a int image
**
*******************************************************************************
**
** const Icomplex_f invfft(const Ifloat &Image)
**
** returns the inverse fft of a float image
**
*******************************************************************************
**
** const Iint conv_fft( const Iint &Im1, const Iint &Im2)
**
** returns the convolution product  Im1 * Im2
**
*******************************************************************************
**
** const Ifloat conv_fft( const Ifloat &Im1, const Ifloat &Im2)
**
** returns the convolution product  Im1 * Im2
**
*******************************************************************************
**
** const Icomplex_f conv_fft( const Icomplex_f &Im1, const Icomplex_f &Im2)
**
** returns the convolution product  Im1 * Im2
**
******************************************************************************/

// static char sccsid[] = "@(#)IM_Fourier.cc 3.1 96/05/02 CEA 1994 @(#)";

#include "IM_Obj.h"


void real (Ifloat& Real, const Icomplex_f &Image) 
{
    int Nl = Image.nl();
    int Nc = Image.nc();

    for (int i=0; i < Nl; i++) 
    for (int j=0; j < Nc; j++) 
         Real(i,j) = Image(i,j).real();

}
void imag (Ifloat& Imag, const Icomplex_f &Image)
{
    int Nl = Image.nl();
    int Nc = Image.nc();  
     
    for (int i=0; i < Nl; i++) 
    for (int j=0; j < Nc; j++) 
         Imag(i,j) = Image(i,j).imag();
}


void real (Ifloat& Real, const Icomplex_d &Image) 
{
    int Nl = Image.nl();
    int Nc = Image.nc();

    for (int i=0; i < Nl; i++) 
    for (int j=0; j < Nc; j++) 
         Real(i,j) = Image(i,j).real();

}
void imag (Ifloat& Imag, const Icomplex_d &Image)
{
    int Nl = Image.nl();
    int Nc = Image.nc();  
     
    for (int i=0; i < Nl; i++) 
    for (int j=0; j < Nc; j++) 
         Imag(i,j) = Image(i,j).imag();
}

// const Icomplex_f invfft(const Iint &Im1)
// {
//     Icomplex_f *Result = new Icomplex_f(Im1.nl(),Im1.nc(), "fft");
// 
//     fft2d (Im1, *Result, -1);
//     return (*Result);
// }
// const Icomplex_f invfft(const Ifloat &Im1)
// {
//     Icomplex_f *Result = new Icomplex_f(Im1.nl(),Im1.nc(), "fft");
// 
//     fft2d (Im1, *Result, -1);
//     return (*Result);
// }
// const Icomplex_f invfft(const Icomplex_f &Image)
// {
//     int Nl = Image.nl();
//     int Nc = Image.nc();
//     Icomplex_f *Result = new Icomplex_f(Nl,Nc, "fft");
// 
//     fft2d (Image, *Result, -1);
//     return (*Result);
// }
// 
//  
