/******************************************************************************
**                   Copyright (C) 2001 by CEA
*******************************************************************************
**
**    UNIT
**
**    Version: 1.0
**
**    Author: Jean-Luc Starck
**
**    Date:  20/03/01
**    
**    File:  FFTN_1D.h
**
*******************************************************************************
**
**    DESCRIPTION  1D FFT Class routine for any size of images
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

#ifndef _CFFTN_1D_H_
#define _CFFTN_1D_H_

#include "GlobalInc.h"
#include "FFTN.h"

class FFTN_1D: public FFTN {
   dblarray BuffW;  // buffer used by the DCT
   intarray BuffIp; // buffer used by the DCT

   inline void pix_swap(complex_f &a, complex_f &b) { complex_f temp=a;a=b;b=temp;}
   inline void pix_swap(double &a, double &b) { double temp=a;a=b;b=temp;}
   inline void pix_swap(complex_d &a, complex_d &b) { complex_d temp=a;a=b;b=temp;}

   public:
     Bool CenterZeroFreq; // if True, the zero frequency is in the middle
                          // of the image, else at the left

     FFTN_1D(){CenterZeroFreq=True;} // Constructor

     void fftn1d(fltarray &Signal, complex_f *Buff, Bool Reverse=False, bool normalize=false);
     // 1D FFT transform of a real signal

     void fftn1d(fltarray &SignalIn, cfarray &SignalOut, Bool Reverse=False, bool normalize=false);
     // 1D FFT transform of a real signal

     void fftn1d (complex_f *Buff, int N, Bool Reverse=False, bool normalize=false);
     // 1D FFT transform of a complex signal

     void fftn1d (cfarray &Signal, Bool Reverse=False, bool normalize=false);
     // 1D FFT transform of a complex signal
     
     void fftn1d (cfarray &SignalIn, cfarray &SignalOut, Bool Reverse=False, bool normalize=false);
     // 1D FFT transform of a complex signal

    void fftn1d (cdarray &Signal, Bool Reverse=False, bool normalize=false);
     // 1D FFT transform of a complex signal
     
     void fftn1d (cdarray &SignalIn, cdarray &SignalOut, Bool Reverse=False, bool normalize=false);
     // 1D FFT transform of a complex signal
          
     void fftn1d(dblarray &Signal, complex_d *Buff, Bool Reverse=False, bool normalize=false);
     // 1D FFT transform of a real signal (single precision)
     
     void fftn1d(fltarray &Signal, complex_d *Buff, Bool Reverse=False, bool normalize=false);
     // 1D FFT transform of a real signal (single precision)
     
     void fftn1d (complex_d *Buff, int N, Bool Reverse=False, bool normalize=false);
     // 1D FFT transform of a complex signal (double precision)
     
      void dct1d (fltarray &Signal, fltarray &SigOut, Bool Reverse=False);
     // 1D DCT transform of a  signal
     // signal size must be a power of 2
     void  dct1d(dblarray &Data,  Bool Reverse=False);
         
     void swap_buff(dblarray &Data,  Bool Reverse=False);
     void swap_buff(complex_f *Buff, int N, Bool Reverse=False);
     void swap_buff(complex_d *Buff, int N, Bool Reverse=False);
     // Swap the complex data (i.e. change the position of the zero 
     // frequency)
    void  center(fltarray &Dat);
    void  center(dblarray &Dat);
    void  uncenter(fltarray &Dat);
    void  uncenter(dblarray &Dat);
    void  center(complex_f *Dat, int N);
    void  uncenter(complex_f *Dat, int N);
    void  center(complex_d *Dat, int N);
    void  uncenter(complex_d *Dat, int N);
    void  convolve(fltarray &D1, fltarray &D2, fltarray &Result);

     ~FFTN_1D() {} // deallocation
};
void im1d_shift(fltarray &Data1, fltarray & Data2, int Dx);
void im1d_shift(dblarray &Data1, dblarray & Data2, int Dx);


#define NBR_STD_WIN 5
enum type_std_win {W_HAMMING, W_HANNING, W_GAUSS, W_BLACKMAN, W_RECT};
#define DEF_STD_WIN W_GAUSS



inline const char * StringSTDWin (type_std_win type)
{
    switch (type)
    {
    case W_HAMMING: 
      return ("Hamming window.");break;
    case W_HANNING: 
      return ("Hanning window.");break;
    case W_GAUSS: 
      return ("Gaussian window.");break;
    case W_BLACKMAN:
      return ("Blackman window.");break;
    case W_RECT:
      return ("Rectangular window.");break;
    }
    return ("Error: bad type of window.");
}
inline void std_win_usage()
{
    fprintf(OUTMAN, "         [-w type_of_window]\n");
    for (int i = 0; i < NBR_STD_WIN ; i++)
    fprintf(OUTMAN, "              %d: %s \n",i+1,
                                            StringSTDWin((type_std_win )i));
    fprintf(OUTMAN, "             default is Gaussian window.\n");
}

class STD_Window {
  public:
  void hamming(fltarray & Data, float W=0.5);
  void hanning(fltarray & Data, float W=0.5);
  void gaussian(fltarray & Data, float W=0.5);
  void blackman(fltarray & Data, float W=0.5);
  void rect(fltarray & Data, float W=0.5);
  void make_win(fltarray & Data, type_std_win type, float W=0.5);
};

#define NBR_SPEC 2
enum type_spec {SPEC_STF, SPEC_WIGNER_VILLE, SPEC_CHOI_WILLIAMS};
#define DEF_SPEC SPEC_STF
inline const char * StringSPEC (type_spec type)
{
    switch (type)
    {
    case SPEC_STF: 
      return ("Short Term Fourier Transform Spectogram.");break;
    case SPEC_WIGNER_VILLE: 
      return ("Wigner-Ville Distribution.");break;
    case SPEC_CHOI_WILLIAMS: 
      return ("Choi-Williams Distribution.");break;
    }
    return ("Error: bad type of window.");
}
inline void spec_usage()
{
    fprintf(OUTMAN, "         [-t TimeFrequency_Distrib]\n");
    for (int i = 0; i < NBR_SPEC ; i++)
    fprintf(OUTMAN, "              %d: %s \n",i+1,
                                            StringSPEC((type_spec )i));
    fprintf(OUTMAN, "             default is 1.\n");
}
#define DEF_SIZE_STD_WIN 512
#define DEF_PARAM_STD_WIN 0.5

class ST_FFTN: public FFTN_1D {
   void recons_direct(complex_f *Trans, fltarray & Data);
        // reconstruction without using the window
	// Data must be allocated to the correct size 
   void recons_win(complex_f *Trans, fltarray & Data);
        // reconstruction using the window
        // Data must be allocated to the correct size 
	
  type_std_win TypeWin; // Type of window to be used in the STF  
  int SizeWin;          // Size of the window. By default SizeWin = Signal Size / 2
  float WinParam;       // Window parameter
  fltarray Win;           // array which contains the window to used in the STF
  int Step;            // Step between two window position.
  int Nlt;             // Number of lines in the transform
  int Nct;             // Number of columns in the transform
  public:
  
  void alloc(int N, type_std_win type, float WinPar=DEF_PARAM_STD_WIN, 
                  int Size=DEF_SIZE_STD_WIN,int StepAna=1);
                  // Set the window to use in the STF
  int nl(int N) {return Nlt; } // return the number of lines in the 
                                   // short term Fourier transform
  int nc(int N) {return Nct; } // return the number of columnsin the 
                                  // short term Fourier transform

  ST_FFTN() {TypeWin = DEF_STD_WIN;SizeWin=0;CenterZeroFreq = False;Step=1;
             Nlt=Nct=1;}
  
  void transform(fltarray & Data, complex_f *Trans);
            // Calculate the STF of Data and store the result in Trans
	    // Trans must be allocated before the call
	    //   Trans = complex_f[N*N], where N is the data size

  void dct_transform(fltarray & Data, fltarray & Trans);
            // Calculate the ST-DCT of Data and store the result in Trans
	    // Trans must be allocated before the call
	    //   Trans = array[N*N], where N is the data size
  void dct_recons(fltarray & Trans, fltarray & Data);
	    	    
  void recons(complex_f *Trans, fltarray & Data, Bool UseWindow=False);
            // inverse short term Fourier transform
	    // if UseWindow == True, all coefficients are used for the
	    // reconstruction, and it may be more stable when the
	    // coefficient  have been modified, It is however not
	    // an exact reconstruction on the border.
	    // if UseWindow == False, it is an exact reconstruction
  void spectogram(complex_f *Buff_STF, int Nx, int Ny, fltarray &Spec);
  void spectogram(fltarray & Data, fltarray &Spec);

  void wigner_wille(fltarray & Data, fltarray & TabWV);
           // Compute the Wigner-Wille distribution
	   // Data(0..N): IN = input data
	   // TabWV(0..SizeWin,0..N/Step): OUT =  Wigner-Wille distribution
	   //                   SizeWin defined the number of frequencies

  void choi_williams(fltarray & Data, fltarray &TabCW,  
                   type_std_win TWinF=DEF_STD_WIN, 
		   int WindowF_Length= -1, 
	           double sigma=1., float WinParamF=DEF_PARAM_STD_WIN);
};

#endif

