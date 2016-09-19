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
**    Date:  29/03/00 
**    
**    File:  WT1D_FFT.cc
**
*******************************************************************************
**
**    DESCRIPTION   routines concerting filters used by the 1D wavelet 
**    -----------   transform which need the FFT   
**
**
******************************************************************************/ 

 
#ifndef _WT1D_FFT_H_
#define _WT1D_FFT_H_

#include <math.h>
#include <stdio.h>
#include <string.h> 

#include "IM_Obj.h"
#include "IM_IO.h"
#include "FFTN_1D.h"

#define WT_PAVE_FFT 1
#define WT_PYR_FFT_DIFF_RESOL 2
#define WT_PYR_FFT_DIFF_SQUARE 3

enum type_wt_filter {WT_PHI,WT_PSI,WT_H,WT_H_TILDE,WT_G,WT_G_TILDE};

/****************************************************************************/

class WT1D_FFT {
   FFTN_1D FFT;
   
   int TypeWT;  // Type of wavelet transform 
                // type available are WT_PAVE_FFT,WT_PYR_FFT_DIFF_RESOL,
                // WT_PYR_FFT_DIFF_SQUARE

   cfarray Ima_h_cf; // buffer for high frequency component
   cfarray Ima_b_cf; // buffer for low frequency component
   cfarray Buff;     // buffer

   double scaling_function (float u, int N);
                        // return the scaling function value at frequency u
                        // N = size of the input data set

   double filter_wavelet (float u, int N);
                       // return the wavelet function value at frequency u

   double filter_h (float u, int N);
                       // return the h filter value at frequency u

   double filter_g (float u, int N);
                       // return the g filter value at frequency u

   double filter_h_tilde (float u, int N);
                       // return the ~h filter value at frequency u

   double filter_g_tilde (float u, int N);
                       // return the ~g filter value at frequency u


   void wave_cf_mult (cfarray &Ima_b, cfarray &Ima_h, int s, int N, Bool Down);
                       // one step transformation or reconstruction
                       // Down: input = if true, wavelet transform
                       //               otherwise wavelet reconstruction
                       // if ( Down == True)
                       //     Ima_b: input/output = input data and out low frequency
                       //     Ima_h: output =  high frequency
                       // else
                       //     Ima_b: input/output = input data and output 
                       //                           reconstructed signal
                       //     Ima_h: input =  high frequency
                       // N: input = data set size
                       // s: input: resolution level

   void under_sampling_2 (cfarray &Imag);
                       // undersample the signal. Signal size is divided by 2.

   void over_sampling_2 (cfarray &Imag, Bool Odd=False);
                       // oversample the signal. Signal size is multiply by 2.

   int pyr_tab_pos(int NbrScale, int N);
                       // initialize TabSize and TabPos for a WT transform 
                       // with NbrScale scales and for a data set with  N points
                       // return the number of pixels of the transform

   intarray TabPos;    // position table: TabPos(s) = starting position of
                       // band s in the transform signal

   intarray TabSize;  // size table: :  TabSize(s) = size of
                      // band s in the transform signal

   int PyrPixNbr;     // Number of pixels of the transform

   public:
    float filter (type_wt_filter Which_Filter, float u,int N);
    // return the value at frequency u of the filter 
    // specified by Which_Filter
    
    float Freq_Coup;  // Cut off frequency of the scaling function

    int pyr_np () { return  PyrPixNbr;} 
            // return the number of pixels of the transform

    int pyr_pos (int s) { return TabPos(s);}
            // return the position of the band s

    int pyr_size (int s) { return TabSize (s);}
            // return the size of the band s

    void create_filter (fltarray &Filter, type_wt_filter Which);
            // create a filter (phi,psi,h,g,~g, or ~h) 

    // Constructor
    WT1D_FFT() {TypeWT=WT_PAVE_FFT;Freq_Coup=0.5; PyrPixNbr=0;}
    WT1D_FFT(int TWT) {TypeWT=TWT;Freq_Coup=0.5; PyrPixNbr=0;}
    WT1D_FFT(int TWT, int NbrScale, int N) 
        {FFT.CenterZeroFreq=True;
	TypeWT=TWT;Freq_Coup=0.5; PyrPixNbr=pyr_tab_pos(NbrScale,N);}

    void transform (fltarray &Signal, fltarray &Trans, int Nbr_Plan);
    // WT transform of a signal with Nbr_Plan scales

    void transform_cf(cfarray & Data, fltarray &Trans, int Nbr_Plan);
    // WT transform of a signal with Nbr_Plan scales
    // input are output are in the Fourier space

    void recons (fltarray  &Trans, fltarray &Signal, int Nbr_plan);
    // WT reconstruction of a signal from its WT

    void recons_cf (fltarray  &Trans, cfarray & Data , int Nbr_plan);
    // WT reconstruction of a signal from its WT.
    // input are output are in the Fourier space

};

/****************************************************************************/

#endif
