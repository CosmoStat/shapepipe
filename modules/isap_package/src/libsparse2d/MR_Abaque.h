/******************************************************************************
**                   Copyright (C) 1994 by CEA
*******************************************************************************
**
**    UNIT
**
**    Version: 14:46:31
**
**    Author: Jean-Luc Starck
**
**    Date:  96/06/13
**    
**    File:  MR_Abaque.h
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


#ifndef __ABAQUE__
#define __ABAQUE__

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
 
#include "Mr_FewEvent.h"

#define NMALLOC(x,n) (x * )malloc((unsigned)(n * sizeof(x)))
#define MALLOC(x) (x *)malloc((unsigned)(sizeof(x)))
#define CALLOC(x,n) (x *)calloc((unsigned)(n*sizeof(x)))

#ifndef MAX
#define	MAX(a,b)	(((a)>(b))?(a):(b))
#endif

#ifndef MIN
#define	MIN(a,b)	(((a)<(b))?(a):(b))
#endif

#define	CUBE_FABS(x)	(fabs(x) * fabs(x) * fabs(x))
#define CUBE(x)		((x) * (x) * (x))

#define MAX_CHAR_NOM	80

// #define DEFAULT_EPSILON 1e-3
#define MIN_EPSILON 1e-6
#define MAX_EPSILON 1.0
// #define SIGMA_BSPLINE_WAVELET 0.040717
#define BIN 512
#define HISTO_SAMPLING_POINTS 2048
#define DEFAULT_NBR_AUTOCONV 25
#define MIN_NBR_AUTOCONV 1
#define MAX_NBR_AUTOCONV 30

#define  Name_Abaque_Default "Abaque.fits"

#define  Name_Bspline           "Aba_bspline"
#define  Name_Wavelet           "Aba_wavelet"
#define  Name_Autoconv_Histo    "Aba_histo"
#define  Name_Distribution      "Aba_distrib"
#define  Name_Log_Distribution  "Aba_log_distrib"
#define  Name_Histo_Mean        "Aba_mean"
#define  Name_Histo_Sigma       "Aba_sigma"

// class StatEventPoisson {
//  int Number_of_AutoConv;
//  Bool InitOK; 
// public:
//  Ifloat Histo_Sigma;  // Histo_Sigma(p) = sigma of the autoconvolved histogram
//  Ifloat Histo_Mean;   // Histo_Mean(p) = mean of the autoconvolved histogram
//  Ifloat Histo_Bound;  // Histo_Bound(p, 0) = minimum reduced coefficient
//                       // Histo_Bound(p, 1) = maximum reduced coefficient
//  Ifloat Histo_Imag;   // Histo_Imag(3*p, 0:Np) = autoconvolved histogram
//                       // Histo_Imag(3*p+1, 0:Np) = Reduced Values
//                       // Histo_Imag(3*p+2, 0:Np) = Normalized Probability
//  Ifloat Histo_Distribution; 
//                   // Histo_Distribution(3*p, 0:Np) = Probability Reduced Values
//                   // Histo_Distribution(3*p+1, 0:Np) = Reduced Values
//                   // Histo_Distribution(3*p+2, 0:Np) = Repartition function
// 
//  Ifloat Threshold;
//                   // Threshold(p, 0) = Negative threshold
//                   // Threshold(p, 1) = Positive threshold
//  StatEventPoisson(int NautoConv=DEFAULT_NBR_AUTOCONV) {
//      Histo_Sigma.alloc(NautoConv+1,1, "Histo_Sigma");
//      Histo_Mean.alloc(NautoConv+1,1, "Histo_Mean");
//      Histo_Bound.alloc(NautoConv+1,2, "Histo_Bound");
//      Histo_Imag.alloc(3*(NautoConv+1),HISTO_SAMPLING_POINTS, "Histo_Imag");
//      Histo_Distribution.alloc(3*(NautoConv+1),
//                              HISTO_SAMPLING_POINTS, "Histo_Distribution");
//      Threshold.alloc(NautoConv+1,2, "Threshold");
//      Number_of_AutoConv = NautoConv;
//      InitOK = False;
//   }
//  void compute_distrib(Bool WriteAllFiles=False);
//                   // Computes the histogram autoconvolutions, and
//                   // the repartition function.
//                   // Histo_Sigma,Histo_Mean,Histo_Bound,
//                   // Histo_Imag, Histo_Distribution are set by this routine
//  void find_threshold(float Epsilon);
//                   // set the Threshold array for a confidence interval
//                   // given by Epsilon
//  float prob(float ValRed, int NEvent);
//                   // Compute the probability to have a reduce wavelet
//                   // coefficient ValRed with NEvent
//  float repartition(float ValRed, int NEvent);
//                   //Compute the integrated probability to have a reduced wavelet
//                   // lower than ValRed with NEvent (Repartition Function)
//  float a_trou_prob(float Coef, int NEvent, int Scale);
//                   // Compute the probability to have a wavelet coefficient
//                   // obtained bu the a trous algorithm at scale Scale
//                   // with NEvent
//  float a_trou_repartition(float Coef, int NEvent, int Scale);
//               //Compute the integrated probability to have a wavelet coefficient
//               // obtained bu the a trous algorithm at scale Scale
//               // lower than Coef with NEvent (Repartition Function)
// ~StatEventPoisson(){Number_of_AutoConv=0;InitOK=False;}
// };

// void bspline_histo_2D(Ifloat &histo, Ifloat &Histo_Bound, Bool WriteAllFiles);
// // Computes the one dimensional bspline  
// // Computes the two dimensional wavelet  
// // Computes the histogram of the wavelet in histo(0,*)  
//  
// void histo_convolution(int N, Ifloat &Imag, Ifloat &Histo_Bound);
// //   performs 2^(N-1) auto-convolution of the 1D-signal Signal1D. 
// //   Actually, it performs only the power-of-2 auto-convolutions. 
// //    But we keep the same number of points, so that means that we loose 
// //   in precision at each auto-convolution.
// //   To achieve this, you can use 2 differents methods : 
// //     "1" : computation is made in the original space ;
// //     "2" : computation is made in the Fourier space.
// //   In INPUT, Imag(0,*) is the original histogram;
// //             N is number of 2-power-auto-convolution ;
// //   In OUTPUT, Imag(3p,*) are the auto-convolued histograms of rank 2^(N-1) ;
// // 	      Imag is an image whose line 3p is the auto-convolued signal 
// //                     of rank 2^p.
// 
// 
// void histo_normalisation(Ifloat &Imag,Ifloat &distrib, Ifloat &Histo_Mean, Ifloat &Histo_Sigma,Ifloat &Histo_Bound,int N);
// // Normalizes the histogram 
// 
// void histo_distribution(Ifloat &distrib ,const Ifloat &Imag,Ifloat &Histo_Bound, Ifloat &Histo_Sigma, int N);
// // Computes the repartition function 
// 
// void histo_threshold(Ifloat &Threshold, const Ifloat &distrib,int N, float Epsilon);
// // computes the threshold max and min of the distribution function
// // to have a probability equal to 1 and 0 with a precision EPSILON. As we only
// // have discrete values for the distribution function, we use a linear
// // interpolation.
// 
// void distrib_log_transf(Ifloat &distrib,int N);
// // Apply the log transformation to the repartition function
#endif
