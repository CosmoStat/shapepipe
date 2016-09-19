

#ifndef _MR1DFEWEVENT_H_
#define _MR1DFEWEVENT_H_


#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "GlobalInc.h"
#include "IM1D_IO.h"
#include "IM_IO.h"

#include <vector>

#define DEF_NBR_AUTOCONV 25
#define SIGMA_BSPLINE_WAVELET 0.040717
#define DEFAULT_EPSILON 1e-3

class cReductHisto {

   int _NbNewHisto;
   int _MaxNbRedHisto;

public:
   fltarray _UsedHistoConv; // number of the two histo used for
                           // compute histogram j
   fltarray _RedHisto;      // Histogram h^3, h^5, h^6 and h^7
   fltarray _RedBound;      // Histogram bounds h^3, h^5, h^6 and h^7
   fltarray _RedBin;        // Histogram bin and nbpoint h^3, h^5, h^6 and h^7
   intarray _TabInd;        // ind of interp histo 
   
   fltarray _HistoConv;       // Histogram information on autoconv  
   
    
public:
   cReductHisto ();
   void setRedHisto (int IndHisto, fltarray& pRedHisto, float Min,
                     float Max, float Bin, int NbPoint);
   void compute_distribution (Bool WriteAllInfo=False);
   
private:   
   void histo_convolution (Bool WriteAllInfo);
};
   
   
enum typedim {E_1D, E_2D, E_3D};






class FewEventPoisson {

private:
   int  _NbAutoConv;         // number max of autoconv used
   bool _InitOk;             // initialisation flag
   
public:
   dblarray _HistoSigma;      // sigma of the autoconvolved histogram
   dblarray _HistoMean;       // mean of the autoconvolved histogram
   dblarray _HistoBound;      // min and max reduced coefficient
   dblarray _HistoBin;        // bin and number of points  
   dblarray _HistoConv;       // Histogram information on autoconv
                             // HistoConv(3*p,0:Np):autoconvolved histogram
			     // HistoConv(3*p+1,0:Np):Reduced Values
			     // HistoConv(3*p+2,0:Np):Normalized Probability
   dblarray _HistoDistrib;    // HistoDistrib(3*p,0:Np):Probability Reduced Values
                             // HistoDistrib(3*p+1,0:Np):Reduced Values
			     // HistoDistrib(3*p+2,0:Np):Repartition function
   Ifloat _Threshold;      // Threshold(p, 0) = Min value threshold
                             // Threshold(p, 1) = Max value threshold
   
   cReductHisto* _RedHisto;   // interpolate histogram
   Bool         _InterpHisto;// idem 
   
   
   Bool _OutForTest;         // file out for test purpose
  
 
public:

   // init attribute
   //
   FewEventPoisson (int NbAutoConv=DEF_NBR_AUTOCONV, Bool OutForTest=False);
   
  
   // create V0 projection ( of DataIn in V0 space)
   // Event_DataIn memorize event in DataIn before projection
   //
   //void building_V0_proj (fltarray &DataIn, intarray &Event_DataIn,
   //                       typedim Dim=E_2D);

   // Computes the histogram autoconvolutions, and the repartition function.
   // Histo_Sigma,Histo_Mean,Histo_Bound,
   // Histo_Imag, Histo_Distribution are set by this routin
   //
   void compute_distribution (Bool WriteAllInfo=False, typedim Dim=E_2D);

   // set the Threshold array for a confidence interval
   // given by Epsilon  
   void find_threshold (double Epsilon, Bool WriteAllInfo=False);
   
   // use interpolate hiso
   void set_interp_histo();
   
   
   
 float prob(float ValRed, int NEvent);
                  // Compute the probability to have a reduce wavelet
                  // coefficient ValRed with NEvent
 float repartition(float ValRed, int NEvent);
                  //Compute the integrated probability to have a reduced wavelet
                  // lower than ValRed with NEvent (Repartition Function)
 float a_trou_prob(float Coef, int NEvent, int Scale);
                  // Compute the probability to have a wavelet coefficient
                  // obtained bu the a trous algorithm at scale Scale
                  // with NEvent
 float a_trou_repartition(float Coef, int NEvent, int Scale);
              //Compute the integrated probability to have a wavelet coefficient
              // obtained bu the a trous algorithm at scale Scale
              // lower than Coef with NEvent (Repartition Function)
 
   
   
   
   
private:

   // Computes the one dimensional bspline  
   // Computes the histogram of the wavelet in histo(0,*)
   //
   void bspline_histo_1D (Bool WriteAllInfo=False);
   
   // Computes the one dimensional bspline  
   // Computes the histogram of the wavelet in histo(0,*)
   //
   void bspline_histo_2D (Bool WriteAllInfo=False);   
   
   // Computes the one dimensional bspline  
   // Computes the histogram of the wavelet in histo(0,*)
   //
   void bspline_histo_3D (Bool WriteAllInfo=False);  
   
   // create V0 projection ( of DataIn in V0 space) 1d
   // Event_DataIn memorize event in DataIn before projection   
   //void building2d_V0_proj (fltarray &DataIn, intarray &Event_DataIn);
     
   // create V0 projection ( of DataIn in V0 space) 2d
   // Event_DataIn memorize event in DataIn before projection   
   //void building1d_V0_proj (fltarray &DataIn, intarray &Event_DataIn);
   
   // performs 2^(N-1) auto-convolution of the 1D-signal Signal1D. 
   // Actually, it performs only the power-of-2 auto-convolutions. 
   // In INPUT, Imag(0,*) is the original histogram;
   //           N is number of 2-power-auto-convolution ;
   // In OUTPUT, Imag(3p,*) are the auto-convolued histograms of rank 2^(N-1) ;
   // 	         Imag is an image whose line 3p is the auto-convolued signal 
   //            of rank 2^p.
   void histo_convolution (Bool WriteAllInfo=False);
   
   // Normalizes the histogram 
   void histo_normalisation (Bool WriteAllInfo=False);
   void new_histo_normalisation (Bool WriteAllInfo=False);
   
   // Computes the repartition function 
   void histo_distribution (Bool WriteAllInfo=False);
   
   void shape_signal (fltarray &Signal1d);
   
   void show_param (char* text, float nbconv, float min, float max, 
                    float bin, float nbechant, dblarray& Histo);
public: //provioire 
   // set the Threshold array for a confidence interval
   // given by Epsilon  
   void histo_threshold (double Epsilon, Bool WriteAllInfo=False);
   
   float event_prob(float ValRed, int NEventReal, dblarray& Tab);
};


#endif

