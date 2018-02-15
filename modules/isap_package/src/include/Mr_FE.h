

#ifndef _MR_FE_H_
#define _MR_FE_H_


#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "GlobalInc.h"
#include "IM1D_IO.h"
#include "IM_IO.h"
#include "Mr_FewEvent.h"

#include <vector>

#define DEF_NBR_AUTOCONV 25
#define SIGMA_BSPLINE_WAVELET 0.040717
#define DEFAULT_EPSILON 1e-3

   



class FewEvent {

private: //--------------------------------------------------------------------

   int  mNbAutoConv;         // number max of autoconv used
   Bool mInitOk;             // initialisation flag
   Bool mWriteAllInfo;    
   Bool _OutForTest;         // file out for test purpose
   Bool mThresholdCoef;
   Bool mWriteThreshold;
   
   unsigned int mIndexBeg;
   unsigned int mIndexEnd;

   dblarray mHisto;
   dblarray mRedHisto;
   dblarray mBin;
   intarray mNbPts;
   dblarray mBound;
 
   dblarray mRedBin;
   intarray mRedNbPts;
   dblarray mRedBound;
   
   dblarray mMean;
   dblarray mStd;
   dblarray mRedWavCoef;
   
   dblarray mDf;
   dblarray mInvDf;
   //dblarray mRedPdf;

public : //--------------------------------------------------------------------
 
    Ifloat mThreshold;

 

public : //--------------------------------------------------------------------

   //-- ctor --
   FewEvent( int NbAutoConv=DEF_NBR_AUTOCONV, Bool OutForTest=False );
   
   //-- manipulator --
   void compute_distribution( Bool WriteAllInfo=False, typedim Dim=E_2D );
  
   enum  TypeFunc { E_Pdf, E_Df } ;        
   
private: //--------------------------------------------------------------------

   // normalisation
   void normalize( int histoNum );
   
   // trace
   void show_param( char* Msg, int histoNum );

   // Computes the one dimensional bspline  
   void bspline_histo_1D( Bool WriteAllInfo=False );
  
public:
   // Computes the two dimensional bspline  
   void bspline_histo_2D( Bool WriteAllInfo=False );   
   
   // performs 2^(N-1) auto-convolution of the 1D-signal Signal1D. 
   // Actually, it performs only the power-of-2 auto-convolutions. 
   void histo_convolution( Bool WriteAllInfo=False );

   // reduced wavelet coef
   void histo_reduced_wavelet_coef( Bool WriteAllInfo=False );

   // Computes the repartition function 
   void histo_distribution (Bool WriteAllInfo=False);
   
   // get prob for event_prob
   float get_prob( float ValRed, int HistoNumber, TypeFunc typeFunc, 
                   Bool ProbSignal=False );

   // set write flag
   void set_write_all( Bool Flag ) ;


public:

   // set the Threshold array for a confidence interval
   // given by Epsilon  
   void find_threshold ( double Epsilon, Bool WriteAllInfo=False );
   void histo_threshold (double Epsilon, Bool WriteAllInfo=False);
   float event_prob( float ValRed, int NEventReal, TypeFunc typeFunc, 
                     Bool ProbSignal=False  );

   // Compute the probability to have a reduce wavelet
   // coefficient ValRed with NEvent  
   float prob( float ValRed, int NEvent );
                 
   // Compute the integrated probability to have a reduced wavelet
   // lower than ValRed with NEvent (Repartition Function)
   float repartition( float ValRed, int NEvent, Bool ProbSignal=False );
               
   // Compute the probability to have a wavelet coefficient
   // obtained bu the a trous algorithm at scale Scale
   // with NEvent
   float a_trou_prob( float Coef, int NEvent, int Scale );
                  
   //Compute the integrated probability to have a wavelet coefficient
   // obtained bu the a trous algorithm at scale Scale
   // lower than Coef with NEvent (Repartition Function)  
   float a_trou_repartition( float Coef, int NEvent, int Scale, 
                             Bool ProbSignal=False );
           
                  
   void write_threshold( Bool Flag );                
 
   
};



#endif

