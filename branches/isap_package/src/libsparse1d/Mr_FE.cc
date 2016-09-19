
#include "Mr_FE.h"
#include <iomanip>
//#include "IM_Noise.h"
//#include "MR1D_NoiseModel.h"
//#include "MR_NoiseModel.h"
//#include "Mr_FewEvent2d.h"


// for wavelet psi
#define FE_BIN1D 12384
#define FE_BIN1D_2 2*FE_BIN1D+1
#define FE_BIN2D 1598
#define FE_BIN2D_2 2*FE_BIN2D+1

// number max of autoconvolution
//#define DEF_NBR_AUTOCONV 25
// length max of autoconvolution


//#define FE_MAXHISTONBPOINT  2048+1
#define FE_MAXHISTONBPOINT  16384+1


// length max of autoconvolution
#define FE_MAXREDHISTO 16

// length of first histo (construct without convolution)
#define FE_LENGTH_FIRST_HISTO 1024+1
//#define FE_LENGTH_FIRST_HISTO 8192+1




// min val for bin*histo(n,i) 
 double FE_MIN_VAL_INTEG_HISTO = 1e-40;

#define FE_Name_Bspline    "_bspline"
#define FE_Name_x_Bspline  "_xbspline"
#define FE_Name_Wavelet    "_wavelet"

#define FE_Name_Histo      "_histo"
#define FE_Name_RedHisto   "_redhisto"
#define FE_Name_FirstHisto "_hfirsthisto"
#define FE_Name_HistoBin   "_hbin"
#define FE_Name_HistoBound "_hbound"
#define FE_Name_HistoNbPts "_hnbpts"

#define FE_Name_RedWavCoef "_hredwav"
#define FE_Name_Pdf        "_hpdf"
#define FE_Name_RedPdf     "_hredpdf"
#define FE_Name_InvPdf     "_hinvpdf"

#define FE_Name_Threshold  "_hthreshold"
#define FE_Name_HistoMean  "_hmean"
#define FE_Name_HistoSigma "_hsigma"

#define	FE_FE_CUBE_FABS(x)	(fabs(x) * fabs(x) * fabs(x))
#define FE_CUBE(x)		((x) * (x) * (x))

//#define DEBUG_FEW 1
#undef DEBUG_FEW   


using namespace std;





 




//---------------------------------------------------------------------------------

FewEvent::FewEvent( int NbAutoConv, Bool OutForTest ) {
//---------------------------------------------------------------------------------

   mNbAutoConv = NbAutoConv;
   mInitOk = False;
   mWriteAllInfo = False;
   mThresholdCoef = True;
   mWriteThreshold = False;
   
   mIndexBeg = 0;
   mIndexEnd = 0;
   
   mHisto.alloc( mNbAutoConv+1, FE_MAXHISTONBPOINT );
   mRedHisto.alloc( mNbAutoConv+1, FE_MAXHISTONBPOINT );
   mBin.alloc( mNbAutoConv+1 );
   mNbPts.alloc( mNbAutoConv+1 );
   mBound.alloc( mNbAutoConv+1, 2 );

   mRedBin.alloc( mNbAutoConv+1 );
   mRedNbPts.alloc( mNbAutoConv+1 );
   mRedBound.alloc( mNbAutoConv+1, 2 );
   
   mThreshold.alloc( mNbAutoConv+1, 2 );
   mMean.alloc( mNbAutoConv+1 );
   mStd.alloc( mNbAutoConv+1 );
   mRedWavCoef.alloc( mNbAutoConv+1, FE_MAXHISTONBPOINT );
   
   mDf.alloc( mNbAutoConv+1, FE_MAXHISTONBPOINT );
   mInvDf.alloc( mNbAutoConv+1, FE_MAXHISTONBPOINT );
   //mRedPdf.alloc( mNbAutoConv+1, FE_MAXHISTONBPOINT );
   
   _OutForTest = OutForTest;
}


 
//---------------------------------------------------------------------------------
void 
FewEvent::compute_distribution( Bool WriteAllInfo, typedim Dim ) {
//---------------------------------------------------------------------------------

   WriteAllInfo = Bool( WriteAllInfo || mWriteAllInfo ); 

   if( mInitOk ) exit (-1);
   
   // compute Bspline Histogramm
   //
   if  (mWriteAllInfo) 
      cout << "Compute the histogram of the Wavelet ... " << endl;
   
   switch( Dim ) {
   case E_1D : bspline_histo_1D( WriteAllInfo ); break;
   case E_2D : bspline_histo_2D( WriteAllInfo ); break;
   case E_3D : exit(-1); break;
   }
       
   // computes the auto-convolued histogramms
   //
   if( WriteAllInfo ) cout << "Compute the autoconvolutions of the histogram ... " << endl;
   histo_convolution( WriteAllInfo );
   
   // computes the normalisation
   //
   if( WriteAllInfo ) cout << "Compute reduced wavelet coef ... " << endl;   
   histo_reduced_wavelet_coef( WriteAllInfo );  
   
   // writing debug result
   if( WriteAllInfo == True ) {
      io_write_ima_float(FE_Name_Threshold, mThreshold);
      fits_write_dblarr( FE_Name_HistoMean, mMean );
      fits_write_dblarr( FE_Name_HistoSigma, mStd );
   } 
   
   // computes the distribution function of histogramms at each leve
   //
   histo_distribution( WriteAllInfo );
    
   mInitOk = True;
}
 
 
 
 
//---------------------------------------------------------------------------------
void
FewEvent::normalize( int HistoNum ) {
//---------------------------------------------------------------------------------

   double sum = 0.;
   for( int i=0; i<mNbPts( HistoNum ); ++i )
      sum += mHisto( HistoNum, i ) * mBin( HistoNum );
   for( int i=0; i<mNbPts( HistoNum ); ++i )
      mHisto( HistoNum, i ) /= sum;
}

//---------------------------------------------------------------------------------
void 
FewEvent::show_param( char* Msg, int HistoNum ) {
//---------------------------------------------------------------------------------

   //cout << " ==> " << Msg << " (histo " << HistoNum+1 << ")" << endl;
   cout << " autoconv for " << pow( 2.0, HistoNum ) 
        << "  events " << endl;
   cout	<< "      -- xmin    :" << mBound( HistoNum, 0 ) 
        << ", -- xmax    :" <<  mBound( HistoNum, 1 )  
        << ", -- bin    :" << mBin( HistoNum ) 
	<< ",   (    NBP:" << mNbPts( HistoNum ) << ")" << endl;
   if( HistoNum > 0 ) {
      cout << " reduction size of new histo" << std::endl;
      cout << "      suppress " << mIndexBeg << " points (histo value < " 
           << FE_MIN_VAL_INTEG_HISTO << ") on lower border" << endl;
      cout << "      suppress " << mNbPts( HistoNum )-1 - mIndexEnd 
           << " points (hist value < " 
           << FE_MIN_VAL_INTEG_HISTO << ") on upper border" << endl;
      cout << "      new size is : " << mNbPts( HistoNum ) - mIndexBeg 
            -( mNbPts( HistoNum )-1 - mIndexEnd ) << std::endl;
      if( mBin( HistoNum ) != mRedBin( HistoNum ) ) {
         cout << " decimation step (if NBP > (FE_MAXHISTONBPOINT-1)/2 = " 
	      << ( FE_MAXHISTONBPOINT-1 ) / 2 << ")" << std::endl; 
         cout << "      bin is increased by 2, bin:" << mRedBin( HistoNum ) << endl;
	 cout << "      new NBP is computed, NBP:" <<  mRedNbPts( HistoNum ) << endl;
      }
      cout	<< "      -- red xmin:" << mRedBound( HistoNum, 0 ) 
           << ", -- red xmax:" <<  mRedBound( HistoNum, 1 )  
           << ", -- red bin:" << mRedBin( HistoNum ) 
	   << ",   (red NBP:" << mRedNbPts( HistoNum ) << ")" << endl;
   }
   cout << std::endl;

   int localSizeSignal = mNbPts( HistoNum );
   dblarray Courbe( localSizeSignal );
   for( int i=0; i<localSizeSignal; ++i ) Courbe(i)=mHisto( HistoNum, i );    
   char FileName[256];
   sprintf( FileName, "_Histo_%d.fits", HistoNum );
   fits_write_dblarr( FileName, Courbe );

}



//---------------------------------------------------------------------------------
void 
FewEvent::bspline_histo_1D( Bool WriteAllInfo ) {
//---------------------------------------------------------------------------------

   WriteAllInfo = Bool( WriteAllInfo || mWriteAllInfo ); 

   // init Phi cubic bspline
   // performs B(x) for x in [-2;+2] :
   // B(x) =1/12 (|x-2|^3 - 4 |x-1|^3 + 6 |x|^3 - 4 |x+1|^3 + |x+2| ^ 3)
   if( WriteAllInfo ) cout << "Compute Bspline..." << endl;  
   dblarray phi_bspline( FE_BIN1D_2 ); 
   int histoNum = 0;
   float x;
   for ( int i=-FE_BIN1D/2; i<=FE_BIN1D/2; ++i ) {
   
      x = 2 * ((float)i)/((float)FE_BIN1D/2);
      phi_bspline(i+FE_BIN1D) = (      FE_CUBE((ABS(x-2))) - 4 * FE_CUBE((ABS(x-1)))
			    + 6 * FE_CUBE((ABS(x)))   - 4 * FE_CUBE((ABS(x+1))) 
			    +     FE_CUBE((ABS(x+2))) ) / 12;
   }      
   if( WriteAllInfo ) fits_write_dblarr( FE_Name_Bspline, phi_bspline );
   
   // performs Psi wavelet 
   // psi(x) = B(x) - 1/4 B(x/2)
   if (WriteAllInfo) cout << "Compute Wavelet..." << endl;
   dblarray psi_wavelet(FE_BIN1D_2);
   double MaxPsi = 0;
   double MinPsi = 0;
   for( int i=-FE_BIN1D; i<=FE_BIN1D; ++i) {
   
      psi_wavelet(i+FE_BIN1D) =         phi_bspline (i+FE_BIN1D)
                           - .5 * phi_bspline (i/2+FE_BIN1D);       
      if (psi_wavelet(i+FE_BIN1D) >= MaxPsi) MaxPsi = psi_wavelet(i+FE_BIN1D);
      if (psi_wavelet(i+FE_BIN1D) <= MinPsi) MinPsi = psi_wavelet(i+FE_BIN1D);
   }   
   if( WriteAllInfo ) fits_write_dblarr( FE_Name_Wavelet , psi_wavelet );      
   
   if( WriteAllInfo ) {
      std::cout << "Histogram : 0"   << endl;
      std::cout << "==============================" << endl;
   }     

   // Fullfills the histogramm (first length = FE_LENGTH_FIRST_HISTO points)
   // 
   if( WriteAllInfo ) cout << "Compute first histogram..." << endl;    
   int index;   
   int LocalSizeSignal = FE_LENGTH_FIRST_HISTO;
   for (int i=-FE_BIN1D;i<=FE_BIN1D;i++) {
   
      // we only use (HSP-1)/2 point in _HistoConv tab
      index = (int) (  (psi_wavelet(i+FE_BIN1D)-MinPsi)*(LocalSizeSignal-1)
                     / (MaxPsi-MinPsi));
      mHisto( histoNum, index )++;  // first histogram, no convolution
   } 
   for( int i=0; i<FE_LENGTH_FIRST_HISTO; ++i )
      mRedHisto( histoNum, i ) = mHisto( histoNum, i );
   mRedBin( histoNum ) = mBin( histoNum );
   mRedNbPts( histoNum ) = mNbPts( histoNum );
   
   // set param for first histogramm
   //
   mBin( histoNum ) = ( MaxPsi-MinPsi ) / ( LocalSizeSignal-1 );  
   mNbPts( histoNum ) = FE_LENGTH_FIRST_HISTO;
   mBound( histoNum, 0 ) = MinPsi;
   mBound( histoNum, 1 ) = MaxPsi;
   mRedBin( histoNum ) = mBin( histoNum );
   mRedNbPts( histoNum ) = mNbPts( histoNum );

   // normalize
   //
   normalize( histoNum );
   
   // debug
   //
   if (WriteAllInfo) show_param( "first histo", histoNum );
    		  
}


//---------------------------------------------------------------------------------
void 
FewEvent::bspline_histo_2D( Bool WriteAllInfo ) {
//---------------------------------------------------------------------------------

   WriteAllInfo = Bool( WriteAllInfo || mWriteAllInfo ); 

   // init Phi cubic bspline
   // performs B(x) for x in [-2;+2] :
   // B(x) =1/12 (|x-2|^3 - 4 |x-1|^3 + 6 |x|^3 - 4 |x+1|^3 + |x+2| ^ 3)
   if( WriteAllInfo ) cout << "Compute Bspline..." << endl;  
   dblarray phi_bspline( FE_BIN2D_2 );
   dblarray x_phi_bspline( FE_BIN2D_2 );
   int histoNum = 0;
   double x;
   for (int i=-FE_BIN2D/2;i<=FE_BIN2D/2;i++) {
   
      x = 2 * ((double)i)/((double)FE_BIN2D/2);
      x_phi_bspline(i+FE_BIN2D) = x;
      phi_bspline(i+FE_BIN2D) = (      FE_CUBE((ABS(x-2))) - 4 * FE_CUBE((ABS(x-1)))
			      + 6 * FE_CUBE((ABS(x)))   - 4 * FE_CUBE((ABS(x+1))) 
			      +     FE_CUBE((ABS(x+2))) ) / 12;
   }      
   if( WriteAllInfo ) {
      fits_write_dblarr( FE_Name_Bspline, phi_bspline );
      fits_write_dblarr( FE_Name_x_Bspline, x_phi_bspline );
   }
   
   // performs Psi wavelet 
   // psi(x,y) = B(x,y) - 1/4 B(x/2,y/2)
   // psi(x,y) = B(x) B(y) - 1/4 B(x/2) B(y/2)
   if( WriteAllInfo ) cout << "Compute Wavelet..." << endl;
   dblarray psi_wavelet( FE_BIN2D_2, FE_BIN2D_2 );
   double MaxPsi = 0;
   double MinPsi = 0;
   for (int i=-FE_BIN2D;i<=FE_BIN2D;i++) {
      for (int j=-FE_BIN2D;j<= FE_BIN2D;j++) {
	 psi_wavelet(i+FE_BIN2D, j+FE_BIN2D) = phi_bspline(i+FE_BIN2D) * phi_bspline(j+FE_BIN2D)
	                      - 0.25 * phi_bspline(i/2+FE_BIN2D) * phi_bspline(j/2+FE_BIN2D);	      
			           
         if (psi_wavelet(i+FE_BIN2D, j+FE_BIN2D) >= MaxPsi) 
	    MaxPsi = psi_wavelet(i+FE_BIN2D, j+FE_BIN2D);
         if (psi_wavelet(i+FE_BIN2D, j+FE_BIN2D) <= MinPsi) 
	    MinPsi = psi_wavelet(i+FE_BIN2D, j+FE_BIN2D);
      }  
   }
   if( WriteAllInfo ) fits_write_dblarr( FE_Name_Wavelet, psi_wavelet );      

   if( WriteAllInfo ) {
      std::cout << "Histogram : 1"   << endl;
      std::cout << "==============================" << endl;
   }     

   // Fullfills the histogramm (first length = FE_LENGTH_FIRST_HISTO points)
   // 
   int index;   
   int LocalSizeSignal = FE_LENGTH_FIRST_HISTO;
   for (int i=-FE_BIN2D;i<=FE_BIN2D;i++) {
      for (int j=-FE_BIN2D;j<=FE_BIN2D;j++) {
   
         // we only use (HSP-1)/2 point in _HistoConv tab
         index = (int) (  (psi_wavelet(i+FE_BIN2D, j+FE_BIN2D)-MinPsi)*(LocalSizeSignal-1)
                        / (MaxPsi-MinPsi));
         mHisto( histoNum, index)++;  // first histogram, no convolution
      } 
   }
   for( int i=0; i<FE_LENGTH_FIRST_HISTO; ++i )
      mRedHisto( histoNum, i ) = mHisto( histoNum, i );
   mRedBin( histoNum ) = mBin( histoNum );
   mRedNbPts( histoNum ) = mNbPts( histoNum );
   
   // set param for first histogramm
   //
   mBin( histoNum ) = ( MaxPsi-MinPsi ) / ( LocalSizeSignal-1 );  
   mNbPts( histoNum ) = FE_LENGTH_FIRST_HISTO;
   mBound( histoNum, 0 ) = MinPsi;
   mBound( histoNum, 1 ) = MaxPsi;
   mRedBin( histoNum ) = mBin( histoNum );
   mRedNbPts( histoNum ) = mNbPts( histoNum );
   mRedBound( histoNum, 0 ) = mBound( histoNum, 0 );
   mRedBound( histoNum, 1 ) = mBound( histoNum, 1 );
  
   // normalize
   //
   normalize( histoNum );
   
   // debug
   //
   if( WriteAllInfo ) show_param( "first histo", histoNum );
      		  
} 


//---------------------------------------------------------------------------------
void 
FewEvent::histo_convolution( Bool WriteAllInfo ) {
//---------------------------------------------------------------------------------

   WriteAllInfo = Bool( WriteAllInfo || mWriteAllInfo ); 
      
   // for all autoconv 1::mNbAutoConv (index 1 to mNbAutoConv+1 in _Histo...)
   //------------------------------------------------------------------------
   for( int nbConv=1; nbConv<mNbAutoConv+1; nbConv++ ) {


      if (WriteAllInfo) {
         cout << "Histogram : "  << nbConv+1 << endl;
	 cout << "==============================" << endl;
      }     
 
      // legth of new histogram : 2 * length of last histo
      int histoNbPts =  2 * ( mRedNbPts( nbConv-1 ) - 1 ) + 1; 
      
      // compute histogram convolution
      // -----------------------------
      // max = HistoNbPoint,  
      // inverse Signal [0->(max-1)]= signal [(max-1)->0]
      // for i in [0, (max-1)/2] do
      //    convol signal in [0, i] with inverse signal in [max-i, max]
      // then for i in [(max-1)/2+1, max-1] do
      //    convol signal in [i-(max-1)/2, (max-1)/2] with inverse 
      //    signal in [0, max-1-i]
      //
      for (long i=0; i<=(histoNbPts-1); i++) mHisto( nbConv, i ) = 0.;
      for (long i=0; i<=(histoNbPts-1)/2; i++) {
         
	 for (long j=0; j<=i; j++)
	    mHisto( nbConv, i ) +=   mRedHisto( nbConv-1, j ) 
	                           * mRedHisto( nbConv-1, i-j );
      }
         
      long jbeg = 1;
      long jend = (histoNbPts-1)/2;  
      for (long i=(histoNbPts-1)/2+1; i<=(histoNbPts-1); i++) {
      
     	 for (long j=jbeg; j<=jend; j++)   	     
	    mHisto( nbConv, i ) +=   mRedHisto( nbConv-1, j ) 
	                           * mRedHisto( nbConv-1, jend+jbeg-j );
	 jbeg ++;      
      }

      // current histogram parameters
      // new min = 2*min
      // new max = 2*max
      // nb pts = HistoNbPoint
      // bin = (new max - mew min)/nb pts-1
      //     = 2*(max-min) / (2*(HistoNbPoint-1)+1 -1)
      //     = (max-min)/HistoNbPoint-1 = old bin
      //
      mBin ( nbConv )     = mRedBin( nbConv-1 );
      mNbPts( nbConv )    = histoNbPts;
      mBound( nbConv, 0 ) = 2 * mRedBound( nbConv-1, 0 );
      mBound( nbConv, 1 ) = 2 * mRedBound( nbConv-1, 1 );
      
      // normalisation for density prob
      // 
      normalize( nbConv );
      
      // threshold pts under FE_MIN_VAL_INTEG_HISTO	
      //
      mIndexBeg = 0;
      mIndexEnd = histoNbPts-1;
      if( mThresholdCoef ) {
         
         while( mHisto( nbConv, mIndexBeg ) < FE_MIN_VAL_INTEG_HISTO ) { mIndexBeg++; }
         while( mHisto( nbConv, mIndexEnd ) < FE_MIN_VAL_INTEG_HISTO ) { mIndexEnd--; }
      
         histoNbPts = mIndexEnd - mIndexBeg + 1;
         if( histoNbPts/2*2 == histoNbPts) {
	    histoNbPts++;
	    if( mIndexBeg > 0 ) mIndexBeg--;
	    else if( mIndexEnd < (unsigned int)histoNbPts-1 ) mIndexEnd++;
	    else 
	       std::cout << "!!!! there's a little pb in the resuduction step...!!!"
	                 << std::endl;
         }
         mRedBound( nbConv, 0 ) = mBound( nbConv, 0 ) + mBin ( nbConv ) * mIndexBeg;
         mRedBound( nbConv, 1 ) = mBound( nbConv, 0 ) + mBin ( nbConv ) * mIndexEnd;
         
      } else {
         
         mRedBound( nbConv, 0 ) = mBound( nbConv, 0 );
         mRedBound( nbConv, 1 ) = mBound( nbConv, 0 );
      
      }
      

      // decimation?
      //
      if (histoNbPts > ( FE_MAXHISTONBPOINT-1 ) / 2 ) { 
      
         int halfSize = ( histoNbPts - 1 ) / 2 + 1;
	 for( int i=0; i<halfSize-1 ; ++i )
	    mRedHisto( nbConv, i ) =   mHisto( nbConv, 2*i+mIndexBeg ) 
	                             + mHisto( nbConv, 2*i+1+mIndexBeg );
         mRedHisto( nbConv, halfSize-1 ) = mHisto( nbConv, histoNbPts - 1+mIndexBeg );
   
	 mRedBin( nbConv ) = mBin( nbConv ) * 2;
	 mRedNbPts( nbConv ) = halfSize;
      
      } else {
      
         for( int i=0; i<histoNbPts; ++i )
            mRedHisto( nbConv, i ) = mHisto( nbConv, i+mIndexBeg );
	 mRedBin( nbConv ) = mBin( nbConv );
	 mRedNbPts( nbConv ) = histoNbPts;
	     
      }
  
      if( WriteAllInfo ) show_param( "", nbConv );
      
   }   
   
   if( WriteAllInfo == True ) {
      fits_write_dblarr( FE_Name_Histo, mHisto );
      fits_write_dblarr( FE_Name_RedHisto, mRedHisto );
      fits_write_dblarr( FE_Name_HistoBound, mBound );  
      fits_write_dblarr( FE_Name_HistoBin, mBin );
      dblarray prov( mNbPts.nx() );
      for( int i=0; i<mNbPts.nx(); ++i )
         prov(i) = mNbPts(i);
      fits_write_dblarr( FE_Name_HistoNbPts, prov );
      
   }
}



//---------------------------------------------------------------------------------
void 
FewEvent::histo_reduced_wavelet_coef( Bool WriteAllInfo ) {
//---------------------------------------------------------------------------------

   WriteAllInfo = Bool( WriteAllInfo || mWriteAllInfo ); 
      
   // for all autoconv 1::mNbAutoConv (index 1 to mNbAutoConv+1 in _Histo...)
   //------------------------------------------------------------------------
   for( int nbConv=0; nbConv<mNbAutoConv+1; nbConv++ ) {
      
      // local var
      // ---------
      double mean  = 0.0;
      double std   = 0.0;
      double bin   = mBin( nbConv );
      double nbPts = mNbPts( nbConv );
      double w_min = mBound( nbConv, 0 );
      
      // compute mean
      // ------------
      for( int j=0; j<nbPts; j++ ) {
         double wav_coef = ((double) j) * bin + w_min;
	 double prov = mHisto( nbConv, j ) * bin;
         mean += prov * wav_coef;
      }   
      mMean( nbConv ) = mean;
      
      for( int j=0; j<nbPts; j++ ) {     
         double wav_coef = ((double) j) * bin + w_min;
	 double prov = mHisto( nbConv, j ) * bin;
         std  += prov * ( wav_coef - mean ) * ( wav_coef - mean );
      }   
      std = sqrt( std );

      mStd( nbConv ) = std;
      
       
      if( WriteAllInfo ) {
     
         cout << "Histogram : "  << nbConv+1 << endl;
	 cout << "==============================" << endl;
         
         //cout << " ==> (histo " << nbConv+1 << ")" << endl;
         cout	<< "      -- xmin:" << mBound( nbConv, 0 ) 
                << ", -- xmax:" <<  mBound( nbConv, 1 )  
                << ", -- bin:" << mBin( nbConv ) 
	        << ",  (NBP:" << mNbPts( nbConv ) << ")" << std::endl;
         cout << " compute wavelet coef " << std::endl;
	 cout << "        -- mean:" << mMean( nbConv )
	        << ", -- sigma:" << mStd( nbConv ) << endl;
      }
      
      double redmean = 0.;
      double redstd = 0.;
      for( int j=0; j<mNbPts( nbConv ); j++ ) {
        
	 /*--- compute real values of the wavelet coefficients ---*/
         mRedWavCoef( nbConv, j ) = 
                  ( ((float) j) * bin + w_min - mean ) / std; 
         double prov = mHisto( nbConv, j ) * bin;
         redmean += prov * mRedWavCoef( nbConv, j );
	 redstd += prov * mRedWavCoef( nbConv, j ) * mRedWavCoef( nbConv, j );
	  
      }
    
      if( WriteAllInfo ) {
      
         cout << " compute reduced wavelet coef " << std::endl;
	 cout << "        -- min red w:" << mRedWavCoef( nbConv, 0 )
	        << ", -- max red w:" << mRedWavCoef( nbConv, mNbPts( nbConv )-1 )
		<< ", -- mean red w:" << redmean 
		<< ", -- sigma red w:" << sqrt( redstd - redmean*redmean )
		<< endl;
	 cout << std::endl;
     
      
      }
   }
   if( WriteAllInfo )       
      fits_write_dblarr( FE_Name_RedWavCoef, mRedWavCoef );

}



//---------------------------------------------------------------------------------
void 
FewEvent::histo_distribution( Bool WriteAllInfo ) {
//---------------------------------------------------------------------------------

   WriteAllInfo = Bool( WriteAllInfo || mWriteAllInfo ); 
      
   // for all histo 
   //------------------------------------------------------------------------
   for( int nbConv=0; nbConv<mNbAutoConv+1; nbConv++ ) {
         
      if( WriteAllInfo ) {   
        cout << "Histogram : "  << nbConv+1 << endl;
        cout << "==============================" << endl;
      }
         
      double bin   = mBin( nbConv );
      int nbPts    = mNbPts( nbConv );
      //double sigma = mStd( nbConv );
      
      //Probability Reduced Values
      //
      // init first values
      mDf( nbConv, 0 ) = mHisto( nbConv, 0 ) * bin;
      //mRedPdf( nbConv, 0 ) = mHisto( nbConv, 0 ) * bin / sigma;
     
      // integ other val (Dist(i) = Dist(i-1) + val(i)*bin)
      for( int nbp=1; nbp<nbPts; nbp++ ) {
         
	 mDf( nbConv, nbp ) =   mDf( nbConv, nbp-1 ) 
	                       + mHisto( nbConv, nbp ) * bin;

         //mRedPdf ( nbConv, nbp ) =   mDf( nbConv, nbp-1 ) 
	 //                          + mHisto( nbConv, nbp ) * bin / sigma;
      }
      
      // inv pdf
      mInvDf( nbConv, nbPts-1 ) = mHisto( nbConv, nbPts-1 ) * bin;
      for( int nbp=1; nbp<nbPts; nbp++ ) {
	 mInvDf( nbConv, nbPts-1-nbp ) =   mInvDf( nbConv, nbPts-1-nbp+1 ) 
	                       + mHisto( nbConv, nbPts-1-nbp ) * bin;
      }
      
       
      
      // normalisation (last value must be one)
      //
      double normInvDf = mInvDf( nbConv, 0 );
      for( int nbp=0; nbp<nbPts; nbp++ ) {
          mDf( nbConv, nbp ) = mDf( nbConv, nbp ) / mDf( nbConv, nbPts-1 );
          //mRedPdf( nbConv, nbp ) =   mRedPdf( nbConv, nbp ) 
	  //                         / mRedPdf( nbConv, nbPts-1 );
          mInvDf( nbConv, nbp ) = mInvDf( nbConv, nbp ) / normInvDf;
      }

      if( WriteAllInfo ) {
         //cout << " ==> (histo " << nbConv+1 << ")" << endl;
         cout << "     -- min val pdf:" << mDf( nbConv, 0 )
              << ", -- max val pdf:" <<  mDf( nbConv, nbPts-1 )
	      << std::endl;
         cout << "     -- min val inv pdf:" << mInvDf( nbConv, nbPts-1 )
              << ", -- max val inv pdf:" <<  mInvDf( nbConv, 0 )
	      << std::endl;
	 cout << std::endl;
      }

      //if( WriteAllInfo ) show_param( "all info", nbConv );
   }
   
   if( WriteAllInfo ) {
        fits_write_dblarr( FE_Name_Pdf, mDf );
        //fits_write_dblarr( FE_Name_RedPdf, mRedPdf );
        fits_write_dblarr( FE_Name_InvPdf, mInvDf );
   }
}


//---------------------------------------------------------------------------------
void 
FewEvent::find_threshold ( double Epsilon, Bool WriteAllInfo ) {
//---------------------------------------------------------------------------------

   WriteAllInfo = Bool( WriteAllInfo || mWriteAllInfo ); 
   
   if( mInitOk == False ) compute_distribution( WriteAllInfo );
   histo_threshold( Epsilon, WriteAllInfo );
}



//---------------------------------------------------------------------------------
void 
FewEvent::histo_threshold ( double Epsilon, Bool WriteAllInfo ) {
//---------------------------------------------------------------------------------

   WriteAllInfo = Bool( WriteAllInfo || mWriteAllInfo ); 

   if( WriteAllInfo ) {
      cout << "FewEvent::histo_threshold" << endl;
      cout << "  epsilon:" << Epsilon << endl;
   }
   
   // for all histo 
   //------------------------------------------------------------------------
   for( int nbConv=0; nbConv<mNbAutoConv+1; nbConv++ ) {
   
      if( WriteAllInfo ) {   
        cout << "Histogram : "  << nbConv+1 << endl;
        cout << "==============================" << endl;
      }     
     
      // search indice of max threshold (RepFunc >= 1-eps) 
      // 
      long indexMax = (long)mNbPts( nbConv ) - 1;
      while(     mInvDf( nbConv, indexMax ) < Epsilon 
              && indexMax > 0 ) { indexMax--; }
	      
      if (WriteAllInfo) {
         std::cout << " compute threshold at the end of the pdf" << std::endl;
         std::cout << "    indexMax for pdf[indexMax->End] = epsilon : " 
	           <<  indexMax <<  std::endl;         
                                 
      }
    
      if( indexMax == mNbPts( nbConv ) - 1 ) {
       
         mThreshold( nbConv, 1 ) = mRedWavCoef( nbConv, indexMax );
// PERHAPS MUST BE ALWAYS TRACED ON SCREEN
         if (WriteAllInfo) {
	    std::cout << "    !! max threshold for histo " << nbConv
                      << " is on last value of histogram !!" << std::endl;
	    std::cout << "    !! epsilon=" << Epsilon 
	              << " is lower than the last integated value !!" << std::endl;
         }
      
      } else {
      
         // inter between IndexMax and IndexMax-1
         double x1 = mInvDf( nbConv, indexMax );
         double y1 = mRedWavCoef( nbConv, indexMax );
         double x2 = mInvDf( nbConv, indexMax+1 );
         double y2 = mRedWavCoef( nbConv, indexMax+1 );
         double a = ( y1 - y2 ) / ( x1 - x2 );
         double b = y1 - a * x1;
         mThreshold ( nbConv, 1 ) = a * ( Epsilon ) + b;

         if (WriteAllInfo) {
            std::cout << "    SeuilMax:=[" << mRedWavCoef( nbConv, indexMax ) << "," 
	              << mRedWavCoef( nbConv, indexMax+1 ) << "]"
	              << ", repartMax:=[" << mInvDf( nbConv, indexMax ) << "," 
	              << mInvDf( nbConv, indexMax+1 ) << "]" << std::endl;
	    std::cout << "    -- " << mInvDf( nbConv, indexMax ) 
	              << "/" << (Epsilon) << "/" 
		      << mInvDf( nbConv, indexMax+1 ) << endl;
         }
     
      }
     
     
     
     
     
     
     
     
      // search indice of min threshold (RepFunc <= eps)
      //
      long indexMin = 0;
      while(     mDf( nbConv, indexMin ) < Epsilon 
              && indexMin < (long)mNbPts( nbConv ) - 1 ) { indexMin++; }
      
      if (WriteAllInfo) {
         std::cout << " compute threshold at the begining of the pdf" << std::endl;
         std::cout << "    indexMin for pdf[0->indexMin] = epsilon " 
	           <<  indexMin << std::endl;         
       }
              
      if( indexMin == 0 ) {
      
         mThreshold( nbConv, 0) = mRedWavCoef( nbConv, indexMin );
// PERHAPS MUST BE ALWAYS TRACED ON SCREEN
         if (WriteAllInfo) {
	    std::cout << "    !! min threshold for histo " << nbConv
                      << " is on first value of histogram !! " << std::endl;
	    std::cout << "    !! epsilon=" << Epsilon 
	              << " is lower than the first integated value !!" << std::endl;
         }
      
      } else {

         // inter between IndexMin and IndexMin+1
         double x1 = mDf( nbConv, indexMin-1 );
         double y1 = mRedWavCoef( nbConv, indexMin-1 );
         double x2 = mDf( nbConv, indexMin );
         double y2 = mRedWavCoef( nbConv, indexMin );
         double a = ( y1 - y2 ) / ( x1 - x2 );
         double b = y1 - a * x1;
         mThreshold ( nbConv, 0)  = a * ( Epsilon ) + b;

         if (WriteAllInfo) {
            std::cout << "    SeuilMin:=[" << mRedWavCoef( nbConv, indexMin-1 ) << "," 
	              << mRedWavCoef( nbConv, indexMin ) << "]"
	              << ", repartMin:=[" << mDf( nbConv, indexMin-1 ) << "," 
	              << mDf( nbConv, indexMin ) << "]" << std::endl;;
            std::cout << "    -- " << mDf( nbConv, indexMin-1 ) << "/" 
	              << (Epsilon) << "/" << mDf( nbConv, indexMin ) 
		      << std::endl;
         }
         
      }
      
    
      
      // trace
      //
      if( WriteAllInfo ) {
         std::cout << " interpolated threshold for epsilon = "
	           << Epsilon << std::endl;
         std::cout << "    ==< " 
                   << " SeuilMin:=" << mThreshold( nbConv, 0 )
                   //<< ", SeuilMax:=" << mThreshold( nbConv, 2 )
                   << ", SeuilMax:=" << mThreshold( nbConv, 1 )
                   << " >==" << std::endl; 
         std::cout << std::endl;
      }   
        
   }
   
   if( WriteAllInfo || mWriteThreshold ) 
      io_write_ima_float( FE_Name_Threshold, mThreshold );
}










 

 
//---------------------------------------------------------------------------------
float 
FewEvent::event_prob( float WavRed, int NEventReal, TypeFunc typeFunc, 
                      Bool ProbSignal ) {
//---------------------------------------------------------------------------------

   if( mWriteAllInfo ) std::cout << "FewEvent::event_prob" << std::endl;
   // compute "power" : 2^(power-1) < NEventReal <= 2^(power)
   //
   int nEvent=1;
   int powerAfter=0;
   int powerBefore=0;
   while( nEvent < NEventReal ) { powerAfter ++; nEvent *= 2; }
   
   // test if the event position must be interpolated
   Bool interpEvent = False;
   if( powerAfter > DEF_NBR_AUTOCONV ) powerAfter=DEF_NBR_AUTOCONV;
   else { 
      if( ( powerAfter > 0 ) && ( nEvent != NEventReal ) ) {
         powerBefore = powerAfter -1;
         interpEvent = True;
      }
   }   

   if( mWriteAllInfo ) {
      std::cout << "    NEventReal : " << NEventReal;
      if( powerAfter == powerBefore ) {
         std::cout << " in [" << pow( 2., powerBefore ) << "]  " 
		   << "(2^" << powerBefore <<  ")" << std::endl;
      
      } else {
         std::cout << " in [" << pow( 2., powerBefore ) << "," 
                   << pow( 2., powerAfter ) << "]  " 
		   << "(2^" << powerBefore << ",2^" << powerAfter << ")" << std::endl;
      }
   }
      
   // Computes the probability related to a reduced wavelet
   // coefficients WavRed, and NEvent events
   //
   // search IndCoef so that  RepFunc (PowerAfter, IndCoef) <  WavRed
   double probAft = get_prob( WavRed, powerAfter, typeFunc, ProbSignal );
    
          
   // if no interpolation, return the "true value" ProbAfter
   if ( !interpEvent ) {
      if( mWriteAllInfo ) {
         std::cout << "    !! No Interp !! reduced coef : " << WavRed << ", event number : " 
                   << NEventReal << ",  prob = " << probAft << std::endl;
      }
      return ( probAft );
   } 
    
   // intrepolate between PowerAfter and PowerBefore
   //    
   // Computes the probability related to a reduced wavelet
   // coefficients WavRed, and n_event/2  events
   double probBef = get_prob( WavRed, powerBefore, typeFunc, ProbSignal );
   
   if( mWriteAllInfo ) {
         std::cout << "    !!    Interp !! reduced coef : " << WavRed 
                   << ", prob_bef:" << probBef
                   << ", prob_aft:" << probAft 
                   << ", n_event:" << nEvent 
                   << ", n_event_real:" << NEventReal << std::endl;
   }
   
   // interpolation
   //
   //double weight = 2.* ( nEvent - NEventReal ) / (float)nEvent;
   //float prob = weight * probBef + ( 1. - weight ) * probAft;   
   double weight = powerAfter - log((double)NEventReal) / log(2.);
   float prob = exp(weight*log(probBef) + (1.-weight)*log(probAft));
   
   if( mWriteAllInfo )
      std::cout << "    fewEvent::prob reduced coef : " << WavRed << ", event number : " 
                << NEventReal << ",  prob = " << prob << std::endl;
    
   return ( prob ); 
 
}






//---------------------------------------------------------------------------------
float 
FewEvent::get_prob( float RedWavCoef, int HistoNumber, TypeFunc typeFunc, 
                    Bool ProbSignal ) {
//---------------------------------------------------------------------------------
       
       
   if( mWriteAllInfo ) 
      std::cout << "search prob that red wav coef (=" << RedWavCoef
                << ") is signifiant" << std::endl;
      
   // Computes the probability related to a reduced wavelet
   // coefficients ValRed, and NEvent events
   //
   // search IndCoef so that  RepFunc (PowerAfter, IndCoef) <  ValRed
   int indCoefUp = 0;
   long nbPts = (long)mNbPts( HistoNumber );
   while (    ( indCoefUp < nbPts-1 ) 
           && ( mRedWavCoef( HistoNumber, indCoefUp) < RedWavCoef ) ) {
      indCoefUp++; 
   }
   int indCoefDown = indCoefUp;
    
   // compute prob RepFunc <  ValRed for PowerAfter  
   //   
   double posProb=0., negProb=0., prob=0.;
   if( indCoefUp==0 ) {
   
      if( typeFunc == E_Pdf ) {
         prob = mHisto( HistoNumber, indCoefUp ) / mStd( HistoNumber );
      }
      if( typeFunc == E_Df ) {
         negProb = mDf( HistoNumber, indCoefUp );
	 posProb = mInvDf( HistoNumber, nbPts - 1 );
	 
	 if( mWriteAllInfo ) {
            std::cout << "    get_prob : indice = " << 0 
		      << ", total indice number = " << nbPts << std::endl;
            std::cout << setprecision( 16 ) << "    get_prob : prob up in ["  
                      << mDf( HistoNumber, 0 ) << "] prob down in ["
		      << 1-mInvDf( HistoNumber, 0 ) << "]" 
	              << setprecision( 6 ) << std::endl;
	 }
      }
      
   } else if( indCoefUp==nbPts ) {
   
      if( typeFunc == E_Pdf ) {
         prob = mHisto( HistoNumber, indCoefUp-1 ) / mStd( HistoNumber );
      }
      if( typeFunc == E_Df ) {
         negProb = mDf( HistoNumber, indCoefUp-1 );
	 posProb = mInvDf( HistoNumber, 0 );
	 
	 if( mWriteAllInfo ) {
            std::cout << "    get_prob : indice = " << nbPts - 1
		      << ", total indice number = " << nbPts << std::endl;
            std::cout << setprecision( 16 ) << "    get_prob : prob up in ["  
                      << 1-mDf( HistoNumber, nbPts - 1 ) << "] prob down in ["
		      << mInvDf( HistoNumber, nbPts - 1 ) << "]" 
	              << setprecision( 6 ) << std::endl;
	 }
      }
      
   } else {
   
      double redWavCoefBef = mRedWavCoef( HistoNumber, indCoefUp-1 );
      double redWavCoefAft = mRedWavCoef( HistoNumber, indCoefUp );
      double valBef=0., valAft=0., weight=0.;
      double posValBef=0., posValAft=0., posWeight=0.;
      double negValBef=0., negValAft=0., negWeight=0.;
      
      if( typeFunc == E_Pdf ) {
         valAft = mHisto( HistoNumber, indCoefUp ) / mStd( HistoNumber );
         valBef = mHisto( HistoNumber, indCoefUp-1 ) / mStd( HistoNumber );
         if ( ABS( valAft - valBef ) > FLOAT_EPSILON ) 
            weight = 1. - ( RedWavCoef - redWavCoefBef ) / ( redWavCoefAft - redWavCoefBef );
         else 
            weight = 0.5;
         prob = weight * valBef + ( 1. - weight ) * valAft;
      }
      if( typeFunc == E_Df ) {
         negValAft = mDf( HistoNumber, indCoefUp );
         negValBef = mDf( HistoNumber, indCoefUp-1 );
         if ( ABS( negValAft - negValBef ) > FLOAT_EPSILON ) 
            negWeight = 1. - ( RedWavCoef - redWavCoefBef ) / ( redWavCoefAft - redWavCoefBef );
         else 
            negWeight = 0.5;
         negProb = negWeight * negValBef + ( 1. - negWeight ) * negValAft;
	 
         posValAft = mInvDf( HistoNumber, indCoefDown+1 );
         posValBef = mInvDf( HistoNumber, indCoefDown );
         if ( ABS( posValAft - posValBef ) > FLOAT_EPSILON ) 
            posWeight = 1. - ( RedWavCoef - redWavCoefBef ) / ( redWavCoefAft - redWavCoefBef );
         else 
            posWeight = 0.5;
         posProb = posWeight * posValBef + ( 1. - posWeight ) * posValAft;
	 
         if( mWriteAllInfo ) {
            std::cout << "   get_prob : red wav coef = " << RedWavCoef << std::endl;
            std::cout << "   get_prob : red wav coef in [" << redWavCoefBef 
	              << "," << redWavCoefAft << "]" <<  std::endl;
            std::cout << "   get_prob : indice up = [" << indCoefUp-1
	              << "," << indCoefUp << "]" << std::endl;
	    std::cout << "   get_prob : indice down = [" << indCoefDown << ","
		      << indCoefDown+1 << "]"<< std::endl;
		      
            std::cout << setprecision( 16 );
            if( RedWavCoef < 0 ) {
               std::cout << "   get_prob : [*] prob up in ["  
		         << mDf( HistoNumber, indCoefUp-1 ) << ","
                         << mDf( HistoNumber, indCoefUp ) << "]" << std::endl;
	       std::cout << "   get_prob : prob down in ["
		         << 1-mInvDf( HistoNumber, indCoefDown+1 ) << ","
		         << 1-mInvDf( HistoNumber, indCoefDown ) << "]" << std::endl;
            }
            if( RedWavCoef > 0 ) {
               std::cout << "   get_prob : prob up in ["  
		         << 1-mDf( HistoNumber, indCoefUp-1 ) << ","
                         << 1-mDf( HistoNumber, indCoefUp ) << "]" << std::endl;
	       std::cout << "   get_prob : [*] prob down in ["
		         << mInvDf( HistoNumber, indCoefDown+1 ) << ","
		         << mInvDf( HistoNumber, indCoefDown ) << "]" << std::endl;
            }
	    std::cout << setprecision( 6 );
	    std::cout << "   get_prob : neg weighted prob " <<  negProb
		      << ", computed in [" << negValBef << "," << negValAft << "]" 
		      << " with weight:" << negWeight  << std::endl;
	    std::cout << "   get_prob : pos weighted prob " <<  posProb
		      << ", computed in [" << posValBef << "," << posValAft << "]" 
		      << " with weight:" << posWeight  << std::endl;
	 }
	 
      }
      
      if( !ProbSignal ) {
         prob = negValAft;
      } else {
         if( RedWavCoef < 0 )   prob = negProb;
         if( RedWavCoef >= 0 )  prob = posProb;
      }
                  
      if( mWriteAllInfo ) {	 
	 std::cout << "   get_prob : weighted prob " << prob  << std::endl;
      }
      
   }     

   return prob;

}





//---------------------------------------------------------------------------------
float 
FewEvent::prob( float ValRed, int NEvent ) {
//---------------------------------------------------------------------------------
   if( mWriteAllInfo ) std::cout << "FewEvent::prob" << std::endl;
   return event_prob( ValRed, NEvent, E_Pdf );
}

//--------------------------------------------------------------------------------- 
float 
FewEvent::repartition( float ValRed, int NEvent, Bool ProbSignal ) {
//---------------------------------------------------------------------------------
   if( mWriteAllInfo ) std::cout << "FewEvent::repartition" << std::endl;
   return event_prob( ValRed, NEvent, E_Df, ProbSignal );
}
 
//---------------------------------------------------------------------------------
float 
FewEvent::a_trou_prob( float Coef, int NEvent, int Scale ) {
//---------------------------------------------------------------------------------
    float RedCoef,alpha,sc;
    if( NEvent != 0 ) {
       for( sc=0,alpha=1.0; sc<Scale; sc++) alpha *= 4.;
       RedCoef = Coef * alpha /  sqrt((float) NEvent) / SIGMA_BSPLINE_WAVELET;
       return  prob( RedCoef, NEvent );
    }
    return 0.;
}
 
//---------------------------------------------------------------------------------
float 
FewEvent::a_trou_repartition(float Coef, int NEvent, int Scale, Bool ProbSignal ) {
//---------------------------------------------------------------------------------
    float RedCoef,alpha,sc;
    if( NEvent != 0 ) {
       for (sc=0,alpha=1.0; sc< Scale; sc++) alpha *= 4.;
       RedCoef = Coef * alpha /  sqrt((float) NEvent) / SIGMA_BSPLINE_WAVELET;
       if( mWriteAllInfo ) 
          std::cout << "!!!!!!!!!!!!!!! Scale:" << Scale << ", NEvent:" << NEvent 
                    << ", w:" << Coef << ", red w:" 
                    << Coef * alpha /  sqrt((float) NEvent) / SIGMA_BSPLINE_WAVELET
                    << std::endl;
        return  repartition(RedCoef, NEvent, ProbSignal);
    }
    return 0.;
}
 
 
//---------------------------------------------------------------------------------
void 
FewEvent::write_threshold( Bool Flag ) {                
//---------------------------------------------------------------------------------
   mWriteThreshold = Flag;
}

 
//---------------------------------------------------------------------------------
void 
FewEvent::set_write_all( Bool Flag ) {                
//---------------------------------------------------------------------------------
   mWriteAllInfo = Flag;
}


