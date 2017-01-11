

#include "Mr_FewEvent.h"
//#include "IM_Noise.h"
//#include "MR1D_NoiseModel.h"
//#include "MR_NoiseModel.h"
//#include "Mr_FewEvent2d.h"


// for wavelet psi
#define BIN1D 12384
#define BIN1D_2 2*BIN1D+1
#define BIN2D 1598
#define BIN2D_2 2*BIN2D+1

// number max of autoconvolution
//#define DEF_NBR_AUTOCONV 25
// length max of autoconvolution


//#define MAXHISTONBPOINT  2048+1
#define MAXHISTONBPOINT  16384+1


// length max of autoconvolution
#define MAXREDHISTO 16

// length of first histo (construct without convolution)
#define LENGTH_FIRST_HISTO 1024+1
//#define LENGTH_FIRST_HISTO 8192+1




// min val for bin*histo(n,i) 
 double MIN_VAL_INTEG_HISTO = 1e-20;

#define Name_Bspline_   "Aba_bspline"
#define Name_Wavelet_   "Aba_wavelet"
#define Name_HistoConv  "Aba_histo"

#define Name_HistoMean  "_histomean"
#define Name_HistoSigma "_histosigma"
#define Name_Threshold "_threshold"
#define Name_FirstHisto "_firsthisto"
#define Name_HistoBin "_histobin"
#define Name_HistoBound "_histobound"

#define	CUBE_FABS(x)	(fabs(x) * fabs(x) * fabs(x))
#define CUBE(x)		((x) * (x) * (x))

//#define DEBUG_FEW 1
#undef DEBUG_FEW   

// for intrep histo, number of hsito used in convol = f(indice)
#define MAXNBNEWHISTO 11
int Tab[MAXNBNEWHISTO]         = {3,5,6,7,9,10,11,12,13,14,15};
int FirstHisto[MAXNBNEWHISTO] =  {1,1,2,3,1, 2, 3, 4, 5, 6, 7};
int SecondHisto[MAXNBNEWHISTO] = {2,4,4,4,8, 8, 8, 8, 8, 8, 8};



FewEventPoisson::FewEventPoisson (int NbAutoConv, Bool OutForTest) {

   _NbAutoConv = NbAutoConv;
   _InitOk = False;
   _HistoConv.alloc (3*(_NbAutoConv+1), MAXHISTONBPOINT);
   _HistoConv.init(0.);
   _HistoBound.alloc (_NbAutoConv+1, 2);
   _HistoBin.alloc (_NbAutoConv+1, 2); 
   _HistoMean.alloc (_NbAutoConv+1);
   _HistoSigma.alloc (_NbAutoConv+1);
   _HistoDistrib.alloc (3*(_NbAutoConv+1), MAXHISTONBPOINT); 
   _Threshold.alloc(_NbAutoConv+1, 2, ""); 
   _InterpHisto = False;
   _RedHisto = (cReductHisto*)NULL;
   
   _OutForTest = OutForTest;
}


void FewEventPoisson::compute_distribution (Bool WriteAllInfo, 
                                              typedim Dim) {


   if (_InitOk) exit (-1);
   
   // compute Bspline Histogramm
   //
   if  (WriteAllInfo) 
      cout << "Compute the histogram of the Wavelet ... " << endl;
   
   switch (Dim) {
   case E_1D : bspline_histo_1D (WriteAllInfo); break;
   case E_2D : bspline_histo_2D (WriteAllInfo); break;
   case E_3D : exit(-1); break;
   }
       
   
   // computes the auto-convolued histogramms
   //
   if  (WriteAllInfo) cout << "Compute the autoconvolutions of the histogram ... " << endl;
   histo_convolution (WriteAllInfo);
   
   // computes the normalisation
   //
   if  (WriteAllInfo) cout << "Compute normalisation ... " << endl;   
   histo_normalisation (WriteAllInfo);  //MODIF
   //new_histo_normalisation (WriteAllInfo);
   
   // writing debug result
   if (WriteAllInfo == True) {
      fits_write_dblarr (Name_HistoConv, _HistoConv);
      fits_write_dblarr (Name_HistoMean, _HistoMean);
      fits_write_dblarr (Name_HistoSigma, _HistoSigma);
      fits_write_dblarr (Name_HistoBound, _HistoBound);      
      fits_write_dblarr (Name_HistoBin, _HistoBin);
   } 
   
   // computes the distribution function of histogramms at each leve
   //
   histo_distribution (WriteAllInfo);
    
   _InitOk = True;
}
 

/*void FewEventPoisson::building_V0_proj (fltarray &DataIn, 
                                        intarray &Event_DataIn,
       				        typedim  Dim) {					  
   switch (Dim) {
   case E_1D : building1d_V0_proj (DataIn, Event_DataIn); break;
   case E_2D : building2d_V0_proj (DataIn, Event_DataIn); break;
   case E_3D : exit(-1); break;
   }					  
					  
}


void FewEventPoisson::building1d_V0_proj (fltarray &DataIn, 
                                          intarray &Event_DataIn) {
   // var 
   float energy_x;
   
   // mem of Data with number of events
   int Nx = DataIn.nx();
   Event_DataIn.reform (Nx);
   for (int i=0;i<Nx;i++) Event_DataIn(i) = (int) (DataIn(i)+.5);
   
   // V0 Proj
   // -------
   DataIn.init();
   // for all point of DataIn
   for (int k=0;k<Nx;k++) {
      
      // for the number of events in DataIn at pos k
      for (int e=0; e<Event_DataIn(k); e++) {
      
        // projection in V0 (scalar product with function phi)
        for (int x_cur=-2;x_cur<=2;x_cur++) {
	   
 	    energy_x = (CUBE_FABS(x_cur- 2) + CUBE_FABS(x_cur + 2)
	          - 4 * (CUBE_FABS(x_cur - 1) + CUBE_FABS(x_cur + 1))
		  + 6 * CUBE_FABS(x_cur)) / 12.; 	      
		  
            if ((k+x_cur >= 0)    &&  (k+x_cur < Nx)) 
	        DataIn(k+x_cur) += energy_x;
         }
      }
   }
} 

   
void FewEventPoisson::building2d_V0_proj (fltarray &DataIn, 
                                          intarray &Event_DataIn) {
   // var 
   float energy_x;
   
   // mem of Data with number of events
   int Nx = DataIn.nx();
   Event_DataIn.reform (Nx);
   for (int i=0;i<Nx;i++) Event_DataIn(i) = (int) (DataIn(i)+.5);
   
   // V0 Proj
   // -------
   DataIn.init();
   // for all point of DataIn
   for (int k=0;k<Nx;k++) {
      
      // for the number of events in DataIn at pos k
      for (int e=0; e<Event_DataIn(k); e++) {
      
        // projection in V0 (scalar product with function phi)
        for (int x_cur=-2;x_cur<=2;x_cur++) {
	   
 	    energy_x = (CUBE_FABS(x_cur- 2) + CUBE_FABS(x_cur + 2)
	          - 4 * (CUBE_FABS(x_cur - 1) + CUBE_FABS(x_cur + 1))
		  + 6 * CUBE_FABS(x_cur)) / 12.; 	      
		  
            if ((k+x_cur >= 0)    &&  (k+x_cur < Nx)) 
	        DataIn(k+x_cur) += energy_x;
         }
      }
   }
} 
*/   


void FewEventPoisson::bspline_histo_1D (Bool WriteAllInfo) {

   if (WriteAllInfo) cout << "Compute Bspline..." << endl;  
   // init Phi cubic bspline
   // performs B(x) for x in [-2;+2] :
   // B(x) =1/12 (|x-2|^3 - 4 |x-1|^3 + 6 |x|^3 - 4 |x+1|^3 + |x+2| ^ 3)
   fltarray phi_bspline(BIN1D_2); 
   float x;
   for (int i=-BIN1D/2;i<=BIN1D/2;i++) {
   
      x = 2 * ((float)i)/((float)BIN1D/2);
      phi_bspline(i+BIN1D) = (      CUBE((ABS(x-2))) - 4 * CUBE((ABS(x-1)))
			    + 6 * CUBE((ABS(x)))   - 4 * CUBE((ABS(x+1))) 
			    +     CUBE((ABS(x+2))) ) / 12;
   }      
//   if (WriteAllInfo) fits_write_fltarr (Name_Bspline_, phi_bspline);
   
   if (WriteAllInfo) cout << "Compute Wavelet..." << endl;
   // performs Psi wavelet 
   // psi(x) = B(x) - 1/4 B(x/2)
   fltarray psi_wavelet(BIN1D_2);
   float MaxPsi = 0;
   float MinPsi = 0;
   for (int i=-BIN1D;i<=BIN1D;i++) {
   
      psi_wavelet(i+BIN1D) =         phi_bspline (i+BIN1D)
                           - .5 * phi_bspline (i/2+BIN1D);       
      if (psi_wavelet(i+BIN1D) >= MaxPsi) MaxPsi = psi_wavelet(i+BIN1D);
      if (psi_wavelet(i+BIN1D) <= MinPsi) MinPsi = psi_wavelet(i+BIN1D);
   }   
//   if (WriteAllInfo) fits_write_fltarr (Name_Wavelet_, psi_wavelet);      
   
   if (WriteAllInfo) cout << "Compute first histogram..." << endl;    
   // Fullfills the histogramm (first length = 1025 points)
   // 
   int index;   
   int LocalSizeSignal = LENGTH_FIRST_HISTO;
   for (int i=-BIN1D;i<=BIN1D;i++) {
   
      // we only use (HSP-1)/2 point in _HistoConv tab
      index = (int) (  (psi_wavelet(i+BIN1D)-MinPsi)*(LocalSizeSignal-1)
                     / (MaxPsi-MinPsi));
      _HistoConv(0,index)++;  // first histogram, no convolution
   } 
   
   // bin
   float LocalHistoBin   = (MaxPsi-MinPsi) / (LocalSizeSignal-1);         
   
   // set min, max, bin and nb points of first level (histogram only)
   //
   _HistoBound(0,0) = MinPsi;
   _HistoBound(0,1) = MaxPsi;
   _HistoBin (0,0) = LocalHistoBin;
   _HistoBin (0,1) = LocalSizeSignal;   
   
   // normalisation of the histogram
   // integral of histogram
//   float sum_check=0.0;
//   for(int i=0;i<LocalSizeSignal;i++)
//      sum_check += _HistoConv(0,i)*LocalHistoBin;
      
   // it must be a probability density -> normalisation by sum_check
//   for(int i=0;i<LocalSizeSignal;i++) _HistoConv(0,i) /= sum_check;

    if (WriteAllInfo) {
//#if DEBUG_FEW      
      dblarray Courbe (LocalSizeSignal);
      for (int i=0;i<LocalSizeSignal;i++) Courbe(i)=_HistoConv(0,i);    
       show_param ("End Convol Compute", 0,
                   _HistoBound(0,0), _HistoBound(0,1),
 		  LocalHistoBin, LocalSizeSignal, Courbe); 
//#endif	
      // char FileName[256];
      // sprintf (FileName, "Histo_%d.fits", 0);
      // fits_write_fltarr (FileName, Courbe);
    
   }
      
   // interp histo 
   if (_InterpHisto) {
      fltarray Courbe (LocalSizeSignal);
      for (int i=0;i<LocalSizeSignal;i++) Courbe(i)=_HistoConv(0,i);
      _RedHisto->setRedHisto (0, Courbe, _HistoBound(0,0), _HistoBound(0,1),
	                      _HistoBin (0,0), (int) _HistoBin (0,1)); 
   } 
   
/*_HistoConv.init(0.);   
fltarray Data;
fits_read_fltarr("Aba_firsthisto.fits", Data);  
for (int i=0;i<Data.n_elem();i++) _HistoConv(0,i) = Data(i);
_HistoBound(0,0) = -0.040418;
_HistoBound(0,1) = 1.456781;
_HistoBin (0,0) = 0.000731;
_HistoBin (0,1) = 511;*/    
    		  
}


void FewEventPoisson::bspline_histo_2D (Bool WriteAllInfo) {


   if (WriteAllInfo) cout << "Compute Bspline..." << endl;  
   // init Phi cubic bspline
   // performs B(x) for x in [-2;+2] :
   // B(x) =1/12 (|x-2|^3 - 4 |x-1|^3 + 6 |x|^3 - 4 |x+1|^3 + |x+2| ^ 3)
   dblarray phi_bspline(BIN2D_2);
   double x;
   for (int i=-BIN2D/2;i<=BIN2D/2;i++) {
   
      x = 2 * ((double)i)/((double)BIN2D/2);
      phi_bspline(i+BIN2D) = (      CUBE((ABS(x-2))) - 4 * CUBE((ABS(x-1)))
			      + 6 * CUBE((ABS(x)))   - 4 * CUBE((ABS(x+1))) 
			      +     CUBE((ABS(x+2))) ) / 12;
   }      
   if (WriteAllInfo) fits_write_dblarr (Name_Bspline_, phi_bspline);
   
   if (WriteAllInfo) cout << "Compute Wavelet..." << endl;
   // performs Psi wavelet 
   // psi(x,y) = B(x,y) - 1/4 B(x/2,y/2)
   // psi(x,y) = B(x) B(y) - 1/4 B(x/2) B(y/2)
   dblarray psi_wavelet(BIN2D_2,BIN2D_2);
   double MaxPsi = 0;
   double MinPsi = 0;
   for (int i=-BIN2D;i<=BIN2D;i++) {
      for (int j=-BIN2D;j<= BIN2D;j++) {
	 psi_wavelet(i+BIN2D, j+BIN2D) = phi_bspline(i+BIN2D) * phi_bspline(j+BIN2D)
	                      - 0.25 * phi_bspline(i/2+BIN2D) * phi_bspline(j/2+BIN2D);	      
			           
         if (psi_wavelet(i+BIN2D, j+BIN2D) >= MaxPsi) 
	    MaxPsi = psi_wavelet(i+BIN2D, j+BIN2D);
         if (psi_wavelet(i+BIN2D, j+BIN2D) <= MinPsi) 
	    MinPsi = psi_wavelet(i+BIN2D, j+BIN2D);
      }  
   }
   if (WriteAllInfo) fits_write_dblarr (Name_Wavelet_, psi_wavelet);      

//#if DEBUG_FEW 
      if (WriteAllInfo) {
         cout << "Histogram : 0"   << endl;
	 cout << "==============================" << endl;
      }     
//#endif       
   // Fullfills the histogramm (first length = 1025 points)
   // 
   int index;   
   int LocalSizeSignal = LENGTH_FIRST_HISTO;
   for (int i=-BIN2D;i<=BIN2D;i++) {
      for (int j=-BIN2D;j<=BIN2D;j++) {
   
         // we only use (HSP-1)/2 point in _HistoConv tab
         index = (int) (  (psi_wavelet(i+BIN2D, j+BIN2D)-MinPsi)*(LocalSizeSignal-1)
                        / (MaxPsi-MinPsi));
         _HistoConv(0,index)++;  // first histogram, no convolution
      } 
   }
   
   // bin
   double LocalHistoBin   = (MaxPsi-MinPsi) / (LocalSizeSignal-1);         
   
   // set min, max, bin and nb points of first level (histogram only)
   //
   _HistoBound(0,0) = MinPsi;
   _HistoBound(0,1) = MaxPsi;
   _HistoBin (0,0) = LocalHistoBin;
   _HistoBin (0,1) = LocalSizeSignal;   
   
   // normalisation of the histogram
   // integral of histogram
   double sum_check=0.0;
   for(int i=0;i<LocalSizeSignal;i++)
      sum_check += _HistoConv(0,i)*LocalHistoBin;
     
   // it must be a probability density -> normalisation by sum_check
   for(int i=0;i<LocalSizeSignal;i++) _HistoConv(0,i) /= sum_check;

//#if DEBUG_FEW
   if (WriteAllInfo) {
   
      dblarray Histo(MAXHISTONBPOINT);
      for (int i=0; i<LENGTH_FIRST_HISTO; i++) 
	  Histo(i) = _HistoConv(0,i);        
      show_param ("End Convol Compute", 0, 
                  _HistoBound(0,0), _HistoBound(0,1),
		  LocalHistoBin, LocalSizeSignal, Histo); 
		
      char FileName[256];
      dblarray Courbe (LocalSizeSignal);
      for (int i=0;i<LocalSizeSignal;i++) Courbe(i)=_HistoConv(0,i);
      sprintf (FileName, "Histo_%d.fits", 0);
      fits_write_dblarr (FileName, Courbe);
   }
//#endif   
      
   // interp histo 
   if (_InterpHisto) {
      fltarray Courbe (LocalSizeSignal);
      for (int i=0;i<LocalSizeSignal;i++) Courbe(i)=_HistoConv(0,i);
      _RedHisto->setRedHisto (0, Courbe, _HistoBound(0,0), _HistoBound(0,1),
	                      _HistoBin (0,0), (int) _HistoBin (0,1)); 
   } 
   
/*_HistoConv.init(0.);   
fltarray Data;
fits_read_fltarr("Aba_firsthisto.fits", Data);  
for (int i=0;i<Data.n_elem();i++) _HistoConv(0,i) = Data(i);
_HistoBound(0,0) = -0.040418;
_HistoBound(0,1) = 1.456781;
_HistoBin (0,0) = 0.000731;
_HistoBin (0,1) = 511;*/    
    		  
}

void FewEventPoisson::histo_convolution (Bool WriteAllInfo) {


   // init var
   //fltarray Signal1d(MAXHISTONBPOINT);
   float LocalBin;
   int NbFiltrage=0;

   // length of first histogram h(0,*)
   int HistoNbPoint = LENGTH_FIRST_HISTO;
   
   // Reduced Histo used for computation
   // init with h0 (_HistoConv(0,*))
   fltarray ReducedHisto (MAXHISTONBPOINT);
   ReducedHisto.init(0.);
   for (int i=0; i<LENGTH_FIRST_HISTO; i++)
      ReducedHisto(i) = _HistoConv(0,i);
   
   // init ReducedMin and ReducedMax value for parameters calculus
   // for the first histo min and _HistoBound(0,0) are the same (id max)
   float ReducedMin = _HistoBound(0,0);
   float ReducedMax = _HistoBound(0,1);
   LocalBin = _HistoBin (0,0);
   
   // for all autoconv 1::_NbAutoConv (index 1 to _NbAutoConv+1 in _Histo...)
   //------------------------------------------------------------------------
   for (int nbConv=1; nbConv<_NbAutoConv+1; nbConv++) {

//#if DEBUG_FEW 
      if (WriteAllInfo) {
         cout << "Histogram : "  << nbConv << endl;
	 cout << "==============================" << endl;
      }     
//#endif   
      //Current ind in _HistoConv tab
      long CurrentConvInd = (long) (3*(nbConv));      
 
      // legth of new histogram : 2 * length of last histo
      HistoNbPoint =  2*(HistoNbPoint-1)+1; //=2*HistoNbPoint-1
      
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
      for (long i=0; i<=(HistoNbPoint-1); i++) _HistoConv(CurrentConvInd, i) = 0.;
      for (long i=0; i<=(HistoNbPoint-1)/2; i++) {
         
	 for (long j=0; j<=i; j++)
	    _HistoConv(CurrentConvInd, i) += (ReducedHisto(j) * ReducedHisto(i-j));
      }
         
      long jbeg = 1;
      long jend = (HistoNbPoint-1)/2;  
      for (long i=(HistoNbPoint-1)/2+1; i<=(HistoNbPoint-1); i++) {
      
     	 for (long j=jbeg; j<=jend; j++)   	     
	     _HistoConv(CurrentConvInd, i) +=  (ReducedHisto(j) * ReducedHisto(jend+jbeg-j));
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
      _HistoBound(nbConv,0) = 2*ReducedMin;
      _HistoBound(nbConv,1) = 2*ReducedMax;
      _HistoBin (nbConv,1) = HistoNbPoint;
  //LocalBin = (_HistoBound(nbConv,1)-_HistoBound(nbConv,0)) / (HistoNbPoint-1); 		 
      _HistoBin (nbConv,0) = LocalBin;
//#if DEBUG_FEW 
      if (WriteAllInfo) { 
         dblarray Histo(HistoNbPoint);
	 for (int i=0; i<HistoNbPoint; i++) 
	    Histo(i) = _HistoConv(CurrentConvInd,i);     
         show_param ("Just after Convolution", nbConv, 
                     _HistoBound(nbConv,0), _HistoBound(nbConv,1),
		     _HistoBin (nbConv,0), _HistoBin (nbConv,1),
		     Histo);  
      }
//#endif   
      
      
      
      // Normalisation for density prob
      // ------------------------------
      float sum=0.;
      for (long i=0; i<HistoNbPoint; i++) 
         sum += _HistoConv(CurrentConvInd, i) * LocalBin;
      for (long i=0; i<HistoNbPoint; i++)
         _HistoConv(CurrentConvInd, i) /= sum;
//#if DEBUG_FEW 
      if (WriteAllInfo) { 
         dblarray Histo(HistoNbPoint);
	 for (int i=0; i<HistoNbPoint; i++) 
	    Histo(i) = _HistoConv(CurrentConvInd,i);     
         show_param ("Just after Normalisation", nbConv, 
                     _HistoBound(nbConv,0), _HistoBound(nbConv,1),
		     _HistoBin (nbConv,0), _HistoBin (nbConv,1),
		     Histo);  
      }
//#endif  	 
	 
	  
      // reduction of histogram 
      // ----------------------
      //
      // if (bin*val(i) < MIN_VAL_INTEG_HISTO) then val(i) is dstroyed
      // only in begin and end of histogram
      // reduced histogram is copied in new histigram table : 
      // _HistoBound (CurrentConvInd+1, *) , not for the last scale!!
      // IndexMin and IndexMax : indice of the new first and last point
      // of reduced histogram
      //
      // search for IndexMin and IndexMax
/*      ReducedMin = _HistoBound(nbConv,0);
      ReducedMax = _HistoBound(nbConv,1);     
      for (long i=0;i<HistoNbPoint;i++)
         ReducedHisto(i) = _HistoConv(CurrentConvInd, i);*/
	 
      
      if( WriteAllInfo ) {  
         std::cout << " ==> Just after integration reduction (histo "
	           << nbConv << ")" << std::endl;
      }    
      long Index = 0;
      while (    (_HistoConv(CurrentConvInd,Index)*_HistoBin (nbConv,0) < MIN_VAL_INTEG_HISTO) 
              && (Index < HistoNbPoint)) {
         Index++;
	 if( WriteAllInfo ) {
	    double prov = (_HistoConv(CurrentConvInd,Index)) * _HistoBin (nbConv,0);
	    std::cout << "       -- histo:" << _HistoConv(CurrentConvInd,Index)
	              << ", bin:" << _HistoBin (nbConv,0) 
		      << ", prod:" << prov << ", admin valmin:" 
		      << MIN_VAL_INTEG_HISTO << std::endl;
	 }
      }
      if( WriteAllInfo ) {
	 std::cout << "       -- nb pts min out : " << Index << std::endl;
      }
	 
      long IndexMin = Index;
      Index = 0;
      while (    (_HistoConv(CurrentConvInd,(HistoNbPoint-1)-Index)*_HistoBin (nbConv,0) < MIN_VAL_INTEG_HISTO) 
              && (Index < HistoNbPoint)) {
         Index++;
	 if( WriteAllInfo ) {
	    double prov = (_HistoConv(CurrentConvInd,(HistoNbPoint-1)-Index))*_HistoBin (nbConv,0);
	    std::cout << "       -- histo:" << _HistoConv(CurrentConvInd,(HistoNbPoint-1)-Index)
	              << ", bin:" << _HistoBin (nbConv,0) 
		      << ", prod:" << prov << ", admin valmin:" 
		      << MIN_VAL_INTEG_HISTO << std::endl;
	 }
      }
      if( WriteAllInfo ) {
	 std::cout << "       -- nb pts max out : " << Index << std::endl;
      }
      long IndexMax = HistoNbPoint-1-Index;
      
      if( WriteAllInfo ) {
         std::cout << "      -- indexMin:" << IndexMin << "(0), -- indexMax:" 
	           << IndexMax << "(" << HistoNbPoint-1 << ")" << std::endl;
      }
      //cout << " ==> IndexMin:" << IndexMin << ", IndexMax:" << IndexMax << endl;
      
      // decal histogram
      //ReducedHisto.init(0.);
      //for (long i=IndexMin; i<IndexMax+1; i++) 
      //   ReducedHisto(i-IndexMin) = _HistoConv(CurrentConvInd, i);
      // compute new number echant of reduced histo
      HistoNbPoint = IndexMax-IndexMin+1;
      //long IndexMax = HistoNbPoint-1;long IndexMin = 0;
      
      // Reduced signal must have an even number of echant for filter section
      // (HistoNbPoint must be a multiple of 4 plus 1 (for simplicity reasons)
      //
      if (HistoNbPoint > (MAXHISTONBPOINT-1)/2) {
         if (HistoNbPoint%4 == 1) {
      
            // OK HistoNbPoint = 4 * K + 1
	 
         } else {
      
            // we have to recup points until cond 1 of while test is ok
            while (    HistoNbPoint%4 != 1 
	           &&  HistoNbPoint<MAXHISTONBPOINT) {
         
	       if (IndexMax+1 <= MAXHISTONBPOINT-1) {
	          IndexMax++;HistoNbPoint++;
	       } else if (IndexMin > 0) {
	          IndexMin--;HistoNbPoint++;
	       }
 	       else {
	          cout << "!!!!!!!HistoNbPoint too large in reduction proc!!!" << endl;
	          break;
	       }
	    }      
         } 
      }
      //cout << " ==> IndexMin:" << IndexMin << ", IndexMax:" << IndexMax << endl;
      
      HistoNbPoint = IndexMax-IndexMin+1;
      // decal histogram
      ReducedHisto.init(0.);
      for (long i=IndexMin; i<IndexMax+1; i++) 
         ReducedHisto(i-IndexMin) = _HistoConv(CurrentConvInd, i);

      // new min, max, bin
      ReducedMin = _HistoBound(nbConv,0)+IndexMin*_HistoBin (nbConv,0);
      ReducedMax = _HistoBound(nbConv,1)-((HistoNbPoint-1)-IndexMax)*_HistoBin (nbConv,0);     
      

//#if DEBUG_FEW 
      if (WriteAllInfo) { 
         dblarray Histo(HistoNbPoint);
	 for (int i=0; i<HistoNbPoint; i++) 
	    Histo(i) = ReducedHisto(i);     
         show_param ("Just after Reduction", nbConv, 
                     ReducedMin, ReducedMax,
		     LocalBin, HistoNbPoint,
		     Histo);  
      }
//#endif  




      // decimation?
      //
      // new histogram at next scale will have 2*HistoNbPoint points
      // if 2*HistoNbPoint > MAXHISTONBPOINT) we have to decimate 
      // we used a filter with kernel [0.25, 0.5, 0.25] (function shape_signal)
      // 
      if (HistoNbPoint > (MAXHISTONBPOINT-1)/2) { 
         //cout << "  =>Filtrage..." << endl;       
         shape_signal (ReducedHisto);
	 HistoNbPoint = (HistoNbPoint-1)/2 + 1;
         LocalBin = (ReducedMax-ReducedMin) / (HistoNbPoint-1);	
	 NbFiltrage++;
      
//#if DEBUG_FEW 
         if (WriteAllInfo) { 
            dblarray Histo(HistoNbPoint);
	    for (int i=0; i<HistoNbPoint; i++) 
	       Histo(i) = ReducedHisto(i);     
            show_param ("Just after SHAPE", nbConv, 
                        ReducedMin, ReducedMax,
		        LocalBin, HistoNbPoint,
		        Histo);  
         }
//#endif       
      
      }
        
	
	
	
      // Normalisation for density prob
      //
//      sum=0;
//      for (long i=0; i<HistoNbPoint; i++) 
//         sum += ReducedHisto(i)*LocalBin;
//      for (long i=0; i<HistoNbPoint; i++)
//         ReducedHisto(i) /= sum;
           

      // interp histo
      if (_InterpHisto)
         _RedHisto->setRedHisto (nbConv, ReducedHisto, ReducedMin, ReducedMax,
	                         LocalBin, HistoNbPoint);
   }
   if (WriteAllInfo) {
      cout << "Nb Filtrage :" << NbFiltrage << endl;   
      cout << "Nb Convol   :" << _NbAutoConv << endl;
   }
}

void FewEventPoisson::histo_normalisation (Bool WriteAllInfo) {

   
   for( int nbConv=0; nbConv<_NbAutoConv+1; nbConv++ ) {

      // compute mean & standard deviation of values of the wavelet coef
      // ---------------------------------------------------------------
      
      // local var
      // ---------
      double mean=0.0;
      double std=0.0;
      double localBin   = _HistoBin(nbConv,0);
      double localNbPts = _HistoBin(nbConv,1);
      double localMinBound = _HistoBound(nbConv,0);
      
      // compute mean
      // ------------
      for( int j=0; j<localNbPts; j++ ) {
         double real_value = ((double) j)*localBin + localMinBound;
	 double surface = _HistoConv(3*nbConv,j)*localBin;
         mean += surface*real_value;
         std  += surface*real_value*real_value;
      }   
      
      _HistoMean(nbConv) = mean;
      std = std - mean*mean;
      if (std < 0.0)  std = 0.0;
      else std = sqrt(std);
      _HistoSigma(nbConv) = std;
      
//# if DEBUG_FEW     
      if (WriteAllInfo) {
     
         cout << "Scale:" << nbConv << endl;
	 cout << "  -min:" << _HistoBound(nbConv,0)
	      << ", -max:" << _HistoBound(nbConv,1)
	      << ", -mean:" << _HistoMean(nbConv)
	      << ", -sigma:" << _HistoSigma(nbConv) << endl;
      }
//#endif
      
      for (int j=0; j<_HistoBin(nbConv,1);j++) {
        
	 /*--- compute real values of the wavelet coefficients ---*/
         _HistoConv(3*nbConv+1,j) = _HistoDistrib(3*nbConv+1,j) =   
                  (((float) j)*localBin + _HistoBound(nbConv,0) 
                   - _HistoMean(nbConv))  / _HistoSigma(nbConv);
 
	
	 /*--- compute histogram associated with the reduced coefficients ---*/
         _HistoConv(3*nbConv+2,j) = _HistoConv(3*nbConv,j)*_HistoSigma(nbConv) ;    
            
	    
	  
      }
      
      /*if (WriteAllInfo && nbConv==20) {
         //char FileName[256];
	 float HistoNbPoint=_HistoBin(nbConv,1);
         //fltarray Courbe (HistoNbPoint);
         //for (int i=0;i<HistoNbPoint;i++) Courbe(i)=_HistoConv(3*nbConv+1,i);
         //sprintf (FileName, "xHisto_%d.fits", nbConv);
	 for (int i=0;i<HistoNbPoint;i++) cout << _HistoConv(3*nbConv+1,i) << "/" ;
         //fits_write_fltarr (FileName, Courbe);
         //for (int i=0;i<HistoNbPoint;i++) Courbe(i)=_HistoConv(3*nbConv+2,i);
         //sprintf (FileName, "yHisto_%d.fits", nbConv);
         //fits_write_fltarr (FileName, Courbe);	    
      }	*/         
   }
}

void FewEventPoisson::new_histo_normalisation (Bool WriteAllInfo) {

   for( int nbConv=0; nbConv<_NbAutoConv+1; nbConv++ ) {

      // compute mean & standard deviation of values of the wavelet coef
      // ---------------------------------------------------------------
      
      // local var
      // ---------
      double mean=0.0;
      double std=0.0;
      double nbPts=0.0;
      double localBin   = _HistoBin(nbConv,0);
      double localNbPts = _HistoBin(nbConv,1);
      double localMinBound = _HistoBound(nbConv,0);
      
      // compute mean
      // ------------
      for( int j=0; j<localNbPts; j++ ) {
         double real_value = ((double) j)*localBin + localMinBound;
         mean += _HistoConv(3*nbConv,j)*real_value;
	 nbPts += _HistoConv(3*nbConv,j);
      }  
      mean = mean / nbPts;
      
      // compute sigma
      // -------------
      for( int j=0; j<localNbPts; j++ ) {
         double center_value = ((double) j)*localBin + localMinBound - mean;
	 double surface = _HistoConv(3*nbConv,j);
         std  += surface*center_value*center_value;
      }   
      
      _HistoMean(nbConv) = mean;
      if (std < 0.0)  std = 0.0;
      else std = sqrt(std);
      _HistoSigma(nbConv) = std;
      
//# if DEBUG_FEW     
      if (WriteAllInfo) {
     
         cout << "Scale:" << nbConv << endl;
	 cout << "  -min:" << _HistoBound(nbConv,0)
	      << ", -max:" << _HistoBound(nbConv,1)
	      << ", -mean:" << _HistoMean(nbConv)
	      << ", -sigma:" << _HistoSigma(nbConv) << endl;
      }
//#endif
      
      for (int j=0; j<_HistoBin(nbConv,1);j++) {
        
	 /*--- compute real values of the wavelet coefficients ---*/
         _HistoConv(3*nbConv+1,j) = _HistoDistrib(3*nbConv+1,j) =   
                  (((float) j)*localBin + _HistoBound(nbConv,0) 
                   - _HistoMean(nbConv))  / _HistoSigma(nbConv);
	
	 /*--- compute histogram associated with the reduced coefficients ---*/
         _HistoConv(3*nbConv+2,j) = _HistoConv(3*nbConv,j)*_HistoSigma(nbConv) ;    
            
	    
	  
      }
      
      /*if (WriteAllInfo && nbConv==20) {
         //char FileName[256];
	 float HistoNbPoint=_HistoBin(nbConv,1);
         //fltarray Courbe (HistoNbPoint);
         //for (int i=0;i<HistoNbPoint;i++) Courbe(i)=_HistoConv(3*nbConv+1,i);
         //sprintf (FileName, "xHisto_%d.fits", nbConv);
	 for (int i=0;i<HistoNbPoint;i++) cout << _HistoConv(3*nbConv+1,i) << "/" ;
         //fits_write_fltarr (FileName, Courbe);
         //for (int i=0;i<HistoNbPoint;i++) Courbe(i)=_HistoConv(3*nbConv+2,i);
         //sprintf (FileName, "yHisto_%d.fits", nbConv);
         //fits_write_fltarr (FileName, Courbe);	    
      }	*/         
   }
}


void FewEventPoisson::histo_distribution (Bool WriteAllInfo) {


   for (int nbConv=0; nbConv<_NbAutoConv+1; nbConv++) {
   
      double LocalBin     = _HistoBin(nbConv,0);
      double LocalNbPoint = _HistoBin(nbConv,1);
      double LocalSigma   = _HistoSigma(nbConv);
      
      //Probability Reduced Values
      //
      // init first values
      _HistoDistrib(3*nbConv,0)   = _HistoConv(3*nbConv,0)*LocalBin;
      _HistoDistrib(3*nbConv+2,0) = _HistoConv(3*nbConv+2,0)*LocalBin/LocalSigma;
      
      // integ other val (Dist(i) = Dist(i-1) + val(i)*bin)
      for (int nbp=1; nbp<LocalNbPoint; nbp++) {
         
	 _HistoDistrib(3*nbConv,nbp)   =   _HistoDistrib(3*nbConv,nbp-1) 
	                          + _HistoConv(3*nbConv,nbp)*LocalBin;

         _HistoDistrib(3*nbConv+2,nbp) =   _HistoDistrib(3*nbConv+2,nbp-1) 
	                          + _HistoConv(3*nbConv+2,nbp)*LocalBin/LocalSigma;
      
      }
      
      // normalisation (last value must be one)
      for( int nbp=0; nbp<LocalNbPoint; nbp++) {
          _HistoDistrib( 3*nbConv, nbp ) =
	     _HistoDistrib( 3*nbConv, nbp ) /  _HistoDistrib( 3*nbConv, LocalNbPoint-1 );
	  _HistoDistrib( 3*nbConv+2, nbp ) =  
	     _HistoDistrib( 3*nbConv+2, nbp ) / _HistoDistrib( 3*nbConv+2, LocalNbPoint-1 );
      }
      
      
   }
//#if DEBUG_FEW   
   if (WriteAllInfo) {
      for (int nbConv=0; nbConv<_NbAutoConv+1; nbConv++) {
         cout << "Scale:" << nbConv << endl;
	 cout << "Min:" << _HistoBound(nbConv,0)
	      << ", Max:" << _HistoBound(nbConv,1)
	      << ", bin:" << _HistoBin(nbConv,0)/_HistoSigma(nbConv) << endl;
	 cout << "d(0):" << _HistoDistrib(3*nbConv,0)
	      << ", d(" << (_HistoBin(nbConv,1)-1)/2 << "):" 
	                << _HistoDistrib(3*nbConv,(int)((_HistoBin(nbConv,1)-1)/2))
	      << ", d(" << _HistoBin(nbConv,1)-1 << "):" 
	                << _HistoDistrib(3*nbConv,(int)(_HistoBin(nbConv,1)-1)) << endl;		
      }
      fits_write_dblarr("NewDistrib.fits", _HistoDistrib);
   }
//#endif
}




void FewEventPoisson::shape_signal (fltarray &Signal1d) {
   // size of signal is reduced by two

   // init var
   long LocalSignalLength = Signal1d.nx();
   fltarray Filter_Signal1d(LocalSignalLength);
  
   // filter signal
   for (int i=1;i<=(LocalSignalLength-1)/2-1;i++)
      Filter_Signal1d(2*i) = 0.5  * Signal1d(2*i) + 
                             0.25 * (Signal1d(2*i-1) + Signal1d(2*i+1));
   
   Filter_Signal1d(0)	= 0.5*(Signal1d(0) + Signal1d(1));		    
   Filter_Signal1d(LocalSignalLength-1) = 0.5*(  Signal1d(LocalSignalLength-2) 
                                               + Signal1d(LocalSignalLength-1));
						
   // Sampling and put zeros
   Signal1d.init(0.);
   for (int i=0;i<=(LocalSignalLength-1)/2;i++) 
      Signal1d(i)= Filter_Signal1d(2*i);
}



void FewEventPoisson::show_param (char* text, float nbconv, float min, 
                                  float max, float bin, float nbechant, 
				  dblarray& Histo) {
				    
				    
   cout << " ==> " << text << " (histo " << nbconv << ")" << endl;
   cout	<< "      -- min:" << min << ", -- max:" << max  
        << ", -- bin:" << bin << ",   (NBP:" << nbechant << ")" << endl;
   //cout << "      -- mean (histogram):" << Histo.mean()
   //     << ", -- sigma (histogram):" << Histo.sigma() << endl;
	
   // compute sigma an mean of realisation.. (wavelet space)	
   float mean=0., std=0.;
   for (int j=0; j<nbechant;j++) {
      
      float real_value = ((float) j)*bin + min;
      mean += Histo(j)*bin*real_value;
   }
   
   for (int j=0; j<nbechant;j++) { 
   
      float real_value = ((float) j)*bin + min;  
      std +=  (   ((Histo(j)*bin*real_value*real_value-mean)-mean)
                * ((Histo(j)*bin*real_value*real_value-mean)-mean));
   } 

   //std = std - mean*mean;
   if (std < 0.0)  std = 0.0;
   else std = sqrt(std);
   //cout << "      -- mean (wavelet space):" << mean 
   //     << ", -- sigma (wavelet space):" << std << endl;
}
 

void FewEventPoisson::histo_threshold (double Epsilon, Bool WriteAllInfo) {


   if( WriteAllInfo ) {
      cout << "FewEventPoisson::histo_threshold" << Epsilon << endl;
      cout << "  epsilon:" << Epsilon << endl;
   }
   
   for (int nbConv=0; nbConv<_NbAutoConv+1; nbConv++) {
   
      // serach indice of max threshold (RepFunc >= 1-eps) 
      // and min threshold (RepFunc <= eps)
      // 
      long IndexMax=0;
      long IndexMin=(long)_HistoBin(nbConv,1)-1;
      // RepFuc = 1-eps in [IndexMax-1:IndexMax]
      while (    (_HistoDistrib(3*nbConv+2,IndexMax) < (1-Epsilon) )
              && (IndexMax<_HistoBin(nbConv,1)-1)) {IndexMax++;}
      // RepFuc = eps in [IndexMin:IndexMin+1]      
      while (    (_HistoDistrib(3*nbConv+2,IndexMin) > Epsilon )
              && (IndexMin>0)) {IndexMin--;}
     // test	      
      if ((IndexMax==0) || (IndexMin==_HistoBin(nbConv,1)-1)) {
	  cout << "ERROR : bad sampling of the distribution function " << endl;
	  exit (-1);
      }   
      if (WriteAllInfo) {
         cout << "  " << nbConv << ", SeuilMin := [" << _HistoDistrib(3*nbConv+1,IndexMin)
	                            << "," <<  _HistoDistrib(3*nbConv+1,IndexMin+1) << "]"
                                    << ", SeuilMax := [" << _HistoDistrib(3*nbConv+1,IndexMax-1) << ","
				    << _HistoDistrib(3*nbConv+1,IndexMax) << "]"
      			            << endl; 
         cout << "  " << nbConv << ", repart := [" << _HistoDistrib(3*nbConv+2,IndexMin)
	                            << "," <<  _HistoDistrib(3*nbConv+2,IndexMin+1) << "]"
                                    << ", XMax := [" << _HistoDistrib(3*nbConv+2,IndexMax-1) << ","
				    << _HistoDistrib(3*nbConv+2,IndexMax) << "]"
      			            << endl; 
      }

      // inter between IndexMax and IndexMax-1
      // X1-eps = Xmax + (1-eps-Ymax)/a,  a=(Ymax-1 - Ymax)/(Xmax-1 - Xmax)
      // 
      double a =   (_HistoDistrib(3*nbConv+1,IndexMax-1)-_HistoDistrib(3*nbConv+1,IndexMax))
	         / (_HistoDistrib(3*nbConv+2,IndexMax-1)-_HistoDistrib(3*nbConv+2,IndexMax));
      double b =        _HistoDistrib(3*nbConv+1,IndexMax-1) 
                 - a * _HistoDistrib(3*nbConv+2,IndexMax-1);
      
         
      //_Threshold (nbConv, 1) =   _HistoDistrib(3*nbConv+1,IndexMax) 
      //                          + (1-Epsilon-_HistoDistrib(3*nbConv+2,IndexMax))/a;
      _Threshold (nbConv, 1) = a*(1-Epsilon) + b;
      
      //_Threshold (1, nbConv) = (1-Epsilon-b)/a;
      if (WriteAllInfo) {
         cout << "  " << "a(max):" << a << ", b(max):" << b <<endl;
	 cout << "  " << "SeuilMax := [" << _HistoDistrib(3*nbConv+1,IndexMax-1) << ","
	                            << _Threshold (nbConv, 1) << "," 
				    << _HistoDistrib(3*nbConv+1,IndexMax) << "]"
      			            << endl;
         cout << "  " << "XMax := [" << _HistoDistrib(3*nbConv+2,IndexMax-1) << ","
	                            << 1-Epsilon << ","
				    << _HistoDistrib(3*nbConv+2,IndexMax) << "]"
      			            << endl; 
      }
      // inter between IndexMin and IndexMin+1
      // Xeps = Xmin + (eps-Ymin)/a,  a=(Ymin+1 - Ymin)/(Xmin+1 - Xmin)
      a =        (_HistoDistrib(3*nbConv+1,IndexMin+1)-_HistoDistrib(3*nbConv+1,IndexMin))
               / (_HistoDistrib(3*nbConv+2,IndexMin+1)-_HistoDistrib(3*nbConv+2,IndexMin));
      b =        _HistoDistrib(3*nbConv+1,IndexMin+1) 
           - a * _HistoDistrib(3*nbConv+2,IndexMin+1);
	       
      //_Threshold (nbConv, 0) =   _HistoDistrib(3*nbConv+1,IndexMin) 
	//                          + (Epsilon-_HistoDistrib(3*nbConv+2,IndexMin))/a; 
      _Threshold (nbConv, 0) = a*(Epsilon) + b;
	
	 
      if (WriteAllInfo) {
         cout << "  " << "a(min):" << a << endl;
      }
      

      if (WriteAllInfo) {
         cout << "Nbr of histogram autoconv.:" << nbConv << ", SeuilMin := " << _Threshold (nbConv, 0)
                                    << ", SeuilMax := " << _Threshold (nbConv, 1)
      			            << endl; 
      }     
   }
    if (WriteAllInfo) 
      io_write_ima_float (Name_Threshold, _Threshold);
}

 
void FewEventPoisson::find_threshold (double Epsilon, Bool WriteAllInfo) {


   if (_InitOk == False) (*this).compute_distribution(WriteAllInfo);
   histo_threshold (Epsilon, WriteAllInfo);
}

 

float FewEventPoisson::event_prob (float ValRed, int NEventReal, dblarray& Tab) {

   // compute "power" : 2^(power-1) < NEventReal <= 2^(power)
   //
   int NEvent=1;
   int PowerAfter=0;
   int PowerBefore=0;
   while (NEvent<NEventReal) {PowerAfter ++; NEvent *= 2;}
   
   // test if the event position must be interpolated
   Bool InterpEvent = False;
   if (PowerAfter > DEF_NBR_AUTOCONV) PowerAfter=DEF_NBR_AUTOCONV;
   else { 
      if ((PowerAfter > 0) && (NEvent != NEventReal)) {
         PowerBefore = PowerAfter -1;
         InterpEvent = True;
      }
   }
      
   // Computes the probability related to a reduced wavelet
   // coefficients ValRed, and NEvent events
   //
   // search IndCoef so that  RepFunc (PowerAfter, IndCoef) <  ValRed
   int IndCoef = 0;
   long LocalNbPoint = (long)_HistoBin(PowerAfter,1);
   while (    (IndCoef < LocalNbPoint-1) 
           && ( Tab(PowerAfter*3+1, IndCoef) < ValRed)) {IndCoef++;}
	   
   // compute prob RepFunc <  ValRed for PowerAfter  
   //   
   float ValRedBef, ValRedAft, ValBef, ValAft, ProbBefore, ProbAfter, Weight;
   if (IndCoef==0) {
      ProbAfter = Tab(3*PowerAfter+2, IndCoef);
   } else if (IndCoef==LocalNbPoint) {
      ProbAfter = Tab(3*PowerAfter+2, IndCoef-1);
   } else {
      ValRedBef = Tab(3*PowerAfter+1, IndCoef-1);
      ValRedAft = Tab(3*PowerAfter+1, IndCoef);
      ValAft = Tab(3*PowerAfter+2, IndCoef);
      ValBef = Tab(3*PowerAfter+2, IndCoef-1);
      if (ABS(ValAft - ValBef) > FLOAT_EPSILON) 
      Weight = 1. - (ValRed - ValRedBef) / (ValRedAft - ValRedBef);
      else Weight = 0.5;
      ProbAfter = Weight * ValBef + (1. - Weight) * ValAft;
   }     
      
   // if no interpolation, return the "true value" ProbAfter
   if (!InterpEvent) return (ProbAfter); 
   
   // intrepolate between PowerAfter and PowerBefore
   //    
   // Computes the probability related to a reduced wavelet
   // coefficients ValRed, and n_event/2  events
   IndCoef = 0;
   while (    (IndCoef < _HistoBin(PowerBefore,1)-1) 
           && ( Tab(PowerBefore*3+1, IndCoef) < ValRed)) {IndCoef++;}
   
   // compute prob RepFunc <  ValRed for PowerAfter  
   //   	     
   if (IndCoef==0) {
      ProbBefore = Tab(3*PowerBefore+2, IndCoef);
   } else if (IndCoef==_HistoBin(PowerBefore,1)) {
      ProbBefore = Tab(3*PowerBefore+2, IndCoef-1);
   } else {             
      ValRedBef = Tab(3*PowerBefore+1, IndCoef-1);
      ValRedAft = Tab(3*PowerBefore+1, IndCoef);
      ValAft = Tab(3*PowerBefore+2, IndCoef);
      ValBef = Tab(3*PowerBefore+2, IndCoef-1);
      if (ABS(ValAft - ValBef) > FLOAT_EPSILON) 
      Weight = 1. - (ValRed - ValRedBef) / (ValRedAft - ValRedBef);
      else Weight = 0.5;
      ProbBefore = Weight * ValBef + (1. - Weight) * ValAft;
              /*cout << "Val red = " << ValRed << endl;
              cout << "ProbBefore = " << ProbBefore << endl;
              cout << "ProbAfter  = " << ProbAfter  << endl;
              cout << "n_event   = " <<  NEvent  << endl;
              cout << "n_event_real  = " <<  NEventReal  << endl;
              cout << " Weight  = " <<   Weight  << endl;*/ 
   }
   Weight = 2.* (NEvent - NEventReal) / (float)NEvent;
   float ValReturn = Weight * ProbBefore + (1. - Weight) * ProbAfter;   
   
   //std::cout << "FewEventPoisson::prob reduced coef : " << ValRed << ", event number : " 
   //          << NEventReal << ",  prob = " << ValReturn << std::endl;
      
   return ( ValReturn); 
      
}

void FewEventPoisson::set_interp_histo () {
   _InterpHisto = True;
   _RedHisto = new cReductHisto;
}  
 
 
float FewEventPoisson::prob(float ValRed, int NEvent) {
   return event_prob(ValRed, NEvent, _HistoConv);
}
 
float FewEventPoisson::repartition(float ValRed, int NEvent) {
   return event_prob(ValRed, NEvent, _HistoDistrib);
}
 
float FewEventPoisson::a_trou_prob(float Coef, int NEvent, int Scale) {
   
    float RedCoef,alpha,sc;
    for (sc=0,alpha=1.0; sc<Scale; sc++) alpha *= 4.;
    RedCoef = Coef * alpha /  sqrt((float) NEvent) / SIGMA_BSPLINE_WAVELET;
//cout << "!!!!!!!!!!!!!!! Scale:" << Scale << ", NEvent:" << NEvent 
//     << ", sigma red:" << alpha /  sqrt((float) NEvent) / SIGMA_BSPLINE_WAVELET
//     << endl;
    return  prob(RedCoef, NEvent);
}
 
float FewEventPoisson::a_trou_repartition(float Coef, int NEvent, int Scale) {
    
    float RedCoef,alpha,sc;
    for (sc=0,alpha=1.0; sc< Scale; sc++) alpha *= 4.;
    RedCoef = Coef * alpha /  sqrt((float) NEvent) / SIGMA_BSPLINE_WAVELET;
//#if DEBUG_FEW  
cout << "!!!!!!!!!!!!!!! Scale:" << Scale << ", NEvent:" << NEvent 
     << ", sigma red:" << alpha /  sqrt((float) NEvent) / SIGMA_BSPLINE_WAVELET
     << endl;
//#endif
     return  repartition(RedCoef, NEvent);
}
 
 
 


cReductHisto::cReductHisto() {

   _NbNewHisto    = 4;    //h^3, h^5, h^6 and h^7
   _MaxNbRedHisto = 8;    // max red histo is 7
   
   _RedHisto.alloc (_MaxNbRedHisto, MAXHISTONBPOINT);
   _RedHisto.init();
   _RedBound.alloc (_MaxNbRedHisto, 2);
   _RedBin.alloc (_MaxNbRedHisto, 2);
   _TabInd.alloc (MAXNBNEWHISTO);
   
   _HistoConv.alloc (3*(_MaxNbRedHisto+1), MAXHISTONBPOINT);
   _HistoConv.init(0.);   
   
   // !! new histo max is 16 !! else you have to modify Tab...
   for (int i=0;i<MAXNBNEWHISTO;i++) _TabInd(i) = Tab[i];
   _UsedHistoConv.alloc (MAXNBNEWHISTO, 2);  
   for (int j=0;j<MAXNBNEWHISTO;j++) {
      _UsedHistoConv(j,0) = FirstHisto[j];
      _UsedHistoConv(j,1) = SecondHisto[j];
   }
}

void cReductHisto::setRedHisto (int IndHisto, fltarray& pRedHisto, float Min,
                               float Max, float Bin, int NbPoint) {
		
   long IndInTab = (long)(pow (2., (double)IndHisto));	       
   if (IndInTab<0 || IndInTab>_MaxNbRedHisto-1) {
   } else {
      for (int i=0;i<pRedHisto.n_elem();i++) 
         _RedHisto(IndInTab,i) = pRedHisto(i);  
      _RedBound(IndInTab,0) = Min;
      _RedBound(IndInTab,1) = Max;
      _RedBin(IndInTab,0) = Bin;
      _RedBin(IndInTab,1) = NbPoint;
   }
}
 

void cReductHisto::compute_distribution (Bool WriteAllInfo) {


   // computes the interp auto-convolued histogramms
   //
   if  (WriteAllInfo) cout << "Compute the interp autoconvolutions of the histogram ... " << endl;
   histo_convolution (WriteAllInfo);
   
   // computes the normalisation
   //
   if  (WriteAllInfo) cout << "Compute normalisation ... " << endl;   
   //histo_normalisation (WriteAllInfo);
   
   // writing debug result
   if (WriteAllInfo == True) {
//      fits_write_fltarr (Name_HistoConv, _HistoConv);
//      fits_write_fltarr (Name_HistoMean, _HistoMean);
//      fits_write_fltarr (Name_HistoSigma, _HistoSigma);
   } 
   
   // computes the distribution function of histogramms at each leve
   //
   //histo_distribution (WriteAllInfo);
}



void cReductHisto::histo_convolution (Bool WriteAllInfo) {
      
   for (int i=0; i<_NbNewHisto; i++) {
   
      // indice of current hisot
      long CurrentConvInd = (long) (_TabInd(i));
      
      // indice of histo used in autoconv
      // the first histo is always shorter than the second
      //
      long IndRedHistoOne = (long)_UsedHistoConv(i,0);
      long IndRedHistoTwo = (long)_UsedHistoConv(i,1);
      long NbPtRedHistoOne;
      long NbPtRedHistoTwo;
      
      if (IndRedHistoOne < IndRedHistoTwo) {
      
         NbPtRedHistoOne = (long)_RedBin(IndRedHistoOne,1);
         NbPtRedHistoTwo = (long)_RedBin(IndRedHistoTwo,1);

      } else {
      
         NbPtRedHistoTwo = (long)_RedBin(IndRedHistoOne,1);
         NbPtRedHistoOne = (long)_RedBin(IndRedHistoTwo,1);      
         IndRedHistoTwo  = (long)_UsedHistoConv(i,0);
	 IndRedHistoOne  = (long)_UsedHistoConv(i,1);
      }
       
      // legth of new histogram : (length 1 histo -1) + (length 2 histo-1) +1
      // --------------------------------------------------------------------
      long HistoNbPoint = NbPtRedHistoOne-1 + NbPtRedHistoTwo-1 + 1;
      
      // compute histogram convolution
      // -----------------------------
      // 
      //
      for (long i=0; i<=(HistoNbPoint-1); i++) _HistoConv(CurrentConvInd, i) = 0;
      
      
      /*for (long i=0; i<=NbPtRedHistoOne-1; i++) {
         fltarray First(i);
	 for (long j=0;j<i;j++) First(j) = RedHisto(IndRedHistoOne,j)
	 fltarray InvSecond(i);
         for (long j=NbPtRedHistoTwo-1;j<NbPtRedHistoTwo-1-i;j--) InvSecond(j) = RedHisto(IndRedHistoTwo,j)
	 
	 for (j=0;j<i;j++)
	    _HistoConv(CurrentConvInd, i) +=  (First(j) * InvSecond(j));
      }
      
      long length = NbPtRedHistoOne;
      fltarray First(length);
      for (long j=0;j<length;j++) First(j) = RedHisto(IndRedHistoOne,j);
      long jend = NbPtRedHistoTwo-1+1;
      long jbeg = jend-length+1;
      for (long i=NbPtRedHistoOne; i<=NbPtRedHistoTwo-1; i++) {
      
	 fltarray InvSecond(length);
	 for (long j=jbeg; j<=jend; j++) InvSecond(j-jbeg) = RedHisto(IndRedHistoTwo,jbeg);
	 
	 for (long j=0; j<length; j++)
	    _HistoConv(CurrentConvInd, i) += (First(j) * InvSecond(j));
	    
	 jbeg--;
	 jend--;      
      }
      
      long jbeg = 1;
      long jend = NbPtRedHistoOne-1;
      for (long i=NbPtRedHistoTwo; i<=HistoNbPoint-1; i++) {
         
	 long length = jend-jbeg+1
	 fltarray First(length);  
	 for (j=jbeg;j<=jend;j++) First(j) = RedHisto(IndRedHistoOne,j);
	 fltarray InvSecond(length);    
         for (j=jbeg;j<=jend;j++) InvSecond(j) = RedHisto(IndRedHistoTwo,j-jbeg);
	 
	 for (long j=0; j<length; j++)
	    _HistoConv(CurrentConvInd, i) += (First(j) * InvSecond(j));
	  
	 jbeg++;
      } */

   }
}


void event_one_scale (intarray&   EventSignal, 
                        int         s, 
		        intarray&   EventCount, 
		        type_border Border,
		        Bool        WriteAllInfo) {

   // support at scale s 
   int Size = (int) (pow((double)2., (double)(s+2)) + 0.5);
   
   // test border
   if (Border != I_MIRROR  && Border != I_CONT)
      cout << "Only Border I_MIRROR and I_CONT are implemented " 
           << "in event count (func event_one_scale)!" << endl;
	   
   int NbPoint = EventSignal.n_elem();
   for (int i=0; i<NbPoint; i++) {
       
      int Total = 0;
      if (Border == I_MIRROR) {
         for (int k =-Size; k <= Size; k++) {
      
            int ind = (i+k < 0 ? -(i+k) : i+k);
	    ind = (ind > (NbPoint-1) ? (NbPoint-1)-(ind - (NbPoint-1)) : ind);
            Total += EventSignal(ind);
         }
      }
      if (Border == I_CONT) {
         for (int k =-Size; k <= Size; k++) {
      
            int ind = (i+k < 0 ? 0 : i+k);
	    ind = (ind > (NbPoint-1) ? (NbPoint-1) : ind);
            Total += EventSignal(ind);
         }
      }     
      

      EventCount(i) = Total;
   }
    
   if (WriteAllInfo) {     
      cout << "Scale:" << s << ", Size:" << 2*Size+1 << endl;
      char FileName[256];
      sprintf (FileName, "EventCount_%d", s);
      fltarray prov(EventSignal.n_elem());
      for (int i=0;i<EventSignal.n_elem();i++) prov(i)=(float)EventCount(i);
      fits_write_fltarr(FileName,prov); 
   }  
}

