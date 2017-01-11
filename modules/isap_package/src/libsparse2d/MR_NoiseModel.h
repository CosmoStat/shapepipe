/******************************************************************************
**                   Copyright (C) 1995 by CEA
*******************************************************************************
**
**    UNIT
**
**    Version: 3.2
**
**    Author: Jean-Luc Starck
**
**    Date:  96/06/13 
**    
**    File:  MR_NoiseModel.h
**
*******************************************************************************
**
**    DESCRIPTION  Noise Modelling Class
**    ----------- 
**                 
**        MRNoiseModel NoiseModel(TNoise, MR_Data);
**     or MRNoiseModel NoiseModel(TNoise, Nl, Nc, ScaleNumber, Transform)
**        allocates space for Noise Modelling
**
**     NoiseModel.model(Imag); 
**        Estimates the detection level at each scale and creates the
**        multiresolution support
**     
**     NoiseModel.model(MR_Data);
**        Estimates the detection level at each scale and creates the
**        multiresolution support. MR_Data is the multiresolution of the
**        image already transformed.
**        
**     MRNoiseModel.support(s,i,j) return True is something has been detected
**        at this scale and at this position
**
**     MRNoiseModel.signif(s,i,j, Val) return True is this value is upper
**        the detection level at this scale and at this position.
**
**     Bool DilateSupport;   // if true, then dilate the support
**     Bool SupIsol;         // if true, then supress isolated pixel
**     Bool MinEvent;        // if true, then detection with a minimum 
**                           //  number of events
**     int MinEventNumber;   // minimum number of event for a detection
**     Bool BadPixel;        //  if true, then take into account bad pixels
**     float BadPixelVal;    //  bad pixel value
**
**     int Type_SigmaDetectionMethod; // method for sigma noise detection
**     float CCD_Gain;         // CCD gain 
**     float CCD_ReadOutSigma; // CCD read-out noise standard deviation 
**     float CCD_ReadOutMean;  // CCD read-out noise mean 
**     float SigmaNoise;       // noise standard deviation (Gaussian Noise) 
**     float NSigma;           // Detection at NSigma * Noise 
**
******************************************************************************/

#ifndef _CNOISEMODEL_H_
#define _CNOISEMODEL_H_
#include "IM_Noise.h"
#include "MR_Noise.h"
#include "MR_Abaque.h"
#include "IM_Deconv.h"
#include "MR_Rayleigh.h"
#include "MR_CorrNoise.h"
#include "Mr_FewEvent.h"
#include "Mr_FE.h"
#include "MR_Threshold.h"  

#define DEF_MINEVENT_NUMBER 4 
#define DEF_NSIGMA_DETECT 3.
#define DEF_FIRST_DETECT_SCALE 0


#define NBR_SIGMA_METHOD 4
enum type_sigma_method {SIGMA_MEDIAN, SIGMA_BSPLINE, SIGMA_CLIPIMA,
                        SIGMA_SUPPORT,SIGM_UNDEFINED=-1};
#define DEF_SIGMA_METHOD SIGMA_MEDIAN

#define VAL_SupNull 0
#define VAL_SupOK 1
#define VAL_SupDill 2
#define VAL_SupEdge 3
#define VAL_SupLastOK 9
#define VAL_SupKill 10
#define VAL_SupMinEv 11
#define VAL_SupFirstScale 12

#define MAXBAND MAX_SCALE*3

const type_noise_speckle TSPECKLE = NOISE_LOG_RAYLEIGH;
 
class MRNoiseModel {
         int NbrScale, Nl, Nc, NbrBand;
         intarray  TabNl;
         intarray  TabNc;
         intarray  TabPos;
	 intarray  TabBandScale;
         int Size;
         unsigned char *TabSupport;
         // multiresolution support of an image
         //    TabSupport(s,i,j) = 0  no detection
         //    TabSupport(s,i,j) = 1  detection
         //    TabSupport(s,i,j) = 2  detection by dilation
         //    TabSupport(s,i,j) = 10 detection but isolated pixel
         //    TabSupport(s,i,j) = 11 detection but not enough events

         float *TabLevel;          // Noise standard deviation array
	                           // for the coefficients
         type_noise TypeNoise;     // Kind of noise in the data
         type_noise NewStatNoise;  // Kind of noise in the data after transform
         set_transform Set_Transform;
         type_transform Transform;
         type_border Border;
	 // internal routine: return the band and the position of a coefficient
	 void pos_mrcoeff(int NumCoef, int &b, int &i, int &j);
         // internal routine: return the index in TabLevel and TabSupport
         int index(int b, int i, int j)
	          const {return (TabPos(b)+TabNc(b)*i + j);}
         // internal routine: initialize TabLevel in case of
         //                   non white gaussian noise
         void max_sigma_in_scale(Ifloat &Ima_Sigma);
         // internal routine: modified to calculate the true sigma at 
         //                   each scale from the RMS map
	 void get_sigma_from_rms_map(Ifloat &Ima_Sigma);
         // internal routine: controle if a coefficient value 
         //                   should be thresholded
	Bool kill_coef(int b,int i,int j, float Val, Bool SetSupport);
         // internal routine: set the number of events used per
         //                   wavelet coefficient calcalution 
         // void set_mr_tabevent(Ifloat & ImaEvent);
         void init_param();
	 // internal routine: set the support in case of
	 // speckle noise
	 void speckle_support(MultiResol &MR_Data);
	 // internal routine: set the support in case of
	 // corellated noise
	 void correl_support(MultiResol &MR_Data);

         FilterAnaSynt *FilterBank; // pointer to the filter bank to use
                                    // in the orthogonal wavelet transform   
         int NbrUndecimatedScale; // Number of undecimated scale used
                                  // in the undecimated wavelet transform
                                  // By default (when NbrUnDecimatedScale=-1), 
                                  // all scale are undecimated.
                                  // if NbrUnDecimatedScale = 0.
                                  // then it is an Orthogonal Wavelet Transform
         float multi_sure_estimation(MultiResol & MR_Data, int b);
         float sure_estimation(MultiResol &MR_Data);
	 

	 Bool mOldPoisson;
	 Bool mWriteThreshold;
	 
// **********************  Public part *************************************   
       public:
         type_threshold TypeThreshold; 
         // type of threshold to be used for the detection of the significant coefficents

         type_undec_filter U_Filter;    // Filter used in a non orthogonal 3 directional
	                                // and undecimated WT
					
         // return a pointer to the filter bank
         sb_type_norm TypeNorm;         // type of normalization
	                                // (in L1 or L2)
         FilterAnaSynt * filter_bank() {return FilterBank;}

         int nbr_scale() { return NbrScale; }
         int nbr_band() { return NbrBand; }
         int nbr_undec_scale() {return NbrUndecimatedScale;}

         int nbr_coeff_in_band(int b) { return TabNl(b)*TabNc(b);}
         int nl() { return Nl; }
         int nc() { return Nc; }
         type_noise which_noise() { return TypeNoise; }
         type_transform type_trans() { return Transform; }
 
         Bool TransImag;  // True is the image is transformed before
                          // computing its multiresolution transform
         int NiterSigmaClip;
         int SizeBlockSigmaNoise;	 
	 
	 Bool NoCompDistrib;           
	 
	 Iint EventImage;    // Original signal: used for poisson event
                                  // only.  	 
	 
         MRNoiseModel();
         void alloc(type_noise TNoise, int Nl_Imag, 
                    int Nc_Imag, int ScaleNumber, type_transform Trans, 
                    FilterAnaSynt *FAS=NULL, sb_type_norm Norm=NORM_L1, 
                    int NbrUndec=-1, int FCT_NDir=DEF_FCT_NDIR);
                          // initialize the Class
         void free();
         MRNoiseModel(type_noise TNoise, int Nl_Imag, 
                      int Nc_Imag, int ScaleNumber, type_transform Trans);
                          // initialize the Class
         MRNoiseModel(type_noise TNoise, MultiResol &MR_Data);
                           // initialize the Class
                           // same as before but took the parameter in MR_Data
         void model(Ifloat & Imag);     
                    // set TabLevel and TabSupport
         void model(Ifloat & Imag, MultiResol &MR_Data);
                    //  compute the multiresolution transform, and
                    //  modelize the noise.
                    // set TabLevel and TabSupport 
         void write_support_mr(char *FileName);
         void write_support_ima(char *FileName);
         void mr_obj(MultiResol &MR_Data);
         Bool operator() (int s,int i,int j, details which_detail);
         Bool operator() (int b,int i,int j);
         Bool operator() (int NumCoef);
         float & sigma(int s,int i,int j, details which_detail);
         float & sigma(int b,int i,int j);
         float & sigma(int NumCoef);
	 
	 void write_in_few_event( Bool Write ) ;
	 
 	 float nsigma(int b);
         // return the NSigma value in band b
	 unsigned char & support(int b,int i,int j);
         unsigned char & support(int s,int i,int j,details which_detail);
         unsigned char & support(int NumCoef);
         Bool signif (float Val, int b, int i, int j, 
                      float LevelMin, float LevelMax);
         Bool signif (float Val,int b,int i,int j);
	 Bool signif (float Val,int b,int i,int j, fltarray &TabNsigma);
         Bool signif (float Val,int s,int i,int j, details which_detail);
         float prob(float Val,int s,int i,int j,details which_detail);
              // return the probility that a noise produce the 
              // wavelet coefficient Val at scale s and at position i,j
         float prob(float Val,int b,int i,int j);
         void prob (MultiResol &MR_Data, Bool Complement=False);
              // return the probility density that a noise produce the 
              // wavelet coefficients of MR_Data
              // Results are stored in MR_Data
              // if Complement == True, then it is the complement to 1
              // which is calculated.
         void prob_noise (MultiResol &MR_Data, Bool Complement=False);
	      // return the probility that a noise produce the 
              // wavelet coefficients of MR_Data
              // Results are stored in MR_Data
              // if Complement == True, then it is the complement to 1
              // which is calculated.
         float prob_noise(float Val, int b, int i, int j);
         float prob_noise(float Val, int s, int i,int j, details which_detail);
              // return the integrated probility (repartition function)
              // that a noise produce the wavelet coefficient
              // Val at scale s and at position i,j
         float prob_signal(float Val, int b, int i, int j);
         float prob_signal(float Val, int s, int i,int j, details which_detail);
              // return the integrated probility (repartition function)
              // that the signal produce the wavelet coefficient
              // Val at scale s and at position i,j
              // prob_signal = 1. - prob_noise
         void im_transform (Ifloat &Image); // apply a transform to an image
         void im_invtransform (Ifloat &Image); 
                                        // apply a inverse transform to an image
         float val_transform(float Val); // apply the transform to a value
         float val_invtransform(float Val); // apply the inverse transform to
                                        //  a value
         void kill_isol(int s);
                    // set to 10 the pixel in the support which are
                    //  detected and isolated
         void dilate_support(int s);
                    //  dilate a scale of the support: dilated values are
                    //  set to 2
         void kill_event(int s, Ifloat &Event_Image, int FirstWin=5);
                    //  suppress detection when the number of events in 
                    //  box are less than MinEventNumber
         void set_support(MultiResol &MR_Data);
                    //  set the support using the noise modelization
                    //  and the multiresolution data
	 void mod_support (MultiResol &MR_Data, fltarray Nsigma);
         void refresh_support(MultiResol &MR_Data);
         void set_sigma(Ifloat  &Image, MultiResol &MR_Data);
         void threshold(MultiResol &MR_Data, Bool SetSupport=False);
                    //  threshold a multiresolution tranform. If SetSupport
                    //  equal to True, then the support is updated
         void weight_snr(MultiResol &MR_Data, Bool SetSupport=False);
         void weight_invsnr(MultiResol &MR_Data, Bool SetSupport=False);
	
	 double prob_signal_few_event( float Val, int b, int i, int j ) ;
	 
         Bool DilateSupport;   // if true, then dilate the support
         Bool GetEgde;         // if true, try to better detect edges
	                       // by applying a 1D wavelet transform on 
			       // columns of the horizontal band and on
			       // lines if the vertical band
			       // Do not use it. Results are only poorly
			       // improved.
         Bool SupIsol;         // if true, then supress isolated pixel
         Bool MinEvent;        // if true, then detection with a minimum 
                               //               number of events
         Bool OnlyPositivDetect; // Detect only positive strucutures
         //StatEventPoisson *CEventPois;
	 FewEventPoisson* CFewEventPoisson2d;
	 FewEvent* CFewEvent2d;
         int MinEventNumber;   // minimum number of event for a detection
         Iint Event_Image;         // Count Image: 
                                   // if TypeNoise != NOISE_EVENT_POISSON
                                   // and MinEvent == True, then Event_Image
                                   // must be initialized before calling 
                                   // the model routine.
                                   // if TypeNoise == NOISE_EVENT_POISSON
                                   // Event_Image is automatically calculated
                                   // only if TransImag == True
         Bool UseRmsMap;           // if UseRmsMap, the RMS map is used
                                   // else it is automatically estimated.
                                   // only used for noise = NOISE_NON_UNI_ADD
                                   // i.e. non uniform additive noise
         Ifloat RmsMap;            // RMS map corresponding to the data
                                   // only used for noise = NOISE_NON_UNI_ADD
                                   // i.e. non uniform additive noise
	 Bool GetRmsCoeffGauss;    // the sigma of a wavelet coefficient is
	                           // calculated assuming independant Gaussian
				   // variable (pixels). Default is True.
				   // only used for noise = NOISE_NON_UNI_ADD
				   // and for the a trous algorithm
				   // with a B3 spline (TO_PAVE_BSPLINE)		   
         MultiResol MR_Data_Event; // number of photons used for the calculation
                                   // of each wavelet coefficients.
         Bool BadPixel;        //  if true, then take into account bad pixels

         float BadPixelVal;    //  bad pixel value
         int FirstDectectScale; // first scale use for the detection
         type_sigma_method SigmaDetectionMethod; 
                                 // method for sigma noise detection
         float CCD_Gain;         // CCD gain 
         float CCD_ReadOutSigma; // CCD read-out noise standard deviation 
         float CCD_ReadOutMean;  // CCD read-out noise mean 
         float SigmaNoise;       // noise standard deviation (Gaussian Noise) 
         float NSigma[MAX_BAND];           // Detection at NSigma * Noise 
         float TabEps[MAX_BAND]; // Epsilon detection level.
                                  // used for poisson noise only
	 StatNoiseMap *CorrelNoiseMap; // class used for correlated noise			  
         StatRayleigh *CSpeckle; // class used for speckle noise
	 Bool SigmaApprox;       // approximate the noise law to a gaussian law
                                 // only used for poisson with few event and
				 // rayley noise
         Bool GradientAnalysis;  // Only used with Transform == TO_DIADIC_MALLAT
                                 // If true, the detection is not done from
                                 // the wavelet coefficient, but from the
                                 // gradients: g_x,y = sqrt( w_x^2 + w_y^2 )
         Bool PoissonFisz;       // if True, then the Anscombe transform is
	                         // replaced by the Fisz transform.
 	 Bool MadEstimFromCenter; // If true, the mad estimation is performed only in the image center
	 void hierarchical_dilate_support();
         void set_old_poisson( Bool Flag );
	 void write_threshold( Bool Flag );
	 Bool old_poisson();
	 void trace();
         ~MRNoiseModel();
};

// void find_typenoise(Ifloat & Imag);
Bool one_level_per_pos_2d(type_noise TNoise); // True if one sigma for each (i,j)


#endif
