/******************************************************************************
**                   Copyright (C) 1996 by CEA
*******************************************************************************
**
**    UNIT
**
**    Version: 0.1
**
**    Author: Rene Gastaud & Jean-Luc Starck
**
**    Date:  22-OCT-1996
**    
**    File:  MR1D_NoiseModel.h
**
*******************************************************************************
**
**    DESCRIPTION  Noise Modelling Class  1 dimension
**    ----------- 
**                 
**        MR1DNoiseModel NoiseModel(TNoise, MR_Data);
**     or MR1DNoiseModel NoiseModel(TNoise, N, ScaleNumber, Transform)
**        allocates space for Noise Modelling
**
**     NoiseModel.model(Imag); 
**        Estimates the detection level at each scale and creates the
**        MR_1Dution support
**     
**     NoiseModel.model(MR_Data);
**        Estimates the detection level at each scale and creates the
**        MR_1Dution support. MR_Data is the MR_1Dution of the
**        image already transformed.
**        
**     NoiseModel.support(s,i) return True is something has been detected
**        at this scale and at this position
**
**     NoiseModel.signif(Val, s, i) return True is this value is upper
**        the detection level at this scale and at this position.
**
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

#ifndef _C1DNOISEMODEL_H_
#define _C1DNOISEMODEL_H_
#include "MR1D_Obj.h"    
#include "MR1D_Sigma.h"
#include "Mr_FewEvent.h"  
#include "MR_Threshold.h"  
#include "IM_Noise.h"  

 
#define DEF_MINEVENT_NUMBER 4 
#define DEF_NSIGMA_DETECT 3.
#define DEF_FIRST_DETECT_SCALE 0

#define VAL_SupNull 0
#define VAL_SupOK 1
#define VAL_SupDill 2
#define VAL_SupLastOK 9
#define VAL_SupKill 10
#define VAL_SupMinEv 11
#define VAL_SupFirstScale 12



class MR1DNoiseModel {
         int NbrScale;  	// number of scales in the wavelet transform
         int Nelem ;       		// size of the input signal
         int TabN[MAX_SCALE_1D]; 	// for pyramidal transform
         int TabPos[MAX_SCALE_1D]; // for pyramidal transform
         int Size;		// size of the transformed signal,
				// size = NbrScale*Nelem if non pyramidal
         unsigned char *TabSupport;
         // multiresolution support of a signal
         //    TabSupport(s,i,j) = 0  no detection
         //    TabSupport(s,i,j) = 1  detection
         //    TabSupport(s,i,j) = 2  detection by dilation
         //    TabSupport(s,i,j) = 10 detection but isolated pixel
         //    TabSupport(s,i,j) = 11 detection but not enough events

         float *TabLevel;	// pointer on the cube of sigma
         type_noise TypeNoise;     // Kind of noise in the data
         type_noise NewStatNoise;  // Kind of noise in the data after transform
         type_trans_1d Type_Transform; // see MR1D_Obj.h
          // enum type_trans_1d { TO1_PAVE_LINEAR, ...,TO1_PYR_LINEAR,...}
    
         void init_param(); // initialization      
                                                            
         set_trans_1d  Set_Transform; // see MR1D_Obj.h
         // enum set_trans_1d {TRANS1_PAVE, TRANS1_PYR, S1_UNDEFINED=-1};

         // internal routine: return the index in TabLevel and TabSupport
         int index(int s, int i) const ;

         // internal routine: initialize TabLevel in case of
         //                   non white gaussian noise
         void max_sigma_in_scale(fltarray & Signal_Sigma);

         // internal routine: controle if a coefficient value 
         //                   should be thresholded
        Bool kill_coef(int s, int i, float Val, Bool SetSupport);

         FilterAnaSynt *FilterBank; // pointer to the filter bank to use
                                    // in the orthogonal wavelet transform   

        void get_sigma_from_rms_map(fltarray  &Dat_Sigma);
        float multi_sure_estimation(MR_1D & MR_Data, int b);
	float sure_estimation(MR_1D &MR_Data);

      public:
         type_threshold TypeThreshold; 
         // type of threshold to be used for the detection of the significant coefficents
 
	 float prob_noise(float Val, int b, int i);
         // return the probility that a noise produce the 
         // wavelet coefficients Val at scale b and position i
	 
         // return a pointer to the filter bank
         sb_type_norm TypeNorm;         // type of normalization
	                                // (in L1 or L2)
         FilterAnaSynt * filter_bank() {return FilterBank;}

         int nbr_scale() { return NbrScale; }
	 int get_size() {return Size;}
	 int get_index (int s, int i) {return index(s,i);}
         type_noise which_noise() {return TypeNoise;}
         int nelem() { return Nelem; }  // Upper case and lower case nelem Nelem
         type_trans_1d type_trans() { return Type_Transform; }
         set_trans_1d set_trans() { return Set_Transform; } // add by RG
         Bool TransImag;  // True is the image is transformed before
                          // computing its MultiResolution transform
         int NiterSigmaClip;
         int SizeBlockSigmaNoise;
	 
	 FewEventPoisson* CFewEventPoisson1d;
	 
	 void alloc (type_noise TNoise, int NbPoint, 
	             int ScaleNumber, type_trans_1d Trans, 
                     FilterAnaSynt *FAS=NULL, sb_type_norm Norm=NORM_L1);

         // initialize the Class
	 MR1DNoiseModel();
         MR1DNoiseModel(type_noise TNoise, int N_signal, 
                      int ScaleNumber, type_trans_1d Trans);
         MR1DNoiseModel(type_noise TNoise, int N_signal, 
                      int ScaleNumber, type_trans_1d Trans,
                      FilterAnaSynt *FAS, sb_type_norm Norm);       
         MR1DNoiseModel(type_noise TNoise, MR_1D &MR_Data);
                           // initialize the Class
                           // same as before but took the parameter in MR_Data
         
         void model(fltarray & Signal);     
                    // set TabLevel and TabSupport
         void model(fltarray &Signal, MR_1D &MR_Data);
                    //  compute the MultiResolution transform, and
                    //  modelize the noise.
                    // set TabLevel and TabSupport 
         void write_support_mr(char *FileName);
         void write_support_ima(char *FileName);
         void mr_obj(MR_1D &MR_Data);
         Bool operator() (int s,int i) const;
         Bool operator() (int NumCoef) const;
         float & sigma(int s, int i) const;
         float & sigma(int NumCoef) const;
         unsigned char & support(int s,int i) const;
         unsigned char & support(int NumCoef) const;

         Bool signif (float Val, int s, int i);         
	 Bool signif (float Val, int s, int i, fltarray& Nsigma);
         void signal_transform (fltarray & Signal); 
                   	// apply a transform to a signal
         void signal_invtransform (fltarray & Signal); 
                        // apply a inverse transform to a signal
         float val_transform(float Val);
                   	// apply a transform to a value
         float val_invtransform(float Val);
                        // apply a inverse transform to a value
         void kill_isol(int s);
                    // set to 10 the pixel in the support which are
                    //  detected and isolated
         void dilate_support(int s);
                    // dilate the multiresolution support
         void kill_event(int s, fltarray &Event_Signal, int FirstWin=5);
                    //  suppress detection when the number of events in 
                    //  box are less than MinEventNumber
         void set_support(MR_1D &MR_Data);
                    //  set the support using the noise modelization
                    //  and the MultiResolution data
         void refresh_support(MR_1D &MR_Data);
         void set_sigma(fltarray & Signal, MR_1D &MR_Data);
         void threshold(MR_1D &MR_Data, Bool Smooth=False, Bool SetSupport=False);
                    //  threshold a MultiResolution tranform. If SetSupport
                    //  equal to True, then the support is updated
         void weight_snr(MR_1D &MR_Data, Bool SetSupport=False);
         void weight_invsnr(MR_1D &MR_Data, Bool SetSupport=False);
         Bool SupIsol;         // if true, then supress isolated pixel
         Bool MinEvent;        // if true, then detection with a minimum 
                               //               number of events
         Bool OnlyPositivDetect; // Detect only positive strucutures
         Bool OnlyNegativDetect; // Detect only negative strucutures
         int MinEventNumber;   // minimum number of event for a detection
	 intarray EventSignal;    // Original signal: used for poisson event
                                  // only.         
	 //Iint Event_Image;         // Original Image: used for poisson event
                                     // only.
	 //intarray Event_Image;         // Original Image: used for poisson event
                                     // only.			     
         MR_1D MR_Data_Event; // number of photons used for the calculation
                                   // of each wavelet coefficients.
         Bool BadPixel;        //  if true, then take into account bad pixels
         Bool DilateSupport;   //  if true, then the support is dilated
         Bool UseRmsMap;           // if UseRmsMap, the RMS map is used
                                   // else it is automatically estimated.
                                   // only used for noise = NOISE_NON_UNI_ADD
                                   // i.e. non uniform additive noise
         fltarray RmsMap;            // RMS map corresponding to the data
                                   // only used for noise = NOISE_NON_UNI_ADD
                                   // i.e. non uniform additive noise
         float BadPixelVal;    //  bad pixel value
         int FirstDectectScale; // first scale use for the detection
         type_sigma_method_1d SigmaDetectionMethod; 
                                 // method for sigma noise detection
         float CCD_Gain;         // CCD gain 
         float CCD_ReadOutSigma; // CCD read-out noise standard deviation 
         float CCD_ReadOutMean;  // CCD read-out noise mean 
         float SigmaNoise;       // noise standard deviation (Gaussian Noise) 
         float NSigma[MAX_SCALE_1D];           // Detection at NSigma * Noise 
         float TabEps[MAX_SCALE_1D]; // Epsilon detection level.
                                     // used for poisson noise only

         void signal_sigma(fltarray &Data, fltarray &Data_Sigma, 
                    int SizeBlock, int Nit); // RG from im_sigma in im_noise.cc
                  //  Data  is an input, I add const for using toto.scale(s)
         ~MR1DNoiseModel();
};

#include "Mr_FewEvent1d.h"


#endif
