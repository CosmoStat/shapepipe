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
**    Date:  5/03/2000 
**    
**    File:  SB_Filter1D.h
**
**    Modification history :
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

#ifndef _SB_FILTER1D_H_
#define _SB_FILTER1D_H_

#include "Border.h"
#include "OptMedian.h"
#include "Filter.h"
#include "Atrou1D.h"

/***********************************************************************/

inline int size_resol(int Scale, int N)
{
       int Val = N;
         for (int s=0; s < Scale; s++) Val = (Val+1) / 2;
         return Val;
}

/***********************************************************************/

// Class for sub-band decomposition with and without decimation
class SubBand1D {
         public:
	   Bool SubSample_H_Even;
           Bool SubSample_G_Odd;
	   int DistPix; // distance between two adjacent pixels 
	                // normally defaulted to 1.
			// only used with the decimated transform
	   type_border Border;    // Border management type
           SubBand1D (){Border=I_MIRROR;DistPix=1;
	                SubSample_H_Even = True;SubSample_G_Odd = True;}
           inline int test_index(int i, int N)
           {
              // if (Border == I_MIRROR) return test_index_mirror(i,N);
              // else return test_index_cont(i,N);
	      return get_index(i, N, Border);
           }           
	   
	   virtual void transform (int , float *, float *, float *){};
	   //                  (int N, float *High, float *Low, float *Det)
	   // sub-band decomposiiton with decimation
	   //  N = number of pixels in the input signal High
	   // High is decomposed in two sub signals Low and Det
	   // Low = low resolution signal: its size is (N+1)/2
	   // Det = Detail information:  its size is N /2 
	   
           virtual void recons (int , float *,  float *,  float *){};
	   //                  (int N, float *Low, float *Det, float *High)
	   // sub-band reconstruction with decimation
	   // inverse transform: the output High (size N) is reconstructed
	   // from Low (size (N+1)/2) and Det (size N/2)
	   
	   
	   virtual void noise_transform (int , float *, float *, float *){};
	   //                  (int N, float *High, float *Low, float *Det)	   
           
 	   
	   virtual void transform (int , float *, float *, float *, int){};
	   //               (int N, float *High, float *Low, float *Det, int Step)
	   // sub-band decomposiiton without decimation
	   //  N = number of pixels in the input signal High
	   // High is decomposed in two sub signals Low and Det
	   // Low = low resolution signal: its size is N
	   // Det = Detail information:  its size is N
	   
           virtual void recons (int , float *,  float *,  float *, int){};
	   //                  (int N, float *Low, float *Det, float *High, int Step)
	   // sub-band reconstruction without decimation
	   // inverse transform: the output High (size N) is reconstructed
	   // from Low (size N) and Det (size N)
	   	   
	   virtual ~SubBand1D(){};
};

/***********************************************************************
* filter bank transform (one step transform) using (bi-) orthogonal
* WT. It can be decimated or an undecimated decomposition
************************************************************************/

class SubBandFilter: public SubBand1D
{
   int NormCoef; // Norm = 2 for L1 normalization 
                  // Norm = 1 for L2 normalization
   type_sb_filter TypeFilter; // type of bands
   float *H0,*G0,*H1,*G1; // filter banks
   int Size_H0, Size_H1, Size_G0, Size_G1; // filter sizes
   int Start_H0, Start_H1, Start_G0, Start_G1; // First position of the 
                                               // filter
   sb_type_norm TypeNorm;
   void reset_param();
   void init(FilterAnaSynt & FILTER, sb_type_norm Norm);
   void init(float *H, int Nh, float *G, int Ng, sb_type_norm FilNorm);
   void init();
    Bool Verbose;
   public:

    SubBandFilter(type_sb_filter T_Filter);
    SubBandFilter(type_sb_filter T_Filter, sb_type_norm Norm);
    SubBandFilter(FilterAnaSynt &FILTER);
    SubBandFilter(FilterAnaSynt &FILTER, sb_type_norm Norm);
    SubBandFilter(char *FileName);
    SubBandFilter(char *FileName, sb_type_norm Norm);
    SubBandFilter(float *H, int Nh, float *G, int Ng, sb_type_norm FilNorm);

    // Decimated SubBand decomposition
    void convol_h0 (int N, float *Input, float *Output);
    // convolved a signal by H0 filter
    void convol_g0 (int N, float *Input, float *Output); 
    // convolved a signal by G0 filter
    void noise_convol_h0 (int N, float *Input, float *Output);
    // convolved a signal by H0^2 filter
    void noise_convol_g0 (int N, float *Input, float *Output);
    // convolved a signal by G0^2 filter
    void convol_h1 (int N, float *Input, float *Output);
    // convolved a signal by H1 filter
    void convol_g1 (int N, float *Input, float *Output);
    // convolved a signal by G1 filter

    void transform (int N, float *SignalIn, float *SignalOut, float *DetailOut);
    // decimated one step transform
    // SignalIn = input signal of size N
    // SignalOut = output smooth signal of size (N+1)/2
    // DetailOut = output wavelet coefficient of size N/2
    // the DistPix field is used in order to know the distance to use
    // between two consecutive pixels. i.e if a given scale 
    // has been undecimated, the pixel to consider must separated by a
    // a distance of 2.

    void recons (int N, float *SignalIn,  float *DetailIn,  float *SignalOut);
    // reconstruct a signal for the low resolution part and the 
    // associated wavelet coefficient
    // N = ouput signal size
    // SignalIn = input signal of size(N+1)/2
    // DetailIn = input wavelet coefficient of size N/2
    // SignalOut = outpout signal of size N

    void noise_transform (int N, float *SignalIn, float *SignalOut, float *DetailOut);
    // SignalIn = input signal of size N
    // SignalOut = output signal of size (N+1)/2 = SignalIn convolved with H0^2
    // DetailOut = output signal coefficient of size N/2 = SignalIn convolved with G0^2

    // Undecimated SubBand decomposition
    void convol_filter(int N, float *Input, float *Output, 
                      float *F, int SizeFilter, int Start_Filter, int Step);
    // convolution routine used for the transformation
    // Input and Output signals have the same size N
    // F = pointer to filter value
    // SizeFilter = filter size
    // Start_Filter = low range index of the filter
    // Step = distance between two cofficients ("trous" size).

    void rec_convol_filter(int N, float *Input, float *Output, 
                       float *F, int SizeFilter, int Start_Filter, int Step);  
    // convolution routine used for the reconstruction
    // Input and Output signals have the same size N
    // F = pointer to filter value
    // SizeFilter = filter size
    // Start_Filter = low range index of the filter
    // Step = distance between two cofficients ("trous" size).

    void transform (int N, float *SignalIn, float *SignalOut, float *DetailOut, int Step);
    // one step undecimated subband transform
    // SignalIn = input signal of size N
    // SignalOut= output smooth signal of size N
    // DetailOut= output wavelet coefficient signal of size N
    // Step = distance between two cofficients ("trous" size).

    void recons (int N, float *SignalIn,  float *DetailIn,  float *SignalOut, int Step);
    // reconstruct a signal for the low undecimated resolution part and the 
    // associated unbdecimated wavelet coefficient
    // N = input-ouput signal size
    // SignalIn = input signal of size N
    // DetailIn = input wavelet coefficient N
    // SignalOut = outpout signal of size N
   
    void transform (fltarray& Sig_in, fltarray& Sig_out, fltarray* Smooth=NULL);
    //decimated subband transform
    // Sig_in = input signal of size N
    // Sig_out= output detail signal of size N
    // Smooth= output smooth signal of size N
   
    void recons (fltarray& Sig_in, fltarray& Sig_out, fltarray& Smooth);
    // reconstruct a signal for the low decimated resolution part and the 
    // associated decimated wavelet coefficient
    // Sig_in = input signal of size N
    // Sig_out = outpout signal of size N
    // Smooth = input smooth wavelet coefficient N
                                
   ~SubBandFilter();
};


/***********************************************************************
* G min-max transform (one step transform)
************************************************************************/

class G_Minmax: public SubBand1D {
  public:
  void transform(int N, float *High, float *Low, float *Det);
  void recons(int N, float *Low, float *Det, float *High);
  
  void transform(int N, float *High, float *Low, float *Det, int Step);
  void recons(int N, float *Low, float *Det, float *High, int Step);
  G_Minmax(){Border = I_CONT;}
  ~G_Minmax(){}
};

/***********************************************************************
*  1D lifting scheme (one step transform)
************************************************************************/

#define NBR_LIFT 7
enum type_lift {TL_UNKNOWN, TL_CDF, TL_MEDIAN, TL_INT_HAAR, TL_INT_CDF, TL_INT_4_2, TL_F79, TL_INT_F79};
#define DEF_LIFT TL_INT_HAAR

const char * StringLSTransform (type_lift  type);
void lifting_usage(type_lift Transform);

class Lifting: public SubBand1D {

   float lift_predict(int i, int N, float *High, int Step=1);
   float lift_update(int i, int N, float *Det, int Step=1);
   void transform_f79(int N, float *High, float *Low, float *Det);
   void transform_int_f79(int N, float *High, float *Low, float *Det);
   void recons_f79(int N, float *Low, float *Det, float *High);
   void recons_int_f79(int N, float *Low, float *Det, float *High);
  public:
  void transform(int N, float *High, float *Low, float *Det);
  void recons(int N, float *Low, float *Det, float *High);
  void transform(int N, float *High, float *Low, float *Det, int Step);
  void recons(int N, float *Low, float *Det, float *High, int Step);
  
  type_lift TypeTrans;
  Lifting(){TypeTrans=DEF_LIFT;}
  Lifting(type_lift TL) {TypeTrans=TL;}
  ~Lifting(){}
};

/***********************************************************************
* Undecimated non-orthgonal filter bank (one step transform)
************************************************************************/

// U_B3SPLINE: H = B3spline, G=Id-H, Ht=H, Gt=Id+H
// U_B3SPLINE_2: H = B3spline, G=Id-H, Gt=G, HT =2*Id-H 

#define NBR_UNDEC_FILTER 4
enum type_undec_filter {U_B3SPLINE,U_B3SPLINE_2,U_B2SPLINE,U_HAAR_B3S_POS,U_HAAR_B3S};
#define DEF_UNDER_FILTER U_B3SPLINE_2

const char * StringUndecFilter (type_undec_filter TypeUnderFilter);
void usb_usage(type_undec_filter Filter);

class UndecSubBandFilter: public SubBand1D
{
   float *FilterH;
   float *FilterG;
   float *FilterTH;
   float *FilterTG;
   int SizeFilterH;
   int SizeFilterG;
   int SizeFilterTH;
   int SizeFilterTG;
   // sb_type_norm TypeNorm;
   type_undec_filter TypeUnderFilter;
   void init();
   void convol_filter(int N, float *Input, float *Output, float *F, int SizeFilter,  int Step);
public:
    Bool Verbose;
    UndecSubBandFilter(){TypeUnderFilter=U_B2SPLINE;init();}
    UndecSubBandFilter(type_undec_filter Type) {TypeUnderFilter=Type;init();}
    // UndecSubBandFilter(){TypeUnderFilter=U_B2SPLINE;TypeNorm=NORM_L1;init();}
    // UndecSubBandFilter(type_undec_filter Type, sb_type_norm TN)  {TypeUnderFilter=Type;TypeNorm= TN;init();}

    void transform(int N, float *High, float *Low, float *Det, int Step);
    // sub-band decomposition without decimation
    // using non-orthogonal filter banks
    //  N = number of pixels in the input signal High
    // High is decomposed in two sub signals Low and Det
    // Low = low resolution signal: its size is N
    // Det = Detail information:  its size is N

    void recons(int N, float *Low, float *Det, float *High, int Step);
    // sub-band reconstruction without decimation
    // inverse transform: the output High (size N) is reconstructed
    // from Low (size N) and Det (size N)

    ~UndecSubBandFilter() { if (FilterH != NULL) delete [] FilterH; FilterH=NULL; 
                            if (FilterG != NULL) delete [] FilterG; FilterG=NULL;
			    if (FilterTG != NULL) delete [] FilterTG; FilterTG=NULL;
                            if ((TypeUnderFilter == U_HAAR_B3S) && (FilterTH  != NULL)) delete [] FilterTH;
			    FilterTH=NULL;
 			  }
};

/***********************************************************************
* (bi-) Othogonal wavelet wavelet transform
***********************************************************************/

class Ortho_1D_WT {
       SubBand1D *Ptr_SB1D;
     public:
        Ortho_1D_WT(SubBand1D &SB1D) {Ptr_SB1D = &SB1D;};
	void transform (fltarray &SignalIn, fltarray &Transf_out, int Nbr_Plan);
        void recons(fltarray &Transf_in, fltarray &Sig_Out, int Nbr_Plan);
        ~Ortho_1D_WT(){Ptr_SB1D = NULL;}
};

/***********************************************************************
* Half undecimated  wavelet wavelet transform
***********************************************************************/

class HALF_1D_WT {
       SubBand1D *Ptr_SB1D;
       int NbrUndecimatedScale;
     public:
        HALF_1D_WT() {Ptr_SB1D = NULL; NbrUndecimatedScale=0;}
        HALF_1D_WT(SubBand1D &SB1D, int NbrUndec) 
	      {Ptr_SB1D = &SB1D;NbrUndecimatedScale=NbrUndec;}
	void alloc(SubBand1D &SB1D, int NbrUndec)
	      {Ptr_SB1D = &SB1D;NbrUndecimatedScale=NbrUndec;}
	void transform (fltarray &SignalIn, fltarray &Transf_out, int Nbr_Plan);
        void recons(fltarray &Transf_in, fltarray &Sig_Out, int Nbr_Plan);
        ~HALF_1D_WT(){Ptr_SB1D = NULL;}
};

/***********************************************************************
* undecimated  wavelet wavelet transform
***********************************************************************/

class PAVE_1D_WT {
       SubBand1D *Ptr_SB1D;
     public:
        PAVE_1D_WT(SubBand1D &SB1D) {Ptr_SB1D = &SB1D;};
        float getAbsMaxTransf (fltarray & WT_Trans, float SigmaNoise, Bool OnlyPositivDetect=False, Bool Verbose=False); 
	void transform (fltarray &SignalIn, fltarray &Transf_out, int Nbr_Plan);
	void transform (fltarray &SignalIn, fltarray &Transf_out, 
                        fltarray &Smooth, int step);
        void recons (fltarray &Transf_in, fltarray &Sig_Out, int Nbr_Plan);
        void recons (fltarray &Transf_in, fltarray &Sig_Out, 
                     fltarray &Smooth, int step);
        ~PAVE_1D_WT(){Ptr_SB1D = NULL;}
};


/***********************************************************************
*  Orthogonal decimated-Undecimated 1D wavelet transform
************************************************************************/

class HALF_DECIMATED_1D_WT {
private:
      SubBandFilter *Ptr_SB1D;
      
      void  set_tabdec(int NumUndec, Bool * & TabDec, int Nbr_Plan);
      void transform (fltarray& Signal, fltarray* TabTrans,  Bool *TabDec);
      float update (float CoefSol, float Threshold, float SoftLevel);
      bool is_on_support (float CoefSol, float Threshold, float SoftLevel);
      void recons (fltarray* TabTrans, fltarray& Signal, Bool *TabDec);
                       
public:
      int NbScale;
      int NbrUndecimatedScale;
      type_border Bord;      
      float SigmaNoise;
      bool OnlyPositivDetect;
      bool SuppressIsolatedPixel;
      int FirstDetectScale;
      bool Verbose;
      bool Write;
      bool UseNormL1;
      float NSigmaSoft;      
            
      HALF_DECIMATED_1D_WT (SubBandFilter &SB1D) {Ptr_SB1D = &SB1D;}
      int alloc (fltarray* & TabTrans, int Nx, int Nbr_Plan, Bool *TabDec);
      int alloc (fltarray* & TabTrans, int Nx, int Nbr_Plan, int NbrUndecimatedScale);
      void free (fltarray* TabTrans, int Nbr_Plan);
      
      void transform (fltarray& Signal, fltarray* TabTrans);
      void KillScaleNotUsed (fltarray* WT_Trans, int FirstDetectScale);
      void KillLastScale (fltarray* WT_Trans);
      void threshold (fltarray* WT_Trans, float NSigma, int IterNumber);
      float getAbsMaxTransf (fltarray* WT_Trans);
      void recons (fltarray* TabTrans, fltarray& Signal);
      ~HALF_DECIMATED_1D_WT() {Ptr_SB1D = NULL;}
};

#endif
