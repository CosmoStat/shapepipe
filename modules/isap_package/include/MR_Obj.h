/******************************************************************************
**                   Copyright (C) 1994 by CEA
*******************************************************************************
**
**    UNIT
**
**    Version: 3.1
**
**    Author: Jean-Luc Starck
**
**    Date:  96/05/02 
**    
**    File:  MR_Obj.h
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

#ifndef _CMULTIRESOL_H_
#define _CMULTIRESOL_H_

#include "IM_Obj.h"
#include "SB_Filter.h"
#include "IM_IO.h"
#include "MeyerWT.h"
#include "LineCol.h"
#include "FCur.h"

#define MAX_BAND 200
#define MAX_SCALE 10
#define MAX_SIZE_LAST_SCALE 4
#define DEF_FCT_NDIR 16 
#define NBR_SET_TRANSFORM 6

enum set_transform {TRANSF_PAVE, TRANSF_PYR,
                    TRANSF_SEMIPYR,
                    TRANSF_MALLAT, 
		    TRANSF_DIADIC_MALLAT,
		    TRANSF_UNDECIMATED_MALLAT,
                    TRANSF_FEAUVEAU,
                   S_UNDEFINED=-1};
const char * StringSetTransform (set_transform type);



#define NBR_TRANSFORM 28 // transform with a noise modelling
#define NBR_TOT_TRANSFORM 31 // transform with or without a noise modelling

// #define NBR_ISOTROP_TRANSFORM 15
#define NBR_ISOTROP_TRANSFORM 13

enum type_transform { TO_PAVE_LINEAR, 
                     TO_PAVE_BSPLINE, 
                     TO_PAVE_FFT,
                     TM_PAVE_MEDIAN,
                     TM_PAVE_MINMAX, 
                     TO_PYR_LINEAR, 
                     TO_PYR_BSPLINE,
                     TO_PYR_FFT_DIFF_RESOL, 
                     TO_PYR_MEYER,
                     TM_PYR_MEDIAN,
                     TM_PYR_LAPLACIAN, 
                     TM_PYR_MINMAX,
                     TM_PYR_SCALING_FUNCTION,
                     TO_MALLAT,
                     TO_FEAUVEAU, 
                     TO_PAVE_FEAUVEAU,
                     TO_LC,  
                     TO_HAAR,
		     TO_SEMI_PYR,
		     TM_TO_SEMI_PYR,
		     TO_DIADIC_MALLAT, 
		     TM_TO_PYR,
                     TO_PAVE_HAAR,
		     TO_UNDECIMATED_MALLAT,
		     TO_UNDECIMATED_NON_ORTHO,
		     TO_PYR_MEYER_ISOTROP,
                     TO_PYR_FFT_DIFF_SQUARE,
 		     TC_FCT,
		     TO_LIFTING,
		     TO_DIV_1, 
		     TO_DIV_2,
		     TO_DIADIC_HAAR,   // this one does not appear in the list 
		     TM_MIN_MAX,  /* G transform */
		     // TLC_DIADIC_HAAR,
		     // TLC_MALLAT,
		     // TO_MAX_DIADIC_MALLAT,
                     T_UNDEFINED=-1};
		     
const char * StringTransform (type_transform type);


inline void transform_usage(type_transform Transform)
{
    fprintf(OUTMAN, "        [-t type_of_multiresolution_transform]\n");
    for (int i = 0; i < NBR_TRANSFORM; i++)
    fprintf(OUTMAN, "              %d: %s \n",i+1,
                                            StringTransform((type_transform)i));
    fprintf(OUTMAN, "             default is %s\n",  StringTransform((type_transform)Transform));
}

inline void all_transform_usage(type_transform Transform)
{
    fprintf(OUTMAN, "        [-t type_of_multiresolution_transform]\n");
    for (int i = 0; i < NBR_TOT_TRANSFORM; i++)
    fprintf(OUTMAN, "              %d: %s \n",i+1,
                                            StringTransform((type_transform)i));
    fprintf(OUTMAN, "             default is %s\n",  StringTransform((type_transform)Transform));
}

inline void isotrop_transform_usage(type_transform Transform)
{
    fprintf(OUTMAN, "        [-t type_of_multiresolution_transform]\n");
    for (int i = 0; i < NBR_ISOTROP_TRANSFORM; i++)
    fprintf(OUTMAN, "              %d: %s \n",i+1,
                                            StringTransform((type_transform)i));
    fprintf(OUTMAN, "             default is %s\n",  StringTransform((type_transform)Transform));
}

inline Bool isotrop(type_transform type)
{
    switch (type)
    {
        case TO_PAVE_LINEAR: 
        case TO_PAVE_BSPLINE: 
        case TO_PAVE_FFT: 
        case TM_PAVE_MEDIAN: 
        case TM_PAVE_MINMAX:
        case TO_PYR_LINEAR: 
        case TO_PYR_BSPLINE: 
        case TO_PYR_FFT_DIFF_RESOL: 
        case TO_PYR_FFT_DIFF_SQUARE: 
        case TM_PYR_MEDIAN: 
        case TM_PYR_MINMAX:
        case TM_PYR_LAPLACIAN: 
        case TM_PYR_SCALING_FUNCTION: 
        case TO_PAVE_FEAUVEAU: 
	case TO_SEMI_PYR:
	case TM_TO_SEMI_PYR:
	case TM_TO_PYR:
	case TO_PYR_MEYER:
	case TO_PYR_MEYER_ISOTROP:
	// case TO_DIADIC_MALLAT:
               return True;
        default:
              return False;
    }
}

/*************************************************************************/

inline int get_nbr_scale (int N)
{
    int ScaleMax;
    ScaleMax  = iround((float)log((float) (N / 4. * 3.) / log(2.)));
    return ScaleMax;
}

/*************************************************************************/


// image size at resolution s (N = size of the image at resolution 0)
inline int size_ima_resol(int Scale, int N)
{
   int Val = N;
   for (int s=0; s < Scale; s++) Val = (Val+1) / 2;
   return Val;
}

Bool test_scale (int Scale, int Nbr_Plan, char *NameProcCall);
Bool test_type_transform (type_transform Typ, char *NameProcCall) ;
Bool test_set_transform (set_transform Set, char *NameProcCall); 
Bool test_detail (int Scale, int Nbr_Plan, set_transform Set_Transform, 
                  details which_detail);

set_transform SetTransform (type_transform Transform);
int number_band_per_resol(type_transform Transform);

inline void check_scale(int Nl, int Nc, int Np)
{  
   int Nmin = MIN(Nl,Nc);
   int ScaleMax=iround(log((float)Nmin/(float)MAX_SIZE_LAST_SCALE) / log(2.)+ 1.);
   if (Np > ScaleMax)
   {
        fprintf(stdout, "Error: the maximum allowed number of scales for a %dX%d image is %d.\n" , Nl, Nc,ScaleMax);
	exit(-1);
   }  
}

inline int get_correct_scale(int Nl, int Nc, int Np, Bool Mes)
{  
   int Nret = Np;
   int Nmin = MIN(Nl,Nc);
   int ScaleMax=iround(log((float)Nmin/(float)MAX_SIZE_LAST_SCALE) / log(2.)+ 1.);
   if (Np > ScaleMax)
   {
      if (Mes == True)
      {
       fprintf(stdout, "Warning: the maximum allowed number of scales for a %dX%d image is %d, \n" , Nl, Nc,ScaleMax);
       fprintf(stdout, "         and the number of scales is set to this value.\n");
      }
      Nret  = ScaleMax;
   }
   return Nret;  
}




// 2D ORTHOGNAL WT 
//      HORIZONTAL = 1er band = psi(x) phi(y)
//      VERTICAL = 2er band =   phi(x) psi(y)
//      DIAGONAL = 1er band =   psi(x) psi(y)
 
// return the band number of a given scale of a transform 
// with NbrBand scales and of type which_detail
int scale2band(int s, type_transform Transform, 
                  int NbrBand, details which_detail);
void band2scale(const int NumBand, type_transform Type_Transform, 
                   int NbrBand, int & s, details &which_detail);
                  
/* Function definitions for the algorithms using FFT  */
#define SCALING_FUNCTION 1
#define FILTER_H 2
#define FILTER_H_TILDE 3
#define FILTER_G 4
#define FILTER_G_TILDE 5
#define WAVELET 6

#define DEFAULT_N_SIGMA 3.
#define DEFAULT_NBR_SCALE 4
#define DEFAULT_ITER_SCALING 15
#define DEFAULT_TRANSFORM TO_PAVE_BSPLINE
#define DEFAULT_CUT_OFF_FREQ 0.5
#define DEFAULT_ITER_COR_PYR 3
#define DEFAULT_MEDIAN_WINDOW_SIZE 5

class MultiResol {
	 SubBandFilter *SBF_LC;     // SubBand filter for the LC transform
  
         Bool FilterBankAlloc;      // True if the filter bank is allocated
                                    // by the class
         FilterAnaSynt *FilterBank; // pointer to the filter bank to use
                                    // in the orthogonal wavelet transform
         FilterAnaSynt *FilterBank_Line; // pointer to the filter bank to use
                                    // in the orthogonal wavelet transform on line
	 FilterAnaSynt *FilterBank_Column; // pointer to the filter bank to use
                                    // in the orthogonal wavelet transform on column
	 UndecSubBandFilter *UndecFilterBank; // pointer to an undecimated filter bank.
	                           // used in a non orthogonal undecimated tri-directional transform.
	 Ifloat *TabBand;   // array of bands.  
         int Nbr_Plan;      // number of scales
         int Nl, Nc;        // size of the input image
         int NbrBand;       // NumBer of Bands: 
                            //   for isotropic transforms, NbrBand = NBr_Plan
         int NbrBand_per_Resol; // Number of band per resolution 
         int Nbr_MrCoef;    // number of multiresolution coefficients
                            // without taking into account the last scale
         int *TabBandNl;  // Number of lines of each band
         int *TabBandNc;  // Number of columns of each band
         void init_param(); // reset all parameters
         int NbrUndecimatedScale; // Number of undecimated scale used
                                  // in the undecimated wavelet transform
                                  // By default (when NbrUnDecimatedScale=-1), 
                                  // all scale are undecimated.
                                  // if NbrUnDecimatedScale = 0.
                                  // then it is an Orthogonal Wavelet Transform
         
         inline void test_band(int B) const
         {
              if ((B < 0) || (B >= NbrBand))
              {
                 cerr << "error: Bad band number :  " <<  B+1 << endl;
                 exit(-1);
              } 
          } 
      void filter_bank_alloc(); // allocate if filter bank if it is necessary
      int mr_io_fill_header( fitsfile *fptr, int Nl, int Nc, int Nbr_Plan, 
         type_transform Type_Transform,set_transform Set_Transform,
         type_format FormatInputImag, char *Name_MR, type_border Border,
         int MedianWindowSize, float Fc, int Nbr_Iter, 
	 Bool ExactPyrRec, type_lift LiftingTrans, 
	 type_sb_filter SBFilter, sb_type_norm TypeNorm, int NbrUndec, type_undec_filter U_Filter);
	 
      MEYER_WT *MeyerWT;               // Meyer Decomposition class
      FCUR *FCT;                       // Pointer to the Fast Curvelet Transform
      Bool FCT_ExtendWT;               // False by default
      Bool FCT_IsotropWT;              // False by default
      Bool FCT_RealCur;                // True by default     
      intarray TabNbrBandPerResol;     // Number of bands per scale
      intarray FCT_Band2Scale;         // Return the scale number of a given band
      intarray FCT_Band2Dir;           // Return the direction number of a given band
      intarray FCT_Scale2Band;         // Return the first band position of a given scale
    public:
         int FCT_Nbr_Angle;
         FCUR * fast_curvelet_trans() {return FCT;}
	 
         Bool ModifiedATWT;             // It true, use the modified AT WT instead of the
	                                // standard WT
         Bool Verbose;                  // Verbose mode
         type_transform Type_Transform; // type of transformation
         set_transform Set_Transform;   // class of transformation
         type_format FormatInputImag;   // data format of the input image
         char Name_MR[256];             // identifier of the object
         type_border Border;            // Border manadgement type
         int MedianWindowSize;         // window size for MMT and PMT transform 
         float Fc;                      // cut-off frequency for FFT Based
                                        // wavelet transform
         int Nbr_Iter;                  // iteration number for pyramidal 
                                        // wavelet transform. Parameter used
         Bool ExactPyrRec;              // true if the Pyramidal WT must
                                        // have an exact reconstruction 
                                        // if true then Nbr_Iter has no effect
         float SigmaNoise;              // Noise standard deviation in the data
	                                // only used with TM_TO_PYR transform
	 float TMTO_NSigBSpline;        // First parameter for TM_TO_PYR transform
	 float TMTO_NSigMedian;         // Second parameter for TM_TO_PYR transform
	 type_lift LiftingTrans;        // type of lifting (in case of lifting
	                                // transform
         type_sb_filter SBFilter;       // Filter used in the subband 
	                                // decomposition (case of Mallat transform).             
         type_undec_filter U_Filter;    // Filter used in a non orthogonal 3 directional
	                                // and undecimated WT
	 sb_type_norm TypeNorm;         // type of normalization
	                                // (in L1 or L2)

         float NormHaar;                // normalization used for Haar transform
	 Bool EdgeLineTransform;        // Each line and each column are
	                                // transform in a second step in order
					// to better reproduces edges in a
					// image. Only used with TO_MALLAT and
					// TO__DIADIC_HAAR.
         LineCol *LC;                     // Line Col decomposition
        int LC_NbrScaleLine;           // Number of scale in y for Line col decomposition
	 int LC_NbrScaleCol;            // Number of scale in x for Line col decomposition 
	 
	 MultiResol (int Nbr_Line, int Nbr_Col, int Nbr_Scale, 
                     type_transform Transform, char *Name);
         MultiResol ();

         // return a pointer to the filter bank
         FilterAnaSynt * filter_bank() {return FilterBank;}
	 FilterAnaSynt * filter_bank_line() {return FilterBank_Line;}
	 FilterAnaSynt * filter_bank_column () {return FilterBank_Column;}
         UndecSubBandFilter* undec_filter_bank() {return  UndecFilterBank;}
	 MEYER_WT * get_meyer_wt () {return MeyerWT;} 
	 
         // return the number of undecimated scales
         int nbr_undec_scale() {return NbrUndecimatedScale;}

         // return the number of scales
         inline int nbr_scale () const {return Nbr_Plan;}
         
         // return the number of bands
         int nbr_band () const { return NbrBand;}
         
         // return the number of band per resolution
         int nbr_band_per_resol () const { return NbrBand_per_Resol;}
         int nbr_band_per_resol (int b) const { return TabNbrBandPerResol(b);} // FOR FCT, the number of bands is not the same at each scale
	 
         // return the number of  multiresolution coefficients
         // without taking into account the last smooth array
         int nbr_mr_coeff () {return Nbr_MrCoef;}
         
         // return in scale and which_detail, the scale and detail type
         // associated to a band
         void band_to_scale(const int NumBand, int & Scale, details &which_detail)
	  {  band2scale(NumBand, Type_Transform, NbrBand, Scale,which_detail);}
 
         // return the scale associated to a band
         inline int band_to_scale(const int NumBand)
	 {
	    int Scale;
	    details which_detail;
	    band_to_scale(NumBand, Scale, which_detail);
	    return Scale;
	 }
	 
         // return the band associated to a scale and a type of details
         int scale_to_band(int s, details which_detail)  
	        { return scale2band(s, Type_Transform, NbrBand, which_detail);}
		
         // alias to scale_to_band
         inline int stb(int s, details which_detail)  
                         {return scale_to_band(s,which_detail);}
                         
         // return the number of multiresolution coefficients
         int nbr_coeff_in_band(int b) {return TabBandNl[b]*TabBandNc[b];}
         
         // return the number of lines of a band
         int size_band_nl (int NumBand)  const {return TabBandNl[NumBand];}
         
         // return the number of columns of a Band
         int size_band_nc (int NumBand)  const {return TabBandNc[NumBand];}
         
         // return the size of a scale
         int size_scale_nl (int s, details which_detail=D_NULL); 
         int size_scale_nc (int s, details which_detail=D_NULL);
          
         // return the size of the input image
         int size_ima_nl () const {return Nl;}
         int size_ima_nc () const {return Nc;}
         
         // return the scale, the coordinates, and the type of detail
         // of a coefficient
         void pos_coeff(int NumCoef, int &s, 
                        int &i,int &j, details & which_detail);
                        
         // return a scale              
         Ifloat extract_scale (int s, details which_detail=D_NULL);
         // extract a band              
         Ifloat extract_band (int b);
                 
          // insert a scale
         void insert_scale (Ifloat &, int s, details which_detail=D_NULL);
         // insert a band
         void insert_band (Ifloat &, int b);
                  
         // normalize each band to NormVal (==> max(Band) = NormVal)
         // if abs == True, then           (==> max( | Band | ) = NormVal)
         void norm(float NormVal=1, Bool Abs=False);
         
         // allocates data structure of the multiresolution object
         void alloc (int, int, int,  type_transform, char *);

         // allocates the data structure of the multiresolution object
         // A set of filter can be given as a parameter in case
         // of orthogonal transforms.
         void alloc (int Nbr_Line, int Nbr_Col, int Nbr_Scale,
                     type_transform Transform, FilterAnaSynt *FAS=NULL, 
                     sb_type_norm Norm=NORM_L1, int NbrUndec=-1, type_undec_filter U_Filter=DEF_UNDER_FILTER);
	 // computes the number of bands with few and possibly no memory allocation
	 int computeNbBand(int Nbr_Line, int Nbr_Col,int Nbr_Scale,type_transform Transform,int NbrUndec=-1);
         // deallocation of the object
         void free ();
         
         // reference to a coef at scale s, pos i,j and details which_detail
         float & operator() (int s, int i, int j, details which_detail);
         
         // return a coef at scale s, pos i,j and details which_detail with 
         // border manadgement
         float  operator() (int s, int i, int j, details which_detail, type_border Bord);

         // return the  reference to the coeff at position i,j of the band B
         float & operator() (int NumBand, int i, int j) const;
         
          // return the coeff at position i,j of the Band B
         float  operator()  (int NumBand, int i, int j, type_border Bord) const;
         
         // return the coeff at position NumCoef in Band B
         float  operator() (int NumBand, int NumCoef, type_border Bord);
         
         //  return a reference to a  coef at position NumCoef and at band NumBand
         float & operator() (int NumBand, int NumCoef);
 
         // return a reference to a multiresolution coef at position NumCoef
         float & operator() (int NumCoef);
         
         // return ref to band b
         Ifloat & band (int b) const;
	 
	 // return the band array
	 Ifloat *& tabband() {return TabBand;}
	 
         // return ref to scale s
         // obsolete. We keep it only for code compatibility 
         // band should be used instead of scale
         Ifloat & scale (int s) const { return band(s);}
                  
         // wavelet transform of an image
         // if Details = False, and Transform is an isotropic transform
         //   then only the smooth resolutions are calculated (and not
         //   the differences Between two resolutions).
         void transform (Ifloat &Image, type_border Border, Bool Details = True);

         void compute_mod_phase ( Ifloat *&prpo_Modulus,
                                  Ifloat *&prpo_Phase);


         // wavelet transform of an image with default option
         void transform (Ifloat &Image);

         // image reconstruction from the WT
         void recons (Ifloat &Image, type_border Border = DEFAULT_BORDER);

         // image reconstruction from the adjoint operator
         void rec_adjoint (Ifloat &Image, Bool UseLastScale = True, type_border Border = DEFAULT_BORDER);

         // return the normalisation value associated to a band or a scale
         double scale_norm(int s, details which_detail=D_HORIZONTAL);
	 double band_norm(int s);
 
         // read and write a WT file
         void read(char *Name);
         void write(char *Name);
	 
         // read and write a scale from a file
         void read(char *Name, int SelectScale);
         void write(char *Name, int SelectScale);
 
          // read and write a band
	  void read_band (char *Name, int SelectBand);
 	  void write_band(char *Name, int SelectBand);

	  // return the detection level for a given band b, 
	  // in case of Gaussian noise (Sigma = standard deviation)
	  // with NSgima detection level
	  double gauss_detect_level(int b, float Sigma=1., float NSgima=3.);
	  
	  // threshold a band
	  // if (UseAbsVal == False) then all negative values are thresholded
	  void threshold(int b, float Level, Bool UseAbsVal=True);
	  
	  // threshold a band
	  // if (UseAbsVal == False) then all negative values are thresholded
 	  void threshold(float Sigma=1., float NSigma=3, Bool UseAbsVal=True);
	  
          void print_info();
          // print some information about the object
       
          ~MultiResol();
};


/* Multiresolution transform procedures */
void mr_correct_pyr (Ifloat &Imag, MultiResol &MR_Data, 
                     int Max_Iter=DEFAULT_ITER_COR_PYR);
void wave_cf_transform (Ifloat &Imag_in, MultiResol &MR_Data);
void mr_cf_transform (Ifloat &Image, MultiResol &MR_Transf);

void wave_2d_mallat_atrou (Ifloat &pro_ImagIn, Ifloat *&prpo_ImagOut,   
                         int pi_NbrPlan,  type_border ps_Border=DEFAULT_BORDER);
void wave_max_2d_mallat_atrou (Ifloat &pro_ImagIn, Ifloat *&prpo_ImagOut,   
                         int pi_NbrPlan,  type_border ps_Border=DEFAULT_BORDER);
				
void transform_g (Ifloat &Imag_in, Ifloat &Imag_Out, int Iter);
void transform_g (Iint &Imag_in, Iint &Imag_Out, int Iter);

void transform_feauveau (Ifloat &Imag_in, Ifloat &Imag_Out, int Iter);
void transform_feauveau_trou (Ifloat &Imag_in, MultiResol &MR_Data);

void haar_2d_transform (Ifloat &Imag_in, MultiResol &MR_Data);
void haar_2d_atrou_dir2_transform (Ifloat &Data, MultiResol &T_Haar, Bool EdgeLineTransform=False);
void haar_2d_atrou_dir3_transform (Ifloat &Data, MultiResol &T_Haar);

double mr_scale_norm(int s, type_transform Transform, 
                  int MedianWindowSize=DEFAULT_MEDIAN_WINDOW_SIZE, details which_detail=D_HORIZONTAL,
		  type_undec_filter U_Filter=U_B2SPLINE);
		 
/* Multiresolution reconstruction procedures */

void wave_cf_recons (MultiResol &MR_Data, Ifloat &Imag_Out);

// void mr_recons (MultiResol &MR_Transf, Ifloat &Image, 
//                type_border Border=DEFAULT_BORDER);
// void mr_transform (Ifloat &Image, MultiResol &MR_Transf, type_border Border=DEFAULT_BORDER, Bool Details=True);
// void mallat_2d_transform (Ifloat &Imag_in, MultiResol &MR_Data);
// void mallat_2d_reconstruct (MultiResol &MR_Data, Ifloat &Imag_Out);
//void mr_pyr_tmto (Ifloat &Image, MultiResol &MR_Transf, float & SigmaNoise, 
//                 float NSigBSpline=DEFAULT_N_SIGMA, float NSigMedian=5);
// void mr_adjoint_rec (MultiResol &MR_Transf, Ifloat &Image, 
//                       type_border Border=DEFAULT_BORDER);

void rec_wave_2d_mallat_atrou (Ifloat *&prpo_ImagIn,  Ifloat &pro_ImagOut, 
                         int pi_NbrPlan,  type_border ps_Border=DEFAULT_BORDER);
void rec_wave_max_2d_mallat_atrou (Ifloat *&prpo_ImagIn,  Ifloat &pro_ImagOut, 
                         int pi_NbrPlan,  type_border ps_Border=DEFAULT_BORDER);
void haar_2d_reconstruct (MultiResol & T_Haar, Ifloat &Data);
void haar_2d_atrou_dir2_recons (MultiResol & T_Haar, Ifloat &Data, Bool EdgeLineTransform=False);
void haar_2d_atrou_dir3_recons (MultiResol & T_Haar, Ifloat &Data);
// void haar_2d_adjoint (MultiResol &T_Haar, Ifloat &Data, type_border Border=DEFAULT_BORDER);
void inverse_transform_g (Ifloat &Imag_in, Ifloat &Imag_Out, int Iter);
void inverse_transform_g (Iint &Imag_in, Iint &Imag_Out, int Iter);

void recons_feauveau (Ifloat &Imag_in, Ifloat &Imag_Out, int Iter);

void ind_orthog_transf(int s, int i,int j, details which_detail, 
                       int Nl, int Nc, int &ind_i, int &ind_j);
void ortho_trans_to_ima(MultiResol &MR_Data, Ifloat &Ima, int Resol=0);
void ima_to_ortho_trans(MultiResol &MR_Data, Ifloat &Ima, int Resol=0);
                       
float mr_tab_noise (int s);

void mr_ortho_regul_ima_rec(MultiResol &MR_Step, Ifloat &Result, 
                            float *TabLevel, int MaxIter, int SizeTabLevel=3);
// MR_Step = in: orthogonal transform with 2 scales (4 bands)
// Result = out: reconstructed image under laplacian constaint
// MaxIter = in: number of iterations
// TabLevel = in: detection level per scale
// SizeTabLevel = in: size of TabLevel
//                if SizeTabLevel = 1 the level is the same for the 3 bands
//                                        i.e. TabLevel[0]
//                if SizeTabLevel = 3, one level per band
//                if SizeTabLevel > 3, one level per wavelet coefficient

void mr_ortho_regul_rec (MultiResol &MR_Data, Ifloat &Result, 
              float *LevelQ, int MaxIter=10, Bool ModifCoef=False);

// median multiresolution tranform: routine working with integer
#define DEF_IPYRMED_WINDOW_SIZE 3
void mri_pos_ind_medpyr (int *Tab_Nl, int *Tab_Col, int *Tab_Pos, 
                         int Nl, int Nc, int Nbr_Plan);
int mri_size_medpyr (int Nl, int Nc, int Nbr_Plan);
void imI_bilinear_interp (int *INimage, int INlin, int INcol, 
                          int *OUTimage, int OUTlin, int OUTcol);
void mri_pyrmed (int *Ima, int *Pyramid, int *Tab_Nl, 
                 int *Tab_Col, int *Tab_Pos, int Nbr_Etap, 
                 int MedianWindowSize=DEF_IPYRMED_WINDOW_SIZE);
void mri_pyrmedian_transform (int *Data, int Nl, int Nc, int **Median, int Nbr_Plan, int MedianWindowSize);
void mri_pyrmedian_rec (int *Pyramid, int *PictInt, 
                        int Nl, int Nc, int Nbr_Plan, float *Tab_Level);
void morphoi_noise_rec (int *Pyramid, int *PictInt, 
                        int Nl, int Nc, int Nbr_Plan);




void search_number_max (Ifloat*& prpo_ModData, 
                        intarray& pro_NumberMaxInLine,
                        intarray& pro_NumberMaxInRow,
                        int pi_NBrBand, Bool Verbose=False);

void search_max (Ifloat*& prpo_ModData, 
                 intarray& pro_NumberMaxInLine,
                 intarray& pro_NumberMaxInRow,
                 int pi_NBrBand);

void init_max (int pi_NbrPlan,                // number of plan
               int pi_Nl,                     // number of line
               int pi_Nr,                     // number of row
	           intarray& pro_NumberMaxInLine, // Number of max in line
               intarray** pppo_MaxModInLine,  // loc of max mod on line
               intarray& pro_NumberMaxInRow,  // Number of max in line
               intarray** pppo_MaxModInRow,   // loc of max mod on lin
               MultiResol& pro_Mr2dData,      // Data IN
               MultiResol& pro_Mr2dRecData);  // max modulus of Data IN

void init_last_scale (int pi_NbrPlan,                // number of plan
                      int pi_Nl,                     // number of line
                      int pi_Nr,                     // number of row
                      MultiResol& po_Mr2dData,       // Data IN
                      MultiResol& po_Mr2dRecData);   // Rec Data OUT

void interpolate  (intarray& pro_NumberMaxInLine,   // Number of max in line
                   intarray** pppo_MaxModInLine,    // loc of max mod on line
                   intarray& pro_NumberMaxInRow,    // Number of max in line
                   intarray** pppo_MaxModInRow,     // loc of max mod on line
                   MultiResol& pro_Mr2dData,        // Data IN
                   MultiResol& pro_Mr2dRecData);    // Rec Data OUT

// interpolate between two consecutive maxima modulus at row pi_Row
void Ortho_Proj_Operator (int pi_Scale,            // current scale 
			              int pi_Line,             // current line 
                          int pi_Row ,             // current rowe
                          int pi_Begining,         // ind of begining
	                      int pi_End,              // ind of end
		                  MultiResol& pro_Mr2dData,       // Data IN
			              MultiResol& pro_Mr2dRecData);  // Rec Data OUT

void compute_max_modulus (Ifloat*& prpo_ModData, 
                          Ifloat*& prpo_PhaData, 
                          Bool     pe_Verbose,
			  int      pi_NbrPlan);
			  
			  

void init_coord_max (Ifloat*&    prpo_ModData, 
                     int         pi_Nl, 
		     int         pi_Nr,
		     intarray&   pro_NumberMaxInLine, 
		     intarray**  pppo_MaxModInLine,
                     intarray&   pro_NumberMaxInRow, 
		     intarray**  pppo_MaxModInRow,
		     int         pi_NbrPlan) ;
		     
void init_multiresol_out (int         pi_Nl, 
                          int         pi_Nr, 
			  int         pi_NbrPlan, 
			  intarray&   pro_NumberMaxInLine, 
			  intarray**  pppo_MaxModInLine,
                          intarray&   pro_NumberMaxInRow,  
			  intarray**  pppo_MaxModInRow, 
			  MultiResol& pro_Mr2dData, 
			  MultiResol& pro_Mr2dRecData);
			  
			  		     		     
                       
#endif
