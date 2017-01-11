/******************************************************************************
**                   Copyright (C) 1994 by CEA
*******************************************************************************
**
**    UNIT
**
**    Version: 3.3
**
**    Author: Jean-Luc Starck
**
**    Date:  96/06/13 
**    
**    File:  MR1D_Obj.h
**
**    Modification history :
**            23-OCT-1996 R Gastaud add pos_mr1dcoeff
**            24-OCT-1996 R Gastaud add nbr_mr_coeff, SetTransform
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

#ifndef _SI_OBJ_H_
#define _SI_OBJ_H_

#include "IM_Obj.h"  
#include "SB_Filter1D.h"

//enum type_border{I_CONT, I_MIRROR, I_ZERO};
//#define DEFAULT_BORDER I_CONT

//#define FLT_MAX 1e99
#define MAX_SCALE_1D 100
#define MAX_SIZE_LAST_SCALE 4

#define DEFAULT_MORLET_NU_0 0.8 
#define DEFAULT_MEDIAN_1D_WINSIZE 5
#define DEFAULT_BORDER_1D I_MIRROR

#define NBR_TRANS_1D 21
#define NBR_DIADIC_TRANS_1D 11

enum type_trans_1d { TO1_PAVE_LINEAR, 
		     TO1_PAVE_B1SPLINE, 
		     TO1_PAVE_B3SPLINE,
		     TO1_PAVE_B3_DERIV,
		     TO1_PAVE_HAAR,
		     TM1_PAVE_MEDIAN,
		     TU1_MALLAT,
		     TU1_UNDECIMATED_NON_ORTHO,
		     TO1_PAVE_B3SPLINE_GEN2 ,
		     TO1_PYR_B3SPLINE, 
		     TM1_PYR_MEDIAN,
		     TO1_PAVE_MORLET,
		     TO1_PAVE_MEX,
		     TO1_PAVE_FRENCH,
		     TO1_PAVE_DERIV_GAUSS,
		     TO1_MALLAT,
		     TO1_LIFTING,
		     WP1_MALLAT,
		     WP1_LIFTING,
 		     WP1_ATROUS,
		      TO1_PYR_LINEAR,
                     T1_UNDEFINED=-1};

enum set_trans_1d {TRANS1_PAVE, TRANS1_PYR, TRANS1_MALLAT, 
                   TRANS1_WP_MALLAT, TRANS1_WP_ATROUS,
		   S1_UNDEFINED=-1};

set_trans_1d which_set_is_trans1d(type_trans_1d type);
const char * StringTransf1D (type_trans_1d type);

/****************************************************************************/

inline void reform_to_1d(fltarray &Dat)
{
   int Naxis = Dat.naxis();
   int Nelem = Dat.n_elem();
   Bool Reform = False;

   if (Naxis > 1)
   {
      for (int i=1; i <= Naxis; i++)
        if ( (Dat.axis(i) > 1) && (Dat.axis(i) != Nelem)) Reform = True;
      if (Reform == True)
      {
         cerr << "Error: data are not mono-dimensional ... " << endl;
         exit(-1);
      }
      Dat.reform(Dat.n_elem());
   }
}

/****************************************************************************/
 
inline int border_ind_test (int ind, int N, type_border Border)
{
    int Val=ind;
    
    switch (Border)
    {
       case I_MIRROR:
         if (ind < 0) Val = -ind;
         else
         {
             if (ind >= N) Val = 2 * N - 2 - ind;
             else Val = ind;
         }
         break;
       case I_CONT:
         if (ind < 0) Val = 0;
         else if (ind >= N) Val = N-1;
         else Val = ind;
         break;
       default:
         cerr << "This border is not implemented in 1D routines ... " << endl;
         break;
    }
    return (Val);
}

/****************************************************************************/

class MR_1D {
         fltarray Data; // buffer where the wavelet transform is stored
         Bool FilterBankAlloc;      // True if the filter bank is allocated
                                    // by the class
         FilterAnaSynt *FilterBank; // pointer to the filter bank to use
                                    // in the orthogonal wavelet transform
         UndecSubBandFilter *UndecFilterBank; // pointer to an undecimated filter bank.
	                           // used in a non orthogonal undecimated tri-directional transform.
 	 int Nbr_Plan;  // number of scales
	 int Nbr_Band;  // number of band (used by wavelet packets
	                // transform
         int Np;        // number of points of the input signal
	                
         int *TabPos;   // TabPos and TabSize gives the position
         int *TabSize; // and the size of a given band.
	               // used only for non-redundant transform
		       // (i.e. size(signal) = siwe(transform))
	 void wp_pos(int N, int NStep, int & IndTab);
	               // initialize TabPos and TabSize
         void reset ();
      public:
         char Name_MR[256]; // object name
         type_trans_1d Type_Transform; // type of transform
         set_trans_1d  Set_Transform;  // class of transform
         type_border Border;           // type of border to user for the
	                               // border intrpolation
         int MedianWinSize;  // window size used by the median transforms
         int Nbr_Voie;       // number of voices per octaves
         float Scale_0;      // size of the first scale
	                     // when using continuous transform
         float Nu_0;         // Morlet parameter (only used for the 
	                     // morlet transform
	 float BorderSize_0; // border size at the first scale
         Bool Interp;        // if true interpolate the transform
	                     // (only used for pyramidal transform)
	 Bool KillBorder;    // set to zero all coeffients contaminated
	                     // by the border
	 type_sb_filter SB_Filter; // Filter for the (bi-orthogonal transform)
	                           // and wavelet paquet
	 type_undec_filter U_Filter;    // Filter used in a non orthogonal  
	                                // and undecimated WT
         type_lift LiftingTrans;   // type of lifting scheme 
	 sb_type_norm Norm;        // Norm (L1 or L2)
	 
         MR_1D (int N, type_trans_1d T, char *Name, int Nbr_Scale = 0, 
                Bool Interp=False,
                float S0=0., float Nu0=DEFAULT_MORLET_NU_0, int N_V=0);
	 // object constructor 	
	 
         MR_1D (){ reset();}
         // object constructor
	 
         inline int nbr_scale () const {return Nbr_Plan;}
	 // return the number of scales
	 
	 inline int nbr_band () const {return Nbr_Band;}
	 // return the number of bands
	 
         //fltarray scale (int s);
	 // return  a scale of the wavelet transform in a 1D array
	 // PB: a buffer is created containing the scale which is not
	 // deallocated. It is recommended to use the next routine
	 
	 void scale (fltarray & TabScale, int s);
	 // return in TabScale a scale of the wavelet transform
	 
         int size_scale_np (int) const;
	 // return the size of a given scale
	 
	 int pos_band(int s) 
	 { 
	    int Pos=0;
	    if (TabPos == NULL) for (int b=0; b < s-1; b++) Pos += size_scale_np(b);
	    else Pos = TabPos[s];
	    return Pos;
	 }
	 
         int size_ima_np () const {return Np;}
	 // return the size of the input signal
	 
         void alloc (int N, type_trans_1d T, char *Name, int Nbr_Scale, 
                float S0=0., float Nu0=DEFAULT_MORLET_NU_0, int N_V=0);
         // allocate the object

         void  alloc (int N, type_trans_1d T, int Nbr_Scale,  
                   FilterAnaSynt *FAS=NULL, sb_type_norm TNorm=NORM_L1,
                   Bool Inter=False, float S0=0., 
                   float Nu0=DEFAULT_MORLET_NU_0, int N_V=0)	;
         // allocate the object

         void free ();
	 // deallocate the object
	 
         float & operator() (int s, int i) const;
 	 // return one wavelet coefficient at scale s and position i
	 
         fltarray & image() {return Data;}
	 // return the computed wavelet transform in a 2D array
	 
         void transform (fltarray &Image, type_border Border);
	 // apply the wavelet transform to a 1D signal
	 // using a given border interpolation
	 
         void transform (fltarray &Image);
	 // apply the wavelet transform to a 1D signal
	 // using the default border interpolation
	 
         void recons (fltarray &Image, type_border Border);
	 // reconstruction of the 1D signal from its wavelet coefficient
	 // using a given border interpolation

         void recons (fltarray &Image);
	 // reconstruction of the 1D signal from its wavelet coefficient
	 // using the default border interpolation
 
         void rec_adjoint (fltarray &Image, Bool UseLastScale = True, type_border Border = DEFAULT_BORDER);
         // reconstuction by the adjoint operator
	 
         int nbr_mr_coeff (); // return the number of wavelet coefficients
	 int border_size(int Scale); // return the border size  
	  
 	 void kill_border(); 
	 // suppress the wavelet coefficients
	 // computed using an interpolation at the border
	 // border problem

        // return a pointer to the filter bank
         FilterAnaSynt * filter_bank() {return FilterBank;}
         UndecSubBandFilter* undec_filter_bank() {return  UndecFilterBank;}
	 
         void modulus_maxima(); 
	 // set to zero all wavelet coefficients which are not
	 // modulus maximum
	 
         void loc_maxima(); 
	 // local maximum

         void loc_optima(); 
	 // local optima
	 
         ~MR_1D();
};

void wave_1d_mex (fltarray & Signal, fltarray & W_1D, 
                     int N, type_border Border,
             int Nbr_Voie, int Nbr_Plan, float Scale_0);
void wave_1d_mex_rec (fltarray &W_1D, fltarray &Signal, int N, 
                     type_border Border,
                     int Nbr_Voie, int Nbr_Plan, float Scale_0);
void wave_1d_french (fltarray & Signal, fltarray & W_1D, int N, 
                      type_border Border,
                int Nbr_Voie, int Nbr_Plan, float Scale_0);
void wave_1d_der_gauss (fltarray & Signal, fltarray & W_1D, int N,type_border Border, 
             int Nbr_Voie, int Nbr_Plan, float Scale_0);
	     
void wave_1d_french_rec (fltarray & W_1D, fltarray & Signal, int N,
                         type_border Border, 
                         int Nbr_Voie, int Nbr_Plan, float Scale_0);
void wave_1d_linear (fltarray &Signal, fltarray & W_1D, int N, int Nbr_Plan,
                          type_border Border = DEFAULT_BORDER);
void wave_1d_haar (fltarray &Signal, fltarray & W_1D, int N, int Nbr_Plan,
                          type_border Border = DEFAULT_BORDER);
void wave_1d_spline1 (fltarray &Signal, fltarray & W_1D, int N, int Nbr_Plan,
                type_border Border = DEFAULT_BORDER);
void wave_1d_spline3 (fltarray &Signal, fltarray & W_1D, int N, int Nbr_Plan,
                type_border Border = DEFAULT_BORDER);

void wave_1d_spline3_gen2 (fltarray &Signal, fltarray & W_1D, int N, int Nbr_Plan,
                           type_border Border= DEFAULT_BORDER);

void wave_1d_spline3_gen2_rec (fltarray &W_1D, fltarray & Signal, int N, int Nbr_Plan);

void wave_1d_algo_trou_rec (fltarray &W_1D, fltarray & Signal, int N, int Nbr_Plan);
void wave_1d_morlet (fltarray & Signal, fltarray & W_1D, 
                int N, type_border Border, 
                int Nbr_Voie, int Nbr_Plan, float Nu_0, float Scale0);
void morlet_mod (fltarray &Wave, fltarray &Wave_Mod);
void morlet_phase (fltarray &Wave, fltarray &Wave_Phase);
void morlet_re (fltarray &Wave_Mod, fltarray &Wave_Ph, fltarray &Wave);
void morlet_im (fltarray &Wave_Mod, fltarray &Wave_Ph, fltarray &Wave);
void mr1d_median (fltarray &Signal, fltarray &W_1D, int Np, 
                  int Nbr_Plan, int MedianWindowSize=DEFAULT_MEDIAN_1D_WINSIZE,
                  type_border Border=DEFAULT_BORDER);
void wave_1d_B3deriv_atrou (fltarray &Signal, fltarray & W_1D,   
                           int Nbr_Plan, type_border Border=DEFAULT_BORDER);
void wave_1d_rec_B3deriv_atrou (fltarray &W_1D, fltarray & Signal, 
                                int Nbr_Plan, type_border Border=DEFAULT_BORDER);
						  
void mr1d_transform (fltarray & Signal, int Np, type_trans_1d Type_Transform, 
                type_border Border, int MedianWinSize,
                int Nbr_Voie,  fltarray & W_1D, int Nbr_Plan, 
                float Scale_0, float Nu_0, Bool Interp=False);
void pyr_1d_linear (fltarray &Signal, fltarray & W_1D, int Np, int Nbr_Plan,
                type_border Border = DEFAULT_BORDER);
void pyr_1d_spline3 (fltarray &Signal, fltarray & W_1D, int Np, int Nbr_Plan,
                type_border Border = DEFAULT_BORDER);
void mr1d_pyr_median (fltarray &Signal, fltarray &W_1D, int Np, int Nbr_Plan,
                      int WinSize=DEFAULT_MEDIAN_1D_WINSIZE, 
                      type_border Border= DEFAULT_BORDER);
void mr1d_pyr_rec (fltarray &Signal, fltarray &W_1D, int Np, int Nbr_Plan);
void mr1d_recons (fltarray &W_1D, int Np, int Nbr_Plan, type_trans_1d T_Transf,
                 type_border Border, 
                 int Nbr_Voie, fltarray & Signal, float Scale_0);

void filt1d_mediane (const fltarray &S1, fltarray &S2, int N, 
                     int Window_Size=DEFAULT_MEDIAN_1D_WINSIZE, 
                     type_border Border = DEFAULT_BORDER);

void pos_mr1dcoeff(int NumCoef, int &s, int &i, 
                set_trans_1d Set_Transform, int Nelem, int Nbr_Plan);

set_trans_1d SetTransform (type_trans_1d Transform);

void wp1d_atrous_transform (fltarray &Data, fltarray &WP,
                            int Nstep, type_border Border=DEFAULT_BORDER,
			    int StartPos=0, int NumStep=0);
			    
void mallat_1d_transform (fltarray &Data, fltarray &Mallat, int Nbr_Plan, SubBand1D &SBF);
void mallat_1d_reconstruct (fltarray &Mallat, fltarray &Data, int Nbr_Plan, SubBand1D &SBF);
// void lift1_1d_transform (fltarray &Data, fltarray &Lift, int Nbr_Plan);
// void lift1_1d_reconstruct (fltarray &Lift, fltarray &Data, int Nbr_Plan);

void wp1d_mallat_transform (fltarray &Data, fltarray &Mallat,
                            int Nstep, SubBand1D &SBF);
void wp1d_mallat_rec (fltarray &Mallat, fltarray &Data, int Nstep, SubBand1D &SBF);
#endif
