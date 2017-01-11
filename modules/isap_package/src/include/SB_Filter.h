/******************************************************************************
**                   Copyright (C) 1998 by CEA
*******************************************************************************
**
**    UNIT
**
**    Version: 1.0
**
**    Author: Jean-Luc Starck
**
**    Date:  18/09/98 
**    
**    File:  SB_Filter.h
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

#ifndef _SB_FILTER_H_
#define _SB_FILTER_H_

#include "SB_Filter1D.h"
#include "LineCol.h"

/***********************************************************************
* 2D  Decomposition (one step transform)
* 2D sub-band decomposision: an image is transformed in four sub-images
* For a decimated transform, an image is transformed into four images
* of half-size.
* For a undecimated transform, an image is transformed into four images
* with the same size as the input image.
************************************************************************/

class SubBand2D {
        SubBand1D *Ptr_SB1D_LINE; // 1D subband decomposition along the lines
	SubBand1D *Ptr_SB1D_COL;  // 1D subband decomposition along the columns
        void set_sub() {Line_SubSample_H_Even=Line_SubSample_G_Odd=True;
	                Col_SubSample_H_Even=Col_SubSample_G_Odd=True;}
      public:
         Bool Line_SubSample_H_Even;
	 Bool Line_SubSample_G_Odd;
	 Bool Col_SubSample_H_Even;
	 Bool Col_SubSample_G_Odd;
         SubBand1D *get_subband_method() {return Ptr_SB1D_LINE;}
 	 SubBand1D *get_subband_method_line() 
	     { Ptr_SB1D_LINE->SubSample_H_Even = Line_SubSample_H_Even;
	       Ptr_SB1D_LINE->SubSample_G_Odd = Line_SubSample_G_Odd;
 	       return Ptr_SB1D_LINE;}
	 SubBand1D *get_subband_method_col() 
	     { Ptr_SB1D_COL->SubSample_H_Even = Col_SubSample_H_Even;
	       Ptr_SB1D_COL->SubSample_G_Odd = Col_SubSample_G_Odd;
 	       return Ptr_SB1D_COL;}
	 
         SubBand2D(SubBand1D &SB1D) {Ptr_SB1D_LINE = &SB1D; Ptr_SB1D_COL = &SB1D;set_sub();}
	 // allocate the class with the same 1D subband decomposition
	 // for both the lines and the columns
	 
	 SubBand2D(SubBand1D &SB1D_Line, SubBand1D &SB1D_Col) 
	              {Ptr_SB1D_LINE = &SB1D_Line; Ptr_SB1D_COL = &SB1D_Col;set_sub();}
	 // allocate the class with tho 1D subband decompositions
	 // for the lines and the columns
	 
	 void transform2d(Ifloat &Data, Bool UseSameBuff=True,
	                 Ifloat *Horiz=NULL, Ifloat *Vert=NULL, 
			 Ifloat *Diag=NULL, Ifloat *Smooth=NULL);
	 // Decimated sub-band decomposition
	 // if UseSameBuff==True then the result is stored in Data
	 // if UseSameBuff==False then the Image *Horiz,*Vert,*Diag,*Smooth
	 // must be allocated
	 
	 void recons2d(Ifloat &Data, Bool UseSameBuff=True,
	                 Ifloat *Horiz=NULL, Ifloat *Vert=NULL, 
			 Ifloat *Diag=NULL, Ifloat *Smooth=NULL);
	 // Decimated sub-band reconstruction
	 // if UseSameBuff=True then Data contains the input and the output
	 // if UseSameBuff=False, the input is *Horiz,*Vert,*Diag,*Smooth
	 // and the output is Data
	 
	 void transform2d(Ifloat &Data, Ifloat &Horiz, Ifloat &Vert, 
	                  Ifloat &Diag, Ifloat &Smooth, int Step);
	 // Undecimated sub-band decomposition
	 // Step is distance between pixels which much be taken into account
	 // when the 1D filter bank is applied (a trous algorithm).
	 // Data = input image (must be allocated)
	 // Horiz,Vert,	Diag,Smooth = output images of the same size as Data
	 //                           they must be allocated before the call
	 	  
	 void recons2d(Ifloat &Horiz, Ifloat &Vert, 
	               Ifloat &Diag, Ifloat &Smooth, Ifloat &Output, int Step);
	 // Undecimated sub-band reconstruction
	 // Data = output image
	 // Horiz,Vert,	Diag,Smooth = input images 
	 ~SubBand2D(){Ptr_SB1D_LINE = NULL; Ptr_SB1D_COL = NULL;}
};

/***********************************************************************
*   Orthogonal 2D wavelet transform
************************************************************************/
 
class Ortho_2D_WT: public SubBand2D {
      public:
          Ortho_2D_WT(SubBand1D &SB1D):SubBand2D(SB1D) {};
	  Ortho_2D_WT(SubBand1D &SB1D_LINE, SubBand1D &SB1D_COL):SubBand2D(SB1D_LINE,SB1D_COL) {};
	  void transform (Ifloat &Imag, int Nbr_Plan);
	  void transform (Iint &Imag, int Nbr_Plan);
	  void transform(Ifloat &Imag_in, Ifloat  &Transf_Out, int Nbr_Plan);
 	  void transform(Iint &Imag_in, Iint  &Transf_Out, int Nbr_Plan);
	  void recons(Ifloat &Imag, int Nbr_Plan);
	  void recons(Iint &Imag, int Nbr_Plan);
          void recons(Ifloat &Transf_in, Ifloat &Imag_Out, int Nbr_Plan);
          void recons(Iint &Transf_in, Iint &Imag_Out, int Nbr_Plan);	  
};

/***********************************************************************
*   Orthogonal Undecimated 2D wavelet transform
************************************************************************/

class PAVE_2D_WT: public SubBand2D {
      public:
          PAVE_2D_WT(SubBand1D &SB1D):SubBand2D(SB1D) {};
	  PAVE_2D_WT(SubBand1D &SB1D_LINE, SubBand1D &SB1D_COL):SubBand2D(SB1D_LINE,SB1D_COL) {};
	  void one_scale_transform (Ifloat &Imag, Ifloat *TabTrans, int Step, int Pos=0);
	  void one_scale_recons(Ifloat *TabTrans, Ifloat &Imag, int Step, int Pos=0);
	  void transform (Ifloat &Imag, Ifloat *TabTrans, int Nbr_Plan);
	  void recons(Ifloat *TabTrans, Ifloat &Imag, int Nbr_Plan);
          int alloc (Ifloat * & TabBand, int Nl, int Nc, int Nbr_Plan);
          void free(Ifloat *TabBand, int Nbr_Plan);
};

/***********************************************************************
*  Orthogonal decimated-Undecimated 2D wavelet transform
************************************************************************/

class HALF_DECIMATED_2D_WT {
        SubBand1D *Ptr_SB1D_LINE; // 1D subband decomposition along the lines
	SubBand1D *Ptr_SB1D_COL;  // 1D subband decomposition along the columns
        void  set_tabdec(int NumUndec, Bool * & TabDec, int Nbr_Plan);
       public:
          Bool FiszTrans;
	  Bool FiszRec;
          HALF_DECIMATED_2D_WT (SubBand1D &SB1D) 
	    {Ptr_SB1D_LINE = &SB1D; Ptr_SB1D_COL = &SB1D;FiszTrans=FiszRec=False;}
	  HALF_DECIMATED_2D_WT (SubBand1D &SB1D_LINE, SubBand1D &SB1D_COL) 
	    {Ptr_SB1D_LINE = &SB1D_LINE; Ptr_SB1D_COL = &SB1D_COL;
	     FiszTrans=FiszRec=False;}
 	  int alloc(Ifloat * & TabTrans, int Nl, int Nc, int Nbr_Plan, Bool *TabDec);
	  int alloc(Ifloat * & TabTrans, int Nl, int Nc, int Nbr_Plan, int NbrUndecimatedScale);
	  void free(Ifloat *TabTrans, int Nbr_Plan);
 	  void transform (Ifloat &Imag, Ifloat *TabTrans, 
	      int Nbr_Plan, Bool *TabDec);
          void transform(Ifloat &Imag, Ifloat *TabTrans, 
	      int Nbr_Plan, int NumUndec);
	  void recons(Ifloat *TabTrans, Ifloat &Imag, int Nbr_Plan, 
	              Bool *TabDec);
          void recons(Ifloat *TabTrans, Ifloat &Imag, int Nbr_Plan, int NumUndec);
         ~HALF_DECIMATED_2D_WT(){Ptr_SB1D_LINE = NULL; Ptr_SB1D_COL = NULL;}
};

void fisz2d(Ifloat &Data, Ifloat &Recons, Bool Reverse=False, 
                  int NbrScale=-1, type_sb_filter Fil = F_HAAR,  
		  int NumUndec = -1);

inline void fisz2d_trans(Ifloat &Data) {fisz2d(Data,Data);}
inline void fisz2d_inv(Ifloat &Data) {fisz2d(Data,Data,True);}

/************************************************************************/
//    Isotropic a trous wavelet transform
/************************************************************************/

class ATROUS_2D_WT {
    UndecSubBandFilter *USBF;
public:
    Bool ModifiedAWT;  
    Bool AdjointRec;
    type_border Bord;
    ATROUS_2D_WT () {AdjointRec=False;Bord=I_CONT;ModifiedAWT=False;
        USBF = new UndecSubBandFilter(U_B3SPLINE);}
    void alloc(Ifloat * & TabTrans, int Nl, int Nc, int Nbr_Plan);
    void free(Ifloat *TabTrans, int Nbr_Plan);
    void b3spline_filtering(Ifloat & Im_in, Ifloat &Im_out, int Step_trou, Bool GTilde=False, int nb_thr=0);
    inline void b3spline_filtering(Ifloat & Im_in, Ifloat &Im_out, int Step_trou, int nb_thr=0) {b3spline_filtering(Im_in,Im_out,Step_trou,False,nb_thr);};//FCS ADDED
    void transform(Ifloat &Image, Ifloat *TabBand, int Nbr_Plan, int nb_thr=0);
    float norm_band(int s);
    void recons(Ifloat *TabTrans, Ifloat &Imag, int Nbr_Plan, Bool AddLastScale, int nb_thr=0);
    inline void recons(Ifloat *TabTrans, Ifloat &Imag, int Nbr_Plan,int nb_thr) {recons(TabTrans,Imag,Nbr_Plan,True,nb_thr);};
    inline void recons(Ifloat *TabTrans, Ifloat &Imag, int Nbr_Plan) {recons(TabTrans,Imag,Nbr_Plan,True,0);};
    ~ATROUS_2D_WT(){delete USBF;}
};

/************************************************************************/

class PMT_2D  {
        public:
          int MedianWindowSize;
          PMT_2D () {MedianWindowSize=5;}
	  void alloc(Ifloat * & TabTrans, int Nl, int Nc, int Nbr_Plan);
	  void free(Ifloat *TabTrans, int Nbr_Plan);
          void transform(Ifloat &Image, Ifloat *TabBand, int Nbr_Plan, type_border Bord=I_CONT);
          float norm_band(int s);
	  void recons(Ifloat *TabTrans, Ifloat &Imag, int Nbr_Plan, type_border Bord=I_CONT);
         ~PMT_2D(){}
};


/****************************************************************************/
// 2D decimated  Quincunx decomposition
/*************************************************************/

class  Quincunx {
       void ind_quinc(int Band, int i,int j,  
                      int Nl, int Nc, int &ind_i, int &ind_j);
       // return the position in a image of size Nl,Nc of pixel i,j in
       // the band of Band (Band = 0 .. NbrBand-1)
       SubBand1D *Ptr_SB1D;
      public:
        Quincunx (SubBand1D &SB1D) {Ptr_SB1D = &SB1D;}
        // Constructor: SB1D = subband method to use

 	int alloc(Ifloat * & TabTrans, int Nl, int Nc, int Nbr_Plan);
        // return the number of bands and allocate TabTrans
        //    TabTrans[s] = band number s (s = 0..NbrBand-1)
        //    Nl,Nc = image size
        //    Nbr_Plan = Number of scales

 	void free(Ifloat *TabTrans, int Nbr_Plan);
        // deallocates TabTrans 

        SubBand1D *get_subband_method() {return Ptr_SB1D;}
        // return the subband decomposition method

	void transform_one_step (Ifloat & Data, Ifloat &Half, Ifloat & Resol, 
                                 Ifloat & Smooth);
        // transform on image in three subimages by quincunx decomposition

	void recons_one_step(Ifloat &Smooth, Ifloat & Resol, 
	                         Ifloat & Half, Ifloat &ImaRec);
        // reconstruct on image from its three subimages 
        // by quincunx reconstruction

        void transform(Ifloat &Image, Ifloat * & TabBand, int Nbr_Plan);
        // quincunx transform of an image

	void recons(Ifloat *TabTrans, Ifloat &Imag, int Nbr_Plan);
        // quincunx reconstruction of an image
 
	void band_to_ima(Ifloat *TabTrans, Ifloat &Imag, int Nbr_Plan);
        // put all coefficient in Imag
        // Imag must be allocated before calling band_to_ima

	~Quincunx (){Ptr_SB1D = NULL;}
};

/************************************************************************/

#endif
