
/******************************************************************************
**                   Copyright (C) 1999 by CEA
*******************************************************************************
**
**    UNIT
**
**    Version: 1.0
**
**    Author: J.L. Starck
**
**    Date:  19.8.99
**    
**    File:  MR3D_Obj.h
**
*******************************************************************************
**
**    DESCRIPTION  Class for 3D wavelet transform
**    ----------- 
**                 
******************************************************************************/

#ifndef _CMR3D_H_
#define _CMR3D_H_

#include "IM_Obj.h"
#include "IM_IO.h"
#include "SB_Filter.h"
#include "IM3D_IO.h"
#include "Atrou3D.h"
#include "MR3D_Obj.h"

#define NBR_TRANS_3D 3
#define DEFAULT_BORDER_3D I_MIRROR
#define MAX_SCALE_1D 100

enum type_trans_3d {TO3_MALLAT,TO3_LIFTING, TO3_ATROUS, T3_UNDEFINED=-1};

enum set_trans_3d {TRANS3_MALLAT, TRANS3_PAVE, S3_UNDEFINED=-1};
inline char * StringTransf3D (type_trans_3d type)
{
    switch (type)
    {
        case TO3_MALLAT:
               return ((char*)"(bi-) orthogonal transform ");break;
               break;
        case TO3_LIFTING:
              return ((char*)"(bi-) orthogonal transform via lifting sheme"); 
              break;
	case TO3_ATROUS:
	      return ((char*)"A trous wavelet transform"); 
              break;
        case T3_UNDEFINED:
              return ((char*)"undefined transform");break;
    }
    return ((char*)"Error: bad type of transform");
}

inline char * StringSet3D (set_trans_3d type)
{
    switch (type)
    {
         case TRANS3_MALLAT:
              return ((char*)"non-redundant"); 
              break;
	case TRANS3_PAVE:
              return ((char*)"redundant"); 
              break;
        case  S3_UNDEFINED:
              return ((char*)"undefined type transform");break;
    }
    return ((char*)"Error: bad type of transform");
}

inline set_trans_3d which_set_is_trans3d(type_trans_3d type)
{
    switch (type)
    {
        case TO3_MALLAT:
        case TO3_LIFTING:
	      return TRANS3_MALLAT; break;
        case TO3_ATROUS:
	      return TRANS3_PAVE; break;
	case T3_UNDEFINED:
              return S3_UNDEFINED;break;
    }
    return S3_UNDEFINED;
} 

inline set_trans_3d SetTransform (type_trans_3d Transform)
{
    return which_set_is_trans3d(Transform);
}


/************************************************************************/

// 3D sub-band decomposision: an image is transformed into eight sub-cubes
class SubBand3D
{
	private:
		SubBand1D *Ptr_SB1D;
	public:
		SubBand3D(SubBand1D &SB1D) {Ptr_SB1D = &SB1D;}
		void transform3d(fltarray &Data);
		void transform3d(fltarray &Data, fltarray* Trans, int Step);
		void recons3d(fltarray  &Data);
		void recons3d(fltarray* Trans, fltarray &Data, int Step);
		~SubBand3D(){Ptr_SB1D = NULL;}
};

/***********************************************************************
*   Orthogonal 3D wavelet transform (decimated)
************************************************************************/

// Orthogonal 3D wavelet transform
class Ortho_3D_WT: public SubBand3D {
	public:
		Ortho_3D_WT(SubBand1D &SB1D):SubBand3D(SB1D) {};
		~Ortho_3D_WT() {};
		void transform (fltarray &Cube, int Nbr_Plan);
		void transform(fltarray &Cube_in, fltarray &Transf_Out, int Nbr_Plan);
		void recons(fltarray &Cube, int Nbr_Plan);
		void recons(fltarray &Transf_in, fltarray &Cube_Out, int Nbr_Plan);
};

/***********************************************************************
*   Orthogonal Undecimated 3D wavelet transform
************************************************************************/

class PAVE_3D_WT: public SubBand3D
{
	public:
		PAVE_3D_WT(SubBand1D &SB1D):SubBand3D(SB1D) {};
		~PAVE_3D_WT() {};
		int alloc (fltarray * & TabBand, int Nx, int Ny, int Nz, int NbrScale);
		void free(fltarray *TabBand, int NbrScale);
		void one_scale_transform (fltarray &Imag, fltarray *TabTrans, int Step, int Pos=0);
		void transform (fltarray &Cube, fltarray *TabTrans, int NbrScale);
		void recons(fltarray *TabTrans, fltarray &Cube, int NbrScale);
};

/************************************************************************/

class MR_3D {
		ATROUS_3D_WT AT3D_WT; // Wavelet transform class for the
	        		   // 3D wavelet transform
		fltarray Data; // buffer where the wavelet transform is stored
						// only used for non redundant transform
		fltarray *TabBand; // used for the redundant transform.
	        				// otherwise TabBand = &Data
		int Nbr_Plan;  // number of scales
		int Nbr_Band;  // number of band
		int Nx,Ny,Nz;  // number of points of the input signal

	    // for non redundant transform only
		int *TabPosX;  // TabPos and TabSize gives the position
		int *TabPosY;  // and the size of a given band.
		int *TabPosZ;  // used only for non-redundant transform
		int *TabSizeNx; // (i.e. size(signal) = size(transform))
		int *TabSizeNy;
		int *TabSizeNz;

		type_3d_format DataFormat; // input data format
		int  mr_io_fill_header(fitsfile *fptr);
        		 // write to the MR file header information

		FilterAnaSynt *FilterBank; // pointer to the filter bank to use
                        		// in the orthogonal wavelet transform
		Bool FilterBankAlloc;      // True if the filter bank is allocated
                        		// by the class   
		void init();
   public:
		Bool Verbose;
		type_trans_3d Type_Transform; // type of transform
		set_trans_3d  Set_Transform;  // class of transform
		type_border Border;           // type of border to user for the
                        		   // border intrpolation
		type_sb_filter SBFilter;       // Filter used in the subband 
	                    		// decomposition (case of Mallat transform).             
		sb_type_norm TypeNorm;         // type of normalization
	                    		// (in L1 or L2)
		type_lift LiftingTrans;        // type of lifting (in case of lifting
	                    		// transform					

		// return a pointer to the filter bank
		FilterAnaSynt * filter_bank() {return FilterBank;}

		MR_3D (){init();}
		MR_3D (int Nx, int Ny, int Nz, type_trans_3d T, int Nbr_Scale)
	    		   { alloc(Nx,Ny,Nz,T,Nbr_Scale);}

		inline int nbr_scale () const {return Nbr_Plan;}
		inline int nbr_band () const {return Nbr_Band;}
		int size_cube_nx () const {return Nx;}
		int size_cube_ny () const {return Ny;}
		int size_cube_nz () const {return Nz;}
		int size_band_nx (int b) const {return TabSizeNx[b];};
		int size_band_ny (int b) const {return TabSizeNy[b];};
		int size_band_nz (int b) const {return TabSizeNz[b];};
		void alloc (int Nx, int Ny, int Nz, type_trans_3d T, int Nbr_Scale,
			 FilterAnaSynt *FAS=NULL, sb_type_norm Norm=NORM_L1);

		// return the data set for non redundant transform
		fltarray & cube() {return Data;}

		// return the data set for  redundant transform
		fltarray * & TBand() {return TabBand;}

		float & operator() (int b, int x, int y, int z) const;

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

		void get_band(int b, fltarray &Band);
		// return in Band a band of the transform

		void insert_band(int b, fltarray &Band);
		// insert a Band a band of the transform

		void free ();
		// deallocate the object

		void read(char *FileName); // read the MR3D object from the disk
		void write(char *FileName); // write the MR3D object to the disk

		void info_pos_band(); // print information about size and position
	        		   // of each band
		void info_band(int b); // print statistical information
	            		// about the band b
		~MR_3D();
};

#endif
