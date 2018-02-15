/******************************************************************************
**                   Copyright (C) 2007 by CEA
*******************************************************************************
**
**    UNIT
**
**    Version: 1.0
**
**    Author: Jean-Luc Starck
**
**    Date:  14/06/2007
**    
**    File:  MR2D1D.h
**
*******************************************************************************
**
**    DESCRIPTION  multiresolution transform of cube, using 2D wavelet transform
**    -----------  followed (or not) by a 1D wavelet transform
**                 
**
******************************************************************************/


#ifndef _MR2D1D_H_
#define _MR2D1D_H_

#include "IM_Obj.h"
#include "IM_IO.h"
#include "SB_Filter.h"
#include "MR1D_Obj.h"
#include "MR_Obj.h"
#include "IM3D_IO.h"
#include "MR3D_Obj.h"
#include "DefFunc.h"

#define TESTIND 1

/****************************************************************************/

class MR2D1D {
  Bool Apply1DTrans;  // If True, a 1D Wavelet transform is applied on the z-axis
  int Nx, Ny, Nz;     // input cube size
  int NbrScale2D;     // Number of scales of the 2D multiscale transform.
  int NbrScale1D;     // Number of scales of the 1D multiscale transform (NbrScale1D=1 if no WT1D is applied)
  int NbrBand2D;      // Number of bands in the 2D multiscale transform
  int NbrBand1D;      // Number of bands in the 1D multiscale transform
  int NbrCoef2D;      // Number of coefficients in the 2D multiscale transform
  int NbrCoef1D;      // Number of coefficients in the 1D multiscale transform
  
  public: // ATTENTION JE L'AI MIS ICI ET NON PLUS LOIN !!
  
  MultiResol WT2D;    // 2D Multiresolution object
  MR_1D WT1D;         // 1D Multiresolution object
  fltarray *TabBand;  // Array of 3D  bands
  intarray TabFirstPosBandNz; 
           // TabBand[b2](Nx, Ny, TabFirstPosBandNz[b1]:TabFirstPosBandNz[b1]:TabSizeBandNz[b2, b1]-1)
	   // contains the band (b2,b1) b2=band 2D, b1=band 1D
  intarray TabSizeBandNx; // TabSizeBandNx[b2, b1] = size along x-axis of the (b2,b1) b2=band 2D, b1=band 1D
  intarray TabSizeBandNy; // TabSizeBandNy[b2, b1] = size along y-axis of the (b2,b1) b2=band 2D, b1=band 1D
  intarray TabSizeBandNz; // TabSizeBandNz[b2, b1] = size along z-axis of the (b2,b1) b2=band 2D, b1=band 1D
  
  // WT 2D
  sb_type_norm Norm;  // Type or normalisation of the 2D transform
  type_sb_filter SB_Filter; // Type of filter in case of 2D WT
  type_border Bord;         // Parameter for border management
  type_undec_filter U_Filter; // Type of filter in case of undecimated WT
  FilterAnaSynt FAS;          // Filter bank object
  
  int mr_io_fill_header(fitsfile *fptr);
  //public:
       Bool Verbose;
       MR2D1D (){ NbrBand2D=NbrBand1D=0;Verbose=False;}
       
       void alloc(int iNx, int iNy, int iNz, type_transform Trans2D, int Ns2D, int Ns1D=0, Bool NoAlloc=False);
       // Allocate the class for a cube of size (iNx, iNy, iNz) using Ns2D scale in 2D and Ns1D scale in 1D
       // If Ns1D < 2 then no wavelet transform is performed along z axis
       // If NoAlloc=True, then TabBand is not allocated and ONLY routine transform_to_vectarray can be used
       
       void transform (fltarray &Data);
       // Apply the 2D-1D wavelet transform. The transformation is stored in TabBand
       
       void recons (fltarray &Data);
       // Apply the 2D-1D reconstruction

       fltarray get_band(int s2, int s1=0);
       // Extract one band
       
       void put_band(fltarray Band, int s2, int s1=0);
       // Insert one band
       
       void write(char *Name);
       // Write the transformation to a file
       
       void read(char *Name);
       // Read a transformation from a file
       
       void transform_to_vectarray(fltarray &Data, fltarray &TabVect);
       // Apply the 2D-1D transfrom, but the transformation is not stored in TabBand, and
       // reconstruction is possible
       //  get_band, put_band, read, write routines cannot be used
       
       int nbr_band_2d () const { return NbrBand2D;}
       // Return the number of bands of the 2D multiresolution transform
       
       int nbr_band_1d () const { return NbrBand1D;}
       // Return the number of bands of the 1D wavelet transform
       
       int size_band_nx(int s2d, int s1d=0) const { return TabSizeBandNx(s2d, s1d);}
       // Return the size of the band (s2d,s1d) along x-axis
       
       int size_band_ny(int s2d, int s1d=0) const { return TabSizeBandNy(s2d, s1d);}
       // Return the size of the band (s2d,s1d) along y-axis
       
       int size_band_nz(int s2d, int s1d=0) const { return TabSizeBandNz(s2d, s1d);}
       // Return the size of the band (s2d,s1d) along z-axis
       
       int nbr_coef_2d()  const { return NbrCoef2D;}
       // Return the number of coefficients of the 2D decomposition
       
       int nbr_coef_1d()  const { return NbrCoef1D;}
       // Return the number of coefficients of the 1D decomposition
       
       int nbr_pix_nx()  const { return Nx;}
       // Return the number of coefficients of the 1D decomposition
       
       int nbr_pix_ny()  const { return Ny;}
       // Return the number of coefficients of the 1D decomposition

       float & operator() (int s2, int s1, int i, int j, int k) const;
       // Return one coefficient
       
       float & operator() (int s2, int i, int j, int k) const;
       // Return one coefficient, for the case where no 1D WT is applied

       ~MR2D1D() { if (TabBand != NULL) delete [] TabBand; NbrBand1D=NbrBand2D=NbrScale1D=NbrScale2D=0;}
};

/****************************************************************************/
 #endif
