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
**    File:  MR1D1D.h
**
*******************************************************************************
**
**    DESCRIPTION  multiresolution transform of cube, using 2D wavelet transform
**    -----------  followed (or not) by a 1D wavelet transform
**                 
**
******************************************************************************/


#ifndef _MR1D1D_H_
#define _MR1D1D_H_

#include "IM_Obj.h"
#include "IM_IO.h"
#include "SB_Filter.h"
#include "MR1D_Obj.h"
#include "MR_Obj.h"
#include "IM3D_IO.h"
// #include "MR3D_Obj.h"
#include "DefFunc.h"

#define TESTIND 1

/****************************************************************************/

class MR1D1D {
 Bool Apply1DTrans;  // If True, a 1D Wavelet transform is applied on the z-axis
 int Nx, Ny;     // input cube size
 int NbrScaleX;     // Number of scales of the 2D multiscale transform.
 int NbrScaleY;     // Number of scales of the 1D multiscale transform (NbrScaleY=1 if no WT_y is applied)
 int NbrBandX;      // Number of bands in the 2D multiscale transform
 int NbrBandY;      // Number of bands in the 1D multiscale transform
 int NbrCoefX;      // Number of coefficients in the 2D multiscale transform
 int NbrCoefY;      // Number of coefficients in the 1D multiscale transform
  
  public:  
  
  MR_1D WT_x;    // 2D Multiresolution object
  MR_1D WT_y;         // 1D Multiresolution object
  fltarray *TabBand;  // Array of 3D  bands
  intarray TabSizeBandNx; // TabSizeBandNx[b2, b1] = size along x-axis of the (b2,b1) b2=band 2D, b1=band 1D
  intarray TabSizeBandNy; // TabSizeBandNy[b2, b1] = size along y-axis of the (b2,b1) b2=band 2D, b1=band 1D
  intarray TabFirstPosBandNy;
  
  // WT 1D
  sb_type_norm Norm;  // Type or normalisation of the 2D transform
  type_sb_filter SB_Filter; // Type of filter in case of 2D WT
  type_border Bord;         // Parameter for border management
  type_undec_filter U_Filter; // Type of filter in case of undecimated WT
  FilterAnaSynt FAS;          // Filter bank object
  
  int mr_io_fill_header(fitsfile *fptr);
  //public:
       Bool Verbose;
       MR1D1D (){ NbrBandX=NbrBandY=0;Verbose=False;}
       
       void alloc(int iNx, int iNy, type_trans_1d Trans1D, int Nsx, int Nsy=0, Bool NoAlloc=False);
       // Allocate the class for a cube of size (iNx, iNy, iNz) using Ns2D scale in 2D and Ns1D scale in 1D
       // If Ns1D < 2 then no wavelet transform is performed along z axis
       // If NoAlloc=True, then TabBand is not allocated and ONLY routine transform_to_vectarray can be used
       
       void transform (fltarray &Data);
       // Apply the 2D-1D wavelet transform. The transformation is stored in TabBand
       
       void recons (fltarray &Data);
       // Apply the 2D-1D reconstruction

       fltarray get_band(int sx, int sy=0);
       // Extract one band
       
       void put_band(fltarray Band, int sx, int sy=0);
       // Insert one band
       
       void write(char *Name);
       // Write the transformation to a file
       
       void read(char *Name);
       // Read a transformation from a file
       
       void transform_to_vectarray(fltarray &Data, fltarray &TabVect);
       // Apply the 2D-1D transfrom, but the transformation is not stored in TabBand, and
       // reconstruction is possible
       //  get_band, put_band, read, write routines cannot be used
       
       int nbr_band_x () const { return NbrBandX;}
       // Return the number of bands of the 2D multiresolution transform
       
       int nbr_band_y () const { return NbrBandY;}
       // Return the number of bands of the 1D wavelet transform
       
       int size_band_nx(int sx, int sy=0) const { return TabSizeBandNx(sx, sy);}
       // Return the size of the band (s2d,s1d) along x-axis
       
       int size_band_ny(int sx, int sy=0) const { return TabSizeBandNy(sx, sy);}
       // Return the size of the band (s2d,s1d) along y-axis
       
       int nbr_coef_x()  const { return NbrCoefX;}
       // Return the number of coefficients of the 2D decomposition
       
       int nbr_coef_y()  const { return NbrCoefY;}
       // Return the number of coefficients of the 1D decomposition
       
       int nbr_pix_nx()  const { return Nx;}
       // Return the number of coefficients of the 1D decomposition
       
       int nbr_pix_ny()  const { return Ny;}
       // Return the number of coefficients of the 1D decomposition

       float & operator() (int s_x, int s_y, int i, int j) const;
       // Return one coefficient
         // Return one coefficient, for the case where no 1D WT is applied

       ~MR1D1D() { if (TabBand != NULL) delete [] TabBand; NbrBandY=NbrBandX=NbrScaleY=NbrScaleX=0;}
};

/****************************************************************************/
 #endif
