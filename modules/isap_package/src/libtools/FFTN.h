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
**    Date:  28/07/00 
**    
**    File:  FFTN.h
**
*******************************************************************************
**
**    DESCRIPTION: FFT Class for 1D,2D,3D FFT with array of any size.  
**    ----------- 
**                 
**
******************************************************************************/

#ifndef _DEF_FFTN_H_
#define _DEF_FFTN_H_

#include "GlobalInc.h"

// void fft_free (void);

/* double precision routine */
int fftn (int ndim, const int dims[], double Re[], double Im[],
                 int isign, double scaling);

/* float precision routine */
int fftnf (int ndim, const int dims[], float Re[], float Im[],
                  int isign, double scaling);

#define NFACTOR	11

class FFTN
{
		void fft_error() 
		{
		   cout << "Error in FFT computation ... " << endl;
		   exit(-1);
		}
		void scale_reverse(float *Buff, int N)
		{
		   double C = 1. / (double) N;
		   for (int i=0; i < 2*N; i++) Buff[i] *= C;
		}
		void scale_reverse(double *Buff, int N)
		{
		   double C = 1. / (double) N;
		   for (int i=0; i < 2*N; i++) Buff[i] *= C;
		}
	public:
		
	// Utilities
		int factor [NFACTOR];
		/* temp space, (void *) since both float and double routines use it */
		void *Tmp0;  /* temp space for real part */
		void *Tmp1;  /* temp space for imaginary part */
		void *Tmp2;  /* temp space for Cosine values */
		void *Tmp3;  /* temp space for Sine values */
		int  *Perm;  /* Permutation vector */

		// double precision routines
		int fftn (int ndim, const int dims[], double Re[], double Im[], int isign, double scaling);
		int fftradix (double Re[], double Im[], size_t nTotal, size_t nPass, size_t nSpan, int isign, int max_factors, int max_perm);
		// float precision routines
		int fftnf (int ndim, const int dims[], float Re[], float Im[], int isign, double scaling);
		int fftradixf (float Re[], float Im[], size_t nTotal, size_t nPass, size_t nSpan, int isign, int max_factors, int max_perm);
		// parameters for memory management 
		size_t SpaceAlloced;
		size_t MaxPermAlloced;

	// Transform functions
		void transform1d(float *CFBuff, int N, Bool Reverse=False, bool normalize=false);
		// Complex 1D FFT transform
		// N: IN = number of elements
		// CFBuff: IN-OUT = buffer [0..2N-1]
		//                   CFBuff[2i] = real part
		//                   CFBuff[2i+1] = imaginary part
		// Reverse: IN = True for inverse FFT transform

		void transform1d(double *CFBuff, int N, Bool Reverse=False, bool normalize=false);
		// Ditto with double

		void transform2d(float *CFBuff, int Nx, int Ny, Bool Reverse=False, bool normalize=false);
		// Complex 2D FFT transform
		// Nx: IN = number of columns
		// Ny: IN = number of lines
		// CFBuff: IN-OUT = buffer [0..N-1] with N = 2NxNy 
		//                   CFBuff[2i] = real part
		//                   CFBuff[2i+1] = imaginary part
		// Reverse: IN = True for inverse FFT transform

		void transform2d(double *CFBuff, int Nx, int Ny, Bool Reverse=False, bool normalize=false);
		// Ditto with double

		void transform3d(float *CFBuff, int Nx, int Ny, int Nz, Bool Reverse=False, bool normalize=false);
		// Complex 2D FFT transform
		// Nx,Ny,Nz: IN = array dimensions
		// CFBuff: IN-OUT = buffer [0..N-1] with N = 2NxNyNz
		//                   CFBuff[2i] = real part
		//                   CFBuff[2i+1] = imaginary part
		// Reverse: IN = True for inverse FFT transform

		void transform3d(double *CFBuff, int Nx, int Ny, int Nz, Bool Reverse=False, bool normalize=false);
		// Ditto with double

		FFTN();
		~FFTN();
};


#endif
