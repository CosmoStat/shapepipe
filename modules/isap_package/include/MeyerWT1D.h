/******************************************************************************
**                   Copyright (C) 2008 by CEA
*******************************************************************************
**
**    UNIT
**
**    Version: 1.0
**
**    Author: Arnaud Woiselle
**
**    Date:  30/09/2008
**    
**    File:  MeyerWT.h
**
**    Modification history :
**
*******************************************************************************
**
**    DESCRIPTION:  1D MEYER WAVELET TRANSFORM 
**    ----------- 
******************************************************************************/

#ifndef _MEYER_WT1D
#define _MEYER_WT1D

#include "GlobalInc.h"
#include "FFTN.h"
#include "FFTN_1D.h"
#include "IM_IO.h"
// #include "MGA_Inc.h"

/***********************************************************************/

class MEYER_WT1D
{
protected:
	// Function used to compute the low pass filtering
	fltarray H;  // H filter
	inline double lowpass_window(double x)
	{
		double l,r;
		l = exp( 1-1 / (1-exp(1-1/(1-x))));
		r = exp(1-1/(1-exp(1-1/x)));
		return (l /= sqrt(l*l+r*r));
	}

	int Nx_Cube;     // Input cube x size
	int NbrScale;    // Nbr scales in the WT

	cfarray TF_ExtData;	// Fourier transform of the extended cube
						// If Extend = False then TF_ExtData = FFT(input_data)
	cfarray *Tabcf_WT_Band; // Wavelet transform of the cube in Fourier domain
	
// Statistic elements
	double *TabStat;	// Std of each subband
	double *MaxCoef;
	
	intarray TabNx;	// TabNx[s] = x size at scale s

	Bool NeedOddSize;	// If NeedOddSize==True, all scales will have odd sizes
                		// In this case, the input cube must have an odd size

	int ExtNx; 	// size of the signal to be transformed
				// = 2*floor(D_ExtNx/2)+1

	double D_ExtNx; // D_ExtNx = Nx if Extend==False and 4./3. Nx otherwise

	void get_hfilter(fltarray &H, double DNx=0., double DNy=0., double DNz=0.);
	// Calculate the low pass filtering in the Fourier domain to be 
	// applied at each step of the wavelet decomposition.

	void get_extFourier(fltarray &Data, cfarray & TF_ExtData);
	//  Take the FFT of the input cube: Data
	//  Extend it by periodization if Extend == True
	//  TF_ExtData contains the FFT of the [extended] input cube.

	void get_IFFT_Cube(cfarray & TF_ExtData, fltarray &Data);
	// Take the inv. FFT of the input  TF_ExtData,
	// with extraction of part of interest if Extend == True

	void ifft_tabcube(cfarray * & TabCF_Cube, fltarray * & Tab_Cube, Bool Alloc=True);
	// Take the inverse FFT of each wavelet scale
	// If Alloc == True, Tab_Cube is allocated
	//                 otherwise it must have been allocated before  the call
	// IN: TabCF_Cube[s] = Fourier transform of the scale s
	// OUT: Tab_Cube[s] = scale s in direct space

	void fft_tabcube(fltarray * & Tab_Cube, cfarray * & TabCF_Cube); 
	// Take the Forward  FFT of each wavelet scale
	// IN: Tab_Cube[s] = scale s in direct space
	// OUT: TabCF_Cube[s] = Fourier transform of the scale s

	void transform_cf(cfarray & TF_ExtData, cfarray * & TabWT);
	// Decomposition in Fourier space of the input TF_ExtData
	// TabWT: output: = TabWT[s] = Fourier transform of the scale s

	void recons_cf(cfarray * &TabWT,  cfarray & TF_ExtData);
	// Reconstruction the Fourier transform of an cube from
	// the Fourier transforms of its scales.
	// TabW: input: = TabWT[s] = Fourier transform of the scale s
	// TF_ExtData: output reconstructed cube in Fourier space

public:
	MEYER_WT1D();
	~MEYER_WT1D();
	FFTN_1D FFT1D;   // public for test purposes, else protected
    
    bool Verbose;

// Initialize the WT class for a given number of scales and a give cube size
	void init(int Nbr_Scale, int Nx, Bool WTNeedOddSize=False);

// Computes the wavelet transform  of a cube
	void transform(fltarray &Data);
	void transform(fltarray &Data, fltarray * & Tab_Cube, Bool Alloc=True);
	// Data: in = input cube
	// TabWT: out = TabWT[s] is the scale s of the WT
	// If Alloc == True then Tab_Cube is allocated

// Reconstruction of a cube from its WT.
	void recons(fltarray &Data);
	void recons(fltarray * & Tab_Cube, fltarray &Data, Bool Alloc=True);
	// Tab_Cube: in = input wavelet scales
	// Data: out = output reconstructed  cube
	// If Alloc == True then Data is allocated

// Size and position tools
	inline int nxcube() {return Nx_Cube;}
	inline double nxd() {return D_ExtNx;}
	inline int nxs(int s) {return TabNx(s);}
	int nbr_scale() {return NbrScale;} // Number of scales of the WT

// Statistic and information tools
	void extract_stat(fltarray *TabBand, char* Outname);
	void values_at(fltarray *TabBand, char * filename, char* Outname); // use after normalize_self, to extract values at filename's coef list
	void noise_calibration(fltarray *TabBand, char* Outname);
	void normalize(fltarray *TabBand, fltarray *TabBandNorm );
	void normalize_self(fltarray *TabBand, bool inverse=false);

// Filtering methods
	void threshold(fltarray *TabBand, float SigmaNoise, float NSigma, filter_type FilterType=FT_HARD, bool force3sigma=false);
	void wiener(fltarray *TabBand, float noise_lvl, int LocalBS);
	
// IO tools
	void write(char *Name, fltarray * & Tab_WCube, bool Normalize);
	void read(char *Name, fltarray * & Tab_WCube, bool *NormalizeInv);
	
	void test(fltarray *TabBand);
};

/***********************************************************************/
//	void write_mono(char *Name, fltarray * & Tab_WCube);
//	void read_mono(char *Name, fltarray * & Tab_WCube);

#endif
