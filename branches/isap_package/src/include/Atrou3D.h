/******************************************************************************
**                   Copyright (C) 2001 by CEA
*******************************************************************************
**
**    UNIT
**
**    Version: 1.0
**
**    Author: Jean-Luc Starck
**
**    Date:  15/10/01 
**    
**    File:  ATrou3D.h
**
*******************************************************************************
**
**    DESCRIPTION  3D wavelet transform using the a trous algorithm
**    -----------  
**                 
**    
******************************************************************************/

#ifndef _Atrou3D_H_
#define _Atrou3D_H_

#include <vector>
#include "GlobalInc.h"
#include "SB_Filter1D.h"
#include "IM_IO.h"

/************************************************************************/

class ATROUS_3D_WT {
	
private:
	int _Nx;
	int _Ny;
	int _Nz;
	int _NbrScale;
	inline int index(int i, int N) {return get_index(i, N, Bord);}

	bool no_fine;

public:
	// constructor & destructor
	ATROUS_3D_WT ();
	~ATROUS_3D_WT() {}

	Bool ModifiedAWT; // kept public for old functions, use set_use_modAWT(bool) instead
	Bool Adjoint; // kept public for old functions, use set_use_modAWT(bool) instead
	
	type_border Bord; // Type of border

	// get a pixel value, taking into account the border.
	float get_pix(fltarray & Cube, int x, int y, int z)
		{return (Cube(index(x,Cube.nx()), index(y,Cube.ny()), index(z,Cube.nz())));} 

	// allocate the memory space 	for a wavelet transform 
	void alloc(fltarray * & TabBand,int Nx,int Ny,int Nz, int NbrBand)
	{ 
		TabBand = new fltarray [NbrBand];
		for (int s = 0; s < NbrBand; s++) TabBand[s].alloc (Nx, Ny, Nz); 
	}

	// deallocate a wavelet transform.
	void free(fltarray * & TabBand, int Nbr_Plan)
		{ if (Nbr_Plan != 0) delete [] TabBand;}

	// B3 spline spline filtering of a cube, taking into accounts
	// the steps.
	void b3spline_filtering(fltarray & Old, fltarray & New, int s);

	// 3D wavelet transform
	void transform(fltarray & CubeIn, fltarray * & TabBandOut, int Nbr_Plan);

	// 3D wavelet reconstruction
	// Nbr_Plan==0 means it's already saved in the transform 
	void recons(fltarray * & TabBandIn, fltarray & CubeOut, int Nbr_Plan=0,
			Bool AddLastScale=True);

	void set_no_fine(bool nf) {no_fine=nf;}
	
	// added functions for mr3d_atrou
	dblarray TabStat;
	void set_use_modAWT(bool ad=true) {ModifiedAWT=(Bool)ad;}
	void set_use_adjoint(bool ad=true) {Adjoint=(Bool)ad;}
	void set_nbr_sale(int ns) {_NbrScale=ns;}
	void normalize_self(fltarray *TabBand, bool inverse=false);
	void threshold(fltarray * & TabBandIn, float thresh, bool soft=false, bool normalize=false);
	void clean_single(fltarray * & TabBandIn, float Sigma);
	void extract_stat(fltarray * & TabBandIn, char* Outname);
	void read(char *Name, fltarray * &TabBand, bool *NormalizeInv);
	void write(char *Name, fltarray * TabBand, bool Normalize);
};

void iwt3d_clear(std::vector< fltarray* > &C);
void iwt3d_transform(fltarray &Data, std::vector< fltarray* > &vTabBand, int NbrScale3D, bool modAWT);
void iwt3d_recons(std::vector< fltarray* > &vTabBand, fltarray &Data, bool modAWT, bool adjoint);
void iwt3d_threshold(std::vector< fltarray* > &vTabBand, bool modAWT, float threshold, filter_type FilterType);
void iwt3d_filter(fltarray &Data, fltarray &Recons, int NbrScale3D, bool modAWT, bool adjoint, float threshold, filter_type FilterType);

/************************************************************************/

#endif
