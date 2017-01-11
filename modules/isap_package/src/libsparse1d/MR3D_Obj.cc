/******************************************************************************
**                   Copyright (C) 1999 by CEA
*******************************************************************************
**
**    UNIT
**
**    Version: 1.0
**
**    Author: Jean-Luc Starck
**
**    Date:  04/09/99
**    
**    File:  mr3d_trans.cc
**
*******************************************************************************
**
**    DESCRIPTION  multiresolution transform of cube
**    ----------- 
**                 
**    Usage: mr3d_trans options cube output
**		
**		Append :
**			07/2008 : Arnaud Woiselle
**					undecimated transforms : PAVE_3D_WT
**		
******************************************************************************/

#include "MR3D_Obj.h"
 
/************************************************************************/
// 3D sub-band decomposision: an image is transformed into eight sub-cubes

void SubBand3D::transform3d (fltarray &Data)
{
//	cerr<<"SubBand3D::transform..."<<endl;
    int Nx = Data.nx();
    int Ny = Data.ny();
    int Nz = Data.nz();
    int Nz2 = (Nz+1)/2;
    SubBand2D SB2D(*Ptr_SB1D);
    Ifloat ImaFrame;
    int j, i, k;
    float *PtrHigh = new float [Nz];
    float *PtrLow  = new float[Nz];
    float *PtrDet  = new float[Nz];
       
    float *Cube = Data.buffer();
    float *Frame = Cube;
    
    // Transform each frame of the cube 
    for (k=0; k<Nz; k++)
    {
       ImaFrame.alloc(Frame,Ny,Nx);
       SB2D.transform2d(ImaFrame);
       Frame += Nx*Ny;
       ImaFrame.free();
    }
   
    // each z-vector is convolved by h and g
    for (i = 0; i < Nx; i++)
    for (j = 0; j < Ny; j++)
    {
       for (k=0; k<Nz; k++) PtrHigh[k] = Data(i,j,k);
       Ptr_SB1D->transform(Nz, PtrHigh, PtrLow, PtrDet);
       for (k=0; k< Nz2; k++) Data(i,j,k) = PtrLow[k];
       for (k=0; k< Nz/2; k++) Data(i,j,k+Nz2) = PtrDet[k];
    }
    
    delete [] PtrHigh;
    delete [] PtrLow;
    delete [] PtrDet;
//	cerr<<"End SubBand3D::transform"<<endl;
}
			
/************************************************************************/

void SubBand3D::transform3d(fltarray &Data, fltarray* TabTrans, int Step)
{
//	cerr<<"SubBand3D::transform3d_redondant..."<<endl;
	int Nx = Data.nx();
	int Ny = Data.ny();
	int Nz = Data.nz();
	SubBand2D SB2D(*Ptr_SB1D);
	
	// Data
	Ifloat CubeFrame;
	float *CubePtr = Data.buffer();
	
	// 2d transformed version
	fltarray Horiz(Nx,Ny,Nz);
	fltarray Vert(Nx,Ny,Nz);
	fltarray Diag(Nx,Ny,Nz);
	fltarray Smooth(Nx,Ny,Nz);
	
	// Pointers associated
	float *HorizPtr = Horiz.buffer();
	float *VertPtr = Vert.buffer();
	float *DiagPtr = Diag.buffer();
	float *SmoothPtr = Smooth.buffer();
	
	// virtual Ifloat associated
	Ifloat HorizFrame;
	Ifloat VertFrame;
	Ifloat DiagFrame;
	Ifloat SmoothFrame;
	
	int j, i, k;
	float *PtrHigh = new float[Nz];
	float *PtrLow  = new float[Nz];
	float *PtrDet  = new float[Nz];

	// Transform each frame of the cube 
	for (k=0; k<Nz; k++)
	{
		HorizFrame.alloc(HorizPtr,Ny,Nx);
		VertFrame.alloc(VertPtr,Ny,Nx);
		DiagFrame.alloc(DiagPtr,Ny,Nx);
		SmoothFrame.alloc(SmoothPtr,Ny,Nx);
		
		CubeFrame.alloc(CubePtr,Ny,Nx);

		SB2D.transform2d (CubeFrame, HorizFrame, VertFrame, DiagFrame, SmoothFrame, Step);
		
		CubePtr += Nx*Ny;
		CubeFrame.free();
		
		HorizPtr += Nx*Ny;
		VertPtr += Nx*Ny;
		DiagPtr += Nx*Ny;
		SmoothPtr += Nx*Ny;
	}
	
	// Smooth => s z
	// Horiz => x xz
	// Diag => xy xyz
	// Vert => y yz
	// bands : 0   1   2   3   4   5   6   7 
	//  xyz : HHH HHL HLH HLL LHH LHL LLH LLL (high/low)
	//        DH  DL  HH  HL  VH  VL  SH  SL
	
	// each z-vector is convolved by h and g
	for (i = 0; i < Nx; i++)
	for (j = 0; j < Ny; j++)
	{
		for (k=0; k<Nz; k++) PtrHigh[k] = Smooth(i,j,k);
		Ptr_SB1D->transform(Nz, PtrHigh, PtrLow, PtrDet, Step);
		for (k=0; k<Nz; k++) (TabTrans[7])(i,j,k) = PtrLow[k];
		for (k=0; k<Nz; k++) (TabTrans[6])(i,j,k) = PtrDet[k];
		
		for (k=0; k<Nz; k++) PtrHigh[k] = Vert(i,j,k);
		Ptr_SB1D->transform(Nz, PtrHigh, PtrLow, PtrDet, Step);
		for (k=0; k<Nz; k++) (TabTrans[5])(i,j,k) = PtrLow[k];
		for (k=0; k<Nz; k++) (TabTrans[4])(i,j,k) = PtrDet[k];

		for (k=0; k<Nz; k++) PtrHigh[k] = Horiz(i,j,k);
		Ptr_SB1D->transform(Nz, PtrHigh, PtrLow, PtrDet, Step);
		for (k=0; k<Nz; k++) (TabTrans[3])(i,j,k) = PtrLow[k];
		for (k=0; k<Nz; k++) (TabTrans[2])(i,j,k) = PtrDet[k];
		
		for (k=0; k<Nz; k++) PtrHigh[k] = Diag(i,j,k);
		Ptr_SB1D->transform(Nz, PtrHigh, PtrLow, PtrDet, Step);
		for (k=0; k<Nz; k++) (TabTrans[1])(i,j,k) = PtrLow[k];
		for (k=0; k<Nz; k++) (TabTrans[0])(i,j,k) = PtrDet[k];
	}

	delete [] PtrHigh;
	delete [] PtrLow;
	delete [] PtrDet;
//	cerr<<"End SubBand3D::transform3d_redondant"<<endl;
}
/************************************************************************/

void SubBand3D::recons3d (fltarray &Data)
{
    int Nx = Data.nx();
    int Ny = Data.ny();
    int Nz = Data.nz();
    int Nz2 = (Nz+1)/2;
    SubBand2D SB2D(*Ptr_SB1D);
    Ifloat ImaFrame;
    int j, i, k;
    float *PtrHigh = new float[Nz];
    float *PtrLow  = new float[Nz];
    float *PtrDet  = new float[Nz];
       
    float *Cube = Data.buffer();
    float *Frame = Cube;

    // each z-vector is convolved by h and g
    for (i = 0; i < Nx; i++)
    for (j = 0; j < Ny; j++)
    {
       for (k=0; k< Nz2; k++)  PtrLow[k] = Data(i,j,k);
       for (k=0; k< Nz/2; k++) PtrDet[k] = Data(i,j,k+Nz2);
       Ptr_SB1D->recons(Nz,PtrLow, PtrDet, PtrHigh);
       for (k=0; k<Nz; k++) Data(i,j,k) = PtrHigh[k];
    }
        
    // Transform each frame of the cube 
    for (k=0; k<Nz; k++)
    {
       ImaFrame.alloc(Frame,Ny,Nx);
       SB2D.recons2d(ImaFrame);
       Frame += Nx*Ny;
       ImaFrame.free();
    }
   
    delete [] PtrHigh;
    delete [] PtrLow;
    delete [] PtrDet;
}	
/************************************************************************/

void SubBand3D::recons3d (fltarray* TabTrans, fltarray &Data, int Step)
{
//	cerr<<"SubBand3D::recons3d_redondant..."<<endl;
	int Nx = TabTrans[0].nx();
	int Ny = TabTrans[0].ny();
	int Nz = TabTrans[0].nz();
	SubBand2D SB2D(*Ptr_SB1D);
	
	// 2d transformed version
	fltarray Horiz(Nx,Ny,Nz);
	fltarray Vert(Nx,Ny,Nz);
	fltarray Diag(Nx,Ny,Nz);
	fltarray Smooth(Nx,Ny,Nz);
	
	// Pointers associated
	float *HorizPtr = Horiz.buffer();
	float *VertPtr = Vert.buffer();
	float *DiagPtr = Diag.buffer();
	float *SmoothPtr = Smooth.buffer();
	
	// virtual Ifloat associated
	Ifloat HorizFrame;
	Ifloat VertFrame;
	Ifloat DiagFrame;
	Ifloat SmoothFrame;
	
	int j, i, k;
	float *PtrHigh = new float[Nz];
	float *PtrLow  = new float[Nz];
	float *PtrDet  = new float[Nz];


	for (i = 0; i < Nx; i++)
	for (j = 0; j < Ny; j++)
	{
		for (k=0; k<Nz; k++) PtrLow[k] = (TabTrans[7])(i,j,k);
		for (k=0; k<Nz; k++) PtrDet[k] = (TabTrans[6])(i,j,k);
		Ptr_SB1D->recons(Nz, PtrLow, PtrDet, PtrHigh, Step);
		for (k=0; k<Nz; k++) Smooth(i,j,k) = PtrHigh[k];
		
		for (k=0; k<Nz; k++) PtrLow[k] = (TabTrans[5])(i,j,k);
		for (k=0; k<Nz; k++) PtrDet[k] = (TabTrans[4])(i,j,k);
		Ptr_SB1D->recons(Nz, PtrLow, PtrDet, PtrHigh, Step);
		for (k=0; k<Nz; k++) Vert(i,j,k) = PtrHigh[k];

		for (k=0; k<Nz; k++) PtrLow[k] = (TabTrans[3])(i,j,k);
		for (k=0; k<Nz; k++) PtrDet[k] = (TabTrans[2])(i,j,k);
		Ptr_SB1D->recons(Nz, PtrLow, PtrDet, PtrHigh, Step);
		for (k=0; k<Nz; k++) Horiz(i,j,k) = PtrHigh[k];
		
		for (k=0; k<Nz; k++) PtrLow[k] = (TabTrans[1])(i,j,k);
		for (k=0; k<Nz; k++) PtrDet[k] = (TabTrans[0])(i,j,k);
		Ptr_SB1D->recons(Nz, PtrLow, PtrDet, PtrHigh, Step);
		for (k=0; k<Nz; k++) Diag(i,j,k) = PtrHigh[k];
	}

	// Data : erase the previously finest scale, now saved in Diag
	Data.alloc(Nx,Ny,Nz);
	Data.init(0);
	Ifloat CubeFrame;
	float *CubePtr = Data.buffer();
	
	
	for (k=0; k<Nz; k++)
	{
		HorizFrame.alloc(HorizPtr,Ny,Nx);
		VertFrame.alloc(VertPtr,Ny,Nx);
		DiagFrame.alloc(DiagPtr,Ny,Nx);
		SmoothFrame.alloc(SmoothPtr,Ny,Nx);
		
		CubeFrame.alloc(CubePtr,Ny,Nx);

		SB2D.recons2d (HorizFrame, VertFrame, DiagFrame, SmoothFrame, CubeFrame, Step);

		CubePtr += Nx*Ny;
		CubeFrame.free();
		
		HorizPtr += Nx*Ny;
		VertPtr += Nx*Ny;
		DiagPtr += Nx*Ny;
		SmoothPtr += Nx*Ny;
	}

    delete [] PtrHigh;
    delete [] PtrLow;
    delete [] PtrDet;
//	cerr<<"End SubBand3D::recons3d_redondant"<<endl;
}	


//****************************************************************************
//            Ortho_3D_WT OBJ : Orthogonal decimated 3D wavelet transform
//**************************************************************************** 

void Ortho_3D_WT::transform (fltarray  &Cube, fltarray &Cube_Out, int NbrScale)
{
   // copy the image in Ima_Aux
   Cube_Out = Cube;
   transform(Cube_Out, NbrScale);
}

/************************************************************************/

void Ortho_3D_WT::transform (fltarray &Data, int NbrScale)
{
   int i,j,k,s;
   int Nx = Data.nx();
   int Ny = Data.ny();
   int Nz = Data.nz();
    
   // Apply a first transform on the cube itself  
   // it allows to save memory for the Data_Aux size
   transform3d (Data);

   Nx = (Nx+1)/2;
   Ny = (Ny+1)/2;
   Nz = (Nz+1)/2;
 
   for (s = 0; s < NbrScale-2; s++)
   {
       fltarray Data_Aux(Nx,Ny,Nz);
       for (i = 0; i < Nx; i++)
       for (j = 0; j < Ny; j++)
       for (k = 0; k < Nz; k++) Data_Aux(i,j,k) = Data(i,j,k);
          
       // transform Data_Aux
       transform3d (Data_Aux);
       
       for (i = 0; i < Nx; i++)
       for (j = 0; j < Ny; j++)
       for (k = 0; k < Nz; k++) Data(i,j,k) = Data_Aux(i,j,k);
       Nx = (Nx+1)/2;
       Ny = (Ny+1)/2;
       Nz = (Nz+1)/2;
    }
}

/************************************************************************/

void Ortho_3D_WT::recons (fltarray &Transf_in, fltarray &Cube_out, int Nbr_Plan) 
{
   Cube_out = Transf_in;
   recons(Cube_out, Nbr_Plan);
   
}

/************************************************************************/

void Ortho_3D_WT::recons (fltarray &Data, int Nbr_Plan) 
{
   int i,j,k,s;
   int Nx = Data.nx();
   int Ny = Data.ny();
   int Nz = Data.nz();
   
    /* last scale size */
   int  Nx_2, Ny_2, Nz_2, Nxs, Nys, Nzs;
 
   for (s = Nbr_Plan-2; s >= 0; s--)
   {
        Nxs = size_resol(s, Nx);
        Nys = size_resol(s, Ny);
        Nzs = size_resol(s, Nz);
	
        Nx_2 = (Nxs+1)/2;
        Ny_2 = (Nys+1)/2;
	Nz_2 = (Nzs+1)/2;
        fltarray Data_Aux(Nxs,Nys,Nzs);

        // copy the transform part in Imag
        for (i = 0; i < Nxs; i++)
        for (j = 0; j < Nys; j++)
	for (k = 0; k < Nzs; k++) Data_Aux(i,j,k) = Data(i,j,k);
       	
	// inverse transform Ima_Aux
        recons3d(Data_Aux);

        for (i = 0; i < Nxs; i++)
        for (j = 0; j < Nys; j++)
	for (k = 0; k < Nzs; k++) Data(i,j,k) = Data_Aux(i,j,k);         
   }
    
}

/************************************************************************/
//            PAVE_3D_WT OBJ 
/****************************************************************************/

int PAVE_3D_WT::alloc (fltarray * & TabBand, int Nx, int Ny, int Nz, int NbrScale)
{
	int s,NbrBand_per_Resol = 7;
	int NbrBand = NbrBand_per_Resol*(NbrScale-1)+1;
	TabBand = new fltarray [NbrBand];
	for (s = 0; s < NbrBand; s++)
		TabBand[s].alloc(Nx,Ny,Nz);
	return NbrBand;
}

/****************************************************************************/

void PAVE_3D_WT::free(fltarray *TabBand, int NbrScale)
{
    if (NbrScale != 0) delete [] TabBand;
}

/****************************************************************************/

void PAVE_3D_WT::transform (fltarray &Cube, fltarray *TabTrans, int NbrScale)
{
//	cerr<<"PAVE_3D_WT::transform..."<<endl;
	for (int s = 0; s < NbrScale-1; s++)
	{
		int Step = POW2(s);
		if (s == 0)
			transform3d(Cube, TabTrans, Step);
		else
			transform3d(TabTrans[7*s], &(TabTrans[7*s]), Step);
	}
//	cerr<<"End PAVE_3D_WT::transform"<<endl;
}

/****************************************************************************/

void PAVE_3D_WT::one_scale_transform (fltarray &Imag, fltarray *TabTrans, int Step, int Pos)
{
    transform3d (Imag, &(TabTrans[Pos]), Step);
}

/****************************************************************************/

void PAVE_3D_WT::recons (fltarray *TabTrans, fltarray &Cube, int NbrScale)
{
//	cerr<<"PAVE_3D_WT::recons..."<<endl;

	for (int s = NbrScale-2; s >= 0; s--)
	{
		int Step = POW2(s);
		if (s == 0)
			recons3d(TabTrans, Cube, Step);
		else
			recons3d(&(TabTrans[7*s]), TabTrans[7*s], Step);
	}
//	cerr<<"End PAVE_3D_WT::recons"<<endl;
}


/************************************************************************/
//            3D MR OBJ 
/****************************************************************************/

float & MR_3D::operator() (int b, int i, int j, int k) const 
{
   if ((b < 0) || (b >= nbr_band())
       || (i < 0) || (i >= size_band_nx(b))
       || (j < 0) || (j >= size_band_ny(b))
       || (k < 0) || (k >= size_band_nz(b)))
   {
      cout << "Error: band coefficient index ... " << endl;
      cout << "       Band number = " << b << " Nb = " << nbr_band() << endl;
      cout << "       X pos = " << i  << " Nx = " <<  size_band_nx(b) << endl;
      cout << "       Y pos = " << j  << " Ny = " <<  size_band_ny(b) << endl;
      cout << "       Z pos = " << k  << " Nz = " <<  size_band_nz(b) << endl;
      exit(-1);
   }       
   if (Set_Transform == TRANS3_PAVE)
   {
      return TabBand[b](i,j,k);
   }
   else
   {    
   int Indi = TabPosX[b] + i;
   int Indj = TabPosY[b] + j;
   int Indk = TabPosZ[b] + k;
   
   if ((Indi < 0) || (Indi >= Nx) || 
       (Indj < 0) || (Indj >= Ny) || (Indk < 0) || (Indk >= Nz))
   {
      cout << "Error: band coefficient index ... " << endl;
      cout << "       Band number = " << b << " Nb = " << nbr_band() << endl;
      cout << "       Indi = " << Indi  << " Nx = " <<  Nx << endl;
      cout << "       Indj = " << Indj  << " Ny = " <<  Ny << endl;
      cout << "       Indk = " << Indk  << " Nz = " <<  Nz << endl;
      exit(-1);
   }   
   return Data(Indi,Indj,Indk);
   }
}

/****************************************************************************/


void set_size_band(int N, Bool FilterLow, int & P, int &Ns)
{
   int L = (N+1)/2;
    
   if (FilterLow == True)
   {
      Ns = L;
      P = 0;
   }
   else
   {
      Ns = N/2;
      P = L;
   }   
}

/****************************************************************************/

void MR_3D::init()
{
   Nx=0;Ny=0;Nz=0;
   Nbr_Plan=0;
   Type_Transform = T3_UNDEFINED;
   Set_Transform = S3_UNDEFINED;
   Verbose=False;
   FilterBank=NULL;
   FilterBankAlloc=False;
   LiftingTrans = DEF_LIFT;
   SBFilter = DEF_SB_FILTER;
   TypeNorm = DEF_SB_NORM;
   TabBand=NULL;
}

/****************************************************************************/

void MR_3D::alloc (int Npx, int Npy, int Npz, type_trans_3d T, int Nbr_Scale,
	           FilterAnaSynt *FAS, sb_type_norm Norm)
{
    extern type_3d_format IO_3D_Format;
    int b,s;

    Nx = Npx;
    Ny = Npy;
    Nz = Npz;
    Nbr_Plan = Nbr_Scale;
    Type_Transform = T;
    Set_Transform = SetTransform (T);
    Border = DEFAULT_BORDER_3D;   
    
    DataFormat = IO_3D_Format;
    if (Set_Transform == TRANS3_PAVE)
    {
        Nbr_Band = Nbr_Plan;
	TabSizeNx = new int [Nbr_Band];
        TabSizeNy = new int [Nbr_Band];
        TabSizeNz = new int [Nbr_Band];
        TabPosX = new int [Nbr_Band];
        TabPosY = new int [Nbr_Band];
        TabPosZ = new int [Nbr_Band];
	AT3D_WT.alloc(TabBand, Nx, Ny, Nz, Nbr_Plan);
	for (b =0 ; b < Nbr_Band; b++)
	{
	   TabSizeNx[b] = Nx;
	   TabSizeNy[b] = Ny;
	   TabSizeNz[b] = Nz;
	   TabPosX[b] = 0;
	   TabPosY[b] = 0;
	   TabPosZ[b] = 0;
	}
    }
    else
    {
    FilterBank = FAS;
    if (FAS != NULL) SBFilter = FAS->type_filter();
    else if (T == TO3_MALLAT)   
    {
        if (FilterBank == NULL)
        {
             SBFilter = DEF_SB_FILTER;
             TypeNorm = DEF_SB_NORM;
             FilterBank = new FilterAnaSynt;
             FilterBank->alloc(SBFilter);
             FilterBankAlloc=True;
        }
    }
    TypeNorm = Norm;
   
    // LiftingTrans = DEF_LIFT;
    // SB_Filter = F_MALLAT_7_9;
    // Norm = NORM_L1;
    Nbr_Band = 7 * (Nbr_Plan-1) + 1;
    TabSizeNx = new int [Nbr_Band];
    TabSizeNy = new int [Nbr_Band];
    TabSizeNz = new int [Nbr_Band];
    TabPosX = new int [Nbr_Band];
    TabPosY = new int [Nbr_Band];
    TabPosZ = new int [Nbr_Band];
         
    Data.alloc (Nx,Ny,Nz);        
    TabBand = &Data;
    if (Set_Transform == TRANS3_MALLAT)  
    {
      int Lx = Nx;
      int Ly = Ny;
      int Lz = Nz;
      
      for (s=0; s < Nbr_Plan-1; s++)
      {
         int Ns,P;
	  
 	 // Horizontal low
	 b = 7*s;
	 set_size_band(Lx, False, P, Ns);
         TabPosX[b] = P;    
         TabSizeNx[b] = Ns;
	 set_size_band(Ly, True, P, Ns);
         TabPosY[b] = P;    
         TabSizeNy[b] = Ns;
	 set_size_band(Lz, True, P, Ns);
         TabPosZ[b] = P;    
         TabSizeNz[b] = Ns;
	 
         // Vertical low
	 b = 7*s+1;
	 set_size_band(Lx, True, P, Ns);
         TabPosX[b] = P;    
         TabSizeNx[b] = Ns;
	 set_size_band(Ly, False, P, Ns);
         TabPosY[b] = P;    
         TabSizeNy[b] = Ns;
	 set_size_band(Lz, True, P, Ns);
         TabPosZ[b] = P;    
         TabSizeNz[b] = Ns; 
	 
         // Diagonal low
	 b = 7*s+2;
	 set_size_band(Lx, False, P, Ns);
         TabPosX[b] = P;    
         TabSizeNx[b] = Ns;
	 set_size_band(Ly, False, P, Ns);
         TabPosY[b] = P;    
         TabSizeNy[b] = Ns;
	 set_size_band(Lz, True, P, Ns);
         TabPosZ[b] = P;    
         TabSizeNz[b] = Ns; 	 

	 // Horizontal High
	 b = 7*s+3;
	 set_size_band(Lx, False, P, Ns);
         TabPosX[b] = P;    
         TabSizeNx[b] = Ns;
	 set_size_band(Ly, True, P, Ns);
         TabPosY[b] = P;    
         TabSizeNy[b] = Ns;
	 set_size_band(Lz, False, P, Ns);
         TabPosZ[b] = P;    
         TabSizeNz[b] = Ns;
	 
         // Vertical high
	 b = 7*s+4;
	 set_size_band(Lx, True, P, Ns);
         TabPosX[b] = P;    
         TabSizeNx[b] = Ns;
	 set_size_band(Ly, False, P, Ns);
         TabPosY[b] = P;    
         TabSizeNy[b] = Ns;
	 set_size_band(Lz, False, P, Ns);
         TabPosZ[b] = P;    
         TabSizeNz[b] = Ns; 
	 
         // Diagonal high
	 b = 7*s+5;
	 set_size_band(Lx, False, P, Ns);
         TabPosX[b] = P;    
         TabSizeNx[b] = Ns;
	 set_size_band(Ly, False, P, Ns);
         TabPosY[b] = P;    
         TabSizeNy[b] = Ns;
	 set_size_band(Lz, False, P, Ns);
         TabPosZ[b] = P;    
         TabSizeNz[b] = Ns; 	 

         // z-high, and x,y low
	 b = 7*s+6;
	 set_size_band(Lx, True, P, Ns);
         TabPosX[b] = P;    
         TabSizeNx[b] = Ns;
	 set_size_band(Ly, True, P, Ns);
         TabPosY[b] = P;    
         TabSizeNy[b] = Ns;
	 set_size_band(Lz, False, P, Ns);
         TabPosZ[b] = P;    
         TabSizeNz[b] = Ns;

         Lx = (Lx+1) / 2;
	 Ly = (Ly+1) / 2;
	 Lz = (Lz+1) / 2;
      } 
      b = Nbr_Band-1;
      TabPosX[b] = 0;
      TabPosY[b] = 0;
      TabPosZ[b] = 0;
      TabSizeNx[b] = Lx;
      TabSizeNy[b] = Ly;
      TabSizeNz[b] = Lz;
      
   }}
}
 
/****************************************************************************/

void MR_3D::info_pos_band()
{
   int s;
   
   switch (Set_Transform)
   {
     case TRANS3_MALLAT:
       cout << "Number of bands = " <<  Nbr_Band << endl;
       for (s=0; s < Nbr_Band; s++)
       {
       cout << "Band " << s+1;
       cout <<  " Pos = [" << TabPosX[s] << "," << TabPosY[s] << "," <<  TabPosZ[s] << "]";
       cout <<  " Size = [" << TabSizeNx[s] << "," << TabSizeNy[s] << "," << TabSizeNz[s] << "]";
       cout <<  " Band = [" << TabPosX[s] << ":" << TabPosX[s]+TabSizeNx[s]-1 << "," << 
	                          TabPosY[s] << ":" << TabPosY[s]+TabSizeNy[s]-1 << "," <<  
				  TabPosZ[s] << ":" << TabPosZ[s]+TabSizeNz[s]-1 << "]" << endl;
       }
       break;
     case  TRANS3_PAVE:
       cout << "Number of bands = " <<  Nbr_Band << endl;
       for (s=0; s < Nbr_Band; s++)
       {
          cout << "Band " << s+1;
          cout <<  " Size = [" << TabSizeNx[s] << "," << TabSizeNy[s] << "," << TabSizeNz[s] << "]" << endl;
       }
      break;
     default: cerr << "Error: bad transform ... " << endl;
              exit(-1);
   }
}

/****************************************************************************/

void MR_3D::get_band(int b, fltarray &Band)
{
   int i,j,k;
   int Nxb = size_band_nx(b);
   int Nyb = size_band_ny(b);
   int Nzb = size_band_nz(b);
   
   if (Band.n_elem() == 0) Band.alloc(Nxb,Nyb,Nzb);
   else if ((Band.naxis() != 3) || (Band.nx() != Nxb) 
             || (Band.ny() != Nyb) || (Band.nz() != Nzb))
   {
      Band.free();
      Band.alloc(Nxb,Nyb,Nzb);
   }
   for (i=0; i < Nxb; i++)
   for (j=0; j < Nyb; j++)
   for (k=0; k < Nzb; k++) Band(i,j,k) = (*this)(b,i,j,k);
}

/****************************************************************************/

void MR_3D::insert_band(int b, fltarray &Band)
{
   int i,j,k;
   int Nxb = size_band_nx(b);
   int Nyb = size_band_ny(b);
   int Nzb = size_band_nz(b);
   
   if ((Band.n_elem() == 0) || (Band.naxis() != 3) 
       || (Band.nx() != Nxb) 
       || (Band.ny() != Nyb) || (Band.nz() != Nzb))
   {
      cerr << "Error: band to insert has not the correct dimensions ... " << endl;
      exit(-1);
   }
   for (i=0; i < Nxb; i++)
   for (j=0; j < Nyb; j++)
   for (k=0; k < Nzb; k++) (*this)(b,i,j,k) = Band(i,j,k);
}

/****************************************************************************/

void MR_3D::info_band(int b)
{
   int i,j,k;
   double Min = (*this)(b,0,0,0);
   double  Max = Min;
   double Flux = 0;
   double  Mean, Sigma=0;
   int Nxb = size_band_nx(b);
   int Nyb = size_band_ny(b);
   int Nzb = size_band_nz(b);
   long Np = Nxb*Nyb*Nzb;
   
   for (i=0; i < Nxb; i++)
   for (j=0; j < Nyb; j++)
   for (k=0; k < Nzb; k++)
   {
      double Coef = (*this)(b,i,j,k);
      if (Coef < Min) Min = (*this)(b,i,j,k);
      if (Coef > Max) Max = (*this)(b,i,j,k);
      Flux += Coef;
   }
   Mean = Flux / (double) Np;
   for (i=0; i < Nxb; i++)
   for (j=0; j < Nyb; j++)
   for (k=0; k < Nzb; k++)
   {
      double Coef = (*this)(b,i,j,k) - Mean;
      Sigma += Coef*Coef;
   }
   Sigma /= (double) Np;

   cout << "Band " << b+1 << endl;
   cout <<  " Pos = [" << TabPosX[b] << "," << TabPosY[b] << "," <<  TabPosZ[b] << "]";
   cout <<  " Size = [" << Nxb << "," << Nyb << "," << Nzb << "]";
   cout <<  " Band = [" << TabPosX[b] << ":" << TabPosX[b]+Nxb-1 << "," << 
	                    TabPosY[b] << ":" << TabPosY[b]+Nyb-1 << "," <<  
		            TabPosZ[b] << ":" << TabPosZ[b]+Nzb-1 << "]" << endl;   
   cout << "  Min   = " << Min <<  " Max = " << Max << endl;
   cout << "  Mean  = " << Mean << " Flux = " << Flux << endl;
   cout << "  Sigma = " <<  Sigma << endl << endl; 
}


/****************************************************************************/

MR_3D::~MR_3D()
{
   free ();
}

/****************************************************************************/

void MR_3D::free ()
{

    Border = DEFAULT_BORDER_3D;
    switch (Set_Transform)
    {
     case TRANS3_MALLAT:
         Data.free();
        break;
     case  TRANS3_PAVE:
         AT3D_WT.free(TabBand, Nbr_Plan);
       break;
     default: cerr << "Error: bad transform ... " << endl;
              exit(-1);
    }
    Nbr_Plan = 0;
    Nbr_Band = 0;
    Nx = Ny = Nz = 0;
    Set_Transform= S3_UNDEFINED;
    Type_Transform=T3_UNDEFINED;
    if (TabPosX != NULL)  delete [] TabPosX;
    if (TabSizeNx != NULL) delete [] TabSizeNx;
    if (TabPosY != NULL)  delete [] TabPosY;
    if (TabSizeNy != NULL) delete [] TabSizeNy;
    if (TabPosZ != NULL)  delete [] TabPosZ;
    if (TabSizeNz != NULL) delete [] TabSizeNz;    
    if ((FilterBankAlloc == True) && (FilterBank != NULL)) 
    {
        delete FilterBank;
        FilterBank = NULL;
    }
    init();
}

/****************************************************************************/

void MR_3D::transform (fltarray &Cube, type_border Bord)
{
   if (Nbr_Plan < 2)
   {
      cerr << "Error: Object not correctly allocated: Nbr_Plan must be > 2  ";
      cerr << endl;
      exit (-1);
   }
    switch (Type_Transform)
    {
       case TO3_ATROUS:
          AT3D_WT.Bord = Bord;
 	  AT3D_WT.transform(Cube, TabBand, Nbr_Plan);
          break;
       case TO3_MALLAT:
        {
	    Data = Cube;
            SubBandFilter WT1D(*FilterBank, TypeNorm);
            Ortho_3D_WT WT3D(WT1D);
            WT3D.transform(Data, Nbr_Plan);
         }
         break;
       case TO3_LIFTING:
        {
	  Data = Cube;
          Lifting Clift1D(LiftingTrans);
          Ortho_3D_WT WT3D(Clift1D);
          WT3D.transform(Data, Nbr_Plan);
         }
         break; 
       default:
           fprintf (stderr, "Error (proc. MR_3D_transform): Unknown transform\n");
           exit (-1);
           break;
    }    
}

/****************************************************************************/

void MR_3D::transform (fltarray &Cube)
{
   this->transform(Cube, Border);
}

/****************************************************************************/

void MR_3D::recons (fltarray &Cube)
{
   this->recons(Cube, Border);
}

/****************************************************************************/

void MR_3D::recons (fltarray &Cube, type_border Bord)
{
    Cube = Data;
    switch (Type_Transform)
    {
       case TO3_ATROUS:
          AT3D_WT.Bord = Bord;
 	  AT3D_WT.recons(TabBand, Cube, Nbr_Plan);
          break;
	case TO3_MALLAT:
        {
            SubBandFilter WT1D(*FilterBank, TypeNorm);
            Ortho_3D_WT WT3D(WT1D);
            WT3D.recons(Cube, Nbr_Plan);
         }
         break;
       case TO3_LIFTING:
        {
          Lifting Clift1D(LiftingTrans);
          Ortho_3D_WT WT3D(Clift1D);
          WT3D.recons(Cube, Nbr_Plan);
         }
         break; 
       default:
           fprintf (stderr, "Error (proc. MR_3D_transform): Unknown transform\n");
           exit (-1);
           break;
    }     
}

/****************************************************************************/

static void mr3d_io_name (char *File_Name_In, char *File_Name_Out)
{
    int L;

    strcpy (File_Name_Out, File_Name_In);

    L = strlen (File_Name_In);
    if ((L < 3) || (File_Name_In[L-1] != 'r')
                || (File_Name_In[L-2] != 'm')
                || (File_Name_In[L-3] != '.'))
    {
        strcat (File_Name_Out, ".mr");
    }
}

/****************************************************************************/
static void PrintError( int status)
{
    /*****************************************************/
    /* Print out cfitsio error messages and exit program */
    /*****************************************************/

    char status_str[FLEN_STATUS], errmsg[FLEN_ERRMSG];
  
    if (status)
      fprintf(stderr, "\n*** Error occurred during program execution ***\n");

    ffgerr(status, status_str);        /* get the error status description */
    fprintf(stderr, "\nstatus = %d: %s\n", status, status_str);

    if ( ffgmsg(errmsg) )  /* get first message; null if stack is empty */
    {
         fprintf(stderr, "\nError message stack:\n");
         fprintf(stderr, " %s\n", errmsg);

         while ( ffgmsg(errmsg) )  /* get remaining messages */
             fprintf(stderr, " %s\n", errmsg);
    }

    exit( status );       /* terminate the program, returning error status */
}

/******************************************************************/

int MR_3D::mr_io_fill_header(fitsfile *fptr)
{
  int status = 0; // this means OK for cfitsio !!!
    /*****************************************************/
     /* write optional keyword to the header */
    /*****************************************************/


// if ( ffpkyj(fptr, "Nx", (long)Nx,"x-axis size",&status)) PrintError( status );  
// if ( ffpkyj(fptr,"Ny",(long)Ny,"y-axis size",&status)) PrintError( status );  
// if ( ffpkyj(fptr,"Nz",(long)Nz,"z-axis size",&status)) PrintError( status ); 
 if ( ffpkyj(fptr, (char*)"Nbr_Plan", (long)Nbr_Plan, (char*)"Number of scales", &status))
     PrintError( status );  
 if (ffpkyj(fptr, (char*)"Set_Transform", (long)Set_Transform,
                      StringSet3D(Set_Transform), &status))
      PrintError( status );  
 if ( ffpkyj(fptr, (char*)"Type_Transform", (long)Type_Transform, 
                         StringTransf3D(Type_Transform), &status))
     PrintError( status );  

 if (Type_Transform  == TO3_MALLAT)
 {
  if ( ffpkyj(fptr, (char*)"SBFilter", (long) SBFilter, (char*)"Type of filters", &status))
        PrintError( status );  
     if ( ffpkyj(fptr, (char*)"NORM", (long)  TypeNorm, (char*)"normalization", &status))
        PrintError( status );  
  // if ( ffpkyj(fptr, "UNDEC", (long) NbrUndec , "nbr of undec. scales", &status))
  //      PrintError( status ); 
  }	
       
 
  if (Type_Transform  == TO3_LIFTING)
  if ( ffpkyj(fptr, (char*)"LiftingTrans", (long)  LiftingTrans, 
                           (char*)StringLSTransform(LiftingTrans), &status))
        PrintError( status );
     
 if ( ffpkyj(fptr, (char*)"DataFormat",(long) DataFormat,(char*)"Input data format", &status))
     PrintError( status );  

 if ( ffpkyj(fptr, (char*)"Border",(long)Border,(char*)"border type", &status))
     PrintError( status );  

 if (SBFilter == F_USER)
 {
   if (UserFilterFileName != NULL)
   {
     if ( ffpkys(fptr, (char*)"FilBank", UserFilterFileName ,(char*)"Filter", &status))
           PrintError( status );
   }
   else if ( ffpkys(fptr, (char*)"FilBank", (char*)"-" ,(char*)"Filter", &status))
          PrintError( status );
 }
     
 return(status);
} 


/****************************************************************************/

void MR_3D::write (char *Name)
/* new version with fits */
{
 char filename[256];
 Ifloat Ima;
 fitsfile *fptr;    
 int status;
 // float *Ptr;
 int b,simple;
 int bitpix;
 long naxis=0;
 long naxes[4];
 long nelements;
 long group = 1;  /* group to write in the fits file, 1= first group */
 long firstpixel = 1;    /* first pixel to write (begin with 1) */
 Ifloat Aux;
 
/* we keep mr as extension even if its fits ! */
 mr3d_io_name (Name, filename);

#if DEGUG_IO  
    cout << "Write on " << filename << endl;
#endif

 FILE *FEXIST = fopen(filename, "rb");
 if (FEXIST)
 {
    fclose(FEXIST);
    remove(filename);               /* Delete old file if it already exists */
 }

 status = 0;         /* initialize status before calling fitsio routines */

    /* open the file */
 if ( ffinit(&fptr, filename, &status) )     /* create the new FITS file */
     PrintError( status );           /* call PrintError if error occurs */
                                                                              
/* write  the header */
 simple   = True;
 bitpix   =  -32;   /* 32-bit real pixel values      */
 long pcount   =   0;  /* no group parameters */
 long gcount   =   1;  /* only a single image/group */
 int  extend   =   False;
 switch (Set_Transform)
     {
      case TRANS3_MALLAT:
         naxis = 3;
         naxes[0] = Nx;
         naxes[1] = Ny;
	 naxes[2] = Nz;
         break;
     case TRANS3_PAVE:
         naxis = 4;
         naxes[0] = Nx;
         naxes[1] = Ny;
	 naxes[2] = Nz;
	 naxes[3] = Nbr_Plan;
         break;
     default:
         fprintf (stderr, "Error in mr_io_write: bad Set_Transform ... \n");
         exit(-1);
         break; 
     }
 
 // write first header part (parameters)
 if ( ffphpr(fptr,simple,bitpix,naxis,naxes,pcount,gcount,extend,&status) )
     PrintError( status );          /* call PrintError if error occurs */
    
 // write the header of the multiresolution file
 status = mr_io_fill_header(fptr);

 // save the data
 // long fpixels[3];
 // long lpixels[3];
 switch (Set_Transform)
     {
     case TRANS3_MALLAT:
             // this works because everything is put in one image
         nelements = naxes[0] * naxes[1] * naxes[2];
         if ( ffppre(fptr, group, firstpixel, nelements, Data.buffer(), 
               &status) )
              PrintError( status );  
         break;
     case TRANS3_PAVE:
         nelements = naxes[0] * naxes[1] * naxes[2];
         for (b=0; b < Nbr_Plan; b++)
	 {
	   if ( ffppre(fptr, group, firstpixel, nelements, (TabBand[b]).buffer(), 
               &status) )
              PrintError( status );
	    firstpixel += nelements;
	 }
	 break;
     default:
         fprintf (stderr, "Error in mr_io_write: bad Type_Transform ..\n\n");
         break; 
     }

 /* close the FITS file */
 if ( ffclos(fptr, &status) )  PrintError(status);  
// cout << " end of write fits " << endl;

}

/****************************************************************************/

void MR_3D::read (char *Name)
{
    // for fits
    char filename[256];
    fitsfile *fptr;           /* pointer to the FITS file */
    int status=0, hdutype ;
    long hdunum;
    char comment[FLEN_COMMENT];
    int  naxis;
    long naxes[4];
    long mon_long;
    float nulval = 0.;
    int anynul = 0;
    long inc[3];
    void PrintError( int status);
    long nelements = 0 ; // naxes[0] * naxes[1] in the image
    Ifloat Ima;
    long firstpixel = 1;
    
     // for multiresol
    float *Ptr;
 
     mr3d_io_name (Name, filename);

    inc[0]=1;  inc[1]=1; inc[2]=1;
 
#if DEBUG_IO
    cout << "Read in " << filename << endl;
#endif
   
    /* open the file */
    status = 0;         /* initialize status before calling fitsio routines */
    if ( ffopen(&fptr, filename, (int) READONLY, &status) ) 
         PrintError( status );
 // cout << " open the file " << endl;
                                   
    // ******* read the HEADER  *******

    hdunum = 1;  /*read  table */
    if ( ffmahd(fptr, hdunum, &hdutype, &status) ) /* move to the HDU */
           PrintError( status );

    int simple, bitpix, extend;
    long pcount, gcount;
    if ( ffghpr(fptr, 3, &simple, &bitpix, &naxis, naxes, &pcount,
            &gcount, &extend, &status)) /* move to the HDU */
           PrintError( status );

     nelements = naxes[0] * naxes[1] * naxes[2];
// cout << " begin to read " << endl;
     /* read Number of lines, columns, plans */
    //if (ffgkyj(fptr,"Nx", &mon_long, comment, &status)) PrintError( status );
    //int Nxi = (int)mon_long; /* there is no function for reading int */
    //if (ffgkyj(fptr,"Ny", &mon_long, comment, &status)) PrintError( status );
    //int Nyi = (int)mon_long; 
    //if (ffgkyj(fptr,"Nz", &mon_long, comment, &status)) PrintError( status );
    //int Nzi = (int)mon_long;
    int Nxi = (int)naxes[0];
    int Nyi = (int)naxes[1];
    int Nzi = (int)naxes[2];
    
    if (ffgkyj(fptr,(char*)"Nbr_Plan", &mon_long, comment, &status)) PrintError( status );
    int NbrPlan = (int)mon_long; 
    
    type_trans_3d  TT;
    if (ffgkyj(fptr,(char*)"Type_Transform", &mon_long, comment, &status)) 
                    PrintError( status );
    else TT = (type_trans_3d) mon_long;
		    
    FilterAnaSynt *PtrFAS = NULL;
    if (TT == TO3_MALLAT)
    {
        if (ffgkyj(fptr,(char*)"SBFilter", &mon_long, comment, &status)) SBFilter = DEF_SB_FILTER;
        else SBFilter = (type_sb_filter) mon_long;
	
	if (ffgkyj(fptr,(char*)"NORM", &mon_long, comment, &status)) TypeNorm = DEF_SB_NORM;
        else  TypeNorm = (sb_type_norm) mon_long;

      // if (ffgkyj(fptr,"UNDEC", &mon_long, comment, &status))  NbrUndecimatedScale = -1;
      // else NU = (int) mon_long;
 
      if (SBFilter == F_USER)
      {
          char FBName[256];
          if (ffgkys(fptr,(char*)"FilBank", FBName, comment, &status))
                                                      PrintError( status );
          if (FBName[0] == '-') UserFilterFileName = NULL; 
          else UserFilterFileName = strdup(FBName);
      }
    
       PtrFAS = new FilterAnaSynt;
       PtrFAS->Verbose = Verbose;
       PtrFAS->alloc(SBFilter);
       FilterBankAlloc=True;
    }
    alloc(Nxi,Nyi,Nzi, TT, NbrPlan, PtrFAS, TypeNorm);   
     	
    if (Type_Transform  == TO3_LIFTING)
    {
         if (ffgkyj(fptr,(char*)"LiftingTrans", &mon_long, comment, &status)) 
              LiftingTrans = DEF_LIFT;
 	 else LiftingTrans = (type_lift) mon_long;
    }

    if (ffgkyj(fptr,(char*)"DataFormat", &mon_long, comment, &status))
         PrintError( status );
    DataFormat = (type_3d_format) mon_long;

    if (ffgkyj(fptr,(char*)"Border", &mon_long, comment, &status))
         PrintError( status );
    Border = (type_border)mon_long; 

    
#if DEBUG_IO 
    cout << "Read in " << filename << endl;
    cout << "Nx = " << Nx << endl;
    cout << "Ny = " << Ny << endl;
    cout << "Nz = " << Nz << endl;
    cout << "Nbr_Plan = " << nbr_scale () << endl;
    cout << "Type_Transform = " << Type_Transform << " " <<
            StringTransf3D(Type_Transform) << endl;
    cout << "Set_Transform = " << Set_Transform << " " <<
            StringSet3D(Set_Transform) <<  endl;
#endif

    // ******* read the images *******

    long fpixels[3];
    long lpixels[3];
    fpixels[0] = 1;
    fpixels[1] = 1;
    fpixels[2] = 1;
    lpixels[0] = Nx;
    lpixels[1] = Ny;    
    lpixels[2] = Nz;

    switch (Set_Transform)
    {
        case TRANS3_MALLAT:
              Ptr = Data.buffer();
             if ( ffgpve(fptr, 1, 1, nelements, nulval,Ptr, &anynul, &status))
             PrintError( status );
             break;
        case  TRANS3_PAVE:
         nelements = naxes[0] * naxes[1] * naxes[2];
         for (int b=0; b < Nbr_Plan; b++)
	 {
	   Ptr = (TabBand[b]).buffer();
	   if (  ffgpve(fptr, 1, firstpixel, nelements, nulval,Ptr, &anynul, &status))
              PrintError( status );
	    firstpixel += nelements;
	 }
	 break;        
   default:
          fprintf (stderr, "Error in mr_io_read: bad Set_Transform .. \n");
          break; 
    }
 
/* close the FITS file */
 if ( ffclos(fptr, &status) ) PrintError( status );
// cout << " end of read fits file " << endl;
}

/****************************************************************************/


