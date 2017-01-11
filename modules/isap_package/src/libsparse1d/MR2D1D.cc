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
**    File:  MR2D1D.cc
**
*******************************************************************************
**
**    DESCRIPTION  multiresolution transform of cube, using 2D wavelet transform
**    -----------  followed (or not) by a 1D wavelet transform
**                 
**
******************************************************************************/


#include "MR2D1D.h"


/****************************************************************************/

static void mr_io_name (char *File_Name_In, char *File_Name_Out)
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

/*--------------------------------------------------------------------------*/
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

/*--------------------------------------------------------------------------*/


int MR2D1D::mr_io_fill_header( fitsfile *fptr)
{
  int status = 0; // this means OK for cfitsio !!!
    /*****************************************************/
     /* write optional keyword to the header */
    /*****************************************************/

// cout << " DDDD " << (long) NbrBand2D << " " << (long)NbrBand1D << endl;

if ( ffpkyj(fptr, "Nx", (long)Nx,"Nx of the input cube",&status))
     PrintError( status );  
if ( ffpkyj(fptr,"Ny",(long)Ny,"Ny of the input cube",&status))
     PrintError( status );  
if ( ffpkyj(fptr,"Nz",(long)Nz,"Nz of the input cube",&status))
     PrintError( status );  
 if ( ffpkyj(fptr, "NSCALE2D", (long) NbrScale2D, "Number of bands 2D", &status))
     PrintError( status );  
 if ( ffpkyj(fptr, "NSCALE1D", (long)NbrScale1D, "Number of bands 1D", &status))
     PrintError( status );  
 if ( ffpkyj(fptr, "Type_Tra", (long) WT2D.Type_Transform, 
                         (char*)StringTransform(WT2D.Type_Transform), &status))
     PrintError( status );  
 return(status);
} 


/****************************************************************************/

void  MR2D1D::read(char *Name)
{
    char filename[256];
    fitsfile *fptr;           /* pointer to the FITS file */
    int status=0, hdutype ;
    long hdunum;
    char comment[FLEN_COMMENT]; 
    int naxis;
    long naxes[3];
    long mon_long;
    int anynul = 0;
    long nulval = 0;
    long inc[3];
    void PrintError( int status);
    long nelements = 0 ; // naxes[0] * naxes[1] in the image
    //long fpixels[3];
    //long int lpixels[3];
    
     // for multiresol
    float *Ptr;
    //int my_logical; // sais pas...

     mr_io_name (Name, filename);

    inc[0]=1;  inc[1]=1; inc[2]=1;
 
#if DEBUG_IO
    cout << "Read in " << filename << endl;
#endif
   
    /* open the file */
    status = 0;         /* initialize status before calling fitsio routines */
    if ( ffopen(&fptr, filename, (int) READONLY, &status) ) 
         PrintError( status );
                                    
    hdunum = 1;  /*read  table */
    if ( ffmahd(fptr, hdunum, &hdutype, &status) ) /* move to the HDU */
           PrintError( status );

    int simple, bitpix, extend;
    long pcount, gcount;
    if ( ffghpr(fptr, 3, &simple, &bitpix, &naxis, naxes, &pcount,
            &gcount, &extend, &status)) /* move to the HDU */
           PrintError( status );

    nelements = naxes[0];
     
    if (ffgkyj(fptr,"Nx", &mon_long, comment, &status)) PrintError( status );
    int Nxi = (int) mon_long;  
    if (ffgkyj(fptr,"Ny", &mon_long, comment, &status)) PrintError( status );
    int Nyi = (int) mon_long; 
    if (ffgkyj(fptr,"Nz", &mon_long, comment, &status)) PrintError( status );
    int Nzi = (int) mon_long; 
    
    if (ffgkyj(fptr,"NSCALE2D", &mon_long, comment, &status)) PrintError( status );
    int iNbrScale2D  = (int) mon_long;
    if (ffgkyj(fptr,"NSCALE1D", &mon_long, comment, &status)) PrintError( status );
    int iNbrScale1D  = (int) mon_long;

    type_transform  TT;
    if (ffgkyj(fptr,"Type_Tra", &mon_long, comment, &status))   PrintError( status );
    else TT = (type_transform) mon_long;
    
    alloc(Nxi,Nyi,Nzi,TT, iNbrScale2D, iNbrScale1D);
     
    fltarray Tab(nelements);
   //cout << "        Nl = " << (TabCF_Band[0][1]).nl() << " Nc = " << TabCF_Band[0][1].nc() << endl;

    Ptr = Tab.buffer();
    if ( ffgpve(fptr, 1, 1, nelements, nulval, Ptr, &anynul, &status))
             PrintError( status );

    int ind=2;
    for (int s=0; s <NbrBand2D; s++) 
    for (int s1=0; s1 < NbrBand1D; s1++)
    {
    
       int Nxb = (int) Tab(ind++);
       int Nyb = (int) Tab(ind++);
       int Nzb = (int) Tab(ind++);
       //  cout << s << " " << s1 << " " << Nxb << " " << Nyb << " " << Nzb << endl;

        //cout << " s = " << s+1 << " b = " << b+1 << " Nl = " << Ny << " Nc = " << Nx << endl;
       //cout << "        Nl = " << (TabCF_Band[s][b]).nl() << " Nc = " << TabCF_Band[s][b].nc() << endl;
       for (int i=0;  i < Nxb; i++)
       for (int j=0;  j < Nyb; j++)
       for (int k=0;  k < Nzb; k++)  (*this)(s,s1,i,j,k) = Tab(ind++);
   }
      
   if ( ffclos(fptr, &status) ) PrintError( status );
// cout << " end of read fits file " << endl;
}

/****************************************************************************/

void  MR2D1D::write(char *Name)
{

 char filename[256];
 Ifloat Ima;
 fitsfile *fptr;    
 int status;
 int simple;
 int bitpix;
 long naxis=0;
 long naxes[4];
 long nelements;
 long group = 1;  /* group to write in the fits file, 1= first group */
 long firstpixel = 1;    /* first pixel to write (begin with 1) */
 Ifloat Aux;
 
/* we keep mr as extension even if its fits ! */
 mr_io_name (Name, filename);

#if DEBUG_IO  
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
  naxis = 1;

  
  int Nelem=2;
  for (int s=0; s < NbrBand2D; s++)
  for (int s1=0; s1 < NbrBand1D; s1++) 
  {
     Nelem += 3 +  TabSizeBandNx(s,s1)*TabSizeBandNy(s,s1)*TabSizeBandNz(s,s1);
  } 
  // cout << " NELEM = " << Nelem << endl;
  fltarray Data(Nelem);
  int ind=2;
  Data(0) = NbrBand2D;
  Data(1) = NbrBand1D;
 
  // Ifloat Band;
  for (int s=0; s < NbrBand2D; s++)
  for (int s1=0; s1 < NbrBand1D; s1++) 
  {
     Data(ind++) = TabSizeBandNx(s,s1);
     Data(ind++) = TabSizeBandNy(s,s1);
     Data(ind++) = TabSizeBandNz(s,s1);
    
     for (int i=0;  i < TabSizeBandNx(s,s1); i++)
     for (int j=0;  j < TabSizeBandNy(s,s1); j++) 
     for (int k=0;  k < TabSizeBandNz(s,s1); k++) Data(ind++) = (*this)(s,s1,i,j,k);
  }
  // cout << " DATAOK = " <<   endl;
  naxes[0] = ind;
  
 // write first header part (parameters)
 if ( ffphpr(fptr,simple,bitpix,naxis,naxes,pcount,gcount,extend,&status) )
     PrintError( status );          /* call PrintError if error occurs */
  
   // write the header of the multiresolution file
  mr_io_fill_header(fptr);
 
 nelements = ind;
 if ( ffppre(fptr, group, firstpixel, nelements, Data.buffer(), &status) )
              PrintError( status );  
     
  /* close the FITS file */
  if ( ffclos(fptr, &status) )  PrintError( status ); 
}

/****************************************************************************/

inline float & MR2D1D::operator() (int s2, int s1, int i, int j, int k) const
{    
#ifdef TESTIND
   if ( (i < 0) || (i >= size_band_nx(s2,s1)) ||
        (j < 0) || (j >= size_band_ny(s2,s1)) ||
	(k < 0) || (k >= size_band_nz(s2,s1)) ||
	(s2 < 0) || (s2 >= nbr_band_2d()) ||
	(s1 < 0) || (s1 >= nbr_band_1d()))
   {
      printf("Error: (s2,s1,i,j,k)=(%d,%d,%d,%d,%d), size band = (%d,%d,%d) , scale=(%d,%d) \n", s2,s1,i,j,k,size_band_nx(s2,s1), size_band_ny(s2,s1), size_band_nz(s2,s1), nbr_band_2d(), nbr_band_1d());
      exit(-1);
   }
#endif
  
   return TabBand[s2] (i,j,k+TabFirstPosBandNz(s1));
}


/****************************************************************************/

inline float & MR2D1D::operator() (int s2, int i, int j, int k) const
{    
#ifdef TESTIND
   if ( (i < 0) || (i >= size_band_nx(s2)) ||
        (j < 0) || (j >= size_band_ny(s2)) ||
	(k < 0) || (k >= size_band_nz(s2)) ||
	(s2 < 0) || (s2 >= nbr_band_2d())  )
   {
      printf("Error: (s2,i,j,k)=(%d,%d,%d,%d), size band = (%d,%d,%d) , scale=(%d,%d) \n", s2,i,j,k,size_band_nx(s2), size_band_ny(s2), size_band_nz(s2), nbr_band_2d(), nbr_band_1d());
      exit(-1);
   }
#endif
  
   return TabBand[s2] (i,j,k);
}

/****************************************************************************/

fltarray MR2D1D::get_band(int s2, int s1)
{
    int Nxb = size_band_nx(s2,s1);
    int Nyb = size_band_ny(s2,s1);
    int Nzb = size_band_nz(s2,s1);
    fltarray *Cube_Return = NULL;

    Cube_Return = new fltarray(Nxb, Nyb, Nzb);
    for (int i=0; i < Nxb; i++)
    for (int j=0; j < Nyb; j++)
    for (int k=0; k < Nzb; k++) (*Cube_Return)(i,j,k) = (*this)(s2,s1,i,j,k);
    
    return (*Cube_Return);
}

/****************************************************************************/

void  MR2D1D::put_band(fltarray Band, int s2, int s1)
{
    int Nxb = size_band_nx(s2,s1);
    int Nyb = size_band_ny(s2,s1);
    int Nzb = size_band_nz(s2,s1);
    
    for (int i=0; i < Nxb; i++)
    for (int j=0; i < Nyb; j++)
    for (int k=0; k < Nzb; k++) (*this)(s2,s1,i,j,k) = Band(i,j,k);
}

/****************************************************************************/

void MR2D1D::alloc (int iNx, int iNy, int iNz, type_transform Trans2D, int Ns2D, int Ns1D, Bool NoAlloc)
{
   Nx = iNx;
   Ny = iNy;
   Nz = iNz;
   NbrScale2D = Ns2D;
   NbrScale1D = Ns1D;
   if (NbrScale1D < 2)
   {
      NbrScale1D = 1;
      Apply1DTrans = False;
   }
   else Apply1DTrans = True;
   ;
   Norm = NORM_L2;
   SB_Filter = F_MALLAT_7_9;
   Bord = I_CONT;
   U_Filter = DEF_UNDER_FILTER; 
   FilterAnaSynt *PtrFAS = NULL;
    if ((Trans2D == TO_MALLAT) || (Trans2D == TO_UNDECIMATED_MALLAT))
    {
        FAS.Verbose = Verbose;
        FAS.alloc(SB_Filter);
        PtrFAS = &FAS;
    }
    int NbrUndec = -1;                     /*number of undecimated scale */
    type_lift LiftingTrans = DEF_LIFT;
    if (Trans2D == TO_LIFTING) WT2D.LiftingTrans = LiftingTrans;
    WT2D.Border = Bord;
    WT2D.Verbose = Verbose;    
    WT2D.alloc (Ny, Nx, Ns2D, Trans2D, PtrFAS, Norm, NbrUndec, U_Filter);
    NbrBand2D = WT2D.nbr_band();

     	       
    Bool Rebin=False;
    WT1D.U_Filter = U_Filter;
    type_trans_1d Trans1D = TO1_MALLAT;
    if (Apply1DTrans == True)
    {
        WT1D.alloc (Nz, Trans1D, Ns1D, PtrFAS, Norm, Rebin);   
        NbrBand1D = WT1D.nbr_band();
    } 
    else NbrBand1D = 1;
    
   if (NoAlloc == False) TabBand = new fltarray [NbrBand2D];
   else TabBand=NULL;
   TabSizeBandNx.resize(NbrBand2D, NbrBand1D);
   TabSizeBandNy.resize(NbrBand2D, NbrBand1D);
   TabSizeBandNz.resize(NbrBand2D, NbrBand1D);
   TabFirstPosBandNz.resize(NbrBand1D);
   TabFirstPosBandNz(0) =0;
   NbrCoef2D = 0;
   for (int b=0; b < NbrBand2D; b++) 
   {
      int Npz = (Apply1DTrans == True) ? WT1D.size_ima_np (): Nz;
      if (NoAlloc == False)  TabBand[b].alloc(WT2D.size_band_nc(b), WT2D.size_band_nl(b),  Npz);
      NbrCoef2D +=  WT2D.size_band_nc(b) * WT2D.size_band_nl(b);
      for (int b1=0; b1 < NbrBand1D; b1++)
      {
         TabSizeBandNx(b,b1) = WT2D.size_band_nc(b);
	 TabSizeBandNy(b,b1) = WT2D.size_band_nl(b);
 	 if (Apply1DTrans == True)  TabSizeBandNz(b,b1) = WT1D.size_scale_np(b1);
	 else TabSizeBandNz(b,b1) = Nz;
      }
   }
   NbrCoef1D = Nz;
   if (Apply1DTrans == True)  
     for (int b1=1; b1 < NbrBand1D; b1++) TabFirstPosBandNz(b1) = TabFirstPosBandNz(b1-1) + TabSizeBandNz(0,b1-1);
}

/****************************************************************************/

void MR2D1D::transform_to_vectarray(fltarray &Data, fltarray &TabVect)
{
   int Nsx = nbr_coef_2d();
   int Nsy = nbr_coef_1d();
   cout << "ALLOC " << Nsx << " " << Nsy << endl;
   TabVect.alloc(nbr_coef_2d(), nbr_coef_1d());
   int i,j,b,z;
   Nx = Data.nx();
   Ny = Data.ny();
   Nz = Data.nz();
   Ifloat Frame(Ny, Nx);
   fltarray Vect(Nz);
   intarray TabPosBand(nbr_band_2d());
    
   // 2D wt transform per frame
   for (z=0; z < Nz; z++)
   {
      int Pix=0;
      for (i=0; i < Ny; i++)
      for (j=0; j < Nx; j++) Frame(i,j) = Data(j,i,z);
      WT2D.transform(Frame);
      for (b=0; b < nbr_band_2d(); b++)
      {
         for (i=0; i < WT2D.size_band_nl(b); i++)
	 for (j=0; j < WT2D.size_band_nc(b); j++) TabVect(Pix++,z) = WT2D(b,i,j);
	 if (z == 0) TabPosBand(b) = Pix;
       }
   }
 cout << "1D " << NbrBand1D << endl;
   // 1D wt 
   if (NbrBand1D >= 2)
   {
     for (b=0; b <  nbr_band_2d(); b++)
     for (i=0; i < WT2D.size_band_nl(b); i++)
     for (j=0; j < WT2D.size_band_nc(b); j++) 
     {
        for (z=0; z < Nz; z++) Vect(z) = TabVect(TabPosBand(b)+i*WT2D.size_band_nc(b)+j,z);     // TabBand[b](j,i,z);
        WT1D.transform(Vect);
        z = 0;
        for (int b1=0; b1 < NbrBand1D; b1++)
        {
         for (int p=0; p < WT1D.size_scale_np (b1); p++)  TabVect(TabPosBand(b)+i*WT2D.size_band_nc(b)+j,z) = WT1D(b1,p); 
        }
      }
   }    
 cout << "END TRANS " << endl;
}


/****************************************************************************/

void MR2D1D::transform (fltarray &Data)
{
   int i,j,b,z;
   Nx = Data.nx();
   Ny = Data.ny();
   Nz = Data.nz();
   Ifloat Frame(Ny, Nx);
   fltarray Vect(Nz);
   
   // 2D wt transform per frame
   for (z=0; z < Nz; z++)
   {
      for (i=0; i < Ny; i++)
      for (j=0; j < Nx; j++) Frame(i,j) = Data(j,i,z);
      WT2D.transform(Frame);
      for (b=0; b < NbrBand2D; b++)
      {
         for (i=0; i < WT2D.size_band_nl(b); i++)
	 for (j=0; j < WT2D.size_band_nc(b); j++) TabBand[b](j,i,z) = WT2D(b,i,j);
      }
   }
 
   // 1D wt 
   if (NbrBand1D >= 2)
   {
     for (b=0; b < NbrBand2D; b++)
     for (i=0; i < WT2D.size_band_nl(b); i++)
     for (j=0; j < WT2D.size_band_nc(b); j++) 
     {
        for (z=0; z < Nz; z++) Vect(z) = TabBand[b](j,i,z);
        WT1D.transform(Vect);
        z = 0;
        for (int b1=0; b1 < NbrBand1D; b1++)
        {
         for (int p=0; p < WT1D.size_scale_np (b1); p++) TabBand[b](j,i,z++) = WT1D(b1,p); 
        }
      }
   }
}

/****************************************************************************/



void MR2D1D::recons (fltarray &Data)
{
   if ((Data.nx() != Nx) || (Data.ny() != Ny) || (Data.nz() != Nz)) Data.resize(Nx, Ny, Nz); 
   int i,j,b,z;
   Nx = Data.nx();
   Ny = Data.ny();
   Nz = Data.nz();
   Ifloat Frame(Ny, Nx);
   fltarray Vect(Nz);
     
   // 1D wt 
   if (NbrBand1D >= 2)
   {
      for (b=0; b < NbrBand2D; b++)
      for (i=0; i < WT2D.size_band_nl(b); i++)
      for (j=0; j < WT2D.size_band_nc(b); j++) 
      {
         // for (z=0; z < Nz; z++) Vect(z) = TabBand[b](j,i,z);
	z = 0;
        for (int b1=0; b1 < NbrBand1D; b1++)
        for (int p=0; p < WT1D.size_scale_np (b1); p++) WT1D(b1,p) = TabBand[b](j,i,z++); 
         Vect.init();
	 WT1D.recons(Vect);
         for (z=0; z < Nz; z++) TabBand[b](j,i,z) = Vect(z);
     }
   }
   
   // 2D wt 
   for (z=0; z < Nz; z++)
   {
      for (b=0; b < NbrBand2D; b++)
      {
         for (i=0; i < WT2D.size_band_nl(b); i++)
	 for (j=0; j < WT2D.size_band_nc(b); j++) WT2D(b,i,j) = TabBand[b](j,i,z);
      }   
      WT2D.recons(Frame);
      for (i=0; i < Ny; i++)
      for (j=0; j < Nx; j++) Data(j,i,z) = Frame(i,j);
   }
}

 /****************************************************************************/

 
