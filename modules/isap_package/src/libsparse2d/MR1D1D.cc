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
**    File:  MR1D1D.cc
**
*******************************************************************************
**
**    DESCRIPTION  multiresolution transform of cube, using 2D wavelet transform
**    -----------  followed (or not) by a 1D wavelet transform
**                 
**
******************************************************************************/

#include "IM_Obj.h"
#include "MR1D1D.h"

// #define DB_MR1D1D

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


int MR1D1D::mr_io_fill_header( fitsfile *fptr)
{
  int status = 0; // this means OK for cfitsio !!!
    /*****************************************************/
     /* write optional keyword to the header */
    /*****************************************************/

// cout << " DDDD " << (long) NbrBandX << " " << (long)NbrBandY << endl;

if ( ffpkyj(fptr, "Nx", (long)Nx,"Nx of the input cube",&status))
     PrintError( status );  
if ( ffpkyj(fptr,"Ny",(long)Ny,"Ny of the input cube",&status))
     PrintError( status );  
 if ( ffpkyj(fptr, "NSCALEX", (long) NbrScaleX, "Number of bands X", &status))
     PrintError( status );  
 if ( ffpkyj(fptr, "NSCALEY", (long)NbrScaleY, "Number of bands Y", &status))
     PrintError( status );  
 if ( ffpkyj(fptr, "Type_Tra", (long) WT_x.Type_Transform, 
                         (char*)StringTransf1D(WT_x.Type_Transform), &status))
     PrintError( status );  
 return(status);
} 


/****************************************************************************/

void  MR1D1D::read(char *Name)
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

    inc[0]=1;  inc[1]=1; 
 
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
     
    if (ffgkyj(fptr,"NSCALEX", &mon_long, comment, &status)) PrintError( status );
    int iNbrScaleX  = (int) mon_long;
    if (ffgkyj(fptr,"NSCALEY", &mon_long, comment, &status)) PrintError( status );
    int iNbrScaleY  = (int) mon_long;

    type_trans_1d  TT;
    if (ffgkyj(fptr,"Type_Tra", &mon_long, comment, &status))   PrintError( status );
    else TT = (type_trans_1d) mon_long;
    
    alloc(Nxi,Nyi,TT, iNbrScaleX, iNbrScaleY);
     
    fltarray Tab(nelements);
   //cout << "        Nl = " << (TabCF_Band[0][1]).nl() << " Nc = " << TabCF_Band[0][1].nc() << endl;

    Ptr = Tab.buffer();
    if ( ffgpve(fptr, 1, 1, nelements, nulval, Ptr, &anynul, &status))
             PrintError( status );

    int ind=2;
    for (int s=0; s <NbrBandX; s++) 
    for (int s1=0; s1 < NbrBandY; s1++)
    {
    
       int Nxb = (int) Tab(ind++);
       int Nyb = (int) Tab(ind++);
        //  cout << s << " " << s1 << " " << Nxb << " " << Nyb << " " << Nzb << endl;

        //cout << " s = " << s+1 << " b = " << b+1 << " Nl = " << Ny << " Nc = " << Nx << endl;
       //cout << "        Nl = " << (TabCF_Band[s][b]).nl() << " Nc = " << TabCF_Band[s][b].nc() << endl;
       for (int i=0;  i < Nxb; i++)
       for (int j=0;  j < Nyb; j++) (*this)(s,s1,i,j) = Tab(ind++);
   }
      
   if ( ffclos(fptr, &status) ) PrintError( status );
// cout << " end of read fits file " << endl;
}

/****************************************************************************/

void  MR1D1D::write(char *Name)
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
  for (int s=0; s < NbrBandX; s++)
  for (int s1=0; s1 < NbrBandY; s1++) 
  {
     Nelem += 2 +  TabSizeBandNx(s,s1)*TabSizeBandNy(s,s1);
  } 
  // cout << " NELEM = " << Nelem << endl;
  fltarray Data(Nelem);
  int ind=2;
  Data(0) = NbrBandX;
  Data(1) = NbrBandY;
 
  // Ifloat Band;
  for (int s=0; s < NbrBandX; s++)
  for (int s1=0; s1 < NbrBandY; s1++) 
  {
     Data(ind++) = TabSizeBandNx(s,s1);
     Data(ind++) = TabSizeBandNy(s,s1);
     
     for (int i=0;  i < TabSizeBandNx(s,s1); i++)
     for (int j=0;  j < TabSizeBandNy(s,s1); j++)  Data(ind++) = (*this)(s,s1,i,j);
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

inline float & MR1D1D::operator() (int sx, int s1, int i, int j) const
{    
#ifdef TESTIND
   if ( (i < 0) || (i >= size_band_nx(sx,s1)) ||
        (j < 0) || (j >= size_band_ny(sx,s1)) ||
 	(sx < 0) || (sx >= nbr_band_x()) ||
	(s1 < 0) || (s1 >= nbr_band_y()))
   {
      printf("Error: (sx,s1,i,j)=(%d,%d, %d,%d), size band = (%d,%d) , scale=(%d,%d) \n", sx,s1,i,j,size_band_nx(sx,s1), size_band_ny(sx,s1), (sx,s1), nbr_band_x(), nbr_band_y());
      exit(-1);
   }
#endif
  
   return TabBand[sx] (i,j+TabFirstPosBandNy(s1));
}


/****************************************************************************/

fltarray MR1D1D::get_band(int sx, int sy)
{
    int Nxb = size_band_nx(sx,sy);
    int Nyb = size_band_ny(sx,sy);
    fltarray *Dat_Return = NULL;

    Dat_Return = new fltarray(Nxb, Nyb);
    for (int i=0; i < Nxb; i++)
    for (int j=0; j < Nyb; j++) (* Dat_Return)(i,j) = (*this)(sx,sy,i,j);
    
    return (*Dat_Return);
}

/****************************************************************************/

void  MR1D1D::put_band(fltarray Band, int sx, int sy)
{
    int Nxb = size_band_nx(sx,sy);
    int Nyb = size_band_ny(sx,sy);
     
    for (int i=0; i < Nxb; i++)
    for (int j=0; i < Nyb; j++) (*this)(sx,sy,i,j) = Band(i,j);
}

/****************************************************************************/

void MR1D1D::alloc (int iNx, int iNy, type_trans_1d Trans1D, int Nsx, int Nsy, Bool NoAlloc)
{
   Nx = iNx;
   Ny = iNy;
   NbrScaleX = Nsx;
   NbrScaleY = Nsy;
   if (NbrScaleY < 2)
   {
      NbrScaleY = 1;
      Apply1DTrans = False;
   }
   else Apply1DTrans = True;
   ;
   Norm = NORM_L2;
   SB_Filter = F_MALLAT_7_9;
   Bord = I_CONT;
   U_Filter = DEF_UNDER_FILTER; 
   FilterAnaSynt *PtrFAS = NULL;
    if ((Trans1D == TO1_MALLAT) || (Trans1D == TU1_MALLAT))
    {
        FAS.Verbose = Verbose;
        FAS.alloc(SB_Filter);
        PtrFAS = &FAS;
    }
    type_lift LiftingTrans = DEF_LIFT;
    if (Trans1D == TO1_LIFTING)  WT_x.LiftingTrans = LiftingTrans;
    WT_x.Border = Bord;
// cout << "ALLOC " << StringTransf1D(Trans1D) << endl;
    WT_x.alloc (Nx, Trans1D, Nsx, PtrFAS, Norm);
    NbrBandX = WT_x.nbr_band();
     	       
    Bool Rebin=False;
    WT_y.U_Filter = U_Filter;
    type_trans_1d Trans1D_Y = TO1_MALLAT;
    if (Apply1DTrans == True)
    {
        WT_y.alloc (Ny, Trans1D_Y, Nsy, PtrFAS, Norm, Rebin);   
        NbrBandY = WT_y.nbr_band();
    } 
    else NbrBandY = 1;

//cout << "NbrBandX = " << NbrBandX << ", NbrBandY = " << NbrBandY << endl;

    
   if (NoAlloc == False) TabBand = new fltarray [NbrBandX];
   else TabBand=NULL;
   TabSizeBandNx.resize(NbrBandX, NbrBandY);
   TabSizeBandNy.resize(NbrBandX, NbrBandY);
   TabFirstPosBandNy.resize(NbrBandX);
   TabFirstPosBandNy(0) =0;
   NbrCoefX = 0;
   for (int b=0; b < NbrBandX; b++) 
   {
//   cout << " b = " << b+1 << ", NbrCoefX = " << NbrCoefX << ", SizeBand = " << WT_x.size_scale_np(b) << endl;
   
      int Npy = (Apply1DTrans == True) ? WT_y.size_ima_np (): Ny;
//   cout << "Npy = " << Npy << endl;
   
      if (NoAlloc == False)  TabBand[b].alloc(WT_x.size_scale_np(b),  Npy);
      
      NbrCoefX +=  WT_x.size_scale_np(b);
//       cout << "NbrCoefX = " << NbrCoefX << endl;
   
      for (int b1=0; b1 < NbrBandY; b1++)
      {
         TabSizeBandNx(b,b1) = WT_x.size_scale_np(b);
         if (Apply1DTrans == True)  TabSizeBandNy(b,b1) = WT_y.size_scale_np(b1);
	     else TabSizeBandNy(b,b1) = Ny;
      }
   }
   NbrCoefY = Ny;
   if (Apply1DTrans == True)  
     for (int b1=1; b1 < NbrBandY; b1++) TabFirstPosBandNy(b1) = TabFirstPosBandNy(b1-1) + TabSizeBandNy(0,b1-1);

// cout << "END ALLOC " << endl;
}

/****************************************************************************/

void MR1D1D::transform_to_vectarray(fltarray &Data, fltarray &TabVect)
{
   int Nsx = nbr_coef_x();
   int Nsy = nbr_coef_y();
#ifdef DB_MR1D1D
   cout << "ALLOC transform_to_vectarray " << Nsx << " " << Nsy << endl;
#endif
   TabVect.alloc(nbr_coef_x(), nbr_coef_y());
   int x,b,y,i;
   Nx = Data.nx();
   Ny = Data.ny();
   fltarray Frame(Nx);
   fltarray Vect(Ny);
   intarray TabPosBand(NbrBandX);
    
   // 2D wt transform per frame
   for (y=0; y < Ny; y++)
   {
         int Pix=0;
        for (i=0; i < Nx; i++) Frame(i) = Data(i,y);
        WT_x.transform(Frame);
        for (b=0; b < NbrBandX; b++)
        {
           for (i=0; i < WT_x.size_scale_np (b); i++)  TabVect(Pix++,y) = WT_x(b,i);
	       if (y == 0) TabPosBand(b) = Pix;  
        }
    }
#ifdef DB_MR1D1D
 cout << "1D_Y " << NbrBandY << endl;
#endif
   // 1D wt 
   if (NbrBandY >= 2)
   {
   	   for (b=0; b < NbrBandX; b++)
       for (i=0; i < WT_x.size_scale_np (b); i++)
       {
           for (y=0; y< Ny; y++) Vect(y) = TabBand[b](i,y);   TabVect(TabPosBand(b)+i*WT_x.size_scale_np(b)+i, y); 
           WT_y.transform(Vect);
           y = 0;
           for (int b1=0; b1 < NbrBandY; b1++)
           for (int p=0; p < WT_y.size_scale_np (b1); p++)  TabVect(TabPosBand(b)+i* WT_x.size_scale_np(b)+i,y) = WT_y(b1,p);
        } 	
   }    
#ifdef DB_MR1D1D
 cout << "END TRANS " << endl;
#endif
}

/****************************************************************************/

void MR1D1D::transform (fltarray &Data)
{
   int i,j,b,y;
   Nx = Data.nx();
   Ny = Data.ny();
   fltarray Frame(Nx);
   fltarray Vect(Ny);
   
   // 2D wt transform per frame
   for (y=0; y < Ny; y++)
   {
      for (i=0; i < Nx; i++) Frame(i) = Data(i,y);
      WT_x.transform(Frame);
      for (b=0; b < NbrBandX; b++)
      for (i=0; i < WT_x.size_scale_np (b); i++) TabBand[b](i,y) = WT_x(b,i);
    }
 
   // 1D wt 
   if (NbrBandY >= 2)
   {
     for (b=0; b < NbrBandX; b++)
     for (i=0; i < WT_x.size_scale_np (b); i++)
      {
        for (y=0; y< Ny; y++) Vect(y) = TabBand[b](i,y);
        WT_y.transform(Vect);
        y = 0;
        for (int b1=0; b1 < NbrBandY; b1++)
        {
         for (int p=0; p < WT_y.size_scale_np (b1); p++) TabBand[b](i,y++) = WT_y(b1,p); 
        }
      }
   }
}

/****************************************************************************/

void MR1D1D::recons (fltarray &Data)
{
  if ((Data.nx() != Nx) || (Data.ny() != Ny)) Data.resize(Nx, Ny); 
   int i,j,b,y;
   Nx = Data.nx();
   fltarray Frame(Nx);
   fltarray Vect(Ny);
     
   // 1D wt  y
   if (NbrBandY >= 2)
   {
      for (b=0; b < NbrBandX; b++)
      for (i=0; i < WT_x.size_scale_np (b); i++)
      {
    	y = 0;
        for (int b1=0; b1 < NbrBandY; b1++)
        for (int p=0; p < WT_y.size_scale_np (b1); p++)  WT_y(b1,p) = TabBand[b](i, y++); 
        Vect.init();
	    WT_y.recons(Vect);
        for (y=0; y < Ny; y++) TabBand[b](i,y) = Vect(y);
     }
   }
   
   // 1D wt  x
   for (y=0; y < Ny; y++)
   {
      for (b=0; b < NbrBandX; b++)
      for (i=0; i < WT_x.size_scale_np (b); i++) WT_x(b,i) = TabBand[b](i,y);
      WT_x.recons(Frame);
      for (i=0; i < Nx; i++)  Data(i,y) = Frame(i);
   }
}

 /****************************************************************************/

 
