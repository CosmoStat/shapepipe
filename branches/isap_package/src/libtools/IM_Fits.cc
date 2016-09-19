
/******************************************************************************
**                   Copyright (C) 1995 by CEA
*******************************************************************************
**
**    UNIT
**
**    Version: 3.21
**
**    Author: Jean-Luc Starck
**
**    Date:  96/05/07 
**    
**    File:  IM_Fits.cc
**
**    History: 12 March 1997 modified by R Gastaud for fits_read_fltarr
**
*******************************************************************************
**
** void fits_read_float(char *File_Name, Ifloat & Image, fitsstruct *Header)
**
** Read a fits image from a file
**
*******************************************************************************
**
** void fits_write_float(char *File_Name, Ifloat & Image, fitsstruct *Header)
**
** Write a fits image to a file
**
******************************************************************************/


#include"GlobalInc.h"
#include"IM_IOTools.h"
#include"IM_IO.h"
#include"Array.h"

#include "fitsio.h"
#include "fitsio2.h"

extern void readimagehead(fitsstruct *field);
extern void writeimagehead(fitsstruct *field);
extern void readdataf(fitsstruct *field, float *ptr);
extern void readdataf(FILE *file, char *filename, int bitpix, int Nelem, 
                      float *ptr, float bscale=1., float bzero=0.);

extern void writedataf(fitsstruct *field, float *ptr);
extern void writedataf(FILE *file, char *filename, int bitpix, int Nelem, 
                      float *ptr, float bscale=1., float bzero=0.);
extern void FitsPrintError( int status);

int FITS_HDU_Number = 0;

/************************************/

char *fitsname(char * NameStep)
{
    char Sreturn[MAXCHAR];
    char Name[MAXCHAR];
    char *ValRet;

    strcpy(Name,NameStep);
    if (std_inout(NameStep) == False) 
    {
       if (strstr(Name, ".fit") != NULL)  strcpy(Sreturn, Name);
       else if (strstr(Name, ".fts") != NULL)  strcpy(Sreturn, Name);
       else if (strstr(Name, ".mr") != NULL)  strcpy(Sreturn, Name);
       else if (strstr(Name, ".rad") != NULL)  strcpy(Sreturn, Name);
       else if (strstr(Name, ".rid") != NULL)  strcpy(Sreturn, Name);
       else if (strstr(Name, ".cur") != NULL)  strcpy(Sreturn, Name);
       else if (strstr(Name, ".bet") != NULL)  strcpy(Sreturn, Name);
       else if (strstr(Name, ".FIT") != NULL)  strcpy(Sreturn, Name);
       else sprintf(Sreturn, "%s.%s", Name, "fits");
       ValRet = strdup(Sreturn);
    }
    else ValRet = strdup(NameStep);
    return (ValRet);
} 

/************************************/

void init_fits_struct(fitsstruct *Ptr, int Nl, int Nc)
{
   initfield(Ptr);
   Ptr->naxis = 2;
   Ptr->bitpix = -32;
   Ptr->width =  Nc;
   Ptr->height = Nl;
   Ptr->TabAxis[0] = Nc;
   Ptr->TabAxis[1] =  Nl;
   Ptr->npix = (Ptr->TabAxis)[0]*(Ptr->TabAxis)[1];
}

/************************************/

void radec(double Ra, double Dec, int &ihr, int &imin, double &xsec,
                                  int &ideg,int &imn, double &xsc)
{
  double ra = Ra;
  double dec = Dec;
  double xmin,xmn;

  // Compute RA

  ra -= (int) ra / 360;          // Make sure between 0 and 24 hours
  if (ra < 0) ra += 360;
  
  ihr = (int)(ra/15.);
  xmin = ABS(ra*4.0-ihr*60.0);
  imin = (int) xmin;
  xsec = (xmin-imin)*60.0;

  // Compute Dec
  ideg = (int) dec;
  xmn = ABS(dec-ideg)*60.0;
  imn = (int) xmn;
  xsc = (xmn-imn)*60.0;

  // Now test for the special case of zero degrees
  if ((ideg == 0.) && (dec < 0))
  {
     if (imn != 0) 
           imn = imn - 2*imn;
     else  xsc = xsc - 2*xsc;
  }
}
  
/************************************/

void adxy(fitsstruct *FitsHeader, double Ra, double Dec, double *X, double *Y)
{
  int status,stat;
  stat = ffxypx(Ra,Dec, FitsHeader->crvalx,  FitsHeader->crvaly, 
                FitsHeader->crpixx, FitsHeader->crpixy, 
                FitsHeader->cdeltx,  FitsHeader->cdelty, 
                FitsHeader->crotay, FitsHeader->CoordType, 
                X,Y,&status); 
  *X -=  1.;
  *Y -=  1.;
}

/************************************/

void xyad(fitsstruct *FitsHeader, double X, double Y, double *Ra, double *Dec)
{
   int status,stat;
   stat = ffwldp(X+1.,Y+1., FitsHeader->crvalx,  FitsHeader->crvaly, 
                 FitsHeader->crpixx, FitsHeader->crpixy, 
                 FitsHeader->cdeltx,  FitsHeader->cdelty, 
                 FitsHeader->crotay, FitsHeader->CoordType, 
                 Ra, Dec, &status);   
}

/************************************/

int cfitstio_read_celestial_coord(char *File_Name, fitsstruct *FitsHeader)
{
   fitsfile *fptr; 
   int status=0,stat=0;
   double xrval,yrval,xrpix,yrpix,xinc,yinc,rot;
   char ctype[5];
   
  if (io_detect_format(File_Name) != F_FITS)
  {
     cerr << "Error: input file is not a fits format file ... " << endl;
     exit(-1);
  }
  
  if (fits_open_file(&fptr, File_Name, (int) READONLY, &stat)) 
  {
     fprintf(stderr,"Error: CFITSIO package cannot open file %s\n",File_Name);
     fprintf(stderr,"  status = %d\n", stat);
     exit (-1);
  }
    xrval =  0.;
    yrval =  0.;
    xrpix =  0.;
    yrpix =  0.;
    xinc =   0.;
    yinc =   0.;
    rot =    0.;
  
    stat = ffgics(fptr, &xrval, &yrval, &xrpix,
                 &yrpix, &xinc, &yinc, &rot, ctype, &status);
  
    if (status)
    {
       fprintf(stderr,"Warning: CFITSIO package cannot read celestial coordinates ...\n");
       fprintf(stderr,"   status = %d\n", status);
    }
    
    FitsHeader->crpixx = xrpix;
    FitsHeader->crpixy = yrpix;
    FitsHeader->crvalx = xrval;
    FitsHeader->crvaly = yrval;
    FitsHeader->cdeltx = xinc;
    FitsHeader->cdelty = yinc;
    FitsHeader->crotay = rot;
    strcpy(FitsHeader->CoordType,ctype);
    strcpy(FitsHeader->ctypex,ctype);
    strcpy(FitsHeader->ctypey,ctype);
   stat = 0;
   if (fits_close_file(fptr, &stat) )
   {
      fprintf(stderr,"Error closing file %s\n",File_Name);
      exit (-1);
   }   
   return status;
}



/************************************/
/*--------------------------------------------------------------------------*/
void FitsPrintError( int status)
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

void fits_read_header(char *File_Name, fitsstruct *Header)
{
    char Name[MAXCHAR];

    initfield(Header);
    Header->filename = strdup(File_Name);
    if (!(Header->file = fopen(Header->filename, "rb")))
    {
        sprintf (Name, "%s.fits", File_Name);
        if (!(Header->file = fopen(Name, "rb")))
        {
            fprintf(stderr,"Error reading file %s\n",File_Name);
            exit (-1);
        }
        else Header->filename = strdup(Name);
    }
   readimagehead(Header);   
   fclose(Header->file);
}

/************************************/

void fits_read_block(char *filename, Ifloat & Image, int Indi, int Indj, Bool NoBscale)
{
    // for fits
    fitsfile *fptr;           /* pointer to the FITS file */
    int status=0, hdutype ;
    long hdunum;
    int naxis;
    long naxes[3];
    float nulval = 0.;
    int anynul = 0;
    long inc[3];
    long nelements = 0 ; // naxes[0] * naxes[1] in the image
     
    float *Ptr;
    inc[0]=1;  inc[1]=1; inc[2]=1;

  /* open the file */
    status = 0;         /* initialize status before calling fitsio routines */
    if ( fits_open_file(&fptr, fitsname(filename), (int) READONLY, &status) ) 
         FitsPrintError( status );
	 
    hdunum = 1;  /*read  table */
    if ( fits_movabs_hdu(fptr, hdunum, &hdutype, &status) ) /* move to the HDU */
           FitsPrintError( status );

    if (NoBscale == True)
       if (fits_set_bscale(fptr, (double) 1.,  (double) 0., &status) ) 
           FitsPrintError( status );
	   
    int simple, bitpix, extend;
    long pcount, gcount;
    if ( fits_read_imghdr(fptr, 3, &simple, &bitpix, &naxis, naxes, &pcount,
            &gcount, &extend, &status)) /* move to the HDU */
           FitsPrintError( status );

    long fpixels[3];
    long lpixels[3];
    fpixels[0] = 1 + Indj;
    fpixels[1] = 1 + Indi;
    fpixels[2] = 1;
    lpixels[0] = fpixels[0] +  Image.nc() - 1;
    lpixels[1] = fpixels[1] +  Image.nl() - 1;
    lpixels[2] = 1;
    Ptr = Image.buffer();
    nelements =  Image.n_elem();
//    cout << "subset: " << fitsname(filename) << " " << naxis << endl;
//     cout << "naxes[0] " << naxes[0] << endl;
//     cout << "naxes[1] " << naxes[1] << endl;
//     cout << "fpixels[0] " << fpixels[0] << endl;
//     cout << "fpixels[1] " << fpixels[1] << endl;
//     cout << "lpixels[0] " << lpixels[0] << endl;
//     cout << "lpixels[1] " << lpixels[1] << endl;
 
    if ( fits_read_subset_flt(fptr, 0 , naxis, naxes, fpixels, lpixels, inc, nulval,
                  Ptr, &anynul, &status))
                 FitsPrintError( status );
    if ( fits_close_file(fptr, &status) ) FitsPrintError( status );
}

/************************************/

void fits_write_block(char *filename, Ifloat & Image, int Indi, int Indj,
                      Bool NoBscale)
{
    // for fits
    fitsfile *fptr;           /* pointer to the FITS file */
    int status=0;
    int naxis=2;
    long naxes[3];
     long inc[3];
    long nelements = 0 ; // naxes[0] * naxes[1] in the image
     
    float *Ptr;
    inc[0]=1;  inc[1]=1; inc[2]=1;

  /* open the file */
    status = 0;         /* initialize status before calling fitsio routines */
    if ( fits_open_file(&fptr, fitsname(filename), (int) READWRITE, &status) ) 
         FitsPrintError( status );

    int simple, bitpix, extend;
    long pcount, gcount;
    if ( fits_read_imghdr(fptr, 3, &simple, &bitpix, &naxis, naxes, &pcount,
            &gcount, &extend, &status)) /* move to the HDU */
           FitsPrintError( status );
	     
    if (NoBscale == True)
       if (fits_set_bscale(fptr, (double) 1.,  (double) 0., &status) ) 
           FitsPrintError( status );
	 
    long group=0;
    long fpixels[3];
    long lpixels[3];
    fpixels[0] = 1 + Indj;
    fpixels[1] = 1 + Indi;
    fpixels[2] = 1;
    lpixels[0] = fpixels[0] +  Image.nc() - 1;
    lpixels[1] = fpixels[1] +  Image.nl() - 1;
    lpixels[2] = 1;
    Ptr = Image.buffer();
    nelements =  Image.n_elem();

    if ( fits_write_subset_flt(fptr, group, naxis, naxes, fpixels, lpixels,  
                               Ptr,  &status))
                 FitsPrintError( status );
    if ( fits_close_file(fptr, &status) ) FitsPrintError( status );
}

/************************************/

void fits_write_block(char *filename, Iint & Image, int Indi, int Indj,
                      Bool NoBscale)
{
    // for fits
    fitsfile *fptr;           /* pointer to the FITS file */
    int status=0;
    int naxis=2;
    long naxes[3];
    long inc[3];
    long nelements = 0 ; // naxes[0] * naxes[1] in the image
     
    int *Ptr;
    inc[0]=1;  inc[1]=1; inc[2]=1;

  /* open the file */
    status = 0;         /* initialize status before calling fitsio routines */
    if ( fits_open_file(&fptr, fitsname(filename), (int) READWRITE, &status) ) 
         FitsPrintError( status );
  
    int simple, bitpix, extend;
    long pcount, gcount;
    if ( fits_read_imghdr(fptr, 3, &simple, &bitpix, &naxis, naxes, &pcount,
            &gcount, &extend, &status)) /* move to the HDU */
           FitsPrintError( status );
	   
    if (NoBscale == True)
       if (fits_set_bscale(fptr, (double) 1.,  (double) 0., &status) ) 
           FitsPrintError( status );
	 
    long  group=0;
    long fpixels[3];
    long lpixels[3];
    fpixels[0] = 1 + Indj;
    fpixels[1] = 1 + Indi;
    fpixels[2] = 1;
    lpixels[0] = fpixels[0] +  Image.nc() - 1;
    lpixels[1] = fpixels[1] +  Image.nl() - 1;
    lpixels[2] = 1;
    Ptr = Image.buffer();
    nelements =  Image.n_elem();

    if ( fits_write_subset_int(fptr, group, naxis, naxes, fpixels, lpixels,
                  Ptr,  &status))
                 FitsPrintError( status );
    if ( fits_close_file(fptr, &status) ) FitsPrintError( status );
}

/************************************/

void fits_read_block(char *filename, Iint & Image, int Indi, int Indj, Bool NoBscale)
{
    // for fits
    fitsfile *fptr;           /* pointer to the FITS file */
    int status=0, hdutype ;
    long hdunum;
    int naxis;
    long naxes[3];
    long nulval = 0;
    int anynul = 0;
    long inc[3];
    long nelements = 0 ; // naxes[0] * naxes[1] in the image
     
    int *Ptr;
    inc[0]=1;  inc[1]=1; inc[2]=1;

  /* open the file */
    status = 0;         /* initialize status before calling fitsio routines */
    // READONLY or READWRITE
    if ( fits_open_file(&fptr, filename, (int) READONLY, &status) ) 
         FitsPrintError( status );
	 
    hdunum = 1;  /* read  table */
    if ( fits_movabs_hdu(fptr, hdunum, &hdutype, &status) ) /* move to the HDU */
           FitsPrintError( status );

    if (NoBscale == True)
       if (fits_set_bscale(fptr, (double) 1.,  (double) 0., &status) ) 
           FitsPrintError( status );
	   
    int simple, bitpix, extend;
    long pcount, gcount;
    int maxdim = 3;
    if ( fits_read_imghdr(fptr, maxdim, &simple, &bitpix, &naxis, naxes, &pcount,
            &gcount, &extend, &status)) /* move to the HDU */
           FitsPrintError( status );

    long fpixels[3];
    long lpixels[3];
    fpixels[0] = 1 + Indj;
    fpixels[1] = 1 + Indi;
    fpixels[2] = 1;
    lpixels[0] = fpixels[0] +  Image.nc() - 1;
    lpixels[1] = fpixels[1] +  Image.nl() - 1;
    lpixels[2] = 1;
    Ptr = Image.buffer();
    nelements = Image.nl()*Image.nc();
    if ( fits_read_subset_int(fptr, 0 , naxis, naxes, fpixels, lpixels, inc, nulval,
                  Ptr, &anynul, &status))
                 FitsPrintError( status );
    if ( fits_close_file(fptr, &status) ) FitsPrintError( status );
}
 
/************************************/

void fits_read_float(char *File_Name, Ifloat & Image, fitsstruct *Header)
{
    char *Name = strdup(fitsname(File_Name));
    int Nl, Nc;
    fitsfile *fptr;                       
    int status,  nfound, anynull;
    int hdutype;
    long  naxes[3];
    long group, fpixel;
    //int bitpix   =  0;   /* -32 for 32-bit float */
    long Nelem;
    float  nullval;
    int fits_read_com( fitsfile *fptr, char *keyname,       
                      char **value );  /* O - pointer to keyword value  */
//     int fits_get_img_dim (fitsfile *fptr, int *naxis, int *status);
//     int fits_get_img_size (fitsfile *fptr, int maxdim,long *naxes, int *status);
//     int fits_get_img_type(fitsfile *fptr, int *bitpix, int *status);
//     int fits_copy_header(fitsfile *infptr, fitsfile *outfptr, int *status);
    status = 0;
// cout << "READ IN " << File_Name << endl;

    if ( fits_open_file(&fptr, Name, READONLY, &status) ) /* open the image */
    {
       printf("Error: cannot open file %s\n", Name);
       exit(status);
    }
    
//     if ( fits_get_img_type(fptr, &bitpix, &status) )
//     {
//        printf("Error: in fits_get_img_type ...\n");
//        exit(status);
//     }
//     // fits_get_img_param : ffgipr 
//     if ( fits_get_img_dim(fptr, &nfound, &status) )
//     {
//         printf("Error: pb NAXIS ...\n");
//         exit(status);
//     }

//     if (nfound == 0)
    {
        // printf(" primary image look for extension \n");
        if (fits_movabs_hdu(fptr, 1+FITS_HDU_Number, &hdutype, &status) )  
        {
            printf("Error: no image extension ...");
            exit(status);
        }
	// printf(" 1 primary image look for extension \n");
        if (hdutype != IMAGE_HDU)  
        {
            printf("Error: no image in this HDU\n");
            exit(-1);
        }
	// printf(" 2 primary image look for extension \n");
        if ( fits_get_img_dim(fptr, &nfound, &status) )
        {
            printf("Error: pb NAXIS...\n");
            exit(status);
        }
    }       
    naxes[2] = 0; naxes[0] = 0; naxes[1] = 0;
    if (fits_get_img_size(fptr, 3, naxes,  &status) )  
    {
       printf("Error: pb NAXIS  ... \n");
       exit(status);
    }
    Header->naxis = nfound;
   
   if (Header->naxis != 2)
   {
     cerr << "Error: the input image is not a 2D array ... " << endl;
     exit(-1);
   }
   // initfield(Header);
   Nl = naxes[1];
   Nc = naxes[0];

   if (Nl*Nc < 1)
   {
     cerr << "Error: the input image is not a 2D array ... " << endl;
     exit(-1);
   }
   
   if (Nl == 1) 
   {
     cerr << "Error: number of lines must be > 1 ..." << endl;
     exit(-1);
   }
   if (Nc == 1) 
   {
     cerr << "Error: number of columns must be > 1 ..." << endl;
     exit(-1);
   }

   Image.alloc(Nl, Nc, File_Name);
    
//    {
//       printf("file %s naxes %d %d %d %d \n", Name, nfound , 2, Image.nc() , Image.nl() );
//       printf(" nb elem data %d \n", Image.n_elem());
//    }
   group    = 1;
   fpixel   = 1;
   nullval  = 0.;                
   
   Nelem = (long) Image.n_elem();
   if(fits_read_img(fptr, TFLOAT, fpixel, Nelem,(void *) &nullval, 
                 Image.buffer(), &anynull, &status) )
   {
      printf("Error: cannot read in file %s.\n", File_Name);
      exit(status);
   }
#ifndef WINDOWS
// replacing tmpnam by mkstemp
//   char FitsName[L_tmpnam]; // = "xxhdmr.fit";
//   char *Res = tmpnam(FitsName); // replacing tmpnam by mkstemp
	char FitsName[] = "hdmrXXXXXX";
	int fd;   
	fd = mkstemp(FitsName);
#else
  char *FitsName = "xxhdmr.fit";
#endif
    
   fitsfile *fptr1;                       
   remove(FitsName);
   if ( fits_create_file(&fptr1, FitsName, &status))      
   {
          printf("Error: cannot open file %s %d \n", FitsName, status);
          exit(-1);
   }
   fits_copy_header(fptr, fptr1, &status);
   if ( fits_close_file(fptr, &status) )          
   {
        printf("Error: cannot close %s\n",Name);
        exit(status);
   }
   fits_close_file(fptr1, &status);          
//    {
//         printf("Error: cannot close %s\n",FitsName);
//         exit(status);
//    }   
   Header->filename = strdup(FitsName);
   Header->file = fopen(Header->filename, "rb");
   if (Header->file ==NULL)
   {
      printf("Error: cannot open %s\n",FitsName);
      exit(status);
   }
   readimagehead(Header);
   fclose(Header->file);
   remove(FitsName);
}

/*****************************************************************/
// void fits_read_float(char *File_Name, Ifloat & Image, fitsstruct *Header)
// {
//     int Nl, Nc;
//  
//     // initfield(Header);
//     Header->filename = strdup(fitsname(File_Name));
//     Header->file = fits_file_des_in(File_Name);
//    readimagehead(Header);
//    Nl = Header->height;
//    Nc = Header->width;
//    
//    if (Header->naxis != 2)
//    {
//      cerr << "Error: the input image is not a 2D array ... " << endl;
//      exit(-1);
//    }
//    
//    if (Nl*Nc < 1)
//    {
//      cerr << "Error: the input image is not a 2D array ... " << endl;
//      exit(-1);
//    }
//    
//    if (Nl == 1) 
//    {
//      cerr << "Error: number of lines must be > 1 ..." << endl;
//      exit(-1);
//    }
//    if (Nc == 1) 
//    {
//      cerr << "Error: number of columns must be > 1 ..." << endl;
//      exit(-1);
//    }
//    
//    Image.alloc(Nl, Nc, File_Name);
//    readdataf(Header, Image.buffer());
//    if (Header->file != stdin) fclose(Header->file);
// }

/*****************************************************************/
void fits_write_header(char *File_Name, fitsstruct *Header)
{
     Header->filename = fitsname(File_Name);
     if (!(Header->file = fopen(Header->filename, "wb")))
     {
        fprintf(stderr,"Error writing on file %s\n",File_Name);
        exit (-1);
     }
     writeimagehead(Header);
     fclose(Header->file);
}

/*****************************************************************/

FILE *fits_file_des_in(char *fname)
{
   FILE *fp;
   
   if (std_inout(fname) == True)  fp = stdin;
   else 
   {
      fp = fopen(fitsname(fname), "rb");
      if (!fp) 
      {
        cerr << "Unable to open file: " <<  fname  << endl;
        exit(-1);
      }
   }
   return fp;
}
FILE *fits_file_des_out(char *fname)
{
   FILE *fp;
   
   if (std_inout(fname) == True) fp = stdout;
   else 
   {
      fp = fopen(fitsname(fname), "wb");
      if (!fp) 
      {
        cerr << "Unable to open file: " <<  fname  << endl;
        exit(-1);
      }
   }
   return fp;
} 
/*****************************************************************/

void fits_write_float(char *File_Name, Ifloat & Image, fitsstruct *Header)
{
    Header->file = fits_file_des_out(File_Name);
    writeimagehead(Header);
    writedataf(Header, Image.buffer());
    if (Header->file != stdout)  fclose(Header->file);
}

/************************************/

void fits_read_int(char *File_Name, Iint & Image, fitsstruct *Header, Bool NoBscale)
{
    int Nl, Nc;
    extern void  readdatai(fitsstruct *field, int *ptr);
    float bscale;
    float bzero;

    // initfield(Header);
    Header->filename = strdup(File_Name);
    Header->file = fits_file_des_in(File_Name);
   readimagehead(Header);
   Nl = Header->height;
   Nc = Header->width;
   bscale = Header->bscale;
   bzero = Header->bzero;
   if (NoBscale == True) 
   {
      Header->bscale=1.;
      Header->bzero=0.;
   }
   Image.alloc(Nl, Nc, File_Name);
   readdatai(Header, Image.buffer());
   if (NoBscale == True)
   {
       Header->bscale=bscale;
       Header->bzero=bzero;
   }
   if (Header->file != stdin) fclose(Header->file);
}

/*****************************************************************/

void fits_write_int(char *File_Name, Iint & Image, fitsstruct *Header)
{
    extern void  writedatai(fitsstruct *field, int *ptr);

    Header->filename = fitsname(File_Name);
    Header->file = fits_file_des_out(File_Name);

    writeimagehead(Header);
    writedatai(Header, Image.buffer());
    if (Header->file != stdout) fclose(Header->file);
}

/************************************/



