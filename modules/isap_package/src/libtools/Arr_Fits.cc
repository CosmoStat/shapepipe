/******************************************************************************
**                   Copyright (C) 1998 by CEA
*******************************************************************************
**
**    UNIT
**
**    Version: 1.0
**
**    Author: Jean-Luc Starck
**
**    Date:  98/07/17 
**    
**    File:  Arr_Fits.cc
**
**
******************************************************************************/


#include"GlobalInc.h"
#include"IM_IOTools.h"
#include"IM_IO.h"
#include"Array.h"
#include <sys/time.h>

extern void readimagehead(fitsstruct *field);
extern void writeimagehead(fitsstruct *field);
extern void readdataf(fitsstruct *field, float *ptr);
extern void readdataf(FILE *file, char *filename, int bitpix, int Nelem, 
                      float *ptr, float bscale=1., float bzero=0.);

extern void writedataf(fitsstruct *field, float *ptr);
extern void writedataf(FILE *file, char *filename, int bitpix, int Nelem, 
                      float *ptr, float bscale=1., float bzero=0.);

/************************************/

void fits_write_fltarr(char *File_Name, fltarray &Mat)
{
   fitsstruct HD;

   initfield(&HD);
   HD.bitpix = -32;
   HD.width =  Mat.nx();
   HD.height = Mat.ny();
   HD.naxis = Mat.naxis();
   HD.npix = Mat.n_elem();
   for (int i=0; i < HD.naxis; i++) HD.TabAxis[i] = Mat.axis(i+1);
   fits_write_fltarr(File_Name, Mat, &HD);
}

/************************************/

void fits_write_fltarr(char *File_Name, fltarray &Mat, fitsstruct *FitsHeader)
{
    if (FitsHeader != NULL)
    {
        FitsHeader->filename = fitsname(File_Name);
        FitsHeader->file = fits_file_des_out(File_Name);
        writeimagehead(FitsHeader);
        writedataf(FitsHeader, Mat.buffer());
        fclose(FitsHeader->file);
    }
    else
    {
        fitsstruct HD;
        initfield(&HD);
        HD.bitpix = -32;
        HD.width =  Mat.nx();
        HD.height = Mat.ny();
        HD.naxis = Mat.naxis();
        HD.npix = Mat.n_elem();
        for (int i=0; i < HD.naxis; i++) HD.TabAxis[i] = Mat.axis(i+1);
        HD.filename = fitsname(File_Name);
        HD.file = fits_file_des_out(File_Name);
        writeimagehead(&HD);
        writedataf(&HD, Mat.buffer());
        fclose(HD.file);  
   }      
}

/************************************/

void fits_read_fltarr(char *File_Name, fltarray & Data)
{
   fitsstruct FitsHeader;
   fits_read_fltarr(File_Name, Data, &FitsHeader, 1);
}


/************************************/

void fits_read_fltarr(char *File_Name, fltarray & Data, fitsstruct *FitsHeader)
{
   fits_read_fltarr(File_Name, Data, FitsHeader, 1);
}

/************************************/

/************************************/

bool file_exists(const char * filename)
{
    if (FILE * file = fopen(filename, "r"))
    {
        fclose(file);
        return true;
    }
    return false;
}

void fits_read_fltarr(char *File_Name, fltarray & Data, fitsstruct *Header,
 int openflag)
{
    char *Name = fitsname(File_Name);
    fitsfile *fptr;                       
    int status,  nfound, anynull;
    int hdutype;
    long  naxes[3];
    long group, fpixel;
    int bitpix   =  0;   /* -32 for 32-bit float */
    long Nelem;
    float  nullval;
    int fits_read_com( fitsfile *fptr, char *keyname,       
                      char **value );  /* O - pointer to keyword value  */
//     int fits_get_img_dim (fitsfile *fptr, int *naxis, int *status);
//     int fits_get_img_size (fitsfile *fptr, int maxdim,long *naxes, int *status);
//     int fits_get_img_type(fitsfile *fptr, int *bitpix, int *status);
//     int fits_copy_header(fitsfile *infptr, fitsfile *outfptr, int *status);
    status = 0;

    if ( fits_open_file(&fptr, Name,READONLY, &status) ) /* open the image */
    {
       printf("Error: cannot open file %s\n", Name);
       exit(status);
    }
    
    if ( fits_get_img_type(fptr, &bitpix, &status) )
    {
       printf("Error: in fits_get_img_type ...\n");
       exit(status);
    }
    printf("BITPIX=%d\n",bitpix);
    // fits_get_img_param : ffgipr 
    if ( fits_get_img_dim(fptr, &nfound, &status) )
    {
        printf("Error: pb NAXIS ...\n");
        exit(status);
    }
     
    if (nfound == 0)
    {
        // printf(" no primary image look for extension \n");
        if ( fits_movabs_hdu(fptr, 2, &hdutype, &status) ) // move to 2nd HDU 
        {
            printf("Error: no image extension ...\n");
            exit(status);
        }
        if (hdutype != IMAGE_HDU)  
        {
            printf("Error: no image in this HDU\n");
            exit(-1);
        }
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
    printf("Naxis=[%d,%d,%d]\n",naxes[0],naxes[1],naxes[2]);
    switch(nfound)
    {
       case 1 : Data.alloc(naxes[0] ); break;
       case 2 : Data.alloc(naxes[0], naxes[1] ); break;
       case 3 : Data.alloc(naxes[0], naxes[1] , naxes[2] ); break;
       default : printf( "Error: dimension =  %d \n", nfound); exit(-1);break;
    }    
  
   group    = 1;
   fpixel   = 1;
   nullval  = 0.;                
   
   Nelem = (long) Data.n_elem();
   if(fits_read_img(fptr, TFLOAT, fpixel, Nelem,(void *) &nullval, 
                 Data.buffer(), &anynull, &status) )
   {
      printf("Error: cannot read in file %s\n", File_Name);
      exit(status);
   }
   
// replacing tmpnam by mkstemp
//  char FitsName[L_tmpnam]; // = "xxhdmr.fit";
//   char *Res = tmpnam(FitsName);
	//char FitsName[] = "hdmrXXXXXX";
	//int fd;   
	//fd = mkstemp(FitsName);
	
	// mkstemp crashes after 1024 calls
	// so we create the temp name ourself
	char FitsName[] = "hdmrXXXXXXXXXX";
	timeval tim;
	gettimeofday(&tim, NULL);
	long int tt=long(tim.tv_sec*1000000L+tim.tv_usec)%long(1e10L);
	char ttt[16];
	sprintf(ttt,"%10ld",tt);
	
	do
	{
		strcpy(FitsName,"hdmr");
		strcat(FitsName,ttt);
		tt+=1;
	}
	while(file_exists(FitsName));
	
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
	if(fits_close_file(fptr1, &status))
	{
		printf("Error: cannot close %s\n",FitsName);
		exit(status);
	}	
   
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
   free(Name);
}

// void fits_read_fltarr(char *File_Name, fltarray & Data, fitsstruct *FitsHeader,
//  int openflag)
// {
//     // initfield(FitsHeader);
//     FitsHeader->filename = strdup(fitsname(File_Name));
//     FitsHeader->file = fits_file_des_in(File_Name);
// //     if (!(FitsHeader->file = fopen(FitsHeader->filename, "rb")))
// //     {
// //         fprintf(stderr,"Error reading file %s\n",File_Name);
// //         exit (-1);
// //     }
//     readimagehead(FitsHeader);
// 
//     switch(FitsHeader->naxis)
//     {
//        case 0: cerr << "Error: bad NAXIS keywords in " << File_Name << endl;
//                exit(0);
//                break;       
//        case 1: Data.alloc((FitsHeader->TabAxis)[0]);
//                break;       
//        case 2: Data.alloc((FitsHeader->TabAxis)[0],
//                                (FitsHeader->TabAxis)[1]);
//                break;    
//        case 3: Data.alloc((FitsHeader->TabAxis)[0],
//                                (FitsHeader->TabAxis)[1],
//                                (FitsHeader->TabAxis)[2]);
//                break;    
//         default: cerr << "Error: bad NAXIS keywords in " << File_Name << endl;
//                  cerr << "       NAXIS must between 1 and 3" << endl;
//                break;    
//     }
// 
//    readdataf(FitsHeader->file,File_Name, FitsHeader->bitpix, Data.n_elem(),
//               Data.buffer(), FitsHeader->bscale, FitsHeader->bzero);
//    if( openflag) fclose(FitsHeader->file);
// }

/************************************/
 
void fits_write_block(char *filename, fltarray & Image, int Indi, int Indj,
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
    lpixels[0] = fpixels[0] +  Image.nx() - 1;
    lpixels[1] = fpixels[1] +  Image.ny() - 1;
    lpixels[2] = fpixels[2] +  Image.nz() - 1;
    Ptr = Image.buffer();
    nelements =  Image.n_elem();

    if ( fits_write_subset_flt(fptr, group, naxis, naxes, fpixels, lpixels,  
                               Ptr,  &status))
                 FitsPrintError( status );
    if ( fits_close_file(fptr, &status) ) FitsPrintError( status );
}

/************************************/

void fits_write_block(char *filename, intarray & Image, int Indi, int Indj,
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
    lpixels[0] = fpixels[0] +  Image.nx() - 1;
    lpixels[1] = fpixels[1] +  Image.ny() - 1;
    lpixels[2] = fpixels[2] +  Image.nz() - 1;
    Ptr = Image.buffer();
    nelements =  Image.n_elem();

    if ( fits_write_subset_int(fptr, group, naxis, naxes, fpixels, lpixels,
                  Ptr,  &status))
                 FitsPrintError( status );
    if ( fits_close_file(fptr, &status) ) FitsPrintError( status );
}

/************************************/


void fits_read_block(char *filename, intarray & Image, int Indi, int Indj, Bool NoBscale)
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
    lpixels[0] = fpixels[0] +  Image.nx() - 1;
    lpixels[1] = fpixels[1] +  Image.ny() - 1;
    lpixels[2] = fpixels[2] +  Image.nz() - 1 ;
    Ptr = Image.buffer();
    nelements = Image.n_elem();
    if ( fits_read_subset_int(fptr, 0 , naxis, naxes, fpixels, lpixels, inc, nulval,
                  Ptr, &anynul, &status))
                 FitsPrintError( status );
    if ( fits_close_file(fptr, &status) ) FitsPrintError( status );
}
 
/************************************/
 
void fits_read_block(char *filename, fltarray & Image, int Indi, int Indj, Bool NoBscale)
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
    lpixels[0] = fpixels[0] +  Image.nx() - 1;
    lpixels[1] = fpixels[1] +  Image.ny() - 1;
    lpixels[2] = fpixels[2] +  Image.nz() - 1 ;
    Ptr = Image.buffer();
    nelements = Image.n_elem();

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

/************************************/

void fits_read_dblarr(char *File_Name, dblarray & Data)
{
    fitsfile *fptr; 
    int status,  nfound, anynull;
    long naxes[3], group, fpixel; // , npixels, ii;
    double nullval;
    status = 0;
    char FitsName[MAXCHAR];
    strcpy(FitsName, fitsname(File_Name));

    if (fits_open_file(&fptr,  FitsName, (int) READONLY, &status) )
    {
          printf("Error: cannot open file %s status=%d ", FitsName, status);
          exit(-1);
    }
 
    naxes[2] = 0; 
    naxes[0] = 0; 
    naxes[1] = 0;
    if (fits_read_keys_lng(fptr, (char*)"NAXIS", 1, 3, naxes, &nfound, &status) )
    {
        printf(" Error: cannot read NAXIS keyword");
        exit(-1);
    }

    switch(nfound)
    {
       case 1 : Data.alloc(naxes[0] ); break;
       case 2 : Data.alloc(naxes[0], naxes[1] ); break;
       case 3 : Data.alloc(naxes[0], naxes[1] , naxes[2] ); break;
       default : printf( " erreur nfound %d \n", nfound); break;
    }
    group  = 1;
    fpixel = 1;
    nullval= 0; 

    if(fits_read_img(fptr, TDOUBLE, fpixel, Data.n_elem(),
                     (void *) &nullval, Data.buffer(),
                     &anynull, &status) )
    {
       printf("\n error in fits_read_img %s", FitsName);
       exit(-1);
    }

    if (fits_close_file(fptr, &status) )    
    {
        printf("\n erreur fermeture %s status = %d",FitsName, status);
        exit(-1);
    }
}

/************************************/

void fits_write_dblarr(char *File_Name, dblarray & Data)
{
    fitsfile *fptr;   // pointer to the FITS file  
    long naxes[3], group, firstpixel, npixels;

    int status = 0;
    int bitpix = -64;   
    int simple = 1;   // TRUE 
    long pcount = 0;  // no group parameters 
    long gcount = 1;  // only a single image/group  
    int extend = 0   ;  // no extension 
    char FitsName[MAXCHAR];
    strcpy(FitsName, fitsname(File_Name));
 
     remove(FitsName);
     if ( fits_create_file(&fptr, FitsName, &status) )      
     {
          printf("Error: cannot open file %s %d ", FitsName, status);
          exit(-1);
     }
     naxes[0] = (long) Data.axis(1);  // int converted to long  
     naxes[1] = (long) Data.axis(2);
     naxes[2] = (long) Data.axis(3);
     // it can be a vector, an image or a cube  
     if (fits_write_grphdr(fptr, simple, bitpix, Data.naxis(), 
                           naxes, pcount, gcount, extend, &status) )
     {
        printf("Error: cannot write header %s %d", FitsName, status);
        exit(-1);
     }
     group  = 1;  
     firstpixel = 1;                  
     npixels = (long) Data.n_elem();
 
     if ( fits_write_img( fptr, TDOUBLE, firstpixel, npixels, Data.buffer(), 
             &status) )
     {
          printf("Error: cannot write %s %d ", FitsName, status);
          exit(-1);
     }
     if ( fits_close_file(fptr, &status) )                    
     {
          printf("Error: cannot close %s %d ", FitsName, status);
          exit(-1);
     }
}

/************************************/

void fits_read_intarr(char * FileName, intarray & Tab)
{
   fltarray Tab_Flt;
   fits_read_fltarr(FileName, Tab_Flt);
   Tab.alloc(Tab_Flt.nx(),Tab_Flt.ny(),Tab_Flt.nz());
   for(int i=0;i<Tab_Flt.nx();i++)
   for(int j=0;j<Tab_Flt.ny();j++)
   for(int k=0;k<Tab_Flt.nz();k++)
		Tab(i,j,k)=iround(Tab_Flt(i,j,k));
}

/************************************/

void fits_write_intarr(char * FileName, intarray & Tab)
{
 	fltarray Tab_Flt(Tab.nx(),Tab.ny(),Tab.nz());
	for(int i=0;i<Tab.nx();i++)
	for(int j=0;j<Tab.ny();j++)
	for(int k=0;k<Tab.nz();k++)
		Tab_Flt(i,j,k)=(float)(Tab(i,j,k));
	fits_write_fltarr(FileName, Tab_Flt);
}

/************************************/

void fits_read_cfarr(char *File_Name, cfarray & DataCF, fitsstruct *Header,
 int openflag)
{
    char *Name = fitsname(File_Name);
    fitsfile *fptr;                       
    int status,  nfound, anynull;
    int hdutype;
    long  naxes[3];
    long group, fpixel;
    int bitpix   =  -32;   /* -32 for 32-bit float */
    long Nelem;
    float  nullval;
    int fits_read_com( fitsfile *fptr, char *keyname,       
                      char **value );  /* O - pointer to keyword value  */
//     int fits_get_img_dim (fitsfile *fptr, int *naxis, int *status);
//     int fits_get_img_size (fitsfile *fptr, int maxdim,long *naxes, int *status);
//     int fits_get_img_type(fitsfile *fptr, int *bitpix, int *status);
//     int fits_copy_header(fitsfile *infptr, fitsfile *outfptr, int *status);
    status = 0;

    if ( fits_open_file(&fptr, Name,READONLY, &status) ) /* open the image */
    {
       printf("Error: cannot open file %s\n", Name);
       exit(status);
    }
    
    if ( fits_get_img_type(fptr, &bitpix, &status) )
    {
       printf("Error: in fits_get_img_type ...\n");
       exit(status);
    }
    // fits_get_img_param : ffgipr 
    if ( fits_get_img_dim(fptr, &nfound, &status) )
    {
        printf("Error: pb NAXIS ...\n");
        exit(status);
    }
     
    if (nfound == 0)
    {
        // printf(" no primary image look for extension \n");
        if ( fits_movabs_hdu(fptr, 2, &hdutype, &status) ) // move to 2nd HDU 
        {
            printf("Error: no image extension ...\n");
            exit(status);
        }
        if (hdutype != IMAGE_HDU)  
        {
            printf("Error: no image in this HDU\n");
            exit(-1);
        }
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
    
fltarray *Data;
    switch(nfound)
    {
		case 2 : 
			Data->alloc(naxes[0], naxes[1] ); 
			DataCF.alloc((std::complex<float>*)Data->buffer(),naxes[0]/2, naxes[1],"toto");
			break;
		case 3 : 
			Data->alloc(naxes[0], naxes[1] , naxes[2] ); 
			DataCF.alloc((std::complex<float>*)Data->buffer(),naxes[0]/2, naxes[1], naxes[2],"toto");
			break;
		default : 
			printf( "Error: dimension =  %d \n", nfound);
			exit(-1);
			break;
    }

  
   group    = 1;
   fpixel   = 1;
   nullval  = 0.;                
   
   Nelem = (long) Data->n_elem();
   if(fits_read_img(fptr, TFLOAT, fpixel, Nelem,(void *) &nullval, 
                 Data->buffer(), &anynull, &status) )
   {
      printf("Error: cannot read in file %s\n", File_Name);
      exit(status);
   }
   
	//replacing tmpnam by mkstemp
	//char FitsName[L_tmpnam]; // = "xxhdmr.fit";
	//char *Res = tmpnam(FitsName);
	//char FitsName[] = "hdmrXXXXXX";
	//int fd;   
	//fd = mkstemp(FitsName);
	
	// mkstemp crashes after 1024 calls
	// so we create the temp name ourself
	char FitsName[] = "hdmrXXXXXXXXXX";
	timeval tim;
	gettimeofday(&tim, NULL);
	long int tt=long(tim.tv_sec*1000000L+tim.tv_usec)%long(1e10L);
	char ttt[16];
	sprintf(ttt,"%10ld",tt);
	
	do
	{
		strcpy(FitsName,"hdmr");
		strcat(FitsName,ttt);
		tt+=1;
	}
	while(file_exists(FitsName));
	
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
	if(fits_close_file(fptr1, &status))
	{
		printf("Error: cannot close %s\n",FitsName);
		exit(status);
	}	
   
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
   free(Name);
}

/************************************/

void fits_read_cfarr(char *File_Name, cfarray & Data)
{
   fitsstruct FitsHeader;
   fits_read_cfarr(File_Name, Data, &FitsHeader, 1);
}


/************************************/

void fits_read_cfarr(char *File_Name, cfarray & Data, fitsstruct *FitsHeader)
{
   fits_read_cfarr(File_Name, Data, FitsHeader, 1);
}

