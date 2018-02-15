/*******************************************************************************
**
**    UNIT
**
**    Version: 3.3
**
**    Author: Jean-Luc Starck
**
**    Date:  96/06/13 
**    
**    File:  IM_IO.h
**
*******************************************************************************
**
**    DESCRIPTION  FITS Include
**    ----------- 
**                 
******************************************************************************/

#ifndef _IM_IO_H_
#define _IM_IO_H_


#include"Array.h"
#include"IM_IOTools.h"
#include"writefits3d.h"
#include<unistd.h>

typedef unsigned char byte;
typedef unsigned long u_long;

#include"IM_Lut.h"

#define IO_FITS_OK 1
#ifdef MIDAS
#define IO_MIDAS_OK 1
#else 
#define IO_MIDAS_OK 0
#endif

#ifdef NO_DISP_IO
#define IO_DISP_OK 0
#else 
#define IO_DISP_OK 1
#endif

#define IO_JPEG_OK 0
#define IO_GIF_OK 0
#define IO_PGM_OK 0
#define IO_TIFF_OK 0

// Available format
enum type_format {F_UNKNOWN, F_DISP, F_MIDAS, F_FITS, F_GIF, F_PGM, F_JPEG};

enum type_data {T_BYTE, T_SHORT, T_INT, T_FLOAT, T_DOUBLE,
		     T_COMPLEX_F, T_COMPLEX_D, UNKNOWN};

void io_strcpy_prefix(char *Prefix, char *Filename);
// copy  Filename into Prefix without the suffix
	     
void io_set_format(type_format Format=F_DISP);
// set the format image to FITS 

type_data which_data_type();
// return the type of the last read data set

void io_set_format(const char *Format);
// determine and set  the format from a file name

type_format io_which_format(const char *File_Name);
// determine the format associated to the file name
// and return it. If no format is detected, return
// F_UNKNOWN

type_format io_detect_format(const char *File_Name);
// determine the format associated to the file name
// and return it. If no format is detected, return
// the default format

type_format which_format();
// return the actual format

Bool std_inout(char *Filename);
// return true if standard input-output are specified.
// std in-out are specified if the filename start with a "-" caracter.

inline const char *StringFormat (type_format type)
{
    switch (type)
    {
        case F_DISP: 
              return ("DISP format");break;
        case F_MIDAS: 
              return ("MIDAS format");break;
        case F_FITS: 
              return ("FITS format");break;
        case F_GIF: 
              return ("GIF format");break;
        case F_PGM: 
              return ("PGM format");break;
	case F_JPEG:
	      return ("JPEG format");break;
        default: 
              return ("Unknown format");break;
    }
}

#define DEFAULT_FORMAT_IMAGE  F_FITS


 /*--------------------------- FITS BitPix coding ----------------------------*/

#define         BP_BYTE         8
#define         BP_SHORT        16
#define         BP_INT          32
#define         BP_FLOAT        (-32)
#define         BP_DOUBLE       (-64)


/*----------------------------- Fits image parameters ----------------------------*/

// typedef struct
class fitsstruct {
 void fitsinit()
 {
    file=NULL;	
    fitsheadsize= 0;
    bitpix = 0;
    bytepix = 0;
    width = 0;
    height = 0;
    npix = 0;
    bscale = 1.;
    bzero = 0.;
    crpixx = 0.;
    crpixy = 0.;
    crvalx = 0.;
    crvaly = 0.;
    cdeltx= 0.;
    cdelty = 0.;
    ngamma=0.;
    pixscale=0.;
    nlevels = 0;
    pixmin  = 0.;
    pixmax  = 0.;
    epoch = 0.;
    crotax = 0.;
    crotay = 0.;
    fitshead = NULL;
    origin = NULL;
    strcpy(ctypex, "");
    strcpy(ctypey, "");
    strcpy(rident, "");

    /* HISTORY & COMMENT fields */
    history = NULL;
    hist_size = 0;
    comment = NULL;
    com_size = 0;

    naxis=0;
    for (int i=0; i < MAX_NBR_AXIS;i++)
    {
       TabAxis[i]=0;
       TabStep[i]=0.;
       TabRef[i]=0.;
       TabValRef[i]=0.;
    }
    filename = NULL;
    origin = NULL;
    fitshead = NULL;
   }
 public:
  
  fitsstruct ()  
  { 
     fitsinit();
  }
  void hd_fltarray(fltarray &Mat, char *History=NULL)
  {
     char *creafitsheader();

     fitshead = creafitsheader();
     bitpix = -32;
     width =  Mat.nx();
     height = Mat.ny();
     naxis = Mat.naxis();
     npix = Mat.n_elem();
     for (int i=0; i < naxis; i++) TabAxis[i] = Mat.axis(i+1);
     if (History != NULL) origin = History;
     else origin=NULL;
  }  
  void hd_init(int Nx, int Ny=0, int Nz=0, char *History=NULL)
  {
     fitsheadsize= 2880;
     char *creafitsheader();
     fitshead = creafitsheader();
     fitsinit();
     bitpix = -32;
     TabAxis[0] = width = Nx;
     TabAxis[1] = height = Ny;
     TabAxis[2] = Nz;
     npix = width;
     naxis = 1;
     if (Ny > 0) 
     {
        naxis++;
        npix *= Ny;
     }
     if (Nz > 0) 
     {
        naxis++;
        npix *= Nz;
     }
     if (History != NULL) origin = History;
     else origin = NULL;
  }  
  ~fitsstruct()  
  { 
      if (filename != NULL) free (filename);
//      if (origin != NULL) delete origin;
//      if (fitshead != NULL) delete[] ((char *) fitshead);
      if (fitshead != NULL) free ((char *) fitshead);
      if (history != NULL)  free (history);
      if (comment != NULL)  free (comment);
     fitsinit();
  }
  
  char		*filename;		/* pointer to the image filename */
  char          *origin;                /* pointer to the origin */
  char		ident[512];		/* field identifier (read from FITS)*/
  char		rident[512];	        /* field identifier (relative) */
  FILE		*file;			/* pointer the image file structure */
  char		*fitshead;		/* pointer to the FITS header */
  int		fitsheadsize;		/* FITS header size */
/* ---- main image parameters */
  int		bitpix, bytepix;	/* nb of bits and bytes per pixel */
  int		width, height;		/* x,y size of the field */
  int		npix;			/* total number of pixels */
  double	bscale, bzero;		/* FITS scale and offset */
  double	ngamma;			/* normalized photo gamma */
  int		nlevels;		/* nb of quantification levels */
  float		pixmin, pixmax;		/* min and max values in frame */
/* ---- basic astrometric parameters */
  double	epoch;			/* epoch for coordinates */
  double	pixscale;		/* pixel size in arcsec.pix-1 */
					/* */
/* ---- astrometric parameters */
  double	crpixx,crpixy;		/* FITS CRPIXn */
  double	crvalx,crvaly;		/* FITS CRVALn */
  double	cdeltx,cdelty;		/* FITS CDELTn */
  double	crotax,crotay;		/* FITS CROTAn */
  char          ctypex[256];            /* FITS CTYPE1 */
  char          ctypey[256];            /* FITS CTYPE2 */
  char          CoordType[256];         
  
/* ---- HISTORY & COMMENT parameters --- */
	 char *history;
	 int hist_size;
	 char *comment;
	 int com_size;

/* ---- for non image use */
  int naxis;
  int TabAxis[MAX_NBR_AXIS];
  double TabStep[MAX_NBR_AXIS];
  double TabRef[MAX_NBR_AXIS];
  double TabValRef[MAX_NBR_AXIS];
  };
 //	} fitsstruct;


// Read the celestial coordinates from a fits file
int cfitstio_read_celestial_coord(char *File_Name, fitsstruct *FitsHeader);

// calculate pixel coordinates from celestial position
//   X,Y coordinates are defined from origin (0,0) (and not (1,1) as in fits
//   format)
void adxy(fitsstruct *FitsHeader, double Ra, double Dec, double *X, double *Y);

// convert Ra, Dec on decimal degrees to HH,MN,SEC and DEG,MN,SEC
void radec(double Ra, double Dec, int &ihr, int &imin, double &xsec,
                                  int &ideg,int &imn, double &xsc);
                                  
// calculate celestial position from pixels coordinates
// output X,Y are given in a referential starting at origin (0,0)
void xyad(fitsstruct *FitsHeader, double X, double Y, double *Ra, double *Dec);

void io_read_ima_int  (char *File_Name, Iint &, fitsstruct *hd = NULL, Bool NoBscale = False);
void io_write_ima_int (char *, Iint &, fitsstruct *hd = NULL);
void io_read_ima_float  (char *, Ifloat &, fitsstruct *hd = NULL);
void io_write_ima_float (char *, Ifloat &, fitsstruct *hd = NULL, Bool NormData=False);
void io_read_ima_complex_f  (char *, Icomplex_f &);
void io_write_ima_complex_f (char *, Icomplex_f &);

char *fitsname(char * NameStep);
FILE *fits_file_des_in(char *fname);
FILE *fits_file_des_out(char *fname);
void fits_read_header(char *File_Name, fitsstruct *Header);
void fits_read_fltarr(char *File_Name, fltarray &Mat);
void fits_read_fltarr(char *File_Name, fltarray &Mat, fitsstruct * FitsHeader);
void fits_read_fltarr(char *File_Name, fltarray &Mat, fitsstruct * FitsHeader,
                      int openflag);



void fits_read_cfarr(char *File_Name, cfarray &Mat);
void fits_read_cfarr(char *File_Name, cfarray &Mat, fitsstruct * FitsHeader);
void fits_read_cfarr(char *File_Name, cfarray & DataCF, fitsstruct *Header,
					  int openflag);



void fits_read_block(char *filename, Ifloat & Image, int Indi, int Indj, Bool NoBscale=False);
void fits_read_block(char *filename, Iint & Image, int Indi, int Indj, Bool NoBscale=False);
void fits_write_block(char *filename, Ifloat & Image, int Indi, int Indj,
                      Bool NoBscale=False);
void fits_write_block(char *filename, Iint & Image, int Indi, int Indj,
                      Bool NoBscale=False);
void fits_read_float(char *File_Name,Ifloat & Image,fitsstruct *Header);
void fits_read_int(char *File_Name,Iint & Image, 
                     fitsstruct *Header, Bool NoBscale=False);

void fits_write_float(char *File_Name,Ifloat & Image,fitsstruct *Header);
void fits_write_int(char *File_Name, Iint & Image,fitsstruct *Header);
void fits_write_header(char *File_Name, fitsstruct *Header);

void fits_write_fltarr(char *File_Name, fltarray &Mat);
void fits_write_fltarr(char *File_Name, fltarray &Mat, fitsstruct *FitsHeader);
void fits_read_dblarr(char *File_Name, dblarray & Data);
void fits_write_dblarr(char *File_Name, dblarray & Data);
void fits_write_intarr(char * FileName, intarray & Tab);
void fits_read_intarr(char * FileName, intarray & Tab);

void makehistory(char *mystring, char *myproc, char *myalgo, char *myargs);
int fitsaddhist_com(fitsstruct *pfitsbuf, char *comment, char *type_com);
int fitsread(char *fitsbuf, char *keyword, void *ptr,
             h_type type, type_data t_type);
char *readfitshead(FILE *file, char *filename, int *nblock);
char *creafitsheader();
void initfield(fitsstruct *Header); /* initialize the structure FITSSTRUCT */
void init_fits_struct(fitsstruct *Ptr, int Nl, int Nc);
void FitsPrintError( int status);

void fits_read_block(char *filename, fltarray & Image, int Indi, int Indj, Bool NoBscale=False);
void fits_read_block(char *filename, intarray & Image, int Indi, int Indj, Bool NoBscale=False);
void fits_write_block(char *filename, intarray & Image, int Indi, int Indj,
                      Bool NoBscale=False);
void fits_write_block(char *filename, fltarray & Image, int Indi, int Indj,
                      Bool NoBscale=False);

class io_disp;
class io_midas;
class io_fits;
class IOInfoData;

class io_image {
        friend class io_disp;
        friend class io_midas;
        friend class io_fits;
          union {
           int *Data_int;
           short int *Data_short;
           unsigned char *Data_byte;
           float *Data_float;
           double *Data_double;
           complex_f *Data_complex_f;
           complex_d *Data_complex_d;
                 };
        public:
           int Nl_io;
           int Nc_io;
           type_data Type;
           io_image() {};
           virtual void read_ima_int  (char *, Iint &){};
           virtual void write_ima_int (char *, Iint &){};
           virtual void read_ima_float  (char *, Ifloat &){};
           virtual void write_ima_float (char *, Ifloat &){};
           virtual void read_ima_complex_f  (char *, Icomplex_f &){};
           virtual void write_ima_complex_f (char *, Icomplex_f &){};
           virtual ~io_image(){};
             };

class io_disp: public virtual io_image {
          void read_type(char *File_Name);   
          void read_size (char *);
          void read_ima (char *);
          void write_ima (char *);
         public:
          io_disp() {Nl_io = 0; Nc_io = 0; Type = UNKNOWN;}
          ~io_disp(){};
          void read_ima_int  (char *, Iint &);
          void write_ima_int (char *, Iint &);
          void read_ima_float  (char *, Ifloat &);
          void write_ima_float (char *,  Ifloat &);
          void read_ima_complex_f  (char *, Icomplex_f &);
          void write_ima_complex_f (char *, Icomplex_f &);
	  void read_info_ima(char *File_Name, int &Nl, int &Nc, 
                            type_data & TypeDat);
	 void read_block_ima(char *File_Name,Ifloat & Image,int Indi, int Indj);
	 void read_block_ima(char *File_Name,Iint & Image,int Indi, int Indj);
	 void write_block_ima(char *File_Name,Ifloat & Image,int Indi, int Indj);
	 void write_block_ima(char *File_Name,Iint & Image,int Indi, int Indj);
	 
         void create_file (char *File_Step_Name); 
              };



/***** Definition of the Midas image structure  */

#define SIZE_PIC_CUNIT 63      /* Max of the comment size on the axes */
#define SIZE_PIC_INDENT 71     /* Max of the comment size on the image */
#define SIZE_PIC_TAB_NAXIS 1   /* Dimensions array size */
#define SIZE_PIC_TAB_NPIX 2    /* Array size for the dimensions on each axis */
#define SIZE_PIC_START 2       /* Origine coordinates array size */
#define SIZE_PIC_STEP 2        /* Step array size  */
#define SIZE_PIC_TAB_CUTS 4    /* Cuts array size  */

typedef struct  {
        float *Rbuf;
           /* Buffer where the data are ranged */
        float Cuts[SIZE_PIC_TAB_CUTS];
           /* Cuts [0] ... Cuts [1] for the visualisation 
              Cuts [2] ... Cuts [3] = min ... max of the image
           */
        char Cunit[SIZE_PIC_CUNIT];
           /* Axis comments  */
        char Ident[SIZE_PIC_INDENT];
           /* Image comments */
        int Naxis [SIZE_PIC_TAB_NAXIS];
           /* Dimension Number */
        int Npix [SIZE_PIC_TAB_NPIX];
           /* Npix [0] = pixels number on x
              Npix [1] = pixels number on y  */
        double Start [SIZE_PIC_START];
           /* Origine coordinates x = Start [0], y = Start */
        double Step [SIZE_PIC_STEP];
           /* Step on the axis */
                 } midas_pic_des;
#if IO_MIDAS_OK

class io_midas: public virtual io_image {
         public:
          io_midas();
          ~io_midas();
          void read_ima_int  (char *, Iint &);
          void write_ima_int (char *,const Iint &);
          void read_ima_float  (char *, Ifloat &);
          void write_ima_float (char *,const Ifloat &);
          void read_ima_complex_f  (char *, Icomplex_f &);
          void write_ima_complex_f (char *,const Icomplex_f &);
	  void read_info (char *File_Name, IOInfoData &InfoDat);
             };
#endif

#include "IM_IOCol.h"

void readgif(char *fname, Ifloat &Data);
void writegif(char *fname, Ifloat & Data, Bool NormData=False);
void readgif(char *fname, Iint & Data);
void readgif(char *fname, PICINFO &pinfo);
void writegif(char *fname, Iint & Data, Bool NormData=False);
int readgif(FILE *fp, PICINFO *pinfo);
int writegif(FILE *fp, PICINFO pinfo);
void writegif(char *fname, PICINFO & pinfo);

void writejpeg(char *fname, Ifloat & Data, Bool NormData=False);
void writejpeg(char *fname, Iint & Data, Bool NormData=False);

// gif block access routines
void readhdgif(char *filename, int &Nl, int &Nc);
void readhdgif(char *File_Name, int &Nl, int &Nc, PICINFO * & GifIma);
void read_block_gif (Ifloat & Image, int Indi, int Indj, PICINFO *GifIma);
void read_block_gif (Ifloat & Image, int Indi, int Indj);
void read_block_gif (Iint & Image, int Indi, int Indj);
void read_block_gif (Iint & Image, int Indi, int Indj, PICINFO *GifIma);
void gif_block_set_ptr_in(PICINFO *GifIma);

void writehdgif (const char *filename, int Nl, int Nc);
void writehdgif (const char *filename, int Nl, int Nc, PICINFO *GifIma);
void write_block_gif (Ifloat & Image, int Indi, int Indj);
void write_block_gif (Iint & Image, int Indi, int Indj);
void closegif(char *FileName); // write the byte array to the disk
void closegif(char *File_Name, PICINFO *GifIma);
void write_block_gif (Iint & Image, int Indi, int Indj, PICINFO *GifIma);
void write_block_gif (Ifloat & Image, int Indi, int Indj, PICINFO *GifIma);
void gif_block_set_ptr_out(PICINFO *GifIma);

// PGM format
char *pgmname(const char * NameStep);
FILE *pgm_file_des_in(char *fname);
FILE *pgm_file_des_out(char *fname);
void readpgm (char *filename, int &Nl, int &Nc, unsigned char *  &buffer);
void readpgm(char *fname, Iint &Data);
void readpgm(char *fname, Ifloat &Data);
void writepgm (char *filename, int Nl, int Nc, unsigned char *buffer);
void writepgm(char *fname, Iint &Data, Bool NormData=False);
void writepgm(char *fname, Ifloat &Data, Bool NormData=False);
// block access
void read_block_pgm (const char *filename, Ifloat & Image, int Indi, int Indj);
void read_block_pgm (const char *filename, Iint & Image, int Indi, int Indj);
void readhdpgm(const char *filename, int &Nl, int &Nc);
void writehdpgm(char *filename, int Nl, int Nc);
void write_block_pgm (char *filename, int Nl, int Nc,
                      Ifloat & Image, int Indi, int Indj);
void write_block_pgm (char *filename, int Nl, int Nc,
                      Iint & Image, int Indi, int Indj);

#endif
