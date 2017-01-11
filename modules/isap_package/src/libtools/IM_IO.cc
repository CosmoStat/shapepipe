/*******************************************************************************
**
**    UNIT
**
**    Version: 3.2
**
**    Author: Jean-Luc Starck
**
**    Date:  96/05/07 
**    
**    File:  IM_IO.cc
**
*******************************************************************************
**
**    DESCRIPTION  input output routines
**    ----------- 
**                 
**    PARAMETRES    
**    ----------    
** 
**    RESULTS      
**    -------  
**
**
** 
*******************************************************************************
**
** void io_set_format(type_format Format)
**
** set Format_Imag to Format
**
*******************************************************************************
**
** void io_set_format(char *Format)
**
** set Format_Imag to Format
**
*******************************************************************************
**
** type_format io_detect_format(char *File_Name) 
**
** seach format from name image
**
*******************************************************************************
**
** void io_read_ima_int  (char *File_Name, Iint &)
** 
** read int image from file
**
*******************************************************************************
**
** void io_write_ima_int (char *, Iint &)
**
** write int image to file
**
*******************************************************************************
**
** void io_read_ima_float  (char *, Ifloat &)
**
** read float image from file
**
*******************************************************************************
**
** void io_write_ima_float (char *, Ifloat &)
**
** write float image to file
** 
*******************************************************************************
**
** void io_read_ima_complex_f  (char *, Icomplex_f &)
**
** read complex float image from file
**
*******************************************************************************
**
** void io_write_ima_complex_f (char *, Icomplex_f &)
**
** write complex float image to file
** 
******************************************************************************/

// static char sccsid[] = "@(#)IM_IO.cc 3.2 96/05/07 CEA 1994 @(#)";
 
#include"GlobalInc.h"
#include "IM_Obj.h"
#include "IM_IO.h"

#if IO_DISP_OK
io_disp IO_DISP;
#endif

#if IO_MIDAS_OK
io_midas IO_MIDAS;
#endif

type_format Format_Imag = F_UNKNOWN;
type_data TypeInputData = UNKNOWN;
C_IO_RGB IO_RGB;
softinfo Soft;

#if IO_JPEG_OK
#include "IM_JPEG.h"
#endif

#ifdef USELM
extern DemoLic MRDEMO;
#endif

/**************************************************************************/

type_data which_data_type()
{
   return TypeInputData;
}

/**************************************************************************/

type_format which_format()
{
  return Format_Imag;
}


/**************************************************************************/

void io_strcpy_prefix(char *Prefix, char *Filename)
{
   char *ptr;
   strcpy(Prefix, Filename);
   if (strrchr(Filename, '.') != NULL)
   {
       if (Format_Imag != F_DISP) 
       {
           ptr = strrchr(Prefix, '.');
           *ptr = '\0';
       }
       else
       {
           ptr = strchr(Prefix, '_');
           *ptr = '\0';
       }
   }
}

/**************************************************************************/

Bool std_inout(char *Filename)
{
   char Pre[MAXCHAR];
   Bool Stdinout = False;
   
   io_strcpy_prefix(Pre, Filename);
   // fprintf(stderr, "pre = %s \n", Pre);
   if ((strlen(Pre) == 1) && (Pre[0] == '-')) Stdinout = True;
   else if ((strlen(Pre) == 2) && (Pre[0] == '-') && (Pre[1] == 'b')) Stdinout = True;

  // if ((strlen(Pre) == 2) && (Pre[0] == 'i') && (Pre[1] == 'n'))   Stdinout = True;
  // else if ((strlen(Pre) == 3) && (strstr(Pre, "inb") != NULL))  Stdinout = True;
  // else if ((strlen(Pre) == 3) &&(strstr(Pre, "out") != NULL))  Stdinout = True;
   return Stdinout;
}


/**************************************************************************/

void io_set_format(type_format Format) 
{
    Format_Imag = Format;
};

/**************************************************************************/

void io_set_format(const char *Format) 
{
    char Name_Up[MAXCHAR];
    char Name[MAXCHAR];
    const char * Disp = "dis";
    const char * Midas = "mid";
    const char * Fits = "fit";
    const char * Gif = "gif";
    const char *Pgm = "pgm";
    const char *Jpg = "jpg";
    unsigned int i;

    strcpy(Name_Up, Format);
    Format_Imag = F_UNKNOWN;
    strcpy(Name, Name_Up);
    for (i=0; i < strlen(Name); i++) Name[i] = tolower(Name[i]);
    if (strstr(Name, Disp) != NULL) Format_Imag = F_DISP;
    else if (strstr(Name, Midas) != NULL) Format_Imag = F_MIDAS;
    else if (strstr(Name, Gif) != NULL) Format_Imag = F_GIF;
    else if (strstr(Name, Pgm) != NULL) Format_Imag = F_PGM;
    else if (strstr(Name, Fits) != NULL) Format_Imag = F_FITS;
    else if (strstr(Name, Jpg) != NULL) Format_Imag = F_JPEG;
};

/**************************************************************************/

type_format io_which_format(const char *Format) 
{
    char Name_Up[MAXCHAR];
    char Name[MAXCHAR];
    const char * Disp = ".d";
    const char * Midas = ".bdf";
    const char * Fits = ".fit";
    const char * Gif = ".gif";
    const char *Pgm = "pgm";
    const char *Jpg = "jpg";
    const char *mr = "mr";
    const char *rad = "rad";
    const char *rid = "rid";
    const char *cur = "cur";
    const char *bet = "bet";
    const char *fts = "fts";
    const char *fit = "fit";
    const char *FIT = "FIT";
    unsigned int i;

    strcpy(Name_Up, Format);
    Format_Imag = F_UNKNOWN;
    strcpy(Name, Name_Up);
    for (i=0; i < strlen(Name); i++) Name[i] = tolower(Name[i]);
    if (strstr(Name, Disp) != NULL) Format_Imag = F_DISP;
    else if (strstr(Name, Midas) != NULL) Format_Imag = F_MIDAS;
    else if (strstr(Name, Gif) != NULL) Format_Imag = F_GIF;
    else if (strstr(Name, Pgm) != NULL) Format_Imag = F_PGM;
    else if (strstr(Name, Fits) != NULL) Format_Imag = F_FITS;
    else if (strstr(Name, Jpg) != NULL) Format_Imag = F_JPEG;
    else if (strstr(Name, mr) != NULL) Format_Imag = F_FITS;
    else if (strstr(Name, rad) != NULL) Format_Imag = F_FITS;
    else if (strstr(Name, rid) != NULL) Format_Imag = F_FITS;
    else if (strstr(Name, cur) != NULL) Format_Imag = F_FITS;
    else if (strstr(Name, bet) != NULL) Format_Imag = F_FITS;
    else if (strstr(Name, fts) != NULL) Format_Imag = F_FITS;
    else if (strstr(Name, fit) != NULL) Format_Imag = F_FITS;
    else if (strstr(Name, FIT) != NULL) Format_Imag = F_FITS;
    return Format_Imag;
};
    
/**************************************************************************/

type_format io_detect_format(const char *Format) 
{
    Format_Imag = F_UNKNOWN;
    Format_Imag = io_which_format(Format);
    if (Format_Imag == F_UNKNOWN) Format_Imag = DEFAULT_FORMAT_IMAGE;
    return Format_Imag;
};
    
/**************************************************************************/

void io_read_ima_int  (char *File_Name, Iint &Image, fitsstruct *Header, Bool NoBscale)
{
    fitsstruct HD,*PtrHD;

    if (Format_Imag == F_UNKNOWN)
                     Format_Imag = io_detect_format (File_Name);
    switch (Format_Imag)
    {
        case F_DISP: 
#if IO_DISP_OK
              IO_DISP.read_ima_int (File_Name, Image);
              TypeInputData = IO_DISP.Type;
#else
              fprintf (stderr, "Error: DISP is not active\n");
              exit (-1); 
#endif
              break;
        case F_MIDAS:
#if IO_MIDAS_OK
              IO_MIDAS.read_ima_int (File_Name, Image);
              TypeInputData = T_FLOAT;
#else
              fprintf (stderr, "Error: MIDAS is not active\n");
              exit (-1); 
#endif
              break;
        case F_FITS:
#if IO_FITS_OK
              if (Header != NULL) PtrHD = Header;
              else PtrHD = &HD;              
              fits_read_int(File_Name, Image, PtrHD, NoBscale);
              switch(PtrHD->bitpix)
              {
                 case  BP_BYTE: TypeInputData = T_BYTE; break;
                 case  BP_SHORT: TypeInputData = T_SHORT; break;
                 case  BP_INT: TypeInputData = T_INT; break;
                 case  BP_FLOAT: TypeInputData = T_FLOAT; break;
                 case  BP_DOUBLE: TypeInputData = T_DOUBLE; break;
                 default: break;
              }
#else
              fprintf (stderr, "Error: FITS is not active\n");
              exit (-1); 
#endif
              break;
        case F_GIF:
#if IO_GIF_OK
               readgif (File_Name, Image);
               TypeInputData = T_BYTE;
#else
              fprintf (stderr, "Error:GIF is not active\n");
              exit (-1); 
#endif
               break;
         case F_PGM: 
#if IO_PGM_OK
               readpgm (File_Name, Image);
               TypeInputData = T_BYTE;
#else
              fprintf (stderr, "Error:PGM is not active\n");
              exit (-1); 
#endif
               break;      
         case F_JPEG:
#if IO_JPEG_OK
               io_read2d_jpeg(File_Name, Image);
               TypeInputData = T_BYTE;
#else
              fprintf (stderr, "Error:JPEG is not active\n");
              exit (-1); 
#endif
               break;
        default:
              fprintf (stderr, "Error: bad image format. cannot read ...\n");
              exit (-1);
              break;
    }
    
#ifdef USELM
    MRDEMO.test(Image.nc(), Image.nl());
#endif
}

/**************************************************************************/

void io_write_ima_int  (char *File_Name, Iint &Image, fitsstruct *Header)
{
     fitsstruct HD;

    if (Format_Imag == F_UNKNOWN)
                     Format_Imag = io_detect_format (File_Name);

    switch (Format_Imag)
    {
        case F_DISP: 
#if IO_DISP_OK
              IO_DISP.write_ima_int (File_Name, Image);
#else
              fprintf (stderr, "Error: DISP is not active\n");
              exit (-1); 
#endif
              break;
        case F_MIDAS:
#if IO_MIDAS_OK
              IO_MIDAS.write_ima_int (File_Name, Image);
#else
              fprintf (stderr, "Error: MIDAS is not active\n");
              exit (-1); 
#endif
              break;
        case F_FITS:
#if IO_FITS_OK
              if (Header != NULL) 
              {
                 Header->naxis = 2;
                 Header->width =  Image.nc();
                 Header->height = Image.nl();
                 (Header->TabAxis)[0] = Image.nc();
                 (Header->TabAxis)[1] = Image.nl();
                 Header->npix =  Image.nc()*Image.nl();
                 fits_write_int(File_Name, Image, Header);
              }
              else 
              {
                  initfield(&HD);
                  HD.naxis = 2;
                  HD.bitpix = 32;
                  HD.width =  Image.nc();
                  HD.height = Image.nl();
                  HD.TabAxis[0] = Image.nc();
                  HD.TabAxis[1] = Image.nl();
                  HD.npix = HD.TabAxis[0]*HD.TabAxis[1];
                  fits_write_int(File_Name, Image, &HD);
              }
#else
              fprintf (stderr, "Error: FITS is not active\n");
              exit (-1); 
#endif
              break;
        case F_GIF:  
#if IO_GIF_OK
             writegif(File_Name, Image); 
#else
              fprintf (stderr, "Error: GIF is not active\n");
              exit (-1); 
#endif
             break;
        case F_PGM:
#if IO_PGM_OK
              writepgm(File_Name, Image); 
#else
              fprintf (stderr, "Error: PGM is not active\n");
              exit (-1); 
#endif
             break;
        case F_JPEG: 
#if IO_JPEG_OK
              io_write2d_jpeg(File_Name, Image); 
#else
              fprintf (stderr, "Error: JPEG is not active\n");
              exit (-1); 
#endif
             break;
          default:
              fprintf (stderr, "Error: bad image format. cannot write ...\n");
              exit (-1);
              break;
    }
}

/**************************************************************************/

void io_read_ima_float (char *File_Name, Ifloat &Image, fitsstruct *Header)
{
    fitsstruct HD,*PtrHD;
  
    if (Format_Imag == F_UNKNOWN)
                     Format_Imag = io_detect_format (File_Name);
    switch (Format_Imag)
    {
        case F_DISP: 
#if IO_DISP_OK
              IO_DISP.read_ima_float (File_Name, Image);              
              TypeInputData = IO_DISP.Type;
#else
              fprintf (stderr, "Error: DISP is not active\n");
              exit (-1); 
#endif
              break;
        case F_MIDAS:
#if IO_MIDAS_OK
              IO_MIDAS.read_ima_float (File_Name, Image);
              TypeInputData = T_FLOAT;
#else
              fprintf (stderr, "Error: MIDAS is not active\n");
              exit (-1); 
#endif
              break;
        case F_FITS:
#if IO_FITS_OK
              if (Header != NULL) PtrHD = Header;
              else PtrHD = &HD;              
              fits_read_float(File_Name, Image, PtrHD);
              switch(PtrHD->bitpix)
              {
                 case BP_BYTE: TypeInputData = T_BYTE; break;
                 case BP_SHORT: TypeInputData = T_SHORT; break;
                 case BP_INT: TypeInputData = T_INT; break;
                 case BP_FLOAT: TypeInputData = T_FLOAT; break;
                 case BP_DOUBLE: TypeInputData = T_DOUBLE; break;
                 default: break;
              }
#else
              fprintf (stderr, "Error: FITS is not active\n");
              exit (-1); 
#endif
              break;
         case F_GIF: 
#if IO_GIF_OK
               readgif (File_Name, Image);
#else
              fprintf (stderr, "Error:GIF is not active\n");
              exit (-1); 
#endif
               break;
         case F_PGM: 
#if IO_PGM_OK
               readpgm (File_Name, Image);
#else
              fprintf (stderr, "Error:PGM is not active\n");
              exit (-1); 
#endif
               break;      
         case F_JPEG:
#if IO_JPEG_OK
              io_read2d_jpeg(File_Name, Image);
#else
              fprintf (stderr, "Error:JPEG is not active\n");
              exit (-1); 
#endif
               break;
        default:
              fprintf (stderr, "Error: bad image format. cannot read ...\n");
              exit (-1);
              break;
    }
#ifdef USELM
    MRDEMO.test(Image.nc(), Image.nl());
#endif
}

/**************************************************************************/

void io_write_ima_float  (char *File_Name, Ifloat &Image, fitsstruct *Header, Bool NormData)
{
    fitsstruct HD;

    if (Format_Imag == F_UNKNOWN)
                     Format_Imag = io_detect_format (File_Name);

    switch (Format_Imag)
    {
        case F_DISP: 
#if IO_DISP_OK
              IO_DISP.write_ima_float (File_Name, Image);
#else
              fprintf (stderr, "Error: DISP is not active\n");
              exit (-1); 
#endif
              break;
        case F_MIDAS:
#if IO_MIDAS_OK
              IO_MIDAS.write_ima_float (File_Name, Image);
#else
              fprintf (stderr, "Error: MIDAS is not active\n");
              exit (-1); 
#endif
              break;
        case F_FITS:
#if IO_FITS_OK
              if (Header != NULL) 
              {
                 Header->naxis = 2;
                 Header->width =  Image.nc();
                 Header->height = Image.nl();
		 Header->filename = strdup(File_Name);
                 (Header->TabAxis)[0] = Image.nc();
                 (Header->TabAxis)[1] = Image.nl();
                 Header->npix =  Image.nc()*Image.nl();
                 fits_write_float(File_Name, Image, Header);
              }
              else 
              {
                  initfield(&HD);
                  HD.naxis = 2;
                  HD.bitpix = -32;
                  HD.width =  Image.nc();
                  HD.height = Image.nl();
                  HD.filename = strdup(File_Name);
                  HD.TabAxis[0] = Image.nc();
                  HD.TabAxis[1] = Image.nl();
                  HD.npix = HD.TabAxis[0]*HD.TabAxis[1];
                  fits_write_float(File_Name, Image, &HD);
              }
#else
              fprintf (stderr, "Error: FITS is not active\n");
              exit (-1); 
#endif
              break;        
        case F_GIF: 
#if IO_GIF_OK
              writegif(File_Name, Image, NormData); 
#else
              fprintf (stderr, "Error: GIF is not active\n");
              exit (-1); 
#endif
              break;
        case F_PGM: 
#if IO_PGM_OK
              writepgm(File_Name, Image, NormData);
#else
              fprintf (stderr, "Error: PGM is not active\n");
              exit (-1); 
#endif
              break;
        case F_JPEG:
#if IO_JPEG_OK
              io_write2d_jpeg(File_Name, Image, NormData); 
#else
              fprintf (stderr, "Error: JPEG is not active\n");
              exit (-1); 
#endif
              break;
        default:
              fprintf (stderr, "Error: bad image format. cannot write ...\n");
              exit (-1);
              break;
    }
}

/**************************************************************************/

static void  read_datacomplex_f (char *File_Step_Name, Icomplex_f & Image)
{
    int i,j;
    char File_Name[100];
    Ifloat Buff_re, Buff_im;
    int Nline, Ncol;
    char Ext[10];

    if (Format_Imag == F_FITS) strcpy(Ext, ".fits");
    else if (Format_Imag == F_GIF) strcpy(Ext, ".gif");
    else if (Format_Imag == F_PGM) strcpy(Ext, ".pgm");
    else if (Format_Imag == F_JPEG) strcpy(Ext, ".jpg");

    /* read the real part  of the image  */
    strcpy (File_Name, File_Step_Name);
    strcat (File_Name, "_re");
    strcat (File_Name, Ext);
    io_read_ima_float (File_Name, Buff_re);

    /* read the imaginary part of the image  */
    strcpy (File_Name, File_Step_Name);
    strcat (File_Name, "_im");
    strcat (File_Name, Ext);
    io_read_ima_float (File_Name, Buff_im);

    Nline = Buff_re.nl();
    Ncol = Buff_re.nc();
    Image.alloc(Nline, Ncol, File_Step_Name);

    for ( i=0; i < Nline; i++)
    for ( j=0; j < Ncol; j++) 
          Image(i,j) = complex_f(Buff_re(i,j), Buff_im(i,j));


}

/****************************************************************/

void io_read_ima_complex_f  (char *File_Name, Icomplex_f &Image)
{
    if (Format_Imag == F_UNKNOWN)
                     Format_Imag = io_detect_format (File_Name);
    switch (Format_Imag)
    {
        case F_DISP: 
#if IO_DISP_OK
              IO_DISP.read_ima_complex_f (File_Name, Image);
              TypeInputData = IO_DISP.Type;
#else
              fprintf (stderr, "Error: DISP is not active\n");
              exit (-1); 
#endif
              break;
        case F_MIDAS:
#if IO_MIDAS_OK
              IO_MIDAS.read_ima_complex_f (File_Name, Image);
              TypeInputData = T_FLOAT;
#else
              fprintf (stderr, "Error: MIDAS is not active\n");
              exit (-1); 
#endif
              break;
        case F_FITS:
        case F_GIF:
        case F_PGM:
        case F_JPEG:
              read_datacomplex_f(File_Name, Image);
              TypeInputData = T_FLOAT;
               break;
        default:
              fprintf (stderr, "Error: bad image format. cannot read ...\n");
              exit (-1);
              break;
    }
#ifdef USELM
    MRDEMO.test(Image.nc(), Image.nl());
#endif
}

/**************************************************************************/

void io_write_ima_complex_f  (char *File_Name, Icomplex_f &Image)
{
    char File_Name_cf[100];
    if (Format_Imag == F_UNKNOWN)
                     Format_Imag = io_detect_format (File_Name);
    switch (Format_Imag)
    {
        case F_DISP: 
#if IO_DISP_OK
              IO_DISP.write_ima_complex_f (File_Name, Image);
#else
              fprintf (stderr, "Error: DISP is not active\n");
              exit (-1); 
#endif
              break;
        case F_MIDAS:
#if IO_MIDAS_OK
              IO_MIDAS.write_ima_complex_f (File_Name, Image);
#else
              fprintf (stderr, "Error: MIDAS is not active\n");
              exit (-1); 
#endif
              break;
        case F_FITS:
        case F_GIF:
        case F_PGM:
        case F_JPEG:
             { 
              char Ext[10];
              if (Format_Imag == F_FITS) strcpy(Ext, ".fits");
              else if (Format_Imag == F_GIF) strcpy(Ext, ".gif");
              else if (Format_Imag == F_PGM) strcpy(Ext, ".pgm");
              else if (Format_Imag == F_JPEG) strcpy(Ext, ".jpg");
              Ifloat Buff(Image.nl(), Image.nc(), "buff");
              strcpy (File_Name_cf, File_Name);
              strcat (File_Name_cf, "_re");
              strcat (File_Name_cf, Ext);
              real(Buff, Image);
              io_write_ima_float(File_Name_cf, Buff);
              strcpy (File_Name_cf, File_Name);
              strcat (File_Name_cf, "_im");
              strcat (File_Name_cf, Ext);
              imag(Buff, Image);
              io_write_ima_float(File_Name_cf, Buff);
             }
              break;
        default:
              fprintf (stderr, "Error: bad image format. cannot write ...\n");
              exit (-1);
              break;
    }
}

