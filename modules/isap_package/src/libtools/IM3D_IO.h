/*******************************************************************************
**
**    UNIT
**
**    Version: 1.0
**
**    Author: Jean-Luc Starck
**
**    Date:  06/08/99
**    
**    File:  IM3D_IO.h
**
*******************************************************************************
**
**    DESCRIPTION  IO 3D routine
**    ----------- 
**                 
******************************************************************************/

#ifndef _IM3D_IO_H_
#define _IM3D_IO_H_

#include"GlobalInc.h"
#include "IM_IO.h"
#include"writefits3d.h"

#if IO_JPEG_OK
#include "IM_JPEG.h"
#endif

typedef unsigned short u_short;
typedef unsigned char  u_char;
typedef unsigned int   u_int;

enum type_3d_format {F3D_UNKNOWN, F3D_TIFF, F3D_FITS, F3D_GIF, F3D_JPEG};
#define DEFAULT_3D_FORMAT  F3D_FITS

inline const char *String3DFormat (type_3d_format  type)
{
    switch (type)
    {
        case F3D_FITS: 
              return ("FITS format");break;
        case F3D_GIF: 
              return ("GIF format");break;
 	case F3D_JPEG:
	      return ("JPEG format");break;
 	case F3D_TIFF:
	      return ("TIFF format");break;
        default: 
              return ("Unknown format");break;
    }
}

void io_3d_read_data(char *File_Name, fltarray & Dat, fitsstruct *FitsHeader=NULL);
void io_3d_write_data(char *File_Name, fltarray & Dat,fitsstruct *FitsHeader=NULL);

type_3d_format  io_which_3d_format(const char *Format);
type_3d_format  io_detect_3dformat(const char *Format);
void io_3d_set_format(const char *Format);

int io_read3d_tiff(char *name, fltarray & Data);
int io_read2d_tiff(char *name, fltarray & Data);
int io_write3d_tiff(  char *name, fltarray & Data);

void io_read3d_gif(char *File_Name, fltarray & Dat);
void io_write3d_gif(char *File_Name, fltarray & Dat, Bool Norm=False);

class IO3DInfoData {
  public:
    Bool PseudoColRead;
    unsigned int AllocBufferSize;
    int NbrReadLine;      // For block per block reading, we need to know
    int LastBlockLineNbr; // which block is in memory (JPEG format)
    int NbrWriteLine;     // Number of written lines
    char File_Name[256];
    fitsstruct *PtrFits;
    PICINFO Pinfo;
#if IO_JPEG_OK
    READ_JPEG RJEP;
    WRITE_JPEG WJEP;
#endif
    Bool PseudoCol;
    Bool FullCol;
    int Nl;    
    int Nc;
    int Nima;
    type_3d_format Format;
    type_data Type;
    IO3DInfoData()  { 
       AllocBufferSize=0;
       PseudoColRead=False;
       PtrFits = NULL;
       Pinfo.pic = NULL;
       FullCol = PseudoCol = False;
       Nl = Nc = 0;
       Nima = 1;
       Format = F3D_UNKNOWN;
       Type = UNKNOWN;
       File_Name[0]='\0';
       Pinfo.type = 0;
       Pinfo.comment = new char [1];
       Pinfo.comment = '\0';
       Pinfo.colType = F_FULLCOLOR;
       NbrReadLine = LastBlockLineNbr = NbrWriteLine = 0;
    }
 

    void init_writing(char *fname, int Nlima, int Ncima, int NbrIma=3);
    void init_reading(char *FileName, Bool ReadInPseudo=False);
    void end_writing();
    void end_reading();
    void read_pseudo_block(intarray &Data, int Indi, int Indj);
    void read_pseudo_block(fltarray &Data, int Indi, int Indj);
    void read_col_block(intarray &Data, int Indi, int Indj);
    void read_col_block(fltarray &Data, int Indi, int Indj);
    void write_pseudo_block(intarray &Data, int Indi, int Indj);
    void write_pseudo_block(fltarray &Data, int Indi, int Indj);
    void write_col_block(intarray &Data, int Indi, int Indj);
    void write_col_block(fltarray &Data, int Indi, int Indj);

    ~IO3DInfoData() 
    { 
	if (Pinfo.pic  != NULL) delete [] Pinfo.pic;
	PtrFits = NULL;
        Pinfo.pic = NULL;
    }
};

void io_3d_read_block_ima(char *File_Name, fltarray &Image, 
                       int Indi, int Indj, IO3DInfoData &InfoDat, Bool NoBscale=False);
void io_3d_write_block_ima(char *File_Name, fltarray &Image, 
                    int Indi, int Indj, IO3DInfoData &InfoDat, Bool NoBscale=False);

void close3dgif(char *File_Name);
void close3dgif(char *File_Name, PICINFO *GifIma);
void write_block3d_gif (fltarray  & Data, int Indi, int Indj);
void io_write3d_block_gif (fltarray & Image, int Indi, int Indj, PICINFO *GifIma);

#endif
