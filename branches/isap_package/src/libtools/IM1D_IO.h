/*******************************************************************************
**
**    UNIT
**
**    Version: 1.0
**
**    Author: Jean-Luc Starck
**
**    Date:  22/9/98 
**    
**    File:  IM1D_IO.h
**
*******************************************************************************
**
**    DESCRIPTION  FITS Include
**    ----------- 
**                 
******************************************************************************/

#ifndef _IM1D_IO_H_
#define _IM1D_IO_H_
#include "IM_IO.h"

enum type_1d_format {F1D_UNKNOWN, F1D_FITS, F_1DASCII, F_1DEXEL};
#define DEFAULT_1D_FORMAT  F1D_FITS

inline const char *String1DFormat (type_1d_format  type)
{
    switch (type)
    {
        case  F1D_FITS: 
              return ("FITS format");break;
        case  F_1DASCII: 
              return ("ASCII format");break;
        case  F_1DEXEL: 
              return ("EXEL format");break;
        default: 
              return ("Unknown format");break;
    }
}

void io_1d_read_data(char *File_Name, fltarray & Dat, int NbrSkip=0);
void io_1d_read_data(char *File_Name, fltarray & Dat, fitsstruct *Header);
void io_1d_read_data(char *File_Name, fltarray & Dat, int NbrSkip,  fitsstruct *Header);

void io_1d_write_data(char *File_Name, fltarray & Dat,fitsstruct *Header=NULL);

void io_read1d_ascii(char *Name_Dat_In,  fltarray & Dat, 
                     int NbrSkip=0, Bool Exell = False);
void io_write1d_ascii(char *Name_Dat_Out,  fltarray & Dat, Bool Exell= False);
void io_read2d_ascii(char *Name_Dat_In,  fltarray &Dat);
void io_write2d_ascii(char *Name_Dat_Out,  fltarray &Dat);
void io_write_ascii(char *Name_Dat_Out,  fltarray &Dat, Bool Exell= False);

type_1d_format  io_which_1d_format(const char *Format);
type_1d_format  io_detect_1dformat(const char *Format); 
void io_1d_set_format(const char *Format);

#endif
