/*******************************************************************************
**
**    UNIT
**
**    Version: 1.0
**
**    Author: Jean-Luc Starck
**
**    Date:  11/04/00 
**    
**    File:  IM_IOCol.h
**
*******************************************************************************
**
**    DESCRIPTION  Tiff, Gif, and JPEG image definition
**    ----------- 
**                 
******************************************************************************/

#ifndef _IM_IOCOL_H_
#define _IM_IOCOL_H_

typedef unsigned char byte;
typedef unsigned short u_short;
typedef unsigned char u_char;
typedef unsigned int u_int;

//  Gif format
# define PARM(a) a
# define PIC8  0
# define PIC24 1

#define F_FULLCOLOR 0
#define F_BWDITHER  2
#define F_GREYSCALE 1

/* MONO returns total intensity of r,g,bfor (i=0; i<256; i++) pinfo.g[i] = (byte) i;
   for (i=0; i<256; i++) pinfo.r[i] = (byte) i;
   for (i=0; i<256; i++) pinfo.b[i] = (byte) i;
     (i = .33R + .5G + .17B) */
#define MONO(rd,gn,bl) ( ((int)(rd)*11 + (int)(gn)*16 + (int)(bl)*5) >> 5)


/* info structure filled in by the LoadXXX() image reading routines */
typedef struct { byte *pic;                  /* image data */
	         int   w, h;                 /* pic size */
		 int   type;                 /* PIC8 or PIC24 */

		 byte  r[256],g[256],b[256];
		                             /* colormap, if PIC8 */

		 int   normw, normh;         /* 'normal size' of image file
					        (normally eq. w,h, except when
						doing 'quick' load for icons */

		 int   frmType;              /* def. Format type to save in */
		 int   colType;              /* def. Color type to save in  */
                    /* also called colorType value F_FULLCOLOR, F_GREYSCALE, */
		 char  fullInfo[128];        /* Format: field in info box */
		 char  shrtInfo[128];        /* short format info */
		 char *comment;              /* comment text */

		 int   numpages;             /* # of page files, if >1 */
		 char  pagebname[64];        /* basename of page files */
	       } PICINFO;


#define xvbzero(s,size) memset(s,0,size)

// Gif format
inline byte float_to_byte(float V)
{
   byte Vb;
   if (V > 255) Vb = 255;
   else if (V < 0) Vb = 0;
   else Vb = (byte) V; 
   return Vb;
}
inline byte int_to_byte(int V)
{
   byte Vb;
   if (V > 255) Vb = 255;
   else if (V < 0) Vb = 0;
   else Vb = (byte) V; 
   return Vb;
} 

 void gifWarning(char *st);
char *gifname(char * NameStep);
char *tiffname(char * NameStep);
int io_read3d_tiff(char *Name,  PICINFO *pinfo);
void io_write3d_tiff( char *Name, PICINFO *pinfo);
#endif
