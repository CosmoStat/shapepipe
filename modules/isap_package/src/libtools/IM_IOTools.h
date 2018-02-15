
/******************************************************************************
**                   Copyright (C) 1994 by CEA
*******************************************************************************
**
**    UNIT
**
**    Version: 3.2
**
**    Author: Jean-Luc Starck
**
**    Date:  96/05/07 
**    
**    File:  IM_IOTools.h
**
*******************************************************************************
**
**    DESCRIPTION  
**    ----------- 
**                 
**    PARAMETRES    
**    ----------    
** 
**    RESULTS      
**    -------  
**
**
******************************************************************************/

#ifndef _IOTOOLS_
#define _IOTOOLS_

#include<ctype.h>
#include<math.h>
#include<stdio.h>
#include<stdlib.h>
#include<string.h>

#define MAXCHAR 256
#define RETURN_OK 0
#define RETURN_ERROR (-1)
#define RETURN_FATAL_ERROR (-2)
#ifdef  NOSMALLHUGE
#define BIG 1e+30   /* a huge number */
#else
#define BIG HUGE_VAL
#endif

#ifndef SEEK_SET
#define SEEK_SET 0
#endif
#ifndef SEEK_CUR
#define SEEK_CUR 1
#endif

#ifndef EXIT_SUCCESS
#define EXIT_SUCCESS 0
#endif
#ifndef EXIT_FAILURE
#define EXIT_FAILURE -1
#endif

/*------------------- a few definitions to read FITS parameters ------------*/

#define FBSIZE  2880L   /* size (in bytes) of one FITS block */

#define FITSTOF(k, def) \
                        ((point = fitsnfind(buf, k, n))? \
                                 atof(strncpy(st, &point[10], 70)) \
                                :(def))

#define FITSTOI(k, def) \
                        ((point = fitsnfind(buf, k, n))? \
                                 atoi(strncpy(st, &point[10], 70)) \
                                :(def))
#define FITSTOS(k, str, def) \
                        { point = fitsnfind(buf, k, n); \
                          if (point != NULL) \
                                { \
                                for (i=0,point+=11; (*point)!='\'' && i < 69;) \
                                        (str)[i++] = *(point++); \
                                (str)[i] = '\0'; \
                                } \
                          else\
                                strcpy(str, def); \
                        }
#define QFREAD(ptr, size, file, fname) \
                if (fread(ptr, (size_t)(size), (size_t)1, file)!=1) \
                  error(EXIT_FAILURE, "*Error* while reading ", fname)

#define QFWRITE(ptr, size, file, fname) \
                if (fwrite(ptr, (size_t)(size), (size_t)1, file)!=1) \
                  error(EXIT_FAILURE, "*Error* while writing ", fname)

#define QFSEEK(file, offset, pos, fname) \
                if (fseek(file, (offset), pos)) \
                  error(EXIT_FAILURE,"*Error*: file positioning failed in ", \
                        fname)
#define QFTELL(pos, file, fname) \
                if ((pos=ftell(file))==-1) \
                  error(EXIT_FAILURE,"*Error*: file position unknown in ", \
                        fname)

/* int     t_size[] = {1, 2, 4, 4, 8}; */
typedef enum {H_INT, H_FLOAT, H_EXPO, H_BOOL, H_STRING, H_COMMENT,
                        H_KEY}  h_type;         /* type of FITS-header data */


extern void swapbytes(void *ptr, int nb, int n);
extern void    error(int num, char *msg1, char *msg2);

#endif
