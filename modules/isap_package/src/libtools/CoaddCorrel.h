/***********************************************************
**      Copyright (C) 1999 CEA
************************************************************
**
**    UNIT
**
**    Version:  1.0
**
**    Author:   N. Devillars
**
**    Date:     27/08/99
**    
**    File:     CoaddCorrel.h
**
************************************************************
**
**  Tools for coaddition from correlation
**  
************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#ifndef _COADDCORR_H_
#define _COADDCORR_H_

typedef float pixelvalue ;

#define NBR_CORR_INTERP 5
enum type_interp_corr {ICF_SINC2,ICF_LANCZOS,ICF_HAMMING,ICF_HANN,ICF_TANH};

inline const char *StringCorrInterp (type_interp_corr type)
{
    switch (type)
    {
        case ICF_SINC2:
              return ("sinc2"); break;
        case ICF_LANCZOS:
              return ("lanczos");break;
	case ICF_HAMMING:
               return ("hamming");break;
        case ICF_HANN:
              return ("hann"); break;
        case ICF_TANH:
              return ("tanh");break;
    }
    return ("Error: bad type of interpolation");
}

		
/*
 * Cross-correlation result structure.
 * Contains everything needed for a cross-correlation definition
 */
 
		
		
typedef struct _XCORR_DEF_ {

    /* Number of points in the structure */
    int         xc_np ;

    /* Stores the offsets in x and y */
    double  *   xc_x ;
    double  *   xc_y ;
    /* Stores x-correlation levels */
    double  *   xc_level ;

    /* Half-sizes of the measure domain */
    int         xc_hx ;
    int         xc_hy ;

    /* Half-sizes of the search domain */
    int         xc_dx ;
    int         xc_dy ;

    /* Flag: 1 if xc_x and xc_y contain sensible values, 0 else */
    int         xc_init ;

    /* Name of the estimation method */
    /* not used at all yet */
    int         xc_method ;

} xcorr_def ;



/*---------------------------------------------------------------------------
   Function :   cube_get_offset()
   In       :   input list of planes
                input pattern image
                size of the input cube: lx, ly, lz
                filled up xcorr_def structure
   Out      :   int 0 if everything went Ok, -1 otherwise
   Job      :   determine the offsets between planes in a given cube
   Notice   :
 ---------------------------------------------------------------------------*/




int cube_get_offset(
    pixelvalue ** cube,
    pixelvalue  * pattern,
    int           lx,
    int           ly,
    int           lz,
    int           PosX,
    int           PosY,
    xcorr_def   * xc
);


/*---------------------------------------------------------------------------
   Function :   cube_resample()
   In       :   input list of planes
                (allocated) output list of planes
                (allocated and meaningful) xcorr_def structure
                input & output cube sizes
                interpolation method name
   Out      :   int 0 if Ok, -1 otherwise
   Job      :   resampling consecutive planes of a cube according to an
                xcorr_def contents.
   Notice   :
                interp_method can take one of the following values:
                "default" -> "tanh"
                "sinc2"
                "lanczos"
                "hamming"
                "hann"
                "tanh"


 ---------------------------------------------------------------------------*/



int cube_resample(
    pixelvalue ** cube_in,
    pixelvalue ** cube_out,
    xcorr_def   * xc,
    int           lx,
    int           ly,
    int           lz,
    int           lx_out,
    int           ly_out,
    char        * interp_method
);


/*---------------------------------------------------------------------------
   Function :   warp_image()
   In       :   pointer to input image
                pointer to (allocated) output image
                size of both images
                required offset and zoom
   Out      :   0 if Ok, -1 otherwise
   Job      :   resample an image with a given transformation
   Notice   :
 ---------------------------------------------------------------------------*/

int
warp_image(
    pixelvalue  *   image_in,
    int             lx,
    int             ly,
    pixelvalue  *   image_out,
    int             lx_out,
    int             ly_out,
    double          dx,
    double          dy,
    double          zoom,
    char        *   interp_method
);

#endif
