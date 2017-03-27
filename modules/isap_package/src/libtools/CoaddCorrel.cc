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
**    File:     CoaddCorrel.cc
**
************************************************************
**
**  Tools for coaddition from correlation
**  
************************************************************/


#include <stdio.h>
#include <stdlib.h>
#include "GlobalInc.h"
#include "CoaddCorrel.h"

/*
 * Declarations private to this module
 */


/* #defines */
#define MAX_PIX_VALUE       (pixelvalue)2147483647L
#define APO_PARABOLE(a,b,c) (0.5*(((a)-(c))/((a)-(2.0*(b))+(c))))


/* Number of tabulations in kernel  */
#define TABSPERPIX      (10000)
#define KERNEL_WIDTH    (2.0)
#define KERNEL_SAMPLES  (1+(int)(TABSPERPIX * KERNEL_WIDTH))
#define TANH_STEEPNESS  (5.0)


/* PRIVATE function prototypes */

/*
 * xcorrelation module
 */
static double
private_correlation(pixelvalue *, pixelvalue *, int, int, int, int,
                    int, int, int, int, int, int, double *, double *) ;
static double apodisation(double, double, double) ;
static double sinc(double) ;



/*
 * warping module
 */
static double * generate_interpolation_kernel(char *);
static double * generate_tanh_kernel(double) ;
static void reverse_tanh_kernel(double *, int) ;
int
warp_image(
    pixelvalue *, int, int,
    pixelvalue *, int, int,
    double, double, double,
    char *);


/*
 * PUBLIC FUNCTIONS
 */
 
/*---------------------------------------------------------------------------
   Function :   cube_get_offset()
   In       :   input list of planes
                input pattern image
                size of the input cube: lx, ly, lz
                filled up xcorr_def structure
   Out      :   int 0 if everything went Ok, -1 otherwise
   Job      :   determine the offsets between planes in a given cube
   Notice   :   xcorr_def.xc_method defined in .h file.
 ---------------------------------------------------------------------------*/
 
int cube_get_offset(pixelvalue ** cube, pixelvalue  * pattern,
                    int lx, int ly, int lz, int PosX, int PosY, xcorr_def  *xc)
{
    int     i ;
    double  dx, dy, level ;

    if (xc==NULL || cube==NULL || pattern==NULL) return -1 ;
    switch (xc->xc_method) {
        default:
        for (i=0 ; i<lz ; i++) {
            level = private_correlation( pattern, cube[i],
                    lx, lx, PosX, PosY,  PosX, PosY,
                    xc->xc_dx, xc->xc_dy, xc->xc_hx, xc->xc_hy,
                    &dx, &dy) ;
            xc->xc_x[i] = dx ;
            xc->xc_y[i] = dy ;
            xc->xc_level[i] = level ;
        }
        break ;
    }
    return 0 ;
}


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
   Notice   :   interp_method values to be found in .h file
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
)
{
    int     i ;

    if (cube_in==NULL || cube_out==NULL || xc==NULL) return -1 ;
    for (i=0 ; i<lz ; i++) {
       printf("RESAMPLE %d: x = %f and y = %f\n",i,(xc->xc_x)[i],(xc->xc_y)[i]);
        warp_image(cube_in[i], lx, ly,
                   cube_out[i], lx_out, ly_out,
                   (xc->xc_x)[i],
                   (xc->xc_y)[i],
                   (double)lx_out / (double)lx,
                   interp_method) ;
    }
    return 0 ;
}




/*
 * PRIVATE FUNCTIONS
 */


/*----------------------------------------------------------------------------
 * Function :   private_correlation()
 * In       :   2 image buffers + search information
 * Out      :   best estimated position, image difference on this pixel
 * Job      :   estimate the place of minimum of squared differences
 *              between 2 image buffers
 * Notice   :   estimation is made up to subpixel accuracy
 *--------------------------------------------------------------------------*/

static double private_correlation(
    pixelvalue    *buff_in1,    /* Search is made on this buffer          */
    pixelvalue    *buff_in2,    /* This is the pattern we are looking for */
    int     lx1,            /* nb pixels/row buffer 1                       */
    int     lx2,            /* nb pixels/row buffer 2                       */
    int     at_x1,          /* Search center in buffer 1                    */
    int     at_y1,
    int     at_x2,          /* Search center in buffer 2                    */
    int     at_y2,
    int     dx_max,         /* Search area size                             */
    int     dy_max,
    int     hx,             /* Measure area size                            */
    int     hy,
    double   *dx,            /* Returned apodized position in x              */
    double   *dy )           /* Returned apodized position in y              */
{
    pixelvalue    *buffer_in1 = buff_in1;
    pixelvalue    *buffer_in2 = buff_in2;
    int                     i,j,k,l;
    double                    *distances;
    int                     inc1,inc2;
    double                   inv_surface ;
    int                     k_min,l_min ;
    double                   somme_min;
    register pixelvalue     *reg1,*reg2;
    register double           value,somme;
    double                   inc_x, inc_y ;
    int                     pos_min ;
    double                   best_distance ;

    /* Error handling: test entries */

    if ((buffer_in1 == (pixelvalue*)NULL)||(buffer_in2 == (pixelvalue*)NULL))
    {
        fprintf(stderr,"Error: NULL buffer in entry: aborting correlation") ;
	exit(-1);
    }

    if ((at_x1<=0) || (at_x1>=lx1) || (at_x2<=0) || (at_x2>=lx2)) 
    {
        fprintf(stderr,"Error: value out of bounds for requested correlation center") ;
        exit(-1) ;
    }

    /* More error handling should be made there to check the whole
     * area is within bounds.
     */

    distances = (double *) malloc((2*dy_max+1)*(2*dx_max+1)*sizeof(double));
    somme_min = (double)MAX_PIX_VALUE*(double)MAX_PIX_VALUE*
                (double)((2*hy+1)*(2*hx+1)) ;
    k_min = l_min = 0 ;

    /* Move into the buffers, up to the requested searching place   */
    buffer_in1 += at_x1 + at_y1*lx1 ;
    buffer_in2 += at_x2 + at_y2*lx2 ;
    inc1 = lx1-hx-hx-1;
    inc2 = lx2-hx-hx-1;

    inv_surface = 1.0 / ((double)(2*hx+1)*(double)(2*hy+1)) ;
    // cout << "hx = " << hx << " hy = " << hy << endl;
    // cout << "lx1 = " << lx1 << " lx2 = " << lx2 << endl;
    // cout << "dy_max = " << dy_max << " dx_max = " << dx_max << endl;
    for (l=-dy_max;l<=dy_max;l++)
    for (k=-dx_max;k<=dx_max;k++) 
    {
        somme = 0;
        reg1 = buffer_in1+k-hx+(l-hy)*lx1;
        reg2 = buffer_in2-hx-hy*lx2;
	
	// cout << l << " " << k << " V1 = " << *reg1 << " V2 = " << *reg2 << endl;
        for (j=-hy;j<=hy;j++) 
	{ 
           for (i=-hx;i<=hx;i++) 
	   {
               value = (double)(*reg1++)-(double)(*reg2++);
               value *= value;
               somme += value;
           }
           reg1+=inc1;
           reg2+=inc2;
        }
	// cout << "S = " << somme << " Sm = " << somme_min << endl;
        if (somme<somme_min) 
	{
            l_min = l;
            k_min = k;
            somme_min = somme;
        }
        distances[dx_max+k+(2*dx_max+1)*(dy_max+l)] = somme * inv_surface ;
    }

    /* apodisation : sub pixel correlation */
    pos_min = dx_max+k_min+(2*dx_max+1)*(dy_max+l_min) ;
    best_distance = distances[pos_min] ;

    /* Take care of edge effects in measure */
    if ((k_min == -dx_max)||(k_min == dx_max)) {
        inc_x = 0.0 ;
    } else {
        inc_x = apodisation(    distances[pos_min-1],
                                distances[pos_min],
                                distances[pos_min+1]) ;
    }

    if ((l_min == -dy_max)||(l_min == dy_max)) {
        inc_y = 0.0 ;
    } else {
        inc_y = apodisation(    distances[pos_min-(2*dx_max+1)],
                                distances[pos_min],
                                distances[pos_min+(2*dx_max+1)]) ;
    }

    *dx = (double)k_min + inc_x  ;
    *dy = (double)l_min + inc_y ;

    /* For debug purposes only  */
    // cout << "DX DY = " << *dx << " " << *dy << endl;
    // for (j=0 ; j<2*dy_max+1 ; j++) {
    //     for (i=0 ; i<2*dx_max+1 ; i++) {
    //        printf("%d %d %g\n",
    //                i-dx_max, j-dy_max, distances[i+j*(2*dx_max+1)]) ;
    //    }
    //    printf("\n") ;
    // }
    

    free(distances) ;
    return(best_distance) ;
}


/*----------------------------------------------------------------------------
 * Function :   apodisation()
 * In       :   3 best correlation distances
 * Out      :   subpixel accurate position
 * Job      :   subpixel accuracy in correlation
 * Notice   :
 *--------------------------------------------------------------------------*/

static double apodisation(
    double   d1,
    double   d2,
    double   d3 )
{
    double   value ;

    if (d1==d2) return (-0.5) ;
    if (d2==d3) return (0.5) ;
    if (d1==d3) return (0.0) ;
    value = (double)APO_PARABOLE(d1, d2, d3) ;
    if (value>0.5) value = 0.5 ;
    if (value<-0.5) value = -0.5 ;
    return value ;
}



/*---------------------------------------------------------------------------
 * Function :   generate_interpolation_kernel()
 * Job      :   Computes tabulated values of the kernel
 * In       :   Kernel type
 * Out      :   Kernel values
 * Notice   :
 *--------------------------------------------------------------------------*/

static double   *
generate_interpolation_kernel(char * kernel_type)
{
    double  *   tab ;
    int         i ;
    double      x ;
    double      alpha ;
    double      inv_norm ;
    int         samples = KERNEL_SAMPLES ;

    if (!strcmp(kernel_type, "default")) {
        tab = generate_interpolation_kernel("tanh") ;
    } else if (!strcmp(kernel_type, "sinc")) {
        tab = (double*)malloc(samples * sizeof(double)) ;
        tab[0] = 1.0 ;
        tab[samples-1] = 0.0 ;
        for (i=1 ; i<samples ; i++) {
            x = (double)KERNEL_WIDTH * (double)i/(double)(samples-1) ;
            tab[i] = sinc(x) ;
        }
    } else if (!strcmp(kernel_type, "sinc2")) {
        tab = (double*)malloc(samples * sizeof(double)) ;
        tab[0] = 1.0 ;
        tab[samples-1] = 0.0 ;
        for (i=1 ; i<samples ; i++) {
            x = 2.0 * (double)i/(double)(samples-1) ;
            tab[i] = sinc(x) ;
            tab[i] *= tab[i] ;
        }
    } else if (!strcmp(kernel_type, "lanczos")) {
        tab = (double*)malloc(samples * sizeof(double)) ;
        for (i=0 ; i<samples ; i++) {
            x = (double)KERNEL_WIDTH * (double)i/(double)(samples-1) ;
            if (fabs(x)<2) {
                tab[i] = sinc(x) * sinc(x/2) ;
            } else {
                tab[i] = 0.00 ;
            }
        }
    } else if (!strcmp(kernel_type, "hamming")) {
        tab = (double*)malloc(samples * sizeof(double)) ;
        alpha = 0.54 ;
        inv_norm  = 1.00 / (double)(samples - 1) ;
        for (i=0 ; i<samples ; i++) {
            x = (double)i ;
            if (i<(samples-1)/2) {
                tab[i] = alpha + (1-alpha) * cos(2.0*M_PI*x*inv_norm) ;
            } else {
                tab[i] = 0.0 ;
            }
        }
    } else if (!strcmp(kernel_type, "hann")) {
        tab = (double*)malloc(samples * sizeof(double)) ;
        alpha = 0.50 ;
        inv_norm  = 1.00 / (double)(samples - 1) ;
        for (i=0 ; i<samples ; i++) {
            x = (double)i ;
            if (i<(samples-1)/2) {
                tab[i] = alpha + (1-alpha) * cos(2.0*M_PI*x*inv_norm) ;
            } else {
                tab[i] = 0.0 ;
            }
        }
    } else if (!strcmp(kernel_type, "tanh")) {
        tab = generate_tanh_kernel(TANH_STEEPNESS) ;
    } else {
        fprintf(stderr, "unrecognized kernel type [%s]: aborting generation",
                kernel_type) ;
        return (double*)NULL ;
    }

    return(tab) ;
}

/*---------------------------------------------------------------------------
 * Function :   sinc()
 * In       :   double
 * Out      :   double
 * Job      :   cardinal sine
 * Notice   :   rescaled by PI
 *--------------------------------------------------------------------------*/

static double
sinc(double x)
{
    if (fabs(x)<1e-4)
        return (double)1.00 ;
    else
        return ((sin(x * (double)M_PI)) / (x * (double)M_PI)) ;
}

/*
 * The following function is a good approximation of a box filter,
 * built up from a product of hyperbolic tangents. It has the following
 * properties:
 * 1. It converges very quickly towards +/- 1
 * 2. The converging transition is sharp
 * 3. It is infinitely differentiable everywhere (smooth)
 * 4. The transition sharpness is scalable
 *
 * It is defined as a macro here (yes, it's dirty) to optimize the code
 */

#define hk_gen(x,s) (((tanh(s*(x+0.5))+1)/2)*((tanh(s*(-x+0.5))+1)/2))


/*---------------------------------------------------------------------------
   Function :   generate_tanh_kernel()
   In       :   double: steepness of the hyperbolic tangent
   Out      :   pointer to (samples) doubles
   Job      :   generate an hyperbolic tangent kernel
   Notice   :   don't touch
 ---------------------------------------------------------------------------*/


static double * generate_tanh_kernel(double steep)
{
    double  *   kernel ;
    double  *   x ;
    double      width ;
    double      inv_np ;
    double      ind ;
    int         i ;
    int         np ;
    int         samples ;

    width   = (double)TABSPERPIX / 2.0 ;
    samples = KERNEL_SAMPLES ;
    np      = 32768 ; /* Hardcoded: should never be changed */
    inv_np  = 1.00 / (double)np ;

    /*
     * Generate the kernel expression in Fourier space
     * with a correct frequency ordering to allow standard FT
     */
    x = (double*)malloc((2*np+1)*sizeof(double)) ;
    for (i=0 ; i<np/2 ; i++) {
        ind      = (double)i * 2.0 * width * inv_np ;
        x[2*i]   = hk_gen(ind, steep) ;
        x[2*i+1] = 0.00 ;
    }
    for (i=np/2 ; i<np ; i++) {
        ind      = (double)(i-np) * 2.0 * width * inv_np ;
        x[2*i]   = hk_gen(ind, steep) ;
        x[2*i+1] = 0.00 ;
    }

    /* 
     * Reverse Fourier to come back to image space
     */
    reverse_tanh_kernel(x, np) ;

    /*
     * Allocate and fill in returned array
     */
    kernel = (double*)malloc(samples * sizeof(double)) ;
    for (i=0 ; i<samples ; i++) {
        kernel[i] = 2.0 * width * x[2*i] * inv_np ;
    }
    free(x) ;
    return kernel ;
}

/*---------------------------------------------------------------------------
   Function :   reverse_tanh_kernel()
   In       :   a tanh generated kernel in Fourier space
   Out      :   a reversed kernel in image space
   Job      :   transforms from Fourier to image space a tanh kernel
   Notice   :   optimized, only useful for generate_tanh_kernel()
 ---------------------------------------------------------------------------*/


#define KERNEL_SW(a,b) tempr=(a);(a)=(b);(b)=tempr
static void reverse_tanh_kernel(double * data, int nn)
{
    unsigned long   n,
                    mmax,
                    m,
                    i, j,
                    istep ;
    double  wtemp,
            wr,
            wpr,
            wpi,
            wi,
            theta;
    double  tempr,
            tempi;

    n = (unsigned long)nn << 1;
    j = 1;
    for (i=1 ; i<n ; i+=2) {
        if (j > i) {
            KERNEL_SW(data[j-1],data[i-1]);
            KERNEL_SW(data[j],data[i]);
        }
        m = n >> 1;
        while (m>=2 && j>m) {
            j -= m;
            m >>= 1;
        }
        j += m;
    }
    mmax = 2;
    while (n > mmax) {
        istep = mmax << 1;
        theta = 2 * M_PI / mmax;
        wtemp = sin(0.5 * theta);
        wpr = -2.0 * wtemp * wtemp;
        wpi = sin(theta);
        wr  = 1.0;
        wi  = 0.0;
        for (m=1 ; m<mmax ; m+=2) {
            for (i=m ; i<=n ; i+=istep) {
                j = i + mmax;
                tempr = wr * data[j-1] - wi * data[j];
                tempi = wr * data[j]   + wi * data[j-1];
                data[j-1] = data[i-1] - tempr;
                data[j]   = data[i]   - tempi;
                data[i-1] += tempr;
                data[i]   += tempi;
            }
            wr = (wtemp = wr) * wpr - wi * wpi + wr;
            wi = wi * wpr + wtemp * wpi + wi;
        }
        mmax = istep;
    }
}
#undef KERNEL_SW



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
)
{
    int         i, j, k ;
    double  *   kernel ;

    int         leaps[16] ;
    double      neighbors[16] ;

    double      inv_zoom ;
    double      x, y ;
    int         px, py ;
    int         pos ;
    int         tabx, taby ;
    double      rsc[8], sumrs ;
    double      cur ;

    if (image_in==NULL || image_out==NULL) return -1 ;
    if ((kernel=generate_interpolation_kernel(interp_method))==NULL)
        return -1 ;

    /* precompute leaps for 16 closest neighbor pixel positions */

    leaps[0] = -1 - lx ;
    leaps[1] =    - lx ;
    leaps[2] =  1 - lx ;
    leaps[3] =  2 - lx ;

    leaps[4] = -1 ;
    leaps[5] =  0 ;
    leaps[6] =  1 ;
    leaps[7] =  2 ;

    leaps[8] = -1 + lx ;
    leaps[9] =      lx ;
    leaps[10]=  1 + lx ;
    leaps[11]=  2 + lx ;

    leaps[12]= -1 + 2*lx ;
    leaps[13]=      2*lx ;
    leaps[14]=  1 + 2*lx ;
    leaps[15]=  2 + 2*lx ;

    inv_zoom = 1.0 / zoom ;

    /* Double loop on the output image  */
    for (j=0 ; j < ly_out ; j++) {
        for (i=0 ; i< lx_out ; i++) {
            /* Compute the original source for this pixel   */
            x = (double)i * inv_zoom - dx ;
            y = (double)j * inv_zoom - dy ;

            /* Which is the closest integer positioned neighbor?    */
            px = (int)x ;
            py = (int)y ;

            if ((px < 1) || (px > (lx-3)) || (py < 1) || (py > (ly-3))) {
                image_out[i+j*lx_out] = (pixelvalue)0.0 ;
            } else {
                /* Now feed the positions for the closest 16 neighbors  */
                pos = px + py * lx ;
                for (k=0 ; k<16 ; k++) {
                    neighbors[k] = (double)image_in[pos+leaps[k]] ;
 		    //if((i == lx_out/2) && (j==ly_out/2)) {
                    //   printf("[%d,%d],[%f,%f]->[%d,%d] [%d], Neighbors[%d]=%f\n",i,j,x,y,px,py,pos,k,neighbors[k]); 
                    //}
                }

                /* Which tabulated value index shall we use?    */
                tabx = (int)((x - (double)px) * (double)(TABSPERPIX)) ;
                taby = (int)((y - (double)py) * (double)(TABSPERPIX)) ;

                /* Compute resampling coefficients  */
                /* rsc[0..3] in x, rsc[4..7] in y   */

                rsc[0] = kernel[TABSPERPIX + tabx] ;
                rsc[1] = kernel[tabx] ;
                rsc[2] = kernel[TABSPERPIX - tabx] ;
                rsc[3] = kernel[2 * TABSPERPIX - tabx] ;
                rsc[4] = kernel[TABSPERPIX + taby] ;
                rsc[5] = kernel[taby] ;
                rsc[6] = kernel[TABSPERPIX - taby] ;
                rsc[7] = kernel[2 * TABSPERPIX - taby] ;
                //if ((i==lx_out/2)&&(j==ly_out/2)) {
                //   printf("Tabx=%d, Taby=%d, rsc = %f,%f,%f,%f,%f,%f,%f,%f\n",tabx,taby,rsc[0],rsc[1],rsc[2],rsc[3],rsc[4],rsc[5],rsc[6],rsc[7]);
                //}             

                sumrs = (rsc[0]+rsc[1]+rsc[2]+rsc[3]) *
                        (rsc[4]+rsc[5]+rsc[6]+rsc[7]) ;

                /* Compute interpolated pixel now   */
                cur =   rsc[4] * (  rsc[0]*neighbors[0] +
                                    rsc[1]*neighbors[1] +
                                    rsc[2]*neighbors[2] +
                                    rsc[3]*neighbors[3] ) +
                        rsc[5] * (  rsc[0]*neighbors[4] +
                                    rsc[1]*neighbors[5] +
                                    rsc[2]*neighbors[6] +
                                    rsc[3]*neighbors[7] ) +
                        rsc[6] * (  rsc[0]*neighbors[8] +
                                    rsc[1]*neighbors[9] +
                                    rsc[2]*neighbors[10] +
                                    rsc[3]*neighbors[11] ) +
                        rsc[7] * (  rsc[0]*neighbors[12] +
                                    rsc[1]*neighbors[13] +
                                    rsc[2]*neighbors[14] +
                                    rsc[3]*neighbors[15] ) ;

                /* Affect the value to the output image */
                image_out[i+j*lx_out] = (pixelvalue)(cur/sumrs) ;
                /* done ! */
            }
        }
    }
    free(kernel) ;
    return 0 ;
}




#ifdef MAINTEST
int main(int argc, char *argv[])
{
}
#endif
