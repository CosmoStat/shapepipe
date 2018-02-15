#ifndef SR_UTIL_H
#define SR_UTIL_H
#include <stdio.h>
#include <stdlib.h>
#include "TempArray.h"
#include <cmath>
#include "IM_Obj.h"
#include "IM_IO.h"
#include "MR_OBJ.h"
#include "MR_Sigma.h"
#include "CoaddCorrel.h"
#include "IM_Deconv.h"
//#include <random>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

const double Pi = 3.14159265358979323846;
const int LR_DENOISING = 0;
const int MR_SUPP_DENOISING=1;
//gsl_rng *rng;
//const gsl_rng_type * T;
//T =  gsl_rng_default;
//rng = gsl_rng_alloc (T);

typedef struct{
  type_transform transf;
  int nb_sc;
  int nb_usc;
  sb_type_norm Norm;
  type_undec_filter U_Filter;
  type_sb_filter SB_Filter;
  type_border Bord;
  FilterAnaSynt *FAS ;
}mr_opt;
int compute_centroid(Ifloat *img,double sig,double *cent,int niter=10);
int compute_centroid_arr(fltarray *data,double sig,double **cent,int niter=10);
void thresholding(Ifloat *data,Ifloat *data_out,Ifloat *thresh,bool thresh_type=0);
void thresholding(fltarray *data_in,fltarray *data_out,fltarray *thresh,bool thresh_type=0);
void circ_thresh(Ifloat *data_in,Ifloat *data_out,double r,double *cen);
void circ_thresh(fltarray *data_in,fltarray *data_out,double r,double **cen);
int sign_num (double a);
int wl_trans(Ifloat *Dat,mr_opt opt,MultiResol* MR_Data);
void wl_thresholding(MultiResol*wav_coeff,fltarray *thresh,bool thresh_type=0);
void wl_filter(Ifloat *img,Ifloat*img_filt,mr_opt opt,double nsig=4,bool thresh_type=0,fltarray* noise_map=NULL,fltarray *coeff_thresh=NULL,Ifloat *coarse_scale=NULL);
int wl_gnoise_est(MultiResol* wav_coeff,fltarray* noise_arr);
mr_opt mr_opt_init(int Nx,type_transform transf = TO_UNDECIMATED_MALLAT,int nb_sc=-1,int nb_usc=-1,sb_type_norm Norm=NORM_L1,type_undec_filter U_Filter=DEF_UNDER_FILTER,type_sb_filter SB_Filter=F_MALLAT_7_9,type_border Bord=I_CONT,Bool Verbose =False);
int noise_est(fltarray* data,fltarray* noise_arr);
void wl_arr_filter(fltarray *data,fltarray*data_filt,mr_opt opt,double nsig=4,bool thresh_type=0);
void coadd_all_frames(fltarray Data, Ifloat &Ima,Bool GetMedian=False);
void coadd_all_frames(fltarray Dat, Ifloat &Ima,float Zoom,
		      Bool GetMedian=True,int MaxDist=4, 
		      int Surface=10, type_interp_corr TypeInterp=ICF_TANH,Bool OptNoOffset=False,Bool MeanSub=False,Bool Verbose=False);
void flux_est(Ifloat *data, double r,double *cent,double*flux);
void flux_est(fltarray *data, double r,double **cent,double*flux,int ref_im_ind=0);
int gnoise_est(Ifloat* data,Ifloat* noise_arr);
int gnoise_est(fltarray* data,fltarray* noise_arr);
void shift_est(fltarray *data, double sig_gfit, double *sig_vec,double **cent,double **shift,int ref_im_ind=0,double nsig=4);
double sinc(double x);
void lanczos(double u,double*mask,int mask_rad=10);
void lanczos(double u[],double**mask,int mask_rad=10);
void decim(Ifloat *img,Ifloat *img_filt,Ifloat *img_dec,int D,Ifloat* mask,Bool filt_en=True,Bool fft_en=True ); 
void rotate(Ifloat *input,Ifloat *output,int dir);
void transpose_decim(Ifloat*im_in,Ifloat*im_out,int D);
void power_meth(void(*opname)(double*,double*),int size_input,int *nb_iter,double *spec_rad,int nb_iter_max=100,double tol=0.01);
double norm2 (double * x,int siz);
void scale (double * x,double siz,double a);
void ineq_cons(Ifloat &A, Ifloat B,double b=1);
void reverse(double*u,int length,double*ur);
void holes(double*h,int sc_nb,int length,double*hout);
void convol1D(double*in,double*out,double*h,int length,int filt_length);
void sep_filt2d(Ifloat* im_in,Ifloat*im_out,double *h,double *g,int lh,int lg);
void transp_sep_filt2d(Ifloat* im_in,Ifloat*im_out,double *h,double *g,int lh,int lg);
//void randomn(double *x,double sig=1,int length=1,double mean=0);
void randomngsl(double *x,double sig,int length,double mean=0);
//void randomn(Ifloat *mat,double sig=1,double mean=0);
void randomngsl(Ifloat *mat,double sig=1,double mean=0);
//void randomn(fltarray *mat,double *sig,double mean=0);
void randomngsl(fltarray *mat,double *sig,double mean=0);
//void randomn(fltarray *mat,double *sig,double*mean);
void randomngsl(fltarray *mat,double *sig,double*mean);
double var_iter(double x,int n,double *mean,double*M);
void check_ineq(Ifloat dat,double thresh,Bool abs_en,Ifloat &flag);
void mr_support_filt(Ifloat img,Ifloat &img_filt, mr_opt opt,double nsig=7,int nb_iter=1,double*mse=NULL,Bool Pos_coeff=True,double lambda=1,Bool coarse_cons=False,Bool pos_cons=True, Bool drop_coarse=False,Bool iso_cons=False, double sig=89.108911);

inline void mod_denoise_usage(void)
{
    fprintf(OUTMAN, "         [-m type_of_model_denoising]\n");
    
    fprintf(OUTMAN, "              %d: %s \n",0,
                                            "Wavelet denoising of LR images before coaddition");
    fprintf(OUTMAN, "              %d: %s \n",1,
                                            "Iterative multiresolution support denoising after coaddition");
    fprintf(OUTMAN, "             default is %s.\n",  "Iterative multiresolution support denoising after coaddition");
}

inline void noise_est_usage(void)
{
    fprintf(OUTMAN, "         [-e Noise_est_meth]\n");
    
    fprintf(OUTMAN, "              %d: %s \n",0,
                                            "Estimation on backprojected residual using a MAD");
    fprintf(OUTMAN, "              %d: %s \n",1,
                                            "Noise simulation");
    //fprintf(OUTMAN, "              %d: %s \n",2,
    //                                      "Direct calculation of noise standard deviation per pixel (option under development)");
	    
    fprintf(OUTMAN, "             default is %s.\n",  "Estimation on backprojected residual using a MAD");
}

inline void thresh_type_usage(void)
{
    fprintf(OUTMAN, "         [-c Thresh_type]\n");
    
    fprintf(OUTMAN, "              %d: %s \n",0,
                                            "Hard Thresholding");
    fprintf(OUTMAN, "              %d: %s \n",1,
                                            "Soft Thresholding");
    //fprintf(OUTMAN, "              %d: %s \n",2,
    //                                      "Direct calculation of noise standard deviation per pixel (option under development)");
	    
    fprintf(OUTMAN, "             default is %s.\n",  "Soft Thresholding");
}

void decim_conv_mat(Ifloat conv_mat,int im_siz[],int decim_fact,Ifloat &output,double flux=1,double sig=1);
void convmask_shift_extend(Ifloat mask,int im_siz[],int shift[],Ifloat &output);
void noise_map_comp(double **shifts, double *sig, double *flux,int nb_im,int lancrad, int im_siz[],mr_opt opt, int decim_fact,fltarray &output);
int renewSeed();
#endif


