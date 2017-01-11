/******************************************************************************
**                   Copyright (C) 1997 by CEA
*******************************************************************************
**
**    UNIT
**
**    Version: 1.0
**
**    Author: Jean-Luc Starck
**
**    Date:  97/10/18 
**    
**    File:  Usage.h
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



#ifndef _USAGE_H_
#define _USAGE_H_

inline void vm_usage()
{
#ifdef LARGE_BUFF
   extern char emem_tmpdirname[1024];
   extern int emem_ramlimit;
    fprintf(OUTMAN, "         [-z]\n");
    fprintf(OUTMAN, "             Use virtual memory.\n");
    fprintf(OUTMAN, "                default limit size: %d\n",  emem_ramlimit);
    fprintf(OUTMAN, "                default directory: %s\n",  emem_tmpdirname); 
    manline();
    fprintf(OUTMAN, "         [-Z VMSize:VMDIR]\n");  
    fprintf(OUTMAN, "             Use virtual memory.\n");
    fprintf(OUTMAN, "                VMSize = limit size (megabytes) \n");
    fprintf(OUTMAN, "                VMDIR = directory name \n");
#endif
}

// ******************************
// Option in mr_transform
// ******************************

inline void nbr_nbr_undec_usage(int N=-1)
{
    fprintf(OUTMAN, "         [-u number_of_undecimated_scales]\n");
    fprintf(OUTMAN, "             Number of undecimated scales used in the Undecimated Wavelet Transform\n");
    if (N < 0) fprintf(OUTMAN, "             Default is all scale.\n");
    else fprintf(OUTMAN, "             Default is %d.\n", N);
}

inline void nbr_scale_usage(int Nbr_Plan)
{
    fprintf(OUTMAN, "         [-n number_of_scales]\n");
    fprintf(OUTMAN, "             Number of scales used in the multiresolution transform\n");
    fprintf(OUTMAN, "             Default is %d.\n", Nbr_Plan);
}

inline void write_scales_x_band_usage()
{
   fprintf(OUTMAN, "         [-x]\n");
   fprintf(OUTMAN, "             Write all bands separately as images with prefix 'band_j' (j being the band number)\n");
}

inline void write_scales_x_usage()
{
   fprintf(OUTMAN, "         [-x]\n");
   fprintf(OUTMAN, "             Write all scales separately as images with prefix 'scale_j' (j being the scale number)\n");
}

inline void write_band_usage()
{
   fprintf(OUTMAN, "         [-b BandNumber]\n");
   fprintf(OUTMAN, "             Extract a band.\n");
}

inline void read_band_usage()
{
   fprintf(OUTMAN, "         [-b BandNumber]\n");
   fprintf(OUTMAN, "             Insert a band.\n");
}

inline void write_scales_b_usage()
{
   fprintf(OUTMAN, "         [-B]\n");
   fprintf(OUTMAN, "             Same as x option, but interpolate by block the bands.\n");
}

inline void write_scales_i_usage()
{
    fprintf(OUTMAN, "         [-i]\n");
    fprintf(OUTMAN, "             Same as B option, but interpolate by a B3 spline the bands.\n");
    fprintf(OUTMAN, "             This option is valid only if the chosen multiresolution \n");
    fprintf(OUTMAN, "             transform is pyramidal (6,7,8,9,10,11,12). \n");
}

inline void iter_transform_usage()
{
    fprintf(OUTMAN, "         [-c iter]\n");
    fprintf(OUTMAN, "             Iterative transformation. Iter = number of iterations. \n");
    fprintf(OUTMAN, "             This option is valid only if the chosen multiresolution  \n");
    fprintf(OUTMAN, "             transform is pyramidal (6,7,8,9,10,11). The reconstruction \n");
    fprintf(OUTMAN, "             is not exact and we need few iterations. Generally, we take 3. \n");
}

// ******************************
// Option in mr_extract
// ******************************

inline void scale_number_usage()
{
   fprintf(OUTMAN, "         [-s scale_number]\n");
   fprintf(OUTMAN, "             Scale number to extract.\n");
}

// ******************************
// Option in mr_insert
// ******************************

inline void scale_number_insert_usage()
{
   fprintf(OUTMAN, "         [-s scale_number]\n");
   fprintf(OUTMAN, "             Scale number to insert.\n");
   fprintf(OUTMAN, "             By default, the first scale is used.\n");
}

// ******************************
// Option in mr_info
// ******************************

inline void analyse_struct_usage()
{
    fprintf(OUTMAN, "         [-a]\n");
    fprintf(OUTMAN, "              Significant structures analysis.\n");
    fprintf(OUTMAN, "              default is no.\n");
}

// ******************************
// Option in mr_filter
// ******************************


inline void nbr_scalep_usage(int DefNp)
{
    fprintf(OUTMAN, "             default is %d in case of poisson noise with few events.\n", DefNp);
}  

inline void nsigma_usage(float Sigma)
{ 
    fprintf(OUTMAN, "         [-s nsigma]\n");
    fprintf(OUTMAN, "             Thresholding at nsigma * SigmaNoise\n");
    fprintf(OUTMAN, "             default is %2.0f.\n", Sigma);
}

inline void gauss_usage()
{
    fprintf(OUTMAN, "         [-g sigma]\n");
    fprintf(OUTMAN, "             sigma = noise standard deviation\n");
    fprintf(OUTMAN, "             default is automatically estimated.\n");
}

inline void ccd_usage()
{
    fprintf(OUTMAN, "         [-c gain,sigma,mean]\n");
    fprintf(OUTMAN, "             Poisson + readout noise, with: \n");
    fprintf(OUTMAN, "                 gain = gain of the CCD\n");
    fprintf(OUTMAN, "                 sigma = read-out noise standard deviation\n");
    fprintf(OUTMAN, "                 mean = read-out noise mean\n");
    fprintf(OUTMAN, "             default is no (Gaussian).\n");
}

inline void max_iter_usage(int MaxIter)
{
    fprintf(OUTMAN, "         [-i number_of_iterations]\n");
    fprintf(OUTMAN, "             Maximum number of iterations\n");
    fprintf(OUTMAN, "             default is %d.\n", MaxIter);
}

inline void converg_param_usage(float Eps)
{
    fprintf(OUTMAN, "         [-e epsilon]\n");
    fprintf(OUTMAN, "             Convergence parameter\n");
    fprintf(OUTMAN, "             default is %f.\n",Eps);
}

inline void convergp_param_usage(float Eps)
{
    fprintf(OUTMAN, "             default is %f in case of poisson noise with few events.\n", Eps);
}

inline void support_file_usage()
{
    fprintf(OUTMAN, "         [-w support_file_name]\n");
    fprintf(OUTMAN, "             Creates an image from the multiresolution support \n");
    fprintf(OUTMAN, "             and save to disk.\n");
}

inline void kill_isol_pix_usage()
{
    fprintf(OUTMAN, "         [-k]\n");
    fprintf(OUTMAN, "             Suppress isolated pixels in the support. Default is no.\n");
}

inline void kill_last_scale_usage()
{
    fprintf(OUTMAN, "         [-K]\n");
    fprintf(OUTMAN, "             Suppress the last scale. Default is no.\n");
}

inline void detect_pos_usage()
{
    fprintf(OUTMAN, "         [-p]\n");
    fprintf(OUTMAN, "             Detect only positive structure. Default is no.\n");
 
}

inline void prec_eps_poisson_usage(float Eps)
{
    fprintf(OUTMAN, "         [-E Epsilon]\n");
    fprintf(OUTMAN, "             Epsilon = precision for computing thresholds\n");
    fprintf(OUTMAN, "                       (only used in case of poisson noise with few events)\n");
    fprintf(OUTMAN, "             default is %5.2e \n", Eps);
}

inline void size_block_usage(int SizeBlock)
{
    fprintf(OUTMAN, "         [-S SizeBlock]\n");
    fprintf(OUTMAN, "             Size of the  blocks used for local variance estimation.\n");
    fprintf(OUTMAN, "             default is %d.\n", SizeBlock);
}

inline void sigma_clip_block_usage(int NiterClip)
{
    fprintf(OUTMAN, "         [-N NiterSigmaClip]\n");
    fprintf(OUTMAN, "             Iteration number used for local variance estimation.\n");
    fprintf(OUTMAN, "             default is %d.\n", NiterClip);
}

inline void first_detect_scale_usage()
{
    fprintf(OUTMAN, "         [-F first_detection_scale]\n");
    fprintf(OUTMAN, "             First scale used for the detection \n");
    fprintf(OUTMAN, "             default is 1.\n");
}

inline void window_size_usage(int SWindowSize)
{
    fprintf(OUTMAN, "         [-W WindowSize]\n");
    fprintf(OUTMAN, "             Window size for median and average filtering.\n");
    fprintf(OUTMAN, "             default is %d.\n", SWindowSize);
}


// ******************************
// Option in mr_deconv
// ******************************


inline void poisson_noise_usage()
{
    fprintf(OUTMAN, "         [-p]\n");
    fprintf(OUTMAN, "             Poisson Noise\n");
    fprintf(OUTMAN, "             default is no (Gaussian).\n");
}

inline void dilate_sup_usage()
{
    fprintf(OUTMAN, "         [-l]\n");
    fprintf(OUTMAN, "             Dilate the support\n");
}

inline void write_residual_usage()
{
    fprintf(OUTMAN, "         [-r residual_file_name]\n");
    fprintf(OUTMAN, "             Residual_file_name = file name\n");
    fprintf(OUTMAN, "             write the residual to the disk \n");
}

inline void fwhm_usage(float Fwhm)
{
    fprintf(OUTMAN, "         [-f Fwhm]\n");
    fprintf(OUTMAN, "             Full width at half maximum.\n");
    fprintf(OUTMAN, "             Default value is %f\n", Fwhm);
}

inline void gain_clean_usage(float Gain)
{
    fprintf(OUTMAN, "         [-G gamma_parameter]\n");
    fprintf(OUTMAN, "             gamma parameter. Only used by CLEAN method.\n"); 
    fprintf(OUTMAN, "             Default value is %f\n", Gain); 
}

inline void psf_not_center_usage()
{
   fprintf(OUTMAN, "         [-S]\n");
   fprintf(OUTMAN, "             Do not shift automatically the maximum  \n");
   fprintf(OUTMAN, "             of the PSF at the center.\n");
}

// ******************************
// Option in mr_psupport
// ******************************

inline void input_poisson_usage()
{
    fprintf(OUTMAN, "         [-a ascii_file]\n");
    manline();
    fprintf(OUTMAN, "         [-I image_file]\n");  
    fprintf(OUTMAN, "         a & I options can't be used together, \n");
    fprintf(OUTMAN, "         and one must be set. \n");
}

inline void min_event_usage (int MinEvent)
{
    fprintf(OUTMAN, "         [-e minimum_of_events]\n");
    fprintf(OUTMAN, "             Minimum number of events for a detection.\n");
    fprintf(OUTMAN, "             default is %d\n", MinEvent);
}

inline void write_wave_mr_usage()
{
    fprintf(OUTMAN, "         [-w]\n");
    fprintf(OUTMAN, "             Write the following file:\n");
    fprintf(OUTMAN, "              xx_Wavelet.mr : contains the wavelet transform\n");
    fprintf(OUTMAN, "              of the image.\n");
}

inline void signif_ana_usage()
{
    fprintf(OUTMAN, "         [-s SignifStructureAnalysis_FileName]\n");
    fprintf(OUTMAN, "             Write in xx_Segment.mr the segmented scales.\n");
    fprintf(OUTMAN, "             Analyse the detected wavelet coefficients,\n");
    fprintf(OUTMAN, "             and write in the file:\n");
    fprintf(OUTMAN, "               Number of detected structures per scale\n");
    fprintf(OUTMAN, "               Percentage of significant wavelet coefficents\n");
    fprintf(OUTMAN, "               Mean deviation of shape from sphericity\n");
    fprintf(OUTMAN, "               For each detected structure, its surface aera, its perimeter, and\n");
    fprintf(OUTMAN, "               its deviation of shape from sphericity, \n");
    fprintf(OUTMAN, "               its angle, its elongation in both axis directions.\n");
}

inline void ascii_signif_ana_usage()
{
    fprintf(OUTMAN, "         [-t SignifStructureAnalysis_FileName]\n");
    fprintf(OUTMAN, "             Same as -s option, but results are stored\n");
    fprintf(OUTMAN, "             in an ascii table format.\n");
    fprintf(OUTMAN, "             The table contains: scale number, structure number, \n");
    fprintf(OUTMAN, "             Max_x, Max_y, Surface, Perimeter, Morpho, \n");
    fprintf(OUTMAN, "             Angle, Sigma_X, Sigma_Y. \n\n"); 
}

inline void abaque_file_usage()
{
    fprintf(OUTMAN, "         [-q abaque_file]\n");
    fprintf(OUTMAN, "              default is Abaque.fits.\n\n");
}
// ******************************
// Option in mr_abaque
// ******************************

inline void abaque_option_usage(char *Name_Abaque_Default)
{
    fprintf(OUTMAN, "         [-n Number]\n");
    fprintf(OUTMAN, "             Number = Number of scales as a power of 2\n");
    fprintf(OUTMAN, "             default is 25\n");
manline();

    fprintf(OUTMAN, "         [-w]\n");
    fprintf(OUTMAN, "             Write the following files:\n");
    fprintf(OUTMAN, "               Aba_histo.fits: contains all histograms\n");
    fprintf(OUTMAN, "                    h(3*i)   = histogram values\n");
    fprintf(OUTMAN, "                    h(3*i+1) = reduced coordinates\n");
    fprintf(OUTMAN, "                    h(3*i+2) = normalized histogram\n");
//     fprintf(OUTMAN, "               Aba_distrib.fits: distribution functions\n");
//     fprintf(OUTMAN, "                    F(3*i) = function values\n");
//     fprintf(OUTMAN, "                    F(3*i+1) = reduced coordinates\n");
//     fprintf(OUTMAN, "                    F(3*i+2) = reduced values\n");
//    fprintf(OUTMAN, "               Aba_log_distrib.fits: log transformation of F\n");
//    fprintf(OUTMAN, "                    L(3*i) = function values\n");
//    fprintf(OUTMAN, "                    L(3*i+1) = real coordinates\n");
//    fprintf(OUTMAN, "                    L(3*i+2) = reduced coordinates\n");
//     fprintf(OUTMAN, "               Aba_mean.fits: contains the mean real values of the histograms\n");
//     fprintf(OUTMAN, "               Aba_sigma.fits: contains the sigma real values of the histograms\n");
    fprintf(OUTMAN, "               Aba_bspline.fits: contains the used Bspline\n");
    fprintf(OUTMAN, "               Aba_wavelet.fits: contains the used wavelet\n");
manline();

    fprintf(OUTMAN, "          [-d]\n");
    fprintf(OUTMAN, "             Use all default parameters\n");
    fprintf(OUTMAN, "                default Number of scales\n\n");
    fprintf(OUTMAN, "                default precision\n\n");
    fprintf(OUTMAN, "                default abaque file name is %s\n\n", Name_Abaque_Default);
}

// ******************************
// Option in mr_pfilter
// ******************************

inline void write_pfilter_usage()
{
    fprintf(OUTMAN, "         [-w]\n\n");
    fprintf(OUTMAN, "           write the following files\n");
    fprintf(OUTMAN, "             xx_Wavelet.mr : contains the wavelet transform\n");
    fprintf(OUTMAN, "              of the image.\n");

    fprintf(OUTMAN, "             xx_Support.mr : contains the thresholded  wavelet transform.\n\n");
}

// ******************************
// Option in mr_sigma
// ******************************

inline void sigma_gain()
{
  fprintf(OUTMAN, "\n");
  fprintf(OUTMAN, "         [-p gain]\n");
  fprintf(OUTMAN, "              performs the standard deviation of the \n");
  fprintf(OUTMAN, "              Gaussian part of the noise.\n");
  fprintf(OUTMAN, "              gain used in the model of image \n");
  fprintf(OUTMAN, "              Input image must of the form :\n");
  fprintf(OUTMAN, "              Gain * Poisson_Noise + Zero_Mean_Gaussian_Noise :\n");
  fprintf(OUTMAN, "              default is no\n");
}

// ******************************
// Option in mr_fusion
// ******************************

inline void fusion_option_usage(int ResMin)
{
    fprintf(OUTMAN, "         [-r res_min]\n");
    fprintf(OUTMAN, "             Miminum resolution for reconstruction\n");
    fprintf(OUTMAN, "             default is %d\n", ResMin);
manline();
    fprintf(OUTMAN, "         [-D dist_max]\n");
    fprintf(OUTMAN, "             Maximum estimated distance between \n");
    fprintf(OUTMAN, "             two identical points in both images.\n");
manline();
//     fprintf(OUTMAN, "         [-l]\n");
//     fprintf(OUTMAN, "              Registration choice: \n");
//     fprintf(OUTMAN, "                0: Sub-scene registration \n");
//     fprintf(OUTMAN, "                1: Sub-scene and scene registration \n");
//     fprintf(OUTMAN, "              Default is Sub-scene registration.\n");
// manline();

    fprintf(OUTMAN, "         [-o]\n");
    fprintf(OUTMAN, "              Manual Options specifications:\n");
    fprintf(OUTMAN, "                 -Matching distance\n");
    fprintf(OUTMAN, "                 -Threshold level\n");
    fprintf(OUTMAN, "                 -Type of deformation model\n");
    fprintf(OUTMAN, "                 -Type of interpolation \n");
    fprintf(OUTMAN, "                 0: None\n");
    fprintf(OUTMAN, "                 1:  Matching distance\n");
    fprintf(OUTMAN, "                 2:  Threshold level\n");
    fprintf(OUTMAN, "                 3:  Type of deformation model\n");
    fprintf(OUTMAN, "                 4:  Matching distance\n");
    fprintf(OUTMAN, "                     Threshold level\n");
    fprintf(OUTMAN, "                     Type of deformation model\n");
    fprintf(OUTMAN, "                 5:  Matching distance\n");
    fprintf(OUTMAN, "                     Threshold level\n");
    fprintf(OUTMAN, "                     Type of deformation model\n");
    fprintf(OUTMAN, "                     Type of interpolation\n");
    fprintf(OUTMAN, "              Default is None. \n");
manline();

    fprintf(OUTMAN, "         [-i InterpolType]\n");
    fprintf(OUTMAN, "                Type of interpolation: \n");
    fprintf(OUTMAN, "                 0: zero order interpolation - nearest neighour\n");
    fprintf(OUTMAN, "                 1: First order interpolation - bilinear \n");
    fprintf(OUTMAN, "                 2: Second order interpolation - cubic \n");
    fprintf(OUTMAN, "              Default is 2. \n");   
manline();
    
    fprintf(OUTMAN, "         [-d DeforModel]\n");
    fprintf(OUTMAN, "                Deformation model: \n");
    fprintf(OUTMAN, "                 0: Polynomial of the first order of type I.\n");
    fprintf(OUTMAN, "                 1: Polynomial of the first order of type II. \n");
    fprintf(OUTMAN, "                 2: Polynomial of the second order. \n");
    fprintf(OUTMAN, "                 3: Polynomial of the third order. \n");
    fprintf(OUTMAN, "              Default is 1. \n");   
manline();
    fprintf(OUTMAN, "         [-w]\n");
    fprintf(OUTMAN, "              write the following files:\n");
    fprintf(OUTMAN, "              - Deformation model parameters for each scale\n");
    fprintf(OUTMAN, "              - control points for each scale\n");
    fprintf(OUTMAN, "              default in None.\n");
}

// ******************************
// Option in mr_visu
// ******************************

inline void mrvisu_option_usage(float NSigma)
{
    fprintf(OUTMAN, "         [-V Type_Visu]\n");
    fprintf(OUTMAN, "              1: Gray level\n");
    fprintf(OUTMAN, "              2: Contour\n");
    fprintf(OUTMAN, "              3: Perspective\n");
    fprintf(OUTMAN, "              Default is 1.\n");
 manline();   
    fprintf(OUTMAN, "         [-b]\n");
    fprintf(OUTMAN, "             Save output image in bi-level.\n");
    fprintf(OUTMAN, "             Only used if Type_Visu equal to 2 or 3.\n");
 manline();
    fprintf(OUTMAN, "         [-c]\n");
    fprintf(OUTMAN, "             Do not apply a normalization on the multiresolution coefficient \n");
    fprintf(OUTMAN, "             Only used if Type_Visu equal to 1.\n");
 manline();
 
//     fprintf(OUTMAN, "         [-w PS_FileName]\n");
//     fprintf(OUTMAN, "             Save also the result in a postscript file\n");
// manline();

    fprintf(OUTMAN, "         [-i Increment]\n");
    fprintf(OUTMAN, "             Number of lines of the image which will be used.\n");
    fprintf(OUTMAN, "             If Increment = 3, only on line in 3 is used.\n");
    fprintf(OUTMAN, "             Only used if Type_Visu equal to 3.\n");
    fprintf(OUTMAN, "             The default value is 1.\n");
manline();
    fprintf(OUTMAN, "         [-s nsigma]\n");    
    fprintf(OUTMAN, "             Plot contour at nsigma*Sigma if Type_Visu equal to 2.\n");
    fprintf(OUTMAN, "             Threshold value upper nsigma*Sigma if Type_Visu equal to 3.\n");
    fprintf(OUTMAN, "             default is %f.\n", NSigma);
}

// ******************************
// Option in mr_detect
// ******************************
inline void last_detect_scale_usage()
{
    fprintf(OUTMAN, "         [-L last_detection_scale]\n");
    fprintf(OUTMAN, "             Last scale used for the detection.\n");
}

inline void verbose_usage()
{   
    fprintf(OUTMAN, "         [-v]\n");
    fprintf(OUTMAN, "             Verbose. Default is no.\n");  
}

inline void mrdetect_option_usage(int Nb_iter_rec, float Eps_ErrorRec)
{
    fprintf(OUTMAN, "         [-i number_of_iterations]\n");
    fprintf(OUTMAN, "             Iteration number per object reconstruction\n");
    fprintf(OUTMAN, "             default is %d\n", Nb_iter_rec);
manline();
    fprintf(OUTMAN, "         [-u object_reconstruction_error]\n");
    fprintf(OUTMAN, "             default is: %f\n", Eps_ErrorRec);
manline();
    fprintf(OUTMAN, "         [-k]\n");
    fprintf(OUTMAN, "             Keep isolated objects\n");
    fprintf(OUTMAN, "             default is no.\n");
manline();
    fprintf(OUTMAN, "         [-K]\n");
    fprintf(OUTMAN, "              Keep objects at the border.\n");
    fprintf(OUTMAN, "              default is no.\n");
manline();
    fprintf(OUTMAN, "         [-A FluxMult]\n");
    fprintf(OUTMAN, "              Flux in tex table are multiplied by FluxMul.\n");
    fprintf(OUTMAN, "              default is 1.\n");    
manline();
    fprintf(OUTMAN, "         [-w writing_parameter]\n");
    fprintf(OUTMAN, "              1: write each object separately in an image. \n");
    fprintf(OUTMAN, "                 The image file name of the object will be: \n");
    fprintf(OUTMAN, "                      ima_obj_xx_yy.fits \n");
    fprintf(OUTMAN, "              2: simulated two images \n");
    fprintf(OUTMAN, "                   xx_ellips.fits: an ellipse is drawn around each object \n");
    fprintf(OUTMAN, "                   xx_simu.fits: image created only from the morphological parameters \n");
    fprintf(OUTMAN, "              3: equivalent to 1 and 2 together \n");
manline();
    fprintf(OUTMAN, "         [-U]\n");
    fprintf(OUTMAN, "             Sub Segmentation.\n");   
manline();
    fprintf(OUTMAN, "         [-p]\n");
    fprintf(OUTMAN, "             Detect also negative structures \n");
    fprintf(OUTMAN, "             default is no.\n");
manline();
    fprintf(OUTMAN, "         [-q]\n");
    fprintf(OUTMAN, "              Define the root of an object from the maximum position and its value \n");
    manline();
    
    fprintf(OUTMAN, "         [-d DistMax]\n");
    fprintf(OUTMAN, "              Maximum distance between two max positions\n");
    fprintf(OUTMAN, "              of the same object at two successive scales.\n");
    fprintf(OUTMAN, "              Default is 1.\n");
    manline();  
}

// ******************************
// Option in mr_comp
// ******************************

inline void mrcomp_option_usage(int Elstr_Size, float SignalQuantif, 
                               float NoiseQuantif, int WindowSize)
{
    fprintf(OUTMAN, "         [-r]\n");
    fprintf(OUTMAN, "              Compress the noise. Default is no.\n");
    manline();    

    fprintf(OUTMAN, "         [-k]\n");
    fprintf(OUTMAN, "              Keep isolated pixel in the support \n");
    fprintf(OUTMAN, "              at the first scale. Default is no.\n");
    fprintf(OUTMAN, "              If the PSF is on only one pixel, this\n");
    fprintf(OUTMAN, "              option should be set\n");
    manline();    

    fprintf(OUTMAN, "         [-l]\n");
    fprintf(OUTMAN, "              Save the noise (for lossless compression)\n");
    fprintf(OUTMAN, "              default is no\n");
    manline();    

    fprintf(OUTMAN, "         [-q signal_quantif]\n");
    fprintf(OUTMAN, "              Signal quantification\n");
    fprintf(OUTMAN, "              default is %5.2f\n", SignalQuantif);
    manline();    

    fprintf(OUTMAN, "         [-e noise_quantif]\n");
    fprintf(OUTMAN, "              Noise quantification\n");
    fprintf(OUTMAN, "              default is %5.2f.\n", NoiseQuantif);
    manline();    

    fprintf(OUTMAN, "         [-f ]\n");
    fprintf(OUTMAN, "              Keep all the fits header. Default is no.\n");
    manline();    

    fprintf(OUTMAN, "         [-b] bad_pixel_value\n");
    fprintf(OUTMAN, "              all pixels with this value will be\n");
    fprintf(OUTMAN, "              considered as bad pixels, and not\n");
    fprintf(OUTMAN, "              used for the noise modeling.\n");
    manline();    

    fprintf(OUTMAN, "         [-O]\n");
    fprintf(OUTMAN, "              optimization. If set, the program\n");
    fprintf(OUTMAN, "              works with integer instead of float.\n");
    manline();    

    fprintf(OUTMAN, "         [-B]\n");
    fprintf(OUTMAN, "              optimization without BSCALE operation\n");
    fprintf(OUTMAN, "              in case of fits images.\n");
    manline();    

    fprintf(OUTMAN, "         [-P]\n");
    fprintf(OUTMAN, "              Keep only positive coefficients.\n");
    manline();    

    fprintf(OUTMAN, "         [-W]\n");
    fprintf(OUTMAN, "              Median window size equal to 3. Default is %d\n", WindowSize);
    manline();    

    fprintf(OUTMAN, "         [-S]\n");
    fprintf(OUTMAN, "              Use a square structural element.\n");
    fprintf(OUTMAN, "              (Only for math. morphology compresssion.\n");
    fprintf(OUTMAN, "               Default structural element is a circle\n");
     manline();    
   
    fprintf(OUTMAN, "         [-D Dim]\n");
    fprintf(OUTMAN, "             Dimension of the structural element.\n");
    fprintf(OUTMAN, "             (Only for math. morphology compresssion.\n");
    fprintf(OUTMAN, "              Default is %d.\n", Elstr_Size);
     manline();    
   
    fprintf(OUTMAN, "         [-R Compression_Ratio]\n");
    fprintf(OUTMAN, "             Fixes the compression ratio.\n");    
    fprintf(OUTMAN, "              Default is no.\n");   
       manline();    
 
    fprintf(OUTMAN, "         [-C BlockSize]\n");
    fprintf(OUTMAN, "              Compress by block. \n");
    fprintf(OUTMAN, "              BlockSize = size of each block.\n");
    fprintf(OUTMAN, "              Default is No.\n");
    manline();  
      
    fprintf(OUTMAN, "       [-i NbrIter]\n");
    fprintf(OUTMAN, "              Apply an iterative compression.\n");
    fprintf(OUTMAN, "              NbrIter = Number of iterations. \n");
    fprintf(OUTMAN, "              Only used with orthogonal transform.\n");
    fprintf(OUTMAN, "              Default is no iteration.\n");
    manline();
    
    fprintf(OUTMAN, "         [-N]\n");
    fprintf(OUTMAN, "             Do not use noise modeling.\n");
}

// ******************************
// Option in mr_decomp
// ******************************

inline void mrdecomptool_option_usage()
{
    fprintf(OUTMAN, "        [-B BlockNbr]\n");
    fprintf(OUTMAN, "              Decompress only one block. \n");
    fprintf(OUTMAN, "              BlockNbr is the block number to decompress.\n");
    fprintf(OUTMAN, "              Default is no.\n");
}

// **********************

inline void mrdecomp_option_usage()
{
    mrdecomptool_option_usage();
    manline();    
  
    fprintf(OUTMAN, "        [-r resolution]\n");
    fprintf(OUTMAN, "          resol = 0..nbr_of_scale-1 \n");
    fprintf(OUTMAN, "          resol = 0 for full resolution (default) \n");
    fprintf(OUTMAN, "          resol = nbr_of_scale-1 for the worse resol.\n");
  
    manline();    
  
    fprintf(OUTMAN, "        [-t] output type\n");
    fprintf(OUTMAN, "              if the input image was a fits image, \n");
    fprintf(OUTMAN, "              the image output type can be fixed by the user \n");
    fprintf(OUTMAN, "              to 'f' for float, to 'i' for integer, or 's' \n");
    fprintf(OUTMAN, "              for short. By default, the output type \n");
    fprintf(OUTMAN, "              is the same as the type of the original image. \n");
       manline();    

    fprintf(OUTMAN, "        [-g] \n");
    fprintf(OUTMAN, "              add a simulated noise to the decompressed\n");
    fprintf(OUTMAN, "              image with the same properties as in the\n");
    fprintf(OUTMAN, "              original image. So they look very similar.\n");
       manline();    
    
    fprintf(OUTMAN, "        [-I IterRecNbr]\n");
    fprintf(OUTMAN, "              Use an iterative reconstruction. \n");
    fprintf(OUTMAN, "              Only used with orthogonal transform.\n");
    fprintf(OUTMAN, "              Default is no iteration.\n");
      
}
 
// ******************************
// Option in mr_background
// ******************************

inline void mrbgr_option_usage(int Npix)
{
    fprintf(OUTMAN, "         [-n number_of_pixels]\n");
    fprintf(OUTMAN, "             Number of pixels of the last scale.\n");
    fprintf(OUTMAN, "             Default is %d.\n", Npix);
       manline();    

    fprintf(OUTMAN, "         [-w background_file_name]\n");
    fprintf(OUTMAN, "             backgroung_file_name = file name\n");
    fprintf(OUTMAN, "             creates the backgroung   \n");
    fprintf(OUTMAN, "             and write it on the disk\n");
}

// ******************************
// Option in mr1d_detect
// ******************************

inline void mr1ddetectr_option_usage(int RecIterNumber)
{
    fprintf(OUTMAN, "         [-a]\n");
    fprintf(OUTMAN, "              detection of Absorption lines. Default is no. \n");
       manline();    

    fprintf(OUTMAN, "         [-e]\n");
    fprintf(OUTMAN, "              detection of Emission lines. Default is no. \n");
       manline();    

    fprintf(OUTMAN, "         [-f FirstScale]\n");
    fprintf(OUTMAN, "             first scale. Default is 1.\n\n");
       manline();    

    fprintf(OUTMAN, "         [-l LastScale]\n");
    fprintf(OUTMAN, "             Last scale. Default is number_of_scales-2.\n");
       manline();    

    fprintf(OUTMAN, "         [-i IterNumber]\n");
    fprintf(OUTMAN, "             Number of iteration for the reconstruction. \n");
    fprintf(OUTMAN, "             Default is %d\n", RecIterNumber);
       manline();    

    fprintf(OUTMAN, "         [-M]\n");
    fprintf(OUTMAN, "             Use the multiresolution median transform \n");
    fprintf(OUTMAN, "             instead of the a-trous algorithm. \n");
       manline();    

    fprintf(OUTMAN, "         [-A]\n");
    fprintf(OUTMAN, "              detect only negative  multiresolution coefficients. Default is no. \n");
       manline();    

    fprintf(OUTMAN, "         [-E]\n");
    fprintf(OUTMAN, "              detect only positive multiresolution coefficients. Default is no. \n");
       manline();    

    fprintf(OUTMAN, "         [-w ]\n");
    fprintf(OUTMAN, "              write other results:\n");
    fprintf(OUTMAN, "                tabadd.fits = sum of the reconstructed objects  \n");
    fprintf(OUTMAN, "                tabseg.fits = segmented wavelet transform  \n");
}
// ******************************
// Option in im_simu
// ******************************

inline void imsimu_option_usage()
{
    fprintf(OUTMAN, "         [-p ]\n");
    fprintf(OUTMAN, "             Poisson Noise. Default is no. \n\n");

    fprintf(OUTMAN, "         [-g sigma]\n");
    fprintf(OUTMAN, "             sigma = noise standard deviation\n");
    fprintf(OUTMAN, "             default is no. \n\n");

    fprintf(OUTMAN, "         [-c sigma]\n");
    fprintf(OUTMAN, "             Poisson Noise + gaussian noise\n");
    fprintf(OUTMAN, "             sigma = gaussian noise standard deviation\n");
    fprintf(OUTMAN, "             default is no.\n\n");

    fprintf(OUTMAN, "         [-r psf_image]\n");
    fprintf(OUTMAN, "             psf_image = instrumental response (PSF)\n");
    fprintf(OUTMAN, "             default is no. \n\n");

    fprintf(OUTMAN, "         [-f width]\n");
    fprintf(OUTMAN, "             width = full width at half-maximum of the\n");
    fprintf(OUTMAN, "                     gaussian instrumental response (FWHM)\n");
    fprintf(OUTMAN, "             default is no.\n\n");

    fprintf(OUTMAN, "         [-I InitRandomVal]\n");
    fprintf(OUTMAN, "             Value used for random value generator initialization.\n");
    fprintf(OUTMAN, "             default is 100. \n\n");
    
    fprintf(OUTMAN, "         [-w psf_file_name]\n");
    fprintf(OUTMAN, "             write the PSF to the disk. Valid only if -f is set.\n");
    fprintf(OUTMAN, "             default is no.\n");
    

}
// ******************************
// Option in im_segment
// ******************************

inline void imsegment_option_usage()
{
    fprintf(OUTMAN, "         [-b]\n");
    fprintf(OUTMAN, "             Eliminates regions at the border.\n");
    fprintf(OUTMAN, "             default is no.\n");
}

// ******************************
// Option in im_opening
// ******************************

inline void imopen_usage()
{
    fprintf(OUTMAN, "         [-n opening_number]\n");
    fprintf(OUTMAN, "             opening number.\n");
    fprintf(OUTMAN, "             default is 1.\n");
}

inline void immorpho_usage(int Elstr_Size)
{
    fprintf(OUTMAN, "         [-s structural_element]\n");
    fprintf(OUTMAN, "             1 => sqare \n");
    fprintf(OUTMAN, "             2 => cross \n");
    fprintf(OUTMAN, "             3 => circle \n");
    fprintf(OUTMAN, "             default is 3.\n\n");
    
    fprintf(OUTMAN, "         [-d Dim]\n");
    fprintf(OUTMAN, "              Dimension of the structural element.\n");
    fprintf(OUTMAN, "              Only for square and circle.\n");
    fprintf(OUTMAN, "              Default is %d\n", Elstr_Size);
}

// ******************************
// Option in im_erode
// ******************************

inline void imerode_usage()
{
   fprintf(OUTMAN, "         [-n erosion_number]\n");
   fprintf(OUTMAN, "             erosion number.\n");
   fprintf(OUTMAN, "             default is 1.\n\n");
} 

// ******************************
// Option in im_dilate
// ******************************

inline void imdilate_usage()
{
   fprintf(OUTMAN, "         [-n dilation_number]\n");
   fprintf(OUTMAN, "             dilation number.\n");
   fprintf(OUTMAN, "             default is 1.\n");  
} 

// ******************************
// Option in im_closing
// ******************************

inline void imclosing_usage()
{
    fprintf(OUTMAN, "         [-n closing_number]\n");
    fprintf(OUTMAN, "             closing number.\n");
    fprintf(OUTMAN, "             default is 1.\n"); 
} 
#endif


