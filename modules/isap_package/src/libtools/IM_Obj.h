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
**    File:  IM_Obj.h
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

#ifndef _IMAGE_H_
#define _IMAGE_H_

#include "GlobalInc.h"
#include "Border.h"
#include <cmath>

#include "TempArray.h"

enum details { D_NULL, D_HORIZONTAL, D_DIAGONAL, D_VERTICAL, I_SMOOTH, D_HALF_RESOL, D_RESOL};

enum typtrans { DATA, SIGMA};

void img_error (char *);

int test_nl (int N); 
int test_nc (int N); 

Bool test_same_size (char *Im1, int Nl1, int Nc1, char *Im2, int Nl2, int Nc2);
Bool test_indice_ij (const char *Name_Img, int i, int j, int Nl, int Nc);


#define DEFAULT_NL 256
#define DEFAULT_NC 256
#define DEFAULT_MED_WINSIZE 3

//#define NOT_USED 1e99 
#define NOT_USED 10000 

const Iint iint (const Ifloat &);
const Iint iint (const Icomplex_f &);
//const Iint operator + (const Iint &, const Iint &);
//const Iint operator - (const Iint &, const Iint &);
//const Iint operator * (const Iint &, const Iint &);
//const Iint operator / (const Iint &, const Iint &);

double sigma (const Iint &);
double average (const Iint &);
int min (const Iint &);
int max (const Iint &);
int flux (const Iint &);
int energy (const Iint &);

const Icomplex_f fft(const Iint &);
const Icomplex_f invfft(const Iint &);
void fft2d (const Iint &, const Icomplex_f &, int Dir = 1, bool normalize=false);
const Iint conv_fft (  const Iint &, const Iint &);
//const Iint conv_direct (const Iint &, const Iint &, type_border Type_Border = I_ZERO);

float correlation (const Iint &, const Iint &);

void smooth_average (const Iint & Im_in, Iint &Im_out, 
         type_border Type = DEFAULT_BORDER, int Step_trou=0, int Size_Window=3);
void smooth_linear (const Iint & Im_in, Iint &Im_out,
         type_border Type = DEFAULT_BORDER, int Step_trou=0);
void smooth_bspline (const Iint & Im_in, Iint &Im_out,
         type_border Type = DEFAULT_BORDER, int Step_trou=0);
void smooth_mediane (const Iint &Imag1, Iint &Imag2, 
         type_border Type=DEFAULT_BORDER, int Step_trou = 0,
         int Window_Size=DEFAULT_MED_WINSIZE);
void fmedian(int *in_image, int Nl, int Nc, int *out_image,
            type_border Border=DEFAULT_BORDER,
            int Window_Size=DEFAULT_MED_WINSIZE);
void fmedian(Iint & in_image, Iint & out_image,
              type_border Border=DEFAULT_BORDER,
              int Window_Size=DEFAULT_MED_WINSIZE);

void im_reduce_size_2 (const Iint &Imag, Iint &Imag_Out);
void im_increase_size_2 (const Iint &Imag, Iint &Pict_out, 
         type_border Border=DEFAULT_BORDER);
void im_extend (const Iint &Imag, Iint &Imag_Out);

void morphoi_erosion (int *Imag1, int *Imag2, int Nl, int Nc, int Window_Size = 3);
void morphoi_dilation (int *Imag1, int *Imag2, int Nl, int Nc, int Window_Size = 3);
void morphoi_cercle_erosion (int *Imag1, int *Imag2, int Nl, int Nc, int Window_Size=5);
void morphoi_cercle_dilation (int *Imag1, int *Imag2, int Nl, int Nc, int Window_Size=5);
 
//void morphoi_dilation_class (Iint &Imag1, Iint &Imag2, int Window_Size=3);
//void cerclei_erosion_class (Iint &Imag1, Iint &Imag2, int Window_Size=5);
//void cerclei_dilation_class (Iint &Imag1, Iint &Imag2, int Window_Size=5);

#define MIN_SIZE_FOR_MEM_ALLOC_CALL 50000
// we use the memory class alocator only for image
// size greater than MIN_SIZE_FOR_MEM_ALLOC_CALL

inline void INFO (Ifloat &Data, char *Mes)  
{ 
   int i,j;
   float Min=Data(0,0),Max=Data(0,0);
   double Flux=0., Mean, Sigma=0.;

   for (i = 0; i < Data.nl(); i++) 
   for (j = 0; j < Data.nc(); j++)
   {
      if (Min > Data(i,j)) Min = Data(i,j);
      if (Max < Data(i,j)) Max = Data(i,j);
      Flux += Data(i,j);
      Sigma += Data(i,j)*Data(i,j);
   }
   Mean = Flux / (Data.nl()*Data.nc());
   Sigma = sqrt(Sigma/(Data.nl()*Data.nc()) - Mean*Mean);
   cout << Mes << ": Min = " << Min << "    Max = " << Max << "    Flux = " << Flux;
   cout << "    Mean = " << Mean << "     Sigma = " << Sigma << endl;
}

inline Bool isolated_pixel(Ifloat &Ima, int i, int j, type_border Bord=I_MIRROR)
{
   Bool ValRet = False;
   
   if ((ABS(Ima(i,j)) > FLOAT_EPSILON) && ( (ABS(Ima(i,j-1,Bord)) < FLOAT_EPSILON)
         && (ABS(Ima(i-1,j-1,Bord)) < FLOAT_EPSILON) && (ABS(Ima(i+1,j-1,Bord)) < FLOAT_EPSILON) 
	 && (ABS(Ima(i-1,j,Bord))  < FLOAT_EPSILON) &&  (ABS(Ima(i+1,j,Bord))  < FLOAT_EPSILON) 
	 &&  (ABS(Ima(i,j+1,Bord))  < FLOAT_EPSILON) && (ABS(Ima(i-1,j+1,Bord))  < FLOAT_EPSILON) 
	 && (ABS(Ima(i+1,j+1,Bord))  < FLOAT_EPSILON))) ValRet = True;
   return ValRet;
}



double sigma (const Ifloat &);
double average (const Ifloat &);
float min (const Ifloat &);
float max (const Ifloat &);
float flux (const Ifloat &);
double energy (const Ifloat &);
double total (const Ifloat &);
void norm_flux (Ifloat &Image, float Val_Norm=1.);
float total_in_block(Ifloat & Image, int i, int j, int Neighbour);

void fft2d (const Ifloat &, const Icomplex_f &, int Dir = 1, bool normalize=false);
// const Icomplex_f fft(const Ifloat &);
// const Icomplex_f invfft(const Ifloat &);

const Ifloat conv_fft (const Ifloat &, const Ifloat &);
//const Ifloat conv_direct (const Ifloat &, const Ifloat &, type_border Type_Border = I_ZERO);
void fft2d_conv(const Ifloat &Im_in1, const Ifloat &Im_in2, Ifloat &Im_out);

float correlation (const Ifloat &, const Ifloat &);
float correlation (const Ifloat & I1 , const Ifloat & I2, float ValMin);
float regression (const Ifloat &Im1, const Ifloat & Im2);
float regression (const Ifloat &Im1, const Ifloat & Im2, float ValMin);
// Search the regression coefficient such that:
//   Im1 = Im2 * CoefRegression
// takes into account only value with absolute value > ValMin

void ifloat (Ifloat& Flt, Iint& Intt) ;

//ifloat (const Icomplex_f &);

void threshold (Ifloat &Image, float T=0.);

void smooth_average (const Ifloat & Im_in, Ifloat &Im_out, 
              type_border Type = DEFAULT_BORDER, 
              int Step_trou=0, int Size_Window=3);
void smooth_linear (const Ifloat & Im_in, Ifloat &Im_out,
                                 type_border Type = DEFAULT_BORDER, 
                                 int Step_trou=0);
void smooth_bspline (const Ifloat & Im_in, Ifloat &Im_out,
                                 type_border Type = DEFAULT_BORDER,
                                 int Step_trou=0);
void smooth_haar (const Ifloat & Im_in, Ifloat &Im_out,
                                   type_border Type= DEFAULT_BORDER, 
				   int Step_trou=0, float Norm=4.);

void smooth_mediane (const Ifloat &Imag1, Ifloat &Imag2, 
             type_border Type = DEFAULT_BORDER,
             int Step_trou = 0, int Window_Size=3);

// void rmedian( Ifloat & Imag1, Ifloat  & Imag2,
//             type_border Type = DEFAULT_BORDER, int Window_Size=3);

void im_reduce_size_2 (const Ifloat &Imag, Ifloat &Imag_Out);
void im_increase_size_2 (const Ifloat &Pict_in, Ifloat &Pict_out, 
             type_border Border=DEFAULT_BORDER);
void im_extract (const Ifloat &Imag, Ifloat &Imag_Out);
void im_extract (const Icomplex_f &Imag, Ifloat &Imag_Out, float Norm=1.);
void im_extend (const Ifloat &Imag, Ifloat &Imag_Out);
void im_extend (const Ifloat &Imag, Icomplex_f &Imag_Out);
void im_zero_padding (const Ifloat &Imag, Ifloat &Imag_Out);
void im_zero_padding (const Ifloat &Imag,  Icomplex_f &Imag_Out);

void morpho_erosion (Ifloat &Imag1, Ifloat& Imag2, int Window_Size = 3);
void morpho_dilation (Ifloat &Imag1, Ifloat &Imag2, int Window_Size = 3);
void morpho_cercle_erosion (Ifloat &Imag1, Ifloat& Imag2, int Window_Size=5);
void morpho_cercle_dilation (Ifloat &Imag1, Ifloat &Imag2, int Window_Size=5);

void im_detect(Ifloat &Imag_in, Ifloat &Imag_Detect, 
               float **x_ef, float **y_ef,  float **f_ef, int *nmax);
void im_segment (Ifloat &Imag, Ifloat &Segment, int &NLabel, 
                float Level=0.5, Bool CleanBord = False, int FirstLabel=1);
void im_thin(Ifloat &Imag,float Level);
void im_edge_kirsh(Ifloat &Ima, Iint &Angle, Ifloat &ImaEdge, int Step=0);
void im_threshold_kirsh_edge(Ifloat &Edge, Iint &Angle, 
                             float Noise, float Nsigma, int Step);

void im_histo (Ifloat &Pict, int **Tab_Histo, int &Nbr_Val);

void morpho4_dilation (Ifloat &Imag1, Ifloat &Imag2, int Step_trou=0);
void morpho4_erosion (Ifloat &Imag1, Ifloat& Imag2, int Step_trou=0);
void im_block_extend(const Ifloat &Imag, Ifloat &Imag_Out);
void im_bilinear_interp (Ifloat &INimage, Ifloat &OUTimage);
//void im_cf_interp (Ifloat &Image, int Nls, int Ncs);

void real (Ifloat& Real, const Icomplex_f &);
void imag (Ifloat& Imag, const Icomplex_f &);
void real (Ifloat& Real, const Icomplex_d &);
void imag (Ifloat& Imag, const Icomplex_d &);

const Ifloat real(const Icomplex_f &);
const Ifloat imag(const Icomplex_f &);

void fft2d (const Icomplex_f &, const Icomplex_f &, int Dir=1, bool normalize=false);
void fft2d (Icomplex_f &, int Dir=1, bool normalize=false);
	    
	    
void im_extract (const Icomplex_f &Imag, Icomplex_f &Imag_Out);
void im_extend (const Icomplex_f &Imag, Icomplex_f &Imag_Out);

const Icomplex_f invfft(const Icomplex_f &);
 
Icomplex_f icomplex_f (const Iint &);
Icomplex_f icomplex_f (const Ifloat &);

void fft1d (const Ifloat &Signal, const Icomplex_f & Signal_cf, int Dir=1);
void fft1d (const Icomplex_f &Signal, const Icomplex_f & Signal_cf, int Dir=1);	    


#endif

