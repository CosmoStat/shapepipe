/*******************************************************************************
**
**    UNIT
**
**    Version: 3.3
**
**    Author: Jean-Luc Starck
**
**    Date:  96/07/02 
**    
**    File:  IM_Deconv.cc
**
*******************************************************************************
**
**    DESCRIPTION  Image Deconvolution routines
**    -----------  
**                 
**
*******************************************************************************
**
** void dec_pos_max (Ifloat &Tab, int &Ind_i, int &ind_j , float &Val_Max)
**
** seek the maximum in an image
** Tab = image
** Ind_i, ind_j = position of the max
** Val_Max = max values
**
*****************************************************************************
**
** void dec_line_column (int N, int &N_Out)
**
** N_Out =  the power of 2 superior or equal at N
**
*****************************************************************************
**
** void dec_center_psf (Ifloat &D_Beam, Ifloat &Psf)
**
** move the psf to the center of a bigger image
** 
*****************************************************************************
**
** void dec_convol_conj (Ifloat &Resi, Icomplex_f &Psf_cf)
**
** convolve the image Resi with the conjugate of the Fourier
** transform of the PSF
**
*****************************************************************************
** 
** void dec_convol (Ifloat &Resi, Icomplex_f &Psf_cf, Ifloat &Result)
**
** convolve the image Resi with the PSF
**
*******************************************************************************
**
** void dec_inverse (Ifloat &Resi, Icomplex_f &Psf_cf, Ifloat &Result)
**
** divides the fourier transform of Resi by Psf_cf
** the inverse fourier transform is put in Result
**
*******************************************************************************
**
** void dec_im_standart (Ifloat &Imag, Ifloat &Obj, Icomplex_f &Psf_cf,
**                     float Eps_cv DEFAULT_EPSILON_DECONV,
**                     int Nbr_Iter = DEFAULT_MAX_ITER_DECONV, 
**                     type_deconv Deconv = DEFAULT_DECONV, float TotalFlux=0.)
**
** deconvolution by standart algorithm
** Deconv = DEC_CITTERT, DEC_GRADIENT, or DEC_LUCY
** Obj contains the result
** a positivity constraint is applied
** if TotalFlux is set, then the obj is renormalized at each iteration
**
*******************************************************************************
**
** void dec_im_support (Ifloat &Imag, Ifloat &Obj, Icomplex_f &Psf_cf,
**             Ifloat &Support, float Eps_cv, int Nbr_Iter, type_deconv Deconv)
**
** Deconvolution similar to dec_im_standart by a support constraint is applied
** 
******************************************************************************/

// static char sccsid[] = "@(#)IM_Deconv.cc 3.3 96/07/02CEA 1994 @(#)";

#include "IM_Obj.h"
#include "IM_Math.h"
#include "IM_Deconv.h"
#include "IM_Noise.h"
#include "IM_Sigma.h"
#include "FFTN_2D.h"

#define INT_POW(x,y,z) { int l,xx,yy; xx = (x) ; yy = (y);  for (l=0,(z)=1;l<yy;++ l,z *= xx); }

#define DEBUG 0
#define VISU_DATA 0

Bool DecVerbose=False;

/****************************************************************************/
    
void dec_pos_max (Ifloat &Tab, int &Ind_i, int &ind_j , float &Val_Max, 
                  Bool SearchPositiv)
{
   int i,j;
   float Val_Abs_Max;
   int Bande = 1;
   int Nl = Tab.nl();
   int Nc = Tab.nc();

   Val_Max = 0.;
   Val_Abs_Max = 0.;

   for (i = Bande; i < Nl-Bande; i++)
   for (j = Bande; j < Nc-Bande; j++)
   {
      float Val_Abs = (SearchPositiv == False) ? ABS(Tab(i,j)) : Tab(i,j);
      if (Val_Abs > Val_Abs_Max)
      {
          Val_Abs_Max = Val_Abs;
          Val_Max = Tab(i,j);
          Ind_i = i;
          ind_j = j;
      }
   }
}

/***************************************************************************/

void dec_line_column (int N, int &N_Out)
/* We must check Nl to make sure it is a power of 2 */
{
    int len_exp,temp;

    len_exp = (int)(0.3+log((double)(N))/(log(2.0)));
    INT_POW(2,len_exp,temp);
    if (temp < N) temp *= 2; 
    N_Out = temp;
}

/**************************************************************************/

void dec_center_psf (Ifloat &D_Beam, Ifloat &Psf)
{
   int Ind_i, Ind_j,i,j,Dep_i, Dep_j;
   float Val_Max;
   int Nl_Beam = D_Beam.nl();
   int Nc_Beam = D_Beam.nc();
   int Nl = Psf.nl();
   int Nc = Psf.nc();
   double Flux = 0.;

   dec_pos_max (D_Beam, Ind_i, Ind_j , Val_Max, True);

   Psf.init();
   for (i = 0; i < Nl_Beam; i++)
   for (j = 0; j < Nc_Beam; j++)
   {
      Dep_i = i - Ind_i  + Nl / 2;
      Dep_j = j - Ind_j  + Nc / 2;
      if ((Dep_i >= 0) && (Dep_i < Nl) && (Dep_j >= 0) && (Dep_j < Nc))
      {
         Psf(Dep_i, Dep_j) = D_Beam(i,j);
	 Flux += Psf(Dep_i, Dep_j);   
      }
    }
    for (i = 0; i < Nl; i++)
    for (j = 0; j < Nc; j++) Psf (i,j) = (float) (Psf(i,j) / Flux);     
}

/**************************************************************************/

static void psf_extend (const Ifloat &Imag, Ifloat &Imag_Out)
{
    int Nl0 = Imag.nl();
    int Nc0 = Imag.nc();
    int Nl1 = Imag_Out.nl();
    int Nc1 = Imag_Out.nc();
    int i1,j1,i0,j0, Depi, Depj;
    double Flux = 0.;
    
    Depi = ceil((double)(Nl1 - Nl0) / 2);
    Depj = ceil((double)(Nc1 - Nc0) / 2);
    
    for (i1 = 0; i1 < Nl1; i1++)
    for (j1 = 0; j1 < Nc1; j1++)
    {
       i0 = i1 - Depi;
       j0 = j1 - Depj;
       if ((i0 < 0) || (j0 < 0) || (i0 >= Nl0) || (j0 >= Nc0)) 
            Imag_Out (i1,j1) = 0.;
       else 
       {
          Imag_Out (i1,j1) = Imag (i0,j0);
	  Flux += Imag_Out (i1,j1);
       }
    }
    for (i1 = 0; i1 < Nl1; i1++)
    for (j1 = 0; j1 < Nc1; j1++) 
                Imag_Out (i1,j1) = (float) (Imag_Out (i1,j1) / Flux);  
}

/**************************************************************************/

void psf_get(Ifloat &InPsf, Icomplex_f &Psf_cf, int ImaNl, int ImaNc, Bool PsfMaxShift)
{
   FFTN_2D FFT2D;
   int Nl = InPsf.nl();
   int Nc = InPsf.nc();
   if (Nl < ImaNl) Nl = ImaNl;
   if (Nc < ImaNc) Nc = ImaNc;
   
   int Nl1=Nl, Nc1=Nc;
   // Bool ModifSize=False;
   // dec_line_column (Nl, Nl1);
   // dec_line_column (Nc, Nc1);
   // if (Nl1 < Nc1) Nl1 = Nc1;
   // else Nc1 = Nl1;
   // if ((Nl != Nl1) || (Nc != Nc1)) ModifSize = True;
    
   Ifloat OutPsf(Nl1, Nc1, "Psf");
   if (PsfMaxShift == True) dec_center_psf (InPsf, OutPsf);
   else psf_extend (InPsf, OutPsf);
   // INFO(InPsf,"IN PSF");
   // INFO(OutPsf,"OUT PSF");
   Psf_cf.resize(Nl1, Nc1);
   // fft2d(OutPsf, Psf_cf);
   // io_write_ima_float("xx_psf.fits", OutPsf);
   FFT2D.fftn2d(OutPsf, Psf_cf);
}   

/****************************************************************************/

void psf_convol (Ifloat &Imag, Icomplex_f &Psf_cf, Ifloat &Imag_out, Bool FluxNorm)
{
   FFTN_2D FFT2D;
   int i,j,Nl = Imag.nl();
   int Nc = Imag.nc();
   int Nl1 = Psf_cf.nl();
   int Nc1 = Psf_cf.nc();
   Icomplex_f O_cf (Nl1, Nc1, "Buffer conv");
   float FluxPsf = (FluxNorm == False) ? 1.: Psf_cf(Nl1/2,Nc1/2).real();
   double Flux;
   
   if ((Nl != Nl1) || (Nc != Nc1))
   {
       im_extend (Imag, O_cf);
   }
   else 
   {
       for (i = 0; i < Nl; i++)
       for (j = 0; j < Nc; j++) O_cf(i,j) = complex_f(Imag(i,j),0);
   }
   // fft2d (O_cf);      
   FFT2D.fftn2d(O_cf);
   O_cf *= Psf_cf;
   // fft2d (O_cf, -1);
   FFT2D.ifftn2d(O_cf);
   if ((Nl != Nl1) || (Nc != Nc1)) 
   {
      im_extract (O_cf, Imag_out, FluxPsf);
      Flux = total(Imag) / total(Imag_out);
      for (i = 0; i < Nl; i++)
      for (j = 0; j < Nc; j++) Imag_out(i,j) = (float) (Imag_out(i,j) * Flux);
   }
   else for (i = 0; i < Nl; i++) 
        for (j = 0; j < Nc; j++)  Imag_out(i,j) = O_cf(i,j).real() / FluxPsf; 
   
   // INFO(Imag,"IN ima");
   // INFO(Imag_out,"OUT ima");
}

/***************************************************************************/

void psf_convol (Ifloat &Imag, Icomplex_f &Psf_cf, Bool FluxNorm)
{
   psf_convol (Imag, Psf_cf, Imag, FluxNorm);
}

/****************************************************************************/

void psf_convol (Ifloat &Imag, Ifloat &Psf, Bool FluxNorm)
{
   Icomplex_f Psf_cf;
   psf_get(Psf, Psf_cf, Imag.nl(), Imag.nc(), True);
   psf_convol (Imag, Psf_cf, Imag, FluxNorm);
}

/****************************************************************************/

void psf_convol (Ifloat &Imag, Ifloat &Psf, Ifloat &ImagOut, Bool FluxNorm, Bool PsfMaxShift)
{
   Icomplex_f Psf_cf;
   psf_get(Psf, Psf_cf, Imag.nl(), Imag.nc(), PsfMaxShift);
   psf_convol (Imag, Psf_cf, ImagOut, FluxNorm);
}

/***************************************************************************/

void psf_convol_conj (Ifloat &Imag, Icomplex_f &Psf_cf, Ifloat &Imag_out)
{
   FFTN_2D FFT2D;
   int i,j;
   int Nl = Imag.nl();
   int Nc = Imag.nc();
   int Nl1 = Psf_cf.nl();
   int Nc1 = Psf_cf.nc();
   Icomplex_f O_cf (Nl1, Nc1, "Buffer conv");
   float FluxPsf = Psf_cf(Nl1/2,Nc1/2).real();

   if ((Nl != Nl1) || (Nc != Nc1)) im_extend (Imag, O_cf);
   else 
   {
       for (i = 0; i < Nl; i++)
       for (j = 0; j < Nc; j++) O_cf(i,j) =  complex_f(Imag(i,j),0);
   }
   // fft2d (O_cf);
   FFT2D.fftn2d(O_cf);
   for (i = 0; i < Nl1; i++)
   for (j = 0; j < Nc1; j++) O_cf(i,j) *= conj(Psf_cf(i,j));
   // fft2d (O_cf, -1);
   FFT2D.ifftn2d(O_cf);
   if ((Nl != Nl1) || (Nc != Nc1)) im_extract (O_cf, Imag_out, FluxPsf);
   else for (i = 0; i < Nl*Nc; i++)  Imag_out(i) = O_cf(i).real() / FluxPsf;
 }

/****************************************************************************/

void psf_convol_conj (Ifloat &Imag, Icomplex_f &Psf_cf)
{
   psf_convol_conj (Imag, Psf_cf, Imag);
}

/****************************************************************************/

void dec_convol_conj (Ifloat &Resi, Icomplex_f &Psf_cf)
{
   FFTN_2D FFT2D;
   int i,j;
   int Nl = Resi.nl();
   int Nc = Resi.nc();
   Icomplex_f O_cf (Nl, Nc, "Buffer conv");
   float FluxPsf = Psf_cf(Nl/2,Nc/2).real();

   // fft2d (Resi, O_cf);
   FFT2D.fftn2d(Resi, O_cf);
   for (i = 0; i < Nl; i++)
   for (j = 0; j < Nc; j++) O_cf(i,j) *= conj(Psf_cf(i,j));
   
   // fft2d (O_cf, -1);
   FFT2D.ifftn2d(O_cf);
   for (i = 0; i < Nl*Nc; i++)  Resi(i) = O_cf(i).real() / FluxPsf; 
}

/***************************************************************************/

void dec_convol (Ifloat &Resi, Icomplex_f &Psf_cf, Ifloat &Result)
{
   FFTN_2D FFT2D;
   int i,Nl = Resi.nl();
   int Nc = Resi.nc();
   Icomplex_f O_cf (Nl, Nc, "Buffer conv");
   float FluxPsf = Psf_cf(Nl/2,Nc/2).real();
   
   // fft2d (Resi, O_cf);
   FFT2D.fftn2d(Resi, O_cf);
   O_cf *= Psf_cf;
   FFT2D.ifftn2d(O_cf);
   // fft2d (O_cf, -1);
   for (i = 0; i < Nl*Nc; i++)   Result(i) = O_cf(i).real() / FluxPsf; 
}

/***************************************************************************/

/***************************************************************************/

void dec_inverse (Ifloat &Imag, Icomplex_f &Psf_cf, Ifloat &Result, float Eps)
{
   FFTN_2D FFT2D;
   int i,j;
   int Nl = Imag.nl();
   int Nc = Imag.nc();
   int Nl1 = Psf_cf.nl();
   int Nc1 = Psf_cf.nc();
   Icomplex_f O_cf (Nl1, Nc1, "Buffer conv");
   // float FluxPsf = (FluxNorm == False) ? 1.: Psf_cf(Nl1/2,Nc1/2).real();
   float FluxPsf = 1. / Psf_cf(Nl1/2,Nc1/2).real();
   
   if ((Nl != Nl1) || (Nc != Nc1)) im_extend (Imag, O_cf);
   else 
   {
       for (i = 0; i < Nl; i++)
       for (j = 0; j < Nc; j++) O_cf(i,j) = complex_f(Imag(i,j),0);
   }
   FFT2D.fftn2d(O_cf);
   // fft2d (O_cf);
   float Den;

   for (i = 0; i < Nl1; i++) 
   for (j = 0; j < Nc1; j++) 
   {
       Den = norm (Psf_cf(i,j));
       O_cf(i,j) *= conj(Psf_cf(i,j));
       if (Den > Eps) O_cf(i,j) /= Den;
       else O_cf(i,j) = complex_f(0.,0.);
   }
   FFT2D.ifftn2d(O_cf);
   // fft2d (O_cf, -1);
   if ((Nl != Nl1) || (Nc != Nc1)) im_extract (O_cf, Result, FluxPsf);
   else for (i = 0; i < Nl; i++)
        for (j = 0; j < Nc; j++)  Result(i,j) = O_cf(i,j).real() / FluxPsf;
}

/************************************************************************/

void dec_im_standart (Ifloat &Imag, Ifloat &Obj, Icomplex_f &Psf_cf,
                     float Eps_cv, int Nbr_Iter, type_deconv Deconv,
                     float TotalFlux, float Noise_Ima)
{
    FFTN_2D FFT2D;
    int Nl = Imag.nl();
    int Nc = Imag.nc();
    Ifloat Imag_n(Nl, Nc, "Imag_n");
    Ifloat Resi(Nl, Nc, "Resi");
    int i,j,Iter=0;
    Bool Stop = False;
    float Xi, Sigma2, OldXi=1e10, Sigma;

    /* Noise estimation in the Data */
    if (Noise_Ima < FLOAT_EPSILON) 
    {
        Noise_Ima = detect_noise_from_med (Imag);
    } 

    /* Initialization */
    Sigma2 = Noise_Ima*Noise_Ima;

    /* Residual estimation */ 
    // dec_convol (Obj, Psf_cf, Imag_n);
    FFT2D.convolve(Obj, Psf_cf, Imag_n);
    for (i = 0; i < Nl; i++) 
    for (j = 0; j < Nc; j++) Resi(i,j) = Imag(i,j) - Imag_n(i,j);

    do
    {

       switch (Deconv)
       {
           case DEC_LUCY:
           case DEC_MR_LUCY:
                /* Calculate the multiplicate term in Lucy's iteration */
                for (i = 0; i < Nl; i++) 
                for (j = 0; j < Nc; j++) 
                {
                    if (ABS(Imag_n(i,j)) > FLOAT_EPSILON)
                                Resi(i,j) = Imag(i,j) / Imag_n(i,j);
                    else Resi(i,j) = 1.;
                }
                // dec_convol_conj (Resi, Psf_cf);
		FFT2D.convolve_conj(Resi, Psf_cf); 
                Obj *= Resi;
                break;
	   case DEC_GRADIENT:
           case DEC_MR_GRADIENT:
                /* Calculate the multiplicate term in Gradient's iteration */
                FFT2D.convolve_conj(Resi, Psf_cf);
		// dec_convol_conj (Resi, Psf_cf); 
           case DEC_CITTERT:
           case DEC_MR_CITTERT:
                /* calculate the next estimation of the object with positivity 
                   constraint 
                */
                for (i = 0; i < Nl; i++) 
                for (j = 0; j < Nc; j++) 
                {
                    if (Resi(i,j) >= - Obj(i,j)) Obj(i,j) += Resi(i,j);
                    else Obj(i,j) = 0.;
                }
                break;
           default:
                cerr << "mr_deconv: Not implemented in this procedure ... ";
                cerr << endl;
                exit (0);
                break;
       }

       if (TotalFlux > FLOAT_EPSILON) norm_flux (Obj, TotalFlux);

       /* Residual estimation */ 
       // dec_convol (Obj, Psf_cf, Imag_n);
       FFT2D.convolve(Obj, Psf_cf, Imag_n);
       Xi = 0.;
       for (i = 0; i < Nl; i++) 
       for (j = 0; j < Nc; j++)
       {
           Resi(i,j) = Imag(i,j) - Imag_n(i,j);
           Xi += Resi(i,j)*Resi(i,j);
       }
       Sigma = Xi / (Nl*Nc);
       Xi = Sigma / Sigma2;

#if VISU_DATA
       if (Iter % 1 == 0)
       {
          printf("%d: Sigma = %f, Xi2 = %f\n", Iter, sqrt(Sigma), Xi);
       }
#endif
       Iter ++;
       if (Iter > Nbr_Iter) Stop = True;
       else if (Eps_cv > FLOAT_EPSILON)
       {
          if (OldXi < Xi) Stop = True;
          else if ((OldXi - Xi) < Eps_cv) Stop = True;
       }
       OldXi = Xi;
    } while (Stop == False);
    if (DecVerbose == True)
            printf ("Iteration Number: %d, Xi2 = %f \n", Iter, Xi);
}

/************************************************************************/

void dec_im_support (Ifloat &Imag, Ifloat &Obj, Icomplex_f &Psf_cf,
               Ifloat &Support, float Eps_cv, int Nbr_Iter, type_deconv Deconv)
{
    FFTN_2D FFT2D;
    int Nl = Imag.nl();
    int Nc = Imag.nc();
    Ifloat Imag_n(Nl, Nc, "Imag_n");
    Ifloat Resi(Nl, Nc, "Imag_n");
    float Sigma, Old_Sigma, Delta;
    int i,j,Iter=0;
 
    float SigmaImag = sigma (Imag);

    /* Initialization */
    Delta = Sigma = 1.e20;
    do
    {
       Old_Sigma = Sigma; 

       /* Residual estimation */ 
       // dec_convol (Obj, Psf_cf, Imag_n);
       FFT2D.convolve(Obj, Psf_cf, Imag_n);
       for (i = 0; i < Nl; i++) 
       for (j = 0; j < Nc; j++)
          if (Support(i,j) > FLOAT_EPSILON)  
                   Resi(i,j) = Support(i,j) * (Imag(i,j) - Imag_n(i,j));
          else Resi(i,j) = 0.;

       /* Standard deviation of the residual */
       Sigma = sigma(Resi);
 
       switch (Deconv)
       {
           case DEC_LUCY:
           case DEC_MR_LUCY:
                /* Calculate the multiplicate term in Lucy's iteration */
                for (i = 0; i < Nl; i++) 
                for (j = 0; j < Nc; j++) 
                {
                    if ( (Support(i,j) > FLOAT_EPSILON) &&
                         (ABS(Imag_n(i,j)) > FLOAT_EPSILON)) 
                                           Resi(i,j) /= Imag_n(i,j);
                    else Resi(i,j) = 0.;
                }
                // dec_convol_conj (Resi, Psf_cf); 
                FFT2D.convolve_conj(Resi, Psf_cf);
               Resi *= Obj;
                break;
           case DEC_GRADIENT:
           case DEC_MR_GRADIENT:
                /* Calculate the multiplicate term in Gradient's iteration */
                // dec_convol_conj (Resi, Psf_cf); 
		FFT2D.convolve_conj(Resi, Psf_cf);
                break;
	   case DEC_CITTERT:
           case DEC_MR_CITTERT:
                break;
           default:
                cerr << "mr_deconv: Not implemented in this procedure ... ";
                cerr << endl;
                exit (0);
                break;
       }
       /* calculate the next estimation of the object with positivity 
          constraint 
       */
       for (i = 0; i < Nl; i++) 
       for (j = 0; j < Nc; j++) 
       {
           if (Resi(i,j) >= - Obj(i,j)) Obj(i,j) += Resi(i,j);
           else Obj(i,j) = 0.;
       }
       Delta = (Old_Sigma - Sigma) / SigmaImag;

       if (DecVerbose == True)
       {
          if ((Iter > 0) && (Iter % 1 == 0))
          {
          printf("%d: Sigma(residual) = %f, Delta = %f\n", Iter, Sigma, Delta );
          }
       }
       Iter ++;
    } while ((Iter <= Nbr_Iter) && (Delta > Eps_cv));
}

/***************************************************************************/
