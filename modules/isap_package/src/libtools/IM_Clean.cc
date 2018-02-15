/*******************************************************************************
**
**    UNIT
**
**    Version: 1.0
**
**    Author: Jean-Luc Starck
**
**    Date:  08/07/99 
**    
**    File:  IM_Clean.cc
**
*******************************************************************************
**
**    DESCRIPTION  Image Deconvolution: CLEAN METHOD
**    -----------  
**                 
**
*******************************************************************************
**
** void dec_clean (Ifloat & Imag, Ifloat & Psf, Ifloat & Imag_Out, 
**                 float Fwhm, float Noise, float Gamma, int Niter, 
**                 int K_Aver, Bool AddResi, Bool SearchPositiv, Bool UseEnergy,
** 		   Bool Verbose)
**
** Deconvolution by the CLEAN method
** Imag = in: dirty map
** Psf = in: dirty beam
** Imag_Out = out: clean map
** Fwhm = in: FWHM of the clean beam (if Fwhm > 0) the clean map is convolved
**            by the cleab beam
** Noise: in: detection level
** Gamma: in: loop gain
** Niter: in: maximum number of iterations
** K_Aver: in: end loop parameter: stop if Max < K_Aver*mean(residual map)
** AddResi: in: if true, the residual map is added to the solution
** SearchPositiv: in: if true, only positive peaks are searched
** UseEnergy: in: if true, maxima are search on the energy map
** Verbose: in: verbose mode
** 
******************************************************************************/
		
#include "IM_Obj.h"
#include "IM_Math.h"
#include "IM_Deconv.h"
#include "IM_Noise.h"
#include "IM_Sigma.h"
 
/***************************************************************************/

static float substract_imag (Ifloat & Dirty_Map, Ifloat & Dirty_Beam, 
                      int Max_Beam_i, int Max_Beam_j, 
                      int Max_i, int Max_j, float  Coef_Mul)
 
// Soustraction de la Dirty_Beam centree en Max_Beam_i, Max_Beam_j
//    a la dirty map. 
//    Dirty_Map[i,j] = Dirty_Map[i,j] - 
//      Coef_Mul*Dirty_Beam [ Max_Beam_i + i - Max_i, Max_Beam_j + j - Max_Beam_j ]
// 
//    La procedure renvoie la moyenne en valeur absolue de l'image .
// 
{
    register int i,j;
    float Aver_Abs = 0.;
    int Dep_Beam_i, Dep_Beam_j;
    int Nl = Dirty_Map.nl();
    int Nc = Dirty_Map.nc();
    
    for (i = 0; i < Nl; i++)
    {
        Dep_Beam_i = i - Max_i + Max_Beam_i;
        for (j = 0; j < Nc; j++)
	{
            Aver_Abs += ABS(Dirty_Map(i,j));  
            Dep_Beam_j = j  - Max_j + Max_Beam_j;
            if ((Dep_Beam_i >= 0) && (Dep_Beam_i < Nl) &&
                    (Dep_Beam_j >= 0) && (Dep_Beam_j < Nc))
	    {
 	        Dirty_Map(i,j) -= Dirty_Beam(Dep_Beam_i,Dep_Beam_j)*Coef_Mul;
	    }
	}
    }
    Aver_Abs /= (Nl*Nc);

    return (Aver_Abs);
}

/************************************************************************/

static void dec_hogbom (Ifloat & Dirty_Map, Ifloat & Dirty_Beam, 
                 float Gamma, float Noise, Ifloat & Clean_Map,
                 int Niter, int K_Aver, Bool SearchPositiv, Bool Verbose)
{
    int iteration=0;
    int Max_i, Max_j;
    float Old_Max, Max_Ima = 0., Coef_Mul, Gain;
    float Aver_Abs = 0.;
    int Max_Beam_i, Max_Beam_j;
    float Max_Beam;

    Clean_Map.init();

    dec_pos_max (Dirty_Beam, Max_Beam_i, Max_Beam_j , Max_Beam, SearchPositiv);

    /* begin of the loop */
    Gain = Gamma / Max_Beam;
    
    /* search the max in the dirty map  */

    dec_pos_max (Dirty_Map, Max_i, Max_j , Max_Ima, SearchPositiv);
    Old_Max = ABS(Max_Ima) + 1.;

   while ( (ABS(Max_Ima) >= K_Aver*Aver_Abs) 
             && (ABS(Max_Ima) >= Noise) 
             && (iteration < Niter)
             && ( (ABS(Old_Max) > ABS(Max_Ima)) || (Max_Ima*Old_Max > 0.))
          )
    {
        Coef_Mul = Max_Ima * Gain;
        Clean_Map(Max_i, Max_j) += Coef_Mul;
  
	/* Subtract the dirty beam to dirty map. the dirty beam is centered
           on Imax, and weighted by Gamma */
         Aver_Abs = substract_imag (Dirty_Map, Dirty_Beam,  Max_Beam_i,  
	                            Max_Beam_j, Max_i, Max_j, Coef_Mul);
		      
         iteration ++;
 
         if ((Verbose == True) && (iteration % 10 == 0))
            printf ("%d:  Max_Ima = (%d,%d): %f\n", iteration,Max_i, Max_j, Max_Ima);
 
         // search the max in the dirty map 
         Old_Max = Max_Ima;
         dec_pos_max (Dirty_Map, Max_i, Max_j , Max_Ima, SearchPositiv);
    }
    
    if (Verbose == True)
    {
        printf ("End loop: number of iterations = %d: Max_Ima = (%d,%d): %f\n", iteration,
                                                  Max_i, Max_j, Max_Ima);
        printf("end loop criterion is :\n");                                      
        if (ABS(Max_Ima) < K_Aver*Aver_Abs)
              cout << "  Max_Ima < " << K_Aver << " * MeanResidu" << endl;
        if (iteration >= Niter) cout << "Nbr_of_iter = MAX_ITER"<< endl;
        if (ABS(Max_Ima) < Noise ) cout << "Max_Ima < Noise Level"   << endl;
        if (ABS(Old_Max) < ABS(Max_Ima) )  cout << "Old_Max > New_Max" << endl;
    }
}

/***************************************************************************/

static void map_energy (Ifloat & Imag, Ifloat &Imag_x, Ifloat &Imag_y, 
                        Ifloat &Plan_Energ)
// energy of an  image 
//  E = som_i( I_i^2 + Der_x(I)_i^2 +  Der_y(I)_i^2)
{
    int i,j;
    int Nl = Imag.nl();
    int Nc = Imag.nc();
    
    for (i = 0; i < Nl; i++)
    for (j = 0; j < Nc; j++)
     Plan_Energ(i,j)= Imag_x(i,j)*Imag_x(i,j) 
                         + Imag_y(i,j)*Imag_y(i,j) + Imag(i,j)*Imag(i,j);
}


/***************************************************************************/

static void fft_derive_image (Ifloat & Imag, Ifloat &Imag_x, Ifloat & Imag_y)
//  Calcul de la derive d'une image dans 
//    les directions x et y par utilisation de
//    la transormee de Fourier:
//         Der_x (I) = TF-1 [ TF(I)*2*i*PI*u ]
//         Der_y (I) = TF-1 [ TF(I)*2*i*PI*v ]

{
    int i, j;
    float u,v;
    int Nl = Imag.nl();
    int Nc = Imag.nc();
    Icomplex_f Pict (Nl, Nc, "Buffer");
    Icomplex_f Pict_x (Nl, Nc, "Buffer");
    Icomplex_f Pict_y (Nl, Nc, "Buffer");
 
    fft2d(Imag,Pict);
    // derivative calculation
    for (i = 0; i < Nl; i++) 
    for (j = 0; j < Nc; j++)
    {
        u = (float)(i - Nl/2) / (float) Nl;
        v = (float)(j - Nl/2) / (float) Nc;
	Pict_x(i,j) = Pict(i,j)*complex_f(0,2.*PI*u);
	Pict_y(i,j) = Pict(i,j)*complex_f(0,2.*PI*v);
    }

    fft2d(Pict_x, -1);
    fft2d(Pict_y, -1);
 
    for (i = 0; i < Nl; i++)
    for (j = 0; j < Nc; j++)
    {
        Imag_x(i,j) = Pict_x(i,j).real();
        Imag_y(i,j) = Pict_y(i,j).real();
    }
}

/***************************************************************************/

void clean_hogbom_energy (Ifloat & Dirty_Map, Ifloat & Dirty_Beam, 
                 float Gamma, float Noise, Ifloat & Clean_Map,
		 Ifloat & Tab_Dirty_X, Ifloat & Tab_Dirty_Y,
 		 Ifloat &Dirty_Beam_X, Ifloat &Dirty_Beam_Y,
                 int Niter, int K_Aver, Bool SearchPositiv, Bool Verbose)
{
    int Nl = Dirty_Map.nl();
    int Nc = Dirty_Map.nc();
    int iteration=0;
    int Max_i, Max_j;
    float Old_Max, Max_Ima = 0.;
    float Gain,Gain_X,Gain_Y;
    float Coef_Mul;
    float Aver_Abs = 0., Aver_Abs_X, Aver_Abs_Y;
    int Max_Beam_i, Max_Beam_j;
    float Max_Beam, Max_Dirty_Map;
    Ifloat Map_Energy(Nl,Nc,"map energy");

    Clean_Map.init();
  
    dec_pos_max (Dirty_Beam, Max_Beam_i, Max_Beam_j , Max_Beam, SearchPositiv);
    Gain = Gamma / Max_Beam;

    // energy map creation
    map_energy (Dirty_Map, Tab_Dirty_X, Tab_Dirty_Y, Map_Energy);

    /* Detection du maximum de la map d'energie */
    dec_pos_max (Map_Energy, Max_i, Max_j, Max_Ima, True);
    Max_Dirty_Map = Dirty_Map(Max_i,Max_j);

    /* Debut de la boucle */
    Gain = Gamma / Max_Beam;
    Gain_X = Gamma / Max_Beam;
    Gain_Y = Gamma / Max_Beam;
    Old_Max = ABS(Max_Ima) + 1;

    while (  (ABS(Max_Dirty_Map) >= Noise) 
             && (iteration < Niter)
             && (Old_Max - Max_Ima > FLOAT_EPSILON)
           )
    {
        Coef_Mul = Dirty_Map(Max_i,Max_j) * Gain;
        Clean_Map(Max_i,Max_j)  += Coef_Mul;

 	// Subtract the dirty beam to dirty map. the dirty beam is centered
        //   on Imax, and weighted by Gamma  
        Aver_Abs_X = substract_imag (Tab_Dirty_X, Dirty_Beam_X, Max_Beam_i, 
                                  Max_Beam_j, Max_i, Max_j, Coef_Mul);
        Aver_Abs_Y = substract_imag (Tab_Dirty_Y, Dirty_Beam_Y, Max_Beam_i, 
                                  Max_Beam_j, Max_i, Max_j, Coef_Mul);
        Aver_Abs = substract_imag (Dirty_Map, Dirty_Beam, Max_Beam_i, 
                                  Max_Beam_j, Max_i, Max_j, Coef_Mul);

        iteration ++;
        Old_Max = Max_Ima;
  	
        map_energy (Dirty_Map, Tab_Dirty_X, Tab_Dirty_Y, Map_Energy);
  	dec_pos_max (Map_Energy, Max_i, Max_j, Max_Ima, True);
	Max_Dirty_Map = Dirty_Map(Max_i,Max_j);
    }

    if (Verbose == True)
    {
        printf ("End loop: number of iterations = %d: Max_Ima = (%d,%d): %f\n", iteration,
                                                  Max_i, Max_j, Max_Dirty_Map);
        printf("end loop criterion is :\n");                                      
        if (iteration >= Niter) cout << "Nbr_of_iter = MAX_ITER"<< endl;
        if (ABS(Max_Dirty_Map) < Noise ) cout << "Max_Ima < Noise Level"   << endl;
        if (ABS(Old_Max) < ABS(Max_Ima) )  cout << "Old_Max > New_Max" << endl;
    }
}


/***************************************************************************/

void dec_clean (Ifloat & Imag, Ifloat & Psf, Ifloat & Imag_Out, 
                float Fwhm, float Noise, float Gamma, int Niter, 
                int K_Aver, Bool AddResi, Bool SearchPositiv, Bool UseEnergy,
		Bool Verbose)
{
    int Nl = Imag.nl();
    int Nc = Imag.nc();
    Ifloat Clean_Beam(Nl, Nc, "Clean_Beam");

    // find the dirac list position
    if (UseEnergy == False)
     dec_hogbom (Imag, Psf, Gamma, Noise, Imag_Out, Niter, K_Aver, SearchPositiv, Verbose);
    else
    {
      Ifloat Imag_x(Nl,Nc,"Imag_x");
      Ifloat Imag_y(Nl,Nc,"Imag_x");
      Ifloat Dirty_Beam_X(Nl,Nc,"Imag_x");
      Ifloat Dirty_Beam_Y(Nl,Nc,"Imag_x");
      fft_derive_image (Imag,  Imag_x,  Imag_y);
      fft_derive_image (Psf,  Dirty_Beam_X, Dirty_Beam_Y);
      clean_hogbom_energy (Imag, Psf, Gamma, Noise,Imag_Out, Imag_x, Imag_y,
		      Dirty_Beam_X,  Dirty_Beam_Y, Niter, K_Aver, SearchPositiv,
		      Verbose);
    }
    
    // create the clean beam 
    if (Fwhm > FLOAT_EPSILON)
    {
        Clean_Beam = im_gaussian (Nl, Nc, Fwhm);

       // Convolution of the result by the clean beam 
       fft2d_conv(Imag_Out, Clean_Beam, Imag_Out);
    }

    if (AddResi == True) Imag_Out += Imag;
}

/***************************************************************************/

