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
**    Date:  96/05/02
**    
**    File:  MR1D_Sigma.cc
**    
**    Modification history:
**           25-OCT-1996 R Gastaud add mr1d_tab_noise
**
*******************************************************************************
**
**    DESCRIPTION 
**    -----------    
**
**    this module contains the routines for 1D filtering
**
**************************************************************************/ 

#include <string>  
#include <cmath>  
#include <iostream>
#include "Array.h"
#include "IM_Noise.h"
#include "MR1D_Obj.h"
#include "MR1D_Sigma.h"


float Tab1DSignifLevel[MAX_SCALE_1D];

#define NOISE_STUDY 0
// see also  NbrScale and Nx

/*******************************************************************/

void mr1d_noise_compute (MR_1D & MR_Data)
{
    int Nbr_Plan = MR_Data.nbr_scale();
    int MedianWindowSize = MR_Data.MedianWinSize;
    type_trans_1d Transform = MR_Data.Type_Transform;
    mr1d_noise_compute (Nbr_Plan, Transform, MR_Data.Border, MedianWindowSize);
}

/*******************************************************************/

void mr1d_noise_compute (int Nbr_Plan, type_trans_1d Transform, 
                         type_border Border, int WinSize, sb_type_norm Norm)
{
    int i;
    //static int Last_Nbr_Plan=0;
    //static type_trans_1d Last_Transform;
    //static int Pas = 0;
    double s2;
    switch(Transform)
    {
        case TO1_MALLAT:
	case WP1_MALLAT:
	case TU1_MALLAT:
             if (Norm == NORM_L1)
             {
                s2 = 1. / sqrt(2.);
 	        for (i = 0; i < MAX_SCALE_1D; i++) 
	            Tab1DSignifLevel[i] =  pow(s2, (double) (i+1));
             }
             else  
               for (i = 0; i < MAX_SCALE_1D; i++)  Tab1DSignifLevel[i] = 1.;
             break;
        case TO1_PAVE_HAAR:
	     s2 = 1. / sqrt(2.);
 	     for (i = 0; i < MAX_SCALE_1D; i++) 
	         Tab1DSignifLevel[i] =  pow(s2, (double) (i+1));
             break;
        case TO1_PAVE_LINEAR:
             Tab1DSignifLevel[0] = 0.605698;
             Tab1DSignifLevel[1] = 0.329922;
	     Tab1DSignifLevel[2] = 0.210199;
	     Tab1DSignifLevel[3] = 0.144264;
	     Tab1DSignifLevel[4] = 0.0974801;
	     Tab1DSignifLevel[5] = 0.0685790;
	     Tab1DSignifLevel[6] = 0.0534222;
	     Tab1DSignifLevel[7] = 0.0370909;
	     Tab1DSignifLevel[8] = 0.0287061;
             for (i = 9; i < MAX_SCALE_1D; i++) Tab1DSignifLevel[i] = 0.0177307;
             break;
 	case TO1_PAVE_B1SPLINE:
             Tab1DSignifLevel[0] = 0.403798;
             Tab1DSignifLevel[1] = 0.276628;
	     Tab1DSignifLevel[2] = 0.196208;
	     Tab1DSignifLevel[3] = 0.142940;
	     Tab1DSignifLevel[4] = 0.101959;
	     Tab1DSignifLevel[5] = 0.0746539;
	     Tab1DSignifLevel[6] = 0.0574192;
	     Tab1DSignifLevel[7] = 0.0428343;
	     Tab1DSignifLevel[8] = 0.0319696;
             for (i = 9; i < MAX_SCALE_1D; i++) Tab1DSignifLevel[i] = 0.0215028;
             break;
        case TO1_PAVE_B3SPLINE:
             Tab1DSignifLevel[0] = 0.721618;
             Tab1DSignifLevel[1] = 0.286209;
	     Tab1DSignifLevel[2] = 0.178281;
	     Tab1DSignifLevel[3] = 0.121010;
	     Tab1DSignifLevel[4] = 0.0831040;
	     Tab1DSignifLevel[5] = 0.0588543;
	     Tab1DSignifLevel[6] = 0.0446171;
	     Tab1DSignifLevel[7] = 0.0289142;
	     Tab1DSignifLevel[8] = 0.0177307;
             for (i = 9; i < MAX_SCALE_1D; i++) Tab1DSignifLevel[i] = 0.0177307;
             break;     
        case TO1_PAVE_B3SPLINE_GEN2:
         Tab1DSignifLevel[0] = 0.649506;
          Tab1DSignifLevel[1] = 0.119327;
	     Tab1DSignifLevel[2] = 0.0485716;
	     Tab1DSignifLevel[3] = 0.0231648;
	     Tab1DSignifLevel[4] = 0.0114491;
	     Tab1DSignifLevel[5] = 0.00570809;
	     Tab1DSignifLevel[6] = 0.00285199;
	     Tab1DSignifLevel[7] = 0.00142574;
	     Tab1DSignifLevel[8] = 0.000712838;
             for (i = 9; i < MAX_SCALE_1D; i++) Tab1DSignifLevel[i] =0.000712838;
             break;     
          case TU1_UNDECIMATED_NON_ORTHO:
             Tab1DSignifLevel[0] = 0.723490;
             Tab1DSignifLevel[1] = 0.285450;
	     Tab1DSignifLevel[2] = 0.177948 ;
	     Tab1DSignifLevel[3] = 0.122223;
	     Tab1DSignifLevel[4] = 0.0858113;
	     Tab1DSignifLevel[5] = 0.0605703;
	     Tab1DSignifLevel[6] = 0.0428107;
	     Tab1DSignifLevel[7] = 0.23;
	     Tab1DSignifLevel[8] = 0.017;
             for (i = 9; i < MAX_SCALE_1D; i++) Tab1DSignifLevel[i] = 0.017;
             break;
 	 case TO1_PAVE_B3_DERIV:
             Tab1DSignifLevel[0] = 0.699926;
             Tab1DSignifLevel[1] = 0.329224;
	     Tab1DSignifLevel[2] = 0.209089;
	     Tab1DSignifLevel[3] = 0.141389;
	     Tab1DSignifLevel[4] = 0.0978073;
	     Tab1DSignifLevel[5] = 0.0709058;
	     Tab1DSignifLevel[6] = 0.0542333;
	     Tab1DSignifLevel[7] = 0.0378617;
	     Tab1DSignifLevel[8] = 0.0278816;
             for (i = 9; i < MAX_SCALE_1D; i++) Tab1DSignifLevel[i] = 0.0216781;
             break;
        case TM1_PAVE_MEDIAN:
             if (WinSize == 5)
             {
                Tab1DSignifLevel[0] = 0.944437;
                Tab1DSignifLevel[1] = 0.334682;
                Tab1DSignifLevel[2] = 0.230710;
                Tab1DSignifLevel[3] = 0.195319;
                Tab1DSignifLevel[4] = 0.135050;
                Tab1DSignifLevel[5] = 0.100175;
                Tab1DSignifLevel[6] = 0.083291;
                Tab1DSignifLevel[7] = 0.095585;
                Tab1DSignifLevel[8] = 0.1;
             }
             else 
             {
                Tab1DSignifLevel[0] = 0.887662;
                Tab1DSignifLevel[1] = 0.387685;
                Tab1DSignifLevel[2] = 0.333953;
                Tab1DSignifLevel[3] = 0.263442;
                Tab1DSignifLevel[4] = 0.199089;
                Tab1DSignifLevel[5] = 0.131932;
                Tab1DSignifLevel[6] = 0.107235;
                Tab1DSignifLevel[7] = 0.140406;
                Tab1DSignifLevel[8] = 0.1;
             }
             for (i = 9; i < MAX_SCALE_1D; i++) Tab1DSignifLevel[i] = 0.1;
             break;
        case TO1_PYR_LINEAR:
             Tab1DSignifLevel[0] = 0.605698;
             Tab1DSignifLevel[1] = 0.326176;
	     Tab1DSignifLevel[2] = 0.206345;
	     Tab1DSignifLevel[3] = 0.143773;
	     Tab1DSignifLevel[4] = 0.0957876;
	     Tab1DSignifLevel[5] = 0.0696705;
	     Tab1DSignifLevel[6] = 0.0558049;
	     Tab1DSignifLevel[7] = 0.0354839;
	     Tab1DSignifLevel[8] = 0.0262708;
             for (i = 9; i < MAX_SCALE_1D; i++) Tab1DSignifLevel[i] = 0.0148184;
             break;
        case TO1_PYR_B3SPLINE:
             Tab1DSignifLevel[0] = 0.715785;
             Tab1DSignifLevel[1] = 0.282886;
	     Tab1DSignifLevel[2] = 0.176579;
	     Tab1DSignifLevel[3] = 0.118032;
	     Tab1DSignifLevel[4] = 0.0814612;
	     Tab1DSignifLevel[5] = 0.0593129;
	     Tab1DSignifLevel[6] = 0.0463765;
	     Tab1DSignifLevel[7] = 0.0317115;
	     Tab1DSignifLevel[8] = 0.0216181;
             for (i = 9; i < MAX_SCALE_1D; i++) Tab1DSignifLevel[i] = 0.0169869;
             break;

        case TM1_PYR_MEDIAN:
             if (WinSize == 5)
             {
                Tab1DSignifLevel[0] = 0.906316;
                Tab1DSignifLevel[1] = 0.351941;
                Tab1DSignifLevel[2] = 0.259543;
                Tab1DSignifLevel[3] = 0.191981;
                Tab1DSignifLevel[4] = 0.145102;
                Tab1DSignifLevel[5] = 0.101198;
                Tab1DSignifLevel[6] = 0.088637;
                Tab1DSignifLevel[7] = 0.110414;
                Tab1DSignifLevel[8] = 0.1;
             }
             else 
             {
                Tab1DSignifLevel[0] = 0.839884;
                Tab1DSignifLevel[1] = 0.463646;
                Tab1DSignifLevel[2] = 0.384971;
                Tab1DSignifLevel[3] = 0.288793;
                Tab1DSignifLevel[4] = 0.216289;
                Tab1DSignifLevel[5] = 0.159931;
                Tab1DSignifLevel[6] = 0.153861;
                Tab1DSignifLevel[7] = 0.169053;
                Tab1DSignifLevel[8] = 0.1;
             }
             for (i = 9; i < MAX_SCALE_1D; i++) Tab1DSignifLevel[i] = 0.1;
             break;
        default:
#if NOISE_STUDY
	     if ((Pas == 0) || (Nbr_Plan > Last_Nbr_Plan) || (Transform != Last_Transform))
             {
	        int NSimu = 10000;
	        MR_1D PyrNoise(NSimu, Transform,"Noise Estimation", Nbr_Plan);
	        fltarray NoiseImage (NSimu);
	        PyrNoise.MedianWinSize = WinSize;
	        PyrNoise.Border = Border;

	        Last_Nbr_Plan = Nbr_Plan;
	        Last_Transform = Transform;
	        Pas = 1;

	        randomn (NoiseImage, NSimu);
		
		
	        PyrNoise.transform (NoiseImage);
	        for (i = 0; i < Nbr_Plan; i++) {
//                 Tab1DSignifLevel[i] = (PyrNoise.scale(i)).sigma ();
//                 Tab1DSignifLevel[i] = sigma (PyrNoise.scale(i));
//                 Phil 08/09/99 Modif function scale
                   fltarray prov;
		   PyrNoise.scale(prov, i);
		   Tab1DSignifLevel[i] = prov.sigma ();
                }

	        for (i = Nbr_Plan; i < MAX_SCALE_1D; i++) Tab1DSignifLevel[i] = 0.;
	     }
#else
             cout << "Error: bad transform ... " << endl;
             exit(-1);   
#endif
             break;
     }
#if NOISE_STUDY 
fprintf (stderr, "\n\nNOISE ESTIMATION\n");
for (i = 0; i < 10; i++)
   fprintf (stderr, "Scale(%d) = %f\n", i+1, Tab1DSignifLevel[i]); 
#endif
}

/****************************************************************************/

float mr1d_tab_noise(int s)
{
    return Tab1DSignifLevel[s];
}
   
/****************************************************************************/

float mr1d_level_noise (MR_1D &MR_Data, int s, float N_Sigma, float N_Sigma0)
{
    float L;

    switch (MR_Data.Set_Transform)
    {
        case TRANS1_PAVE:
        case TRANS1_PYR:
           if (s == 0) L = N_Sigma0 * Tab1DSignifLevel[s];
           else L = N_Sigma * Tab1DSignifLevel[s];
           break;
        default: L = 0.;break;
    }
    return (L);
}

/************************************************************************/

float mr1d_detect_noise(fltarray &Sig, int NbrScale, float N_Sigma, 
                                       int Niter, int NiterClip)
{
    float Noise, OldEst;
    int Nx = Sig.axis(1);
    unsigned char *PixNoise = new unsigned char [Nx];
    float *TabLevel = new float [NbrScale];  
    float std, Val, Cvg;
    int s,i,nb,k=0;

    mr1d_noise_compute (NbrScale, TO1_PAVE_B3SPLINE);

    Noise = detect1d_noise_from_med(Sig); 
    fltarray Buff(Nx);
    filt1d_mediane  (Sig, Buff, Nx, 3);
    Buff = Sig - Buff;
//    Noise = sigma_clip(Buff, NiterClip) / 0.972463;
    Noise = Buff.sigma_clip (NiterClip) / 0.972463;
    
    // Kill too large value in the input signal for noise estimation
    for (i=0;i<Nx;i++) 
    {
       if (ABS(Buff(i)) > 10*Noise)  
       {
          Buff(i) = Sig(i) - Buff(i);  // takes the median value
          PixNoise[i] = 0;
       }
       else Buff(i) = Sig(i);
    }

    MR_1D MR_Data_Pave (Nx, TO1_PAVE_B3SPLINE, "MR_Deglitch", NbrScale);
    MR_Data_Pave.transform (Buff);
    do
    {
        OldEst = Noise;

        for (i=0;i<NbrScale-1;i++) 
           TabLevel[i] = Noise*mr1d_level_noise(MR_Data_Pave,i,N_Sigma,N_Sigma);
        for (i=0;i<Nx;i++) 
        {
           if (ABS( Buff(i)-Sig(i) ) < 5*Noise)
           {
              s = 0;
              while ((s < NbrScale-1) && (ABS(MR_Data_Pave(s,i)) < TabLevel[s])) s++;
              if (s == NbrScale-1) PixNoise[i] = 1;
              else PixNoise[i] = 0;
           }
        }

        /* noise standard deviation calculation with the mask*/
        std = 0.0;
        nb=0;
        for (i=0;i<Nx;i++)
        if (PixNoise[i] == 1)
        {
	   Val = Sig(i) - MR_Data_Pave(NbrScale-1,i);
	   std = std + Val*Val;
	   nb ++;
        }

        if (nb == 0)  std = 0.0;
        else std = sqrt(std/nb) / 0.973703; 
        Noise = std;
        k++;
        if (Noise > FLOAT_EPSILON) Cvg = ABS(OldEst-Noise)/Noise;
#if PRINT_DATA
      cout << "Iter " <<k <<" Old Noise = "<<OldEst<< "  cvg = " << Cvg << endl;
        cout << "       nb = " << nb << " New noise = " << Noise << endl;
#endif
        } while (k <= Niter);

 //      } while ( (k <= Niter) && (Cvg > 1.e-5));

    delete[] TabLevel ; 
    delete[] PixNoise ; 
    return (Noise);
}


/************************************************************************/

float detect1d_noise_from_med (fltarray &Sig, int NiterClip)
{
    float Noise;
    fltarray Buff(Sig.axis(1));

    filt1d_mediane  (Sig, Buff, Sig.axis(1), 3);
    Buff = Sig - Buff;
    Noise = Buff.sigma_clip (NiterClip);
//    Noise = sigma_clip (Buff, NiterClip);
    return (Noise/0.972463);
}
