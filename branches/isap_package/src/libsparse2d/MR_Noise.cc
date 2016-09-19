/*******************************************************************************
**
**    UNIT
**
**    Version: 3.2
**
**    Author: Jean-Luc Starck
**
**    Date:  96/06/13 
**    
**    File:  MR_Noise.cc
**
*******************************************************************************
**
**    DESCRIPTION routines for noise treatment 
**    -----------  
**                 
**
*******************************************************************************
**
**  void noise_compute (int Nbr_Plan, type_transform Transform, int Nl, int Nc)
** 
** fill the array  TabSignificantLevel with  the good values for the transform
** Transform
**           
*******************************************************************************
**
** float mr_noise_estimation (MultiResol &MR_Data)
**
** return the noise estimation computed from the first scale of the
** transformation MR_Data
**
*******************************************************************************
**
** void mr_threshold (MultiResol &MR_Data, float &Noise_Ima, float N_Sigma)
**
** Threshold a multiresolution transform:
** if Noise_Ima = 0, compute the noise
** the Threshold level is at N_Sigma*Noise
**
*******************************************************************************
**
** void noise_support_threshold (MultiResol &MR_Data, float Noise_Ima, 
**                            float N_Sigma, Bool SetSupport)
**
** Threshold a multiresolution transform:
** if Noise_Ima = 0, compute the noise
** the Threshold level is at N_Sigma*Noise
** if SetSupport is True, then if MR_Data(s,i,j) > DataSupport (s,i,j) = 1
**
*******************************************************************************
**
** Bool mallat_significant (MultiResol &MR_Data, int s, int i, int j,
**                         float Level, details Det)
**
** return True if the wavelet coefficient is significant
**
*******************************************************************************
**
** float mr_level_noise (MultiResol &MR_Data, int s, float N_Sigma)
**
** return the significant level of the multiresolution MR_Data at the 
** scale s, for an assumtion of gaussian noise with a sigma equal to 1.
** the level is multiplied by N_Sigma 
**  
******************************************************************************/

// static char sccsid[] = "@(#)MR_Noise.cc 3.2 96/06/13 CEA 1994 @(#)";


#include "IM_Obj.h"
#include "MR_Obj.h"
#include "IM_IO.h"
#include "IM_Math.h"
#include "IM_Noise.h"
#include "MR_Noise.h"
#include "MR_Sigma.h"

float TabSignificantLevel[MAX_BAND];

extern MultiResol DataSupport;
 
#define PRINT_DATA 0

/****************************************************************************/
 
void noise_compute (MultiResol &MR_Data)
{
    int Nbr_Plan = MR_Data.nbr_scale();
    int Nl = MR_Data.size_ima_nl();
    int Nc = MR_Data.size_ima_nc();
    int MedianWindowSize = MR_Data.MedianWindowSize;
    type_transform Transform = MR_Data.Type_Transform;
    noise_compute (Nbr_Plan, Transform, Nl, Nc, MedianWindowSize);
}

/****************************************************************************/

void noise_compute (int Nbr_Plan, type_transform Transform,
                    int Nl, int Nc, int MedianWindowSize)
{
    int i;
    for (i = 0; i < MAX_SCALE; i++)  
          TabSignificantLevel[i] = mr_scale_norm(i, Transform, MedianWindowSize); 
    
    if (0)
    {
        MultiResol PyrNoise(Nl,Nc,Nbr_Plan,Transform,"Noise Estimation");
        Ifloat NoiseImage (Nl, Nc, "IMAGE noise_compute");
        PyrNoise.Border = I_CONT;

        NoiseImage = im_noise (1., Nl, Nc);
        PyrNoise.transform (NoiseImage);
        for (i = 0; i < Nbr_Plan; i++) 
	                 TabSignificantLevel[i] = sigma (PyrNoise.scale(i));
        for (i = Nbr_Plan; i < MAX_SCALE; i++)  TabSignificantLevel[i] = 0.;
    }
}

/****************************************************************************/

float mr_noise_estimation (MultiResol &MR_Data)
{
    float Noise=0.;

    if (TabSignificantLevel[0] < FLOAT_EPSILON) noise_compute (MR_Data);
    switch (MR_Data.Set_Transform)
    {
        case TRANSF_PAVE:
        case TRANSF_PYR:
          if (   (MR_Data.Type_Transform == TM_PYR_MINMAX)
               || (MR_Data.Type_Transform == TM_PAVE_MINMAX)
               || (MR_Data.Type_Transform == TM_PYR_SCALING_FUNCTION))
               Noise = detect_noise_sigma (MR_Data.scale(0), True);
          else Noise = detect_noise_sigma (MR_Data.scale(0), False);
          Noise /= TabSignificantLevel [0];
          break;
        case TRANSF_MALLAT:
          Noise=detect_noise_sigma(MR_Data.extract_scale(0,D_DIAGONAL),False);
          Noise /= TabSignificantLevel [0];
          break;
        case TRANSF_FEAUVEAU:
          Noise = detect_noise_sigma (MR_Data.extract_scale(0, D_RESOL), False);
          Noise /= TabSignificantLevel [1];
          break;
        default:
          fprintf (stderr, "Error in noise_mr_estimation: bad Set_Transform");
          break; 
    }
#if PRINT_DATA
fprintf (stderr, "\n\nNOISE Standard Deviation Estimated = %f\n", Noise);
#endif
     return (Noise);
}

/****************************************************************************/

Bool mallat_significant (MultiResol &MR_Data, int s, int i, int j,
                         float Level, details Det)
{
    if (ABS (MR_Data(s,i,j, Det)) < Level)  return False;
    else return True;
}

/****************************************************************************/

float mr_tab_noise(int s)
// ?? Kept for compatibity: see MultiResol.scale_norm
{
    return TabSignificantLevel[s];
}

/****************************************************************************/

float mr_level_noise (MultiResol &MR_Data, int s, float N_Sigma)
// ?? Kept for compatibity: see MultiResol.gauss_detect_level 
{
    float L;

    switch (MR_Data.Set_Transform)
    {
        case TRANSF_PAVE:
        case TRANSF_PYR:
           if (s == 0) L = (N_Sigma+1) * TabSignificantLevel[s];
           else L = N_Sigma * TabSignificantLevel[s];
           if ( (MR_Data.Type_Transform == TM_PAVE_MINMAX) 
                  || (MR_Data.Type_Transform == TM_PYR_MINMAX)
                  || (MR_Data.Type_Transform == TM_PYR_SCALING_FUNCTION))
                                                        L *= 2.;
           break;
        case TRANSF_MALLAT:
           if (s == 0)  L = (N_Sigma+1) * TabSignificantLevel[s];
           else L = N_Sigma * TabSignificantLevel[s];
           break;
        case TRANSF_FEAUVEAU:
           if (s == 0) L = (N_Sigma+1) * TabSignificantLevel[s];
           else L = N_Sigma * TabSignificantLevel[s];
           break;
        default: L = 0.;break;
    }
    return (L);
}

/****************************************************************************/

void mr_threshold (MultiResol &MR_Data, float &Noise_Ima, float N_Sigma)
{
    int s;
    float Level;
    int i,j;
    int Nbr_Plan = MR_Data.nbr_scale()-1;
    int Nl,Nc;

    /* Noise behaviour in the multiresolution space */
    noise_compute (MR_Data);

    /* Noise estimation in the Data */
    if (Noise_Ima < FLOAT_EPSILON) Noise_Ima = mr_noise_estimation (MR_Data);

    switch (MR_Data.Set_Transform)
    {
        case TRANSF_PAVE:
        case TRANSF_PYR:
          for (s = 0; s < Nbr_Plan; s++)
          {
              Level = Noise_Ima * mr_level_noise (MR_Data,  s,  N_Sigma);
              noise_threshold_imag (MR_Data.scale(s),  Level);
          }
          break;
        case TRANSF_MALLAT:
          for (s = 0; s < Nbr_Plan; s++)
          {
              Level = Noise_Ima * mr_level_noise (MR_Data,  s,  N_Sigma);
              Nl = MR_Data.size_scale_nl(s,D_HORIZONTAL);
              Nc = MR_Data.size_scale_nc(s,D_HORIZONTAL);
              for (i = 0; i < Nl; i++)
              for (j = 0; j < Nc; j++)
              {
                 if (!mallat_significant(MR_Data, s, i, j, Level, D_NULL))
                 {
                     MR_Data(s,i,j, D_HORIZONTAL) = 0.;
                     MR_Data(s,i,j, D_DIAGONAL) = 0.;
                     MR_Data(s,i,j, D_VERTICAL) = 0.;
                 }
              }
          }
          break;
        case TRANSF_FEAUVEAU:
          for (s = 0; s < Nbr_Plan; s++)
          {
              Level = Noise_Ima * mr_level_noise (MR_Data,  2*s,  N_Sigma);
              for (i = 0; i < MR_Data.size_scale_nl(s,D_HALF_RESOL); i++)
              for (j = 0; j < MR_Data.size_scale_nc(s,D_HALF_RESOL); j++)
                 if (fabs (MR_Data(s,i,j, D_HALF_RESOL)) < Level)
                                      MR_Data(s,i,j, D_HALF_RESOL) = 0.;
              Level = Noise_Ima  * mr_level_noise (MR_Data,  2*s+1,  N_Sigma);
              for (i = 0; i < MR_Data.size_scale_nl(s,D_RESOL); i++)
              for (j = 0; j < MR_Data.size_scale_nc(s,D_HALF_RESOL); j++)
                 if (fabs (MR_Data(s,i,j, D_RESOL)) < Level)
                                      MR_Data(s,i,j, D_RESOL) = 0.;
          }
          break;
        default:
          fprintf (stderr, "Error in noise_mr_threshold: bad Set_Transform");
          break; 
    }
}

/****************************************************************************/

void mr_support_threshold (MultiResol &MR_Data, float Noise_Ima, 
                            float N_Sigma, Bool SetSupport)
{
    int s;
    float Level, Val;
    int i,j;
    int Nbr_Plan = MR_Data.nbr_scale()-1;
    int Nl = MR_Data.scale(0).nl();
    int Nc = MR_Data.scale(0).nc();

    switch (MR_Data.Set_Transform)
    {
        case TRANSF_PAVE:
        case TRANSF_PYR:
          for (s = 0; s < Nbr_Plan; s++)
          {
              Nl = MR_Data.scale(s).nl();
              Nc = MR_Data.scale(s).nc();
              Level = Noise_Ima * mr_level_noise (MR_Data,  s,  N_Sigma);
              Val = mr_level_noise (MR_Data,  s,  N_Sigma);
              Level = Noise_Ima * Val;
#if PRINT_DATA
        cout << endl;
        cout << "Scale " << s+1 << " Noise_Ima = " << Noise_Ima << endl;
        cout << " mr_level_noise = " << Val << endl;
        cout << " N_Sigma = " << N_Sigma << endl;
        cout << " Level = " << Level << endl;
        cout << " Sigma = " << sigma(MR_Data.scale(s)) << endl;
#endif
              noise_support_threshold_imag (MR_Data.scale(s), Level, 
                                            DataSupport.scale(s), SetSupport);
          }
          break;
        case TRANSF_MALLAT:
          for (s = 0; s < Nbr_Plan; s++)
          {
              Level = Noise_Ima * mr_level_noise (MR_Data,  s,  N_Sigma);
              Nl = MR_Data.size_scale_nl(s,D_HORIZONTAL);
              Nc = MR_Data.size_scale_nc(s,D_HORIZONTAL);
              for (i = 0; i < Nl; i++)
              for (j = 0; j < Nc; j++)
              {
                 if (!mallat_significant(MR_Data, s, i, j, Level,  D_HORIZONTAL))
                 {
                     if (DataSupport(s,i,j, D_HORIZONTAL) < FLOAT_EPSILON)
                                           MR_Data(s,i,j, D_HORIZONTAL) = 0.;
                 }
                 else  if (SetSupport) DataSupport(s,i,j, D_HORIZONTAL) = 1.;
                 
                 
                 if (!mallat_significant(MR_Data, s, i, j, Level,  D_VERTICAL))
                 {
                     if (DataSupport(s,i,j, D_VERTICAL) < FLOAT_EPSILON)
                                           MR_Data(s,i,j, D_VERTICAL) = 0.;
                 }
                 else  if (SetSupport) DataSupport(s,i,j, D_VERTICAL) = 1.;
                 
                 
                 if (!mallat_significant(MR_Data, s, i, j, Level,  D_DIAGONAL))
                 {
                     if (DataSupport(s,i,j, D_DIAGONAL ) < FLOAT_EPSILON)
                                           MR_Data(s,i,j,  D_DIAGONAL) = 0.;
                 }
                 else  if (SetSupport) DataSupport(s,i,j, D_DIAGONAL) = 1.;
              }
          }
          break;
        case TRANSF_FEAUVEAU:
          for (s = 0; s < Nbr_Plan; s++)
          {
              Level = Noise_Ima * mr_level_noise (MR_Data,  2*s,  N_Sigma);
              Nl = MR_Data.size_scale_nl(s,D_HALF_RESOL);
              Nc = MR_Data.size_scale_nc(s,D_HALF_RESOL);
              for (i = 0; i < Nl; i++)
              for (j = 0; j < Nc; j++)
              {
                 if (fabs (MR_Data(s,i,j, D_HALF_RESOL)) < Level)
                 {
                     if (DataSupport(s,i,j, D_HALF_RESOL) < FLOAT_EPSILON)
                         MR_Data(s,i,j, D_HALF_RESOL) = 0.;
                 }
                 else   if (SetSupport)
                             DataSupport(s,i,j, D_HALF_RESOL) = 1.;
              }
              Level = Noise_Ima * mr_level_noise (MR_Data,  2*s+1,  N_Sigma);
              Nl = MR_Data.size_scale_nl(s,D_RESOL);
              Nc = MR_Data.size_scale_nc(s,D_RESOL);
              for (i = 0; i < Nl; i++)
              for (j = 0; j < Nc; j++)
              {
                 if (fabs (MR_Data(s,i,j, D_RESOL)) < Level)
                 {
                     if (DataSupport(s,i,j, D_RESOL) < FLOAT_EPSILON)
                     {
                         MR_Data(s,i,j, D_RESOL) = 0.;
                     }
                 }
                 else  if (SetSupport) DataSupport(s,i,j, D_RESOL) = 1.;
              }
          }
          break;
        default:
         fprintf (stderr,"Error in noise_mr_iter_threshold: bad Set_Transform");
         exit (-1);
         break; 
    }
}


/****************************************************************************/
