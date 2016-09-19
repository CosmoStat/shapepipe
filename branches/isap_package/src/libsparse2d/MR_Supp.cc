/*******************************************************************************
**
**    UNIT
**
**    Version: 3.1
**
**    Author: Jean-Luc Starck
**
**    Date:  96/05/02 
**    
**    File:  MR_Supp.cc
**
*******************************************************************************
**
**    DESCRIPTION  Multiresolution support procedures
**    -----------  
**                 
**
*******************************************************************************
**
** void mr_create_support (MultiResol &MR_Data, float &Noise_Ima, float N_Sigma)
**
** create a multiresolution support associate to a multiresolution
** transform. The multiresolution support is contained in the global
** varaible DataSupport which is a multiresolution object
** if Noise_Ima <= 0, estimate the noise from the first scale
**
*******************************************************************************
**
** void mr_dilate_support (Bool SmoothSupport)
** 
** dilate each scale of the support DataSupport
** if SmoothSupport, smooth the dilated scale
**
*******************************************************************************
**
** void mr_smooth_support ()
**
** smooth each scale of the multiresolution support DataSupport
**
*******************************************************************************
** 
** void mr_hierar_support ()
**
** Hierarchical model of the first order
** put zero at positions where the model is not OK
**
*******************************************************************************
**
** void mr_add_map_support (Ifloat &MapPos)
**
** Add a Boolean image to the first scale of the multiresolution support
** The global variable DataSupport (type MultiResol) is modified
** 
*******************************************************************************
** void mr_support (Ifloat &Imag, float &Noise_Ima, type_transform Transform, 
**                 int Nbr_Plan, float N_Sigma, type_noise Stat_Noise)
**
** Multiresolution support creation from an image
** If Noise_Ima <= 0, then the noise in the image is estimated
**
** Imag = IN: image
** Noise_Ima: Standard deviation of the noise.
**              if Sigma_Noise = 0, the standart deviation of the noise
**              is estimated in this routine from an automatically way
** Type_Transform: multiresolution transform algorithm choosen
** Nbr_Plan: number of scales 
**           it is not necessary to take values for Nbr_Plan superior to 5
** Stat_Noise: defines the statistic of the noise
**           Stat_Noise = NOISE_GAUSSIAN, NOISE_POISSON, ...
** N_Sigma: Used for detection level, the level is estimated by:
**          Level = N_Sigma * Noise_Standart_Deviation
**          N_Sigma = 3 is a standart value
**
**
** The global variable DataSupport (type MultiResol) is set
**
******************************************************************************/

// static char sccsid[] = "@(#)MR_Supp.cc 3.1 96/05/02 CEA 1995 @(#)";

#include "MR_Support.h"
#include "MR_Noise.h"
#include "MR_Obj.h"

MultiResol DataSupport;
extern float TabSignificantLevel[MAX_BAND];

/****************************************************************************/

static void support_set (MultiResol &MR_Data, float Noise_Ima, float N_Sigma)
{
    int s;
    float Level;
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
              for (i = 0; i < MR_Data.scale(s).nl(); i++)
              for (j = 0; j < MR_Data.scale(s).nc(); j++)
                 if (ABS (MR_Data(s,i,j)) >= Level) DataSupport(s,i,j) = 1.;
                 else DataSupport(s,i,j) = 0.;
          }
          DataSupport.scale(Nbr_Plan).init(1);
          break;
        case TRANSF_MALLAT:
          for (s = 0; s < Nbr_Plan; s++)
          {
              Level = Noise_Ima * mr_level_noise (MR_Data,  s,  N_Sigma);
              Nl = MR_Data.size_scale_nl(s,D_HORIZONTAL);
              Nc = MR_Data.size_scale_nc(s);
              for (i = 0; i < Nl; i++)
              for (j = 0; j < Nc; j++)
              {
                 if (ABS (MR_Data(s, i, j, D_HORIZONTAL)) >= Level) 
                                MR_Data(s, i, j, D_HORIZONTAL) = 1;
                 else MR_Data(s, i, j, D_HORIZONTAL) = 0;
                 if (ABS (MR_Data(s, i, j, D_DIAGONAL)) >= Level) 
                                MR_Data(s, i, j, D_DIAGONAL) = 1;
                 else MR_Data(s, i, j, D_DIAGONAL) = 0;
                 if (ABS (MR_Data(s, i, j, D_VERTICAL)) >= Level) 
                                MR_Data(s, i, j, D_VERTICAL) = 1;
                 else MR_Data(s, i, j, D_VERTICAL) = 0;
              }
          }
          break;
        case TRANSF_FEAUVEAU:
          for (s = 0; s < Nbr_Plan; s++)
          {
              Level = Noise_Ima * mr_level_noise (MR_Data,  2*s,  N_Sigma);
              Nl = MR_Data.size_scale_nl(s,D_HALF_RESOL);
              Nc = MR_Data.size_scale_nc(s);
              for (i = 0; i < Nl; i++)
              for (j = 0; j < Nc; j++)
                 if (ABS (MR_Data(s,i,j, D_HALF_RESOL)) >= Level)
                                  DataSupport(s,i,j, D_HALF_RESOL) = 1;
              Level = Noise_Ima * mr_level_noise (MR_Data,  2*s+1,  N_Sigma);
              Nl = MR_Data.size_scale_nl(s,D_RESOL);
              Nc = MR_Data.size_scale_nc(s);
              for (i = 0; i < Nl; i++)
              for (j = 0; j < Nc; j++)
                 if (ABS (MR_Data(s,i,j, D_RESOL)) >= Level)
                                  DataSupport(s,i,j, D_RESOL) = 1;
          }
          break;
        default:
         fprintf (stderr,"Error in support_set: bad Set_Transform");
         exit (-1);
         break; 
    }
}

/****************************************************************************/

void mr_create_support (MultiResol &MR_Data, float &Noise_Ima, float N_Sigma)
{
    int Nl = MR_Data.scale(0).nl();
    int Nc = MR_Data.scale(0).nc();
    int Nbr_Plan = MR_Data.nbr_scale();

    /* Allocate space for the support */
    DataSupport.alloc (Nl, Nc, Nbr_Plan, MR_Data.Type_Transform, "Support");

    /* Noise behaviour in the multiresolution space */
    noise_compute (MR_Data);

    /* Noise estimation in the Data */
    if (Noise_Ima < FLOAT_EPSILON) 
                                Noise_Ima = mr_noise_estimation (MR_Data);

    /* initialize the support */
    support_set (MR_Data, Noise_Ima, N_Sigma);
}

/****************************************************************************/


void mr_dilate_support (Bool SmoothSupport)
{
    int Window, Size=1;
    int s=0;

    int Nl = DataSupport.scale(0).nl();
    int Nc = DataSupport.scale(0).nc();
    int Nbr_Plan = DataSupport.nbr_scale();
    Ifloat Buff (Nl, Nc, "Buff Dilate Support");

    switch (DataSupport.Set_Transform)
    {
        case TRANSF_PAVE:
        case TRANSF_PYR:
          for (s = 0; s < Nbr_Plan-1; s++)
          {
              Window = 2*Size+1;
              Nl = DataSupport.scale(s).nl();
              Nc = DataSupport.scale(s).nc();
              Buff.resize (Nl,Nc);
              if (!SmoothSupport)
              {
                  morpho_dilation (DataSupport.scale(s), Buff, Size);
                  DataSupport.scale(s) = Buff; 
              }
              else
              {
                  morpho_dilation (DataSupport.scale(s), Buff, Size); 
                  smooth_bspline (Buff, DataSupport.scale(s));
              }
              if (DataSupport.Set_Transform == TRANSF_PAVE) Size = 2*Size;
          }
          break;
        case TRANSF_MALLAT:
        case TRANSF_FEAUVEAU:
          morpho_dilation (DataSupport.scale(0), Buff);
          DataSupport.scale(0) = Buff;
          break; 
        default:
         fprintf (stderr,"Error in mr_dilate_support: bad Set_Transform");
         exit (-1);
         break; 
    }
}

/****************************************************************************/

void mr_smooth_support ()
{
    int Nl = DataSupport.scale(0).nl();
    int Nc = DataSupport.scale(0).nc();
    int Nbr_Plan = DataSupport.nbr_scale();
    Ifloat Buff (Nl, Nc, "Buff Dilate Support");
    int s = 0;

    switch (DataSupport.Set_Transform)
    {
        case TRANSF_PAVE:
        case TRANSF_PYR:
          for (s = 0; s < Nbr_Plan-1; s++)
          {
              Nl = DataSupport.scale(s).nl();
              Nc = DataSupport.scale(s).nc();
              Buff.resize (Nl,Nc);
              smooth_bspline (DataSupport.scale(s), Buff);
              DataSupport.scale(s) = Buff;
          }
          break;
        case TRANSF_MALLAT:
        case TRANSF_FEAUVEAU:
          smooth_bspline (DataSupport.scale(0), Buff);
          DataSupport.scale(0) = Buff;
          break; 
        default:
         fprintf (stderr,"Error in mr_smooth_support: bad Set_Transform");
         exit (-1);
         break; 
    }
}

/****************************************************************************/

void mr_hierar_support ()
{
    int i,j;
    Bool PutZero;

    for (int s = 0; s < DataSupport.nbr_scale()-2; s++)
     switch (DataSupport.Set_Transform)
     {
         case TRANSF_PAVE:
              for (i = 1; i < DataSupport.scale(s).nl()-1; i++)
              for (j = 1; j < DataSupport.scale(s).nc()-1; j++)
              {
                  PutZero = True;
                  if (DataSupport.scale(s)(i,j) == 1)
                  {
                    if (DataSupport.scale(s)(i-1,j) == 1) PutZero = False;
                    else if (DataSupport.scale(s)(i+1,j) == 1) PutZero = False;
                    else if (DataSupport.scale(s)(i,j+1) == 1) PutZero = False;
                    else if (DataSupport.scale(s)(i,j-1) == 1) PutZero = False;
                    else if (DataSupport.scale(s+1)(i,j) == 1) PutZero = False;
                  }
                  if (PutZero) DataSupport.scale(s)(i,j) = 0;
              }
              break;
         case TRANSF_PYR:
              for (i = 1; i < DataSupport.scale(s).nl()-1; i++)
              for (j = 1; j < DataSupport.scale(s).nc()-1; j++)
              {
                  PutZero = True;
                  if (DataSupport.scale(s)(i,j) == 1)
                  {
                    if (DataSupport.scale(s)(i-1,j) == 1) PutZero = False;
                    else if (DataSupport.scale(s)(i+1,j) == 1) PutZero = False;
                    else if (DataSupport.scale(s)(i,j+1) == 1) PutZero = False;
                    else if (DataSupport.scale(s)(i,j-1) == 1) PutZero = False;
                    else if (DataSupport.scale(s+1)(i/2,j/2)==1) PutZero= False;
                  }
                  if (PutZero) DataSupport.scale(s)(i,j) = 0;
              }
              break;
         default: break;
     }
}

/************************************************************************/

void mr_add_map_support (Ifloat &MapPos)
{
    int i,j;
    int Nl = MapPos.nl();
    int Nc = MapPos.nc();

    for (i = 0; i < Nl; i++)
    for (j = 1; j < Nc; j++)   
    { 
       /* Modification of the first scale of the
       multiresolution support */
       if (MapPos(i,j) > 0.5) DataSupport(0,i,j) = 1.;

       /* Modification of others scales */
/*       for (s = 1; s < DataSupport.nbr_scale()-1; s++)
          switch (DataSupport.Set_Transform)
          {
            case TRANSF_PAVE:
               Step = pow2(s);
               for (is = Step; is < Nl-Step; is++)
               for (js = Step; js < Nc-Step; js++)
                  if (MapPos(is,js) > 0.5) DataSupport(s,i,j) = 1.;
               break;
            case TRANSF_PYR:
               Step = pow2(s);
               for (is = Step; is < Nl-Step; is++)
               for (js = Step; js < Nc-Step; js++)
                  if (MapPos(is,js) > 0.5) DataSupport(s,i/Step,j/Step) = 1.;
               break;
          }
*/
    }     
}

/************************************************************************/

void mr_support (Ifloat &Imag, float &Noise_Ima, type_transform Transform, 
                 int Nbr_Plan, float N_Sigma, type_noise Stat_Noise)
{
    int Nl = Imag.nl();
    int Nc = Imag.nc();
    Ifloat ImagGauss (Nl, Nc, "ImagGauss");


    /* Allocate space for the support */
    DataSupport.alloc (Nl, Nc, Nbr_Plan, Transform, "Support");

    /* Noise behaviour in the multiresolution space */
    noise_compute (Nbr_Plan, Transform, Nl, Nc);

    switch (Stat_Noise)
    {
        case NOISE_POISSON: 
        case NOISE_GAUSS_POISSON:
                noise_poisson_transform (Imag, ImagGauss);
                Noise_Ima = 1.;

                /* multiresolution transform */
                DataSupport.transform(ImagGauss);

               /* Noise estimation in the Data */
               if (Noise_Ima < FLOAT_EPSILON) 
                                Noise_Ima = mr_noise_estimation (DataSupport);

               /* initialize the support */
               support_set (DataSupport, Noise_Ima, N_Sigma);
               break;
        case NOISE_GAUSSIAN: 
                /* multiresolution transform */
                DataSupport.transform(Imag);

               /* Noise estimation in the Data */
               if (Noise_Ima < FLOAT_EPSILON) 
                                Noise_Ima = mr_noise_estimation (DataSupport);

               /* initialize the support */
               support_set (DataSupport, Noise_Ima, N_Sigma);
               break;
        default:
               cout << "Error: not implemented in this routine... " << endl;
               exit(0);
               break;
    }
}

/****************************************************************************/
