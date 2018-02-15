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
**    File:  MR1D_Filtering.cc
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

#include "Array.h"
#include "IM_Obj.h"
#include "IM_Math.h"
#include "IM_Noise.h"
#include "MR1D_Obj.h"
#include "NR.h"

// RG #include "MR_Sigma.h"
#include "MR1D_Filter.h"   // RG for mr1d_support_threshold 

MR_1D Sup1d;

#define PRINT_DATA 0
#define DEBUG_DATA 0

// see also  NbrScale and Nx


void mr1d_support_set (MR_1D &MR_Data, float Noise_Ima, float N_Sigma, float N_Sigma0)
{
    int s;
    float Level;
    int i;
    int Nbr_Plan = MR_Data.nbr_scale()-1;

    switch (MR_Data.Set_Transform)
    {
        case TRANS1_PAVE:
        case TRANS1_PYR:
          for (s = 0; s < Nbr_Plan; s++)
          {
              Level = Noise_Ima * mr1d_level_noise (MR_Data,  s,  N_Sigma, N_Sigma0);
              for (i = 0; i < MR_Data.size_scale_np(s); i++)
              {
                 if (ABS (MR_Data(s,i)) >= Level) Sup1d(s,i) = 1.;
                 else Sup1d(s,i) = 0.;
              }
          }
          break;
        default:
         fprintf (stderr,"Error in support_set: bad Set_Transform");
         exit (-1);
         break; 
    }
}

/****************************************************************************/

void mr1d_create_support (MR_1D &MR_Data, float Noise_Ima, float N_Sigma, float N_Sigma0)
{
    int Nc = MR_Data.size_ima_np();
    int Nbr_Plan = MR_Data.nbr_scale();

    /* Allocate space for the support */
    Sup1d.alloc (Nc, MR_Data.Type_Transform, "Support", Nbr_Plan);

    /* Noise behaviour in the multiresolution space */
    mr1d_noise_compute ( Nbr_Plan, MR_Data.Type_Transform);

    /* initialize the support */
    mr1d_support_set (MR_Data, Noise_Ima, N_Sigma, N_Sigma0);
}

/****************************************************************************/

void mr1d_support_threshold (MR_1D &MR_Data, float Noise_Ima, float N_Sigma, 
                             float N_Sigma0)
// in MR1D_Filter.h  N_Sigma0=DEF_NSIG0=4.0
{
    int s,i,cpt;
    int Nbr_Plan = MR_Data.nbr_scale()-1;
    float Level;

    switch (MR_Data.Set_Transform)
    {
        case TRANS1_PAVE:
        case TRANS1_PYR:
          for (s = 0; s < Nbr_Plan; s++)
          {
              Level = Noise_Ima * mr1d_level_noise (MR_Data,s,N_Sigma,N_Sigma0);
              cpt =0;
              for (i = 0; i < MR_Data.size_scale_np(s); i++)
              {
                 if (Sup1d (s,i) < FLOAT_EPSILON) 
                 {
                    if (ABS (MR_Data(s,i)) >= Level) Sup1d(s,i) = 1.;
                    else 
                    {
                       MR_Data(s,i) = 0.;
#if DEBUG_DATA
                       cpt++;
#endif
                    }
                 }
              }
#if DEBUG_DATA
cout << "Scale " << s+1 << " Level = " << Level << " T = " << cpt / (float) MR_Data.size_scale_np(s) * 100. << endl;
#endif

          }
          break;
        default:
         fprintf (stderr,"Error in mr1d_support_threshold: bad Set_Transform");
         exit (-1);
         break; 
    }
}


/****************************************************************************/

void mr1d_iter_threshold (fltarray &Data, fltarray &Result, float &Noise_Ima,
                   type_trans_1d Transform, int Nbr_Plan, 
                   float N_Sigma, float Epsilon, int Max_Iter)
{
   int i,Iter = 0;
   int N = Data.axis(1);
   MR_1D MR_Data (N, Transform, "MR_Filtering", Nbr_Plan);
   fltarray Resi (N);
   float Sigma;
   float Old_Sigma, Delta;
   float Energ;

   Energ = Data.sigma();

    /* Initialization */
    Delta = Sigma = 1.e10;

   while (True)
   {
       Old_Sigma = Sigma; 
       Resi = Data - Result;

       Sigma = Resi.sigma();
       MR_Data.transform (Resi);
       mr1d_support_threshold (MR_Data, Noise_Ima, N_Sigma);
       MR_Data.recons (Resi);
       Delta = (Old_Sigma - Sigma)  / Energ;

       Result += Resi;

       /* Threshold negative values */
       for (i = 0; i < N; i++) if (Result(i) < 0.) Result(i) = 0.;

#if PRINT_DATA
       cout << "Iter " << Iter+1 << " ==> Delta  = " << ABS(Delta) << endl;
#endif

       Iter ++;
     if ((Iter >= Max_Iter) || (Delta < Epsilon)) break;
   }
}

/****************************************************************************/

void mr1d_filter (fltarray &Data, fltarray &Result, float &Sigma_Noise,
                     type_trans_1d Transform, int Nbr_Plan,
                     type_noise Stat_Noise, float N_Sigma, 
                     float Epsilon, int Max_Iter, Bool WriteFilterMr)
{
    int N = Data.axis(1);
    fltarray DataGauss (N);
    float MinDat;
    MR_1D MR_Data (N, Transform, "MR_Filtering", Nbr_Plan);

//    MinDat = min(Data);
    MinDat = Data.min(); 

    if (Stat_Noise == NOISE_POISSON)
    {
        noise_poisson_transform (Data, DataGauss);
        Sigma_Noise = 1.;
    }
    else 
    {
        /* we want to use the positivity constraint, so we set 
           the min at zero
        */
        /* for ( j = 0; j < N; j++) DataGauss(j) = Data(j) - MinDat;*/
        DataGauss = Data;
    }

    /* Noise estimation in the Data */
    if (Sigma_Noise < FLOAT_EPSILON) 
                            Sigma_Noise = detect1d_noise_from_med (DataGauss);
#if PRINT_DATA
    cout << "Noise estimation : " << Sigma_Noise << endl;
#endif

   MR_Data.transform (DataGauss);
   mr1d_create_support (MR_Data, Sigma_Noise, N_Sigma);
   mr1d_support_threshold (MR_Data, Sigma_Noise, N_Sigma);
   
   if (WriteFilterMr) {
      fltarray Dat = MR_Data.image();
      fits_write_fltarr("FiterMr1d.mr", Dat);
   }
   
   
   MR_Data.recons (Result);

    mr1d_iter_threshold (DataGauss, Result, Sigma_Noise, Transform, 
                         Nbr_Plan,  N_Sigma,  Epsilon,  Max_Iter);
    if (Stat_Noise == NOISE_POISSON)
                     noise_inverse_poisson_transform (Result, Result);
 /*   else
      for ( j = 0; j < N; j++)    Result(j) += MinDat; */
}


/***************************************************************************/
 
void MR1DFiltering::threshold(fltarray &Data, float T)
{
   for (int i=0; i < Data.n_elem(); i++) if (Data(i) < T) Data(i) = 0.;
}

/********************************************************************/

void MR1DFiltering::MaxThreshold(fltarray &Data, float T)
{
   for (int i=0; i < Data.n_elem(); i++) if (Data(i) > T) Data(i) = T;
}
 
/********************************************************************/ 

void MR1DFiltering::kill_last_scale(MR_1D & MR_Data)
{
   int s = MR_Data.nbr_band()-1;
   for(int w=0; w < MR_Data.size_scale_np(s); w++) MR_Data(s,w) = 0;
}

/********************************************************************/ 
 
void MR1DFiltering::filter(fltarray &Imag, fltarray &Result)
{
    int b;

    if (NoiseModel == NULL)
    {
       cerr << "Error: Noise model not initialized ... " << endl;
       exit(-1);
    }
   
   int Np = Imag.nx();
   int NbrPlan =  NoiseModel->nbr_scale(); 
   MR_1D MR_Data;
   MR_Data.alloc(Np,NoiseModel->type_trans(), NbrPlan,
                 NoiseModel->filter_bank(), NoiseModel->TypeNorm);
              //   cout << "SigmaNoise1 = " <<  NoiseModel->SigmaNoise << endl;

   if (IsDataModelled == False)  NoiseModel->model(Imag, MR_Data);
              // cout << "SigmaNoise2 = " <<  NoiseModel->SigmaNoise << endl;

   if (NoiseModel->TransImag==True)  NoiseModel->signal_transform(Imag);

   // in this case, WT of image has not yet been computed                                              
   if (IsDataModelled == True)  MR_Data.transform(Imag);
   
   if (WriteFilterMR == True)
   {
      MR_1D MR_Sup;
      MR_Sup.alloc(Np,NoiseModel->type_trans(), NbrPlan,
                 NoiseModel->filter_bank(), NoiseModel->TypeNorm);
      if (NameSupport == NULL)
      {
         cout << "Error: cannot write the multiresolution support ... " << endl;
         exit(-1);
      }
      MR_Sup.transform(Imag);
      NoiseModel->threshold(MR_Sup);
      fltarray Dat = MR_Sup.image();
      fits_write_fltarr(NameSupport, Dat);
   }
    
      
   if ((Verbose == True)  && (StatNoise  == NOISE_GAUSSIAN))
                 cout << "SigmaNoise = " <<  NoiseModel->SigmaNoise << endl;

   // if  TransImag = True ==> wavelet coefficient are not those of the 
   // input image (but the transformed image). We have transform again the image.       
   // if (NoiseModel->TransImag==True)  NoiseModel->signal_transform(Imag);
                                                 
   // in this case, WT of image has not yet been computed                                              
   //if (IsDataModelled == True)  MR_Data.transform(Imag);

   if (Verbose == True)
   {
      for (b=0; b < MR_Data.nbr_band()-1; b++)
      {
         float NSig = (NoiseModel->NSigma)[b];
         cout << "Band " << b+1 << " Nsig = " << NSig << " Sigma = " <<  NoiseModel->sigma(b,0)  << " T  = " <<  NoiseModel->sigma(b,0) * NSig  <<  endl;
      }
   }
      
   switch (Filter)
   {
        case FIL_1D_ITER_THRESHOLD:
	case FIL_1D_TV_CONSTRAINT:
               mr_support_filter(Imag, MR_Data, Result);   
              break;
        case FIL_1D_THRESHOLD:
              NoiseModel->threshold(MR_Data);
              if (KillLastScale == True) kill_last_scale(MR_Data);
              MR_Data.recons(Result);
              break;
        case FIL_1D_SOFT_THRESHOLD:
	     for (b=0; b < MR_Data.nbr_band()-1; b++)
             {        
                for (int i=0;i<MR_Data.size_scale_np(b);i++)
                {
                    float NoiseLevel = (NoiseModel->NSigma)[b]*NoiseModel->sigma(b,i);
                    MR_Data(b,i) = soft_threshold(MR_Data(b,i),NoiseLevel);
                 }
              }
              if (KillLastScale == True) kill_last_scale(MR_Data);
              MR_Data.recons(Result);
              break;
        default:            
          cerr << "Error: this filtering method is not implemented in this routine ... " << endl;
          exit(-1);
     }
   
   // if the data have been transform, we have to take the inverse transform
   // of the solution
   if (NoiseModel->TransImag == True)
   { 
       if (KillLastScale != True)
       {
          NoiseModel->signal_invtransform(Imag);
          NoiseModel->signal_invtransform(Result);
       }
       else
       {
           Result = Imag - Result;
           NoiseModel->signal_invtransform(Imag);
           NoiseModel->signal_invtransform(Result);
           Result = Imag - Result;
       }
    }
    if (PositivIma == True) threshold(Result);
    if (MaxIma == True) MaxThreshold(Result);
}

/********************************************************************/ 
 
inline int sgn(float Val) {return ( (Val >= 0) ? 1: -1);}

inline float markov_val2(fltarray &Sig, int i)
{
      float MarkovPowerParam=1.1;
      double p = MarkovPowerParam;
      double dp=p-1;
      int N = 2;
      dblarray Tab(N);
      int ip = (i == Sig.nx() -1) ? Sig.nx() -2: i+1;
      int im = (i == 0) ? 1: i-1;

      Tab(0) = Sig(i) - Sig(ip);
      Tab(1) = Sig(i) - Sig(im);
      double  Val = 0;
      for (int i=0; i < N; i++) Val += p*sgn(Tab(i))*pow(ABS(Tab(i)),dp);
      return (float) Val / (float) N;
}

/********************************************************************/ 

static void regul1d(fltarray &Result, fltarray &Resi, float RegulParam)
{
   int N = Result.nx();
   int i;
   for (i = 0; i < N; i ++) Resi(i) -=  RegulParam*markov_val2(Result,i); 
}

/********************************************************************/ 

void MR1DFiltering::mr_support_filter(fltarray &Imag, MR_1D & MR_Data, fltarray &Result)
{
   int Iter = 0;
   int i,Np = Imag.nx();
   fltarray Resi (Np);
   float Sigma, Noise_Ima;
   float Old_Sigma, Delta;
   float  Conv = 1.;
   if (RegulParam > 0.5) Conv = 1. / (2.*RegulParam);
   if  (MR_Data.Type_Transform == TO1_PAVE_B3SPLINE_GEN2)  UseAdjointRec = False;

   Delta = Sigma = 1.e9;
   Noise_Ima=NoiseModel->SigmaNoise;
   if (Noise_Ima < FLOAT_EPSILON) Noise_Ima = 1.;
   // if ((NoiseModel->TransImag == True) && (NoiseModel->which_noise() != NOISE_EVENT_POISSON)) Sup_Set = False;
   Result.init();
   Resi = Imag;

   Sigma = Resi.energy();
   Old_Sigma = Sigma; 
   while (True)
   {
        MR_Data.transform (Resi, Border);
        NoiseModel->threshold(MR_Data, False, Sup_Set);
         
        if (KillLastScale == True) kill_last_scale(MR_Data);
        if (UseAdjointRec == True) MR_Data.rec_adjoint (Resi);
        else MR_Data.recons (Resi);
	if (RegulParam > 0) 
	{
	   regul1d(Result, Resi, RegulParam);
           for (i = 0; i < Np; i ++) 
	            Result(i) += Conv*Resi(i);
	}
	else Result += Resi;	   
        if (PositivIma == True) threshold(Result);
        if (MaxIma == True) MaxThreshold(Result);

        Old_Sigma = Sigma; 
        Resi = Imag - Result;
        Sigma = Resi.energy();
        Delta = (Old_Sigma - Sigma)  / (Np*Noise_Ima*Noise_Ima);

        Iter ++;
        float Xi = Sigma / (Np*Noise_Ima*Noise_Ima);
        if (Verbose == True)
        {
           if (((NoiseModel->which_noise() == NOISE_GAUSSIAN) ||
               (NoiseModel->which_noise() == NOISE_POISSON)) &&
               (KillLastScale == False))
             cout << "Iter "<< Iter << ": Xi = "<< Xi <<" ==> Delta  = " <<Delta << endl;
           else cout << "Iter "<< Iter  <<" ==> Delta  = " <<Delta<< ", Regul =  " << RegulParam <<  endl;
        }
        if ((Iter >= Max_Iter) || (Delta < Epsilon)) break;
    }
}

/********************************************************************/ 

// void MR1DFiltering::grad_adj_filter(MR_1D &MR_Data, fltarray &Result)
// {
//    int i,Iter = 0;
//    int Np =  Result.nx();
//    Bool UseLastScale = (KillLastScale == True) ? False: True;
//  
//    fltarray Resi (Np);
//    fltarray Sol0 (Np);
//    fltarray temp (Np);
//    float Sigma, Noise_Ima;
//    float Delta;
//    float num,den,w;
// 
//    Delta = Sigma = 1.e9;
//    Noise_Ima=NoiseModel->SigmaNoise;
//    if (Noise_Ima < FLOAT_EPSILON) Noise_Ima = 1.;
// 
//    NoiseModel->threshold(MR_Data);
//    MR_Data.rec_adjoint(Sol0, UseLastScale, Border);
// 
//    Result = Sol0;
//    Resi = Result;
//    Sigma = Resi.sigma();
//    for (i=0; i< Np; i++)
//    {
//       if (PositivIma == True)
//             temp(i) = (Result(i) > 0.) ? Result(i) : 0.;
//        else temp(i) = Result(i); 
//        if (MaxIma == True)
// 	    temp(i) = (Result(i) < DEX_MAX_IMA) ? Result(i) : DEX_MAX_IMA;
//        else temp(i) = Result(i);      
//    }
// 
//    while (True)
//    {
//         MR_Data.transform (temp, Border);
//         NoiseModel->threshold(MR_Data);
//         MR_Data.rec_adjoint (temp);
//         for (i=0; i< Np; i++) Resi(i) = Sol0(i) - temp(i);
// 
//         MR_Data.transform (Resi, Border);
//         NoiseModel->threshold(MR_Data);
//         MR_Data.rec_adjoint (temp);
//         num = 0.;
//         den = 0.;
//         for (i=0; i< Np; i++)
//         {
//             num += Resi(i)*temp(i);
//             den += temp(i)*temp(i);
//         }
//         if(!den) w = 1.0;
//         else w = MAX (1.0, num/den); 
// 
//         for (i=0; i< Np; i++)
//         {
//            Result(i) +=  w*Resi(i);
//            if (PositivIma == True)
//               temp(i) = (Result(i) > 0.) ? Result(i) : 0.;
//            else temp(i) = Result(i);
// 	   if (MaxIma == True)
// 	      temp(i) = (Result(i) < DEF_1D_MAX_SIGNAL) ? Result(i) : DEF_1D_MAX_SIGNAL;
// 	   else temp(i) = Result(i);
//         }
// 
//         Sigma = sigma(Resi);
//         Iter ++;
//         if (Verbose == True)
//              cout << "Iter "<< Iter  << ": sigma(resi) = "<< Sigma <<"   convergence coeff  = " <<w<< endl;
//         if (Iter >= Max_Iter)  break;
//    }
//    if (PositivIma == True) threshold(Result);
//    if (MaxIma == True) MaxThreshold(Result);
// }
// 
// /********************************************************************/ 
