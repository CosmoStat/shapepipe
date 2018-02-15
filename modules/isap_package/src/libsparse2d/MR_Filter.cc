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
**    File:  MR_Filtering.cc
**
*******************************************************************************
**
**    DESCRIPTION  Image filtering using multiresolution
**    -----------  
**                 
**
** 
******************************************************************************/

#include "IM_Obj.h"
#include "MR_Obj.h"
#include "IM_IO.h"
#include "MR_NoiseModel.h"
#include "MR_Filter.h"
#include "MR_Sigma.h"
#include "NR.h"
#include "IM_Regul.h"
#include "MR_SoftRegul.h"
#include "MR_Psupport.h"

#define PRINT_DATA 0
#define WRITE_DATA 0

void MaxThreshold (Ifloat &Image);
void bessel_estimate_bayes(float *, int , float, Bool);
void im_wave_variance(Ifloat &Data, Ifloat &Ima_Variance, int SizeBlock, int Step=1);
  
  /********************************************************************/ 
/*
** Operation
**      . A = | Plan [i] | / Noise
**      . If A > 1 then Seuil = 0
**      . Else 
**      .      B = | Plan_Next[i] | / Sigma_Next 
**      .      If B > 1 Then Coef_B = 0     ** linear
**      .      Else B = 1. - 1. / N_Sigma * B     ** interpolation
**      .      Seuil = Noise * B 
**
** Notes:
**     When the wavelet coefficient is superior to the noise 
**   then we don't threshold (Seuil = 0).
**   In the contrary, we look the next plane and we compute
**     B = | Plan_Next[i] | / Sigma_Next. t
**   If B is significant (B > 1), then we don't threshold  (Seuil = 0).
**   If B is under the statistical significant level 
**   then we wheigh the thresholding level by the distance of B
**   from the significant level
*/
static float hierarchical_thresh(float CoefWave, float CoefNext, 
                                 float Sigma, float  SigmaNext)
{
    float Coef_Alpha, Coef_Beta, Seuil;
     	    
    if ((Sigma < FLOAT_EPSILON) || (SigmaNext <  FLOAT_EPSILON))
    {
       Seuil = 0;
    }
    else
    {
       Coef_Alpha = ABS(CoefWave) / Sigma;
       if (Coef_Alpha > 1) Seuil = 0.;
       else
       {
         Coef_Beta = ABS(CoefNext) / SigmaNext;
         if (Coef_Beta > 1) Coef_Beta = 0.;
         else Coef_Beta = 1. - Coef_Beta;
         Seuil =  Sigma * Coef_Beta;
       }
       Seuil = (ABS(CoefWave) < Seuil) ? 0 : CoefWave;
    }
    return Seuil;
}

/********************************************************************/ 
/*
** Plan[i] = x1 * Plan[i]
** with:
**      x1 = S^2 / (B^2 + S^2)
** and
**      B^2 = noise variance= Sigma_Noise * Sigma_Noise
**      S^2 = signal variance = Standart_Deviation(Plan)^2 - B^2
*/
static float wiener(float CoefWave, float SigmaNoise, float SigmaSignal)
{
   float Coef_B2, Coef_S2, Coef_Alpha;

    Coef_B2 = pow ((double)(SigmaNoise),2.);
    Coef_S2 = pow ((double)(SigmaSignal),2.) - Coef_B2;
    if (Coef_S2 < 0) Coef_S2 = 0.; 
    Coef_Alpha = Coef_S2 / (Coef_B2 + Coef_S2);
    return (CoefWave*Coef_Alpha);
}

/********************************************************************/ 
/*
** Hierarchical Wiener Filtering:
**   Plan[i] = x1 * Plan[i] + x2 * Plan_Next[i]
**   with
**      x1 = T^2 / (B^2 + T^2 + Q^2)
**      x2 = B^2 / (B^2 + T^2 + Q^2)
**      Q^2 = T^2 * B^2 / S^2
**    and
**      B^2 = noise variance = Sigma_Noise * Sigma_Noise
**      S^2 = signal variance = Standart_Deviation(Plan)^2 - B^2
**      T^2 = inter-scale dependance variance
**          = Standart_Deviation(Plan - Plan_Next)^2
**
**  These parameters are computed from the hypothesis that
**  signal and noise are gaussian.
*/

static float hierar_wiener(float CoefWave, float CoefNext, float SigmaSignal,
                           float SigmaDiff, float SigmaNoise)
{
    double Coef_T2,Coef_B2,Coef_S2,Coef_Q2,Coef_x1,Coef_x1_Bis;
    float ValRet=0.;
    
    /* Computes the filtering coefficients */
    Coef_T2 = pow ((double)(SigmaDiff),2.);
    Coef_B2 = pow ((double)(SigmaNoise),2.);
    Coef_S2 = pow ((double)(SigmaSignal),2.) - Coef_B2; 

    if (Coef_S2 > 0)
    {
       Coef_Q2 = Coef_T2 * Coef_B2 / Coef_S2;
       Coef_x1 = Coef_T2 / (Coef_B2 + Coef_T2 + Coef_Q2);
       Coef_x1_Bis = Coef_B2 / (Coef_B2 + Coef_T2 + Coef_Q2);
       ValRet = (float) (Coef_x1 * CoefWave  + Coef_x1_Bis * CoefNext);
    }
    return (ValRet);
}                                 


/********************************************************************/ 

void MRFiltering::kill_last_scale(MultiResol & MR_Data)
{
   int s = MR_Data.nbr_band()-1;
   for(int w=0; w < MR_Data.nbr_coeff_in_band(s); w++) MR_Data(s,w) = 0;
}

/********************************************************************/

static void apply_flat_field(Ifloat &Imag, Ifloat Flat)
{
   for (int i=0; i < Imag.nl(); i++)
   for (int j=0; j < Imag.nc(); j++)  
     if (Flat(i,j) > 0) Imag(i,j) /= Flat(i,j);
}
 /********************************************************************/ 

void MRFiltering::filter(Ifloat &Imag, Ifloat &Result,Bool ModifiedATWT)
{
    int b,i,j;
    if (NoiseModel == NULL)
    {
       cerr << "Error: Noise model not initialized ... " << endl;
       exit(-1);
    }
    
   // cout << "INPUT FILTER" << endl;
   // cout << "INPUT FILTER 2" << endl;
  
   int Nl = Imag.nl();
   int Nc = Imag.nc(); 
   int NbrPlan =  NoiseModel->nbr_scale(); 
   MultiResol MR_Data;
   MR_Data.alloc(Nl,Nc,NbrPlan,NoiseModel->type_trans(), 
                 NoiseModel->filter_bank(), NoiseModel->TypeNorm, NoiseModel->nbr_undec_scale(), NoiseModel->U_Filter);

   MR_Data.ModifiedATWT=ModifiedATWT;
   if ((NoiseModel->type_trans() == TO_DIADIC_MALLAT) && (Filter == FILTER_THRESHOLD))
       NoiseModel->GradientAnalysis = True;
         
   if (IsDataModelled == False)  
   {
      NoiseModel->model(Imag, MR_Data);
      IsDataModelled = True;
   }
 
   if ((Verbose == True)  && (StatNoise  == NOISE_GAUSSIAN))
                 cout << "SigmaNoise= " <<  NoiseModel->SigmaNoise << endl;

   // if  TransImag = True ==> wavelet coefficient are not those of the 
   // input image (but the transformed image). We have to transform again the image.       
   if (NoiseModel->TransImag==True)  NoiseModel->im_transform(Imag);
                                                 
   // in this case, WT of image has not yet been computed                                              
   if (IsDataModelled == True)  MR_Data.transform(Imag);
   
   if (BackgroundImage == True)
   {
      if (NoiseModel->TransImag==True) 
         NoiseModel->im_transform(BackgroundData);
      Imag -= BackgroundData;
      MR_Data.transform(Imag);
      if (StatNoise  !=  NOISE_EVENT_POISSON)  NoiseModel->set_support(MR_Data);
      else 
      {
          mr_psupport(MR_Data, *NoiseModel, MR_Data.Border);
          if (NoiseModel->SupIsol == True) for (int b = 0; b < NbrPlan-2; b++) NoiseModel->kill_isol(b);
      }
   }
   
   switch (Filter)
   {
        case FILTER_ITER_THRESHOLD:
             mr_support_filter(Imag, MR_Data, Result);
             break;
        case FILTER_WAVELET_CONSTRAINT:
             mr_wc_filter(Imag, MR_Data, Result);
             break;
	case FILTER_TV_CONSTRAINT:
             mr_support_filter(Imag, MR_Data, Result);
             break;
	case FILTER_ITER_ADJOINT_COEFF:
               grad_adj_filter(MR_Data, Result);
 	       if (UseFlatField == True) apply_flat_field(Result, FlatField);
            break;
        case FILTER_THRESHOLD:
              if (Verbose == True)
              {
                 for (int b=0; b < MR_Data.nbr_band()-1; b++)
                 {
                    float NSig = (NoiseModel->NSigma)[b];
                    cout << "Band " << b+1 << " Nsig = " << NSig << " Sigma = " <<  NoiseModel->sigma(b,0,0)  << " T  = " <<  NoiseModel->sigma(b,0,0) * NSig  <<  endl;
                 }
              }
              NoiseModel->threshold(MR_Data); 
              if (KillLastScale == True) kill_last_scale(MR_Data);
              MR_Data.recons(Result);
	      if (UseFlatField == True) apply_flat_field(Result, FlatField);
               // if (MR_Data.Type_Transform != TO_HAAR) MR_Data.recons(Result);
              // else mr_ortho_regul_rec (MR_Data, Result, 10);
              // mr_ortho_regul_ima_rec(MR_Data, Result,10);
 
              break;
        case FILTER_HIERARCHICAL_TRESHOLD:
        case FILTER_HIERARCHICAL_WIENER:
	case FILTER_BIVARIATE_SHRINKAGE:
              hierar_filter (MR_Data, Result);
              if (KillLastScale == True) kill_last_scale(MR_Data);
              MR_Data.recons(Result);
 	      if (UseFlatField == True) apply_flat_field(Result, FlatField);
              break;
        case FILTER_MULTI_RES_WIENER:
	    {
	      int BlockSize =  NoiseModel->SizeBlockSigmaNoise;
              int Step = 1;
 	      for (b=0; b < MR_Data.nbr_band()-1; b++)
              {
 	        int B2 = BlockSize / 2;
		if (b < MR_Data.nbr_band_per_resol()) 
		{
		  B2 = NoiseModel->SizeBlockSigmaNoise/2+1;
		  BlockSize = 2*B2+1;
		}
		else 
		{
		  BlockSize =  NoiseModel->SizeBlockSigmaNoise;
		  B2 = BlockSize / 2;
		}
 	        float NSig = 1. + sqrt(2./(BlockSize*BlockSize));
		double Sigma;
 //  		for (i=0;i < MR_Data.size_band_nl(b); i++) 
// 		for (j=0;j < MR_Data.size_band_nc(b); j++)   
// 		{
// 		   Sigma = 0.;
// 		   for (int k=i-B2*Step; k <= i+B2*Step; k+=Step) 
// 		   for (int l=j-B2*Step; l <= j+B2*Step; l+=Step)  
// 		      Sigma += (MR_Data.band(b))(k,l,I_MIRROR)*(MR_Data.band(b))(k,l,I_MIRROR);
// 		   Sigma /= (float)(BlockSize*BlockSize); 
//                    float NoiseLevel = NoiseModel->sigma(b,i,j)*NoiseModel->sigma(b,i,j);
// 		   if (Sigma < NSig*NoiseLevel) MR_Data(b,i,j) =0.;
//   		}
 		for (i=0;i < MR_Data.size_band_nl(b); i++) 
		for (j=0;j < MR_Data.size_band_nc(b); j++)   
		{
		   if (MR_Data(b,i,j) != 0)
		   {
		      Sigma = 0.;
		      for (int k=i-B2*Step; k <= i+B2*Step; k+=Step) 
		      for (int l=j-B2*Step; l <= j+B2*Step; l+=Step)  
		         Sigma += (MR_Data.band(b))(k,l,I_MIRROR)*(MR_Data.band(b))(k,l,I_MIRROR);
		      Sigma /= (float)(BlockSize*BlockSize); 
                     float Coef = MR_Data(b,i,j);
                     float NoiseLevel = NoiseModel->sigma(b,i,j)*NoiseModel->sigma(b,i,j);
                     float SigmaSignal = MAX(0, Sigma - NoiseLevel);
                     if (Sigma > 0) MR_Data(b,i,j) *= SigmaSignal / Sigma;
                     else MR_Data(b,i,j) = 0.;
		   }
 		}
		
                // if ((MR_Data.nbr_undec_scale() > b/MR_Data.nbr_band_per_resol()) && ((b+1) % MR_Data.nbr_band_per_resol() == 0)) Step *= 2;
  	    }
	    if (KillLastScale == True) kill_last_scale(MR_Data);
            MR_Data.recons(Result);
 	      if (UseFlatField == True) apply_flat_field(Result, FlatField);
	   }
	    break;
        case FILTER_SOFT_THRESHOLD:
	       for (b=0; b < MR_Data.nbr_band()-1; b++)
	       {
                  float NSig = (NoiseModel->NSigma)[b];
		   if (Verbose == True)
                     cout << "Band " << b+1 << " Nsig = " << (NoiseModel->NSigma)[b] << " Sigma = " <<  NoiseModel->sigma(b,0,0)  << " T  = " <<  NoiseModel->sigma(b,0,0) * (NoiseModel->NSigma)[b]  <<  endl;
  		  for (i=0;i < MR_Data.size_band_nl(b); i++) 
		  for (j=0;j < MR_Data.size_band_nc(b); j++)
 		       MR_Data(b,i,j)= soft_threshold(MR_Data(b,i,j), NSig*NoiseModel->sigma(b,i,j));
 	       }
	       if (KillLastScale == True) kill_last_scale(MR_Data);
                MR_Data.recons(Result);
 	       if (UseFlatField == True) apply_flat_field(Result, FlatField);
	       break;
        case FILTER_MULTI_HARD_MAD:
	       for (b=0; b < MR_Data.nbr_band()-1; b++)
	       {
	          float NoiseMad = detect_noise_from_mad(MR_Data.band(b));
                  float NoiseLevel = (NoiseModel->NSigma)[b] * NoiseMad;
		  if (Verbose == True)
		   cout << "Band " << b+1 << " Nsig = " << (NoiseModel->NSigma)[b] << " MAD sigma = " << NoiseMad << " T = " <<  NoiseLevel <<  endl;

		  for (i=0;i < MR_Data.size_band_nl(b); i++) 
		  for (j=0;j < MR_Data.size_band_nc(b); j++)
 		      if (ABS(MR_Data(b,i,j)) < NoiseLevel) MR_Data(b,i,j) = 0.;
 	       }
	       if (KillLastScale == True) kill_last_scale(MR_Data);
               MR_Data.recons(Result);
 	       if (UseFlatField == True) apply_flat_field(Result, FlatField);
	       break;
        case FILTER_MULTI_SOFT_MAD:
	       for (b=0; b < MR_Data.nbr_band()-1; b++)
	       {
	          float NoiseMad = detect_noise_from_mad(MR_Data.band(b));
 		  float NoiseLevel = (NoiseModel->NSigma)[b] * NoiseMad;
		  if (Verbose == True)
		   cout << "Band " << b+1 << " Nsig = " << (NoiseModel->NSigma)[b] << " MAD sigma = " << NoiseMad << " T = " <<  NoiseLevel <<  endl;

		  for (i=0;i < MR_Data.size_band_nl(b); i++) 
		  for (j=0;j < MR_Data.size_band_nc(b); j++)
 		       MR_Data(b,i,j) = soft_threshold(MR_Data(b,i,j),NoiseLevel);
  	       }
	       if (KillLastScale == True) kill_last_scale(MR_Data);
               MR_Data.recons(Result);
 	       if (UseFlatField == True) apply_flat_field(Result, FlatField);
 	       break;
 	case FILTER_BAYES_BESSEL:
	        for (b=0; b < MR_Data.nbr_band()-1; b++)
	        {
	          float *Ptr = MR_Data.band(b).buffer();
	          int N = MR_Data.band(b).nl()*MR_Data.band(b).nc();
	          bessel_estimate_bayes(Ptr, N,  NoiseModel->sigma(b,0,0), True);
                }
		if (KillLastScale == True) kill_last_scale(MR_Data);
                MR_Data.recons(Result);
 	        if (UseFlatField == True) apply_flat_field(Result, FlatField);
                break;
         default:           
           cerr << "Error: unknown filtering method  ... " << endl;
          exit(-1);
     }
   
   // if the data have been transform, we have to take the inverse transform
   // of the solution
   if (NoiseModel->TransImag == True)
   { 
      if (KillLastScale != True)
      {
         NoiseModel->im_invtransform(Imag);
         NoiseModel->im_invtransform(Result);
      }
      else
      {
          Result = Imag - Result;
          NoiseModel->im_invtransform(Imag);
          NoiseModel->im_invtransform(Result);
          Result = Imag - Result;
      }
   }       
   
   if (BackgroundImage == True) Imag += BackgroundData;
   if (PositivIma == True) threshold(Result);
   if (MaxIma == True) MaxThreshold(Result);
} 
 
/********************************************************************/ 
/****************************************************************************/

void im_wave_variance(Ifloat &Data, Ifloat &Ima_Variance, int SizeBlock, int Step)
{
  int i,j,k,l;
  int Dsize = SizeBlock / 2 * Step;
  float  S0 = SizeBlock*SizeBlock;
 
  for (i=0; i< Data.nl(); i++)
  for (j=0; j< Data.nc(); j++)
  {
     Ima_Variance(i,j)=0;
     for (k=i-Dsize; k<= i+Dsize; k+=Step)
     for (l=j-Dsize; l<= j+Dsize; l+=Step)
     {
         Ima_Variance(i,j) +=Data(k,l,I_MIRROR)*Data(k,l,I_MIRROR);
      }
     Ima_Variance(i,j) /= S0;
  }
}

/****************************************************************************/
void MRFiltering::hierar_filter(MultiResol & MR_Data, Ifloat &Result)
{
    int Nl = Result.nl ();
    int Nc = Result.nc ();
    int i,j,s;
    float Coef,w,Nsig;
    float NoiseLevel, SigmaNext=0., SigmaSignal=0., SigmaDiff=0., SigmaData=0.;
    Ifloat Plan_Next(Nl, Nc, (char *) "Plan_Next");
    Ifloat Ima_DifVariance;
    int NbPerRes = MR_Data.nbr_band_per_resol();
    Ifloat Ima_Variance;
//     if (isotrop(MR_Data.Type_Transform) != True)
//     {
//        cerr << "Error: this filtering method need an isotropic transform ... " << endl;
//        exit(-1);
//     }
    
    /* Filtering scale by scale */
    int ind=1;
    int Step=1;
    intarray TabStep(MR_Data.nbr_band());
    for (s = 0; s < MR_Data.nbr_band() -1; s++)
    {
       TabStep(s) = Step;
       if ((MR_Data.nbr_undec_scale() > s/NbPerRes) && (ind % NbPerRes == 0)) Step *= 2;
       ind = ind + 1;
    }
    int BlockSize =  NoiseModel->SizeBlockSigmaNoise;
    for (s =  MR_Data.nbr_band() - 2 - NbPerRes; s >= 0; s--)
    {
        Nl = MR_Data.band(s).nl ();
        Nc = MR_Data.band(s).nc ();
	
	Ima_Variance.resize(Nl,Nc);
        Plan_Next.resize (Nl, Nc);
        im_block_extend(MR_Data.band(s+NbPerRes), Plan_Next);
 	
  	Nsig = NoiseModel->NSigma[s];
        if (Filter == FILTER_HIERARCHICAL_WIENER)
        {
	    Ima_Variance = MR_Data.band(s);
	    Ima_Variance  -= Plan_Next;
	    Ima_DifVariance.resize (Nl, Nc);
            im_wave_variance(Ima_Variance, Ima_DifVariance, BlockSize, TabStep(s));
        }
	else if (Filter == FILTER_HIERARCHICAL_TRESHOLD)
	{
	   Ima_DifVariance.resize (Nl, Nc);
	   im_wave_variance(Plan_Next, Ima_DifVariance, BlockSize, TabStep(s)*2);
	}
	im_wave_variance(MR_Data.band(s), Ima_Variance, BlockSize, TabStep(s));
        
        for (i=0;i<Nl;i++)
        for (j=0;j<Nc;j++)
        {
           float C1,C2;
	   SigmaData  = sqrt(Ima_Variance(i,j));
           NoiseLevel = NoiseModel->sigma(s,i,j);
           SigmaSignal = MAX(0, SigmaData - NoiseLevel);
           Coef = MR_Data(s,i,j);
           switch (Filter)
           {
             case FILTER_HIERARCHICAL_TRESHOLD:
                 w = linear_weight(Coef, NoiseLevel, Nsig);
		 SigmaNext = Ima_DifVariance(i,j); 
                 MR_Data(s,i,j) = w*Coef+(1-w)* hierarchical_thresh(Coef, 
                                          Plan_Next(i,j),
					  NoiseLevel*Nsig,
					  SigmaNext*NoiseModel->NSigma[s+NbPerRes]);
                 break;
            case FILTER_HIERARCHICAL_WIENER:
	         SigmaDiff = sqrt(Ima_DifVariance(i,j));
                 // w = linear_weight(Coef, NoiseLevel, Nsig);
                 // MR_Data(s,i,j) =  w*Coef+(1-w)*hierar_wiener (Coef,  
                 //                      Plan_Next(i,j), SigmaSignal, SigmaDiff,
                 //                      NoiseLevel);
		 MR_Data(s,i,j) = hierar_wiener(MR_Data(s,i,j), Plan_Next(i,j), SigmaData, SigmaDiff, NoiseLevel); 
		 break;
           case FILTER_BIVARIATE_SHRINKAGE:
 		 SigmaSignal = Ima_Variance(i,j) - NoiseLevel*NoiseLevel;
		 if (SigmaSignal < 0) SigmaSignal = 0.;
		 else SigmaSignal = sqrt(SigmaSignal);
		 C2 = sqrt(Coef*Coef+Plan_Next(i,j)*Plan_Next(i,j));
		 C1 = C2 - sqrt(3.)*NoiseLevel*NoiseLevel/SigmaSignal;
		 if (C1 < 0) C1 = 0.;
		 if (C2 == 0) MR_Data(s,i,j) = 0.;
		 else  MR_Data(s,i,j) *= C1 / C2;
                 break;
             default:
                  fprintf (stderr, "Bad Type Filtering\n");
                  exit (-1);
           }
       }
    }
}

/********************************************************************/ 

// void MRFiltering::mr_wc_filter (Ifloat &Imag, MultiResol & MR_Data, Ifloat &Result)
// {
//    int b,i,j,Iter = 0;
//    int Nl = Imag.nl();
//    int Nc = Imag.nc();
//    float Resi;
//    float Sigma, Noise_Ima;
//    float Old_Sigma, Delta;
//    int NbrPlan =  NoiseModel->nbr_scale(); 
//    MR_Regul RI;
//    RI.ExpDecreasingLambda = True;
//    RI.NbrScale = MR_Data.nbr_scale();
//    MultiResol MR_Sol;
//    
//    MR_Sol.alloc(Nl,Nc,NbrPlan,NoiseModel->type_trans(), 
//                  NoiseModel->filter_bank(), NoiseModel->TypeNorm, NoiseModel->nbr_undec_scale());
// 
//    MR_Data.transform (Imag, Border); 
//    NoiseModel->threshold(MR_Data);
//    MR_Data.recons (Result);
//    
//    Delta = Sigma = 1.e9;
//    Noise_Ima=NoiseModel->SigmaNoise;
//    if (Noise_Ima < FLOAT_EPSILON) Noise_Ima = 1.;
//   
//    Sigma = energy(Imag);
//    Old_Sigma = Sigma; 
//    float SoftThreshold, Lambda = 1.;
//    float StepL = 1. / (float) Max_Iter;
//    float DetectCoefTol = 0.5;
//    while (True)
//    {
//         //cout << "----------------------------------------------------------" << endl;
// 	//cout << "Iter number " << Iter << endl;
//         //cout << "Begin loop Sigma="<<Sigma<<", OldSigma="<<Old_Sigma<< endl;
//         Lambda -= StepL; 
//         if (Lambda < 0) Lambda = 0.;
//           
//         MR_Sol.transform (Result, Border); 
// 	for (b=0; b < MR_Data.nbr_band()-1; b++)
// 	for (i=0; i < MR_Data.size_band_nl(b); i++)
// 	for (j=0; j < MR_Data.size_band_nc(b); j++)
// 	{
//            //float NSig = (NoiseModel->NSigma)[b];
//  	   float Noise =  (*NoiseModel).sigma(b,i,j);	   
// 	   float Interval = Noise * DetectCoefTol;
//            Resi = MR_Data(b,i,j) - MR_Sol(b,i,j);
//            SoftThreshold = Lambda*Noise;
// 	   if ((*NoiseModel)(b,i,j) == True)
// 	   {
//  	      if (ABS(Resi) > Interval)
// 	      {
// 	         MR_Sol(b,i,j) =  MR_Data(b,i,j);
// 	      }
// 	   }
// 	   else 
// 	   {
// 	      if (MR_Sol(b,i,j) > Interval) MR_Sol(b,i,j) =  Interval;
// 	      else if (MR_Sol(b,i,j) < -Interval) MR_Sol(b,i,j) = -Interval;
//  	   }
// 	   MR_Sol(b,i,j) = soft_threshold(MR_Sol(b,i,j), SoftThreshold);
// 	}
//         b = MR_Data.nbr_band()-1;
// 	MR_Sol.band(b) = MR_Data.band(b);
// 	
//         MR_Sol.recons (Result);
//         if (RegulParam > 0)      
// 	{
// 	   float LambdaTV = RegulParam*Noise_Ima;
//            RI.im_soft_threshold(Result, Result, LambdaTV);
//         }
// 	if ((PositivIma == True) && (NoiseModel->PoissonFisz == False))
// 	                             threshold(Result);
//         if (MaxIma == True) MaxThreshold(Result);
// 
//         Old_Sigma = Sigma; 
//         Sigma = energy(Result);
//         Delta = (Old_Sigma - Sigma)  / (Nl*Nc*Noise_Ima*Noise_Ima);
// 
//         Iter ++;
//         float Xi = Sigma / (Nl*Nc*Noise_Ima*Noise_Ima);
//         if (Verbose == True)
//         {
//            if (((NoiseModel->which_noise() == NOISE_GAUSSIAN) ||
//                (NoiseModel->which_noise() == NOISE_POISSON)) &&
//                (KillLastScale == False))
//              cout << "Iter "<< Iter+1 << ": Xi = "<< Xi <<" ==> Delta  = " <<Delta<< endl;
//            else cout << "Iter "<< Iter+1  <<" ==> Delta  = " <<Delta<< endl;
//         }
//         if ((Iter >= Max_Iter) || (ABS(Delta) < Epsilon)) break;
// 	
//         //cout << "Iter numer " << Iter-1 << endl;
//         //cout << "End   Loop Sigma="<<Sigma<<", OldSigma="<<Old_Sigma<< endl;
//         //cout << "Nrj (imag to filter [IN])  = " << energy(Resi) << endl;
//         //cout << "Nrj (filtered imag  [OUT]) = " << energy(Result) << endl;   	
//     }
// }

void MRFiltering::mr_wc_filter (Ifloat &Imag, MultiResol & MR_Data, Ifloat &Result)
{
   int b,i,j,Iter = 0;
   int Nl = Imag.nl();
   int Nc = Imag.nc();
   Ifloat Resi (Nl, Nc, (char *) "Residual");
   float Sigma, Noise_Ima;
   float Old_Sigma, Delta, LambdaTV;
   MR_Regul RI;
   RI.ExpDecreasingLambda = True;
   

   Delta = Sigma = 1.e9;
   Noise_Ima=NoiseModel->SigmaNoise;
   if (Noise_Ima < FLOAT_EPSILON) Noise_Ima = 1.;
   Result.init();
   Resi = Imag;
   if (UseFlatField == True)
   {
       if ((FlatField.nl() != Imag.nl()) || (FlatField.nc() != Imag.nc()))
       {
           cout << "Error: flat field image has not the same size as the input image ... " << endl;
	   exit(-1);
       }
    }
       
   LambdaTV = RegulParam*Noise_Ima;
   if (Filter == FILTER_TV_CONSTRAINT)
   {
      RI.NbrScale = 2;
      if (RegulParam <= 0) LambdaTV = 0.1*Noise_Ima;
   }
   else RI.NbrScale = MR_Data.nbr_scale();
    
   Sigma = energy(Resi);
   Old_Sigma = Sigma; 
   Sup_Set=False;
   MR_Data.ModifiedATWT = True;
   UseAdjointRec = False;
   if (Verbose == True)
   {
      cout << "  Nbr Max iter = " << Max_Iter<< endl;
      cout << "  Epsilon cvg  = " << Epsilon<< endl;
      if (RegulParam > 0) cout << "  TV regul param = " << RegulParam<< endl;
      if (UseAdjointRec == True) cout << "  Use Adjoint reconstruction " << endl;
      if (Sup_Set == True) cout << "  Update the support at each iter" << endl;
      if (PositivIma == True)  cout << "  Positivity " << endl;
      if (NoiseModel->OnlyPositivDetect == True)  cout << "  Detect only positive coeff" << endl;
   }
   
   float SoftThreshold, Lambda = 1.;
   float StepL = RegulParam / (float) (Max_Iter);
   float DetectCoefTol = 1.;
   float MaxFlat = (UseFlatField == False) ? 1: FlatField.max();
   while (True)
   {
        Lambda -= StepL; 
        if (Lambda < 0) Lambda = 0.;
        //cout << "----------------------------------------------------------" << endl;
	//cout << "Iter numer " << Iter << endl;
        //cout << "Begin loop Sigma="<<Sigma<<", OldSigma="<<Old_Sigma<< endl;
        //cout << "Nrj (imag to filter [IN])  = " << energy(Resi) << endl;
        //cout << "Nrj (filtered imag  [OUT]) = " << energy(Result) << endl;   
    
        MR_Data.transform (Resi, Border); 
        // NoiseModel->threshold(MR_Data, Sup_Set);
 	for (b=0; b < MR_Data.nbr_band()-1; b++)
	{
 	   for (i=0; i < MR_Data.size_band_nl(b); i++)
	   for (j=0; j < MR_Data.size_band_nc(b); j++)
	   {
	       float Noise = (*NoiseModel).sigma(b,i,j);
	       float NSig = (*NoiseModel).nsigma(b);	   
 	       float Interval = Noise*DetectCoefTol;
	       // if (ABS(MR_Data(b,i,j)) < Interval) MR_Data(b,i,j) = 0;
	       if (((*NoiseModel)(b,i,j) == False) || (ABS(MR_Data(b,i,j)) < Interval)) MR_Data(b,i,j) = 0;
  	   }
	}          
 	   
        if (KillLastScale == True) kill_last_scale(MR_Data);
        if (UseAdjointRec == True) MR_Data.rec_adjoint (Resi);
        else MR_Data.recons (Resi);
	// INFO(Resi,"resi");
        if (UseFlatField == False) Result += Resi; 
	else 
	{
	   for (i=0; i < Resi.nl(); i++)
	   for (j=0; j < Resi.nc(); j++)  
	   {
	      float Val =  (ABS(FlatField(i,j)) > 0) ? Resi(i,j) / FlatField(i,j) : 0.;
	      Result(i,j) += Val;
	   }
        }
	
    if (Lambda*RegulParam> 0)
	{
	   MR_Data.transform (Result, Border); 
  	   for (b=0; b < MR_Data.nbr_band()-1; b++)
	   {
 	      for (i=0; i < MR_Data.size_band_nl(b); i++)
	      for (j=0; j < MR_Data.size_band_nc(b); j++)
	      {
	         float Noise = (*NoiseModel).sigma(b,i,j) * MaxFlat;	  
 	         float Interval = Noise * (*NoiseModel).nsigma(b);
                 SoftThreshold = Lambda*Noise*RegulParam;
    	         MR_Data(b,i,j) = soft_threshold(MR_Data(b,i,j), SoftThreshold);
	       // if (((*NoiseModel)(b,i,j) == False) && (MR_Data(b,i,j) > Interval)) MR_Data(b,i,j) = Interval;
	       // else if (((*NoiseModel)(b,i,j) == False) && (MR_Data(b,i,j) < -Interval)) MR_Data(b,i,j) = -Interval;
 	      }
	   }          
	   MR_Data.recons(Result);
	}
	if ((PositivIma == True) && (NoiseModel->PoissonFisz == False))
	                             threshold(Result);
        if (MaxIma == True) MaxThreshold(Result);
        Old_Sigma = Sigma; 
        if (UseFlatField == False) 
	{
	   for (i=0; i < Resi.nl(); i++)
	   for (j=0; j < Resi.nc(); j++)  Resi(i,j) = Imag(i,j)  - Result(i,j);
	}
	else
	{
	   for (i=0; i < Resi.nl(); i++)
	   for (j=0; j < Resi.nc(); j++) Resi(i,j) = Imag(i,j)  - Result(i,j)*FlatField(i,j);
	}
    if (MissingData == True)
    {
        for (i=0; i < Resi.nl(); i++)
        for (j=0; j < Resi.nc(); j++) Resi(i,j) *=  MaskIma(i,j);
    }
       
	
        Sigma = energy(Resi);
        Delta = (Old_Sigma - Sigma)  / (Nl*Nc*Noise_Ima*Noise_Ima);

        Iter ++;
        float Xi = Sigma / (Nl*Nc*Noise_Ima*Noise_Ima);
        if (Verbose == True)
        {
           if (((NoiseModel->which_noise() == NOISE_GAUSSIAN) ||
               (NoiseModel->which_noise() == NOISE_POISSON)) &&
               (KillLastScale == False))
             cout << "Iter "<< Iter+1 << ": Xi = "<< Xi <<" ==> Delta  = " <<Delta<< " Lambda = " << Lambda << endl;
           else cout << "Iter "<< Iter+1  <<" ==> Delta  = " << Delta <<  ", Lambda = " << Lambda << endl;
        }
        // if ((Iter >= Max_Iter) || (Delta < Epsilon)) break;
	 if (Iter >= Max_Iter) break;
	
        //cout << "Iter numer " << Iter-1 << endl;
        //cout << "End   Loop Sigma="<<Sigma<<", OldSigma="<<Old_Sigma<< endl;
        //cout << "Nrj (imag to filter [IN])  = " << energy(Resi) << endl;
        //cout << "Nrj (filtered imag  [OUT]) = " << energy(Result) << endl;   	
    }
}

/********************************************************************/

void MRFiltering::mr_support_filter(Ifloat &Imag, MultiResol & MR_Data, Ifloat &Result)
{
   int b,i,j,Iter = 0;
   int Nl = Imag.nl();
   int Nc = Imag.nc();
   Ifloat Resi (Nl, Nc, (char *) "Residual");
   float Sigma, Noise_Ima;
   float Old_Sigma, Delta, LambdaL1, LambdaTV;
   MR_Regul RI;
   RI.ExpDecreasingLambda = True;
   

   Delta = Sigma = 1.e9;
   Noise_Ima=NoiseModel->SigmaNoise;
   if (Noise_Ima < FLOAT_EPSILON) Noise_Ima = 1.;
   Result.init();
   if (UseFlatField == True)
   {
       if ((FlatField.nl() != Imag.nl()) || (FlatField.nc() != Imag.nc()))
       {
           cout << "Error: flat field image has not the same size as the input image ... " << endl;
	   exit(-1);
       }
    }
   Resi = Imag;
   
   LambdaTV = RegulParam*Noise_Ima;
   if (Filter == FILTER_TV_CONSTRAINT)
   {
      RI.NbrScale = 2;
      if (RegulParam <= 0) LambdaTV = 0.1*Noise_Ima;
   }
   else 
   {
      RI.NbrScale = MR_Data.nbr_scale();
      if (RegulParam == 0.1) 
      {
         LambdaL1 = 0.;  // 0.1 is the default value in mr_filter, and by default no l_1 regul
      }
      else LambdaL1 = RegulParam;
   }
   Sigma = energy(Resi);
   Old_Sigma = Sigma; 
   
   if (MissingData == True) 
   {
      Sup_Set = False;
      Border= I_MIRROR;
   }
   if (Sup_Set == False)
   {
      MR_Data.ModifiedATWT = True;
      UseAdjointRec = False;
   }
    
   if (Verbose == True)
   {
      cout << "  Nbr Max iter = " << Max_Iter<< endl;
      cout << "  Epsilon cvg  = " << Epsilon<< endl;
      if (RegulParam > 0) cout << "  Regul param = " << RegulParam<< endl;
      if (UseAdjointRec == True) cout << "  Use Adjoint reconstruction " << endl;
      if (Sup_Set == True) cout << "  Update the support at each iter" << endl;
      if (PositivIma == True)  cout << "  Positivity " << endl;
      if (NoiseModel->OnlyPositivDetect == True)  cout << "  Detect only positive coeff" << endl;
   }

     
  
   while (True)
   {
        //cout << "----------------------------------------------------------" << endl;
	//cout << "Iter numer " << Iter << endl;
        //cout << "Begin loop Sigma="<<Sigma<<", OldSigma="<<Old_Sigma<< endl;
        //cout << "Nrj (imag to filter [IN])  = " << energy(Resi) << endl;
        //cout << "Nrj (filtered imag  [OUT]) = " << energy(Result) << endl;   
    
        MR_Data.transform (Resi, Border); 
        NoiseModel->threshold(MR_Data, Sup_Set);
        
        if (KillLastScale == True) kill_last_scale(MR_Data);
        if ((Iter ==0) && (MissingData == True) && (LambdaL1 == 0))
        {
            
            for (b=0; b < MR_Data.nbr_band()-1; b++)
            for (i=0;i < MR_Data.size_band_nl(b); i++) 
            for (j=0;j < MR_Data.size_band_nc(b); j++)
            {
                if (ABS(MR_Data(b,i,j) / MR_Data.band_norm(b)) > LambdaL1) LambdaL1 = ABS(MR_Data(b,i,j) ) / MR_Data.band_norm(b) ;
                LambdaL1 /= Noise_Ima;
                RegulParam = LambdaL1;
            }
        }
	       
        if (UseAdjointRec == True) MR_Data.rec_adjoint (Resi);
        else MR_Data.recons (Resi);
	// INFO(Resi,"resi");
	 
        if (UseFlatField == False) Result += Resi; 
	else 
	{
	   for (i=0; i < Resi.nl(); i++)
	   for (j=0; j < Resi.nc(); j++)  
	   {
	      float Val =  (ABS(FlatField(i,j)) > 0) ? Resi(i,j) / FlatField(i,j) : 0.;
	      Result(i,j) += Val;
	   }
        }
        
	if ((RegulParam > 0) && (Filter == FILTER_TV_CONSTRAINT))  RI.im_soft_threshold(Result, Result, LambdaTV);
    else if (LambdaL1 > 0)
	{
	   MR_Data.transform (Result, Border); 
       for (b=0; b < MR_Data.nbr_band(); b++)
	   {
              for (i=0;i < MR_Data.size_band_nl(b); i++) 
              for (j=0;j < MR_Data.size_band_nc(b); j++)
              {
                  if (b != MR_Data.nbr_band()-1) NoiseModel->support(b,i,j) = (NoiseModel->signif(MR_Data(b,i,j), b,i,j)) ? VAL_SupOK: VAL_SupNull;
                  MR_Data(b,i,j) = soft_threshold(MR_Data(b,i,j), LambdaL1*MR_Data.band_norm(b));
              }
  	   }
	   LambdaL1 -= RegulParam / (Max_Iter+1.);
	   if (LambdaL1 < 0) LambdaL1  = 0.;
	   MR_Data.recons (Result);
    }       
	
 	if ((PositivIma == True) && (NoiseModel->PoissonFisz == False))  threshold(Result);
        if (MaxIma == True) MaxThreshold(Result);
	
        Old_Sigma = Sigma; 
	if (UseFlatField == False) 
	{
	   for (i=0; i < Resi.nl(); i++)
	   for (j=0; j < Resi.nc(); j++)  Resi(i,j) = Imag(i,j)  - Result(i,j);
	}
	else
	{
	   for (i=0; i < Resi.nl(); i++)
	   for (j=0; j < Resi.nc(); j++) Resi(i,j) = Imag(i,j)  - Result(i,j)*FlatField(i,j);
	}
    if (MissingData == True)
    {
        for (i=0; i < Resi.nl(); i++)
        for (j=0; j < Resi.nc(); j++) Resi(i,j) *=  MaskIma(i,j);
    }
    
     	
        Sigma = energy(Resi);
        Delta = (Old_Sigma - Sigma)  / (Nl*Nc*Noise_Ima*Noise_Ima);

        Iter ++;
        float Xi = Sigma / (Nl*Nc*Noise_Ima*Noise_Ima);
        if (Verbose == True)
        {
           if (((NoiseModel->which_noise() == NOISE_GAUSSIAN) ||
               (NoiseModel->which_noise() == NOISE_POISSON)) &&
               (KillLastScale == False))
             cout << "Iter "<< Iter+1 << ": Xi = "<< Xi <<" ==> Delta  = " << Delta <<  ", Lambda = " << LambdaL1 <<endl;
           else cout << "Iter "<< Iter+1  <<" ==> Delta  = " <<Delta << " Lambda = " << LambdaL1 << endl;
        }
        if ((LambdaL1 == 0) && ((Iter >= Max_Iter) || (ABS(Delta) < Epsilon))) break;
//         if (Iter > 3)
//         {
//            MR_Data.ModifiedATWT = True;
//            Sup_Set = False;
//         }
        //cout << "Iter numer " << Iter-1 << endl;
        //cout << "End   Loop Sigma="<<Sigma<<", OldSigma="<<Old_Sigma<< endl;
        //cout << "Nrj (imag to filter [IN])  = " << energy(Resi) << endl;
        //cout << "Nrj (filtered imag  [OUT]) = " << energy(Result) << endl;   	
    }
}

/********************************************************************/ 

void MRFiltering::grad_adj_filter(MultiResol &MR_Data, Ifloat &Result)
{
   int i,j,Iter = 0;
   int Nl =  Result.nl();
   int Nc =  Result.nc();
   Bool UseLastScale = (KillLastScale == True) ? False: True;
 
   Ifloat Resi (Nl, Nc, (char *) "Residual");
   Ifloat Sol0 (Nl, Nc, (char *) "Sol0");
   Ifloat temp (Nl, Nc, (char *) "temp");
   float Sigma, Noise_Ima;
   float Delta;
   float num,den,w;

   Delta = Sigma = 1.e9;
   Noise_Ima=NoiseModel->SigmaNoise;
   if (Noise_Ima < FLOAT_EPSILON) Noise_Ima = 1.;

   NoiseModel->threshold(MR_Data);
   MR_Data.rec_adjoint(Sol0, UseLastScale, Border);

   Result = Sol0;
   Resi = Result;
   Sigma = sigma(Resi);
   for (i=0; i< Nl; i++)
   for (j=0; j< Nc; j++)
   {
      if (PositivIma == True)
            temp(i,j) = (Result(i,j) > 0.) ? Result(i,j) : 0.;
       else temp(i,j) = Result(i,j); 
       if (MaxIma == True)
	    temp(i,j) = (Result(i,j) < DEX_MAX_IMA) ? Result(i,j) : DEX_MAX_IMA;
       else temp(i,j) = Result(i,j);      
   }

   while (True)
   {
        MR_Data.transform (temp, Border);
        NoiseModel->threshold(MR_Data);
        MR_Data.rec_adjoint (temp);
        Resi = Sol0 - temp;

        MR_Data.transform (Resi, Border);
        NoiseModel->threshold(MR_Data);
        MR_Data.rec_adjoint (temp);
        num = 0.;
        den = 0.;
        for (i=0; i< Nl; i++)
        for (j=0; j< Nc; j++) 
        {
            num += Resi(i,j)*temp(i,j);
            den += temp(i,j)*temp(i,j);
        }
        if(!den) w = 1.0;
        else w = MAX (1.0, num/den); 

        for (i=0; i< Nl; i++)
        for (j=0; j< Nc; j++) 
        {
           Result(i,j) +=  w*Resi(i,j);
           if (PositivIma == True)
              temp(i,j) = (Result(i,j) > 0.) ? Result(i,j) : 0.;
           else temp(i,j) = Result(i,j);
	   if (MaxIma == True)
	      temp(i,j) = (Result(i,j) < DEX_MAX_IMA) ? Result(i,j) : DEX_MAX_IMA;
	   else temp(i,j) = Result(i,j);
        }

        Sigma = sigma(Resi);
        Iter ++;
        if (Verbose == True)
             cout << "Iter "<< Iter  << ": sigma(resi) = "<< Sigma <<"   convergence coeff  = " <<w<< endl;
        if (Iter >= Max_Iter)  break;
   }
   if (PositivIma == True) threshold(Result);
   if (MaxIma == True) MaxThreshold(Result);
}


/********************************************************************/ 
void MaxThreshold (Ifloat &Image) {
    int i,j;
    for (i = 0; i < Image.nl(); i++)
    for (j = 0; j < Image.nc(); j++) if (Image(i,j) > DEX_MAX_IMA) Image(i,j) = DEX_MAX_IMA;
}
