/******************************************************************************
**                   Copyright (C) 1998 by CEA
*******************************************************************************
**
**    UNIT
**
**    Version: 1.0
**
**    Author: Jean-Luc Starck
**
**    Date:  98/01/12 
**    
**    File:  MW_Filter.cc
**
*******************************************************************************
**
**    DESCRIPTION  Image filtering using the multiscale entropy
**    ----------- 
**                 
**
******************************************************************************/

#include "IM_Obj.h"
// #include "IM_IO.h"
#include "MR_Obj.h"
#include "MR_NoiseModel.h"
#include "CErf.h"
#include "CMem.h"
#include "MW_Filter.h"


/****************************************************************************/ 

void mw_filter(MultiResol &MR_Data, MRNoiseModel & NoiseModel, 
        int TypeFilter, fltarray &TabAlpha, Bool DataSNR, 
	Bool DataModel, MultiResol *MR_Model,
	float ConvgParam, int MaxIter, Bool PositivConst, Bool Verbose, Bool UseSupport)
/* Apply the multiscale entropy to the multiresolution data set MR_Data
   MR_Data = in-out: multiresolution transform of an image
   NoiseModel = in: noise modeling of the image in the wavelet space
   TypeFilter = in: type of filtering
             = DEF_MEM_ALPHA_CST  ==> alpha is fixed 
	     = DEF_MEM_ALPHA_OPT  ==> alpha is globally optimized
	     = DEF_MEM_ALPHA_OPT_SCALE ==> alpha is optimized at each scale
	     
   TabAlpha = in-out: TabAlpha(b) gives the alpha regularization parameter for
                  the band b
   DataModel = in: true if we have a model
   MR_Model = in: pointer to a multiresolution model or NULL (if none)
   DataSNR  = in: true if alpha parameter are modified using the SNR
                    of the data
   UseSupport= in: true if alpha parameter are modified using the multiresolution
                        support
  Warning: works well, when each wavelet can be modelized by a Gaussian		    
*/	
{
   CMemWave MemWObj;

   switch (TypeFilter)
   {
      case DEF_MEM_ALPHA_CST:
        mw_fm1(MR_Data, NoiseModel,TabAlpha, MemWObj, DataModel, MR_Model, DataSNR, UseSupport);
 	break;
      case DEF_MEM_ALPHA_OPT:
        mw_fm2( MR_Data, NoiseModel,  TabAlpha, MemWObj, DataSNR, 
	          DataModel, MR_Model, ConvgParam, MaxIter, PositivConst, Verbose, UseSupport);
        break;  
       case DEF_MEM_ALPHA_OPT_SCALE:
        mw_fm3( MR_Data, NoiseModel,  TabAlpha, MemWObj, DataSNR, 
	          DataModel, MR_Model, ConvgParam, MaxIter, Verbose, 
		  True, PositivConst, UseSupport);
        break; 
      default: 
        cerr << "Error: unknow type of filtering ... " << endl;
	exit(-1);
	break;
   }
}
 
/****************************************************************************/

void mw_fm1(MultiResol & MR_Data, MRNoiseModel & NoiseModel, 
            fltarray &TabAlpha,  CMemWave & MemWObj, 
	    Bool UseModel,  MultiResol *MR_Model, 
	    Bool UseDataSNR, Bool UseSupport)
/* Apply the multiscale entropy to the multiresolution data set MR_Data
   MR_Data = in-out: multiresolution transform of an image
   NoiseModel = in: noise modeling of the image in the wavelet space
   TabAlpha = in: TabAlpha(b) gives the alpha regularization parameter for
                  the band b
   MemWObj = in: class for entropy filtering
   UseModel = in: true if we have a model
   MR_Model = in: pointer to a multiresolution model or NULL (if none)
   UseDataSNR = in: true if alpha parameter are modified using the SNR
                    of the data
   UseSupport= in: true if alpha parameter are modified using the multiresolution
                        support
*/		    		  
    
{
   float AlphaP,Val;
   int b,i,j;
   float Model=0.;

   // apply the multiscale entropy on each band
   for (b = 0; b < MR_Data.nbr_band()-1; b++)
   {
      int s = MR_Data.band_to_scale(b);
		    
      // for all wavelet coefficients
      for (i = 0; i < MR_Data.size_band_nl(b); i++)
      for (j = 0; j < MR_Data.size_band_nc(b); j++)
      {
          if (s <  NoiseModel.FirstDectectScale) MR_Data(b, i,j) = 0;
	  else
	  {
 	    float Sigma = NoiseModel.sigma(b,i,j);
            float RegulParam = TabAlpha(b);
	    Val = MR_Data(b,i,j);
	  
	    // use the input model for the wavelet coefficients
	    if ((UseModel != True) || (MR_Model == NULL)) Model = 0.;
	    else  Model = (*MR_Model)(b,i,j);
 	  
	    // Modify Alpha using the SNR of the data
            if ((UseDataSNR == True) 
	    && ((NoiseModel.OnlyPositivDetect == False) || (Val >= 0)))
            {
	     // SNR = MAX [ ABS(data) / (N*Sigma), 1 ]
             AlphaP = ABS(Val / (Sigma* NoiseModel.NSigma[s]));
             if (AlphaP > 1) AlphaP = 1;
	     
	     // if AlphaP = 0, J = h_n ==> Solution = model
	     // when RegulParam < 0, the filtering routine set to
	     // the model the solution
             if (AlphaP < FLOAT_EPSILON) RegulParam = -1;
	     else RegulParam *= (1-AlphaP)/AlphaP;
            }
	    else if ((UseSupport == True)  
	            &&  (NoiseModel(b,i,j) == True)) RegulParam = 0;
	  
	    // correct the wavelet coefficient
 	    if (RegulParam > 0) MR_Data(b, i,j) = MemWObj.filter(Val, RegulParam, Sigma, Model);
	    else MR_Data(b, i,j) = Val;
	  
// 	  if ((b==0) && (i==59) && (j==47))
// 	  {
// 	     if (NoiseModel(b,i,j) == True) cout << "Support coef ... " << endl;
// 	     if (UseSupport == True) cout << "Use Support" << endl;
//  	     cout << "Data = " << MR_Data(b, i,j) << " Corr = " << Val << endl;
//           }
         }
       }
    }                                
}
 
/***************************************/


void mw_fm2(MultiResol &MR_Data, MRNoiseModel & NoiseModel, 
        fltarray &TabAlpha, CMemWave & MemWObj, Bool DataSNR, Bool DataModel, MultiResol *MR_Model,
	float ConvgParam, int MaxIter, Bool PositivConst, Bool Verbose, Bool UseSupport)
/* Apply the multiscale entropy to the multiresolution data set MR_Data
   The Parameter Alpha is estimated iteratively (by dichotomy) in order to have
   a residu compatible with the noise.
   
   MR_Data = in-out: multiresolution transform of an image
   NoiseModel = in: noise modeling of the image in the wavelet space
   TabAlpha = in: TabAlpha(b) gives the alpha regularization parameter for
                  the band b
   MemWObj = in: class for entropy filtering
   UseModel = in: true if we have a model
   MR_Model = in: pointer to a multiresolution model or NULL (if none)
   UseDataSNR = in: true if alpha parameter are modified using the SNR
                    of the data
   ConvgParam = in: Convergence parameter
   MaxIter = in: maximum number of iterations
   PositivConst = in: apply the positivy constraint
   UseSupport= in: true if alpha parameter are modified using the multiresolution
                        support
*/	
{
   int Nbr_Band = NoiseModel.nbr_band()-1;
   int Nbr_Plan = NoiseModel.nbr_scale();
   int Nl = MR_Data.size_ima_nl();
   int Nc = MR_Data.size_ima_nc();
   type_transform Transform = NoiseModel.type_trans();
   MultiResol MR_Sol(Nl,Nc,Nbr_Plan,Transform,"MR_Sol");
   Ifloat Data(Nl,Nc, "data");
   Ifloat Ima(Nl,Nc, "Ima");
   double RegulMin, RegulMax, SigmaNoise, Delta;
   int b,Iter=0;
   float RegulParam;
   float AlphaUser = TabAlpha(0);
   
   // we need the raw data in order to be able to calculate the residual
   // image: resi = data - solution
   MR_Data.recons(Data);
   
   // initialization
   MR_Sol.band(Nbr_Band) = MR_Data.band(Nbr_Band);
   RegulMin=0.0;                // minimum alpha value
   RegulMax=DEF_MAX_MEM_ALPHA;  // maximum alpha value
   
   cout <<  "Sigma Noise= " <<  NoiseModel.SigmaNoise <<  endl;
   do 
   {
      Iter ++;
      // copy the raw wavelet coef. into the solution array
      for (b = 0; b < Nbr_Band; b++)  MR_Sol.band(b) = MR_Data.band(b);
      
      // RegulParam = alpha parameter new calculation
      RegulParam = (RegulMin + RegulMax)/2;
      
      // put it in TabAlpha
      for (b = 0; b <  Nbr_Band; b++) TabAlpha(b) = RegulParam;
      
      // correct the wavelet coefficients
      mw_fm1(MR_Sol, NoiseModel,TabAlpha, MemWObj, DataModel, MR_Model, False, False);
      
      // reconstruct the filtered image
      MR_Sol.recons(Ima);
      
      // Apply positivuty constraint
      if (PositivConst == True) threshold(Ima);
      
      // SigmaNoise = sigma(residual)
      SigmaNoise = sigma(Data-Ima);  
      
      // modify RegulMin,RegulMax by dichotomy 
      if (SigmaNoise > NoiseModel.SigmaNoise) RegulMax=RegulParam; 
      else RegulMin = RegulParam;
      
      // distance between  sigma(residual) and sigma(noise)
      Delta = ABS(SigmaNoise - NoiseModel.SigmaNoise);

      if (Verbose == True) 
         cout << Iter << ": Sigma Resi = " << SigmaNoise << " Regul. Param = " << RegulParam << endl;

     // convergence test
   } while( (Delta>ConvgParam) && (Iter < MaxIter));
   
   // the calculate alpha is 
   RegulParam = (RegulMin + RegulMax)/2;
   for (b = 0; b <  Nbr_Band; b++) TabAlpha(b) = RegulParam*AlphaUser;
   
   if (Verbose == True)  
      cout << "Optimal Regul. Param = " << RegulParam << endl;
   
   // final correction of the wavelet coefficients
   mw_fm1(MR_Data, NoiseModel,TabAlpha, MemWObj, DataModel, MR_Model,DataSNR,UseSupport);
}
	       
/***************************************/

void mw_fm3(MultiResol &MR_Data, MRNoiseModel & NoiseModel, 
        fltarray &TabAlpha, CMemWave & MemWObj, Bool DataSNR, Bool DataModel, MultiResol *MR_Model,
	float ConvgParam, int MaxIter, Bool Verbose, 
	Bool RecEachIter, Bool PosConstraint, Bool UseSupport)
/* Apply the multiscale entropy to the multiresolution data set MR_Data
   The Parameter Alpha is estimated iteratively at scale separately
   (by dichotomy) in order to have a residu compatible with the noise
   at all the scales
   
   MR_Data = in-out: multiresolution transform of an image
   NoiseModel = in: noise modeling of the image in the wavelet space
   TabAlpha = in: TabAlpha(b) gives the alpha regularization parameter for
                  the band b
   MemWObj = in: class for entropy filtering
   UseModel = in: true if we have a model
   MR_Model = in: pointer to a multiresolution model or NULL (if none)
   UseDataSNR = in: true if alpha parameter are modified using the SNR
                    of the data
   ConvgParam = in: Convergence parameter
   MaxIter = in: maximum number of iterations
   PositivConst = in: apply the positivy constraint
   UseSupport= in: true if alpha parameter are modified using the multiresolution
                        support
*/	
{
   int Nbr_Band = NoiseModel.nbr_band()-1;
   int Nbr_Plan = NoiseModel.nbr_scale();
   int Nl = MR_Data.size_ima_nl();
   int Nc = MR_Data.size_ima_nc();
   type_transform Transform = NoiseModel.type_trans();
   MultiResol MR_Sol(Nl,Nc,Nbr_Plan,Transform,"MR_Sol");
   Ifloat ImaSol(Nl,Nc, "Ima");
   double  SigmaNoise;
   int b,i,j, Iter=0;
   fltarray  RegulMin( Nbr_Band ); 
   fltarray  RegulMax( Nbr_Band );    
   fltarray  TabDelta( Nbr_Band );    
   float AlphaUser = TabAlpha(0);

   // initialization
   MR_Sol.band(Nbr_Band) = MR_Data.band(Nbr_Band);
 
  // initialization
   for (b = 0; b < Nbr_Band; b++)
   {
      RegulMin(b) = 0;                 // minimum alpha value
      RegulMax(b) = DEF_MAX_MEM_ALPHA; // maximum alpha value
      if (MR_Data.band_to_scale(b) == 0) TabAlpha(b)= 5;
      else TabAlpha(b)= 1.;
      TabDelta(b) = 0.;
   }   
   if (Verbose == True) 
   {
         if (RecEachIter == True) cout << "Reconstruction at each iteration ... " << endl;
	 if (UseSupport == True) cout << "Use the multiresolution support ... " << endl;
	 if (DataSNR == True) cout << "Protect high coeff from the regularization ... " << endl;
	 cout << " AlphaUser = " << AlphaUser << endl;
   }

   do 
   {
      Iter ++;
      if (Verbose == True) cout << "Iter " << Iter << endl;
      
      // copy the raw wavelet coef. into the solution array
      for (b = 0; b <= Nbr_Band; b++)  MR_Sol.band(b) = MR_Data.band(b);
      
      // RegulParam = alpha parameter new calculation
      for (b = 0; b <  Nbr_Band; b++) 
                 TabAlpha(b) = (RegulMin(b) + RegulMax(b))/2;
      
      // correct the wavelet coefficients
      mw_fm1(MR_Sol, NoiseModel,TabAlpha, MemWObj, DataModel, MR_Model, DataSNR, UseSupport);

      // transform the result, apply the positivity contrain
      // and reconstruct
      if (RecEachIter == True)
      {
         MR_Sol.recons(ImaSol);
	 INFO(ImaSol, "RES");
	 if (PosConstraint == True) threshold(ImaSol);
	 MR_Sol.transform(ImaSol);
	 // MR_Sol.band(Nbr_Band) = MR_Data.band(Nbr_Band);
      }
      
      // new estimation of the Alpha paramters
      for (b = 0; b <  Nbr_Band; b++) 
      {
         int Nlb = MR_Data.size_band_nl(b);
	 int Ncb = MR_Data.size_band_nc(b);
         SigmaNoise = 0.;
	 long Nc = 0;
         for (i = 0; i < Nlb; i++)
         for (j = 0; j < Ncb; j++)
         {
	    if ((UseSupport == False ) || 
	        ((UseSupport == True) && (NoiseModel(b,i,j) == False)))
	    {
  	       double Resi = MR_Data(b,i,j) - MR_Sol(b,i,j);
 	      double SigmaCoef = NoiseModel.sigma(b,i,j);
	      double ExpectVariance =  SigmaCoef*SigmaCoef;
	      if (ExpectVariance > FLOAT_EPSILON)
 	          SigmaNoise += (Resi*Resi) / ExpectVariance;
              else SigmaNoise += 1 + Resi*Resi;
	      Nc++;
	   }
 	 }
	 if (Nc > 0) SigmaNoise = sqrt(SigmaNoise/(double) (Nc));
	 else  SigmaNoise = 1.;
	 if (SigmaNoise >= 1)  RegulMax (b) =  TabAlpha(b);
	 else  RegulMin(b) = TabAlpha(b);      
  
	 TabDelta(b) = ABS(SigmaNoise - 1.);
	 if (Verbose == True)
 	   cout << "   band " << b+1 << " Normalized sigma " << SigmaNoise << " Alpha = " <<  TabAlpha(b) << "  Delta = " << TabDelta(b) << endl;      
      }
     // convergence test
   } while( (TabDelta.max() > ConvgParam) && (Iter < MaxIter));

   // the calculate alpha is 
   for (b = 0; b <  Nbr_Band; b++) 
                 TabAlpha(b) = (RegulMin(b) + RegulMax(b))/2*AlphaUser;    

   if (Verbose == True)
     for (b = 0; b <  Nbr_Band; b++) 
      cout << "band " << b+1 << " Optimal Regul. Param = " <<  TabAlpha(b) << endl;
    
   // final correction of the wavelet coefficients
   mw_fm1(MR_Data, NoiseModel, TabAlpha, MemWObj, DataModel, MR_Model,DataSNR,UseSupport);
   MR_Data.recons(ImaSol);
}
	       
/***************************************/

// same as fm3 but alpha parameter is not estimated by dichotomy, but
// fix step gradient
// void mw_fm3_bis(MultiResol &MR_Data, MRNoiseModel & NoiseModel, 
//                fltarray &TabAlpha, CMemWave &MemWObj, Bool DataSNR,  
// 	       Bool DataModel, MultiResol *MR_Model,
// 	       float ConvgParam, int MaxIter)
// /* Apply the multiscale entropy to the multiresolution data set MR_Data.
//    The alpha parameter is different at each scale, and all alpha are
//    estimated iteratively in order to have a residu compatible with the noise
//    at all the scale.
//    
//    MR_Data = in-out: multiresolution transform of an image
//    NoiseModel = in: noise modeling of the image in the wavelet space
//    TabAlpha = in: TabAlpha(b) gives the alpha regularization parameter for
//                   the band b
//    MemWObj = in: class for entropy filtering
//    UseModel = in: true if we have a model
//    MR_Model = in: pointer to a multiresolution model or NULL (if none)
//    UseDataSNR = in: true if alpha parameter are modified using the SNR
//                     of the data
//    ConvgParam = in: Convergence parameter
//    MaxIter = in: maximum number of iterations
//    
//    minimization of alpha by fixed step gradient
//    ccl: dichotomy seem better adapted.
// */		       
// {
//    int Nbr_Band =  MR_Data.nbr_band()-1;
//    int Nl = MR_Data.size_ima_nl();
//    int Nc = MR_Data.size_ima_nc();
//    type_transform Transform = NoiseModel.type_trans();
//    int Nbr_Plan = MR_Data.nbr_scale();
//    MultiResol MR_Sol(Nl,Nc,Nbr_Plan,Transform,"MR_Sol");
//    MultiResol MR_Resi(Nl,Nc,Nbr_Plan,Transform,"MR_Resi");
//    
//    Ifloat Data(Nl,Nc, "data");
//    Ifloat Ima(Nl,Nc, "data");
//    double SigmaNoise;
//    int b,i,j,Iter=0;
//    fltarray TabDelta( Nbr_Band ); 
//    fltarray TabStep( Nbr_Band ); 
//    
//    // we need the raw data in order to be able to calculate the residual
//    // image: resi = data - solution   
//    MR_Data.recons(Data);
//    
//    // initialization
//    for (b = 0; b < Nbr_Band; b++)
//    {
//        if (MR_Data.band_to_scale(b) == 0) TabAlpha(b)= 5;
//        else TabAlpha(b)= 1.;
//        TabDelta(b) = 0.;
//    }
//    
//  
// 
//    do 
//    {
//       Iter ++;
//       cout << "Iter " << Iter << endl;
// 
//       // copy the raw wavelet coef. into the solution array   
//       for (b = 0; b <= Nbr_Band; b++)  MR_Sol.band(b) = MR_Data.band(b); 
//             
//       // correct the wavelet coefficients
//       mw_fm1(MR_Sol, NoiseModel,TabAlpha, MemWObj, DataModel, MR_Model,DataSNR);
//       
//       // calculate the residual image
//       //MR_Sol.recons(Ima);
//       //Ima = Data - Ima;
//       
//       // wavelet transform of the residual
//       //MR_Resi.transform(Ima);
//       
//       // new estimation of the Alpha paramters
//       for (b = 0; b <  Nbr_Band; b++) 
//       {
//          int Nlb = MR_Data.size_band_nl(b);
// 	 int Ncb = MR_Data.size_band_nc(b);
//          SigmaNoise = 0.;
// 	 
// 	 // standard deviation of the residual in the band b
// 	 // and Alpha paramter new estimation
// 	 float SumNum = 0.;
// 	 // float SumDen = 0.;
//          for (i = 0; i < Nlb; i++)
//          for (j = 0; j < Ncb; j++)
//          {
// 	    float dhs,dhn,Coef = MR_Sol(b,i,j);
// 	    float Resi= MR_Data(b,i,j) - Coef;
// 	    //float Resi= MR_Resi(b,i,j);
// 	    float SigmaCoef = NoiseModel.sigma(b,i,j);
// 	    float ExpectVariance =  SigmaCoef*SigmaCoef;
//  	    SigmaNoise += (Resi*Resi) / ExpectVariance;
// 	    
// 	    dhn = MemWObj.grad_hn(Coef);
// 	    dhs = MemWObj.grad_hs(Resi);
// 	    SumNum += dhn/ExpectVariance * Resi;
// 	    //SumNum += dhn/ExpectVariance*(-Resi-dhs);
// 	    //SumDen += dhn*dhn/ExpectVariance;
// 	 }
// 	 SigmaNoise = sqrt(SigmaNoise/(Nlb*Ncb));
// 	 if (SigmaNoise >= 1)  TabAlpha(b) -= 0.2*SumNum/(Nlb*Ncb);
// 	 else TabAlpha(b) += 0.2*SumNum/(Nlb*Ncb);
// 	 if ( TabAlpha(b)  <= 0.) TabAlpha(b) = 0.;
// 	 //if (SumDen > FLOAT_EPSILON)  SumNum /= SumDen;
// 	 //else TabAlpha(b) = 1.;
//  	 //if (SumNum  <= 0.) TabAlpha(b) /= 2.;
// 	 //else TabAlpha(b) = SumNum; 
//   	 
// 	 // distance between  sigma(residual) and sigma(noise)
//  	 TabDelta(b) = ABS(SigmaNoise - 1.);
// 	 
// 	 cout << "   band " << b+1 << " Normalized sigma " << SigmaNoise << " Alpha = " <<  TabAlpha(b) << "  Delta = " << TabDelta(b) << endl;
//       }  
//       
//    // convergence test
//    } while( (TabDelta.max() > ConvgParam) && (Iter < MaxIter));
//    
//       
// 
//    for (b = 0; b <  Nbr_Band; b++) 
//       cout << "band " << b+1 << " Optimal Regul. Param = " <<  TabAlpha(b) << endl;
// 
//    // final correction of the wavelet coefficients
//    mw_fm1(MR_Data, NoiseModel,TabAlpha, MemWObj, DataModel, MR_Model,DataSNR);
// }      

/****************************************************************************/
