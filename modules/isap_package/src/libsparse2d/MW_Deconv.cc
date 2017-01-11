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
**    Date:  23/12/98 
**    
**    File:  MW_Deconv.cc
**
*******************************************************************************
**
**    DESCRIPTION  Image Deconvolution using the multiscale entropy
**    ----------- 
**                 
**
******************************************************************************/

#include "IM_Obj.h"
#include "IM_Math.h"
// #include "IM_IO.h"
#include "MR_Obj.h"
#include "MR_NoiseModel.h"
#include "CErf.h"
#include "CMem.h"
#include "IM_Deconv.h"
#include "MR_Deconv.h"
#include "MW_Deconv.h"
#include "MW_Filter.h"


/*****************************************************************************/
 
void MWDeconv::mw_grad(Ifloat & Gradient)
// Apply the multiscale entropy to the multiresolution data set MR_Obj  
{
   float Alpha=0.,Alphab;
   float Sigma;
   int Nbr_band = MR_Obj.nbr_band();
   int b,i,j;
   
   // apply the multiscale entropy on each band
   MR_Obj.transform(Obj); 
    
   for (b=0; b < Nbr_band-1; b++)
   {
       // for all wavelet coefficients
      for (i=0; i < MR_Obj.size_band_nl(b); i++)
      for (j=0; j < MR_Obj.size_band_nc(b); j++)
      {      	  
	  float Val = MR_Obj(b,i,j);
          Sigma =  ModelData->sigma(b,i,j);
  	  Alphab = RegulParam*TabAlpha(b);
 	  //Modify Lambda using the SNR of the data
	  switch(TypeRegul)
	  {
	     case DEC_NO_REGUL: 
	            Alpha = 0.; 
		    break;
	     case DEC_REGUL_NO_PROTECT: 
	            Alpha = Alphab; 
		    break;
	     case DEC_REGUL_PROBA: 
	           Alpha = Alphab*MR_Data(b,i,j); 
		   break;
	     case DEC_REGUL_SUPPORT: 
	           //if ((*ModelData)(b,i,j) == True) Alpha = 0.;
		   //else 
		   Alpha = Alphab* MR_RegulObj(b,i,j); 
		   break;
	     case DEC_REGUL_PROBA_SUPPORT:
	           //if ((*ModelData)(b,i,j) == True) Alpha = 0.;
		   //else Alpha =  Alphab*MR_Data(b,i,j);
		   Alpha =  Alphab * MR_Data(b,i,j) * MR_RegulObj(b,i,j);
		   break;
	  }
          // correct the wavelet coefficients 
	  if (DecMethod == DEC_MR_MEM_NOISE)
	  {
	    int Sign = (Val >= 0) ? 1: -1;	     
	    MR_Obj(b,i,j) = MemWObj.grad_hn(ABS(Val)/Sigma)*Sigma*Sign*Alpha;
 	  }	  
          else MR_Obj(b,i,j) = Val*Alpha;
      }
      
   }
   MR_Obj.band(Nbr_band-1).init();
   MR_Obj.recons(Buff);
   Gradient -= Buff;
    

}

/*********************************************************************/

float MWDeconv::mw_find_optim(Ifloat &Gradient)
// MR_Obj contains the wavelet transform of the current solution
{
    int b,i,j;
    int Nbr_band = MR_Obj.nbr_band();
    float hcroise,hgrad, coef_num,coef_den, Wg, ValRet;
    
    MultiResol MR_Grad(Resi.nl(), Resi.nc(), MR_Obj.nbr_scale(), 
                       MR_Obj.Type_Transform, "MR_Grad");
    MR_Grad.transform(Gradient);
    MR_Obj.transform(Obj); 
    psf_convol (Gradient, Psf_cf, Buff);
    coef_den = energy(Buff);
    Buff *= Resi;
    coef_num = flux(Buff);
    
    hcroise=0;
    hgrad=0;
    for (b = 0; b < Nbr_band-1; b++)
    for (i = 0; i <  MR_Grad.size_band_nl(b); i++)
    for (j = 0; j <  MR_Grad.size_band_nc(b); j++)
    {
	Wg = (TypeRegul == DEC_REGUL_PROBA) ? MR_Data(b,i,j): 1.;
	Wg *= TabAlpha(b)*RegulParam;
	hcroise += MR_Obj(b,i,j) * MR_Grad(b,i,j) * Wg;
	hgrad   += MR_Grad(b,i,j) * MR_Grad(b,i,j)* Wg;
    }
    coef_num -=  hcroise;
    coef_den +=  hgrad;
    cout  << "Coef num = " <<  coef_num << " Den = " << coef_den;
     cout << " param = " << coef_num /  coef_den << endl; 
    ValRet  =  coef_num /  coef_den;
    return ValRet;
}


/***************************************************************************/

void MWDeconv::compute_mem_resi()
{
   int b,i,j;
   
    psf_convol (Obj, Psf_cf, Resi);
    for (i=0; i< Nl*Nc; i++)  Resi(i) = Imag(i) - Resi(i);
    SigmaResi = sigma(Resi);
    
    switch (TypeWeight)
    {
       case DEC_WEIGHT_PROBA:
           MR_Obj.transform(Resi);
           for (b = 0; b < NbrBand-1; b++)
           for (i = 0; i < MR_Obj.size_band_nl(b); i++)
           for (j = 0; j < MR_Obj.size_band_nc(b); j++) 
                       MR_Obj(b,i,j) *= (1. - MR_Data(b,i,j));
           MR_Obj.recons(Resi);
	   break;
      case DEC_WEIGHT_SUPPORT:
           MR_Obj.transform(Resi);
	   ModelData->threshold(MR_Obj);
           MR_Obj.recons(Resi);
	   break;
     case DEC_NO_WEIGHT:
           break;
     default:
           cerr << "Error: unknown data weighting method ... " << endl;
	   exit(-1);
	   break;
     }
}

/****************************************************************************/

float MWDeconv::sigma_obj(float Sigma)
//    The routine estimate the regularization paramter per band
//    1) Create a random Gaussian noise: N
//    2) Convolve it by the PSF: I = P*N
//    3) Take its wavelet transform: W = WT(I)
//    4) Calculate the standard deviation in each band j: TabSigma[j]
//    5) MaxSigma = max(TabSigma)
//    6) The regularisation parameter at each band j is given by
//         TabAlpha[j] = Sigma_j / TabSigma[i]
//    7) Renormalisation: TabAlpha[j] = TabAlpha[j] / MaxSigma
{
   int b;
   Ifloat Dat(Nl, Nc, "noise");
   TabSigmaObj.alloc(NbrBand);
   im_noise_gaussian (Dat, 1.);
   psf_convol_conj (Dat, Psf_cf);
   MR_Obj.transform(Dat);
   
   for (b = 0; b < NbrBand-1; b++)
      TabSigmaObj (b) = Sigma * MR_Obj.band_norm(b) / sigma( MR_Obj.band(b));
      
   float Max =  TabSigmaObj.max();
   for (b = 0; b < NbrBand-1; b++) TabAlpha(b) = TabSigmaObj(b) / Max;
   
   if (Verbose == True)
   {
      for (b = 0; b < NbrBand-1; b++)
      {
        cout << "Band " << b << "  alpha = " << TabAlpha(b) << endl;
        // cout << "       Sigma obj  = " << TabSigmaObj(b) << "  P2 = " 
	//        << TabSigmaObj(b)*TabSigmaObj(b)<< endl;
      }
   }
   return (Sigma / sigma(Dat));
}

/****************************************************************************/

void MWDeconv::get_mr_support_first_guess(type_deconv Deconv, int Niter)
{
    // save real parameters
    type_deconv RealDecMethod = DecMethod;
    int RealMaxIter = MaxIter;
    float RealIterCvg = IterCvg;
    Bool RealVerb = Verbose;
    Bool RealOptim = OptimParam;
    // int b,i,j;
    //fltarray TabNsig[NbrBand];
    
    // Change some paramters 
    IterCvg = 1.;
    DecMethod = Deconv;
    MaxIter = Niter;
    Verbose = False;
    OptimParam = False;
    //float NSig = 5.;
    //for (b=0; b < NbrBand; b++) TabNsig[b] =  (ModelData->NSigma)[b];
    //for (b=0; b < NbrBand; b++)  (ModelData->NSigma)[b] = NSig;
    // run the deconvolution
    im_iter_deconv();
       
    // restore real parameters
    DecMethod = RealDecMethod;
    MaxIter = RealMaxIter;
    IterCvg = RealIterCvg;
    Verbose = RealVerb;
    OptimParam = RealOptim;
    //for (b=0; b < NbrBand; b++) (ModelData->NSigma)[b] = TabNsig[b];

    if (Verbose == True)
    {
        INFO(Obj, "MR Van-Citter Sol: ");
    }
}

/****************************************************************************/

void MWDeconv::find_support_obj()
{
   int b,j,i;

   float Sig = sigma(Obj);
    
   MR_Obj.transform(Obj);
   // io_write_ima_float("sol1.fits", Obj);
   for (b = 0; b < NbrBand-1; b++)
   {
       int Scale = MR_Obj.band_to_scale(b);
       int Band = (int) pow ((double) 2., (double)(Scale+1));
       
       Buff.resize(MR_Obj.size_band_nl(b), MR_Obj.size_band_nc(b));
       float Level = NSigmaObj*Sig*MR_Obj.band_norm(b);
       for (i = 0; i < MR_Obj.size_band_nl(b); i++)
       for (j = 0; j < MR_Obj.size_band_nc(b); j++)
       {   
          if (ABS(MR_Obj(b,i,j)) > Level) Buff(i,j) = 0.; 	  
          else Buff(i,j) =  1;
  	}
	smooth_bspline(Buff, MR_RegulObj.band(b), I_ZERO, Scale);
	for (i = Band; i < MR_Obj.size_band_nl(b)-Band; i++)
        for (j = Band; j < MR_Obj.size_band_nc(b)-Band; j++)
                 if (ABS(MR_Obj(b,i,j))  > Level) MR_RegulObj(b,i,j) = 0.; 	  
    }
    Buff.resize(Nl,Nc);
    // MR_RegulObj.write("sup.mr");
}
	  
/****************************************************************************/

void MWDeconv::mw_vaguelet(Bool UseSupport, Ifloat *FirstGuess, Ifloat *ICF)
{
   int b,i,j;
   CMemWave MemWObj;
   Bool RecEachIter = False;
   Bool DataSNR = False;
   
   WaveFilterResi = True;
   NbrBand = ModelData->nbr_band();
   Nscale = ModelData->nbr_scale();
   Transform= ModelData->type_trans();
   TabAlpha.alloc(NbrBand);
   init_deconv(FirstGuess,  ICF);
    
   Ifloat NoiseSimu(Nl,Nc,"simu");
   va_inverse_data(NoiseSimu);
   // io_write_ima_float("xx_noise.fits", NoiseSimu);
   
   Resi = Obj;
//    if (UseSupport == True)
//    {
//       MR_Data.transform(Obj);
//       if (MR_Data.Set_Transform != TRANSF_MALLAT)
//       {
//           type_filter Filter = (MR_Data.Set_Transform == TRANSF_MALLAT) ?
//                            FILTER_THRESHOLD: FILTER_ITER_THRESHOLD;
//           MRFiltering CFilter(*ModelData, Filter);
//           NoiseSimu=Obj;
//           CFilter.Sup_Set = False;
//           CFilter.Verbose = Verbose;
//           CFilter.KillLastScale = KillLastScale;
//           CFilter.mr_support_filter(NoiseSimu, MR_Data, Resi);
//       }
//       else
//       {
//           ModelData->threshold(MR_Data);
// 	  MR_Data.recons(Resi);
//       }
//       Obj -= Resi;
//    }
    
   // io_write_ima_float("xx_inv.fits", Obj);
   for (b = 0; b < NbrBand-1; b++) 
                            TabAlpha(b) = (RegulParam <= 0) ? 1.: RegulParam;

    // MR_Data.transform(Obj);   
    // apply the filtering
    Bool DataEdge = False;
    MultiResol  *MR_Model = NULL;
    float FilCvgParam= DEF_MEM_FILT_CONGER;
    int FilNbrIter = DEF_MEM_FILT_MAX_ITER;
    
    MRNoiseModel MD(NOISE_CORREL, Nl, Nc, Nscale, Transform); 
    MD.RmsMap = NoiseSimu;
    for (b=0; b < NbrBand; b++) MD.NSigma[b]= ModelData->NSigma[b];
    for (b=0; b < NbrBand; b++) MD.TabEps[b]= ModelData->TabEps[b];
    MD.OnlyPositivDetect = ModelData->OnlyPositivDetect;
    MD.SupIsol = ModelData->SupIsol;
     
    MD.model(Obj, MR_Data);
    // MD.threshold(MR_Data);
    mw_fm3(MR_Data, MD,  TabAlpha, MemWObj, DataSNR, 
	          DataEdge, MR_Model, FilCvgParam, FilNbrIter, Verbose,
		  RecEachIter, PositivConstraint, UseSupport);
    // MR_Data.write("xx.mr");
    int Lb = MR_Data.nbr_band()-1;
    if (KillLastScale == True) MR_Data.band(Lb).init();
    
    MR_Data.recons(Obj);
    // Resi -= Obj;
    // io_write_ima_float("xx_resi.fits", Resi);
    // if (UseSupport == True) Obj += Resi;
    
    if (PositivConstraint == True) threshold(Obj);
    if (KillLastScale == True)
    {
       MultiResol MR_Resi(Nl, Nc, MD.nbr_scale(), MD.type_trans(), "MR_Data");		 
       for (int Iter=0; Iter < 10; Iter++)	
       {
           MR_Resi.transform(Obj);
           for (b = 0; b < MR_Data.nbr_band()-1; b++)
           for (i = 0; i < MR_Data.size_band_nl(b); i++)
           for (j = 0; j < MR_Data.size_band_nc(b); j++)
                 MR_Resi(b,i,j) = MR_Data(b,i,j) - MR_Resi(b,i,j);
	   MR_Resi.band(Lb).init();	 
           MR_Resi.recons(Resi);
           Obj += Resi;
           if (PositivConstraint == True) threshold(Obj);
       }
    }
    
    compute_resi();
    if (GaussConv == True)  convol_gauss(Obj, Fwhm);
    else if (UseICF == True) psf_convol(Obj, Ima_ICF, Obj);
}

/****************************************************************************/

void MWDeconv::mw_deconv(Ifloat *FirstGuess, Ifloat *ICF)
{
    int b,i,Iter=0;
    Bool Stop = False;
    float  OldFonc=1e+15, Fonc, Delta;
    float Cvgparam = IterCvg;
//     if ((OptimParam == True) && (DecMethod == DEC_MR_MEM_NOISE))
//     {
//        cerr << "Optimization is not implemented when using noise entropy ... " << endl;
//        exit(-1);
//     }
 
    
    // initialization
    // ----------------
    WaveFilterResi = True;
    NbrBand = ModelData->nbr_band();
    Nscale = ModelData->nbr_scale();
    Transform= ModelData->type_trans();
    TabAlpha.alloc(NbrBand);
    // Border = I_CONT;
    init_deconv(FirstGuess,ICF);
    NbIter_SupportCalc = DEF_NB_ITER_SUP_CAL;
    if ((RegulParam <= 0.) || (TypeRegul == DEC_NO_REGUL))
    {
       RegulParam = 0.;
       TypeRegul = DEC_NO_REGUL;
    }
  
   
    int Np = Nl*Nc;
    Ifloat Gradient(Nl, Nc, "Info");
    Buff.alloc(Nl, Nc, "Buffer");
    MR_Obj.alloc(Nl, Nc, Nscale, Transform, "MR_Obj");
    MR_Obj.Border = Border; 
    if ((TypeRegul == DEC_REGUL_PROBA_SUPPORT) ||
        (TypeRegul == DEC_REGUL_SUPPORT)) 
    {
         MR_RegulObj.alloc(Nl, Nc, Nscale, Transform, "MR_Obj");
	 MR_RegulObj.Border = Border;
    }
    Noise_Ima = ModelData->SigmaNoise;
     
    // Noise_Obj = sigma_obj(Noise_Ima);
    for (b = 0; b < NbrBand-1; b++) TabAlpha(b) = 1.;
    float FluxImag=0.;
    
    if (NormFlux == True) FluxImag = flux(Imag);
    
   // Find the first guess using the regularized Van Cittert Method  
    if (GetMR_FirstGuess == True) get_mr_support_first_guess(DEC_MR_GRADIENT,50);
    
    if ((TypeWeight == DEC_WEIGHT_PROBA) || (TypeRegul == DEC_REGUL_PROBA)
          || (TypeRegul == DEC_REGUL_PROBA_SUPPORT))
    {
       ModelData->prob_noise(MR_Data);
    }
    
    
    if (Verbose == True)
    {
       cout << "Image size = " << Nl << " " << Nc << endl;
       cout << "Noise_Ima = " << Noise_Ima << endl;
       // cout << "Sigma obj = " << Noise_Obj << endl;
    }     
      
    // Residual estimation 
    compute_mem_resi();

    // if (Verbose == True) cout << "start iterate" << endl;
    do
    {
        // Calculate the multiplicate term in Gradient's iteration 
        // io_write_ima_float("xx_resi.fits",  Resi);
        psf_convol_conj (Resi, Psf_cf, Gradient);

        // calculates the gradient
   	if (TypeRegul != DEC_NO_REGUL) mw_grad(Gradient);
	    
  	    // Functional gradient
	    //io_write_ima_float("xx_grad.fits", Gradient);
	    //if ((Verbose == True) && (Iter == MaxIter-1))
            //{
            //  i = 670;
            //  cout << "P(10,10) = " << Gradient(i) << endl;
            //}
 	    //io_write_ima_float("xx_fgrad.fits",  Gradient);
 	
	// Supress the last scale in the Gradient
        if (KillLastScale == True)
        {
            int Nbr_band = MR_Obj.nbr_band();
            MR_Obj.transform(Gradient);
            MR_Obj.band(Nbr_band-1).init();
            MR_Obj.recons(Gradient);
        }
	  	
        // convergence parameter
	if (OptimParam == True) 
	{
	   // Cvgparam = mw_find_optim(Gradient);
           // PB: mw_find_optim does not work correctly
	   //     it replaces by a X2 optimization
	   
	   // Xi2 optimization
           psf_convol(Gradient, Psf_cf, Buff);
           Cvgparam = flux(Buff*Resi)/energy(Buff);
           if (Cvgparam > 10) Cvgparam = 10.;
     	}

        // calculate the next object estimation with positivity  constraint
        for (i = 0; i < Np; i++) Obj(i) += Cvgparam * Gradient(i);		    

        if (PositivConstraint == True) threshold(Obj);
        if (NormFlux == True)
        {
          float FluxObj = flux(Obj);
	  FluxObj = FluxImag / FluxObj;
	  for (i= 0; i < Np; i++) Obj(i) *= FluxObj;
        }	
       
       // New Residual estimation  
       compute_mem_resi();
       
       Iter ++;
       Fonc = SigmaResi / Noise_Ima;
       Delta = (Iter > 1) ? OldFonc - Fonc : Fonc;
       if (Iter > MaxIter) Stop = True;
       else if ((ABS(Delta) < EpsCvg) && (EpsCvg > 0)) Stop = True;
       OldFonc = Fonc; 
       
       if ((Verbose == True) || (Stop == True))
       {
          if (Iter % 2 == 0)
          {   
//              cout.width(8);  // see p343 c++ stroustrup
//              cout.fill(' ');
//              cout.setf(ios::right,ios::adjustfield);
//              cout.setf(ios::scientific,ios::floatfield);
//              cout.precision(4);
    	     if (OptimParam == True)            	 
 	       printf("%d:  Sigma(grad) = %f, Delta = %f, Sigma(resi) = %f, Cvgparam = %f \n", 
	           Iter,  Fonc, Delta, SigmaResi, Cvgparam);
             else printf("%d:  Sigma(grad) = %f, Delta = %f, Sigma(resi) = %f \n", 
	           Iter,  Fonc, Delta, SigmaResi);
          }    
       }
       
       // Get new multiresolution support
       if ((Iter % NbIter_SupportCalc == 0) &&
           ((TypeRegul == DEC_REGUL_SUPPORT) || 
	    (TypeRegul == DEC_REGUL_PROBA_SUPPORT)))
       {
          find_support_obj ();
       }
       
    } while (Stop == False);   

    psf_convol (Obj, Psf_cf, Resi);
    for (i=0; i< Nl*Nc; i++)  Resi(i) = Imag(i) - Resi(i);

    if (GaussConv == True)  convol_gauss(Obj, Fwhm);
    else if (UseICF == True) psf_convol(Obj, Ima_ICF, Obj);
       
    if (Verbose == True)
    {
       INFO(Obj, "Solution: ");
       INFO(Resi, "Resi: ");
    }
}

/*****************************************************************************/

inline float laplacian_val(Ifloat &Obj, int i, int j)
{
  type_border type=I_CONT;
  float Val = - Obj(i,j) + 0.25*( Obj(i+1,j,type)
		    +  Obj(i-1,j,type) +  Obj(i,j+1,type) +  Obj(i,j-1,type));
  return Val;
}		    

/*****************************************************************************/

// static void laplacien_regul(Ifloat &Obj, Ifloat &Grad, float Regul)
// {
//    int Nl = Obj.nl();
//    int Nc = Obj.nc();
//    int i,j;
//    Ifloat Aux(Nl,Nc,"Aux");
//    	
//     for (i = 0; i < Nl; i ++)
//     for (j = 0; j < Nc; j ++) Aux(i,j) =  laplacian_val(Obj,i,j);
//      
//     for (i = 0; i < Nl; i ++)
//     for (j = 0; j < Nc; j ++)
//     {
//         Grad(i,j) -=  Regul*laplacian_val(Aux,i,j);
//     }
// }

/****************************************************************************/
/*
// TRY to build two solutions: one free of noise and the second
// containing the part contminated by the noise
void MWDeconv::mw_deconv(Ifloat *FirstGuess, Ifloat *ICF)
{
    Ifloat Gradient(Nl, Nc, "Info");
    float Func, OldFunc=1e10, Sigma, Delta;
    int i,Iter=0;
    Bool Stop = False;
    float Cvgparam = IterCvg;
      
    // initialization
    // ----------------
    WaveFilterResi = True;
    NbrBand = ModelData->nbr_band();
    Nscale = ModelData->nbr_scale();
    Transform= ModelData->type_trans();
    TabAlpha.alloc(NbrBand);
    // Border = I_CONT;
    init_deconv(FirstGuess,ICF);
    int Np = Nl*Nc;
    NbIter_SupportCalc = DEF_NB_ITER_SUP_CAL;
    if ((RegulParam <= 0.) || (TypeRegul == DEC_NO_REGUL))
    {
       RegulParam = 0.;
       TypeRegul = DEC_NO_REGUL;
    }
  
    Buff.alloc(Nl, Nc, "Buffer");
    Noise_Ima = ModelData->SigmaNoise;
    float Sigma2N = Noise_Ima*Noise_Ima*Np;
    if (NormFlux == True) FluxImag = flux(Imag);
    else FluxImag=0.;
        
    // Residual estimation
    compute_resi();
    // INFO(Resi, "RESI");         
    // start iterations
    // cout << "Start iterate " << endl;
    // INFO(Obj, "obj");
    ObjN.alloc(Nl,Nc,"objn");
    ObjS.alloc(Nl,Nc,"objn");
    cout << "Cvgparam = " << Cvgparam<< endl;
    do
    {
       // compute the gradient
       Gradient = Resi;
       MR_Data.transform (Resi);
       ModelData->threshold(MR_Data);
       MR_Data.recons (Resi);
       Gradient -= Resi;
       //INFO(Resi, "Resi sig");
       //INFO(Gradient, "Resi non sig");
       psf_convol_conj (Resi, Psf_cf);
       //INFO(Resi, "Resi Conv");
       for (i = 0; i < Np; i++) ObjS(i) += Cvgparam *Resi(i);
       //INFO(ObjS, "ObjS");
       
       psf_convol_conj (Gradient, Psf_cf, Resi);
       if (RegulParam > 0) 
       {
          laplacien_regul(ObjN, Resi, RegulParam);
          for (i = 0; i < Np; i++) ObjN(i) += Cvgparam*Resi(i);
          // INFO(ObjN, "ObjN");
          for (i = 0; i < Np; i++) Obj(i) = ObjS(i) + ObjN(i);
       }
       else  Obj = ObjS;
       
       if (PositivConstraint == True) threshold(Obj);
       if (NormFlux == True)
       {
          float FluxObj = flux(Obj);
	  FluxObj = FluxImag / FluxObj;
	  for (i= 0; i < Np; i++) Obj(i) *= FluxObj;
       }

       // next residual
       compute_resi();
       Sigma = sigma(Resi) / Sigma2N;
       Func = Sigma / Noise_Ima;
       Delta = (Iter > 1) ? OldFunc - Func : Func;
       if (Iter > MaxIter) Stop = True;
       else if ((ABS(Delta) < EpsCvg) && (EpsCvg > 0)) Stop = True;
       OldFunc = Func; 
              
       if ((Verbose == True) && ((Iter >= 0) && (Iter % 1 == 0)))
       {
          Func = fonctional() / (float) Np;
          printf("%d: Sigma = %f, Func = %f\n", Iter, Sigma,  Func);
       }
       
       Iter ++;
       if (Iter > MaxIter) Stop = True;
       else if (EpsCvg > FLOAT_EPSILON)
       {       
          if (OldFunc < Sigma) Stop = True;
          if ((ABS(OldFunc - Sigma)) < EpsCvg) Stop = True;
       }
       OldFunc = Sigma;
    } while (Stop == False);
    
    if (Verbose == True)
    {
       INFO(Obj, "Solution: ");
       INFO(Resi, "Resi: ");
    }
}
*/
