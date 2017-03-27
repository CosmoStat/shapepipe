/*******************************************************************************
**
**    UNIT
**
**    Version: 1.0
**
**    Author: Jean-Luc Starck
**
**    Date:  20/12/98 
**    
**    File:  MR_Deconv.cc
**
*******************************************************************************
**
**    DESCRIPTION  Image deconvolution using multiresolution
**    -----------  
**
** 
******************************************************************************/

#include "MR_Deconv.h"

/*****************************************************************************/

void MRDeconv::init_param()
{
    StatNoise = NOISE_GAUSSIAN;
    KeepImagn = False;
    Nl= Nc = 0;
    PsfMaxShift = True;
    GaussConv = False;
    Fwhm = 0.;
    Verbose = False;
    RegulParam = 0;
    MaxIter = DEFAULT_MAX_ITER_DECONV;
    DecMethod = DEFAULT_STD_DECONV;
    EpsCvg = 0.0000;
    IterCvg = 0.1;
    PositivConstraint = True;
    NormFlux = False;
    OptimParam = False;
    KillLastScale = False;
    UseICF = False;
    Noise_Ima = 0.;
    N_Sigma = 3.;
    FluxImag = 0.;
    TypeInit = DEC_INIT_FLAT;
    ModelData = NULL;
    WaveFilterResi = False;
    UseModel = False;
    Border = DEFAULT_BORDER;
    KeepImagn = True;   
    AdjointRec = False;
    MaxMRCleanIter = DEFAULT_CLEAN_NITER;
    SupportDilate = False;
    UseMRCEnergy = True;
    CleanLastScale=True;
    CleanFirstScale=0;
    TabCB_cf=NULL;
}

/*****************************************************************************/

void MRDeconv::convol_gauss(Ifloat &Data, float FWHM)
{
   Ifloat Gauss;
   Gauss = im_gaussian(Data.nl(), Data.nc(), FWHM);
   norm_flux(Gauss);
   psf_convol (Data, Gauss);
}

/*****************************************************************************/

void MRDeconv::init_deconv(Ifloat *FirstGuess, Ifloat *ICF)
{
   float Flux;
   Nl = Imag.nl();
   Nc = Imag.nc();
   Bool NonAddMethod = False;
   int i;
     
   // set NonAddMethod to True for LUCY and MAP
   if ((DecMethod == DEC_LUCY) || (DecMethod == DEC_MR_LUCY) 
       || (DecMethod == DEC_MAP) || (DecMethod == DEC_MR_MAP))
	                                         NonAddMethod = True;     

   if (FirstGuess != NULL) TypeInit = DEC_INIT_GUESS;
   else Obj.alloc(Nl, Nc, "object");
   
   // Multiplyer image allocation MEM in memthod
   if ((DecMethod == DEC_MEM) || (DecMethod == DEC_MEM_MODEL)) 
                                              Mult.alloc(Nl,Nc,"mult");

   // if an ICF is introduced, it can be either a Gaussian
   // or an image given by the user
   // the PSF must be convolved by the ICF
   if (GaussConv == True) convol_gauss(Psf, Fwhm);
   else if (ICF != NULL)
   {
      Ima_ICF.alloc(Psf_cf.nl(), Psf_cf.nc(), "ICF");
      dec_center_psf (*ICF, Ima_ICF);
      UseICF = True;
      FFT2D.convolve(Psf, Ima_ICF, Psf);
      // fft2d_conv(Psf, Ima_ICF, Psf);
   }
    
   // centering the PSF (normally the max of PSF is at the center of the image)
   // the size of Psf_cf is a power of 2, which is not necessary the size of
   // Psf, or the image size
   psf_get(Psf, Psf_cf, Nl, Nc, PsfMaxShift); 
            
   switch (TypeInit)
   {
      case DEC_INIT_ZERO: Obj.init();  break;
      case DEC_INIT_FLAT: 
             Flux = flux(Imag) / (float) (Nl*Nc);
  	     Obj.init(Flux);
	     break;
      case DEC_INIT_IMA:
             Obj = Imag;
	     break;
      case DEC_INIT_GUESS:
             Obj = *FirstGuess;
             break;
   }
        
   Resi.alloc(Nl,Nc,"resi");
// if ((DecMethod == DEC_LUCY) || 
//      ((OptimParam == True) && (StatNoise == NOISE_POISSON))) 
//         KeepImagn = True;
   if (KeepImagn == True) Imag_n.alloc(Nl,Nc,"Imag_n");

   /* Noise estimation in the Data */
   if (StatNoise == NOISE_GAUSSIAN)
   {
        if (Noise_Ima < FLOAT_EPSILON)  
	{
	   Noise_Ima = detect_noise_from_med (Imag);
	   if (Verbose == True)
	      cout << "Sigma Noise = " << Noise_Ima << endl;
	}
   }
   if (Noise_Ima < FLOAT_EPSILON) Noise_Ima = 1.;
   
   if ((WaveFilterResi == True) || (KillLastScale == True))
   {
      WaveFilterResi = True;
      // MR_Data.alloc(Nl, Nc, ModelData->nbr_scale(), 
      //              ModelData->type_trans(), "MR_Data");
      MR_Data.alloc(Nl,Nc, ModelData->nbr_scale(), ModelData->type_trans(), 
                  ModelData->filter_bank(),  ModelData->TypeNorm,  ModelData->nbr_undec_scale());

      MR_Data.Border = Border;
      ModelData->model(Imag, MR_Data);
      MR_Data.ModifiedATWT = True;
 
      if (KillLastScale == True) 
      {
         int Nb = ModelData->nbr_band()-1;
         MR_Data.transform(Obj);
         if (TypeInit == DEC_INIT_IMA) ModelData->threshold(MR_Data);
         MR_Data.band(Nb).init();
         MR_Data.recons(Obj);
      }
   }
   
   // Zero value in the first guess leads to problem with
   // Lucy and MAP methods
   if (NonAddMethod == True)
   {
      for (i=0; i < Nl*Nc; i++)
      {
 	  if (Obj(i) < FLOAT_EPSILON) Obj(i) = FLOAT_EPSILON;
	  if (Imag(i) < FLOAT_EPSILON) Imag(i) = FLOAT_EPSILON;
      }
   }
   
   switch (DecMethod)
   {
      case DEC_MARKOV_LUCY:
      case DEC_MARKOV:
           RIm.GradOperType = OPER_MARKOV_8;
	   MRIm.NbrScale = ModelData->nbr_scale();
           break;
      case DEC_TIKHONOV:
           RIm.GradOperType = OPER_LAPLACIAN;
           break;
      case DEC_MR_MEM:
      case DEC_MR_MEM_NOISE:
      case DEC_MEM: 
      case DEC_MEM_MODEL:
      case DEC_MAP: 
      case DEC_CITTERT: 
      case DEC_GRADIENT: 
      case DEC_INVERSE: 
      case DEC_LUCY: 
      case DEC_CLEAN:
      case DEC_MR_CITTERT: 
      case DEC_MR_GRADIENT:
      case DEC_MR_LUCY:
      case DEC_MR_MAP:
      case DEC_MR_CLEAN:
      case DEC_MR_VAGUELET:
      case DEC_MR_INVERSE:
           break;
   }
}

/****************************************************************************/

void MRDeconv::compute_resi()
{
   int i;
   psf_convol (Obj, Psf_cf, Resi);
   if (KeepImagn == True) Imag_n = Resi;
   for (i=0; i< Nl*Nc; i++)  Resi(i) = Imag(i) - Resi(i);
}

/****************************************************************************/

void MRDeconv::im_grad(Ifloat &Gradient)
{
   // int Np = Nl*Nc;
   int i,j;
   float Scale;
   double Logmin = 1e-9;
   
   
     // Regularization by residual filtering
    // INFO(Resi, "RESI");
 
    if ((WaveFilterResi == True) || (KillLastScale == True))
                                            MR_Data.transform (Resi);
    if (WaveFilterResi == True) ModelData->threshold(MR_Data);
  
    if (KillLastScale == True)
    {
       if (AdjointRec == False) 
       {
	     MR_Data.band(ModelData->nbr_band()-1).init();
	     MR_Data.recons (Resi);
       }
       else MR_Data.rec_adjoint (Resi, False);
    }
    else if (WaveFilterResi == True)
    {
       if (AdjointRec == False) MR_Data.recons (Resi);
       else MR_Data.rec_adjoint (Resi);
    }
    
   switch (DecMethod)
   {
            case DEC_MEM:
  	       Scale = flux(Imag_n) / FluxImag;
	       for (i = 0; i < Nl; i++)
	       for (j = 0; j < Nc; j++)
	       {
		  // Mult(i) += Imag(i)*Scale - Imag_n(i);
                  float ValNum = Imag(i,j) * Scale;
		  if (ValNum < Logmin) ValNum = Logmin;
		  float ValDen = Imag_n(i,j);
		  if (ValDen < Logmin) ValDen = Logmin;
		  Mult(i,j) += log(ValNum/ValDen);
	       }
               Gradient = Mult; 		  
  	       psf_convol(Gradient, Psf_cf);
	       for (i = 0; i < Nl; i++)
	       for (j = 0; j < Nc; j++)
	            Gradient(i,j) = exp(Gradient(i,j));
	       break;
	   case DEC_MEM_MODEL:
	       Gradient = Resi;
               psf_convol_conj (Gradient, Psf_cf);
	       if (UseModel == True)
	       {
	          for (i = 0; i < Nl; i++) 
		  for (j = 0; j < Nc; j++) 
	            if (Obj(i,j) > MemModel(i,j)) 
		       Gradient(i,j) -= RegulParam*log(Obj(i,j)/MemModel(i,j));
	       }
	       else 
	       {
	          for (i = 0; i < Nl; i++)
		  for (j = 0; j < Nc; j++)
		    if (Obj(i,j) > FLOAT_EPSILON) 
		       Gradient(i,j) -= RegulParam*(1+log(Obj(i,j)));
	       }
	       break;
	   case DEC_MR_MAP:   
 	   case DEC_MAP:
	        for (i = 0; i < Nl; i++)  
		for (j = 0; j < Nc; j++)
		   if (Imag_n(i,j) > FLOAT_EPSILON)  
		        Gradient(i,j) = Resi(i,j) / Imag_n(i,j);
		   else Gradient(i,j) =  FLOAT_EPSILON;
 		psf_convol_conj (Gradient, Psf_cf);
		for (i = 0; i < Nl; i++)
		for (j = 0; j < Nc; j++) Gradient(i,j) = exp( Gradient(i,j));
 		break;
	   case DEC_MR_LUCY:
           case DEC_LUCY:
	   case DEC_MARKOV_LUCY:		
                for (i = 0; i < Nl; i++)
		for (j = 0; j < Nc; j++) 
                {
                   if (Imag_n(i,j) > FLOAT_EPSILON)  Gradient(i,j) = Resi(i,j) / Imag_n(i,j);
                   else  Gradient(i,j) = 0.;
                }
 		psf_convol_conj (Gradient, Psf_cf);
                Gradient *= Obj;
		if (RegulParam > 0) MRIm.obj_regul(Obj, Gradient, RegulParam*Noise_Ima);
                break;
	   case DEC_MARKOV:
           case DEC_MR_GRADIENT:
	   case DEC_GRADIENT:
 		Gradient = Resi;
 		psf_convol_conj (Gradient, Psf_cf);
		if (RegulParam > 0) MRIm.obj_regul(Obj, Gradient, RegulParam*Noise_Ima);
                break;
	   case DEC_TIKHONOV:
 		Gradient = Resi;
 		psf_convol_conj (Gradient, Psf_cf);
		if (RegulParam > 0) RIm.obj_regul(Obj, Gradient, RegulParam);
                break;
	   case DEC_MR_CLEAN:
	   case DEC_MR_CITTERT:	
	   case DEC_CITTERT:
	         Gradient = Resi;
		 if (RegulParam > 0) MRIm.obj_regul(Obj, Gradient, RegulParam*Noise_Ima);
                 break;
           default:
                cerr << "mr_deconv: Not implemented in this procedure ... ";
                cerr << endl;
                exit (0);
                break;
   }
}

/****************************************************************************/

float MRDeconv::find_optim_poisson(Ifloat &Gradient)
{
  float Cvgparam = 1., OldCvgparam;
  Ifloat Delta(Nl,Nc,"Buff");
  int i,Iter=0;
  int Np = Nl*Nc;
  float MinVal=1., MaxVal=100.;
  float Val, Deriv1=0., Deriv2=0.;
  if (RegulParam > 0.5) MinVal = 1. / (2.*RegulParam);
  
  // sear the maximum value
  for (i=0; i < Np; i++) 
    if (Gradient(i) < - FLOAT_EPSILON)
    {
       Val = - Obj(i) / Gradient(i);
       if (Val < MaxVal) MaxVal = Val;
    }
  // cout << "interval: Min = " << MinVal << " Max = " << MaxVal << endl;
  
  if (MaxVal > 0.)
  {
     // dec_convol (Gradient, Psf_cf, Delta);
     psf_convol(Gradient, Psf_cf, Delta);
     // Cvgparam = (MinVal + MaxVal) / 2.;
     do
     {
        Iter++;
        OldCvgparam = Cvgparam;
        // Calculate first and second derivative
        for (i=0; i < Np; i++)
	{
	   // First derivative
	   Val = Imag_n(i) + Cvgparam*Delta(i);
	   if (ABS(Val) > FLOAT_EPSILON)
	         Deriv1 += Imag(i)*Delta(i) / Val ;

           // Second derivative 
	   Val *= Val;
	   if (ABS(Val) > FLOAT_EPSILON)
	          Deriv2 -= Imag(i)*Delta(i)*Delta(i) / Val;
	}	  
        Cvgparam -= Deriv1/ Deriv2;
	
	// cout << "Iter " << Iter << " Cvgparam = " << Cvgparam << endl;
     } while ((Iter < 10) && ((ABS(Cvgparam-OldCvgparam) > 0.1*OldCvgparam)));
     if (Cvgparam > MaxVal) Cvgparam = MaxVal;
  }
  if (Cvgparam < MinVal) Cvgparam = MinVal;
  return Cvgparam;
}
/****************************************************************************/

float MRDeconv::find_optim_tikhonov(Ifloat &Gradient)
{
   float Cvgparam = IterCvg;
   Ifloat Buff(Nl,Nc,"Buff");
   float LapG,LapO, Num, Den;
   int i,j;
   
   if (RegulParam > 0.)
   {
      // dec_convol (Gradient, Psf_cf, Buff);
      psf_convol(Gradient, Psf_cf, Buff);
      Num = flux(Buff*Resi);
      Den = energy(Buff);
      for (i = 0; i < Nl; i ++)
      for (j = 0; j < Nc; j ++)
      {
         LapG = RIm.laplacian_val(Gradient,i,j);
	 LapO = RIm.laplacian_val(Obj,i,j);
	 Num -= RegulParam*LapG*LapO;
         Den += RegulParam*LapG*LapG;
      }
      Cvgparam =  Num / Den;
   }
   else Cvgparam = find_optim_xi2(Gradient);
   
   return Cvgparam;  
}

/****************************************************************************/

float MRDeconv::find_optim_xi2(Ifloat &Gradient)
{
   float Cvgparam = IterCvg;
   Ifloat Buff(Nl,Nc,"Buff");
   float MinVal = 1.;
   if (RegulParam > 0.5) MinVal = 1. / (2.*RegulParam);
   // dec_convol (Gradient, Psf_cf, Buff);
   psf_convol(Gradient, Psf_cf, Buff);
   Cvgparam = flux(Buff*Resi)/energy(Buff);
   if (Cvgparam < MinVal) Cvgparam = MinVal;
   else if (Cvgparam > 10) Cvgparam = 10.;
   return Cvgparam;  
}

/****************************************************************************/

float MRDeconv::im_find_optim(Ifloat &Gradient)
{
   float Cvgparam = 1.;
   if (DecMethod == DEC_TIKHONOV) Cvgparam = find_optim_tikhonov(Gradient);
   else if (StatNoise == NOISE_GAUSSIAN)  Cvgparam = find_optim_xi2(Gradient);
   else if (StatNoise == NOISE_POISSON) Cvgparam = find_optim_poisson(Gradient);
   if (Verbose == True)
   {
      cout << "Optim Val = " << Cvgparam << endl;
   }
   
   if (Cvgparam < 0) Cvgparam = 0.;
   return Cvgparam;  
}

/****************************************************************************/

float MRDeconv::fonctional()
{
   float Val,Func = 0.;
   int i,j;
     
   // Maxumum likehood part of the functional
   if (StatNoise == NOISE_GAUSSIAN) Func = energy(Resi); 
   else  // POISSON noise
   {
      Func = energy(Imag);  // we add a constant in order to have 
                            // a positive functional
      for (i = 0; i < Nl; i ++)
      for (j = 0; j < Nc; j ++)
        if (Imag_n(i,j) > 0)  Func -= Imag(i,j)*log(Imag_n(i,j));
   }
   
   // regularized part of the functional
   if (DecMethod == DEC_TIKHONOV)  // case of Tikhonov regularization
   { 	
      for (i = 0; i < Nl; i ++)
      for (j = 0; j < Nc; j ++) 
      {
           Val = RIm.laplacian_val(Obj,i,j);
           Func += Val*Val*RegulParam;
      }
   }
   else if ((DecMethod == DEC_MEM) || 
         ((DecMethod == DEC_MEM_MODEL) && (UseModel == False))) // case of MEM regularization
   { 	
      for (i = 0; i < Nl; i ++) 
      for (j = 0; j < Nc; j ++)
          if (Obj(i,j) > 0)  Func += RegulParam*Obj(i,j)*log(Obj(i,j));
   }
   else if (DecMethod == DEC_MEM_MODEL)
   { 	
      for (i = 0; i < Nl; i ++) 
      for (j = 0; j < Nc; j ++)
               if (Obj(i,j) > 0)  
	        Func +=  RegulParam*((-Obj(i,j) + MemModel(i,j)) + Obj(i,j)*log(Obj(i,j)/MemModel(i,j)));
   }
   return Func;
}

/****************************************************************************/

void MRDeconv::im_iter_deconv()
{
    Ifloat Gradient(Nl, Nc, "Info");
    float Func, OldFunc=1e10, Sigma;
    int i,j,Iter=0;
    Bool Stop = False;
    float Cvgparam = IterCvg;
    float Sigma2N = Noise_Ima*Noise_Ima*Nl*Nc;
    if (NormFlux == True) FluxImag = flux(Imag);
    else FluxImag=0.;

    // Residual estimation
    compute_resi();
    // INFO(Resi, "RESI");         
    // start iterations
    if ((Verbose == True) && (RegulParam > 0))
           cout << "Start iterate: Regul =  " <<  RegulParam << endl;
    do
    {         
       // compute the gradient
       im_grad(Gradient);
  
       if (OptimParam == True) Cvgparam = im_find_optim(Gradient);
        
       // calculate the next object estimation with positivity  constraint
       if ((DecMethod == DEC_MAP) || (DecMethod == DEC_MR_MAP))
       {
          for (i = 0; i < Nl; i++)
	  for (j = 0; j < Nc; j++) Obj(i,j) *= Gradient(i,j);
       }
       else if (DecMethod == DEC_MEM) Obj = Gradient;
       else 
       {
          for (i = 0; i < Nl; i++) 
	  for (j = 0; j < Nc; j++) Obj(i,j) += Cvgparam * Gradient(i,j);
       }
       if (PositivConstraint == True) threshold(Obj);
       if ((NormFlux == True) && (KillLastScale == False))
       {
          float FluxObj = flux(Obj);
	  FluxObj = FluxImag / FluxObj;
	  for (i= 0; i < Nl; i++)
	  for (j= 0; j < Nc; j++) Obj(i,j) *= FluxObj;
       }

       // next residual
       compute_resi();
       Sigma = energy(Resi) / Sigma2N;
       
       if ((Verbose == True) && ((Iter >= 0) && (Iter % 1 == 0)))
       {
          Func = fonctional() / (float) (Nl*Nc);
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

/****************************************************************************/

void MRDeconv::init_mrc()
{
   Nl = Imag.nl();
   Nc = Imag.nc();
   Ifloat D_Beam(Nl, Nc, "PSF");

   TypeInit = DEC_INIT_ZERO;
   AdjointRec=False;
   Obj.alloc(Nl, Nc, "object");
   Resi.alloc(Nl,Nc,"resi");
   if (KeepImagn == True) Imag_n.alloc(Nl,Nc,"Imag_n");

   /* Noise estimation in the Data */
   if (StatNoise == NOISE_GAUSSIAN)
   {
        if (Noise_Ima < FLOAT_EPSILON)  
	{
	   Noise_Ima = detect_noise_from_med (Imag);
	   if (Verbose == True)
	      cout << "Sigma Noise = " << Noise_Ima << endl;
	}
   }
   if (Noise_Ima < FLOAT_EPSILON) Noise_Ima = 1.;
   
   WaveFilterResi = True;
   // MR_Data.alloc(Nl, Nc, ModelData->nbr_scale(), 
   //                 ModelData->type_trans(), "MR_Data");
   MR_Data.alloc(Nl,Nc, ModelData->nbr_scale(),ModelData->type_trans(), 
                 ModelData->filter_bank(), ModelData->TypeNorm, 
                 ModelData->nbr_undec_scale());
   MR_Data.Border = Border;
   ModelData->model(Imag, MR_Data);
   dec_center_psf (Psf, D_Beam);
   norm_flux(D_Beam);
   Psf.resize(Nl,Nc);
   Psf = D_Beam;
}

/*****************************************************************************/

void MRDeconv::mrc_rec_obj(MultiResol &MR_SOL, Icomplex_f *TabCB_cf)
{
   int b;
   int Nb = (KillLastScale == False) ? MR_Data.nbr_band(): MR_Data.nbr_band()-1;

   for (b=0; b < CleanFirstScale; b++) MR_Data.band(b).init();
   for (b=CleanFirstScale; b < Nb; b++)
   {  
      psf_convol (MR_SOL.band(b), TabCB_cf[b], MR_Data.band(b), False);
   }
   if (KillLastScale == True) MR_Data.band(Nb).init();
   if (AdjointRec == False) MR_Data.recons (Obj);
   else MR_Data.rec_adjoint (Obj);
   threshold(Obj);
}

/*****************************************************************************/

void MRDeconv::mrc(Ifloat *ICF)
{
   int Iter=0;
   int b,i,j;
   int Kaver = 1;
   float DetectLevel,GammaClean = RegulParam;
   extern int DecVerbose;
   float OldFunc=1e10, Sigma;
   Bool Stop = False;
   float Sigma2N = Noise_Ima*Noise_Ima;
   int Last = MR_Data.nbr_band()-1;

   DecVerbose = Verbose;
   MultiResol MR_SOL(Nl,Nc,ModelData->nbr_scale(), 
                         ModelData->type_trans(), "MR_SOL");
   MultiResol MR_PSF(Psf.nl(),Psf.nc(),ModelData->nbr_scale(), 
                         ModelData->type_trans(), "MR_PSF");
   MR_PSF.transform(Psf);
   //norm_flux(MR_PSF.band(MR_Data.nbr_band()-1));
   MR_Data.transform(Imag);
   TabCB_cf = new Icomplex_f [MR_Data.nbr_band()];
   int Nb = (CleanLastScale == True) ? MR_Data.nbr_band(): MR_Data.nbr_band()-1;
   // int Nb = MR_Data.nbr_band()-1;
   // Apply CLEAN multiresolution
   for (b=0; b < Nb; b++)
   {
      // initialization
      DetectLevel = ModelData->sigma(b,0,0)*ModelData->NSigma[b];
      if (b == MR_Data.nbr_band()-1)
      {
         Kaver = 4;
	 DetectLevel = 0.;
      }
      if (Verbose == True) cout << "CLEAN Band " << b << endl;
      
      // run CLEAN
      Bool AddResi=False;
      Bool SearchPositiv=True;
      if (b >= CleanFirstScale)
             dec_clean (MR_Data.band(b), MR_PSF.band(b), MR_SOL.band(b), 
                  0.,  DetectLevel, GammaClean, MaxMRCleanIter, Kaver, AddResi,
		  SearchPositiv,UseMRCEnergy);

      int Cpt=0;
      for (i=0; i < MR_Data.size_band_nl(b); i++)
      for (j=0; j < MR_Data.size_band_nc(b); j++)
      {
	  if (MR_SOL(b,i,j) != 0) 
	  {
	     Cpt++;
	     if (b != Last)  ModelData->support(b,i,j)= VAL_SupOK;
	  }
	  else if (b != Last) ModelData->support(b,i,j)= VAL_SupKill; 
      }
      if (Verbose == True) cout << "  number of peaks = " << Cpt << endl;
      if ((b != Last) && (SupportDilate == True)) ModelData->dilate_support(b);
   }
   // MR_Data.write("xx_resi.mr");
   
   // Create the wavelet CLEAN BEAM
   Resi.init();
   Resi(Nl/2,Nc/2) = 1.;
   MR_PSF.transform(Resi);
   float *TabMax = new float [MR_Data.nbr_band()];
   for (b=0; b < MR_Data.nbr_band(); b++)
   {
      if (MR_PSF.Set_Transform == TRANSF_PYR)  TabMax[b] = pow((double) 4.,(double) b);
      else TabMax[b] = 1;
      TabCB_cf[b].alloc(MR_PSF.size_band_nl(b), MR_PSF.size_band_nc(b), "band"); 
      fft2d(MR_PSF.band(b), TabCB_cf[b]);
    }
    
   // Reconstruct the solution
   mrc_rec_obj(MR_SOL, TabCB_cf);
//   MR_PSF.write("xx_psf.mr");
//   MR_Data.write("xx_sol0.mr");
               
   // Create the Multiresolution support
//    for (b=0; b < MR_Data.nbr_band()-1; b++)
//    for (i=0; i < MR_Data.size_band_nl(b); i++)
//    for (j=0; j < MR_Data.size_band_nc(b); j++)
//      if (ModelData->support(b,i,j) < 9) MR_Data(b,i,j) = 1;
//      else MR_Data(b,i,j) = 0;
//    MR_Data.write("xx_sol.mr");
//    io_write_ima_float("xx_sol.fits", Obj);
   
    // Iterative refinment of the solution
    if (GaussConv == True) 
    {
       GaussConv = True;
       convol_gauss(Psf, Fwhm);
    }
    else if (ICF != NULL)
    {
      Ima_ICF.alloc(Psf_cf.nl(), Psf_cf.nc(), "ICF");
      dec_center_psf (*ICF, Ima_ICF);
      UseICF = True;
      fft2d_conv(Psf, Ima_ICF, Psf);
    }
    psf_get(Psf, Psf_cf, Nl, Nc, PsfMaxShift); 

    compute_resi();
    do
    {       
       // WT of the residual
       MR_Data.transform(Resi);
       
       // Update the dirac list
       for (b=CleanFirstScale; b < Nb; b++)
       {
           for (i=0; i < MR_Data.size_band_nl(b); i++)
           for (j=0; j < MR_Data.size_band_nc(b); j++)
           if (((b !=  MR_Data.nbr_band()-1) && ((*ModelData)(b,i,j) == True))
	      || ((b == MR_Data.nbr_band()-1) && (MR_SOL(b,i,j) != 0)))
	   {
	     MR_SOL(b,i,j) += MR_Data(b,i,j) / TabMax[b];
	     if (MR_SOL(b,i,j) < FLOAT_EPSILON) MR_SOL(b,i,j) = FLOAT_EPSILON;
	   }
        }
       if ((CleanLastScale == False) && (KillLastScale == False))
                                         MR_SOL.band(Nb) += MR_Data.band(Nb);
       // New solution
       mrc_rec_obj(MR_SOL, TabCB_cf);
       
       // New residual
       compute_resi();
       
       Sigma = sigma(Resi);
       Sigma *= Sigma / Sigma2N;
       if ((Verbose == True) && ((Iter >= 0) && (Iter % 1 == 0)))
       {
          printf("%d: Sigma = %f\n", Iter+1, Sigma);
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
   //MR_SOL.write("xx_sol1.mr");
   //MR_Data.write("xx_sol2.mr");  
   delete [] TabCB_cf;
}

/****************************************************************************/

void dec_std_inverse (Ifloat &Imag, Icomplex_f &Psf_cf, Ifloat &Result)
{
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
   fft2d (O_cf);
   // O_cf /= Psf_cf;
    for (i = 0; i < Nl1; i++) 
    for (j = 0; j < Nc1; j++) 
    {
        if (Psf_cf(i,j).real() > 0.01)
                        O_cf(i,j) /= Psf_cf(i,j).real();
        else O_cf(i,j) = complex_f(0,0);
    }

   fft2d (O_cf, -1);
   if ((Nl != Nl1) || (Nc != Nc1)) im_extract (O_cf, Result, FluxPsf);
   else for (i = 0; i < Nl; i++)
        for (j = 0; j < Nc; j++)  Result(i,j) = O_cf(i,j).real() / FluxPsf;
}

/************************************************************************/

void MRDeconv::va_inverse_data(Ifloat &NoiseSimu)
{
    if (StatNoise == NOISE_GAUSSIAN) im_noise_gaussian (NoiseSimu, Noise_Ima);
   else
   {
      // Find the noise map by filtering the data
      type_filter Filter = (MR_Data.Set_Transform == TRANSF_MALLAT) ?
                           FILTER_THRESHOLD: FILTER_ITER_THRESHOLD;
      MRFiltering CFilter(*ModelData, Filter);
      CFilter.Verbose = Verbose;
      CFilter.filter(Imag, Obj);
   
      for (int i=0; i < Nl; i++)
      for (int j=0; j < Nc; j++) NoiseSimu(i,j) = Imag(i,j) - Obj(i,j);
   }
    
   // inverse deconvolution
   // dec_std_inverse (Imag, Psf_cf, Obj);
   // dec_std_inverse (NoiseSimu, Psf_cf, NoiseSimu);
   dec_inverse(Imag, Psf_cf, Obj);
   dec_inverse(NoiseSimu, Psf_cf, NoiseSimu);
}

/************************************************************************/

void MRDeconv::vaguelet()
{
   int b,i,j;
   int Nb = MR_Data.nbr_band()-1;
   Ifloat NoiseSimu(Nl,Nc,"simu");
   // cout << "NoiseSimu = " << "Nl = " << Nl << " Nc = " << Nc << " Nelem = " << NoiseSimu.n_elem() << endl;
   // inverse deconvolution on both imag and noise map
   va_inverse_data(NoiseSimu);
   // io_write_ima_float("xx_inv.fits", Obj);
   // io_write_ima_float("xx_noise.fits", NoiseSimu);
   
   // fltarray TabSigma(Nb);
   //MR_Data.transform(NoiseSimu);
   //for (b=0; b < Nb; b++) TabSigma(b) = sigma(MR_Data.band(b));
   //for (b=0; b < Nb; b++) cout << "Band " << b << " sigma = " << TabSigma(b)*ModelData->NSigma[b]  << endl;
   
    // Find the threshold considering a correlated noise map
    // INFO(NoiseSimu, "NN");
    StatNoiseMap STM;
    STM.alloc(NoiseSimu, MR_Data.nbr_scale(), MR_Data.Type_Transform);
    fltarray Tab_Proba(Nb+1);
    for (b=0; b < Nb; b++) 
     Tab_Proba(b) =  (1. - erf((double) ModelData->NSigma[b] / sqrt((double) 2.)));   
    fltarray ThresholdMin(Nb+1);
    fltarray ThresholdMax(Nb+1);
    STM.find_threshold(MR_Data.nbr_band(), Tab_Proba, ThresholdMin, ThresholdMax);
    //if (Verbose == True)
    //  for (b=0; b < Nb; b++) 
    //   cout << "Band " << b << " Eps = " << Tab_Proba(b) << " Level Min = " << ThresholdMin(b) <<  " Level Max = " << ThresholdMax(b) << endl;
 
    // Update the multiresolution support
    MR_Data.transform(Obj);   
    for (b=0; b < Nb; b++) 
    for (i=0; i < MR_Data.size_band_nl(b); i++)
    for (j=0; j < MR_Data.size_band_nc(b); j++)
    {
        if (ModelData->OnlyPositivDetect == True)
	{
	  if (MR_Data(b,i,j) > ThresholdMax(b)) ModelData->support(b,i,j) = 1;
        }
        else 
	{
	   if ((MR_Data(b,i,j) < ThresholdMin(b)) || (MR_Data(b,i,j) > ThresholdMax(b))) 
                                                ModelData->support(b,i,j) = 1;
	   //if (ABS(MR_Data(b,i,j)) < ModelData->NSigma[b]*TabSigma(b)) 
	   //ModelData->support(b,i,j) = 0;
	   //else  ModelData->support(b,i,j) = 1;
	}
   }
    
    // Iterative filtering of the solution
    if (MR_Data.Set_Transform != TRANSF_MALLAT)
    {
       type_filter Filter = (MR_Data.Set_Transform == TRANSF_MALLAT) ?
                           FILTER_THRESHOLD: FILTER_ITER_THRESHOLD;
       MRFiltering CFilter(*ModelData, Filter);
       NoiseSimu=Obj;
       CFilter.Max_Iter = MaxIter;
       CFilter.Epsilon = EpsCvg;
       CFilter.Sup_Set = False;
       CFilter.Verbose = Verbose;
       CFilter.KillLastScale = KillLastScale;
       CFilter.mr_support_filter(NoiseSimu, MR_Data, Obj);
    }
    else
    {
        fltarray TabLevel(Nb+1);
        for (b=0; b < Nb; b++)  
           TabLevel(b) = MAX(ABS(ThresholdMin(b)), ThresholdMax(b));
        ModelData->threshold(MR_Data);
	if (MaxIter <= 1) MR_Data.recons(Obj);
        else mr_ortho_regul_rec(MR_Data, Obj, TabLevel.buffer(), MaxIter);
    }
    if (PositivConstraint == True) threshold(Obj);
}
 
/****************************************************************************/

void MRDeconv::im_deconv(Ifloat *FirstGuess, Ifloat *ICF)
{
   float Mean,Sigma;
   float GammaClean = RegulParam;
   Bool IterMethod = True;
   
   switch (DecMethod)
   {
      case DEC_CLEAN:
            IterMethod = False;
            StatNoise = NOISE_GAUSSIAN;
 	    Nl = Imag.nl();
	    Nc = Imag.nc();
	    Obj.alloc(Nl,Nc, "Obj");
	    Resi.alloc(Nl,Nc, "Resi");
            Resi = Imag;
            norm_flux (Psf, 1.);
            if (Noise_Ima < FLOAT_EPSILON)  
                            Noise_Ima = detect_noise_from_med (Resi);
            sigma_clip(Resi,Mean,Sigma);
	    cout << "Mean+N_Sigma*Noise_Ima = " << Mean+N_Sigma*Noise_Ima << endl;
	    cout << "GammaClean = " << GammaClean << endl;
	    if (GaussConv == True) cout << "GaussConv = " << Fwhm << endl;
	    
            dec_clean (Resi, Psf, Obj, 0., Mean+N_Sigma*Noise_Ima, 
                       GammaClean, MaxIter, 
                       DEFAULT_CLEAN_KAVER, DEFAULT_CLEAN_ADDRESI);
            if (ICF != NULL)
            {
               Ima_ICF.alloc(Nl,Nc, "Psf");
               dec_center_psf (*ICF, Ima_ICF); 
	       UseICF = True;           
	    }   
	    break;
       case DEC_INVERSE:
            IterMethod = False;
            StatNoise = NOISE_GAUSSIAN;
            init_deconv(FirstGuess,  ICF);
	    dec_inverse (Imag, Psf_cf, Obj);
	    compute_resi();
            break;      
      case DEC_MAP:
             TypeInit = DEC_INIT_FLAT;
             OptimParam = False;
	     StatNoise = NOISE_POISSON;
 	     break;
      case DEC_MARKOV_LUCY:
      case DEC_LUCY:
             TypeInit = DEC_INIT_FLAT;
             StatNoise = NOISE_POISSON;
 	     break;	
      case DEC_CITTERT:
      case DEC_MEM:
      case DEC_MEM_MODEL:
             TypeInit = DEC_INIT_FLAT;
             NormFlux = True;
             OptimParam = False;
	     StatNoise = NOISE_GAUSSIAN;	
 	     break;
      case DEC_MARKOV:
      case DEC_GRADIENT:
      case DEC_TIKHONOV:
             StatNoise = NOISE_GAUSSIAN;	
 	     break;
      case DEC_MR_CITTERT:
      case DEC_MR_GRADIENT:
             WaveFilterResi = True; 
             break;
      case DEC_MR_LUCY:
      case DEC_MR_MAP:
             TypeInit = DEC_INIT_FLAT;
             WaveFilterResi = True;
  	    break;
      case DEC_MR_CLEAN:
            IterMethod = False;
            if (FirstGuess != NULL)  
            {
              cerr << "Error: a first guess cannot be used by Multiresolution CLEAN ... " << endl;
              exit(-1);
            }            
            init_mrc();
            mrc(ICF);
	    break;  
      case DEC_MR_VAGUELET:
            IterMethod = False;
	    WaveFilterResi = True;
            init_deconv(FirstGuess,  ICF);
	    vaguelet();
	    compute_resi();
            break;
      case DEC_MR_MEM:
      case DEC_MR_MEM_NOISE:
      case DEC_MR_INVERSE:
      default:
                cerr << "mr_deconv: Not implemented in this procedure ... ";
                cerr << endl;
                exit (0);
                break;
   }
   if (IterMethod == True)
   {
      init_deconv(FirstGuess,ICF);
      im_iter_deconv();
   }
   
   if (GaussConv == True)  convol_gauss(Obj, Fwhm);
   else if (UseICF == True) psf_convol(Obj, Ima_ICF, Obj);
}
