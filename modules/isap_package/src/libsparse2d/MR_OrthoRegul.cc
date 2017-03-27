/*******************************************************************************
**
**    UNIT
**
**    Version: 1.0
**
**    Author: Jean-Luc Starck
**
**    Date:  97/09/03 
**    
**    File:  MR_OrthoRegul.cc
**
*******************************************************************************
**
**    DESCRIPTION  Modify the wavelet coefficient of an orthogonal
**    -----------  wavelet transform (type MALLAT) in order to sastisfy
**                 a regularization constraint. 
**
**
*******************************************************************************/

#include "IM_Obj.h"
#include "MR_Obj.h"
#include "IM_IO.h"
// #include "MR_NoiseModel.h"
// #include "MR_Filter.h"
// #include "MR_Sigma.h"

/****************************************************************************/

// static void mr_ortho_regul_ima_rec(MultiResol &MR_Step, Ifloat &Result, float *Quant, int MaxIter)
// {
//    int Nl = Result.nl();
//    int Nc = Result.nc();
//    int w,Iter=0;
//    Ifloat Lap(Nl,Nc, "Dat");
//    MultiResol MR_Iter (Nl, Nc, MR_Step.nbr_scale(), MR_Step.Type_Transform, "MR_Iter");
//    int Nw = MR_Step.nbr_mr_coeff();
//    unsigned char *TabSup = new unsigned char [Nw];
//    float F;
//    float F_old=0; 
//    float Delta_F=10;;
// 
//    type_border Bord=I_CONT;
//    float LevelT=INFINITY;
//    float *TabSave = new float [Nw];
//    float Delta_w,LevelQ;
//    details Det;
//    int Scale,Indi,Indj;
//    
//    for (w=0; w < MR_Step.nbr_mr_coeff(); w++)
//      {
//        // Save initial value of dequantized coefficients for range constraint
//        TabSave[w]=MR_Step(w);
//        
//        // compute the multiresolution support
//        TabSup[w] = (ABS(MR_Step(w)) > FLOAT_EPSILON) ? 1 : 0; 
//        
//        // Surch the thresholding level
//        if ((ABS(MR_Step(w)) > FLOAT_EPSILON)&&(ABS(MR_Step(w)) <LevelT ))
// 	 LevelT=ABS(MR_Step(w));
//      }
// 	  
//    while (Iter<MaxIter)
//    // while (((Delta_F*Delta_F) > .001 )&&(Iter<MaxIter))
//    {
//        MR_Step.recons(Result);
//        
//        // compute the Laplacian
//        for (int i=0;i<Nl;i++)
//        for (int j=0;j<Nc;j++)
// 	   Lap(i,j) =  4*Result(i,j) - (Result(i-1,j,Bord) +
// 					Result(i,j-1, Bord) + 
// 					Result(i,j+1, Bord) +
// 					Result(i+1,j, Bord));
// 
//        // positivity constraint on reconstructed image
//        threshold(Result);
//        
//        // compute the multiresolution transform of the laplacian 
//        MR_Iter.transform(Lap); 
//        F=0.;   
//        // compute the  new solution
//        for (w=0; w < MR_Step.nbr_mr_coeff(); w++)
// 	 {
// 	   //main constraint
// 	   MR_Step(w) -= 0.2 * MR_Iter(w);
// 	   F += MR_Iter(w);
// 
// 	   // range constraint on thresholded-restored coefficients 
// 	   if ((ABS(MR_Step(w))>LevelT) && (TabSup[w] == 0)) 
// 	             MR_Step(w) = (MR_Step(w)>0.0)? LevelT:-LevelT;
// 	   
// 	   // range constraint on quantized-restored coefficients 
// 	   MR_Step.pos_coeff(w, Scale, Indi, Indj, Det);
// 	   switch (Det)
// 	   {
// 	      case D_HORIZONTAL: LevelQ =  Quant[0]; break;
// 	      case D_VERTICAL: LevelQ =  Quant[1]; break;
// 	      case D_DIAGONAL: LevelQ =  Quant[2]; break;
// 	      default: LevelQ = 1;
// 	   }
// 	   LevelQ = LevelQ / 2.;
// 	   Delta_w =  MR_Step(w) - TabSave[w];
// 	   if ((ABS(Delta_w) > LevelQ/2.)&&(TabSup[w]==1)) 
// 	     MR_Step(w) = (Delta_w>0.0)? TabSave[w]+LevelQ:TabSave[w]-LevelQ;
// 	 }
//        if (MR_Step.nbr_mr_coeff()!=0)  
// 	 Delta_F = (F-F_old)/ MR_Step.nbr_mr_coeff();
//        else  Delta_F = 0.;
// 
//        // cerr << "Iter = " << Iter+1<< " F = " << F << "  Delta_F =" <<  Delta_F*Delta_F <<  endl;
//        F_old=F;
//        Iter++;
//    }
//    MR_Step.recons(Result);
//    // positivity constraint on reconstructed image
//    threshold(Result);
// 
//    delete TabSup;
//    delete TabSave;  
// }


void mr_ortho_regul_ima_rec(MultiResol &MR_Step, Ifloat &Result, 
                            float *TabLevel, int MaxIter, int SizeTabLevel)
// MR_Step = in: orthogonal transform with 2 scales (4 bands)
// Result = out: reconstructed image under laplacian constaint
// MaxIter = in: number of iterations
// TabLevel = in: detection level per scale
// SizeTabLevel = in: size of TabLevel
//                if SizeTabLevel = 1 the level is the same for the 3 bands
//                                        i.e. TabLevel[0]
//                if SizeTabLevel = 3, one level per band
//                if SizeTabLevel > 3, one level per wavelet coefficient
{
   int Nl = Result.nl();
   int Nc = Result.nc();
   int b,i,j,w,Iter=0;
   Ifloat Lap(Nl,Nc, "Dat");
   MultiResol MR_Iter (Nl, Nc, MR_Step.nbr_scale(), MR_Step.Type_Transform, "MR_Iter");
   MR_Iter.SBFilter = MR_Step.SBFilter;
   MR_Iter.TypeNorm =  MR_Step.TypeNorm;
   MR_Iter.LiftingTrans = MR_Step.LiftingTrans;
   MR_Iter.NormHaar = MR_Step.NormHaar;
   MR_Iter.Border = MR_Step.Border;

   int Nw = MR_Step.nbr_mr_coeff();
   unsigned char *TabSup = new unsigned char [Nw];
   float F;
   float F_old=0; 
   float Delta_F=10;;
   float Level;
   
   type_border Bord=I_CONT;
   float *TabSave = new float [Nw];
   float Delta_w,LevelQ,LevelT=INFINITY;
     
   w =0;
   for (b=0; b < MR_Step.nbr_band()-1; b++)
   for (i=0; i < MR_Step.size_band_nl(b); i++)
   for (j=0; j < MR_Step.size_band_nc(b); j++)
   {
      // Save initial value of dequantized coefficients for range constraint
      TabSave[w]=MR_Step(b,i,j);
      float Val = ABS(TabSave[w]);
      
      // compute the multiresolution support
      TabSup[w] = (Val > FLOAT_EPSILON) ? 1 : 0; 
      w++;
      
      // Surch the thresholding level
      if ((Val > FLOAT_EPSILON) && (Val < LevelT)) LevelT = Val;
   }
	  
   while (Iter <  MaxIter)
   // while (((Delta_F*Delta_F) > .001 )&&(Iter<MaxIter))
   {
       // cerr << "Iter = " << Iter+1<<  Nl << " " << Nc << endl;

       MR_Step.recons(Result);
       
       // compute the Laplacian
       for (i=0;i<Nl;i++)
       for (j=0;j<Nc;j++)
	   Lap(i,j) =  4*Result(i,j) - (Result(i-1,j,Bord) +
					Result(i,j-1, Bord) + 
					Result(i,j+1, Bord) +
					Result(i+1,j, Bord));

       // positivity constraint on reconstructed image
       threshold(Result);
       
       // compute the multiresolution transform of the laplacian 
       MR_Iter.transform(Lap); 
       F=0.;   
       // compute the  new solution
       w=0;
       for (b=0; b < MR_Step.nbr_band()-1; b++)
       for (i=0; i < MR_Step.size_band_nl(b); i++)
       for (j=0; j < MR_Step.size_band_nc(b); j++)
       {
          if (SizeTabLevel == 1) Level = TabLevel[0];
	  else if (SizeTabLevel == MR_Step.nbr_band()-1) Level = TabLevel[b];
          else Level = TabLevel[w];
	  
	  //main constraint
	  MR_Step(b,i,j) -= 0.2 * MR_Iter(b,i,j);
	  F += MR_Iter(b,i,j);

	  // range constraint on thresholded-restored coefficients 
	  if ((ABS(MR_Step(b,i,j))>LevelT) && (TabSup[w] == 0)) 
	             MR_Step(b,i,j) = (MR_Step(b,i,j)>0.0)? LevelT:-LevelT;
	   
	   // range constraint on quantized-restored coefficients 
	   LevelQ =  Level / 2.;
  	   Delta_w =  MR_Step(b,i,j) - TabSave[w];
	   if ((ABS(Delta_w) > LevelQ/2.)&&(TabSup[w]==1)) 
	     MR_Step(b,i,j) = (Delta_w>0.0)? TabSave[w]+LevelQ:TabSave[w]-LevelQ;
           w++;
       }
       if (MR_Step.nbr_mr_coeff()!=0)  
	 Delta_F = (F-F_old)/ MR_Step.nbr_mr_coeff();
       else  Delta_F = 0.;

       // cerr << "Iter = " << Iter+1<< " F = " << F << "  Delta_F =" <<  Delta_F*Delta_F <<  endl;
       F_old=F;
       Iter++;
   }
   MR_Step.recons(Result);
   // positivity constraint on reconstructed image
   threshold(Result);

   delete [] TabSup;
   delete [] TabSave;  
}
 
/****************************************************************************/

void mr_ortho_regul_rec (MultiResol &MR_Data, Ifloat &Result, 
              float *LevelQ, int MaxIter, Bool ModifCoef)
{
   Ifloat Imag;
   int NbrPlan = MR_Data.nbr_scale();
   int Nb = MR_Data.nbr_band();
   int b,s,Nls,Ncs;
   MultiResol MR_Step;
   float B[3];
 
   s = NbrPlan-2;
   Nls = MR_Data.size_scale_nl(s, I_SMOOTH);
   Ncs = MR_Data.size_scale_nc(s, I_SMOOTH); 
   Result.resize(Nls,Ncs);
   if ((Nls != MR_Data.band(Nb-1).nl()) || (Ncs != MR_Data.band(Nb-1).nc()))
   {
       cerr << "Error: band size are not the same in mr_ortho_regul_rec ... " << endl;
       cerr  << "   Nls = " << Nls << " Ncs = " << Ncs << endl;
       exit(-1);
   }   
   Result = MR_Data.band(Nb-1);
 
   for (s = NbrPlan-2; s >= 0; s--)
   {
      Nls = (s > 0) ? MR_Data.size_scale_nl(s-1, I_SMOOTH): MR_Data.size_ima_nl();
      Ncs = (s > 0) ? MR_Data.size_scale_nc(s-1, I_SMOOTH): MR_Data.size_ima_nc(); 
 
      MR_Step.alloc(Nls, Ncs, 2, MR_Data.Type_Transform, "MR_Step");
      MR_Step.band(3) = Result;
      MR_Step.SBFilter = MR_Data.SBFilter;
      MR_Step.TypeNorm =  MR_Data.TypeNorm;
      MR_Step.LiftingTrans =  MR_Data.LiftingTrans;
      MR_Step.NormHaar =  MR_Data.NormHaar;
      MR_Step.Border =  MR_Data.Border;
      for (b=0; b < MR_Step.nbr_band()-1; b++) 
      {
          if ((MR_Data.band(3*s+b).nl() != MR_Step.band(b).nl()) ||
	      (MR_Data.band(3*s+b).nc() != MR_Step.band(b).nc()))
	  {
	     cerr << "Error: band size are not the same in filter_noiter ... " << endl;
	     exit(-1);
	  }
	 MR_Step.band(b) = MR_Data.band(3*s+b);
	 B[b] = LevelQ[3*s+b];
      }
      
      Result.resize(Nls,Ncs);
      if (MaxIter > 1) mr_ortho_regul_ima_rec(MR_Step, Result, B, MaxIter, 3);
      else MR_Step.recons(Result); 
      MR_Step.free(); 
   }   
   if (ModifCoef == True) MR_Data.transform(Result);
}

/****************************************************************************/
