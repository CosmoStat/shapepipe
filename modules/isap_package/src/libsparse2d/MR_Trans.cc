/******************************************************************************
**                   Copyright (C) 1994 by CEA
*******************************************************************************
**
**    UNIT
**
**    Version: 3.3
**
**    Author: Jean-Luc Starck
**
**    Date:  96/06/13 
**    
**    File:  MR_Trans.cc
**
*******************************************************************************
**
**    DESCRIPTION : multiresolution transform procedures
**    -----------  
**
*******************************************************************************
**
**  void mr_transform (Ifloat &Image, MultiResol &MR_Transf, type_border Border, **  Bool Details)
**
** compute a multiresolution transform of the image Image
** 
*******************************************************************************
**
**  void MultiResol::transform (Ifloat &Image, type_border Border, Bool Details)
**
**  multiresoluiton class menber:
** compute a multiresolution transform of the image Image
** by calling mr_transform
**
*******************************************************************************
**
** void mr_scaling_function (Ifloat &Imag, MultiResol &MR_Data, int Max_Iter)
**
** compute a positive multiresolution transform
** MR_Data should be a pyramidal transform of Imag
** MR_Data is modified in order to have only positive coefficients.
** the reconstruction of MR_Data gives a good approximation of Imag
** !!! no transform will give again MR_Data
** 
** the positive coefficient are found by iterating. Max_Iter is the number of
** iteration
**
******************************************************************************/

// static char sccsid[] = "@(#)MR_Trans.cc 3.3 96/06/13 CEA 1994 @(#)";

#include "MR_Obj.h"
#include "IM_IO.h"
#include "IM_Sigma.h"
// #include "SB_Spline.h"

#define PRINT_DATA 0
#define DEBUG_TROBJ 0


static double TabNormPyrB3Spline[MAX_SCALE] = 
    {0.88937,  0.200924, 0.08544, 0.0410743, 0.0208055, 0.0108577, 0.00528016,
     0., 0., 0.};     
     
static double TabNormPyrMedian5[MAX_SCALE] = 
    {0.970291,  0.338748, 0.178373, 0.0992649, 0.0539517, 0.0317369, 0.0188998,
      0.013650, 0., 0.}; 
            
// static double TabNormPyrMedian3[MAX_SCALE] = 
//     {0.999234,  0.438327, 0.2347, 0.132756, 0.0757592, 0.0466552, 0.0326198,
//      0.0401325, 0., 0};      

static void mr_transform (Ifloat &Image, MultiResol &MR_Transf, Bool EdgeLineTransform, type_border Border=DEFAULT_BORDER, Bool Details=True);
     
/****************************************************************************/

void mr_correct_pyr (Ifloat &Imag, MultiResol &MR_Data, int Max_Iter)
{
   int Iter = 0, s,i,j;
   int Nl = Imag.nl();
   int Nc = Imag.nc();
   Ifloat Resi (Nl, Nc, "Residual");
   int Nbr_Plan = MR_Data.nbr_scale();
   type_transform Transform = TM_PYR_MINMAX;
   MultiResol MR_Oper (Nl, Nc, Nbr_Plan, Transform, "MR_Oper");
   float Dist;

   cout << "Max_Iter = " << Max_Iter << endl;

   for (Iter = 0; Iter < Max_Iter; Iter ++)
   {
       MR_Data.recons (Resi);
       Resi = Imag - Resi;
#if DEBUG_TROBJ
       cout << "  ITER " << Iter << endl;
       cout << "  Sigma Resi =  " << sigma(Resi) << endl;
       cout << "  Min Resi =  " << min(Resi) << endl;
       cout << "  Max Resi =  " << max(Resi) << endl;
#endif

       MR_Oper.transform (Resi);
       for (s = 0; s < Nbr_Plan - 1; s++)
       {
           Dist = 0;

#if DEBUG_TROBJ
cout << "  Min Scale  " << s+1 << " = " << min (MR_Oper.scale(s)) << endl;
cout << "  Max Scale  " << s+1 << " = " << max (MR_Oper.scale(s)) << endl;
#endif
           for (i = 0; i < MR_Data.scale(s).nl(); i ++)
           for (j = 0; j < MR_Data.scale(s).nc(); j ++)
           {
               Dist += ABS(MR_Oper(s,i,j));
               MR_Data(s,i,j) += MR_Oper(s,i,j);
           }
#if PRINT_DATA
cout << "Iter: " << Iter+1 << "=> DIST Scale  " << s+1 << " = " << Dist << endl;
#endif
       }
#if PRINT_DATA
cout << endl;
#endif
       s = Nbr_Plan-1;
       MR_Data.scale(s) += MR_Oper.scale(s);
   }

#if PRINT_DATA
   MR_Data.recons (Resi);
   Resi = Imag - Resi;
   cout << "  Min Resi =  " << min(Resi) << endl;
   cout << "  Max Resi =  " << max(Resi) << endl;
#endif

}

/****************************************************************************/

static void mr_scaling_function (Ifloat &Imag, MultiResol &MR_Data, int Max_Iter)
{
   int Iter = 0, s,i,j;
   int Nl = Imag.nl();
   int Nc = Imag.nc();
   Ifloat Resi (Nl, Nc, "Residual");
   int Nbr_Plan = MR_Data.nbr_scale();
   type_transform Transform = TM_PYR_MINMAX;
   MultiResol MR_Oper (Nl, Nc, Nbr_Plan, Transform, "MR_Oper");
   float Dist;

   // cout << "Max_Iter = " << Max_Iter << endl;

   for (Iter = 0; Iter < Max_Iter; Iter ++)
   {
       MR_Data.recons (Resi);
       Resi = Imag - Resi;
#if DEBUG_TROBJ
       cout << "  ITER " << Iter << endl;
       cout << "  Sigma Resi =  " << sigma(Resi) << endl;
       cout << "  Min Resi =  " << min(Resi) << endl;
       cout << "  Max Resi =  " << max(Resi) << endl;
#endif

       MR_Oper.transform (Resi);
       for (s = 0; s < Nbr_Plan; s++)
       {
           Dist = 0;

#if DEBUG_TROBJ
cout << "  Min Scale  " << s+1 << " = " << min (MR_Oper.scale(s)) << endl;
cout << "  Max Scale  " << s+1 << " = " << max (MR_Oper.scale(s)) << endl;
#endif
           for (i = 0; i < MR_Data.scale(s).nl(); i ++)
           for (j = 0; j < MR_Data.scale(s).nc(); j ++)
           {
               Dist += ABS(MR_Oper(s,i,j));
               MR_Data(s,i,j) += MR_Oper(s,i,j);
               if (MR_Data(s,i,j) < 0) MR_Data(s,i,j) = 0;
           }
#if PRINT_DATA
cout << "Iter: " << Iter+1 << "=> DIST Scale  " << s+1 << " = " << Dist << endl;
#endif
       }
#if PRINT_DATA
cout << endl;
#endif
       // s = Nbr_Plan-1;
       // MR_Data.scale(s) += MR_Oper.scale(s);
   }

#if PRINT_DATA
   MR_Data.recons (Resi);
   io_write_ima_float("s_pos_imag", Resi);
   Resi = Imag - Resi;
   io_write_ima_float("s_pos_resi", Resi);
   cout << "  Min Resi =  " << min(Resi) << endl;
   cout << "  Max Resi =  " << max(Resi) << endl;
#endif

}

/****************************************************************************/

static void rmedian (Ifloat &Data, Ifloat & Med, type_border Border, int Size)
{
   smooth_mediane (Data, Med,  Border, 0, Size);
}

/****************************************************************************/

static void mr_pyr_tmto (Ifloat &Image, MultiResol &MR_Transf, float & SigmaNoise, 
                 float NSigBSpline, float NSigMedian)
{
    int i,j,s,Nbr_Plan = MR_Transf.nbr_band();
    int Nl = Image.nl();
    int Nc = Image.nc();
    float Sigmaj,LevelScaleBspline,LevelScaleMedian;
    type_border Border = MR_Transf.Border;
    Ifloat BuffMed(Nl,Nc,"BuffMed");
    Ifloat BuffB3(Nl,Nc,"BuffMed");
    float Coeff_h0 = 3. / 8.;
    float Coeff_h1 = 1. / 4.;
    float Coeff_h2 = 1. / 16.;
    int Step = 1;   
       
    // find the noise level
    if ( SigmaNoise < FLOAT_EPSILON) SigmaNoise = detect_noise_from_med (Image);

    // intialization
    MR_Transf.band(0) = Image;

    for (s = 0; s < Nbr_Plan -1; s++)
    {
       int Nls = MR_Transf.size_band_nl(s);
       int Ncs = MR_Transf.size_band_nc(s);
       BuffMed.resize(Nls, Ncs);
       BuffB3.resize(Nls, Ncs);
       
       // detection level in case of B3 spline wavelet transform
       Sigmaj = SigmaNoise*TabNormPyrB3Spline[s];  
       if (s == 0) LevelScaleBspline = (NSigBSpline+1)*Sigmaj ;
       else LevelScaleBspline =  NSigBSpline*Sigmaj;

       // detection level in case of  pyramidal median tranform
       LevelScaleMedian =  NSigMedian*SigmaNoise*TabNormPyrMedian5[s];
        
       // median filtering
       rmedian (MR_Transf.band(s), BuffMed, Border, MR_Transf.MedianWindowSize);
        
       // Replace in BuffMed unsignificant pixels    
       //  (i.e.  | where scale(w) - median(w) | < LevelScaleMedian)
       // by their non-filtered value ==> only high values are
       // filtered bu the median filtering
       for (i = 0; i < Nls; i++) 
       for (j = 0; j < Ncs; j++) 
            if (ABS(MR_Transf(s,i,j) - BuffMed(i,j)) < LevelScaleMedian) 
	                                   BuffMed(i,j) = MR_Transf(s,i,j);
        
  
       // smooth BuffMed by a B3 spline
       // high peaks have removed by the median filtering
       // result is store in BuffMed and coefficient in MR_Transf(s,i,j)
       for (i = 0; i < Nls; i ++)
       for (j = 0; j < Ncs; j ++)
          BuffB3(i,j) = Coeff_h0 *  BuffMed(i,j)
                 + Coeff_h1 * (   BuffMed(i,j-Step, Border) 
                               +  BuffMed(i,j+Step, Border)) 
                 + Coeff_h2 * (   BuffMed(i,j-2*Step, Border) 
                               +  BuffMed(i,j+2*Step, Border));
      for (i = 0; i < Nls; i ++)
      for (j = 0; j < Ncs; j ++)
      {
          BuffMed(i,j) = Coeff_h0 * BuffB3(i,j)
                 + Coeff_h1 * (  BuffB3 (i-Step, j,  Border) 
                               + BuffB3 (i+Step, j,  Border)) 
                 + Coeff_h2 * (  BuffB3 (i-2*Step, j,  Border) 
                               + BuffB3 (i+2*Step, j,  Border));        
         MR_Transf(s,i,j) -= BuffMed(i,j); 
      }
      
      // reduce the size of the smoothed array by 2
      im_reduce_size_2 (BuffMed, MR_Transf.scale(s+1));
   }
}


/****************************************************************************/

static void mr_semi_pyr_tmto (Ifloat &Image, MultiResol &MR_Transf, float & SigmaNoise, 
                 float NSigBSpline, float NSigMedian)
{
    int i,j,s,Nbr_Plan = MR_Transf.nbr_band();
    int Nl = Image.nl();
    int Nc = Image.nc();
    float Sigmaj,LevelScaleBspline,LevelScaleMedian;
    type_border Border = MR_Transf.Border;
    Ifloat BuffMed(Nl,Nc,"BuffMed");
    Ifloat BuffB3(Nl,Nc,"BuffMed");
        
    // find the noise level
    if ( SigmaNoise < FLOAT_EPSILON) SigmaNoise = detect_noise_from_med (Image);
    // cout << "SigmaNoise = " << SigmaNoise << endl;
    if (SigmaNoise == 0) SigmaNoise = FLOAT_EPSILON;
    
    MR_Transf.band(0) = Image;
    for (s = 0; s < Nbr_Plan -1; s++)
    {
       int Nls = MR_Transf.size_band_nl(s);
       int Ncs = MR_Transf.size_band_nc(s);
       BuffMed.resize(Nls, Ncs);
       BuffB3.resize(Nls, Ncs);
       
       // detection level in case of B3 spline wavelet transform
       Sigmaj = SigmaNoise*MR_Transf.band_norm(s);   //  TabNormPyrB3Spline[s];  
       if (s == 0) LevelScaleBspline = (NSigBSpline+1)*Sigmaj ;
       else LevelScaleBspline =  NSigBSpline*Sigmaj;
        
       // median filtering
       if (s == 0)
       {
        // LevelScaleMedian =  NSigMedian*SigmaNoise*TabNormPyrMedian5[s];
	LevelScaleMedian =  NSigMedian*SigmaNoise;
        rmedian (Image, BuffMed, I_MIRROR, 5);
       }
       else
       {
        //LevelScaleMedian =  NSigMedian*SigmaNoise*TabNormPyrMedian5[s];
	LevelScaleMedian =  NSigMedian*SigmaNoise*MR_Transf.band_norm(s-1);
        rmedian (MR_Transf.band(s), BuffMed, I_MIRROR, 9);
       }
	
       // Replace in BuffMed unsignificant pixels    
       //  (i.e.  | where scale(w) - median(w) | < LevelScaleMedian)
       // by their non-filtered value ==> only high values are
       // filtered bu the median filtering
       for (i = 0; i < Nls; i++) 
       for (j = 0; j < Ncs; j++) 
            if (ABS(MR_Transf(s,i,j) - BuffMed(i,j)) < LevelScaleMedian) 
	                                   BuffMed(i,j) = MR_Transf(s,i,j);
        
       // smooth BuffMed by a B3 spline
       // high peaks have removed by the median filtering
       // result is store in BuffMed and coefficient in MR_Transf(s,i,j)
       if (s > 0)  
       {
          smooth_bspline (BuffMed , BuffB3, Border, 1);
          // smooth_bspline (BuffB3 , BuffMed, Border, 1);
          // MR_Transf.band(s) -= BuffMed;
	  // im_reduce_size_2 (BuffMed, MR_Transf.band(s+1));
	  MR_Transf.band(s) -= BuffB3;
	  im_reduce_size_2 (BuffB3, MR_Transf.band(s+1));
       }
       else 
       {
           smooth_bspline (BuffMed , BuffB3, Border, 0);
           MR_Transf.band(s)  -=  BuffB3;
           MR_Transf.band(s+1) = BuffB3;
       }
    }
}

/*********************************************************************/

static void mr_transform (Ifloat &Image, MultiResol &MR_Transf, Bool EdgeLineTransform, type_border Border, Bool Details)
{
    int Nbr_Plan = MR_Transf.nbr_scale();
    void mr_tf_transform(Ifloat &, MultiResol &);
    int s;
    void im_min_reduce_size_2 (const Ifloat &, Ifloat &);
    int MedianWindowSize = MR_Transf.MedianWindowSize;
    int Nl = MR_Transf.size_ima_nl();
    int Nc = MR_Transf.size_ima_nc();
    Ifloat Buff, Buff_Dilat;
    
    MR_Transf.Border = Border;
    switch (MR_Transf.Type_Transform)
    {
       case TO_PAVE_LINEAR:
           MR_Transf.band(0) = Image;
           for (s = 0; s < Nbr_Plan -1; s++)
           {
               smooth_linear (MR_Transf.band(s),MR_Transf.band(s+1),Border,s); 
               if (Details == True) 
	                MR_Transf.band(s) -= MR_Transf.band(s+1);
           }
           break; 
       case TO_PAVE_BSPLINE:
            MR_Transf.band(0) = Image;
            for (s = 0; s < Nbr_Plan -1; s++)
            {
               smooth_bspline (MR_Transf.band(s),MR_Transf.band(s+1),Border,s);
               if  (Details == True) MR_Transf.band(s) -= MR_Transf.band(s+1);
            }
            break;
       case TO_PAVE_HAAR:
           MR_Transf.band(0) = Image;
           for (s = 0; s < Nbr_Plan -1; s++)
           {
                smooth_haar (MR_Transf.band(s),MR_Transf.band(s+1),Border,s,MR_Transf.NormHaar);
                if  (Details == True) MR_Transf.band(s) -= MR_Transf.band(s+1);
           }
           break;
       case TM_PAVE_MEDIAN:
           MR_Transf.band(0) = Image;
           //printf ("MedianWindowSize = %d \n", MedianWindowSize);
           for (s = 0; s < Nbr_Plan -1; s++)
           {
                int Size, Demi=MedianWindowSize/2;
                Size = POW2(s)*2*Demi+1;
                rmedian(MR_Transf.band(s), MR_Transf.band(s+1),Border, Size);
                if  (Details == True) MR_Transf.band(s) -= MR_Transf.band(s+1);
           } 
           break;
       case TM_PAVE_MINMAX:
 	   MR_Transf.scale(0) = Image;
           Buff.alloc(Nl, Nc, "Buff: TM_PAVE_MINMAX");
           for (s = 0; s < Nbr_Plan -1; s++)
           {
               Nl = MR_Transf.size_band_nl(s);  
               Nc = MR_Transf.size_band_nc(s);   
               Buff.resize(Nl, Nc);
               morpho_erosion (MR_Transf.band(s), Buff, 2*s+3);
               morpho_dilation (Buff, MR_Transf.band(s+1), 2*s+3);
               if   (Details == True) MR_Transf.band(s) -= MR_Transf.band(s+1);
           }
           break;
//        case TO_PYR_LINEAR: // TO_PYR_UNSER: Unser Spline test
//           {
// 	    SPLINE_SUBBAND SSB;
//             // FIL_SPLINE, FIL_SPLINE_L2, FIL_SPLINE_CENT, FIL_SPLINE_L2_CENT 
// 	    type_spline_filter Tsf = FIL_SPLINE_CENT;
// 	    int SplineOrder = 1;
// 	    SSB.alloc(Tsf, SplineOrder);
// 	    MR_Transf.scale(0) = Image;
//             Buff.alloc(Nl, Nc, "Buff: TM_PYR_MEDIAN");
//             for (s = 0; s < Nbr_Plan -1; s++)
//             {
// 	       cout << endl << "SSB Band " << s+1 << endl;
//                Nl = MR_Transf.size_band_nl(s);  
//                Nc = MR_Transf.size_band_nc(s);
// 	       int Nls = MR_Transf.size_band_nl(s+1);
// 	       int Ncs = MR_Transf.size_band_nc(s+1);
//                Buff.resize(Nl, Nc);
// 	       float *InBuff = MR_Transf.band(s).buffer();
//                float *OutBuff = MR_Transf.band(s+1).buffer();
//  	       SSB.Reduce_2D(InBuff, Nc, Nl, OutBuff);
// 	       INFO(MR_Transf.band(s+1), "REDUCE");
// 	       SSB.Expand_2D(OutBuff, Ncs, Nls, Buff.buffer());
// 	       INFO(Buff, "EXPAND");
//                MR_Transf.band(s) -= Buff;
// 	       INFO(MR_Transf.band(s), "DIFF");
//             }  
// 	   }
//            break;
       case TO_PYR_LINEAR:
 	   MR_Transf.scale(0) = Image;
           Buff.alloc(Nl, Nc, "Buff: TO_PYR_LINEAR");
           for (s = 0; s < Nbr_Plan -1; s++)
           {
	       Nl = MR_Transf.size_band_nl(s);  
               Nc = MR_Transf.size_band_nc(s);   
               Buff.resize(Nl, Nc);
               smooth_linear (MR_Transf.band(s), Buff, Border);
               if   (Details == True) MR_Transf.band(s) -= Buff;
               im_reduce_size_2 (Buff, MR_Transf.band(s+1));
           }
           break;
       case TO_PYR_BSPLINE:
 	   MR_Transf.band(0) = Image;
	   Buff.alloc(Nl, Nc, "Buff: TO_PYR_BSPLINE");
           for (s = 0; s < Nbr_Plan -1; s++)
           {
               Nl = MR_Transf.size_band_nl(s);  
               Nc = MR_Transf.size_band_nc(s);  
               Buff.resize(Nl, Nc);
               smooth_bspline (MR_Transf.band(s), Buff, Border);
               if   (Details == True) MR_Transf.band(s) -= Buff;
               im_reduce_size_2 (Buff, MR_Transf.band(s+1));
           }
           break;
       case TO_SEMI_PYR:
           MR_Transf.band(0) = Image;
	   Buff.alloc(Nl, Nc, "Buff: TO_PYR_BSPLINE");
           for (s = 0; s < Nbr_Plan -1; s++)
           {
	       int ValTrou;
	       ValTrou = (s == 0) ? 0 : 1;
               Nl = MR_Transf.size_band_nl(s);  
               Nc = MR_Transf.size_band_nc(s);  
               Buff.resize(Nl, Nc);
               smooth_bspline (MR_Transf.band(s), Buff, Border, ValTrou);
               if   (Details == True) MR_Transf.band(s) -= Buff;
               if (s > 0) 
	            im_reduce_size_2 (Buff, MR_Transf.band(s+1));
	       else MR_Transf.band(s+1) = Buff;
           }
           break;
       case TM_TO_SEMI_PYR:
          mr_semi_pyr_tmto (Image, MR_Transf, MR_Transf.SigmaNoise,
	               MR_Transf.TMTO_NSigBSpline, MR_Transf.TMTO_NSigMedian);
           break;
       case TM_TO_PYR:
          mr_pyr_tmto (Image, MR_Transf, MR_Transf.SigmaNoise,
	               MR_Transf.TMTO_NSigBSpline,MR_Transf.TMTO_NSigMedian);
          break;
       case TM_PYR_MINMAX:
 	   MR_Transf.band(0) = Image;
           Buff.alloc(Nl, Nc, "Buff: TO_PYR_MINMAX");
           Buff_Dilat.alloc(Nl, Nc, "Buff_Dilat: TO_PYR_MINMAX");
           for (s = 0; s < Nbr_Plan -1; s++)
           {
	       Nl = MR_Transf.size_band_nl(s);  
               Nc = MR_Transf.size_band_nc(s);  
               Buff.resize(Nl, Nc);
               Buff_Dilat.resize(Nl, Nc);

               morpho_erosion (MR_Transf.band(s), Buff);
               morpho_dilation (Buff, Buff_Dilat);
               if   (Details == True) MR_Transf.band(s) -= Buff_Dilat;
               im_min_reduce_size_2  (Buff_Dilat, MR_Transf.band(s+1));
           }
           break;
       case TM_PYR_SCALING_FUNCTION:
 	   MR_Transf.band(0) = Image;
           Buff.alloc(Nl, Nc, "Buff: TO_PYR_MINMAX");
           Buff_Dilat.alloc(Nl, Nc, "Buff_Dilat: TO_PYR_MINMAX");

           for (s = 0; s < Nbr_Plan -1; s++)
           {
	       Nl = MR_Transf.size_band_nl(s);  
               Nc = MR_Transf.size_band_nc(s);
               Buff.resize(Nl, Nc);
               Buff_Dilat.resize(Nl, Nc);
               morpho_erosion (MR_Transf.scale(s), Buff);
               morpho_dilation (Buff, Buff_Dilat);
               im_min_reduce_size_2  (Buff_Dilat, MR_Transf.scale(s+1));
               //im_increase_size_2 (MR_Transf.scale(s+1), Buff_Dilat, Border);
               if   (Details == True) MR_Transf.scale(s) -= Buff_Dilat;
               //   threshold(MR_Transf.scale(s));
           }
           mr_scaling_function (Image, MR_Transf, MR_Transf.Nbr_Iter);
           break;
       case TM_PYR_MEDIAN:
 	   MR_Transf.band(0) = Image;
           Buff.alloc(Nl, Nc, "Buff: TM_PYR_MEDIAN");
           for (s = 0; s < Nbr_Plan -1; s++)
           {
               Nl = MR_Transf.size_band_nl(s);  
               Nc = MR_Transf.size_band_nc(s);
	       if (MR_Transf.ExactPyrRec == False)
               {
                   Buff.resize(Nl, Nc);
                   smooth_mediane (MR_Transf.band(s), Buff, Border, 0,
                                   MedianWindowSize);
                   if   (Details == True) MR_Transf.band(s) -= Buff;
                   im_reduce_size_2 (Buff, MR_Transf.band(s+1));
               }
               else
               {
                   Buff.resize(Nl, Nc);

                   smooth_mediane (MR_Transf.band(s), Buff, Border, 0,
                                   MedianWindowSize);
                   im_reduce_size_2 (Buff, MR_Transf.band(s+1));
                   im_increase_size_2 (MR_Transf.band(s+1), Buff, Border);
                   MR_Transf.band(s) -= Buff;
               }
           }
           break;
       case TM_PYR_LAPLACIAN:
 	   MR_Transf.band(0) = Image;
           Buff.alloc(Nl, Nc, "Buff: TM_PYR_LAPLACIAN");

           for (s = 0; s < Nbr_Plan -1; s++)
           {
 	       Nl = MR_Transf.size_band_nl(s);  
               Nc = MR_Transf.size_band_nc(s);
               Buff.resize(Nl, Nc);
               smooth_bspline (MR_Transf.band(s), Buff, Border, s);
               im_reduce_size_2 (Buff, MR_Transf.band(s+1));
               im_increase_size_2 (MR_Transf.band(s+1), Buff, Border);
               MR_Transf.band(s) -= Buff;
           }
           break;
       case TO_PAVE_FFT:
       case TO_PYR_FFT_DIFF_RESOL:
       case TO_PYR_FFT_DIFF_SQUARE:
           if (Details) wave_cf_transform (Image, MR_Transf);
           else mr_cf_transform (Image, MR_Transf);
           break;
       case TO_DIADIC_MALLAT:
           wave_2d_mallat_atrou (Image, MR_Transf.tabband(), Nbr_Plan, Border);
	   break;
       case TO_MALLAT:
       case TO_UNDECIMATED_MALLAT:
       case TO_DIV_1:
       case TO_DIV_2:
           {
	      // SubBandFilter WT1D(MR_Transf.SBFilter, MR_Transf.TypeNorm);
	      // cout << "TO_UNDECIMATED_MALLAT " << endl;
              FilterAnaSynt *FAS_line = MR_Transf.filter_bank_line();
              if (FAS_line == NULL)
              {
                  cout << "Error: filter FAS_line bank is not defined ... " << endl;
                  exit(-1);
              }
	      // cout << "TO " << endl;
	      FilterAnaSynt *FAS_col = MR_Transf.filter_bank_column();
              if (FAS_col == NULL)
              {
                  cout << "Error: filter bank FAS_col is not defined ... " << endl;
                  exit(-1);
              }
	      //   cout << "WT1D_L " << endl;
	      SubBandFilter WT1D_L(*FAS_line, MR_Transf.TypeNorm);
	      //    cout << "WT1D_C " << endl;
	      SubBandFilter WT1D_C(*FAS_col, MR_Transf.TypeNorm);
              WT1D_L.Border= MR_Transf.Border;
              WT1D_C.Border= MR_Transf.Border;
              HALF_DECIMATED_2D_WT P2WT(WT1D_L, WT1D_C);
	      //   cout << "RUN TR " << endl;
              P2WT.transform(Image, MR_Transf.tabband(), Nbr_Plan,  MR_Transf.nbr_undec_scale() );
	      //   cout << "END " << endl;
           }
           // haar_2d_atrou_dir3_transform (Image, MR_Transf);
	   break;
       case TO_LC:
              (MR_Transf.LC)->undec_transform (Image, MR_Transf.tabband(), MR_Transf.LC_NbrScaleCol, MR_Transf.LC_NbrScaleLine);
	       // Image.init();
 	      //(MR_Transf.LC)->undec_recons (MR_Transf.tabband(), Image, MR_Transf.LC_NbrScaleCol, MR_Transf.LC_NbrScaleLine);
	      //io_write_ima_float("xx.fits", Image);
	       break;
        case TC_FCT:
	    {
	       // cout << "TRANS " << endl;
	      
	       FCUR *FCT = MR_Transf.fast_curvelet_trans();
	       FCT->cur_trans(Image);
	       
	      // cout << "TRANSOK " << endl;
	       int IndBand=0;
	       for (s = 0; s < FCT->nbr_scale(); s++)
	       {
               for (int b=0; b < FCT->nbr_band(s); b++) 
	       {
	          for (int i=0; i < FCT->size_band_nl(s,b); i++)
		  for (int j=0; j < FCT->size_band_nc(s,b); j++) 
		  {
// 		      if (i >= MR_Transf.size_band_nl(IndBand)) cout << "ERROR: i = " << i << " SIZE_BAND_NL = " <<  MR_Transf.size_band_nl(IndBand) << endl;
// 		      if (j >= MR_Transf.size_band_nc(IndBand)) cout << "ERROR: j = " << j << " SIZE_BAND_NC = " <<  MR_Transf.size_band_nc(IndBand) << endl;
// 		      if (i >= (FCT->TabCF_Band)[s][b].nl()) cout << "ERROR: i = " << i << " SIZE_BAND_CUR_NL = " << (FCT->TabCF_Band)[s][b].nl() << endl;
// 		      if (j >= (FCT->TabCF_Band)[s][b].nc()) cout << "ERROR: j = " << j << " SIZE_BAND_CUR_NC = " << (FCT->TabCF_Band)[s][b].nc() << endl;
		      MR_Transf(IndBand, i, j) = (*FCT)(s,b,i,j); 
		  }
		  // (FCT->TabCF_Band)[s][b].init();
		 //  cout << "(" << s << "," << b << ") sigma = " << MR_Transf.band(IndBand).sigma() << " Nl = " << FCT->size_band_nl(s,b) << " " <<  FCT->size_band_nc(s,b) <<  " Nc = " << FCT->size_band_nc(s,b) << " " <<(FCT->TabCF_Band)[s][b].nl() << " " << (FCT->TabCF_Band)[s][b].nc() <<    endl;
	          IndBand ++;
	       }
	       }
// 	       IndBand=0;
// 	       for (s = 0; s < FCT->nbr_scale(); s++)
// 	       {
//                for (int b=0; b < FCT->nbr_band(s); b++) 
// 	       {
// 	          for (int i=0; i < FCT->size_band_nl(s,b); i++)
// 		  for (int j=0; j < FCT->size_band_nc(s,b); j++) (*FCT)(s,b,i,j) = MR_Transf(IndBand, i, j);
//  	          IndBand ++;
// 	       }
// 	       }
// 	       Ifloat D(Image.nl(), Image.nc(), "D");
// 	       FCT->cur_recons(D);
// 	       D -= Image;
// 	       INFO(D, "ERROR REC = ");
	       // cout << "ERROR REC = " << D.info("ERR") << endl;
 	     //  cout << "ETRANS CP " << endl;
	    }
	    break;
       case TO_PYR_MEYER:
       case TO_PYR_MEYER_ISOTROP:
           {
	      MEYER_WT *WT = MR_Transf.get_meyer_wt ();
//		  WT->init(4,64,64, True);
//		  cerr<<"plop";
	      WT->transform(Image, MR_Transf.tabband(), False);
	   }
	   break;
       case TO_UNDECIMATED_NON_ORTHO:
           {
	      UndecSubBandFilter *USF = MR_Transf.undec_filter_bank();
              if (USF == NULL)
              {
                  cout << "Error: undecimated filter bank is not defined ... " << endl;
                  exit(-1);
              }
              PAVE_2D_WT WT(*USF);  
              WT.transform(Image, MR_Transf.tabband(), Nbr_Plan);
           }
 	   break;
       case TO_DIADIC_HAAR:
           haar_2d_atrou_dir2_transform (Image, MR_Transf, EdgeLineTransform);
	   break;
       //case TO_MAX_DIADIC_MALLAT:
           //wave_max_2d_mallat_atrou (Image,MR_Transf.tabband(),Nbr_Plan,Border);
	   //break;
//        case TO_MALLAT:
//            {
// 	     // SubBandFilter WT1D(MR_Transf.SBFilter, MR_Transf.TypeNorm);
//              FilterAnaSynt *FAS = MR_Transf.filter_bank();
//              SubBandFilter WT1D(*FAS, MR_Transf.TypeNorm);
//              Ortho_2D_WT WT2D(WT1D);
//  	     Buff.alloc(Image.nl(), Image.nc(), "aux");
//              Buff = Image;
//              WT2D.transform(Buff, MR_Transf.nbr_scale());
//              ima_to_ortho_trans(MR_Transf, Buff);
//            }
//            break;
       case TO_LIFTING:
           {
	     Lifting Clift1D(MR_Transf.LiftingTrans);
             Ortho_2D_WT WT2D(Clift1D);
 	     Buff.alloc(Image.nl(), Image.nc(), "aux");
             Buff = Image;
             WT2D.transform(Buff, MR_Transf.nbr_scale());
             ima_to_ortho_trans(MR_Transf, Buff);
           }
            break;
       case TM_MIN_MAX:
           {
	      G_Minmax C_GTM;
              Ortho_2D_WT WT2D(C_GTM);
 	      Buff.alloc(Image.nl(), Image.nc(), "aux");
              Buff = Image;
              WT2D.transform(Buff, MR_Transf.nbr_scale());
              ima_to_ortho_trans(MR_Transf, Buff);
           }  
           break;
       case TO_FEAUVEAU:
 	   Buff.alloc(Image.nl(), Image.nc(), "aux");
	   transform_feauveau (Image,  Buff, MR_Transf.nbr_scale()-1);
	   ima_to_ortho_trans(MR_Transf,  Buff);
           break;
       case TO_PAVE_FEAUVEAU:
           transform_feauveau_trou (Image, MR_Transf);
           break;
       case TO_HAAR:
           haar_2d_transform (Image, MR_Transf);
           break;
       default:
           fprintf (stderr, "Error (proc. mr_iso_transform):Unknown transform\n");
           exit (-1);
           break;
    }
}


/****************************************************************************/

void MultiResol::filter_bank_alloc()
{
    if ((Type_Transform == TO_MALLAT) || (Type_Transform == TO_UNDECIMATED_MALLAT) || (Type_Transform == TO_LC))
    {
       if (FilterBank == NULL)
       {
           FilterBank = new FilterAnaSynt;
           FilterBank->Verbose = Verbose;
           FilterBank->alloc(SBFilter);
           // FilterBank = new FilterAnaSynt(SBFilter);
           FilterBankAlloc = True;
       }
       else if (SBFilter != (FilterBank->type_filter)())
       {
          cout << "Error: the allocated filter bank is different from the specified one ... " << endl;
          exit(-1);
       }
    }
    else if (Type_Transform == TO_DIV_1)
    {
       if ((FilterBank_Line == NULL) || (FilterBank_Column == NULL))
       {
          FilterBank_Line = new FilterAnaSynt;
          FilterBank_Line->alloc(F_5_3);
	  FilterBank_Column = new FilterAnaSynt;
          FilterBank_Column->alloc(F_4_4);
	}
     }
     else if (Type_Transform == TO_DIV_2)
     {
	 if ((FilterBank_Line == NULL) || (FilterBank_Column == NULL))
	 {
	      TypeNorm = DEF_SB_NORM;
              FilterBank_Line = new FilterAnaSynt;
              FilterBank_Line->alloc(F_4_4);
	      FilterBank_Column = new FilterAnaSynt;
              FilterBank_Column->alloc(F_5_3);
	 }
    }
    
    
   // if ((Type_Transform == TO_DIV_1) || (Type_Transform == TO_DIV_2))
   // {
   //  }
}

/****************************************************************************/

void MultiResol::transform (Ifloat &Image)
{
    filter_bank_alloc();
    if ((ModifiedATWT == True) && (Type_Transform == TO_PAVE_BSPLINE))
    {
       ATROUS_2D_WT AWT;
       AWT.ModifiedAWT = True;
       AWT.transform(Image, TabBand, Nbr_Plan);
    } 
    else mr_transform (Image, (*this), EdgeLineTransform, Border);
}

/****************************************************************************/

void MultiResol::transform (Ifloat &Image, type_border Bord, Bool Details)
{
#if DEBUG_TROBJ
fprintf (stderr, "TRANSFORM: of %s (%d,%d)\n",  
                    Image.Name_Obj, Image.nl(), Image.nc());
fprintf (stderr, "Nl = %d, Nc = %d\n" , Nl, Nc);
fprintf (stderr, "Nbr_Plan = %d\n" , Nbr_Plan);
fprintf (stderr, "Transform = %s\n", StringTransform(Type_Transform));
#endif
    Border = Bord;

    filter_bank_alloc();
    if ((ModifiedATWT == True) && (Type_Transform == TO_PAVE_BSPLINE))
    {
       ATROUS_2D_WT AWT;
       AWT.ModifiedAWT = True;
       AWT.transform(Image, TabBand, Nbr_Plan);
    } 
    else mr_transform (Image, (*this), EdgeLineTransform, Bord, Details);

#if DEBUG_TROBJ
fprintf (stderr, "END TRANSFORM\n");
#endif
}

/****************************************************************************/

