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
**    File:  MR_Recon.cc
**
*******************************************************************************
**
**    DESCRIPTION  multiresolution reconstruction
**    -----------  
**                 
**
*******************************************************************************
**
** void mr_recons (MultiResol &MR_Transf, Ifloat &Image, type_border Border)
**
** reconstruction of an image from its wavelet transform
**
*******************************************************************************
**
** void MultiResol::recons (Ifloat &Image, type_border Border)
**
** class member procedure
** reconstruction of an image from its wavelet transform
**
******************************************************************************/

#include "IM_Obj.h"
#include "MR_Obj.h"
// #include "SB_Spline.h"

/*******************************************************************/

static void mr_adjoint_rec (MultiResol &MR_Transf, Ifloat &Image, Bool UseLastScale, type_border Border)
{
    int	i,j,n;
    int Nl = MR_Transf.size_ima_nl();
    int Nc = MR_Transf.size_ima_nc();
    int NScale = MR_Transf.nbr_band();
    int Nb_ech = NScale - 1;
    Ifloat temp(Nl, Nc, "Rec adjoint pave");
    static int Pas=0;

    switch (MR_Transf.Type_Transform)
    {
       case TO_PAVE_BSPLINE:
 	 if (UseLastScale == True) Image = MR_Transf.band(Nb_ech);
	 else Image.init();
	 for (n=Nb_ech-1; n>= 0 ; n--)
	 {
	     smooth_bspline (Image, temp, Border, n);
	     for (i=0; i < Image.nl(); i++)
	     for (j=0; j < Image.nc(); j++) 
	          Image(i,j) = temp(i,j) +  MR_Transf(n,i,j);
	     //      pb + Image = temp + MR_Transf.band(n);
	     // or for the real adjoint operaot: 
	     //   Image = temp
	     //   smooth_bspline (MR_Transf.band(n), temp, Border, n);
	     //   Image += MR_Transf.band(n) - temp;
         }
	break;
       case  TO_PYR_BSPLINE:
       case  TO_SEMI_PYR:
	    Image.init ();
	    Image.resize (MR_Transf.size_band_nl(Nb_ech),
                      MR_Transf.size_band_nc(Nb_ech));
	    if (UseLastScale == True) Image = MR_Transf.band(Nb_ech);
	 
	    for (n=Nb_ech-1; n >= 0; n--)
	    {
               int Nls = MR_Transf.size_band_nl(n);
               int Ncs = MR_Transf.size_band_nc(n);
               temp.resize (Nls,Ncs);
               if ((n > 0) || (MR_Transf.Type_Transform != TO_SEMI_PYR))
	       {
	         im_increase_size_2 (Image, temp);
                 Image.resize (Nls,Ncs);
	       }
               else Image = temp;
               if (n != Nb_ech-1) 
               {
                 if ((MR_Transf.Set_Transform == TRANSF_PYR)  || (n==0))
                     smooth_bspline (temp, Image, Border, 0);
                else smooth_bspline (temp, Image, Border, 1);
               }
               else Image = temp;
               Image += MR_Transf.band(n);
	     }
 	break;
        default:
             if (Pas == 0)
             {
                cout << "Warning: adjoint reconstruction operator not implemeted " << endl;
	        cout << "         use a simple reconstruction ... " << endl;
                Pas = 1;
             }
	     if (UseLastScale == False) MR_Transf.band(Nb_ech).init();
	     MR_Transf.recons(Image);
 	break;
    }
}


/****************************************************************************/

static void mr_recons (MultiResol &MR_Transf, Ifloat &Image, 
                       type_border Border, Bool EdgeLineTransform)
{
    int Nl, Nc, s, i,j, Nbr_Plan = MR_Transf.nbr_scale();
    int Nbr_Iter;
    Ifloat Buffer;
    
    switch (MR_Transf.Type_Transform)
    {
       case TO_PAVE_LINEAR:
       case TO_PAVE_BSPLINE:
       case TO_PAVE_FEAUVEAU:
       case TM_PAVE_MEDIAN:
       case TM_PAVE_MINMAX:
       case TO_PAVE_HAAR:
          Image = MR_Transf.band(Nbr_Plan-1);
          for (s = Nbr_Plan -2; s >= 0; s--) 
                          Image += MR_Transf.band(s);
          break;
//        case TO_PYR_LINEAR:  Unser test Spline
//            {
// 	     SPLINE_SUBBAND SSB;
//              // FIL_SPLINE, FIL_SPLINE_L2, FIL_SPLINE_CENT, FIL_SPLINE_L2_CENT 
// 	     type_spline_filter Tsf = FIL_SPLINE_CENT;
// 	     int SplineOrder = 1;
// 	     SSB.alloc(Tsf, SplineOrder);
//              Buffer.alloc(MR_Transf.size_ima_nl(), MR_Transf.size_ima_nc(),
//                           "Buffer mr_iso_recons");
//  	     Nl = MR_Transf.size_band_nl(Nbr_Plan-1);  
//              Nc = MR_Transf.size_band_nc(Nbr_Plan-1);  
//              Image.resize (Nl, Nc);
//              Image = MR_Transf.band(Nbr_Plan-1);
//              for (s = Nbr_Plan -2; s >= 0; s--)
//              {
//  		 Nl = MR_Transf.size_band_nl(s);  
//                  Nc = MR_Transf.size_band_nc(s);   
// 	         int Nls = MR_Transf.size_band_nl(s+1);
// 	         int Ncs = MR_Transf.size_band_nc(s+1);
//                  Buffer.resize (Nl, Nc);
//  	         float *InBuff = Image.buffer();
//                  float *OutBuff = Buffer.buffer();   
// 		 SSB.Expand_2D(InBuff, Ncs, Nls, OutBuff);           
//                  Image.resize (Nl, Nc);
// 		 for (i = 0; i < Nl; i++)
//                  for (j = 0; j < Nc; j++)
//                            Image(i,j) = Buffer(i,j) + MR_Transf(s,i,j);
//               }
// 	    }
//            break;
       case TO_PYR_LINEAR:
       case TO_PYR_BSPLINE:
       case TM_PYR_MEDIAN:
       case TM_PYR_LAPLACIAN:       
       case TM_TO_PYR:       
       case TO_SEMI_PYR:
       case TM_TO_SEMI_PYR:
              Buffer.alloc(MR_Transf.size_ima_nl(), MR_Transf.size_ima_nc(),
                         "Buffer mr_iso_recons");
 	     Nl = MR_Transf.size_band_nl(Nbr_Plan-1);  
             Nc = MR_Transf.size_band_nc(Nbr_Plan-1);  
             Image.resize (Nl, Nc);
             Image = MR_Transf.band(Nbr_Plan-1);
             for (s = Nbr_Plan -2; s >= 0; s--)
             {
 		 Nl = MR_Transf.size_band_nl(s);  
                 Nc = MR_Transf.size_band_nc(s);   
                 Buffer.resize (Nl, Nc);
                 if ((s > 0) || (MR_Transf.Set_Transform != TRANSF_SEMIPYR))
		 {
		    im_increase_size_2 (Image, Buffer, Border);
		    // im_bilinear_interp(Image, Buffer);
                    Image.resize (Nl, Nc);
		    for (i = 0; i < Nl; i++)
                    for (j = 0; j < Nc; j++)
                           Image(i,j) = Buffer(i,j) + MR_Transf(s,i,j);
		 }
		 else Image += MR_Transf.band(s);
              }
           break;
       case TM_PYR_MINMAX:
       case TM_PYR_SCALING_FUNCTION:
              Buffer.alloc( MR_Transf.size_ima_nl(), MR_Transf.size_ima_nc(),
                          "Buffer mr_iso_recons");

             Nl = MR_Transf.size_band_nl(Nbr_Plan-1);  
             Nc = MR_Transf.size_band_nc(Nbr_Plan-1);  
             Image.resize (Nl, Nc);
             Image = MR_Transf.band(Nbr_Plan-1);
             for (s = Nbr_Plan -2; s >= 0; s--)
             {
                 Nl = MR_Transf.size_band_nl(s);  
                 Nc = MR_Transf.size_band_nc(s); 
                 Buffer.resize (Nl, Nc);
                 im_bilinear_interp (Image, Buffer);
                 Image.resize (Nl, Nc);
                 for (i = 0; i < Nl; i++)
                 for (j = 0; j < Nc; j++)
                           Image(i,j) = Buffer(i,j) + MR_Transf(s,i,j);
             }
           break;
       case TO_PAVE_FFT:
       case TO_PYR_FFT_DIFF_RESOL:
       case TO_PYR_FFT_DIFF_SQUARE:
          wave_cf_recons (MR_Transf, Image);
          break;
       case TO_LIFTING:
           {
 	     ortho_trans_to_ima(MR_Transf, Image);
	     Lifting Clift1D(MR_Transf.LiftingTrans);
             Ortho_2D_WT WT2D(Clift1D);
             WT2D.recons(Image, MR_Transf.nbr_scale());
            }        
	    break;  
//        case TO_MALLAT:
//            {
//  	     ortho_trans_to_ima(MR_Transf, Image);
// 	     // SubBandFilter WT1D(MR_Transf.SBFilter, MR_Transf.TypeNorm);
//              FilterAnaSynt *FAS = MR_Transf.filter_bank();
// 	     SubBandFilter WT1D(*FAS, MR_Transf.TypeNorm);
//              Ortho_2D_WT WT2D(WT1D);
//              WT2D.recons(Image, MR_Transf.nbr_scale());
//            }       
//            break;
       case TO_DIADIC_MALLAT:
          rec_wave_2d_mallat_atrou (MR_Transf.tabband(),  Image,   
                     Nbr_Plan, Border);   
            break;
       case TO_DIADIC_HAAR:
          haar_2d_atrou_dir2_recons (MR_Transf, Image, EdgeLineTransform);
 	  break; 
       case TO_MALLAT:
       case TO_UNDECIMATED_MALLAT:
       case TO_DIV_1:
       case TO_DIV_2:
           {	
              FilterAnaSynt *FAS_line = MR_Transf.filter_bank_line();
              if (FAS_line == NULL)
              {
                  cout << "Error: filter bank is not defined ... " << endl;
                  exit(-1);
              }
	      FilterAnaSynt *FAS_col = MR_Transf.filter_bank_column();
              if (FAS_col == NULL)
              {
                  cout << "Error: filter bank is not defined ... " << endl;
                  exit(-1);
              }
	      SubBandFilter WT1D_L(*FAS_line, MR_Transf.TypeNorm);
	      SubBandFilter WT1D_C(*FAS_col, MR_Transf.TypeNorm);
              HALF_DECIMATED_2D_WT P2WT(WT1D_L, WT1D_C);
   	      P2WT.recons(MR_Transf.tabband(), Image,Nbr_Plan, MR_Transf.nbr_undec_scale());
           }
  	  break; 
      case TO_LC:
           // (MR_Transf.LC)->undec_transform (Image, MR_Transf.tabband(), MR_Transf.LC_NbrScaleCol, MR_Transf.LC_NbrScaleLine);
            (MR_Transf.LC)->undec_recons (MR_Transf.tabband(), Image, MR_Transf.LC_NbrScaleCol, MR_Transf.LC_NbrScaleLine);
	     break;
      case TC_FCT:
	    {
	       FCUR *FCT = MR_Transf.fast_curvelet_trans();
 	       int IndBand=0;
	       for (s = 0; s < FCT->nbr_scale(); s++)
               for (int b=0; b < FCT->nbr_band(s); b++) 
	       {
	          for (int i=0; i < FCT->size_band_nl(s,b); i++)
		  for (int j=0; j < FCT->size_band_nc(s,b); j++)
		       (*FCT)(s,b,i,j) = MR_Transf(IndBand, i, j);
	          IndBand ++;
	       }
               FCT->cur_recons(Image);
	    }
	    break;
      case TO_PYR_MEYER:
      case TO_PYR_MEYER_ISOTROP:
           {
	      MEYER_WT *WT = MR_Transf.get_meyer_wt ();
	      WT->recons(MR_Transf.tabband(), Image, True);
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
	     WT.recons(MR_Transf.tabband(), Image,Nbr_Plan);
           }
  	  break; 
       case TM_MIN_MAX:
          {
 	      ortho_trans_to_ima(MR_Transf, Image);
	      G_Minmax C_GTM;
              Ortho_2D_WT WT2D(C_GTM);
              WT2D.recons(Image, MR_Transf.nbr_scale());
            }  
           break;
       case TO_FEAUVEAU:
 	    Buffer.alloc(MR_Transf.size_ima_nl(),MR_Transf.size_ima_nc(),"aux");
	    ortho_trans_to_ima(MR_Transf,  Buffer);
            Nbr_Iter = MR_Transf.nbr_scale()-1;
            recons_feauveau(Buffer, Image, Nbr_Iter);
           break;
       case TO_HAAR:
          haar_2d_reconstruct (MR_Transf, Image);
          break;
       default:
          fprintf (stderr, "Error (proc. mr_recons):Unknown transform\n");
          exit (-1);
          break;
    }
} 
  
/****************************************************************************/

void MultiResol::rec_adjoint (Ifloat &Image, Bool UseLastScale, type_border Bord)
{ 
    filter_bank_alloc();
    if (Type_Transform == TO_PAVE_BSPLINE)
    {
       ATROUS_2D_WT AWT;
       AWT.Bord = Bord;
       AWT.ModifiedAWT = True;
       // AWT.AdjointRec=True;
       Image.resize( TabBand[0].nl(), TabBand[0].nc());
       AWT.recons(TabBand, Image, Nbr_Plan, UseLastScale);
    }
    else mr_adjoint_rec (*this, Image, UseLastScale, Bord);
}

/****************************************************************************/

void MultiResol::recons (Ifloat &Image, type_border Bord)
{
    filter_bank_alloc();
    if (Type_Transform == TO_PAVE_BSPLINE)
    {
       ATROUS_2D_WT AWT;
       AWT.Bord = Bord;
       AWT.ModifiedAWT = ModifiedATWT;
       Image.resize( TabBand[0].nl(), TabBand[0].nc());
       AWT.recons(TabBand, Image, Nbr_Plan);
    }
    else  mr_recons (*this, Image, Bord, EdgeLineTransform);
}

/****************************************************************************/

