/******************************************************************************
**                   Copyright (C) 1996 by CEA
*******************************************************************************
**
**    UNIT
**
**    Version: 18:28:57
**
**    Author: Jean-Luc Starck
**
**    Date:  97/05/06
**    
**    File:  MR_NoiseModel.cc
**
*******************************************************************************
**
**    DESCRIPTION  noise modelling
**    -----------  
**                 
*******************************************************************************
**
** MRNoiseModel::MRNoiseModel(type_noise TNoise, int Nl_Imag, 
**             int Nc_Imag, int ScaleNumber, type_transform Trans)
**
** Class builder
** type_noise TNoise is the type of noise in the data
** Nl_Imag,Nc_Imag size of the image
** ScaleNumber: number of scales of the transform
** type_transform Trans: type of multiresoltuion transform
**
******************************************************************************/


#include "MR_Obj.h"
#include "IM_Noise.h"
#include "MR_Noise.h"
#include "MR_NoiseModel.h"
#include "MR_Sigma.h"
#include "MR_Abaque.h"
#include "MR_Psupport.h"
#include "NR.h"
#include "MR1D_Obj.h"
#include "FFTN_2D.h"
#include "Mr_FewEvent2d.h"

#define WRITE_DATA 0

/************************************************************************/

Bool one_level_per_pos_2d(type_noise TNoise)
// return True if we need only one level per position
// in the multiresolution space
{
    Bool ValRet= False;

    if ((TNoise == NOISE_EVENT_POISSON) || (TNoise == NOISE_NON_UNI_ADD)
         || (TNoise == NOISE_NON_UNI_MULT)  || (TNoise == NOISE_UNDEFINED))  ValRet = True;

    return ValRet;
}

/************************************************************************/

void MRNoiseModel::init_param()
{
   TabLevel = NULL;
   TabSupport = NULL;
   NbrScale = 0;   
   NbrBand = 0;
   Nl = 0;
   Nc = 0;
   Size = 0;
   //CEventPois = NULL;
   CFewEventPoisson2d = NULL;
   CFewEvent2d = NULL;
   CSpeckle = NULL;
   CorrelNoiseMap = NULL;
   Transform=T_UNDEFINED;
   Set_Transform=S_UNDEFINED;
   SigmaApprox = False;
   FilterBank = NULL;
   TypeNorm = DEF_SB_NORM;
   NbrUndecimatedScale = -1;
   GradientAnalysis = False;
   NoCompDistrib = False;
   PoissonFisz = False;
   MadEstimFromCenter = False;
   TypeThreshold = DEF_THRESHOLD;
   U_Filter = DEF_UNDER_FILTER;
   mOldPoisson = False;
   mWriteThreshold = False;
}

/************************************************************************/

void MRNoiseModel::alloc(type_noise TNoise, int Nl_Imag, 
                      int Nc_Imag, int ScaleNumber, type_transform Trans, 
                      FilterAnaSynt *FAS, sb_type_norm Norm, int NbrUndec, int FCT_NDir)
{
    int s, Nl_s=Nl_Imag, Nc_s=Nc_Imag;
    Nl = Nl_Imag;
    Nc = Nc_Imag;
    NbrScale = ScaleNumber;
    TypeNoise = TNoise;
    Transform = Trans;
    Set_Transform = SetTransform(Trans);
    NbrBand = NbrScale;
    FilterBank = FAS;
    TypeNorm = Norm;
    NbrUndecimatedScale = (NbrUndec >= 0) ? NbrUndec: NbrScale;

    if (Set_Transform == TRANSF_DIADIC_MALLAT) NbrBand = NbrScale*2-1;
    else if ((Set_Transform == TRANSF_UNDECIMATED_MALLAT)|| (Set_Transform == TRANSF_MALLAT))
            NbrBand = (NbrScale-1)*3+1;

   if (Set_Transform != TRANSF_UNDECIMATED_MALLAT)
   {
       TabNl.alloc(NbrBand);
       TabNc.alloc(NbrBand);
       TabPos.alloc(NbrBand);
       TabBandScale.alloc(NbrBand);
   }	
       
    switch (Set_Transform)
    {
        case TRANSF_UNDECIMATED_MALLAT:
              {
 	        if (Trans == TC_FCT)
		{
 		   NbrBand = fct_real_get_band_size(NbrScale, Nl, Nc, FCT_NDir, TabNl,  TabNc);        
                   // cout << NbrScale << " " << NbrScale << " " <<   NbrBand << " " <<  FCT_NDir << endl;
                    TabPos.alloc(NbrBand);
	           TabBandScale.alloc(NbrBand);
                   Size=0;
		   for (s = 0; s < NbrBand-1; s++) 
		   {
		      TabPos(s) = Size;
		      Size += TabNl(s)*TabNc(s);
		   }
 		}
		else if  (Trans == TO_LC) 
		{
 		    NbrScale = get_nbr_scale(Nc_Imag);
		    int NbrBandPerResol = get_nbr_scale(Nl_Imag);
		    NbrBand = NbrScale * NbrBandPerResol;
		    cout << NbrScale << " " << NbrBandPerResol << " " <<   NbrBand << endl;
                    TabNl.alloc(NbrBand);
                    TabNc.alloc(NbrBand);
                    TabPos.alloc(NbrBand);
	            TabBandScale.alloc(NbrBand);
 		    Size=0;
                    for (int i=0; i < NbrBand; i++) 
	            {
			    TabNl(i) =Nl_s;
			    TabNc(i) = Nc_s;
			    TabPos(i) = Size;
			    Size += (Nl_s*Nc_s);
		    }
		 }
		 else
		 {
                   Size=0;
                   int NbrBandPerResol = 3;
		   TabNl.alloc(NbrBand);
                   TabNc.alloc(NbrBand);
                   TabPos.alloc(NbrBand);
	           TabBandScale.alloc(NbrBand);
                   for (s = 0; s < NbrBand-1; s+=NbrBandPerResol)
                   {
                      if (s/NbrBandPerResol >= NbrUndecimatedScale)
                      { 
                       TabNl(s) = (Nl_s+1)/2;
                       TabNc(s) = Nc_s/2;  
                       TabPos(s) = Size;
                       Size +=  TabNl(s)*TabNc(s);
                       TabNl(s+1) = Nl_s/2;
                       TabNc(s+1) = (Nc_s+1)/2;
                       TabPos(s+1) = Size;
                       Size +=  TabNl(s+1)*TabNc(s+1);
                       TabNl(s+2) = Nl_s/2;
                       TabNc(s+2) = Nc_s/2;
                       TabPos(s+2) = Size;
                       Size +=  TabNl(s+2)*TabNc(s+2);
                       Nl_s = (Nl_s+1)/2;
                       Nc_s = (Nc_s+1)/2;
                      }
                      else
                      {
                        TabNl(s) = TabNl(s+1) =  TabNl(s+2) = Nl_s;
                        TabNc(s) = TabNc(s+1) =  TabNc(s+2) = Nc_s;
                        TabPos(s) = Size;
                        TabPos(s+1) = Size+Nl_s*Nc_s;
                        TabPos(s+2) = Size+2*(Nl_s*Nc_s);
                        Size += 3*(Nl_s*Nc_s);
		      }
                    }
                 }
 	      }
             //for (s = 0; s < NbrBand-1; s++)
             //    cout << "Band " << s+1 << " " << TabNl(s) << " " << TabNc(s) << endl;
              break;
        case TRANSF_DIADIC_MALLAT:
  	case TRANSF_PAVE: 
               Size = Nl*Nc*(NbrBand-1);
               for (s = 0; s < NbrBand-1; s++)
               {
                   TabNl(s) = Nl;
                   TabNc(s) = Nc;
                   TabPos(s) = s*Nl*Nc;
               }
               break;
        case TRANSF_PYR:
 	case TRANSF_SEMIPYR:
               Size=0;
               TabNl(0) = Nl_s;
               TabNc(0) = Nc_s;
               TabPos(0) = 0;
               for (s = 0; s < NbrScale-1; s++)
               {
                   Size += Nl_s*Nc_s;
                   TabPos(s+1) = TabPos(s)+Nl_s*Nc_s;
		   if ((s != 0) || (Set_Transform != TRANSF_SEMIPYR))
	           {
                      Nl_s = (Nl_s+1)/2;
                      Nc_s = (Nc_s+1)/2;
		   }
                   TabNl(s+1) = Nl_s;
                   TabNc(s+1) = Nc_s;
               }
               break;
        case TRANSF_MALLAT:  
               NbrBand = 3*(NbrScale -1)+1;
               TabPos(0) = 0;
	       // cout << Nl_s << " " << Nc_s << endl;
	       Size = 0;
               for (s = 0; s < NbrBand-1; s+=3)
               {
                  TabNl(s) = (Nl_s+1)/2;
                  TabNc(s) = Nc_s/2;
		  TabPos(s) = Size;
		  Size += TabNl(s)*TabNc(s);
                  TabNl(s+1) = Nl_s/2;
                  TabNc(s+1) = (Nc_s+1)/2;
		  TabPos(s+1) = Size;
		  Size += TabNl(s+1)*TabNc(s+1);
                  TabNl(s+2) = Nl_s/2;
                  TabNc(s+2) = Nc_s/2;
		  TabPos(s+2) = Size;
		  Size += TabNl(s+2)*TabNc(s+2);
                  Nl_s = (Nl_s+1)/2;
                  Nc_s = (Nc_s+1)/2;
		  // cout << "Band " << s << ": " << TabNl(s) << " " << TabNc(s) << " " << TabPos(s) << endl;
 		  // cout << "Band " << s+1 << ": " << TabNl(s+1) << " " << TabNc(s+1) << " " << TabPos(s+1) << endl;
		  // cout << "Band " << s+2 << ": " << TabNl(s+2) << " " << TabNc(s+2) << " " << TabPos(s+2) << endl;
              }
               s= NbrBand-1;
               TabNl(s) = Nl_s;
               TabNc(s) = Nc_s;
	       TabPos(s) = Size;
               Size = Nl*Nc;
	       // for (s = 0; s < NbrBand-1; s++) cout << "bb " << s << ": " << TabNl(s) << " " << TabNc(s) << " " << TabPos(s) << endl;
               break;
        case TRANSF_FEAUVEAU:               
               NbrBand = 2*( NbrScale-1)+1;
               TabPos(0) = 0;
               for (s = 0; s < NbrBand-1; s+=2)
               {
                  int Nlw = Nl_s;
                  int Ncw = Nc_s/2;
                  TabNl(s) = Nlw;
                  TabNc(s) = Ncw;
                  TabPos(s+1) = TabPos(s) + Nlw*Ncw;
                  Nlw /= 2;
                  Ncw = (Nc_s+1)/2;
                  TabNl(s+1) = Nlw;
                  TabNc(s+1) = Ncw;
                  TabPos(s+2) = TabPos(s+1) + Nlw*Ncw;
                  Nl_s = (Nl_s+1)/2;
                  Nc_s = (Nc_s+1)/2;
               }
               s= NbrBand-1;
               TabNl(s) = Nl_s;
               TabNc(s) = Nc_s;
               // Size = TabPos(s);
               Size=Nl*Nc;
               break;
        default: cerr << "Not implemented" << endl;
                 exit (0);
                 break;
    }

    // if (one_level_per_pos_2d(TNoise) == True) TabLevel = new float[Size];
     if (one_level_per_pos_2d(TNoise) == True) 
     {
         TabLevel = (float *) alloc_buffer((size_t) (Size*sizeof(float)));
     }
    else TabLevel = new float[NbrBand];

    // TabSupport = new unsigned char [Size];
    TabSupport = (unsigned char *) alloc_buffer((size_t) (Size*sizeof(char)));
    for (s = 0; s < Size; s++) TabSupport[s] = VAL_SupNull;
    // cout << "Size = " << Size << endl;

    CCD_Gain = 1.; 
    CCD_ReadOutSigma=0.; 
    CCD_ReadOutMean=0.; 
    SigmaNoise=0.; 
    for (s = 0; s < NbrBand-1; s++) 
    {
        details which_detail;
        band2scale(s, Transform, NbrBand,  TabBandScale(s), which_detail);
	NSigma[s] = DEF_NSIGMA_DETECT;
	TabEps[s] = DEFAULT_EPSILON;
	if (TabBandScale(s) == 0) NSigma[s] += 1;
	// cout << " band " << s+1 << " scale = " << TabBandScale(s) << " NSigma " << NSigma[s] << endl;
    }
    BadPixelVal=0.;
    BadPixel=False;
    MinEvent=False;
    SupIsol=False;
    FirstDectectScale = DEF_FIRST_DETECT_SCALE;
    OnlyPositivDetect=False;
    DilateSupport=False;
    GetEgde=False;
    MinEventNumber=DEF_MINEVENT_NUMBER;
    SigmaDetectionMethod = DEF_SIGMA_METHOD;
    NiterSigmaClip = 1;
    SizeBlockSigmaNoise = DEFAULT_SIZE_BLOCK_SIG;
    TransImag=False;
    NewStatNoise = TNoise;
    UseRmsMap = False;
    GetRmsCoeffGauss = True;
    //CEventPois = NULL;
    CFewEventPoisson2d = NULL;
    CFewEvent2d = NULL;
    CSpeckle = NULL;
    CorrelNoiseMap = NULL;
    Border = DEFAULT_BORDER;
    SigmaApprox = False;
    switch (TypeNoise)
    {
        case NOISE_GAUSSIAN: 
        case NOISE_NON_UNI_ADD:
        case NOISE_UNDEFINED:
               TransImag=False;    
               break;
        case NOISE_POISSON:
        case NOISE_GAUSS_POISSON:
        case NOISE_MULTI: 
        case NOISE_NON_UNI_MULT:
               TransImag=True;
               SigmaNoise=1.;
               break;
        case NOISE_EVENT_POISSON:
               TransImag=True;
               MinEvent=True;
        //       CEventPois = new StatEventPoisson;
	       if( mOldPoisson ) {
	          CFewEventPoisson2d = new FewEventPoisson;
                  if (!NoCompDistrib) {  // call from mr_pfilter (Abaque exist)
	//          CEventPois->compute_distrib();
		     CFewEventPoisson2d->compute_distribution();
		  }
	       } else{
	          CFewEvent2d = new FewEvent;
		  CFewEvent2d->write_threshold( mWriteThreshold );
                  if (!NoCompDistrib) 
		     CFewEvent2d->compute_distribution();
	       } 
               break;
        case NOISE_UNI_UNDEFINED: 
               NiterSigmaClip = 3;
               break;
        case NOISE_CORREL:
	       break;
        case NOISE_SPECKLE:
	       TransImag=True;
	       CSpeckle = new StatRayleigh(TSPECKLE, NbrScale, 1, Transform);
	       break;
    }
}

/************************************************************************/

MRNoiseModel::MRNoiseModel()
{
 (*this).init_param();
}

/************************************************************************/

MRNoiseModel::MRNoiseModel(type_noise TNoise, int Nl_Imag, 
                      int Nc_Imag, int ScaleNumber, type_transform Trans)
{
   (*this).init_param();
   (*this).alloc(TNoise, Nl_Imag, Nc_Imag, ScaleNumber, Trans);
}

/****************************************************************************/

MRNoiseModel::MRNoiseModel (type_noise TNoise,  MultiResol &MR_Data)
{
    int Nl_s = MR_Data.size_ima_nl();
    int Nc_s = MR_Data.size_ima_nc();
    
   (*this).init_param();
   (*this).alloc(TNoise, Nl_s, Nc_s, MR_Data.nbr_scale(), 
                 MR_Data.Type_Transform, MR_Data.filter_bank(), 
                 MR_Data.TypeNorm, MR_Data.nbr_undec_scale());    

    //MRNoiseModel(TNoise, Nl_s, Nc_s, 
    //             MR_Data.nbr_scale(), MR_Data.Type_Transform);
}

/****************************************************************************/

void MRNoiseModel::pos_mrcoeff(int NumCoef, int &b, int &i, int &j)
{
    int PosCoef = NumCoef;
    b=0;
    while (PosCoef > TabNc(b)*TabNl(b))   
    {
       PosCoef -= TabNc(b)*TabNl(b);
       b++;
    } 
    i = PosCoef /  TabNc(b);
    j = PosCoef - i *  TabNc(b);
}

/****************************************************************************/

Bool MRNoiseModel::operator() (int b,int i,int j)  
{
   Bool ValRet=False; 
   int Ind = index(b,i,j);
//   cout << b << " " << i << " " << j << " ==> " << Ind << endl;
   if ((TabSupport[Ind] > 0) && (TabSupport[Ind] <= VAL_SupLastOK)) ValRet = True;
   return ValRet;
}

// ***************

Bool MRNoiseModel::operator() (int s,int i,int j, details which_detail)
{
   int b = scale2band(s,Transform, NbrBand,which_detail);
   return (*this)(b,i,j);
}

// ***************

Bool MRNoiseModel::operator() (int NumCoef)
{
   int b,i,j;
   pos_mrcoeff(NumCoef, b, i, j);
   return  (*this)(b,i,j);
}

/****************************************************************************/

float & MRNoiseModel::sigma(int b, int i, int j)
{
   int Ind;    
 
   if (one_level_per_pos_2d(TypeNoise) == True) Ind = index(b,i,j);
   else Ind = b;
   return TabLevel[Ind];
}

// ***************

float & MRNoiseModel::sigma(int s,int i,int j, details which_detail)  
{
   int b = scale2band(s,Transform, NbrBand,which_detail);
   return sigma(b,i,j);
}

// ***************

float & MRNoiseModel::sigma(int NumCoef)  
{
   if (one_level_per_pos_2d(TypeNoise) == True) {
     int i,j,b;
     pos_mrcoeff(NumCoef, b, i, j);
     return sigma(b,i,j);
   } else {
	 return TabLevel[NumCoef];
   }
}

/****************************************************************************/

float MRNoiseModel::nsigma(int b)  
{
   return NSigma[b];
}

/****************************************************************************/

unsigned char & MRNoiseModel::support(int b,int i,int j)
{
    int Ind = index(b,i,j);
    // cout << b << " " << i << " " << j << " ==> " << Ind << endl;
    return TabSupport[Ind];
}
// ***************
unsigned char & MRNoiseModel::support(int s,int i,int j, details which_detail)  
{
   int b = scale2band(s,Transform, NbrBand,which_detail);
   return  support(b,i,j);
}
// ***************
unsigned char & MRNoiseModel::support(int NumCoef)  
{
   int i,j,b;
   pos_mrcoeff(NumCoef, b, i, j);
   return support(b,i,j);
}
 

/****************************************************************************/

Bool MRNoiseModel::signif (float Val, int b, int i, int j, float LevelMin, float LevelMax)
{
   Bool ValRet = False;
   float Level;
   Level = sigma(b,i,j)*NSigma[b];
   if (OnlyPositivDetect == True)   
   {
       if (Val > LevelMax) ValRet = True;
   }
   else if ((Val > LevelMax) || (Val < LevelMin)) ValRet = True;
   
   if ((TabBandScale(b) < FirstDectectScale) && (ValRet == True)) ValRet = False;
   return ValRet;
}

/****************************************************************************/

Bool MRNoiseModel::signif (float Val, int b, int i, int j)
{
   Bool ValRet = False;
   float Level;
   Level = sigma(b,i,j)*NSigma[b];
   if (OnlyPositivDetect == True)   
   {
       if (Val > Level) ValRet = True;
   }
   else  if (ABS(Val) > Level) ValRet = True;
   
   if ((TabBandScale(b) < FirstDectectScale) && (ValRet == True)) ValRet = False;
   return ValRet;
}

Bool MRNoiseModel::signif (float Val, int b, int i, int j, fltarray & TNsigma)
{
   Bool ValRet = False;
   float Level;
   Level = sigma(b,i,j) * TNsigma(b);
   if (OnlyPositivDetect == True)   
   {
       if (Val > Level) ValRet = True;
   }
   else  if (ABS(Val) > Level) ValRet = True;
   
   if ((TabBandScale(b) < FirstDectectScale) && (ValRet == True)) ValRet = False;
   return ValRet;
}

// ***************

Bool MRNoiseModel::signif (float Val, int s, int i, int j, details which_detail)
{
   int b = scale2band(s,Transform, NbrBand,which_detail);
   return signif(Val, b,i,j);
}

/****************************************************************************/

float MRNoiseModel::prob(float Val, int b, int i, int j)  
{
    float Sig,P=0.;

    switch (TypeNoise)
    {
       case NOISE_EVENT_POISSON:
         {
            int k,l,Win = (int) (pow((double)2., (double)(b+2)) + 0.5);
            int Nevent = 0;
            for (k =i-Win; k <= i+Win; k++)
            for (l =j-Win; l <= j+Win; l++)
               Nevent += Event_Image(k, l, Border);
            if (NoCompDistrib) {
	       cout << "Error: histogram have to be computed first ..." << endl;
	       exit(-1);
	    }
	    //P = CEventPois->a_trou_prob(Val, Nevent, b);
	    if( mOldPoisson )
	       P = CFewEventPoisson2d->a_trou_prob( Val, Nevent, b );
	    else
	       P = CFewEvent2d->a_trou_prob( Val, Nevent, b );  
          }
	 break;
	case NOISE_CORREL:
	   P = CorrelNoiseMap->prob(b, Val);
	   break;
       case NOISE_SPECKLE:
           P = CSpeckle->prob(b, Val); 
         break;
       default:
         Sig = sigma(b,i,j);
         P = exp(-Val*Val / (2*Sig*Sig) ) / sqrt(2.*PI);
 	break;
    }	
 
  
    return P;
}

// ***************

float MRNoiseModel::prob (float Val, int s, int i, int j, details which_detail)
{
   int b = scale2band(s,Transform, NbrBand,which_detail);
   return prob(Val,b,i,j);    
}


/****************************************************************************/

void MRNoiseModel::prob (MultiResol &MR_Data, Bool Complement)
{
    float Sig,Val;
    int s,i,j;


    switch (TypeNoise)
    {
       case NOISE_EVENT_POISSON:
        {
          Iint EventCount(Nl,Nc,"ImagCount");
          for (s = 0; s < MR_Data.nbr_band()-1; s++)
          {
             event_one_scale(Event_Image, s, EventCount, MR_Data.Border);
             for (i=0; i< MR_Data.size_band_nl(s); i++)
             for (j=0; j< MR_Data.size_band_nc(s); j++)
             {
	     if (NoCompDistrib) {
	        cout << "Error: histogram have to be computed first ... " << endl;
	        exit(-1);
	     }
             //MR_Data(s,i,j) = CEventPois->a_trou_prob(MR_Data(s,i,j), EventCount(i,j), s);
             if( mOldPoisson )
	        MR_Data(s,i,j) = CFewEventPoisson2d->a_trou_prob( MR_Data(s,i,j),
		                                                  EventCount(i,j),
								  s );
             else
                MR_Data(s,i,j) = CFewEvent2d->a_trou_prob( MR_Data(s,i,j), 
		                                           EventCount(i,j), s );
	     
	     if (Complement == True) MR_Data(s,i,j) = 1. - MR_Data(s,i,j) ;
             }
           } 
         }
	 break; 
	case  NOISE_CORREL:
	  for (s = 0; s < MR_Data.nbr_band()-1; s++)
          for (i=0; i< MR_Data.size_band_nl(s); i++)
          for (j=0; j< MR_Data.size_band_nc(s); j++)
          {
	     Val = MR_Data(s,i,j);
	     MR_Data(s,i,j) = CorrelNoiseMap->prob(s, Val);
 	  } 
          break;    
       case NOISE_SPECKLE:
          for (s = 0; s < MR_Data.nbr_band()-1; s++)
          for (i=0; i< MR_Data.size_band_nl(s); i++)
          for (j=0; j< MR_Data.size_band_nc(s); j++)
          {
	     Val = MR_Data(s,i,j);
	     MR_Data(s,i,j) = CSpeckle->prob(s, Val);
 	  } 
          break;
       default:
          for (s = 0; s < MR_Data.nbr_band()-1; s++)
          for (i=0; i< MR_Data.size_band_nl(s); i++)
          for (j=0; j< MR_Data.size_band_nc(s); j++)
          {
             Val = MR_Data(s,i,j);
             Sig = sigma(s,i,j);
	     if (Sig > FLOAT_EPSILON)
                Val =  1. / sqrt(2.*PI)*exp(-Val*Val / (2*Sig*Sig));
             else Val = 0.;
	     if (Complement == True) Val = 1. - Val;
             MR_Data(s,i,j) = Val;
          }
	  break;
    }
}

/****************************************************************************/

void MRNoiseModel::prob_noise (MultiResol &MR_Data, Bool Complement)
{
    int s,i,j;

    switch (TypeNoise)
    {
       case NOISE_EVENT_POISSON:
        {
          Iint EventCount(Nl,Nc,"ImagCount");
          for (s = 0; s < MR_Data.nbr_band()-1; s++)
          {
             event_one_scale(Event_Image, s, EventCount, MR_Data.Border);
             for (i=0; i< MR_Data.size_band_nl(s); i++)
             for (j=0; j< MR_Data.size_band_nc(s); j++)
             {
	         if (NoCompDistrib) {
	           cout << "Error: histogram have to be computed first ..." << endl;
	           exit(-1);
	         }	         
		 //float P = CEventPois->a_trou_repartition(MR_Data(s,i,j), 
		 //                                   EventCount(i,j) , s);
		 double P1, P2;
//std::cout << "s:" << s << ", [x:" << i << ",y:" << j << "]" << std::endl;
		 if( mOldPoisson ) {
                    P1 = CFewEventPoisson2d->a_trou_repartition( MR_Data(s,i,j), 
		                                                EventCount(i,j),
							        s);
		    if (MR_Data(s,i,j) > 0) P1 = 1. - P1;
                    if (Complement == True) MR_Data(s,i,j) = 1. - P1;
		    else MR_Data(s,i,j) = P1;
                 } else {
		    P1 = CFewEvent2d->a_trou_repartition( MR_Data(s,i,j), 
		                                         EventCount(i,j), s);
		    P2 = CFewEvent2d->a_trou_repartition( MR_Data(s,i,j), 
		                                         EventCount(i,j), s, True);
//std::cout << "MRNoiseModel::prob_noise : P1 = " << P1 << ", 1-P1 = " << 1. - P1
//          << ", P2 = " << P2 << std::endl;
                    if (Complement == True) MR_Data(s,i,j) = 1. - P2;
		    else MR_Data(s,i,j) = P2;
		 }
             }
           } 
         }
	 break;      
       default:
          for (s = 0; s < MR_Data.nbr_band()-1; s++)
	  for (i=0; i< MR_Data.size_band_nl(s); i++)
          for (j=0; j< MR_Data.size_band_nc(s); j++)
          {
	     MR_Data(s,i,j) = prob_noise(MR_Data(s,i,j), s, i, j);
	     if (Complement == True) MR_Data(s,i,j) = 1. - MR_Data(s,i,j);
 	  }
          break;
    }
}

/****************************************************************************/
 
float MRNoiseModel::prob_noise(float Val, int b, int i, int j)  
{
    float Sig,P=0.;
    double Vn=0.;
    
    switch (TypeNoise)
    {
       case NOISE_EVENT_POISSON:
        {
	   if (NoCompDistrib) {
	        cout << "Error: histogram have to be computed first ..." << endl;
	        exit(-1);
	   }
	   int k,l,Win = (int) (pow((double)2., (double)(b+2)) + 0.5);
           int Nevent = 0;
           for (k =i-Win; k <= i+Win; k++)
           for (l =j-Win; l <= j+Win; l++)
              Nevent += Event_Image(k, l, Border);
           //P = CEventPois->a_trou_repartition(Val, Nevent, b);
	   if( mOldPoisson ) {
              P = CFewEventPoisson2d->a_trou_repartition(Val, Nevent, b);
	      if (Val > 0) P = 1. - P; 
	   } else {
	      //double P1 = CFewEvent2d->a_trou_repartition(Val, Nevent, b);
	      double P2 = CFewEvent2d->a_trou_repartition( Val, Nevent, b, True );
   // std::cout << "MRNoiseModel::prob_noise : P1 = " << P1 << ", 1-P1 = " << 1. - P1  << ", P2 = " << P2 << std::endl;
              P = P2;
	   }	  
        }  
	break;
       case NOISE_CORREL:
           P = CorrelNoiseMap->repartition(b, Val);
	   if (Val > 0) P = 1. - P;
           break;
       case NOISE_SPECKLE:
           P = CSpeckle-> repartition(b, Val);
	   if (Val > 0) P = 1. - P;
           break;
      default:    
        Sig = sigma(b,i,j);
        if (ABS(Val) < FLOAT_EPSILON) P = 1.;
	else 
	{
	   if (Sig < FLOAT_EPSILON) P = 0;
	   else 
	   {
	      Vn = ABS(Val) / (sqrt(2.)*Sig);
	      if (Vn > 3.5) P = 0;
	      else P = (float) erfc (Vn);
	   }
	}
	break;
    }
 
    return P;
    
}

/****************************************************************************/
double 
MRNoiseModel:: prob_signal_few_event( float Val, int b, int i, int j )  
{

   double P = 0 ;    
   if( NoCompDistrib ) {
      std::cout << "Error: histogram have to be computed first ..." << std::endl ;
      exit( -1 ) ;
   }
   int k,l,Win = (int) ( pow((double)2., (double)( b+2 )) + 0.5) ;
   int Nevent = 0 ;
   for( k =i-Win; k <= i+Win; k++ )
   for( l =j-Win; l <= j+Win; l++ )
      Nevent += Event_Image(k, l, Border) ;
   if( mOldPoisson ) {
      P = CFewEventPoisson2d->a_trou_repartition( Val, Nevent, b ) ;
      if( Val > 0 ) P = 1. - P ; 
   } else {
      //double P1 = CFewEvent2d->a_trou_repartition(Val, Nevent, b);
      P = CFewEvent2d->a_trou_repartition( Val, Nevent, b, True ) ;
      //std::cout << "MRNoiseModel::prob_signal_few_event : P = " << std::endl ;
   }  
   if( kill_coef( b, i, j, Val, False ) == True ) P = 1. ; 
   return P;
}

 
// ***************


float MRNoiseModel::prob_noise(float Val, int s, int i,
                                int j, details which_detail)
{
   int b = scale2band(s,Transform, NbrBand,which_detail);
   return prob_noise(Val,b,i,j);    
}

/****************************************************************************/

float MRNoiseModel::prob_signal(float Val, int b, int i, int j)  
{
   return (1. - prob_noise(Val,b,i,j));  
}

/****************************************************************************/

float MRNoiseModel::prob_signal(float Val, int s, int i,
                                int j, details which_detail)
{
   int b = scale2band(s,Transform, NbrBand,which_detail);
   return  (1. - prob_noise(Val,b,i,j));   
}
 
/****************************************************************************/

float MRNoiseModel::val_transform(float Val)
{
    float ValRet = Val;
    switch (TypeNoise)
    {
        case NOISE_POISSON:
        case NOISE_GAUSS_POISSON:
               ValRet = Val*CCD_Gain + 3. / 8. * CCD_Gain*CCD_Gain + 
                   CCD_ReadOutSigma*CCD_ReadOutSigma - CCD_Gain*CCD_ReadOutMean;
               if (ValRet < 0.) ValRet=0;
               else  ValRet = 2. / CCD_Gain *sqrt(ValRet);
               break;
        case NOISE_GAUSSIAN: 
        case NOISE_NON_UNI_ADD:
        case NOISE_UNDEFINED:
        case NOISE_UNI_UNDEFINED:
	case NOISE_CORREL:
               break;
        case NOISE_MULTI: 
        case NOISE_NON_UNI_MULT:
	case NOISE_SPECKLE:
               if (Val > 0) ValRet = log(Val+1.);
               else ValRet = 0.;
               break;
        case NOISE_EVENT_POISSON:
               break;	
    }
    return ValRet;

}

/****************************************************************************/

float MRNoiseModel::val_invtransform(float Val)
{
    float ValRet = Val;
    switch (TypeNoise)
    {
        case NOISE_POISSON:
        case NOISE_GAUSS_POISSON:
              ValRet = Val*Val/4.*CCD_Gain - ( 3./8.*CCD_Gain*CCD_Gain +
                          CCD_ReadOutSigma*CCD_ReadOutSigma -
                          CCD_Gain*CCD_ReadOutMean ) / CCD_Gain;
               break;
        case NOISE_GAUSSIAN: 
        case NOISE_NON_UNI_ADD:
        case NOISE_UNDEFINED:
        case NOISE_UNI_UNDEFINED:
	case NOISE_CORREL:
               break;
        case NOISE_MULTI: 
        case NOISE_NON_UNI_MULT:
	case NOISE_SPECKLE:
               if (Val > 0) ValRet = exp(Val) - 1.;
               else ValRet = 0.;
               break;
        case NOISE_EVENT_POISSON:
               break;
    }
    return ValRet;
}

/****************************************************************************/

void MRNoiseModel::im_transform(Ifloat &Image)
{
    switch (TypeNoise)
    {
        case NOISE_POISSON:
	       if (PoissonFisz == False) noise_poisson_transform (Image, Image);
               else fisz2d_trans(Image);
	       SigmaNoise = 1.;
               break;
        case NOISE_GAUSS_POISSON:
               noise_poisson_transform (Image, Image);
               SigmaNoise = 1.;
               break;
        case NOISE_GAUSSIAN: 
        case NOISE_NON_UNI_ADD:
        case NOISE_UNDEFINED:
        case NOISE_UNI_UNDEFINED:
	case NOISE_CORREL:
               break;
        case NOISE_MULTI: 
        case NOISE_NON_UNI_MULT:
	case NOISE_SPECKLE:
               noise_log_transform (Image, Image);
               break;
        case NOISE_EVENT_POISSON:
               // event_poisson_transform (Image, Event_Image);
	       building_imag_imag(Image, Event_Image);
               break;
    }
}


/****************************************************************************/

void MRNoiseModel::im_invtransform(Ifloat &Image)
{
    switch (TypeNoise)
    {
        case NOISE_POISSON:
	       if (PoissonFisz == False) noise_inverse_poisson_transform (Image, Image);
               else fisz2d_inv(Image);
	       break;
        case NOISE_GAUSS_POISSON:
               noise_inverse_poisson_transform (Image, Image);
               break;
        case NOISE_GAUSSIAN: 
        case NOISE_NON_UNI_ADD:
        case NOISE_UNDEFINED:
        case NOISE_UNI_UNDEFINED:
	case NOISE_CORREL:
               break;
        case NOISE_MULTI: 
        case NOISE_NON_UNI_MULT:
	case NOISE_SPECKLE:
               noise_inv_log_transform (Image, Image);
               break;
        case NOISE_EVENT_POISSON:
               break;
    }
}

/****************************************************************************/

void MRNoiseModel::get_sigma_from_rms_map(Ifloat &Ima_Sigma)
{
  int RmsNl = Ima_Sigma.nl();
  int RmsNc = Ima_Sigma.nc();
  Ifloat Ptr_Rms;
  Ifloat Dirac;
  int i,j,s;
  int Ind;
  FFTN_2D FFT;
  FFT.CenterZeroFreq = True;
 
  Dirac.alloc(RmsNl,RmsNc,"Dirac");
  Ptr_Rms.alloc(RmsNl,RmsNc,"Ptr_Rms");
  Dirac(RmsNl/2-1,RmsNc/2-1)=1.;

  // First we compute the square image of Rms
  for(i=0; i < RmsNl; i++)
  for(j=0; j < RmsNc; j++) Ptr_Rms(i,j) = Ima_Sigma(i,j)*Ima_Sigma(i,j);
  
  // Then we transform the dirac image :
  MultiResol MR_Dirac(RmsNl, RmsNc, NbrScale, Transform, "Dirac");
  MR_Dirac.transform(Dirac);
 
  // Then we compute the square of the Dirac transform ...
  for (s=0 ; s < MR_Dirac.nbr_band()-1; s++)
  for (i=0 ; i < RmsNl; i++)
  for (j=0 ; j < RmsNc; j++)  MR_Dirac(s,i,j) *= MR_Dirac(s,i,j);
 
  // Then we convolve the RMS by the Dirac (at each scale)
  // The problem is now to take the sqrt of the result image. The FFT may 
  // produce some negative values. Since the result is lower for high number
  // of scales, we cannot use a test based on a comparision of FLOAT_EPSILON
  // or some power of it. We look simply if there was something at the place
  // on the original map ...

  for (s=0 ; s < MR_Dirac.nbr_band()-1 ; s++)
  {
     FFT.convolve(MR_Dirac.band(s), Ptr_Rms);
     for (i=0 ; i < RmsNl ; i++)
     for (j=0 ; j < RmsNc ; j++)
     {
	Ind = index(s,i,j);
	if (MR_Dirac(s,i,j) <= 0) TabLevel[Ind] = 0;  
	else TabLevel[Ind] = sqrt(MR_Dirac(s,i,j));
     }
  }
}

/****************************************************************/

void MRNoiseModel::max_sigma_in_scale(Ifloat &Ima_Sigma)
{
   int s,i,j,Kb,Ke,Lb,Le,k,l;
   int Nls=Nl,Ncs=Nc,Nlp=Nl,Ncp=Nc;
   int Ind, Indp,Step,SizeWave=2;
   float Max;
   extern float mr_tab_noise(int s);

   // FIRST SCALE:
   // find the maximum noise standard deviation in a box
   // put the result in TabLevel
   for (i=0; i< Nl-1; i++)
   for (j=0; j< Nc-1; j++)
   {
           Max = Ima_Sigma(i,j);
           Kb = (i-SizeWave >= 0) ? i-SizeWave : 0;
           Ke = (i+SizeWave <= Nl-1) ? i+SizeWave : Nl-1;
           Lb = (j-SizeWave >= 0) ? j-SizeWave : 0;
           Le = (j+SizeWave <= Nc-1) ? j+SizeWave : Nc-1;

           for (k=Kb; k<= Ke; k++)
           for (l=Lb; l<= Le; l++)
                   if (Max < Ima_Sigma(k,l)) Max =Ima_Sigma(k,l);
           Ind = index(0,i,j);
           TabLevel[Ind] = Max*mr_tab_noise(0) / 0.972463;
    }

   Step=1;
   for (s=1; s < NbrBand-1; s++)
   {
       Nls = TabNl(s);
       Ncs = TabNc(s);
       Nlp = TabNl(s-1);
       Ncp = TabNc(s-1);
       if ((Set_Transform == TRANSF_PAVE) || 
           ((Set_Transform == TRANSF_SEMIPYR) && (s == 1)))
       {
           Step = 2*Step;
       }
       for (i=0; i< Nls-1; i++)
       for (j=0; j< Ncs-1; j++)
       {
           Max = 0.;
           if ((Set_Transform == TRANSF_PAVE) || 
             ((Set_Transform == TRANSF_SEMIPYR) && (s == 1)))
           {
               Kb = (i-Step >= 0) ? i-Step : i;
               Ke = (i+Step <= Nl-1) ? i+Step : i;
               Lb = (j-Step >= 0) ? j-Step : j;
               Le = (j+Step <= Nc-1) ? j+Step : j;
 
           }
           else
           {
               Kb = (2*i-2 >= 0) ? 2*i-2 : 0;
               Ke = (2*i+2 <= Nlp-1) ? 2*i+2 : Nlp-1;
               Lb = (2*j-2 >= 0) ? 2*j-2 : 0;
               Le = (2*j+2 <= Ncp-1) ? 2*j+2 : Ncp-1;
           }
           for (k=Kb; k<= Ke; k+=Step)
           for (l=Lb; l<= Le; l+=Step)
           {
              Indp = index(s-1,k,l);
              if (Max < TabLevel[Indp]) Max = TabLevel[Indp];
           }	   
           Ind = index(s,i,j);
           TabLevel[Ind] = Max*mr_tab_noise(s)/mr_tab_noise(s-1);
       }
   }
}

/****************************************************************************/

void MRNoiseModel::kill_isol(int s)
{
   int i,j;
   Bool PutZero;
   int Nls=Nl;
   int Ncs=Nl;

   // cout << "kill_isol" << endl;
   if (Set_Transform == TRANSF_PYR) 
   {
      Nls = TabNl(s);
      Ncs = TabNc(s);
   }
   if (Set_Transform == TRANSF_PAVE)
   {
       for (i = 1; i < Nl-1; i++)
       for (j = 1; j < Nc-1; j++)
       {
           PutZero = True;
           if (support(s,i,j) == VAL_SupOK)
           {
               if (support(s,i-1,j) == VAL_SupOK) PutZero = False;
               else if (support(s,i+1,j) == VAL_SupOK) PutZero = False;
               else if (support(s,i,j+1) == VAL_SupOK) PutZero = False;
               else if (support(s,i,j-1) == VAL_SupOK) PutZero = False;
               else if (support(s+1,i,j) == VAL_SupOK) PutZero = False;
               if (PutZero == True) support(s,i,j) = VAL_SupKill;
           }
       }
   }
   else  if ((Set_Transform == TRANSF_PYR) || (Set_Transform ==  TRANSF_SEMIPYR)
             || (Set_Transform == TRANSF_DIADIC_MALLAT)
	     || (Set_Transform == TRANSF_UNDECIMATED_MALLAT))
   {
       Nls = TabNl(s);
       Ncs = TabNc(s);
       for (i = 1; i < Nls-1; i++)
       for (j = 1; j < Ncs-1; j++)
       {
           PutZero = True;
           if (support(s,i,j) == VAL_SupOK)
           {
               if (support(s,i-1,j) == VAL_SupOK) PutZero = False;
               else if (support(s,i+1,j) == VAL_SupOK) PutZero = False;
               else if (support(s,i,j+1) == VAL_SupOK) PutZero = False;
               else if (support(s,i,j-1) == VAL_SupOK) PutZero = False;
               // else if (support(s+1,i/2,j/2) == VAL_SupOK) PutZero = False;
               if (PutZero == True) support(s,i,j) = VAL_SupKill;
           }
       }
   }
}

/****************************************************************************/

void MRNoiseModel::dilate_support(int s)
{
    int dim=1;
    int Ind;
    unsigned char *Buff;
    int i,j,k,l,Nls,Ncs;
    int kb,ke,lb,le;
    unsigned char Max;

    if ((Set_Transform == TRANSF_PYR) || (Set_Transform == TRANSF_PAVE)
       || (Set_Transform == TRANSF_SEMIPYR) 
       || (Set_Transform == TRANSF_DIADIC_MALLAT)
       || (Set_Transform == TRANSF_UNDECIMATED_MALLAT))
    {
        if (Set_Transform == TRANSF_PAVE) 
                          dim   = (int)(pow((double)2., (double) (s+1)) + 0.5);

        Nls = TabNl(s);
        Ncs = TabNc(s);
        Buff = new unsigned char [Nls*Ncs];
        for (i = 0; i < Nls*Ncs; i++) Buff[i] =0;

        for (i = 0; i < Nls; i++) 
        for (j = 0; j < Ncs; j++) 
        { 
            Max = ( (*this)(s,i,j) == True ) ? 1 : 0;
            kb = (i-dim >= 0) ? i-dim: 0;
            lb = (j-dim >= 0) ? j-dim: 0;
            ke = (i+dim < Nls) ? i+dim: Nls-1;
            le = (j+dim < Ncs) ? j+dim: Ncs-1;
            k=kb;
            while ( (Max == 0) && (k <= ke))
            {
                l = lb;
                while ( (Max == 0) && (l <= le))
                {
                        if ((*this)(s,k,l) == True) Max=1;
                        l++;
                }
                k++;
            }
            Ind = i*Ncs+j;
            Buff[Ind] = Max;
        }

        for (i = 0; i < Nls; i++) 
        for (j = 0; j < Ncs; j++) 
        { 
            Ind = i*Ncs+j;
            if ((support(s,i,j) != VAL_SupOK) && (Buff[Ind] == 1))
            {
                support(s,i,j) = VAL_SupDill;
            }
        }
        delete [] Buff;
    }
}

/****************************************************************************/

void MRNoiseModel::hierarchical_dilate_support()
{
    int dim=1;
    int i,j,k,l,Nls,Ncs;
    int kb,ke,lb,le;
    unsigned char Max;

    if (Set_Transform == TRANSF_PAVE)
    for (int s = NbrBand-3; s >= 0; s--) 
    {
       int NextScale = s+1;
       dim   = (int)(pow((double)2., (double) s) + 0.5);
       Nls = TabNl(s);
       Ncs = TabNc(s);
       for (i = 0; i < Nls; i++) 
       for (j = 0; j < Ncs; j++) 
       { 
           Max = ( (*this)(NextScale,i,j) == True ) ? 1 : 0;
           kb = (i-dim >= 0) ? i-dim: 0;
           lb = (j-dim >= 0) ? j-dim: 0;
           ke = (i+dim < Nls) ? i+dim: Nls-1;
           le = (j+dim < Ncs) ? j+dim: Ncs-1;
           k=kb;
           if ((*this)(s,i,j) == False)
           {
              while ( (Max == 1) && (k <= ke))
              {
                 l = lb;
                 while ( (Max == 1) && (l <= le))
                 {
                      if ((*this)(NextScale,k,l) == False) Max=0;
                      l++;
                 }
                 k++;
               }
               if (Max == 1) support(s,i,j) = VAL_SupDill;
            }
        }
    }
}

/****************************************************************************/

void MRNoiseModel::kill_event(int s, Ifloat &IEvent_Image, int FirstWin)
{
    int i,j,k,l,kb,ke,Nls,Ncs,sc=0;
    int Step=1;
    int Win = FirstWin/2;
    float Total;

    if (Set_Transform == TRANSF_PAVE)
    {
        Nls = TabNl(s);
        Ncs = TabNc(s);
        while (sc++ < s) Win*=2;

        for (i =0; i < Nl; i++)
        {
            Total = 0;
            kb = (i-Win < 0) ? 0: i-Win;
            ke = (i+Win >= Nls) ? Nls-1: i+Win;

            for (k =kb; k <= ke; k++)
            for (l =0; l <= Win; l++) Total += IEvent_Image(k, l);
            if ((support(s,i,0) == VAL_SupOK) && (Total < MinEventNumber))
                                                support(s,i,0) = VAL_SupMinEv;

            for (j =1; j < Nc; j+=Step)
            {

               for (k = -Win; k <= Win; k++)
               {
                   Total -= (int) IEvent_Image(kb,j-Win-1,I_ZERO);
                   Total += (int) IEvent_Image(kb,j+Win,I_ZERO);
               }
               if ((support(s,i,j) == VAL_SupOK) && (Total < MinEventNumber))
                                                  support(s,i,j) = VAL_SupMinEv;
           }
       }
    }
}

/****************************************************************************/

float  MRNoiseModel::sure_estimation(MultiResol &MR_Data)
{
   int s,i,j;
   int Ind = 1;
   unsigned long N = 0;
   for (s = 0; s < MR_Data.nbr_band()-1; s++)
        N +=  MR_Data.size_band_nl(s) * MR_Data.size_band_nc(s);
   float *Tab = new float [N+1];

   for (s = 0; s < MR_Data.nbr_band()-1; s++)
   for (i = 0; i < MR_Data.size_band_nl(s) ; i++)
   for (j = 0; j < MR_Data.size_band_nc(s) ; j++)
   {
       float Coef = MR_Data(s,i,j) / sigma(s,i,j);
       Tab[Ind++] = Coef*Coef;
   }
 
   sort(N, Tab);
   double MinRisk=0.;
   double CumVal = 0;
   int IndRisk=0;
   for (i = 1; i <= (int) N; i++) 
   {
      CumVal += Tab[i];
      int c = N-i;
      double s = CumVal + c * Tab[i];
      double r = ((double) N - (2.*i) + s) /  (double) N;
      if ((i == 1) || (r < MinRisk)) 
      {
         MinRisk = r;
	 IndRisk = i;
      }
   }
   double T = sqrt(Tab[IndRisk]);
   delete [] Tab;  
   return (float) T; 
}

/********************************************************************/ 

float MRNoiseModel::multi_sure_estimation(MultiResol & MR_Data, int b)
{
   int i,j;
   int Ind = 1;
   int Nl = MR_Data.size_band_nl(b);
   int Nc = MR_Data.size_band_nc(b); 
   unsigned long N = Nl*Nc;
   float *Tab = new float [N+1];
   for (i = 0; i < MR_Data.size_band_nl(b) ; i++)
   for (j = 0; j < MR_Data.size_band_nc(b) ; j++)
   {
       float Coef = MR_Data(b,i,j) / sigma(b,i,j);
       Tab[Ind++] =  Coef*Coef;
   }
   sort(N, Tab);
   double MinRisk=0.;
   double CumVal = 0;
   int IndRisk=0;
   for (i = 1; i <= (int) N; i++) 
   {
      CumVal += Tab[i];
      int c = N-i;
      double s = CumVal + c * Tab[i];
      double r = ((double) N - (2.*i) + s) /  (double) N;
      if ((i == 1) || (r < MinRisk)) 
      {
         MinRisk = r;
	 IndRisk = i;
      }
   }
   double T = sqrt(Tab[IndRisk]);
   delete [] Tab;  
   return (float) T; 
}

/********************************************************************/ 

void MRNoiseModel::set_support(MultiResol &MR_Data)
{
    int b,i,j;
    int Nls, Ncs;
    float SureLevel;
 
     if ((Transform == TO_DIADIC_MALLAT) && (GradientAnalysis == True))
     {
         for (b = 0; b < NbrBand-1; b+=2)
         {
            Nls = MR_Data.size_band_nl(b);
            Ncs = MR_Data.size_band_nc(b);
            for (i = 0; i < Nls;i++)
            for (j = 0; j < Ncs; j++)
            {
               float Coef = sqrt(POW(MR_Data(b,i,j),2.)+POW(MR_Data(b+1,i,j), 2.));
               if (signif(Coef, b,i,j) == True)
               {
                  support(b,i,j)=VAL_SupOK;
                  support(b+1,i,j)=VAL_SupOK;
               }
               else 
               {
                  support(b,i,j)=VAL_SupNull;
                  support(b+1,i,j)=VAL_SupNull;
               }
            }            
         }
     }
     else 
     {
        switch (TypeThreshold)
        {
	   case T_FDR:
            for (b = 0; b < NbrBand-1; b++)
            {
              Nls = MR_Data.size_band_nl(b);
              Ncs = MR_Data.size_band_nc(b);
	      dblarray TabPValue(Ncs,Nls);
	      for (i = 0; i < Nls;i++)
              for (j = 0; j < Ncs; j++) TabPValue(j,i) =  prob_noise(MR_Data(b,i,j), b, i, j);
	      float Alpha = (b==0) ? 1. - erf(NSigma[b] / sqrt(2.)): Alpha*2;
	      if (Alpha > 0.5) Alpha = 0.5;
	      
	      double PDet = fdr_pvalue(TabPValue.buffer(), Nls*Ncs, Alpha);
	      for (i = 0; i < Nls; i++)
              for (j = 0; j < Ncs; j++)
	      {
 		  if (TabPValue(j,i) < PDet)
		  {
		     support(b,i,j) = VAL_SupOK;
		     if ((OnlyPositivDetect == True) && (MR_Data(b,i,j) < 0))
	                                                      support(b,i,j) = VAL_SupNull;   
                     if (TabBandScale(b) < FirstDectectScale) support(b,i,j) = VAL_SupFirstScale;
		  }
		  else support(b,i,j) =  VAL_SupNull;
 	      }
	      TabEps[b] = PDet;
	      //if (PDet < DOUBLE_EPSILON) NSigma[b] = 5.;
	      //else 
	         // NSigma[b] = ABS(xerf(PDet/2.));
	        NSigma[b] = ABS(xerfc(0.5+(1-PDet)/2.));
		if((NSigma[b] < 5) && (NSigma[b] > 0))
			NSigma[b] = ABS(xerfc(0.5+(1-PDet)/2.));
		else
			NSigma[b] = 5;
	         // if (Verbose == True) 
	      printf("FDR: band %d ==> Alpha = %f, NSigma = %f\n", b+1, Alpha, NSigma[b]);
             }
	     break;
	   case T_SURE:
	       SureLevel = sure_estimation(MR_Data);
	       break;
           case T_KSIGMA: 
 	      break;
	   case T_UNIVERSAL: 
 	       for (b = 0; b < NbrBand-1; b++) 
	           NSigma[b] = sqrt(2.*log((float)(MR_Data.size_band_nl(b)*MR_Data.size_band_nc(b))));
	       break;
 	   case T_MRSURE: 
	       for (b = 0; b < NbrBand-1; b++)  
	                 NSigma[b] = multi_sure_estimation(MR_Data,b);	   
               break;
         }
         if (TypeThreshold != T_FDR) 
         {
            for (b = 0; b < NbrBand-1; b++)
            for (i = 0; i < MR_Data.size_band_nl(b); i++)
            for (j = 0; j < MR_Data.size_band_nc(b); j++)
                  support(b,i,j) = (signif(MR_Data(b,i,j), b,i,j) == True) ?  
	                                                    VAL_SupOK: VAL_SupNull;
        }
    }
    Bool SpatialAdaptiv = False;
    if ((Transform == TO_UNDECIMATED_NON_ORTHO) && (SpatialAdaptiv == True))
    {
        // Here we compare the SNR for each coeff in the three bands to the 
	// SNR of the sum of the coeff. If the SNR(SUM) > MAX(SNR_i) i=1,2,3 then we habe an isotropic feature
	// and we set the support to 1 in the three bands. 
        for (b = 0; b < NbrBand-1; b+=3)
        {
	   extern double TabNormPaveB3Spline[MAX_SCALE]; 
	   extern double TabNormUndecNonOrth_B3SPLINE_ISOTROP_2[MAX_SCALE]; 
	   float SumBand;
           Nls = MR_Data.size_band_nl(b);
           Ncs = MR_Data.size_band_nc(b);
 	   for (i = 0; i < Nls;i++)
           for (j = 0; j < Ncs;j++) 
	   {
 	      float Sig1 = sigma(b,i,j);  // std of the noise in horiz band
	      float Sig2 = sigma(b+1,i,j);// std of the noise in vertical band
	      float Sig3 = sigma(b+2,i,j); //std of the noise in diag band
	      float MaxSig = MAX(MR_Data(b,i,j)/Sig1, MR_Data(b+1,i,j)/Sig2);
	      MaxSig = MAX(MaxSig, MR_Data(b+3,i,j)/Sig3); // get the MAX SNR of the three bands
 	      float NSigMIN = MIN(NSigma[b], NSigma[b+1]); // 
	      NSigMIN = MIN(NSigMIN,NSigma[b+2]); // get the detection level parameter 
              float SigSum=0;
	      SumBand = (MR_Data(b,i,j) + MR_Data(b+1,i,j) + MR_Data(b+2,i,j)); // Sum of the three bands
 	      if (NewStatNoise == NOISE_GAUSSIAN)  
	      {
	         switch(U_Filter)
		 {
		    case  U_B3SPLINE:SigSum = TabNormPaveB3Spline[b/3]*SigmaNoise; break;
		    case  U_B3SPLINE_2:SigSum = TabNormUndecNonOrth_B3SPLINE_ISOTROP_2[b/3]*SigmaNoise; break;
		    default: SumBand = 0; SigSum = 1; break;
		 }
	      }
	      else SigSum = sqrt(Sig1*Sig1+Sig2*Sig2+Sig3*Sig3)*1.23; // Correction because of the correlation of the coef.
	      SumBand /= SigSum;
	      if ((SumBand > MaxSig) && (SumBand > NSigMIN))  
	      {
	         support(b,i,j) = support(b+1,i,j) = support(b+2,i,j) = VAL_SupOK;
	      }
	   }
        }
    }
   
/*
     if ((GetEgde == True) && (Transform == TO_PAVE_BSPLINE))
     {
        for (b = 0; b < NbrBand-1; b++)
	{
	   int b1,E_Nscale;
	   Nls = MR_Data.size_band_nl(b);
           Ncs = MR_Data.size_band_nc(b);
	   Ifloat Buff(Nls,Ncs,"Buff");	           
	   E_Nscale =  iround((float)log((float) (MIN(Nls,Ncs) / 4. * 3.) / log(2.)));
 	   if (one_level_per_pos_2d(TypeNoise) == True)
	   {
	       cout << "Error: this noise model cannot be used with this option ... " << endl;
	       exit(-1);
	   } 
           float Level = TabLevel[b] * NSigma[b];
           MultiResol MR_EDGE(Nls, Ncs, E_Nscale, TO_MALLAT, "MR EDGE");
	   MR_EDGE.TypeNorm = NORM_L2;
	   // MR_EDGE.SBFilter=F_HAAR;
	   MR_EDGE.SBFilter=F_MALLAT_7_9;
	   MR_EDGE.transform(MR_Data.band(b));
 	   for (b1= 0; b1 < MR_EDGE.nbr_band(); b1 ++)
	   for (i = 0; i < MR_EDGE.size_band_nl(b1); i++)
           for (j = 0; j < MR_EDGE.size_band_nc(b1); j++)
 	      if (ABS(MR_EDGE(b1,i,j)) < Level) MR_EDGE(b1,i,j) = 0.;
           MR_EDGE.recons(Buff);
	   
           for (i = 0; i < Nls;i++)
           for (j = 0; j < Ncs; j++)
               if ((support(b,i,j)==VAL_SupNull) &&
		   (ABS(Buff(i,j)) > FLOAT_EPSILON)) support(b,i,j)=VAL_SupEdge;	   
	 }    
      }
     
     if ((GetEgde == True)
          && ((Set_Transform == TRANSF_DIADIC_MALLAT) ||
	       (Set_Transform == TRANSF_MALLAT)))
     {
        for (b = 0; b < NbrBand-1; b++)
	{
           int s,b1,w,Np,LC_Nscale;
           details which_detail;
           MR_Data.band_to_scale(b, s, which_detail);
	
 	   switch (which_detail)
	   {
	    case D_HORIZONTAL: 
	           LC_Nscale =  iround((float)log((float) (Ncs / 4. * 3.) / log(2.))) - s ;
		   Np = Ncs;
		   break;
	    case  D_VERTICAL: 
	           LC_Nscale =  iround((float)log((float) (Nls / 4. * 3.) / log(2.))) - s ;
		   Np = Nls;
		   break;      
	    default: LC_Nscale = 0; break;
	   }
   
          if (LC_Nscale > 1)
          {	 
            MR_1D MR_LINE(Np, TO1_MALLAT, "line trans", LC_Nscale);
            MR_LINE.SB_Filter=F_HAAR;
            MR_LINE.Norm = NORM_L2; 
	    MR_LINE.Border = I_CONT;
            fltarray Col(Np);
            float Nsig = 2.*NSigma[b];
 	    if (which_detail == D_HORIZONTAL)
 	    {
             for (j = 0; j < Ncs; j++)
 	     {
	         for (i = 0; i < Np; i++) Col(i) =  MR_Data(b,i,j);
 	         MR_LINE.transform(Col);
	         // thresholding assuming Gaussian noise
		 for (b1= 0; b1 < MR_LINE.nbr_band(); b1 ++)
		 for (w = 0; w < MR_LINE.size_scale_np(b1); w++)
 	             if (ABS(MR_LINE(b1,w)) < Nsig*sigma(b,i,j)) MR_LINE(b1,w) = 0.;
  	         MR_LINE.recons(Col);
  	         for (i = 0; i < Np; i++)  
 		   if ((support(b,i,j)==VAL_SupNull) &&
		        (ABS(Col(i)) > FLOAT_EPSILON)) support(b,i,j)=VAL_SupEdge;
              }	     
	    }
	    else 
	    {  
 	      for (i = 0; i < Nls; i++)
 	      {
	         for (j = 0; j < Np; j++) Col(j) =  MR_Data(b,i,j);
 	         MR_LINE.transform(Col);
	         // thresholding assuming Gaussian noise
	         for (b1= 0; b1 < MR_LINE.nbr_band(); b1 ++)
		 for (w = 0; w < MR_LINE.size_scale_np(b1); w++)
 	             if (ABS(MR_LINE(b1,w)) < Nsig*sigma(b,i,j)) MR_LINE(b1,w) = 0.;
   	         MR_LINE.recons(Col);
 	         for (j = 0; j < Np; j++)
 		    if ((support(b,i,j)==VAL_SupNull) &&
		        (ABS(Col(j)) > FLOAT_EPSILON)) support(b,i,j)=VAL_SupEdge;
              }	   
	   }
        }
      }
    }
*/    
    if (SupIsol == True) for (b = 0; b < NbrBand-2; b++) kill_isol(b);
    if (DilateSupport == True) 
                       for (b = 0; b < NbrBand-1; b++) dilate_support(b);
}

/****************************************************************************/

void MRNoiseModel::mod_support(MultiResol &MR_Data, fltarray Nsigma)
{
    int i,j,k,s;
    int Nls, Ncs;

     for (s = 0; s < NbrBand-1; s++) {
     
        Nls = MR_Data.size_band_nl(s);
        Ncs = MR_Data.size_band_nc(s);
        for (i = 0; i < Nls;i++)
        for (j = 0; j < Ncs; j++) 
	{
 	    if (support(s,i,j) == VAL_SupOK) 
	    {
 	     	for (k=0; k<NbrBand; k++)
            		if (signif(MR_Data(k,i,j), k,i,j, Nsigma) == True)
                             support(k,i,j)=VAL_SupOK;
            		else support(k,i,j)=VAL_SupNull;
	    }
        }
    }
    //if (SupIsol == True) for (s = 0; s < NbrBand-2; s++) kill_isol(s);
    //if (DilateSupport == True) 
    //                   for (s = 0; s < NbrBand-1; s++) dilate_support(s);

}
/****************************************************************************/

Bool MRNoiseModel::kill_coef(int b,int i,int j, float Val, Bool SetSupport)
{
    Bool ValRet=False;

    if ((*this)(b,i,j) == False)
    {
        if (SetSupport == False) ValRet = True;
        else
        {
            if ((support(b,i,j) == VAL_SupNull) &&
                (signif(Val, b,i,j) == True))
            {
               support(b,i,j)=VAL_SupOK;
            }
            else ValRet = True;
        }
    }
    return ValRet;
}

/****************************************************************************/

void MRNoiseModel::threshold(MultiResol &MR_Data, Bool SetSupport)
{
    int i,j,b;
    int Nls, Ncs;

    for (b = 0; b < NbrBand-1; b++)
    {
       // int w, LC_Nscale, Np;
       int b1,s;
       details which_detail;
       MR_Data.band_to_scale(b, s, which_detail);
       Nls = TabNl(b);
       Ncs = TabNc(b);
       
       if ((GetEgde == True) && (Transform == TO_PAVE_BSPLINE))
       {
 	   int E_Nscale;
 	   Ifloat Buff(Nls,Ncs,"Buff");	           
	   E_Nscale =  iround((float)log((float) (MIN(Nls,Ncs) / 4. * 3.) / log(2.)));
  	   if (one_level_per_pos_2d(TypeNoise) == True)
	   {
	       cout << "Error: this noise model cannot be used with this option ... " << endl;
	       exit(-1);
	   } 
           float Level = TabLevel[b] * NSigma[b];
           MultiResol MR_EDGE(Nls, Ncs, E_Nscale, TO_MALLAT, "MR EDGE");
	   MR_EDGE.TypeNorm = NORM_L2;
	   // MR_EDGE.SBFilter=F_HAAR;
	   MR_EDGE.SBFilter=F_MALLAT_7_9;
	   MR_EDGE.transform(MR_Data.band(b));
 	   for (b1= 0; b1 < MR_EDGE.nbr_band(); b1 ++)
	   for (i = 0; i < MR_EDGE.size_band_nl(b1); i++)
           for (j = 0; j < MR_EDGE.size_band_nc(b1); j++)
 	      if (ABS(MR_EDGE(b1,i,j)) < Level) MR_EDGE(b1,i,j) = 0.;
           MR_EDGE.recons(Buff);
	   
           for (i = 0; i < Nls;i++)
           for (j = 0; j < Ncs; j++)  
              if (kill_coef(b,i,j,MR_Data(b,i,j),SetSupport) == True) 
                                              MR_Data(b,i,j) = Buff(i,j);
       } 
         
       
/*       if ((GetEgde == True)
          && ((Set_Transform == TRANSF_DIADIC_MALLAT) ||
	       (Set_Transform == TRANSF_MALLAT)))
       {
 	  switch (which_detail)
	  {
	    case D_HORIZONTAL: 
	           LC_Nscale =  iround((float)log((float) (Ncs / 4. * 3.) / log(2.))) - s ;
		   Np = Ncs;
		   break;
	    case  D_VERTICAL: 
	           LC_Nscale =  iround((float)log((float) (Nls / 4. * 3.) / log(2.))) - s ;
		   Np = Nls;
		   break;      
	    default: LC_Nscale = 0; break;
	  }
      }
      else LC_Nscale = 0;
 LC_Nscale = 0;
      if (LC_Nscale > 1)
      {	 
         MR_1D MR_LINE(Np, TO1_MALLAT, "line trans", LC_Nscale);
         MR_LINE.SB_Filter=F_HAAR;
         MR_LINE.Norm = NORM_L2; 
	 MR_LINE.Border = I_CONT;
         fltarray Col(Np);
         float Nsig = 2.*NSigma[b];
 	 if (which_detail == D_HORIZONTAL)
 	 {
             for (j = 0; j < Ncs; j++)
 	     {
	         for (i = 0; i < Np; i++) Col(i) =  MR_Data(b,i,j);
 	         MR_LINE.transform(Col);
	         // thresholding assuming Gaussian noise
		 for (b1= 0; b1 < MR_LINE.nbr_band(); b1 ++)
		 for (w = 0; w < MR_LINE.size_scale_np(b1); w++)
 	             if (ABS(MR_LINE(b1,w)) < Nsig*sigma(b,i,j)) MR_LINE(b1,w) = 0.;
  	         MR_LINE.recons(Col);
  	         for (i = 0; i < Np; i++)  
 		   if (kill_coef(b,i,j,MR_Data(b,i,j),SetSupport) == True)
                                                         MR_Data(b,i,j)=Col(i);
              }	     
	  }
	  else 
	  {  
 	      for (i = 0; i < Nls; i++)
 	      {
	         for (j = 0; j < Np; j++) Col(j) =  MR_Data(b,i,j);
 	         MR_LINE.transform(Col);
	         // thresholding assuming Gaussian noise
	         for (b1= 0; b1 < MR_LINE.nbr_band(); b1 ++)
		 for (w = 0; w < MR_LINE.size_scale_np(b1); w++)
 	             if (ABS(MR_LINE(b1,w)) < Nsig*sigma(b,i,j)) MR_LINE(b1,w) = 0.;
		     
  	         MR_LINE.recons(Col);
 	         for (j = 0; j < Np; j++)
 		    if (kill_coef(b,i,j,MR_Data(b,i,j),SetSupport) == True) 
                                                        MR_Data(b,i,j)=Col(j);
              }	   
	  }
       } */
       else
       {
          for (i = 0; i < Nls; i++)
          for (j = 0; j < Ncs; j++)
           if (kill_coef(b,i,j,MR_Data(b,i,j),SetSupport) == True) 
                                                     MR_Data(b,i,j)=0.;
       }
   }
}

/****************************************************************************/

void MRNoiseModel::refresh_support(MultiResol &MR_Data)
{
    int i,j,s;
    int ValSup;
    int Nls, Ncs;
    float CoefMr;
    
    for (s = 0; s < NbrBand-1; s++)
    {
       Nls = TabNl(s);
       Ncs = TabNc(s);
       for (i = 0; i < Nls; i++)
       for (j = 0; j < Ncs; j++)
       {
           CoefMr = MR_Data(s,i,j);
	   ValSup = support(s,i,j);
           if ( (ValSup == VAL_SupNull)  && (signif(CoefMr,s,i,j) == True))
                                   support(s,i,j)=VAL_SupOK;
       }
    }
}

/****************************************************************************/

void MRNoiseModel::weight_snr(MultiResol &MR_Data, Bool SetSupport)
{
    int i,j,s;
    int ValSup;
    int Nls, Ncs;
    float CoefMr,Weight;
    
    for (s = 0; s < NbrBand-1; s++)
    {
       Nls = TabNl(s);
       Ncs = TabNc(s);
       for (i = 0; i < Nls; i++)
       for (j = 0; j < Ncs; j++)
       {
           CoefMr = MR_Data(s, i, j);
           ValSup = support(s,i,j);
           if ( (SetSupport == True) && (ValSup == VAL_SupNull) 
                && (signif(CoefMr,s,i,j) == True))
                                   support(s,i,j)=VAL_SupOK;

          Weight = ABS(CoefMr) / (NSigma[s]*sigma(s,i,j));
          if (Weight > 1.) Weight = 1.;
          MR_Data(s, i, j) *= Weight;
       }
    }
}
 
/****************************************************************************/

void MRNoiseModel::weight_invsnr(MultiResol &MR_Data, Bool SetSupport)
{
    int i,j,s;
    float CoefMr,Weight;
    int ValSup;
    int Nls, Ncs;

    for (s = 0; s < NbrBand-1; s++)
    {
       Nls = TabNl(s);
       Ncs = TabNc(s);
       for (i = 0; i < Nls; i++)
       for (j = 0; j < Ncs; j++)
       {
           CoefMr = MR_Data(s, i, j);
           ValSup = support(s,i,j);
	   if (   (SetSupport == True) && (ValSup == VAL_SupNull) 
            && (signif(CoefMr,s,i,j) == True))
                                   support(s,i,j)=VAL_SupOK;

           Weight = ABS(CoefMr) / (NSigma[s]*sigma(s,i,j));
           if (Weight > 1.) Weight = 1.;
           MR_Data(s, i, j) *= (1. - Weight);
       }
    }
}

/****************************************************************************/

void MRNoiseModel::set_sigma(Ifloat & Imag, MultiResol &MR_Data)
{
    extern float mr_tab_noise(int s);
    int b, s,i,j,Ind;

    noise_compute (MR_Data);

//cout << "noise_compute: set_sigma : Noise = " << SigmaNoise << endl;
//Imag.info("MM");
    if ((SigmaNoise < FLOAT_EPSILON) && (NewStatNoise == NOISE_GAUSSIAN))
    {
        switch (SigmaDetectionMethod)
        {
            case SIGMA_MEDIAN: SigmaNoise = detect_noise_from_med (Imag); break;
            case SIGMA_BSPLINE:
                       SigmaNoise = detect_noise_from_bspline (Imag); break;
            case SIGMA_CLIPIMA:SigmaNoise = detect_noise_sigma (Imag); break;
            case SIGMA_SUPPORT: 
                    SigmaNoise=detect_noise_from_support(Imag,MR_Data,(*this));
                    break;
            default: cerr<< "Error: MRNoiseModel:set_sigma, unknown method" << endl;
                     exit(0);break;
        }
    }
 // cout << "After set_sigma : Noise = " << SigmaNoise << endl;

   switch (NewStatNoise)
   {
      case NOISE_SPECKLE:
       {
        // we assume that the level at 3.719sigma (eps = 1e-4) gives
	// a good fit of the gaussian law to the rayley noise distribution
        fltarray R_Eps(NbrBand-1);
	fltarray Tmin(NbrBand-1);
	fltarray Tmax(NbrBand-1);
        for (s = 0; s < NbrBand-1; s++)  R_Eps(NbrBand-1) = 1.e-04;
	CSpeckle->find_threshold(NbrBand,R_Eps, Tmin, Tmax);
        for (s = 0; s < NbrBand-1; s++)
	{
	    float Max = ABS(Tmin(s));
	    if (ABS(Tmin(s)) < ABS(Tmax(s))) Max = Tmax(s);
            TabLevel[s] =  Max / 3.719;
 	 }
        }
        break;
      case NOISE_CORREL:
          for (s = 0; s < NbrBand-1; s++)  
	        TabLevel[s] = CorrelNoiseMap->sigma_band(s);
         break;
      case NOISE_GAUSSIAN:
        if ((Transform == TO_DIADIC_MALLAT) && (GradientAnalysis == True))
        {
           extern double TabNormGradDiadicMallat[MAX_SCALE];
           for (s = 0; s < NbrBand-1; s++)
                          TabLevel[s] = SigmaNoise*TabNormGradDiadicMallat[s/2];
        }
        else  for (s = 0; s < NbrBand-1; s++)
                          TabLevel[s] = SigmaNoise*MR_Data.band_norm(s);  
        break;
      case  NOISE_NON_UNI_ADD:
        if ((Set_Transform == TRANSF_MALLAT) || 
	    (Set_Transform == TRANSF_FEAUVEAU))  
        {
           cerr << "Error: this kind of transform cannot be used " << endl;
           cerr << "       with this noise model " << endl;
           exit(-1);
        }
        if (UseRmsMap == False)
        {
           // Ifloat ImaSigmaNoise(Nl, Nc, "ImaSigmaNoise");
           Ifloat Buff(Nl, Nc, "Buff");
           RmsMap.alloc(Nl, Nc, "ImaSigmaNoise");
           smooth_mediane (Imag, Buff, I_MIRROR, 0, 3);
           Buff = Imag - Buff;
           im_sigma(Buff, RmsMap, SizeBlockSigmaNoise, NiterSigmaClip);
	   // io_write_ima_float("xx_rms.fits", RmsMap);
        } 
        // calculate the true sigma at each scale :
	if ((RmsMap.nl() != Nl) || (RmsMap.nc() != Nc))
        {
           cerr << "Error: RMS map is not correctly initialized ... " << endl;
           exit(0);
        }
	if (((Transform == TO_PAVE_BSPLINE) && (GetRmsCoeffGauss == True))
	          || (Transform == TO_UNDECIMATED_MALLAT)
		  || (Transform == TO_UNDECIMATED_NON_ORTHO)
		  || (Transform == TO_PAVE_FEAUVEAU)
		  || (Transform == TO_DIADIC_MALLAT))
 	        get_sigma_from_rms_map(RmsMap);
	else max_sigma_in_scale(RmsMap);
 	UseRmsMap=True;
       break;
    case NOISE_UNI_UNDEFINED:
        for (b = 0; b < NbrBand-1; b ++)
        {
           //float SigmaScale, MeanScale;
           //sigma_clip(MR_Data.band(b), MeanScale, SigmaScale, NiterSigmaClip);
           //TabLevel[b] = SigmaScale;          
	   TabLevel[b] = detect_noise_from_mad(MR_Data.band(b), MadEstimFromCenter);
        }
       break;
    case NOISE_UNDEFINED:
      {
       Ifloat ImaSigmaNoise(Nl, Nc, "ImaSigmaNoise");
       int BlockSize = SizeBlockSigmaNoise;
       for (b = 0; b < NbrBand-1; b++)
       {
           int Nlb = MR_Data.size_band_nl(b);
           int Ncb = MR_Data.size_band_nc(b);
           float SigmaScale, MeanScale;
           //details which_detail;
           //MR_Data.band_to_scale(b, s, which_detail);
           sigma_clip(MR_Data.band(b), MeanScale, SigmaScale, NiterSigmaClip);
           ImaSigmaNoise.resize(Nlb, Ncb);
	   // If we want the full resolution, we need to uncomment the next line.
	   // Otherwise the noise standard deviation is calculated on for position
	   // at a distance of BlockSize/2
           // im_sigma(MR_Data.band(b), ImaSigmaNoise, BlockSize, NiterSigmaClip);
	   // use the MAD instead of sigma_clipping: im_sigma_block_mad(MR_Data.band(b), ImaSigmaNoise, BlockSize, 32);
           im_sigma_block(MR_Data.band(b), ImaSigmaNoise, BlockSize, NiterSigmaClip, BlockSize/2);
           for (i = 0; i < Nlb;i++)
           for (j = 0; j < Ncb; j++)
           {
               Ind = index(b,i,j);
               TabLevel[Ind] = ImaSigmaNoise(i,j); // MAX(SigmaScale, ImaSigmaNoise(i,j));
           }
           if ((Set_Transform == TRANSF_PAVE) 
	        || ((Set_Transform == TRANSF_UNDECIMATED_MALLAT) && (b%3 == 2))
		|| ((Set_Transform == TRANSF_DIADIC_MALLAT)  && (b%2 == 1)))
		BlockSize *=2;
        }
     }
     break;
   case NOISE_EVENT_POISSON:
     {
        Iint EventCount(Nl,Nc,"ImagCount");
        float alpha;
        for (s = 0; s < MR_Data.nbr_band()-1; s++)
        {
	    alpha=1.;
            for (int sc=0; sc < s; sc++) alpha *= 4.;
            event_one_scale(Event_Image, s, EventCount, MR_Data.Border);
            for (i=0;i<Nl;i++)
            for (j=0;j<Nc;j++)
            {
	     Ind = index(s,i,j);
	     TabLevel[Ind] = 
	        SIGMA_BSPLINE_WAVELET * sqrt((float) EventCount(i,j)) / alpha;
	    }	    
	}
      }	 
     break;
   default: break;
  }
}

/****************************************************************************/

void MRNoiseModel::speckle_support(MultiResol &MR_Data)
{
   int Nlb, Ncb, i,j,b;
   float Max;
   fltarray R_Eps(NbrBand-1);
   fltarray Tmin(NbrBand-1);
   fltarray Tmax(NbrBand-1);
         
   for (b = 0; b < NbrBand-1; b++) R_Eps(b) = TabEps[b];
   CSpeckle->find_threshold(MR_Data.nbr_band(), R_Eps, Tmin, Tmax);

   for (b = 0; b < MR_Data.nbr_band()-1; b++)
   {
      Max = ABS(Tmin(b));
      if (Max < ABS(Tmax(b))) Max = Tmax(b);
      TabLevel[b] =  ABS(Max) / NSigma[b];
	    
      Nlb = MR_Data.size_band_nl(b);
      Ncb = MR_Data.size_band_nc(b);
      for (i = 0; i < Nlb;i++)
      for (j = 0; j < Ncb; j++)
      {
         if (signif(MR_Data(b,i,j), b,i,j,Tmin(b),Tmax(b)) == True)
                                            support(b,i,j)=VAL_SupOK;
         else support(b,i,j)=VAL_SupNull;
      }   
   }
   
    if (SupIsol == True) for (b = 0; b < NbrBand-2; b++) kill_isol(b);
    if (DilateSupport == True) 
                       for (b = 0; b < NbrBand-1; b++) dilate_support(b);
}

/****************************************************************************/

void MRNoiseModel::correl_support(MultiResol &MR_Data)
{
   int Nlb, Ncb, i,j,b;
   fltarray R_Eps(NbrBand-1);
   fltarray Tmin(NbrBand-1);
   fltarray Tmax(NbrBand-1);
   
   for (b = 0; b < NbrBand-1; b++) R_Eps(b) = TabEps[b];
   CorrelNoiseMap->find_threshold(MR_Data.nbr_band(), R_Eps, Tmin, Tmax);

   for (b = 0; b < MR_Data.nbr_band()-1; b++)
   {
      TabLevel[b] = CorrelNoiseMap->sigma_band(b);
      Nlb = MR_Data.size_band_nl(b);
      Ncb = MR_Data.size_band_nc(b);
      for (i = 0; i < Nlb;i++)
      for (j = 0; j < Ncb; j++)
      {
         if (signif(MR_Data(b,i,j), b,i,j,Tmin(b),Tmax(b)) == True)
                                            support(b,i,j)=VAL_SupOK;
         else support(b,i,j)=VAL_SupNull;
      }   
   }
   
    if (SupIsol == True) for (b = 0; b < NbrBand-2; b++) kill_isol(b);
    if (DilateSupport == True) 
                       for (b = 0; b < NbrBand-1; b++) dilate_support(b);
}

/****************************************************************************/

void MRNoiseModel::model(Ifloat & Imag, MultiResol &MR_Data)
{
    NewStatNoise = TypeNoise;
     //int s;

    // if TypeNoise != NOISE_EVENT_POISSON and MinEvent == True
    // the parameter Event_Image must be intialized before the 
    // call to model
    if ((MinEvent == True) &&  (TypeNoise != NOISE_EVENT_POISSON))
    {
       if (  (TypeNoise == NOISE_EVENT_POISSON) 
             && (Imag.nl() != Event_Image.nl())
             && (Imag.nc() != Event_Image.nc()))
       {
        cerr << "Error: MRNoiseModel::model : Event_Image is not correctly Initialized ..." << endl;
        exit(-1);
       }
    }

    // if TypeNoise == NOISE_EVENT_POISSON and TransImag == False
    // the parameter Event_Image must be intialized before the 
    // call to model
    if (  (TypeNoise == NOISE_EVENT_POISSON) && (TransImag == False) &&
          ((Imag.nl() != Event_Image.nl()) || (Imag.nc() != Event_Image.nc())))
    {
        cerr << "Error: MRNoiseModel::model : Event_Image is not correctly Initialized ..." << endl;
        exit(-1);
    }


    switch (TypeNoise)
    {
        case NOISE_POISSON:
        case NOISE_GAUSS_POISSON:
	       SigmaNoise = 1.;
	case NOISE_MULTI:
               if (TransImag == True) im_transform(Imag);
               NewStatNoise = NOISE_GAUSSIAN;
               break;
        case NOISE_GAUSSIAN:   
        case NOISE_NON_UNI_ADD:
        case NOISE_UNI_UNDEFINED:
        case NOISE_UNDEFINED:
	case NOISE_CORREL:
                break;
        case NOISE_NON_UNI_MULT:
               if (TransImag == True) im_transform(Imag);
               NewStatNoise = NOISE_NON_UNI_ADD;
               break;
        case NOISE_EVENT_POISSON:
	        if (TransImag == True) im_transform(Imag);
		// in case of Poisson noise, the MIRROR border must be used
		MR_Data.Border = I_MIRROR;
		TransImag = False;
		SigmaNoise = 1.;
		break;
	case NOISE_SPECKLE:
               if (TransImag == True) im_transform(Imag);
               break;
     }
     
    // transform the input data
    Border = MR_Data.Border;
    U_Filter = MR_Data.U_Filter;
    TypeNorm = MR_Data.TypeNorm;
    FilterBank = MR_Data.filter_bank();
    NbrUndecimatedScale = MR_Data.nbr_undec_scale();
    MR_Data.SigmaNoise = SigmaNoise;
    MR_Data.transform(Imag);
    
    switch (TypeNoise)
    {
        case NOISE_POISSON:
        case NOISE_GAUSS_POISSON:
	case NOISE_MULTI:
	case NOISE_NON_UNI_MULT:
 	       set_sigma(Imag, MR_Data);
               set_support(MR_Data);
               if (TransImag == True) im_invtransform(Imag);
               break;
        case NOISE_GAUSSIAN: 
        case NOISE_NON_UNI_ADD:
        case NOISE_UNDEFINED:
        case NOISE_UNI_UNDEFINED:
	      set_sigma(Imag, MR_Data);
 	      set_support(MR_Data);
	      break;
	case NOISE_CORREL:
	       if (CorrelNoiseMap != NULL) delete CorrelNoiseMap;
	       CorrelNoiseMap = new StatNoiseMap(RmsMap, NbrScale, Transform, 
	                          FilterBank, TypeNorm, NbrUndecimatedScale, U_Filter);
	       if (SigmaApprox == False) correl_support(MR_Data);
	       else
	       {
	          set_sigma(Imag, MR_Data);
                  set_support(MR_Data);
	       }
               break;
	case NOISE_EVENT_POISSON:
	       if (SigmaApprox == False)  // Default
	       {
	          Bool WriteAllInfo = False;
	          mr2d_psupport( MR_Data, (*this), MR_Data.Border, WriteAllInfo );
                  //mr_psupport(MR_Data, (*this), MR_Data.Border);
                  if (SupIsol == True) for (int b = 0; b < NbrBand-2; b++) kill_isol(b);
                }
		else
		{
		  set_sigma(Imag, MR_Data);
                  set_support(MR_Data);
		}	     
		break;
	case NOISE_SPECKLE:
	       if (SigmaApprox == False) speckle_support(MR_Data);
	       else
 	       {
		  set_sigma(Imag, MR_Data);
                  set_support(MR_Data);
	       }
               if (TransImag == True) im_invtransform(Imag);
               break;
     }
 
    // ALREADY done in   MR_Psup.cc 
    // verify if minimum number of count OK
//     if (MinEvent == True) 
//     {
//       //  if (TypeNoise == NOISE_EVENT_POISSON) 
//       //           for (int s=0; s < NbrBand-1; s++) kill_event(s, Imag);
//     }
}

/****************************************************************************/

void MRNoiseModel::model(Ifloat & Imag)
{
   // MultiResol MR_Data(Nl, Nc, NbrScale, Transform, "MRNoiseModel");
   MultiResol MR_Data;
   MR_Data.alloc (Nl, Nc, NbrScale, Transform, FilterBank, TypeNorm, NbrUndecimatedScale, U_Filter);
   model(Imag, MR_Data);
}

/****************************************************************************/

void MRNoiseModel::free()
{
   if (Size > 0)
   {
      if (one_level_per_pos_2d(TypeNoise) == True) free_buffer((char *) TabLevel);
      else delete [] TabLevel;
      free_buffer((char *) TabSupport);
      NbrScale = 0;
      Nl = 0;
      Nc = 0;
      Size = 0;
      //if (CEventPois != NULL) delete CEventPois;
      if (CFewEventPoisson2d != NULL) delete CFewEventPoisson2d;
      if (CFewEvent2d != NULL) delete CFewEvent2d;
      if (CSpeckle != NULL) delete CSpeckle;
      if (CorrelNoiseMap  != NULL) delete CorrelNoiseMap; 
   }
}
/****************************************************************************/

MRNoiseModel::~MRNoiseModel()
{
   (*this).free();
}

/****************************************************************************/

void MRNoiseModel::mr_obj(MultiResol &MR_Data)
{
   int i,j,s;
   
   for (s = 0; s < NbrBand-1; s++)
   for (i = 0; i < TabNl(s);i++)
   for (j = 0; j < TabNc(s); j++)
   {
      if ((*this)(s,i,j)  == True)  MR_Data(s,i,j) = 1.;
      else MR_Data(s,i,j) = 0.;
   }
}

/****************************************************************************/

void MRNoiseModel::write_support_mr(char *FileName)
{
   MultiResol MR_Data;
   MR_Data.alloc (Nl, Nc, NbrScale, Transform, FilterBank, 
                  TypeNorm, NbrUndecimatedScale);
    (*this).mr_obj(MR_Data);
    MR_Data.write(FileName);
}

/****************************************************************************/

void MRNoiseModel::write_support_ima(char *FileName)
{
   MultiResol MR_Data;
   MR_Data.alloc (Nl, Nc, NbrScale, Transform, FilterBank, 
                  TypeNorm, NbrUndecimatedScale);
   (*this).mr_obj(MR_Data);
   float Coef;
   Ifloat Ima(Nl, Nc, "MRNoiseModel write");
   int s,i,j;

   Ima.init();
    switch (Set_Transform)
    {
        case TRANSF_PAVE:
        case TRANSF_PYR:
	case TRANSF_DIADIC_MALLAT:
	case TRANSF_UNDECIMATED_MALLAT:
	case TRANSF_SEMIPYR:
          {
             Ifloat Dat(Nl, Nc, "MRNoiseModel write");
             for (s = NbrBand-2; s >= 0; s--)
             {
                Coef = (float) ( 1 << (NbrBand-1-s));
                im_block_extend(MR_Data.band(s), Dat);
                // cout << "Coef = " << Coef << endl;
                // cout << MR_Data.band(s).nl() << " to " << Dat.nl() << endl;
                // cout << MR_Data.band(s).nc() << " to " << Dat.nc() << endl;
                for (i = 0; i < Dat.nl();i++)
                for (j = 0; j < Dat.nc(); j++)
                   if (Dat(i,j) > FLOAT_EPSILON) Ima(i,j) = Coef;
             }
          }
          break;
        case TRANSF_MALLAT:
        case TRANSF_FEAUVEAU:
           ortho_trans_to_ima(MR_Data, Ima);
           break;
        default:
          cerr << "Error: Unknown set transform ... " << endl;
          exit(0);
          break;

    }

    io_write_ima_float(FileName, Ima);
}

void MRNoiseModel::set_old_poisson( Bool Flag ) {
   mOldPoisson = Flag;
   if( mOldPoisson ) 
      std::cout << "!!! Odl poisson few event class is used !!!"
                << std::endl;
}

Bool MRNoiseModel::old_poisson() {
   return mOldPoisson;
}

void MRNoiseModel::write_threshold( Bool Flag ) {
   mWriteThreshold = Flag;
}


void MRNoiseModel::write_in_few_event( Bool Write ) {
   if( CFewEvent2d != (FewEvent*)NULL )
      CFewEvent2d->set_write_all( Write ) ;
}



/****************************************************************************/

void MRNoiseModel::trace() {

   cout << "class MRNoiseModel" << endl;
   cout << "------------------" << endl;
   
   cout << "NbrScale = " << NbrScale << endl;
   cout << "NbrBand = " << NbrBand  << endl;
   cout << "Nl = " << Nl << endl;
   cout << "Nc = " << Nc << endl;
   cout << "Size = " << Size << endl;
   cout << "Transform = " << StringTransform(Transform) << endl;
   cout << "Set_Transform = " << StringSetTransform(Set_Transform) << endl;
   cout << "TypeNoise = " << StringNoise(TypeNoise) << endl;
   cout << "Border = " << (int)Border << endl;
   cout << "NSigma : " << endl;
   for (int i=0;i<NbrBand;i++)
      cout << " scale[i] : " << NSigma[i] << endl;
   if (SigmaApprox) cout << "SigmaApprox = true" << endl;
   cout << "TypeNorm = " << (int)DEF_SB_NORM << endl;
   cout << "NbrUndecimatedScale = " << (int)NbrUndecimatedScale << endl;
   if (GradientAnalysis) cout << "GradientAnalysis = true" << endl;
   if (NoCompDistrib) cout << "NoCompDistrib = true" << endl;
   cout << "Size  = " << (int)Size << endl;
   cout << "BadPixelVal = " << BadPixelVal << endl;
   if (BadPixel) cout << "BadPixel = true" << endl;
   if (MinEvent) cout << "MinEvent = true" << endl;
   if (SupIsol) cout << "SupIsol = true" << endl; 
   cout << "FirstDectectScale = " << FirstDectectScale << endl;
   if (OnlyPositivDetect) cout << "OnlyPositivDetect = true" << endl; 
   if (DilateSupport) cout << "DilateSupport = true" << endl; 
   if (GetEgde) cout << "GetEgde = true" << endl; 
   cout << "MinEventNumber = " << MinEventNumber << endl;
   cout << "SigmaDetectionMethod = " << (int)SigmaDetectionMethod << endl;
   cout << "NiterSigmaClip = " << (int)NiterSigmaClip << endl;
   cout << "SizeBlockSigmaNoise = " << SizeBlockSigmaNoise << endl;
   if (TransImag) cout << "TransImag = true" << endl; 
   if (UseRmsMap) cout << "UseRmsMap = true" << endl; 
   if (GetRmsCoeffGauss) cout << "GetRmsCoeffGauss = true" << endl; 
   if (SigmaApprox) cout << "SigmaApprox = true" << endl; 

   
   if (TabLevel == NULL)  cout << "TabLevel = NULL" << endl;
   if (TabSupport == NULL) cout << "TabSupport = NULL" << endl;
   //if (CEventPois == NULL) cout << "CEventPois = NULL" << endl;
   if (CFewEventPoisson2d == NULL) cout << "CFewEventPoisson2d = NULL" << endl;
   if (CFewEvent2d == NULL) cout << "CFewEventPoisson2d = NULL" << endl;
   if (CSpeckle == NULL)   cout << "CSpeckle = NULL" << endl; 
   if (CorrelNoiseMap == NULL)   cout << "CorrelNoiseMap = NULL" << endl; 
   if (FilterBank == NULL)   cout << "FilterBank = NULL" << endl; 
 
   cout << "end class MRNoiseModel" << endl;  
 
}
