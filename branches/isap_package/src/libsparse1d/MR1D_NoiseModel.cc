 /******************************************************************************
**                   Copyright (C) 1996 by CEA
*******************************************************************************
**
**    UNIT
**
**    Version: 0.1
**
**    Author: Rene Gastaud & Jean-Luc Starck
**
**    Date:  22-OCT-1996
**    
**    File:  MR1D_NoiseModel.cc
**
*******************************************************************************
**
**    DESCRIPTION  noise modelling for a 1 dimensional signal
**    -----------  
**                 
*******************************************************************************
**
** MR1DNoiseModel::MR1DNoiseModel(type_noise TNoise, int Nelem, 
**             int ScaleNumber, type_trans_1d Trans)
**
** Class builder
** type_noise TNoise is the type of noise in the data
** Nelem size of the input signal
** ScaleNumber: number of scales of the transform
** type_trans_1d Trans: type of multiresoltuion transform
**
******************************************************************************/
 
// static char sccsid[] = "@(#)MR1D_NoiseModel.cc 0.1 96/10/23 CEA 1996 @(#)";

#include "GlobalInc.h" // for Bool
#include "IM_Obj.h"  //rg type_border
#include "IM_Noise.h"  // for noise_poisson_transform
#include "MR1D_NoiseModel.h"  // which include MR1D_Obj.h and MR1D_filter.h
#include "FFTN_1D.h"  // for DEFAULT_EPSILON
#include "NR.h"

#define DEBUG_MODEL 0 
Bool TraceFewEvent=False;
Bool OutForTest=False;

extern void building_signal_signal (fltarray &Signal, intarray &Event_Image);

/************************************************************************/

static Bool one_level_per_pos(type_noise TNoise)
// return True if we need only one level per position
// in the multiresolution space
{
    Bool ValRet= False;

    if ((TNoise == NOISE_EVENT_POISSON) || (TNoise == NOISE_NON_UNI_ADD)
         || (TNoise == NOISE_NON_UNI_MULT)  || (TNoise == NOISE_UNDEFINED))  ValRet = True;

    return ValRet;
}

/************************************************************************/

void MR1DNoiseModel::init_param()
{
   NbrScale = 0;
   Nelem = 0;
   for (int i=0; i<MAX_SCALE_1D; i++) {
      TabN[i] = 0;
      TabPos[i] = 0;
      TabEps[i] = 0;
      NSigma[i] = 0;
   }
   Size = 0;
   TabSupport = NULL;
   TabLevel = NULL;
   TypeNoise = (type_noise)0;
   NewStatNoise = (type_noise)0;
   Type_Transform = (type_trans_1d)0;
   TransImag = False;
   NiterSigmaClip = 0;
   SizeBlockSigmaNoise = 0;
   SupIsol = False;
   MinEvent = False;
   OnlyPositivDetect = False;
   OnlyNegativDetect = False;
   MinEventNumber = DEF_MINEVENT_NUMBER;
   BadPixel = False;
   DilateSupport = False;
   BadPixelVal = 0.;
   FirstDectectScale = 0;
   SigmaDetectionMethod = (type_sigma_method_1d)0;
   CCD_Gain = CCD_ReadOutSigma = CCD_ReadOutMean = SigmaNoise = 0.0;
   FilterBank = NULL;
   TypeNorm = DEF_SB_NORM;
   TypeThreshold = DEF_THRESHOLD;
   UseRmsMap = False;
}

/************************************************************************/

MR1DNoiseModel::MR1DNoiseModel () 
{
  init_param();
}   
   
/************************************************************************/

MR1DNoiseModel::MR1DNoiseModel(type_noise TNoise, int N_signal,
                      int ScaleNumber, type_trans_1d Trans) 
{
   init_param();	      
   alloc (TNoise, N_signal, ScaleNumber, Trans);		      
}

/************************************************************************/

MR1DNoiseModel::MR1DNoiseModel(type_noise TNoise, int N_signal,
                      int ScaleNumber, type_trans_1d Trans,
                      FilterAnaSynt *FAS, sb_type_norm Norm) 
{
   init_param();	      
   alloc (TNoise, N_signal, ScaleNumber, Trans, FAS, Norm);		      
}
 
/************************************************************************/

void MR1DNoiseModel::alloc (type_noise TNoise, int N_signal, 
	                    int ScaleNumber, type_trans_1d Trans,
                            FilterAnaSynt *FAS, sb_type_norm Norm)
{
    int s;
    int index ;

    init_param();	      
    Nelem = N_signal;
    index = Nelem;

    NbrScale = ScaleNumber;
    TypeNoise = TNoise;
    Type_Transform = Trans;
    Set_Transform = SetTransform(Trans); // see MR1D_obj.CC
    FilterBank = FAS;
    TypeNorm = Norm;

    switch (Set_Transform)
    {
        case TRANS1_PAVE:   // see MR1D_obj.h
               Size = Nelem*(NbrScale-1);
               for (s = 0; s < NbrScale-1; s++)
               {
                   TabN[s] = Nelem;
                   TabPos[s] = s*Nelem;
               }
               break;
        case TRANS1_PYR:
               Size=0;
               TabN[0] = index;
               TabPos[0] = 0;
               for (s = 0; s < NbrScale-1; s++)
               {
                   Size += index;
                   TabPos[s+1] = TabPos[s]+index;
                   index = (index+1)/2;
                   TabN[s+1] = index;
               }
               break;
        case TRANS1_MALLAT:
	case TRANS1_WP_MALLAT:
	    Size = Nelem;
 	    for (s = 0; s < NbrScale-1; s++)
            {
                TabN[s] = index/2;
		index = (index+1)/2;
		TabPos[s] = index;
            }
	default: cerr << "Not implemented" << endl;
                 exit (0);
                 break;
    }

    if (one_level_per_pos(TNoise) == True) TabLevel = new float[Size];
    else TabLevel = new float[NbrScale];

    TabSupport = new unsigned char [Size];
    for (s = 0; s < Size; s++) TabSupport[s] = VAL_SupNull;

    CCD_Gain = 1.; 
    CCD_ReadOutSigma=0.; 
    CCD_ReadOutMean=0.; 
    SigmaNoise=0.; 
    for (s = 0; s < NbrScale-1; s++) NSigma[s] = DEF_NSIGMA_DETECT; 
    for (s = 0; s < NbrScale-1; s++) TabEps[s] = DEFAULT_EPSILON; 

    NSigma[0] = NSigma[0] + 1.;
    BadPixelVal=0.;
    BadPixel=False;
    MinEvent=False;
    SupIsol=False;
    FirstDectectScale = DEF_FIRST_DETECT_SCALE;
    OnlyPositivDetect=False;
    OnlyNegativDetect=False;  
    DilateSupport=False;
    MinEventNumber=DEF_MINEVENT_NUMBER;
    SigmaDetectionMethod = SIGMA_SUPPORT_1D;
    NiterSigmaClip = 1;
    SizeBlockSigmaNoise = DEFAULT_SIZE_BLOCK_SIG;
    TransImag=False;
    NewStatNoise = TNoise;
    CFewEventPoisson1d = NULL;
    TransImag=False;
    switch (TypeNoise)
    {
        case NOISE_GAUSSIAN: 
        case NOISE_NON_UNI_ADD:
        case NOISE_UNDEFINED:
        case NOISE_CORREL:    
               break;
        case NOISE_POISSON:
        case NOISE_GAUSS_POISSON:
        case NOISE_MULTI: 
        case NOISE_NON_UNI_MULT:
	       TransImag=True;
	       break;
        case NOISE_EVENT_POISSON:
               TransImag=True;
               MinEvent=True;
	       CFewEventPoisson1d = new FewEventPoisson();
	       CFewEventPoisson1d->compute_distribution(TraceFewEvent, E_1D);
               break;
        case NOISE_UNI_UNDEFINED: 
               NiterSigmaClip = 3;
               break;
        default:
	      cerr << "Error: this kind of noise is not treated in MR1D_NoiseModel ... " << endl;
	      exit(-1);
	      break;
    }
}

/****************************************************************************/

MR1DNoiseModel::MR1DNoiseModel (type_noise TNoise,  MR_1D &MR_Data)
{
    int N_Signal = MR_Data.size_ima_np();  // see MR1D_Obj.h
                // returns NP the number of points of the input signal

    alloc (TNoise, N_Signal, MR_Data.nbr_scale(), MR_Data.Type_Transform,
           MR_Data.filter_bank(), MR_Data.Norm);
    //MR1DNoiseModel(TNoise, N_Signal, 
    //             MR_Data.nbr_scale(), MR_Data.Type_Transform);
	// for nbr_scale() and Type_Transform  see MR1D_Obj.h
        //  returns Nbr_Plan
}

/****************************************************************************/

int MR1DNoiseModel::index(int s, int i) const
{
   int Ind=0;

   switch (Set_Transform)
   {
        case TRANS1_PAVE: // see MR1D_obj.h
		Ind = s*Nelem + i;
		break;
        case TRANS1_PYR: 
		Ind = TabPos[s]+ i; // bug + TabN[s] 
		break;
        default:
               cerr << "Error: MR1DNoiseModel, unknown set transform" << endl;
               exit(0);
        break;
   }
   return Ind;
}

/****************************************************************************/

Bool MR1DNoiseModel::operator() (int s,int i) const
{
   int Ind;
   Bool ValRet=False;

   Ind = index(s,i);
   
//if ((int)TabSupport[Ind] > 0)
//  cout << "TabSup(" << Ind << ")=" << (int)TabSupport[Ind] << " ,   "; 
   
   if ((TabSupport[Ind] > 0) && (TabSupport[Ind] <= VAL_SupLastOK)) 
                                                               ValRet = True;
   return ValRet;
}

/****************************************************************************/

Bool MR1DNoiseModel::operator() (int NumCoef) const
{
   Bool ValRet=False;
   int i,s, Ind;

   pos_mr1dcoeff(NumCoef, s, i, Set_Transform, Nelem,  NbrScale);
   Ind = index(s,i);
   if ((TabSupport[Ind] > 0) && (TabSupport[Ind] <= VAL_SupLastOK)) 
                                                          ValRet = True;
   return ValRet;
}

/****************************************************************************/

float & MR1DNoiseModel::sigma(int s,int i) const
{
   int Ind;

   if (one_level_per_pos(TypeNoise) == True) 
   {
       Ind = index(s, i);
   }
   else Ind = s;

   return TabLevel[Ind];
}

/****************************************************************************/

float & MR1DNoiseModel::sigma(int NumCoef) const
{
   int i, s, Ind;
   
   pos_mr1dcoeff(NumCoef, s, i, Set_Transform, Nelem,  NbrScale);

   if (one_level_per_pos(TypeNoise) == True) 
   {
       Ind = index(s, i);
   }
   else Ind = s;

   return TabLevel[Ind];
}

/****************************************************************************/

unsigned char & MR1DNoiseModel::support(int s,int i) const
{
   int Ind;

   Ind = index(s,i);

   return TabSupport[Ind];
}

/****************************************************************************/

unsigned char & MR1DNoiseModel::support(int NumCoef) const
{
   int i, s, Ind;
   pos_mr1dcoeff(NumCoef, s, i, Set_Transform, Nelem,  NbrScale);
   Ind = index(s,i);
   return TabSupport[Ind];
}

/****************************************************************************/

Bool MR1DNoiseModel::signif (float Val, int s, int i )
{
   Bool ValRet = False;
   float Level;

   Level = sigma(s,i)*NSigma[s];
   if (OnlyPositivDetect == True) 
   {
      if (Val > Level) ValRet = True;
   }
   else 
   {
       if (OnlyNegativDetect == True)
       {
           if (Val < -Level) ValRet = True;
       }
       else if (ABS(Val) > Level) ValRet = True;
   }
   if ((s < FirstDectectScale) && (ValRet == True)) ValRet = False;

   return ValRet;
}

/****************************************************************************/
 
float MR1DNoiseModel::prob_noise(float Val, int b, int i)  
{
    float Sig,P=0.;
    double Vn=0.;
    
    switch (TypeNoise)
    {
       case NOISE_EVENT_POISSON:
        {
// 	   if (NoCompDistrib) {
// 	        cout << "histogram have to be computed first" << endl;
// 	        exit(-1);
// 	   }
// 	   int k,l,Win = (int) (pow((double)2., (double)(b+2)) + 0.5);
//            int Nevent = 0;
//            for (k =i-Win; k <= i+Win; k++)
//                Nevent += Event_Image(k, l, Border);
//            //P = CEventPois->a_trou_repartition(Val, Nevent, b);
//            P = CFewEventPoisson2d->a_trou_repartition(Val, Nevent, b);
// 	   if (Val > 0) P = 1. - P;
           cout << "Error: not valid function for this noise model ... " << endl;
	   exit(-1);
        }  
	break;
      default:    
        Sig = sigma(b,i);
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

Bool MR1DNoiseModel::signif (float Val, int s, int i ,fltarray& Nsigma)
{
   Bool ValRet = False;
   float Level;

   Level = sigma(s,i)*Nsigma(s);
   if (OnlyPositivDetect == True) 
   {
      if (Val > Level) ValRet = True;
   }
   else 
   {
       if (OnlyNegativDetect == True)
       {
           if (Val < -Level) ValRet = True;
       }
       else if (ABS(Val) > Level) ValRet = True;
   }
   if ((s < FirstDectectScale) && (ValRet == True)) ValRet = False;

   return ValRet;
}
/****************************************************************************/

float MR1DNoiseModel::val_transform(float Val)
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
               if (Val > 0) ValRet = log(Val);
               else ValRet = 0.;
               break;
        case NOISE_EVENT_POISSON:
               break;
       default:
	      cerr << "Error: this kind of noise is not treated in MR1D_NoiseModel ... " << endl;
	      exit(-1);
	      break;    }
    return ValRet;
}

/****************************************************************************/

float MR1DNoiseModel::val_invtransform(float Val)
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
               if (Val > 0) ValRet = exp(Val);
               else ValRet = 0.;
               break;
        case NOISE_EVENT_POISSON:
               break;
       default:
	      cerr << "Error: this kind of noise is not treated in MR1D_NoiseModel ... " << endl;
	      exit(-1);
	      break;
    }
    return ValRet;
}

/****************************************************************************/

void MR1DNoiseModel::signal_transform (fltarray & Signal)
                        // apply a transform to a signal
{
    int i;
    switch (TypeNoise)
    {
        case NOISE_POISSON:
        case NOISE_GAUSS_POISSON:
               noise_poisson_transform (Signal, Signal); // see im_noise
               if (SigmaNoise < FLOAT_EPSILON) SigmaNoise = 1.;
               break;
        case NOISE_GAUSSIAN: 
        case NOISE_NON_UNI_ADD:
        case NOISE_UNDEFINED:
        case NOISE_UNI_UNDEFINED:
	case NOISE_CORREL:
               break;
        case NOISE_MULTI: 
        case NOISE_NON_UNI_MULT:
              for (i = 0; i < Signal.n_elem(); i++)
                         Signal(i) = val_transform(Signal(i));
               break;
        case NOISE_EVENT_POISSON:
	       // event_poisson_transform (Image, Event_Image);
	       if (CFewEventPoisson1d != (FewEventPoisson*)NULL)
	          building_signal_signal (Signal, EventSignal);
	       else {
	          cerr << "Pb in building_signal_signal call, NOISE_EVENT_POISSON" << endl;
                  exit(-1);
	       }
	       break; 
        default:
 	      cerr << "Error: this kind of noise is not treated in MR1D_NoiseModel ... " << endl;
	      exit(-1);
	      break;
     }
}


/****************************************************************************/

void MR1DNoiseModel::signal_invtransform (fltarray & Signal)
{
    int i;
    switch (TypeNoise)
    {
        case NOISE_POISSON:
        case NOISE_GAUSS_POISSON:
               noise_inverse_poisson_transform (Signal, Signal); // see im_noise
               break;
        case NOISE_GAUSSIAN: 
        case NOISE_NON_UNI_ADD:
        case NOISE_UNDEFINED:
        case NOISE_UNI_UNDEFINED:
	case NOISE_CORREL:
               break;
        case NOISE_MULTI: 
        case NOISE_NON_UNI_MULT:
              for (i = 0; i < Signal.n_elem(); i++)
                         Signal(i) = val_invtransform(Signal(i));
               break;
        case NOISE_EVENT_POISSON:
               break;
       default:
	      cerr << "Error: this kind of noise is not treated in MR1D_NoiseModel ... " << endl;
	      exit(-1);
	      break;
    }
}

/****************************************************************************/

void MR1DNoiseModel::max_sigma_in_scale(fltarray & Signal_Sigma)
{
   int s, i, Kb, Ke, k;
   int Nls=Nelem, Nlp=Nelem;
   int Ind, Indp,Step,SizeWave=2;
   float Max;
   //extern float mr_tab_noise(int s); // function which returns the noise
           // at scale s for an input signal which noise is 1.0
           //  see MR_Noise.cc
   // see also mr1d_noise_compute in mr1d_filtering.cc
   for (i=0; i< Nelem; i++)
   {
           Max = Signal_Sigma(i);
           Kb = (i-SizeWave >= 0) ? i-SizeWave : 0;
           Ke = (i+SizeWave <= Nls-1) ? i+SizeWave : Nelem-1;

           for (k=Kb; k<= Ke; k++)
                   if (Max < Signal_Sigma(k)) Max = Signal_Sigma(k);
           Ind = index(0,i);
           TabLevel[Ind] = Max*mr1d_tab_noise(0) / 0.972463;
    }

   Step=1;
   for (s=1; s < NbrScale-1; s++)
   {

       if (Set_Transform == TRANS1_PYR)
       {
           Nls = TabN[s];
           Nlp = TabN[s-1];
       }
       else 
       {
           Step = 2*Step;
       }
       for (i=0; i< Nls; i++)
       {
           Max = 0.;
           if (Set_Transform == TRANS1_PAVE) // see MR1D_obj.h
           {
               Kb = (i-Step >= 0) ? i-Step : i;
               Ke = (i+Step <= Nelem-1) ? i+Step : i;
               for (k=Kb; k<= Ke; k+=Step)
               {
                  Indp = index(s-1,k);
                  if (Max < TabLevel[Indp]) Max = TabLevel[Indp];
               }
           }
           else
           {
               Kb = (2*i-2 >= 0) ? 2*i-2 : 0;
               Ke = (2*i+2 <= Nlp-1) ? 2*i+2 : Nlp-1;
               for (k=Kb; k<= Ke; k++)
               {
                  Indp = index(s-1, k);
                  if (Max < TabLevel[Indp]) Max = TabLevel[Indp];
               }
           }
           Ind = index(s, i);
           TabLevel[Ind] = Max*mr1d_tab_noise(s)/mr1d_tab_noise(s-1);
       }
   }
}

/****************************************************************************/

void MR1DNoiseModel::dilate_support(int s)
{
   int i;
   int Nls = Nelem;

   if (Set_Transform == TRANS1_PYR) // see MR1D_obj.h
   {
      Nls = TabN[s];
   }

   for (i = 1; i < Nls-1; i++)
   {
      if (support(s, i) == VAL_SupNull)
      {
         if (support(s, i-1) == VAL_SupOK) support(s, i) = VAL_SupDill;
         if (support(s, i+1) == VAL_SupOK) support(s, i) = VAL_SupDill;
      }
   }
}

/****************************************************************************/

void MR1DNoiseModel::kill_isol(int s)
{
   int i;
   Bool PutZero;
   int Nls = Nelem;

   if (Set_Transform == TRANS1_PYR) // see MR1D_obj.h
   {
      Nls = TabN[s];
   }
   if (Set_Transform == TRANS1_PAVE)
   {
       for (i = 1; i < Nelem-1; i++)
       {
           PutZero = True; 
           if (support(s, i) == VAL_SupOK)
           {
               if (support(s, i-1) == VAL_SupOK) PutZero = False;
               else if (support(s, i+1) == VAL_SupOK) PutZero = False;
               else if (support(s+1,i) == VAL_SupOK) PutZero = False;
               if (PutZero == True) support(s,i) = VAL_SupKill;
           }
       }
   }
   else  if (Set_Transform == TRANS1_PYR) 
   {
       Nls = TabN[s];
       for (i = 1; i < Nls-1; i++)
       {
           PutZero = True;
           if (support(s, i) == VAL_SupOK)
           {
               if (support(s, i-1) == VAL_SupOK) PutZero = False;
               else if (support(s, i+1) == VAL_SupOK) PutZero = False;
               else if (support(s, i) == VAL_SupOK) PutZero = False;
               else if (support(s+1, i/2) == VAL_SupOK) PutZero = False;
               if (PutZero == True) support(s, i) = VAL_SupKill;
           }
       }
   }
}

/****************************************************************************/

void MR1DNoiseModel::kill_event(int s, fltarray &Event_Signal, int FirstWin)
//   FirstWin = 5  in MR1D_NoiseModel.h
{
    int i, k, Nls,sc=0;
    int Step=1;
    int Win = FirstWin/2;
    float Total;

    if ((Set_Transform == TRANS1_PYR) || (Set_Transform == TRANS1_PAVE))
            // see MR1D_obj.h
    {
        Nls = TabN[s];
        while (sc++ < s) Win*=2;

        Total = 0;
        for (k =0; k <= Win; k++) Total += Event_Signal (k);
        if ((support(s, 0) == VAL_SupOK) && (Total < MinEventNumber))
                                                support(s, 0) = VAL_SupMinEv;

        for (i =1; i <= Win; i+=Step)
        {
            Total = 0;
            Total += (int) Event_Signal(i+Win);
            if ((support(s, i) == VAL_SupOK) && (Total < MinEventNumber))
                                              support(s, i) = VAL_SupMinEv;
       }
        for (i = Win+1; i < Nelem - Win; i+=Step)
        {
            Total = 0;
            // kb = (i-Win < 0) ? 0: i-Win;  
            // ke = (i+Win >= Nls) ? Nls-1: i+Win;
            //Total -= (int) Event_Signal(i-Win-1,I_ZERO);
            //Total += (int) Event_Signal(i+Win,I_ZERO);
            Total -= (int) Event_Signal(i-Win-1);
            Total += (int) Event_Signal(i+Win);
            if ((support(s, i) == VAL_SupOK) && (Total < MinEventNumber))
                                              support(s, i) = VAL_SupMinEv;
        }
        for (i = Nelem - Win; i < Nelem; i+=Step)
        {
            Total = 0;
            Total -= (int) Event_Signal(i-Win-1);
            if ((support(s, i) == VAL_SupOK) && (Total < MinEventNumber))
                                              support(s, i) = VAL_SupMinEv;
       }
    }
}

/****************************************************************************/

float MR1DNoiseModel::multi_sure_estimation(MR_1D & MR_Data, int b)
{
   int i;
   int Ind = 1;
   int Np = MR_Data.size_scale_np(b);
   unsigned long N = (unsigned long) Np;
   float *Tab = new float [N+1];
   for (i = 0; i < MR_Data.size_scale_np(b) ; i++)
   {
       float Coef = MR_Data(b,i) /  sigma(b,i);
       Tab[Ind++] =  Coef*Coef;
   }
   sort(N, Tab);
   double MinRisk=0.;
   double CumVal = 0;
   int IndRisk=0;
   for (i = 1; i <= Np; i++) 
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

float MR1DNoiseModel::sure_estimation(MR_1D &MR_Data)
{
   int s,i;
   int Ind = 1;
   unsigned long N = 0;
   for (s = 0; s < MR_Data.nbr_band()-1; s++)
        N +=  MR_Data.size_scale_np(s);
   float *Tab = new float [N+1];

   for (s = 0; s < MR_Data.nbr_band()-1; s++)
   for (i = 0; i < MR_Data.size_scale_np(s) ; i++)
   {
       float Coef = MR_Data(s,i) / sigma(s,i);
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

/****************************************************************************/

void MR1DNoiseModel::set_support(MR_1D &MR_Data)
                    //  set the support using the noise modelization
                    //  and the MultiResolution data
{
    int i , s;
    int Nls;
    float SureLevel;
    
    switch (TypeThreshold)
    {
        case T_KSIGMA: break;
	case T_FDR:
	    for (s = 0; s < NbrScale-1; s++)
	    {
	       int Np = MR_Data.size_scale_np(s);
	       dblarray TabPValue(Np);
	       for (i = 0; i < Np;i++)
                   TabPValue(i) =  prob_noise(MR_Data(s,i), s, i);
	       float Alpha = 1. - erf(NSigma[s] / sqrt(2.));
	       if (Alpha > 0.5) Alpha = 0.5;  // Does not make any sense to tolerate more than 50% of false detections
	       double PDet = fdr_pvalue(TabPValue.buffer(), Np, Alpha);
	       for (i = 0; i < Np; i++)
 	       {
	          if (TabPValue(i) < PDet) 
		  {
		     support(s,i) = VAL_SupOK;
		     if  ((OnlyPositivDetect == True) && (MR_Data(s,i) < 0)) support(s,i) = VAL_SupNull;
		     if  ((OnlyNegativDetect == True) && (MR_Data(s,i) > 0)) support(s,i) = VAL_SupNull;
                     if (s < FirstDectectScale) support(s,i) = VAL_SupNull;
                  }
	          else support(s,i) = VAL_SupNull;
	       }
               TabEps[s] = PDet;
               // if (PDet < DOUBLE_EPSILON) NSigma[s] = 5.;
               // else   
	       NSigma[s] = ABS(xerfc(0.5+(1-PDet)/2.));
  	       // NSigma[s] = ABS(xerf(PDet/2.));
	    }
	    break;
	case T_UNIVERSAL:break;
	    for (s = 0; s < NbrScale-1; s++) 
	       NSigma[s] = sqrt(2.*log((float)(MR_Data.size_scale_np(s))));
	    break;
	case T_SURE:
	   SureLevel = sure_estimation(MR_Data);
	   for (s = 0; s < NbrScale-1; s++) NSigma[s] = SureLevel;
	   break;
	case T_MRSURE: 
	  for (s = 0; s < NbrScale-1; s++)  NSigma[s] = multi_sure_estimation(MR_Data,s);
	break;  
    }
    
    switch (Set_Transform)
    {
        case TRANS1_PAVE:  // see MR1D_obj.h
        case TRANS1_PYR:
          for (s = 0; s < NbrScale-1; s++)
          {
              Nls = MR_Data.size_scale_np(s);  
            // MR_Data.scale is  a fltarray of size Ns
                                       // given by size_scale_np
              for (i = 0; i < Nls;i++)
              {   // MR1DNoiseModel::signif (float Val, int s, int i )
                 if (signif(MR_Data(s, i), s, i) == True) 
                                            support(s, i) = VAL_SupOK;
                 else support(s, i) = VAL_SupNull;
              }
          }
          break;
        default:
         cerr << "Error in noise_mr_iter_threshold: bad Set_Transform"<< endl;
         exit (-1);
         break; 
    }
    if (SupIsol == True) for (s = 0; s < NbrScale-2; s++) kill_isol(s);
    if (DilateSupport == True) 
                    for (s = 0; s < NbrScale-1; s++) dilate_support(s);

}

/****************************************************************************/
Bool MR1DNoiseModel::kill_coef(int s, int i, float Val, Bool SetSupport)
{
    Bool ValRet=False;

    if ((*this)(s, i) == False)
    {
        if (SetSupport == False) ValRet = True;
        else
        {
            if ((support(s, i) == VAL_SupNull) &&
                (signif(Val, s, i) == True))
            {
               support(s, i) = VAL_SupOK;
            }
            else ValRet = True;
        }
    }
    return ValRet;
}    

/****************************************************************************/

void MR1DNoiseModel::threshold(MR_1D &MR_Data, Bool Smooth, Bool SetSupport)
// SetSupport is False
{
    int i,s;
    int Nls;

    switch (Set_Transform)
    {
        case TRANS1_PAVE: // see MR1D_obj.h
        case TRANS1_PYR:
         {
          fltarray Buff (Nelem);
 
          for (s = 0; s < NbrScale-1; s++)
          {
              Nls = TabN[s];
              //Buff.resize (Nls,Ncs);
              Buff.reform(Nls);  // see array.cc

              for (i = 0; i < Nls;i++)
              {                           
                 if (kill_coef(s, i, MR_Data(s, i), SetSupport) == True)
		    MR_Data(s, i) = 0.;
		 //else cout << "MR_Data(" << s << "," << i << ")="
		 //          << MR_Data(s, i) << endl;
              }

              if (Smooth == True)
              {
                  float Coef = 1. / (1. - (3./8.)*(3./8.));
                  float Ener = Coef;

                  if (Set_Transform == TRANS1_PYR)
                  {
                        cout << "Not yet implemented" << endl; 
                       // xxx RG smooth_bspline in IM_smooth.cc
                       // works with Ifloat and not with fltarray !!
                      // smooth_bspline(MR_Data.scale(s), Buff);
                      // Buff = MR_Data.scale(s) - Buff;

                      for (i = 0; i < Nls;i++)
                         if ( (*this)(s, i) == False) MR_Data(s, i) = Buff(i);
                  }
                  else 
                  {
                     fltarray BuffAux (Nelem);
                     // type_border Border = MR_Data.Border;  
                          // Border is public in mr1d_obj
                     //Buff = MR_Data.scale(s); 
		     MR_Data.scale(Buff, s);
                         // scale is a function which returns a fltarr(Ns)
                     for (i=0; i <= s; i++)
                     {
                        cout << "Not yet implemented" << endl;
                       // see above xx RG
                       // smooth_bspline (Buff, BuffAux, Border, i);
                        Buff = BuffAux;
                      //  if (i != s) Buff = BuffAux;
                      //  else Buff -= BuffAux;
                        Ener *= Coef;
                     }
                  }
                  for (i = 0; i < Nls;i++)
                    if ( (*this)(s, i) == False) MR_Data(s, i)=Buff(i)*Ener;
              } 
          }
         }
          break;
        default:
         cerr << "Error in noise_mr_iter_threshold: bad Set_Transform"<< endl;
         exit (-1);
         break; 
    }
}

/****************************************************************************/

void MR1DNoiseModel::refresh_support(MR_1D &MR_Data)
{
    int i, s;
    int ValSup;
    for (int c = 0; c < MR_Data.nbr_mr_coeff(); c++)
    {
        float CoefMr;

        pos_mr1dcoeff(c, s, i, Set_Transform, Nelem,  NbrScale);
        CoefMr = MR_Data(s, i); 
        ValSup = support(s,i);

        if ( (ValSup == VAL_SupNull)  && (signif(CoefMr, s, i) == True))
                                   support(s, i) = VAL_SupOK;
    }
}

/****************************************************************************/

void MR1DNoiseModel::weight_snr(MR_1D &MR_Data, Bool SetSupport)
{
    int i, s;
    float Weight;
    int ValSup;

    for (int c = 0; c < MR_Data.nbr_mr_coeff(); c++)
    {
        float CoefMr;

        pos_mr1dcoeff(c, s, i, Set_Transform, Nelem, NbrScale);
        CoefMr = MR_Data(s, i);
        ValSup = support(s, i);

        if (   (SetSupport == True) && (ValSup == VAL_SupNull) 
            && (signif(CoefMr, s, i) == True))
                                   support(s, i) = VAL_SupOK;

        Weight = ABS(CoefMr) / (NSigma[s]*sigma(s, i));
        if (Weight > 1.) Weight = 1.;
        MR_Data(s, i) *= Weight;
    }
}

/****************************************************************************/

void MR1DNoiseModel::weight_invsnr(MR_1D &MR_Data, Bool SetSupport)
{
    int i, s;
    float Weight;
    int ValSup;

    for (int c = 0; c < MR_Data.nbr_mr_coeff(); c++)
    {
        float CoefMr;
       
        pos_mr1dcoeff(c, s, i, Set_Transform, Nelem, NbrScale);
        // s and i output
        CoefMr = MR_Data(s, i);
        ValSup = support(s, i);

        if (   (SetSupport == True) && (ValSup == VAL_SupNull) 
            && (signif(CoefMr,s, i) == True))
                                   support(s, i) = VAL_SupOK;

        Weight = ABS(CoefMr) / (NSigma[s]*sigma(s, i));
        if (Weight > 1.) Weight = 1.;
        MR_Data(s, i) *= (1. - Weight);
    }
}

/****************************************************************************/

void MR1DNoiseModel::signal_sigma(fltarray &Data,// input : the signal
                         fltarray &Signal_Sigma, // output the sigma
                         int SizeBlock, //input parameter the size of the box
                         int Nit     // input parameter the number of iterations
                           )
// RG from im_sigma in im_noise.cc
// warning it can divide by 0 ! 
/* 
** computes the standard deviation in a block around each pixel.
** if Nit > 1, then compute the standard deviation by a k-sigma clipping,
** and give an estimation of the local gaussian noise.
*/
{
 int i, k, Nl, It;
  Nl = Data.n_elem();  // see Array.h
  int Dsize = SizeBlock / 2;
  float x, S0, S1, S2, Mean = 0., Sigma = 0., Sm = 0.;
  int Kb, Ke;

  for (i=0; i< Nl-1; i++)
  {
      for (It = 0; It < Nit; It++)
      {
         S0 = 0.;
         S1 = 0.;
         S2 = 0.;
         Kb = (i-Dsize >= 0) ? i-Dsize : 0;
         Ke = (i+Dsize <= Nl-1) ? i+Dsize : Nl-1;
         for (k=Kb; k<= Ke; k++)
         {
            x = Data(k);
            if ((It == 0) || (ABS(x - Mean) < Sm))
            {
               S0 ++;
               S1 += x;
               S2 += x*x;
            }
         }
         Mean = S1 / S0;
         Sigma = sqrt(S2/S0- Mean*Mean);
         Sm = 3. * Sigma;
      }
      Signal_Sigma(i) = Sigma;
  }  

}

/****************************************************************************/

void MR1DNoiseModel::get_sigma_from_rms_map(fltarray  &Dat_Sigma)
{
  //Bool ModifSize = False;
  int N = Dat_Sigma.nx();
  fltarray Dat_RMS(N);
  fltarray Dirac(N);
  fltarray Convolved(N);
  fltarray RmsResult(N);
  fltarray Band(N);
  int i,s;
  int Ind;
  
  Dirac(N/2-1)=1.;

  // First we compute the square image of Rms
  for(i=0; i < N; i++) Dat_RMS(i) = Dat_Sigma(i)*Dat_Sigma(i);
 
  // Then we transform the dirac image :
  MR_1D MR_Dirac(N, Type_Transform, "Dirac", NbrScale);
  MR_Dirac.transform(Dirac);

  // Then we compute the square of the Dirac transform ...
  for (s=0 ; s < NbrScale-1; s++)
  for (i=0 ; i < N; i++) MR_Dirac(s,i) *= MR_Dirac(s,i);

  // Then we convolve the RMS by the Dirac (at each scale)
  // The problem is now to take the sqrt of the result image. The FFT may 
  // produce some negative values. Since the result is lower for high number
  // of scales, we cannot use a test based on a comparision of FLOAT_EPSILON
  // or some power of it. We look simply if there was something at the place
  // on the original map ...

  for (s=0 ; s < NbrScale-1 ; s++)
  {
    FFTN_1D CFT;
    MR_Dirac.scale(Band, s);
    CFT.convolve(Dat_RMS, Band, RmsResult);
//     cout << " SCALE " << s+1 << endl;
//     cout << "MIN = " << Dat_RMS.min() <<  " MAX = " << Dat_RMS.max() <<endl;
//     cout << "MIN = " << Band.min() <<  " MAX = " << Band.max() <<endl;
//     cout << "MIN = " << RmsResult.min() <<  " MAX = " << RmsResult.max() <<endl;

    for (i=0 ; i < N ; i++)
    {
	Ind = index(s,i);
	if (Dat_Sigma(i) == 0) TabLevel[Ind] = 0;  
	else TabLevel[Ind] = sqrt(RmsResult(i));
    }
  }
}

/****************************************************************/

void MR1DNoiseModel::set_sigma(fltarray & Signal, MR_1D &MR_Data)
{
    int s,i;
    int Ind;

    // noise_compute (NbrScale, Transform, Nl, Nc);
    mr1d_noise_compute (NbrScale, Type_Transform,
                        MR_Data.Border, MR_Data.MedianWinSize, TypeNorm); 
                                            // RG in MR1D_filtering.cc
                                            // definition in  MR1D_filter.h

#if DEBUG_MODEL
cout << "set_sigma : Noise = " << SigmaNoise << endl;
#endif

    if ((SigmaNoise < FLOAT_EPSILON) && (NewStatNoise == NOISE_GAUSSIAN))
    {
        switch (SigmaDetectionMethod)
        {
            case SIGMA_MEDIAN_1D: 
                       SigmaNoise = detect1d_noise_from_med (Signal); 
                       break; // RG in MR1D_filtering.cc and MR1D_filter.h
            case SIGMA_BSPLINE_1D:
                       // SigmaNoise = detect_noise_from_bspline (Imag);
                        cout << "Not yet implemented" << endl; // xx RG
                       break;
            case SIGMA_CLIPIMA_1D:
                      // SigmaNoise = detect_noise_sigma (Imag); 
                        cout << "Not yet implemented" << endl; // xx RG
                       break;
            case SIGMA_SUPPORT_1D: 
                    SigmaNoise = mr1d_detect_noise(Signal);
                    break;
            default: cerr<< "Error: MR1DNoiseModel:model, unknown method" << endl;
                     exit(0);break;
        }
    }

#if DEBUG_MODEL
 cout << "After set_sigma : Noise = " << SigmaNoise << endl;
#endif

    if (NewStatNoise == NOISE_GAUSSIAN)   
    {
        for (s = 0; s < NbrScale-1; s++)
            TabLevel[s] = SigmaNoise*mr1d_tab_noise(s);
    }

    if (NewStatNoise == NOISE_NON_UNI_ADD)
    {
        if (UseRmsMap == False)
        {
	   RmsMap.alloc(Nelem);
           fltarray Buff(Nelem);
           filt1d_mediane (Signal, Buff, Nelem, 25);  // Window size = 25
           Buff = Signal - Buff;
           signal_sigma(Buff, RmsMap, SizeBlockSigmaNoise, NiterSigmaClip);
	   UseRmsMap=True;
 	}
	if  ((Type_Transform == TO1_PAVE_B3SPLINE) || (Type_Transform == TO1_PAVE_B3SPLINE_GEN2) ||   (Type_Transform == TU1_MALLAT))
	{
	    if (RmsMap.nx() != Signal.nx())
	    {
	       cout << " Error: rms data file has not the same size as the input data .. " << endl;
	       exit(-1);
	    }
	    get_sigma_from_rms_map(RmsMap);
	}
	else  max_sigma_in_scale(RmsMap);  
    }
    if (NewStatNoise == NOISE_UNI_UNDEFINED)
    {
       for (s = 0; s < NbrScale-1; s++)
       {
           // float MeanScale;
           float SigmaScale;
           // sigma_clip(MR_Data.scale(s), MeanScale, SigmaScale, NiterSigmaClip);
	   //MR_Data.scale(s).sigma_clip(MeanScale, SigmaScale, NiterSigmaClip);
	   // Phil 08/09/99 Modif function scale
	   fltarray prov;
	   MR_Data.scale (prov, s);
	   // prov.sigma_clip(MeanScale, SigmaScale, NiterSigmaClip);
	   SigmaScale = get_sigma_mad(prov.buffer(), prov.nx());
	   // sigma_clip in Array.cc and Array.h
           TabLevel[s] = SigmaScale;
       }
    }

    if (NewStatNoise == NOISE_UNDEFINED)
    {
       if ((Set_Transform != TRANS1_PAVE) && (Set_Transform != TRANS1_PYR))
       {
          cerr << "Error: this kind of transform cannot be used " << endl;
          cerr << "       with this noise model " << endl;
          exit(-1);
       }
       fltarray SignalSigmaNoise(Nelem);
       int BlockSize = SizeBlockSigmaNoise;
       for (s = 0; s < NbrScale-1; s++)
       {
           float SigmaScale, MeanScale;

           // sigma_clip(MR_Data.scale(s), MeanScale, SigmaScale, NiterSigmaClip);
	   //MR_Data.scale(s).sigma_clip(MeanScale, SigmaScale, NiterSigmaClip);
           // Phil 08/09/99 Modif function scale
	   fltarray prov;
	   MR_Data.scale(prov, s);
	   prov.sigma_clip(MeanScale, SigmaScale, NiterSigmaClip);
	   
	   SignalSigmaNoise.reform(TabN[s]);
	  
           //signal_sigma(MR_Data.scale(s), SignalSigmaNoise, BlockSize, 
           //                          NiterSigmaClip); // see above
           // Phil 08/09/99 Modif function scale
	   fltarray prov1;
           MR_Data.scale(prov, s);
           signal_sigma(prov, SignalSigmaNoise, BlockSize, NiterSigmaClip);   
	   
           for (i = 0; i < TabN[s];i++)
           {
               Ind = index(s, i);
               TabLevel[Ind] = MAX(SigmaScale, SignalSigmaNoise(i));
           }
           if (Set_Transform == TRANS1_PAVE) BlockSize *=2;  // see MR1D_obj.h
       }
    }
    if (NewStatNoise == NOISE_EVENT_POISSON) {
    
       
       intarray EventCount(Nelem);
       float alpha=1.;
       for (s = 0; s < MR_Data.nbr_band()-1; s++) {
       
          alpha=1.;
          for (int sc=0; sc < s; sc++) alpha *= 2.;
          event1d_one_scale(EventSignal, s, EventCount, MR_Data.Border);
          for (i=0;i<Nelem;i++) {
	     Ind = index(s,i);
	     TabLevel[Ind] = 
	        SIGMA_BSPLINE_WAVELET * sqrt((float) EventCount(i)) / alpha;
          }
      }    
   }
   
   // JLS 26/02/01   sigma = mad(scale)
   if (NewStatNoise == NOISE_CORREL)
   {
       for (s = 0; s < NbrScale-1; s++)
       {
	   fltarray prov;
           MR_Data.scale(prov, s); 
           TabLevel[s] = get_sigma_mad(prov.buffer(), prov.nx());
	   // cout << "Band " << s+1 << " Sig = " << TabLevel[s] << endl;	     
       }
   }
}

/****************************************************************************/

void MR1DNoiseModel::model(fltarray & Signal, MR_1D &MR_Data)
{
    NewStatNoise = TypeNoise;
    int s;
    
    
    // if TypeNoise != NOISE_EVENT_POISSON and MinEvent == True
    // the parameter Event_Image must be intialized before the 
    // call to model
    if ((MinEvent == True) && (TypeNoise != NOISE_EVENT_POISSON)) {
       if (  (TypeNoise == NOISE_EVENT_POISSON) 
             && (Signal.n_elem() != EventSignal.n_elem())) {
        cerr << "Error: MRNoiseModel::model : EventSignal is not correctly Initialized ..." << endl;
        exit(-1);
       }
    }    
    
    // if TypeNoise == NOISE_EVENT_POISSON and TransImag == False
    // the parameter Event_Image must be intialized before the 
    // call to model
    if (  (TypeNoise == NOISE_EVENT_POISSON) && (TransImag == False) &&
          (Signal.n_elem() != EventSignal.n_elem())) {
        cerr << "Error: MR1DNoiseModel::model : EventSignal is not correctly Initialized ..." << endl;
        exit(-1);
    }
    
    
    
    if (TransImag == True)
    {
    switch (TypeNoise)
    {
        case NOISE_POISSON:
        case NOISE_GAUSS_POISSON:
               signal_transform(Signal);  // see above for signal_transform
               NewStatNoise = NOISE_GAUSSIAN;
               break;
        case NOISE_GAUSSIAN: 
	case NOISE_CORREL:
               break;
        case NOISE_MULTI: 
               signal_transform(Signal);  // see above for signal_transform
               NewStatNoise = NOISE_GAUSSIAN;
               break;
        case NOISE_NON_UNI_ADD:
        case NOISE_UNI_UNDEFINED:
        case NOISE_UNDEFINED:
               break;
        case NOISE_NON_UNI_MULT:
               signal_transform(Signal);  // see above for signal_transform
               NewStatNoise = NOISE_NON_UNI_ADD;
               break;
        case NOISE_EVENT_POISSON:
               // building_imag_imag(Imag, Event_Image);
               // in Lib_mr2d/MR_Psup.cc and MR_Psupport.h
               //cout << "Not yet implemented, in progres..." << endl;  
               if (TransImag == True) signal_transform(Signal);
	       // in case of Poisson noise, the MIRROR border must be used
	       MR_Data.Border = I_MIRROR;
	       MinEvent=True;
               TransImag = False;
	       SigmaNoise = 1.;
	       break;
       default:
	      cerr << "Error: this kind of noise is not treated in MR1D_NoiseModel ... " << endl;
	      exit(-1);
	      break;    }
    }
    TypeNorm = MR_Data.Norm;
    FilterBank = MR_Data.filter_bank();
    MR_Data.transform(Signal);
    if (TypeNoise == NOISE_EVENT_POISSON) {
       mr1d_psupport(MR_Data, (*this), MR_Data.Border, TraceFewEvent);
       
       if (SupIsol == True ) 
          for (int s=0;s<NbrScale-2; s++) kill_isol(s);
   
       if (CFewEventPoisson1d->_OutForTest) {
          for (int s=0;s<NbrScale; s++) {
             int NbPoint = MR_Data.size_ima_np();
             fltarray support(NbPoint);
             char FileName[256];
             sprintf (FileName, "Sup_Scale%d", s);
             for (int i=0;i<NbPoint;i++) {
                support(i)=this->support(s,i);}
             fits_write_fltarr(FileName,support);
          }
          fits_write_fltarr("Mr1d_.fits", MR_Data.image());
       }
   
    } else {
       set_sigma(Signal, MR_Data);  // see above
       set_support(MR_Data);        // see above
    }

	 
    if (TransImag == True)
    {
    switch (TypeNoise)
    {
        case NOISE_POISSON:
        case NOISE_GAUSS_POISSON:
               signal_invtransform(Signal);
               break;
        case NOISE_GAUSSIAN:
	case  NOISE_CORREL:
               break;
        case NOISE_MULTI: 
               signal_invtransform(Signal);
               break;
        case NOISE_NON_UNI_MULT: 
               signal_invtransform(Signal);
               break;
        case NOISE_NON_UNI_ADD:
        case NOISE_UNDEFINED:
        case NOISE_UNI_UNDEFINED:
               break;
        case NOISE_EVENT_POISSON:
	            
	       break;               
       default:
	      cerr << "Error: this kind of noise is not treated in MR1D_NoiseModel ... " << endl;
	      exit(-1);
	      break;    }
    }
	 
    // verify if minimum number of count OK
    if (MinEvent == True) 
    {
        if (TypeNoise != NOISE_EVENT_POISSON) 
                 for (s=0; s < NbrScale-1; s++) kill_event(s, Signal);
    }   
}

/****************************************************************************/

void MR1DNoiseModel::model(fltarray & Signal)
{
//   Multiresolution MR_Data(Nelem, NbrScale, Transform, "MR1DNoiseModel");
//   MR_1D (int N, type_trans_1d T, char *Name,
//              int Nbr_Scale, Bool Inter, float S0,
//              float Nu0, int N_V)
   MR_1D MR_Data;
   MR_Data.alloc(Nelem,  Type_Transform, NbrScale, FilterBank, TypeNorm);
   model(Signal, MR_Data); // see above
}


/****************************************************************************/

MR1DNoiseModel::~MR1DNoiseModel()
{
   if (Size > 0)
   {
      delete[] TabLevel;
      delete[] TabSupport;
   }
}

/****************************************************************************/

void MR1DNoiseModel::mr_obj(MR_1D &MR_Data)
{
   int i, s;

    switch (Set_Transform)
    {
        case TRANS1_PAVE:
        case TRANS1_PYR:
          for (s = 0; s < NbrScale-1; s++)
          for (i = 0; i < TabN[s];i++)
              {
                 if ((*this)(s, i)  == True)  MR_Data(s, i) = 1.;
                 else MR_Data(s, i) = 0.;
              }
          break;
        default:
          cerr << "Error: Unknown set transform ... " << endl;
          exit(0);
          break;
    }
}

/****************************************************************************/

void MR1DNoiseModel::write_support_mr(char *FileName)
{
   //MR_1D MR_Data(Nelem,  Type_Transform, FileName, NbrScale);
 
    //(*this).mr_obj(MR_Data);
//    MR_Data.write(FileName);  see MR_io.cc
     cout << "Warning in  MR1DNoiseModel::write_support_mr: Not yet implemented" << endl; // xx RG

    // !!!!!!!!! provisoire !!!!!!!! to do.....
    fltarray prov(Size);
    for (int i=0;i<Size;i++) {
       //if (TabSupport[i] != 0) cout << "WIN !!!!!!!!!!!!" << endl;
       prov(i) = (float)TabSupport[i];
    }
    fits_write_fltarr (FileName, prov);

}

/****************************************************************************/

void MR1DNoiseModel::write_support_ima(char *FileName)
{
   MR_1D MR_Data(Nelem,  Type_Transform, FileName, NbrScale);
//   MultiResol MR_Data(Nl, Nc, NbrScale, Transform, "MR1DNoiseModel");
   (*this).mr_obj(MR_Data);
   float Coef;
   fltarray Signal(Nelem);
   int s, i;

   Signal.init();
    switch (Set_Transform)
    {
        case TRANS1_PAVE:
        case TRANS1_PYR:
          {
             fltarray Dat(Nelem);
             for (s = NbrScale-2; s >= 0; s--)
             {
                Coef = (float) ( 1 << (NbrScale-1-s));
                //im_block_extend(MR_Data.scale(s), Dat);
               // IM_SubImag.cc  IM_obj.h
     cout << "Not yet implemented" << endl; // xx RG

                for (i = 0; i < Dat.n_elem();i++)
                   if (Dat(i) > FLOAT_EPSILON) ;// xx rg Ima(i) = Coef;
             }
          }
          break;
        default:
          cerr << "Error: Unknown set transform ... " << endl;
          exit(0);
          break;

    }

//    io_write_ima_float(FileName, Ima);
     cout << "Not yet implemented" << endl; // xx RG
}

/****************************************************************************/



