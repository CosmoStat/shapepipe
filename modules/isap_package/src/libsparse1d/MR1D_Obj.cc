/*******************************************************************************
**
**    UNIT
** 
**    Version: 3.3
**
**    Author: Jean-Luc Starck
**
**    Date:  96/06/13 
**    
**    File:  MR1D_Obj.cc
**
**    Modification history:
**           23-OCT-1996 R Gastaud add pos_mr1dcoeff
**            24-OCT-1996 R Gastaud add nbr_mr_coeff
**
*******************************************************************************
**
**    DESCRIPTION  1D multiresolution transform
**    -----------  
**                 
*******************************************************************************
**
** int MR_1D::size_scale_np (int s) const 
** 
** return the number of points of  a scale
**
*******************************************************************************
**
** MR_1D::MR_1D (int N, int int Nbr_Scale,
**               type_trans_1d Transform, char *Name)
**
** Constructor: allocate the memory
** initialize the varaibles
**
** N = number of line of the signal
** Nbr_Scale = number of scales of the transform
** Transform = type of transform
** Name = name of the object (usefull for debugging)
** 
*******************************************************************************
**
** MR_1D::~MR_1D()
** 
** Destructor:
**
** deallocates the memory
**
*******************************************************************************
**
** void MR_1D::alloc (int N, int Nbr_Scale,
**                    type_trans_1d Transform, char *Name)
**
** allocate the memory
** initialize the varaibles
**
** N = number of points of the signal
** Nbr_Scale = number of scales of the transform
** Transform = type of transform
** Name = name of the object (usefull for debugging)
** 
*******************************************************************************
**
** void MR_1D::free ()
** 
** deallocates the memory
**
******************************************************************************/
 
#include "MR1D_Obj.h"

int  DEF_MedianWinSize = DEFAULT_MEDIAN_1D_WINSIZE;

/************************************************************************/

set_trans_1d which_set_is_trans1d(type_trans_1d type)
{
    switch (type)
    {
        case TU1_UNDECIMATED_NON_ORTHO:
        case TU1_MALLAT: 
        case TO1_PAVE_LINEAR: 
        case TO1_PAVE_B1SPLINE: 
        case TO1_PAVE_B3SPLINE:
        case TO1_PAVE_B3SPLINE_GEN2:
        case TO1_PAVE_MORLET:
        case TM1_PAVE_MEDIAN: 
        case TO1_PAVE_MEX: 
        case TO1_PAVE_FRENCH: 
	case TO1_PAVE_B3_DERIV:
	case TO1_PAVE_DERIV_GAUSS:
	case TO1_PAVE_HAAR:
              return TRANS1_PAVE;break;
        case TO1_PYR_LINEAR: 
        case TO1_PYR_B3SPLINE: 
        case TM1_PYR_MEDIAN: 
              return TRANS1_PYR;break;
        case TO1_MALLAT:
        case TO1_LIFTING:
	      return TRANS1_MALLAT; break;
	case WP1_MALLAT:
	case WP1_LIFTING:
	      return TRANS1_WP_MALLAT;
	      break;
	case WP1_ATROUS:
	      return TRANS1_WP_ATROUS;
	      break;
        case T1_UNDEFINED:
              return S1_UNDEFINED;break;
    }
    return S1_UNDEFINED;
}

const char * StringTransf1D (type_trans_1d type)
{
    switch (type)
    {
        case TU1_UNDECIMATED_NON_ORTHO:
              return ("Non orthogonal undecimated transform");break;
        case TU1_MALLAT:
	      return ("Undecimated (bi-) orthogonal wavelet transform");break;
        case TO1_PAVE_HAAR:
	      return ("Undecimated Haar wavelet transform: a trous algorithm");break;
        case TO1_PAVE_LINEAR: 
              return ("Linear wavelet transform: a trous algorithm");break;
        case TO1_PAVE_B1SPLINE: 
              return ("B1spline wavelet transform: a trous algorithm");break;
        case TO1_PAVE_B3SPLINE: 
              return ("B3spline wavelet transform: a trous algorithm");break;
        case TO1_PAVE_B3SPLINE_GEN2: 
            return ("Modified positive B3spline wavelet transform: a trous algorithm");break;
        case TO1_PAVE_MORLET:
              return ("Morlet\'s wavelet transform");break;
        case TM1_PAVE_MEDIAN: 
              return ("Morphological median transform");break;
        case TO1_PAVE_MEX: 
              return ("Mexican hat wavelet transform");break;
        case TO1_PAVE_FRENCH: 
              return ("French hat wavelet transform");break;
        case TO1_PYR_LINEAR: 
              return ("Pyramidal linear wavelet transform");break;
        case TO1_PYR_B3SPLINE: 
              return ("Pyramidal b3spline wavelet transform");break;
        case TM1_PYR_MEDIAN: 
              return ("Pyramidal median transform");break;
	case TO1_PAVE_B3_DERIV:
	      return ("Derivative of a b3spline: a trous algorithm");break;
	case TO1_PAVE_DERIV_GAUSS:   
	       return ("Gaussian Derivative wavelet transform");break;
        case TO1_MALLAT:
	       return ("(bi-) orthogonal wavelet transform ");break;
	       break;
        case TO1_LIFTING:
 	      return ("(bi-) orthogonal transform via lifting sheme"); 
	      break;
	case WP1_LIFTING:
 	      return ("Wavelet packets from lifting sheme"); 
	      break;
	case WP1_MALLAT:
	      return  ("Wavelet packets");
	      break;
	case WP1_ATROUS:      
	      return  ("Wavelet packets using the a-trous algorithm)");
	      break;
	case T1_UNDEFINED:
              return ("undefined transform");break;
    }
    return ("Error: bad type of transform");
}

/************************************************************************/


static float set_scale0 (type_trans_1d T, float Nu_0)
{
    float S=0.;

    switch (T)
    {
       case TO1_PAVE_MEX:
       case TO1_PAVE_DERIV_GAUSS:
          /* The minimum scale is such that the step size is just the place
             where the wavelet changes sign.
          */
            S =  1. / sqrt (3.); break;
       case TO1_PAVE_FRENCH:
           /* The minimum scale is such that we have just two points on each
              side of the central point for the computation of the integral.
              Remember that the French hat goes from -3 scale to 3 scale (2/3).
            */
            S =  .66; break;
       case TO1_PAVE_MORLET:
            S = 2. * Nu_0; break;
       case TO1_PAVE_LINEAR:
       case TO1_PAVE_HAAR:
       case TO1_PAVE_B1SPLINE:
       case TO1_PAVE_B3SPLINE:
        case TO1_PAVE_B3SPLINE_GEN2:
       case TM1_PAVE_MEDIAN:
       case TO1_PYR_LINEAR:
       case TO1_PYR_B3SPLINE:
       case TM1_PYR_MEDIAN:
       case TO1_PAVE_B3_DERIV:
       case TO1_MALLAT:
       case TU1_MALLAT:
       case TU1_UNDECIMATED_NON_ORTHO:
       case TO1_LIFTING:
       case WP1_LIFTING:
       case WP1_MALLAT:
       case WP1_ATROUS:
            S = 0.; break;
       case T1_UNDEFINED:
            cerr << "Error: undefined transform ..." << endl;
            exit(-1);
            break;
    }
    return S;
}
/************************************************************************/

static int set_scale_number (int N, type_trans_1d T, float Scale0, int N_V)
{
   int N_Scale=0;

   switch (T)
   {
     case TO1_PAVE_MEX:
     case TO1_PAVE_DERIV_GAUSS:
       /* N_Scale is the number of scales, i.e. the number of voices per
       octave multiplyed by the number of octaves. The number of octaves
       is the integral part of the binary logarithm of the ratio
       scale_max / scale_min. Since the extension of the
       wavelet is $8 scale$, we have $scale_{max} = n/8$.
       */
       N_Scale = iround ((float) N_V * log( (float) N / (8. * Scale0))  / log(2.));
       break;
     case TO1_PAVE_FRENCH:
       /* Since the extension of the
          wavelet is $6 scale$, we have scale_max = n/6.
       */
       N_Scale = iround((float) N_V  * log((float) N / (6. * Scale0)) / log(2.));
       break;
     case TO1_PAVE_MORLET:
       N_Scale = iround((float) N_V * log((float) N / (12. * Scale0)) / log(2.));
       break;
     case TO1_PAVE_LINEAR:
     case TO1_PAVE_HAAR:
     case TO1_PAVE_B1SPLINE:
     case TO1_PAVE_B3SPLINE:
     case TO1_PAVE_B3SPLINE_GEN2:
     case TM1_PAVE_MEDIAN:
     case TM1_PYR_MEDIAN:
     case TO1_PYR_B3SPLINE:
     case TO1_PYR_LINEAR:
     case TO1_PAVE_B3_DERIV:
     case TO1_MALLAT:
     case TU1_MALLAT:
     case TU1_UNDECIMATED_NON_ORTHO:
     case TO1_LIFTING:
     case WP1_MALLAT:
     case WP1_LIFTING:
     case WP1_ATROUS:
        /* The size of the wavelet for the i th scale is 2^{i+1}.
       We want to have a correct computation of the wavelet coefficients
       for at least a quater of the signal at the largest scale.
       Therefore we must have a number \verb|nx| of scales such that
       2^{nx+1} = 3 n / 4
       */
       N_Scale = iround((float)log((float) (N / 4. * 3.) / log(2.)));
       break;
       case T1_UNDEFINED:
            cerr << "Error: undefined transform ..." << endl;
            exit(-1);
            break;

   }
   return N_Scale;
}

/****************************************************************************/

static int set_voice_number (type_trans_1d T)
{
   int Nbr_V=0;
   switch (T)
   {
       case TO1_PAVE_LINEAR:
       case TO1_PAVE_HAAR:
       case TO1_PAVE_B1SPLINE:
       case TO1_PAVE_B3SPLINE:
       case TO1_PAVE_B3SPLINE_GEN2:
       case TM1_PAVE_MEDIAN:
       case TO1_PYR_LINEAR:
       case TO1_PYR_B3SPLINE:
       case TM1_PYR_MEDIAN:
       case TO1_PAVE_B3_DERIV:
       case TO1_MALLAT:
       case TU1_MALLAT:
       case TU1_UNDECIMATED_NON_ORTHO:
       case TO1_LIFTING:
       case WP1_MALLAT:
       case WP1_LIFTING:
       case WP1_ATROUS:
                 Nbr_V = 1; break;
       case TO1_PAVE_MORLET:
       case TO1_PAVE_MEX:
       case TO1_PAVE_FRENCH:
       case TO1_PAVE_DERIV_GAUSS:
              Nbr_V = 12;break;
       case T1_UNDEFINED:
            cerr << "Error: undefined transform ..." << endl;
            exit(-1);
            break;
   }
   return (Nbr_V);
}

/****************************************************************************/

set_trans_1d SetTransform (type_trans_1d Transform)
{
    return which_set_is_trans1d(Transform);
}

/************************************************************************/

int MR_1D::size_scale_np (int s) const {
   int Val_Return=0;
   int N=Data.axis(1);
   int i;

   switch (Set_Transform)
   {
        case TRANS1_PAVE:
	case TRANS1_WP_ATROUS:
               Val_Return = Data.axis(1);
               break;
        case TRANS1_PYR:
               for (i=0; i < s; i++) N = N / 2 + N % 2;
               Val_Return = N;
               break;
        case TRANS1_MALLAT:
	case TRANS1_WP_MALLAT:
 	      Val_Return = TabSize[s];
	      break;
	default:
               fprintf (stderr,"Error: unknown transform\n");
               exit (-1);
               break;
   }
   return (Val_Return);
}

/****************************************************************************/

MR_1D::MR_1D (int N, type_trans_1d T, char *Name,  
              int Nbr_Scale, Bool Inter, float S0, 
              float Nu0, int N_V)
{
    reset();
    Interp=Inter;
    alloc (N, T, Name,  Nbr_Scale, S0, Nu0, N_V);
}
 
/****************************************************************************/

void MR_1D::wp_pos(int N, int NStep, int & IndTab)
{
// cout << "wp_pos: NStep = " << NStep << "  N = " << N << " IndTab = " << IndTab << endl;
   if (NStep > 1)
   {
      wp_pos(N/2, NStep-1, IndTab);
      wp_pos((N+1)/2, NStep-1, IndTab);
   }
   else
   {
      TabSize[IndTab] = N/2;
      if (IndTab > 0)
        TabPos[IndTab] = TabPos[IndTab-1] - TabSize[IndTab];
      else TabPos[IndTab] = Np - N/2;
      IndTab++;
      TabSize[IndTab] = (N+1)/2;
      TabPos[IndTab] = TabPos[IndTab-1] - TabSize[IndTab];
      IndTab++;
   }
}
/****************************************************************************/

void MR_1D::reset()
{
   Np=0;
   Nbr_Plan=0;
   Type_Transform = T1_UNDEFINED;
   Set_Transform = S1_UNDEFINED; 
   Scale_0 = 0.;
   Interp=False;Nu_0 = 0.; 
   Nbr_Voie = 0; 
   FilterBankAlloc=False;
   Name_MR[0]='\0';FilterBank = NULL;
   SB_Filter = F_MALLAT_7_9; 
   Norm = NORM_L1;   
   U_Filter = DEF_UNDER_FILTER;  
   UndecFilterBank=NULL;
}

/****************************************************************************/
	  
void MR_1D::alloc (int N, type_trans_1d T, int Nbr_Scale, 
                   FilterAnaSynt *FAS, sb_type_norm TNorm,
                   Bool Inter, float S0, float Nu0, int N_V)
{
   Interp=Inter;
   char *Name = strdup("multi trans");
   FilterBank = FAS;
   if (FAS != NULL) SB_Filter = FAS->type_filter();
   Norm = TNorm;
   alloc (N, T, Name, Nbr_Scale, S0, Nu0, N_V);
}

/****************************************************************************/

void MR_1D::alloc (int N, type_trans_1d T, char *Name, 
              int Nbr_Scale, float S0, float Nu0, int N_V)
{
    int s;
    if (Nbr_Plan > 0) (*this).free();

    // FilterBank=NULL;
    // Norm = NORM_L1;
    FilterBankAlloc=False;
    Type_Transform = T;
    if ((T == TO1_MALLAT) || (T == TU1_MALLAT) || (T == WP1_MALLAT)) 
    {
       if (FilterBank == NULL)
       {
          SB_Filter = F_MALLAT_7_9;
          FilterBank = new FilterAnaSynt;
          FilterBank->alloc(SB_Filter);
          FilterBankAlloc=True;
       }
     }
     if (Type_Transform == TU1_UNDECIMATED_NON_ORTHO)
     {
         UndecFilterBank = new UndecSubBandFilter(U_Filter);
     }    
    Set_Transform = SetTransform (T);
    Np = N;
    MedianWinSize = DEF_MedianWinSize;
    Border = DEFAULT_BORDER_1D;
    KillBorder= False;
    LiftingTrans = DEF_LIFT;
    if (Nu0 < FLOAT_EPSILON)
         Nu_0 = DEFAULT_MORLET_NU_0;
    else Nu_0 = Nu0;
    if (S0 < FLOAT_EPSILON) 
         Scale_0 = set_scale0 (Type_Transform, Nu_0);
    else Scale_0 = S0;
    if (N_V < 1) 
         Nbr_Voie = set_voice_number (Type_Transform);
    else Nbr_Voie = N_V;
    if (Nbr_Scale < 2) 
         Nbr_Plan = set_scale_number (N, T, Scale_0, Nbr_Voie);
    else Nbr_Plan = Nbr_Scale;
    strcpy (Name_MR, Name);

     Nbr_Band = Nbr_Plan;
     if ((Set_Transform ==  TRANS1_WP_MALLAT) || 
 	 (Set_Transform ==  TRANS1_WP_ATROUS)) 
     {
	  Nbr_Band = (int) pow(2., (double) (Nbr_Plan-1));
     }
    TabPos = (int*)NULL;
    TabSize = (int*)NULL;
    if (Type_Transform == TO1_PAVE_MORLET) Data.alloc (N, 2*Nbr_Plan);
    else if ((Set_Transform ==  TRANS1_MALLAT) || 
             (Set_Transform ==  TRANS1_WP_MALLAT)) 
	 {
	     Data.alloc (N);        
	     TabPos = new int [Nbr_Band];
	     TabSize = new int [Nbr_Band];
	 }
    else Data.alloc (N, Nbr_Band);
    
    if (Set_Transform == TRANS1_MALLAT)  
    {
      int L = N;
      for (s=0; s < Nbr_Plan-1; s++)
      {
	 TabSize[s] = L / 2;
	 L = (L+1) / 2;
	 TabPos[s] = L;
       }
      s = Nbr_Plan-1;
      TabPos[s] = 0;
      TabSize[s] = L;
      //for (s=0; s < Nbr_Band; s++)
      //    cout << "Band " << s+1 << " pos = " << TabPos[s] << " " << TabSize[s] << endl;
    }
    if (Set_Transform == TRANS1_WP_MALLAT)
    {
      int Ind = 0;
       wp_pos(N, Nbr_Plan-1, Ind);
       //for (s=0; s < Nbr_Band; s++)
       //   cout << "Band " << s+1 << TabPos[s] << " " << TabSize[s] << endl;
    }
    
    switch (Type_Transform)
    {
        case TO1_PAVE_HAAR:
        case TO1_PAVE_LINEAR:  
        case TO1_PAVE_B1SPLINE: 
	case TO1_PYR_LINEAR: 
	         BorderSize_0 = 1;break;
        case TO1_PAVE_B3SPLINE: 
    case TO1_PAVE_B3SPLINE_GEN2: 
	case TO1_PYR_B3SPLINE: 
	case WP1_ATROUS:
	case TO1_MALLAT:
	case TU1_MALLAT:
	case TU1_UNDECIMATED_NON_ORTHO:
	case TO1_LIFTING:
	case WP1_MALLAT:
	case WP1_LIFTING:
	         BorderSize_0 = 2;break;
	case TO1_PAVE_B3_DERIV : 
	         BorderSize_0 = 3;break;
        case TM1_PAVE_MEDIAN: 
	case TM1_PYR_MEDIAN:
	         BorderSize_0 = DEFAULT_MEDIAN_1D_WINSIZE/2;break;
        case TO1_PAVE_MORLET:
	         BorderSize_0 =  6*Scale_0;
        case TO1_PAVE_MEX:
 	case TO1_PAVE_DERIV_GAUSS: 
 	         BorderSize_0 = 2*Scale_0;
	         break;
 	case  TO1_PAVE_FRENCH:
	          BorderSize_0 = 3*Scale_0;
                 break;
        default:
               fprintf (stderr, "Bug in SetTransform: bad parameter Transform\n");
	       exit(-1);
               break;
    }
}
 
/****************************************************************************/

MR_1D::~MR_1D()
{
   free ();
}

/****************************************************************************/

void MR_1D::free ()
{
    if ((FilterBankAlloc == True) && (FilterBank != NULL)) 
    {
        delete FilterBank;
        FilterBank = NULL;
    }
    if (UndecFilterBank != NULL)
    {
       delete UndecFilterBank;
       UndecFilterBank = NULL;
    }
    MedianWinSize = DEFAULT_MEDIAN_1D_WINSIZE;
    Border = DEFAULT_BORDER_1D;
    Data.free();
    Nbr_Plan = 0;
    Nbr_Band = 0;
    Np = 0;
    Scale_0 = 0.;
    Nu_0 = 0.; 
    Nbr_Voie = 0;
    Interp = False;
    Type_Transform=T1_UNDEFINED;
    Name_MR[0]='\0';
    if ((Set_Transform ==  TRANS1_MALLAT) || 
        (Set_Transform ==  TRANS1_WP_MALLAT)) {
	delete[] TabPos;
        delete[] TabSize;
    }   
    Set_Transform= S1_UNDEFINED;
}

/****************************************************************************/

float & MR_1D::operator() (int s, int i) const 
{
   if ((Set_Transform == TRANS1_MALLAT) || 
       (Set_Transform == TRANS1_WP_MALLAT))
   {
        return  Data(TabPos[s] + i);
   }
   else return Data(i,s);
}

/****************************************************************************/

//  fltarray MR_1D::scale(int s)
//  {
//     int Ns = size_scale_np(s);
//     fltarray *Result = new fltarray(Ns);
//     int j;
//     for (j=0; j < Ns; j++)  (*Result)(j) = (*this)(s,j);
//     return (*Result);
// } 

/****************************************************************************/

void MR_1D::scale(fltarray & TabScale, int s)
{
   int j,Ns = size_scale_np(s);
   TabScale.alloc(Ns);
   for (j=0; j < Ns; j++)  TabScale(j) = (*this)(s,j);
}

/****************************************************************************/

int MR_1D::border_size(int s)
{
  int ValRet = 0;
  
  if ((Set_Transform == TRANS1_PYR) ||
       (Set_Transform == TRANS1_MALLAT) || 
       (Set_Transform == TRANS1_WP_MALLAT))  ValRet = (int) BorderSize_0;
  else
  {
    if (Nbr_Voie == 1)
       ValRet = (int) (BorderSize_0 * pow((double)2., (double) s) + 0.5);
    else
    {    
     float Scale = pow ((double)2., (double)((double)s / (double) Nbr_Voie));
     ValRet = (int) (Scale*BorderSize_0+1);
    }
  }
  return ValRet;
}

/****************************************************************************/

void MR_1D::kill_border()
{
   int s,i,b,Np;
   for (s = 0; s < Nbr_Plan; s++)
   {
      b = border_size(s);
      Np = size_scale_np(s);
      for (i=0; i < b; i++)
      {
         if (b < Np) (*this)(s,i) =0.;
	 if (Np-i-1 > 0) (*this)(s, Np-i-1) =0.;
      }
   }
}

/****************************************************************************/

void MR_1D::modulus_maxima()
{
   int s,i,Np;
   for (s = 0; s < Nbr_Plan-1; s++)
   {
      Np = size_scale_np(s);
      (*this)(s,0) = (*this)(s,Np-1)=0;
      for (i=1; i < Np-1; i++)
      {
         if ((ABS((*this)(s,i)) < ABS((*this)(s,i-1))) ||
	      (ABS((*this)(s,i)) < ABS((*this)(s,i+1)))) (*this)(s,i) = 0.;
      }
   }
}
/****************************************************************************/

void MR_1D::transform (fltarray &Signal, type_border Bord)
{
   if (Nbr_Plan < 2)
   {
      cerr << "Error: Object not correctly allocated: Nbr_Plan must be > 2  ";
      cerr << endl;
      exit (-1);
   }

    switch (Type_Transform)
    {
       case TU1_UNDECIMATED_NON_ORTHO:
        {
           UndecSubBandFilter *USF = undec_filter_bank();
           if (USF == NULL)
           {
                  cout << "Error: undecimated filter bank is not defined ... " << endl;
                  exit(-1);
           }
	   PAVE_1D_WT UWT(*USF);            
	   UWT.transform(Signal, Data, Nbr_Plan);
	}
	break;
       case TU1_MALLAT:
        {
            FilterAnaSynt *FAS = filter_bank();
             if (FAS == NULL)
              {
                  cout << "Error: filter bank is not defined ... " << endl;
                  exit(-1);
              }
            SubBandFilter  SBF(*FAS, Norm);
            SBF.Border = Border;
            PAVE_1D_WT UWT(SBF);
            UWT.transform(Signal, Data, Nbr_Plan );
 	 }
	 break;
       case TO1_MALLAT:
        {
           FilterAnaSynt *FAS = filter_bank();
             if (FAS == NULL)
              {
                  cout << "Error: filter bank is not defined ... " << endl;
                  exit(-1);
              }
           SubBandFilter  SBF(*FAS, Norm);
           SBF.Border = Border;
           // SubBandFilter SBF(SB_Filter, Norm);
           mallat_1d_transform (Signal, Data, Nbr_Plan, SBF);
 	 }
	 break;
       case TO1_LIFTING:
        {
          Lifting SBF(LiftingTrans);
           SBF.Border = Border;
          mallat_1d_transform (Signal, Data, Nbr_Plan, SBF);
 	 }
         // lift1_1d_transform (Signal, Data, Nbr_Plan);
	 break; 
       case WP1_LIFTING:
        {
          Lifting SBF(LiftingTrans);
           SBF.Border = Border;
          wp1d_mallat_transform (Signal, Data, Nbr_Plan-1, SBF);
 	 }
	 break;
       case WP1_MALLAT:
        {
           FilterAnaSynt *FAS =  filter_bank();
             if (FAS == NULL)
              {
                  cout << "Error: filter bank is not defined ... " << endl;
                  exit(-1);
              }
           SubBandFilter  SBF(*FAS, Norm);
           SBF.Border = Border;
	   // SubBandFilter SBF(SB_Filter, Norm);
           wp1d_mallat_transform (Signal, Data, Nbr_Plan-1, SBF);
	}
	 break;
       case WP1_ATROUS:
         wp1d_atrous_transform (Signal, Data,Nbr_Plan-1, Bord);
	 break;
       case TO1_PAVE_HAAR:
         wave_1d_haar (Signal, Data, Np, Nbr_Plan, Bord); 
         break;
        case TO1_PAVE_LINEAR:
         wave_1d_linear (Signal, Data, Np, Nbr_Plan, Bord); 
         break;
       case TO1_PAVE_B1SPLINE:
         wave_1d_spline1 (Signal, Data, Np, Nbr_Plan, Bord);
         // mr1d_medianspline3_bis(Signal,Data,Nbr_Plan,MedianWinSize, Bord);
         break;
       case TO1_PAVE_B3SPLINE:
         wave_1d_spline3 (Signal, Data, Np, Nbr_Plan,Bord);
         break;
        case TO1_PAVE_B3SPLINE_GEN2:
            wave_1d_spline3_gen2 (Signal, Data, Np, Nbr_Plan,Bord);
            break;
       case TO1_PAVE_B3_DERIV:
         wave_1d_B3deriv_atrou (Signal, Data, Nbr_Plan, Bord);
         break;
       case TO1_PAVE_MORLET:
         wave_1d_morlet (Signal, Data, Np, Bord,
                         Nbr_Voie, Nbr_Plan, Nu_0, Scale_0);
         break;
       case TO1_PAVE_MEX:
         wave_1d_mex (Signal, Data, Np, Bord, Nbr_Voie, Nbr_Plan, Scale_0);
         break;
	case TO1_PAVE_DERIV_GAUSS:
	 wave_1d_der_gauss (Signal, Data, Np, Bord,
	                    Nbr_Voie, Nbr_Plan, Scale_0);
	 break;
       case TO1_PAVE_FRENCH:
         wave_1d_french (Signal, Data, Np, Bord,Nbr_Voie, Nbr_Plan, Scale_0);
         break;
       case TM1_PAVE_MEDIAN:
         mr1d_median (Signal, Data, Np, Nbr_Plan, MedianWinSize, Bord);
         break;
       case TO1_PYR_LINEAR:
         pyr_1d_linear (Signal, Data, Np, Nbr_Plan, Bord);
         break;
       case TO1_PYR_B3SPLINE:
         pyr_1d_spline3 (Signal, Data, Np, Nbr_Plan, Bord);
         break;
       case TM1_PYR_MEDIAN:
         mr1d_pyr_median (Signal, Data, Np, Nbr_Plan, MedianWinSize, Bord);
         break;
       case T1_UNDEFINED:
         cerr << "Error: undefined transform ..." << endl;
         exit(-1);
         break;
   }
   
   if ((Interp==True) && (which_set_is_trans1d(Type_Transform) == TRANS1_PYR))
   {
       int Iter,i,j;
       fltarray Resi (Np);
       fltarray MR_Iter (Np, Nbr_Plan);
       
        for (Iter = 0; Iter < Nbr_Plan; Iter++)
        {  
            for (i = 0; i < Nbr_Plan; i++)
            for (j = 0; j < Np; j++)  MR_Iter(j,i) = 0.;

            for (i = 0; i < Np; i++)  MR_Iter(i,Iter) = Data(i,Iter);
                mr1d_pyr_rec (Resi, MR_Iter, Np, Nbr_Plan);
            for (i = 0; i < Np; i++)  Data(i,Iter) = Resi(i);
        }
  }
   
  if (KillBorder == True) kill_border();
}

/****************************************************************************/

void MR_1D::transform (fltarray &Signal)
{
   this->transform(Signal, Border);
}

/****************************************************************************/

void MR_1D::recons (fltarray &Signal, type_border Bord)
{
    switch (Type_Transform)
    {
         case TU1_UNDECIMATED_NON_ORTHO:
	  {
	     UndecSubBandFilter *USF =  undec_filter_bank();
	     if (USF == NULL)
             {
                  cout << "Error: undecimated filter bank is not defined ... " << endl;
                  exit(-1);
             }
	     PAVE_1D_WT UWT(*USF);
             UWT.recons(Data, Signal, Nbr_Plan);
	  }
	  break;
         case TU1_MALLAT:
         {
           FilterAnaSynt *FAS = filter_bank();
             if (FAS == NULL)
              {
                  cout << "Error: filter bank is not defined ... " << endl;
                  exit(-1);
              }
           SubBandFilter  SBF(*FAS, Norm);
           SBF.Border = Border;
           PAVE_1D_WT UWT(SBF);
           UWT.recons(Data, Signal, Nbr_Plan);
 	 }
	 break;
       case TO1_PAVE_FRENCH: 
         wave_1d_french_rec (Data, Signal, Np, Bord,
                             Nbr_Voie, Nbr_Plan, Scale_0); 
         break;
       case TO1_PAVE_MEX: 
         wave_1d_mex_rec (Data, Signal, Np, Bord,
                             Nbr_Voie, Nbr_Plan, Scale_0);
         break;
       case TO1_PAVE_B3SPLINE_GEN2 :
            wave_1d_spline3_gen2_rec (Data, Signal, Np, Nbr_Plan);
            break;
       case TO1_PAVE_LINEAR:
       case TO1_PAVE_HAAR:
       case TO1_PAVE_B1SPLINE :
       case TO1_PAVE_B3SPLINE :
       case TM1_PAVE_MEDIAN:
         wave_1d_algo_trou_rec (Data, Signal, Np, Nbr_Plan);
	 break;
       case WP1_ATROUS:       
         wave_1d_algo_trou_rec (Data, Signal, Np, Nbr_Band);
         break;
       case TO1_PAVE_B3_DERIV:
         wave_1d_rec_B3deriv_atrou (Data, Signal, Nbr_Plan, Bord);
	 break;
       case TO1_PAVE_MORLET:
       case TO1_PAVE_DERIV_GAUSS:
         cerr << "Error: ";
         cerr << "This reconstruction is not implemented" << endl;
         exit (-1);
         break;
       case TO1_PYR_LINEAR:
       case TO1_PYR_B3SPLINE:
       case TM1_PYR_MEDIAN:
           mr1d_pyr_rec (Signal, Data, Np, Nbr_Plan);
         break;
       case TO1_MALLAT:
        {
           FilterAnaSynt *FAS =  filter_bank();
             if (FAS == NULL)
              {
                  cout << "Error: filter bank is not defined ... " << endl;
                  exit(-1);
              }
           SubBandFilter  SBF(*FAS, Norm);
           SBF.Border = Border;
           // SubBandFilter SBF(SB_Filter, Norm);
 	 mallat_1d_reconstruct (Data, Signal, Nbr_Plan, SBF);
 	 }
         break;
       case TO1_LIFTING:
         {
         Lifting SBF(LiftingTrans);
         SBF.Border = Border;
  	 mallat_1d_reconstruct (Data, Signal, Nbr_Plan, SBF);
 	 }
           // lift1_1d_reconstruct(Data, Signal, Nbr_Plan);
         break;
       case WP1_LIFTING:
         {
         Lifting SBF(LiftingTrans);
          SBF.Border = Border;
 	 wp1d_mallat_rec  (Data, Signal, Nbr_Plan-1, SBF);
 	 }
	 break;
       case WP1_MALLAT:
         {
           FilterAnaSynt *FAS =  filter_bank();
             if (FAS == NULL)
              {
                  cout << "Error: filter bank is not defined ... " << endl;
                  exit(-1);
              }
           SubBandFilter  SBF(*FAS, Norm);
           SBF.Border = Border;
         // SubBandFilter SBF(SB_Filter, Norm);
 	 wp1d_mallat_rec (Data, Signal, Nbr_Plan-1, SBF);
 	 }
         break;
       case T1_UNDEFINED:
         cerr << "Error: undefined transform ..." << endl;
         exit(-1);
         break;
   }

}

/****************************************************************************/

void MR_1D::recons (fltarray &Signal)
{
   this->recons(Signal, Border);
}

/****************************************************************************/

void MR_1D::rec_adjoint(fltarray &Signal,Bool UseLastScale, type_border Bord)
{
  MR_1D TempMR(Np, Type_Transform, "Adjoint rec", Nbr_Plan);
  int s,i;
  TempMR.Border = Bord;
  fltarray Tab;

  if (UseLastScale == False)
  {
     s = Nbr_Band-1;
     for (i=0; i < (*this).size_scale_np(s); i++) (*this)(s,i) = 0.;
  }

  switch (Type_Transform)
  {
       case TO1_PAVE_FRENCH: 
       case TO1_PAVE_MEX: 
       case TO1_PAVE_LINEAR:
       case TO1_PAVE_HAAR:
       case TO1_PAVE_B1SPLINE :
       case TO1_PAVE_B3SPLINE :
      case TO1_PAVE_B3SPLINE_GEN2:
       case TM1_PAVE_MEDIAN:
       case TO1_PAVE_B3_DERIV:
         for (s=0; s < Nbr_Plan-1; s++)
         {
	    int Nps = size_scale_np(s);
            // pb memory de-alloc: Tab = (*this).scale(s);
	    Tab.reform(Nps);
	    for (i=0; i < Nps; i++) Tab(i) = (*this)(s,i);
	                          
            TempMR.transform(Tab);
            for (i=0; i < Nps; i++) (*this)(s,i) = TempMR(s,i);
         } 
         (*this).recons(Signal); 
         break;
       case TO1_PAVE_DERIV_GAUSS:
       case TO1_PAVE_MORLET:
         cerr << "Error: ";
         cerr << "This reconstruction is not implemented" << endl;
         exit (-1);
         break;
       case TO1_PYR_LINEAR:
       case TO1_PYR_B3SPLINE:
       case TM1_PYR_MEDIAN:
       case TO1_MALLAT:
       case TU1_MALLAT:
       case TU1_UNDECIMATED_NON_ORTHO:
       case TO1_LIFTING:
       case WP1_MALLAT:
       case WP1_LIFTING:
       case WP1_ATROUS: (*this).recons(Signal);
         break;
       case T1_UNDEFINED:
         cerr << "Error: undefined transform ..." << endl;
         exit(-1);
         break;
   }
}

/****************************************************************************/

void pos_mr1dcoeff(int NumCoef, // input: coefficient of the wavelet transform
                   int &s, // output:  scale in the transform
                   int &i, // output: index in the image s of the transform
                   set_trans_1d Set_Transform,//input: which wavelet transform
                   int Nelem, // input: number of element of the input signal
                   int Nbr_Plan // input: number of scales in the transform 
                                //       (useful for Haar transform)
                   )

{
   int Val, my_index;
   
   switch (Set_Transform)
   {
        case TRANS1_PAVE:
	case TRANS1_WP_ATROUS:
               s = NumCoef / (Nelem);
               i = NumCoef - s * Nelem;
               break;
	case TRANS1_MALLAT:
	case TRANS1_WP_MALLAT:
	      cerr << "Error in pos_mr1d_coeff: not implemented transform ..." << endl;
	      exit(-1);
	      break;
        case TRANS1_PYR:
               Val = NumCoef ;
               my_index = Nelem;
               s = 0;	
               while (Val > my_index)
               {
                   s++;
                   Val -= Nelem;
                   my_index = (my_index+1)/2;
               }
               i = my_index;
               break;
        default:
               fprintf (stderr,"Error: unknown transform\n");
               exit (-1);
        break;
   }
   if (s >= Nbr_Plan)
   {
      fprintf (stderr,"Error: NumCoef too large ... \n");
      exit (-1);
   }
}
  
/****************************************************************************/

int MR_1D::nbr_mr_coeff ()
{
    int NbrCoef=0;
    int s, Nl_s = Np; // size of the input signal

    switch (Set_Transform)
    {
        case TRANS1_PAVE:
	case TRANS1_WP_ATROUS:
               NbrCoef = (Nbr_Band-1)*Nl_s;
               break;
	case TRANS1_MALLAT:
	case TRANS1_WP_MALLAT:
	       NbrCoef = Np - TabSize[Nbr_Band-1];
	       break;
        case TRANS1_PYR:
               for (s = 0; s < Nbr_Plan-1; s++)
               {
                   NbrCoef += Nl_s;
                   Nl_s = (Nl_s+1)/2;
               }
               break;               
        default:
               fprintf (stderr, "Not implemented\n");
               exit (0);
               break;
    }
    return NbrCoef;
}

/****************************************************************************/

void MR_1D::loc_maxima()
{
   int s,i,Np;
   for (s = 0; s < Nbr_Plan-1; s++)
   {
      Np = size_scale_np(s);

      intarray ao_NonMax (Np);
      int ind=0;

      for (i=1; i < Np-1; i++)
      {
         if (    ((*this)(s,i) > (*this)(s,i-1)) 
              && ((*this)(s,i) > (*this)(s,i+1)) ) {}
         else
			ao_NonMax(ind++) = i;
      }
      (*this)(s,0) = (*this)(s,Np-1) = 0.;
      for (i=0; i<ind; i++) (*this)(s,ao_NonMax(i)) = -0.;
   }
}
 
 
// 0 and Np-1 are not touched
void MR_1D::loc_optima()
{
   int s,i,Np;
   for (s = 0; s < Nbr_Plan-1; s++)
   {
      Np = size_scale_np(s);

      intarray ao_NonOpt (Np);
      int ind=0;

      for (i=1; i < Np-1; i++)
      {
         if (  (   ((*this)(s,i) > (*this)(s,i-1)) 
                && ((*this)(s,i) > (*this)(s,i+1)))
	     ||(   ((*this)(s,i) >= (*this)(s,i-1)) 
                && ((*this)(s,i) > (*this)(s,i+1)))
	     ||(   ((*this)(s,i) > (*this)(s,i-1)) 
                && ((*this)(s,i) >= (*this)(s,i+1)))
		
		
	     ||(   ((*this)(s,i) <= (*this)(s,i-1)) 
                && ((*this)(s,i) < (*this)(s,i+1)))
	     ||(   ((*this)(s,i) < (*this)(s,i-1)) 
                && ((*this)(s,i) <= (*this)(s,i+1)))						
             ||(   ((*this)(s,i) < (*this)(s,i-1)) 
                && ((*this)(s,i) < (*this)(s,i+1))) ) {}
         else
			ao_NonOpt(ind++) = i;
      }
      //(*this)(s,0) = (*this)(s,Np-1) = 0.;
      for (i=0; i<ind; i++) (*this)(s,ao_NonOpt(i)) = -0.;
   }
}
