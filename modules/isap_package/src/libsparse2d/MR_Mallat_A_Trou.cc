/*******************************************************************************
**
**    UNIT
**
**    Version: 
**
**    Author: 
**
**    Date:   
**    
**    File:  
**
*******************************************************************************
**    DESCRIPTION 
*******************************************************************************
*****************************************************************************/ 

// static char sccsid[] = "@(#)MR_Mallat.cc 3.1 96/05/02 CEA 1998 @(#)";

 
#include "MR_Obj.h"
#include "MR1D_Obj.h"   //a enlever, passer border_ind_test en communs a mresol
#include "IM_Edge.h"

#include "macro.h"

#define SIZE_FILTER 7
#define EPS 1e-07

// filter coef for foward 2-D wavelet transform
static double std_Hn[SIZE_FILTER] = {0, 0, 0.125, 0.375, 0.375, 0.125, 0};
static double std_Gn[SIZE_FILTER] = {0, 0, 0, -2, 2, 0, 0};

// filter coef for inverse 2-D wavelet transform
static double std_Kn[SIZE_FILTER] = {0.0078125, 0.0546875, 0.171875, -0.171875,\
                                     -0.0546875, -0.0078125, 0};
static double std_Ln[SIZE_FILTER] = {0.0078125, 0.046875, 0.1171875, 0.65625, \
                                     0.1171875, 0.046875, 0.0078125};
static double std_HnMod[SIZE_FILTER] = {0, 0.125, 0.375, 0.375, 0.125, 0, 0};



/*static void convol1d (fltarray &Signal, 
                      double *Filter, 
                      int SizeFilter,
                      fltarray &Result, 
                      int Step, 
                      type_border Border) {
   int N = Signal.nx();
   int i,j,Ind;
   int DemiSize = SizeFilter/2;
   double Val;
   for (i = 0; i < N; i ++)
   {
      Val = 0;
      for (j= -DemiSize; j <= DemiSize; j++)
      {
         Ind = border_ind_test (i+j*Step, N, Border);
         Val += (double) Signal(Ind) * Filter[j+DemiSize];
      }
      Result(i) = (float) Val;
   }
}*/


static void convolx1d ( Ifloat &pro_ImagIn, 
                        double* ppd_Filter,  
                        int pi_SizeFilter,
                        Ifloat &pro_ImagOut, 
                        int pi_Step, 
                        type_border ps_Border) {

  // local var
  int Nl = pro_ImagIn.nl();              // number of line
  int Nc = pro_ImagIn.nc();              // number or row
  int ai_Ind;                            // curent ind of calcul
  int ai_DemiSize = pi_SizeFilter/2;     // 
  double ad_Val;                         // 

  for (int i=0; i<Nl; i++) {             // each line
    for (int j=0; j<Nc; j++) {           // each point

      ad_Val = 0;
      for (int k=-ai_DemiSize; k <= ai_DemiSize; k++) {
	  ai_Ind = border_ind_test (j+k*pi_Step, Nc, ps_Border);
	  ad_Val += (double) pro_ImagIn(i, ai_Ind) * ppd_Filter[k+ai_DemiSize];
      }
      pro_ImagOut(i, j) = (float) ad_Val;
    }
  }
}


static void convoly1d ( Ifloat &pro_ImagIn, 
                        double* ppd_Filter, 
                        int pi_SizeFilter,
                        Ifloat &pro_ImagOut,
                        int  pi_Step, 
                        type_border ps_Border) {

  // local var
  int Nl = pro_ImagIn.nl();               // number of line
  int Nc = pro_ImagIn.nc();               // number or row
  int ai_Ind;                             // curent ind of calcul
  int ai_DemiSize = pi_SizeFilter/2;    
  double ad_Val;  

  for (int j=0; j<Nc; j++) {              // each row
    for (int i=0; i<Nl; i++) {            // each point

      ad_Val = 0;
      for (int k= -ai_DemiSize; k <= ai_DemiSize; k++) {
	  ai_Ind = border_ind_test (i+k*pi_Step, Nl, ps_Border);
	  ad_Val += (double) pro_ImagIn(ai_Ind, j) * ppd_Filter[k+ai_DemiSize];
      }
      pro_ImagOut(i, j) = (float) ad_Val;
    }
  }
}


static void convolxy1d ( Ifloat &pro_ImagIn, 
                         double* ppd_Filter, 
                         int pi_SizeFilter,
                         Ifloat &pro_ImagOut, int pi_Step, 
                         type_border ps_Border) {

  Ifloat ao_Inter (pro_ImagIn.nl(), pro_ImagIn.nc(), "ao_inter");

  convolx1d (pro_ImagIn, ppd_Filter, pi_SizeFilter, ao_Inter, pi_Step, ps_Border);
  convoly1d (ao_Inter, ppd_Filter, pi_SizeFilter, pro_ImagOut, pi_Step, ps_Border);

}


/************************************************************/
/*
void wave_2d_step_mallat_atrou (Ifloat & ImagIn,  Ifloat & ImagOut_1,  
				Ifloat & ImagOut_2, Ifloat &  ImagOut_S,
 				int Scale, Bool ModPhase, type_border ps_Border) 
{

  // loc var
  int i,s=Scale;                      // ind scale  
  int ai_Step=0;                   // step at scle s
  float Mod, Phase;
     
   // steep et scale s
   ai_Step = iround (pow((double)2., (double)s));

   // convolution over x
    convolx1d (ImagIn, std_Gn, SIZE_FILTER, 
               ImagOut_1, ai_Step, ps_Border);

    // convolution over y
    convoly1d (ImagIn, std_Gn, SIZE_FILTER, 
                ImagOut_2 , ai_Step, ps_Border);

    // convolution over xy
    convolxy1d (ImagIn, std_Hn, SIZE_FILTER, 
                ImagOut_S, ai_Step, ps_Border);

   if (ModPhase == True)
   {
      for (i=0; i< ImagIn.nl()*ImagIn.nc(); i++)  
      {              
          Mod = sqrt(ImagOut_1(i)*ImagOut_1(i) + ImagOut_2(i)*ImagOut_2(i));
 	  ARG (ImagOut_1(i), ImagOut_2(i), Phase);  
	  ImagOut_1(i) = Mod;
	  ImagOut_2(i) = Phase;
      }
   }
 }
*/
/************************************************************/
/*
void mod_phase_to_edge(Ifloat & Mod, Ifloat & Phase, Ifloat &Edge)
{
   int Nl = Edge.nl();
   int Nc = Edge.nc();
   int i,j,i1,j1,i2,j2;
   int Angle,K=1;
   float Pha;
   type_border Border = I_MIRROR;
   
// Normal direction to the edge    
//     Angle = 1 ==> 0 deg
//          = 2 ==> 45 
//          = 3 ==> 90
//          = 4 ==> 135
//          = 5 ==> 180
//          = 6 ==> 225
//          = 7 ==> 270
//          = 8 ==> 315  
    
   
   for (i=0;i< Nl;i++)
   for (j=0;j< Nc;j++)
   {
      // Pha = floatting angle between 1 and 8
      Pha =  ( Phase(i,j) >= 0) Phase(i,j) : 2.*PI - Phase(i,j);
      Pha = Pha / PI * 4. + 1;
      
      // Angle takes the closer integer value
      Angle = (int)(Pha + 0.5);
      if (Angle >= 9) Angle = 1;
      
      // the modulus value must be a maximum in the 
      // direction orthogonal to the edge
     Edge(i,j) =  Mod(i,j);
     i1 = i2 = i;
     j1 = j2 = j;
     switch(Angle)
     {
        case 1: 
	case 5: j1-=K; j2+=K; break;
	case 2:
	case 6: i1-=K; j1-=K; i2+=K; j2+=K; break;
	case 3:
	case 7: i1-=K; i2+=K; break;
	case 4:
	case 8: i1+=K; j1-=K; i2-=K; j2+=K; break;
        default:
	   cerr << "Error: bad orientation ... " << endl;
	   exit(-1);
	   break;
     }
 
     // detect local maximum
     if (    (ABS(Mod(i1,j1,Border)) >  ABS(Mod(i,j)) 
          || (ABS(Mod(i2,j2,Border)) >  ABS(Mod(i,j))) Edge(i,j) = 0;
   }
}
*/                    
/************************************************************/

void wave_2d_mallat_atrou ( Ifloat &pro_ImagIn, 
	                        Ifloat *&prpo_ImagOut,   
                            int pi_NbrPlan, 
	                        type_border ps_Border) {

  // precondition
  _mac_precondition (prpo_ImagOut!=NULL, "Out cube not allocated");
  _mac_precondition (pi_NbrPlan>=1, "Nb Plan >= 1");

  // loc var
  int s=0;                         // ind scale for loop
  int ai_Step=0;                   // step at scle s

  // int end level of prpo_ImagOut[end] with pro_ImagIn
  prpo_ImagOut[2*pi_NbrPlan-1] = pro_ImagIn;

  //Compute the wavelet coefficients
  for (s = 0; s < pi_NbrPlan-1; s++) { 
 
    // steep et scale s
    ai_Step = iround (pow((double)2., (double)s));

    // convolution over x
    convolx1d ( prpo_ImagOut[2*pi_NbrPlan-1], std_Gn, SIZE_FILTER, 
                prpo_ImagOut[2*s], ai_Step, ps_Border);

    // convolution over y
    convoly1d ( prpo_ImagOut[2*pi_NbrPlan-1], std_Gn, SIZE_FILTER, 
                prpo_ImagOut[2*s+1], ai_Step, ps_Border);

    // convolution over xy
    convolxy1d ( prpo_ImagOut[2*pi_NbrPlan-1], std_Hn, SIZE_FILTER, 
                 prpo_ImagOut[2*pi_NbrPlan-2], ai_Step, ps_Border);

    prpo_ImagOut[2*pi_NbrPlan-1] = prpo_ImagOut[2*pi_NbrPlan-2];
    
  }
}


void rec_wave_2d_mallat_atrou ( Ifloat *&prpo_ImagIn, 
	                            Ifloat &pro_ImagOut,   
                                int pi_NbrPlan, 
	                            type_border ps_Border) {

  // precondition
  _mac_precondition (prpo_ImagIn!=NULL, "IN cube not initialised");
  _mac_precondition (pi_NbrPlan>=1, "Nb Plan >= 1");

  // loc var
  int s=0;                         // ind scale for loop
  int ai_Step=0;                   // step at scle s

  // inter rec Imag
  Ifloat ao_S1(pro_ImagOut.nl(), pro_ImagOut.nc(), "ao_S1");
  Ifloat ao_S2(pro_ImagOut.nl(), pro_ImagOut.nc(), "ao_S2");
  Ifloat ao_S3(pro_ImagOut.nl(), pro_ImagOut.nc(), "ao_S3");
  Ifloat ao_Inter(pro_ImagOut.nl(), pro_ImagOut.nc(), "ao_Inter");

  // init pro_ImagOut with end level of prpo_ImagIn
  pro_ImagOut = prpo_ImagIn[2*pi_NbrPlan-2];

  //Compute the wavelet coefficients
  for (s = pi_NbrPlan-2; s >= 0 ; s--) { 
 
    // steep et scale s
    ai_Step = iround (pow((double)2., (double)s));

    // convolution x and y on Dx for S1(s)
    convolx1d ( prpo_ImagIn[2*s], std_Kn, SIZE_FILTER, 
                ao_Inter, ai_Step, ps_Border);

    convoly1d ( ao_Inter, std_Ln, SIZE_FILTER, 
                ao_S1, ai_Step, ps_Border);

    // convolution x and y on Dy for S2(s)
    convolx1d ( prpo_ImagIn[2*s+1], std_Ln, SIZE_FILTER, 
                ao_Inter, ai_Step, ps_Border);

    convoly1d ( ao_Inter, std_Kn, SIZE_FILTER, 
                ao_S2, ai_Step, ps_Border);

    // convolution x and y on Dy for S3(s)
    convolxy1d ( pro_ImagOut, std_HnMod, SIZE_FILTER, 
                 ao_S3, ai_Step, ps_Border);

    pro_ImagOut = ao_S1;
    pro_ImagOut += ao_S2;
    pro_ImagOut += ao_S3;
  }
}



void calc_mod_pha_wave_2d_mallat_atrou (Ifloat *&prpo_W2dTransfIn,
                                        Ifloat *&prpo_Modulus,
                                        Ifloat *&prpo_Phase,
                                        int pi_NbrPlan) {

  // warning : all scales must have same dimension in x and y
  //           dir x is at scale 2*i and Y at 2*i+1

  // local var
  int Nl = prpo_W2dTransfIn[0].nl();         // number of line
  int Nc = prpo_W2dTransfIn[0].nc();         // number of row

  //Compute the wavelet coefficients
  for (int s = pi_NbrPlan-2; s >= 0 ; s--) { // each scale
    for (int i=0; i<Nl; i++) {               // each line
      for (int j=0; j<Nc; j++) {             // each point
	  prpo_Modulus[s](i,j) = 
            sqrt( prpo_W2dTransfIn[2*s](i,j)*prpo_W2dTransfIn[2*s](i,j) + 
                  prpo_W2dTransfIn[2*s+1](i,j)*prpo_W2dTransfIn[2*s+1](i,j));
	  ARG ( prpo_W2dTransfIn[2*s](i,j), prpo_W2dTransfIn[2*s+1](i,j),
            prpo_Phase[s](i,j));
      }
    }
  }
}




/************************************************************************
******* Maxima Wavelet Transform ****************************************
************************************************************************/



// interpolate between two consecutive maxima modulus at row pi_Row
void Ortho_Proj_Operator (int pi_Scale,            // current scale 
			  int pi_Line,             // current line 
                          int pi_Row ,             // current rowe
                          int pi_Begining,         // ind of begining
	                  int pi_End,              // ind of end
                          Ifloat *&prpo_ImagOut,   // Data IN
			  MultiResol& pro_Mr2dRecData) {  // Rec Data OUT


  int pi_NumScale = pi_Scale/2 +1;   // pi_NumScale goes from 1 to Max, 
                                     // and pi_Scale from 0 to max-1
  pi_End = pi_End - pi_Begining;
  int ai_Dec = pi_Begining;
  pi_Begining = 0;

  // coef of interpolate function
  double af_Step1 = pow((double)2., (double)(-pi_NumScale))*pi_Begining;
  double af_Step2 = pow((double)2., (double)(-pi_NumScale))*pi_End;
  double ad_Delta = exp(af_Step1) * exp(-af_Step2) 
                  - exp(af_Step2) * exp(-af_Step1);

  int ai_FirstIndLine  = (pi_Line != -1 ? pi_Line : pi_Begining+ai_Dec);
  int ai_SecondIndLine = (pi_Line != -1 ? pi_Line : pi_End+ai_Dec);
  int ai_FirstIndRow   = (pi_Row  != -1 ? pi_Row : pi_Begining+ai_Dec);
  int ai_SecondIndRow  = (pi_Row  != -1 ? pi_Row : pi_End+ai_Dec);

  double ad_DeltAlpha = 
          (prpo_ImagOut[pi_Scale](ai_FirstIndLine, ai_FirstIndRow) - 
           pro_Mr2dRecData (pi_Scale, ai_FirstIndLine,ai_FirstIndRow))
           * exp(-af_Step2) -  exp(-af_Step1) *
          (prpo_ImagOut[pi_Scale](ai_SecondIndLine, ai_SecondIndRow) - 
           pro_Mr2dRecData(pi_Scale, ai_SecondIndLine, ai_SecondIndRow));


  double ad_DeltBeta  = 
          (prpo_ImagOut[pi_Scale](ai_SecondIndLine, ai_SecondIndRow) - 
           pro_Mr2dRecData (pi_Scale, ai_SecondIndLine, ai_SecondIndRow))
           * exp(af_Step1) - exp(af_Step2) *
          (prpo_ImagOut[pi_Scale](ai_FirstIndLine, ai_FirstIndRow) - 
           pro_Mr2dRecData (pi_Scale, ai_FirstIndLine, ai_FirstIndRow));

  double ad_Alpha = ad_DeltAlpha / ad_Delta;
  double ad_Beta  = ad_DeltBeta  / ad_Delta;

  // interpolate between pi_Begining and pi_End 
  int ai_IndLine=0, ai_IndRow=0;
  for (int i=pi_Begining+1; i<pi_End; i++) {
      double ad_inter = exp (pow((double)2., (double)(-pi_NumScale))*i)
                    * ad_Alpha + ad_Beta *
                    exp (- pow((double)2., (double)(-pi_NumScale*i))*i);
      ai_IndLine = (pi_Line != -1 ? pi_Line : i+ai_Dec);
      ai_IndRow  = (pi_Row  != -1 ? pi_Row : i+ai_Dec);
      pro_Mr2dRecData (pi_Scale, ai_IndLine, ai_IndRow) += ad_inter;
  }
}






void interpolate  (int pi_NbrPlan,               // Number of scale
		   intarray& pro_NumberMaxInLine,// Number of max in line
                   intarray** pppo_MaxModInLine, // loc of maxmod on line
                   intarray& pro_NumberMaxInRow, // Number of max in line
                   intarray** pppo_MaxModInRow,  // loc of maxmod on line
                   Ifloat *&prpo_ImagOut,       // Data IN
                   MultiResol& pro_Mr2dRecData) {// Rec Data OUT



  int ai_Nl = pro_Mr2dRecData.size_ima_nl();   // number of line of image
  int ai_Nr = pro_Mr2dRecData.size_ima_nc();   // number of row of image

  // x plane
  for (int s=0; s<pi_NbrPlan-1; s+=2) {

    // get a line y=cste
    for (int i=0; i<ai_Nl; i++) {

      // interpolate between max of line i at scale s
      for (int l=0; l<pro_NumberMaxInLine(s,i)-1; l++){

	Ortho_Proj_Operator (2*s, i, -1, 
                             (**(pppo_MaxModInLine+s*ai_Nl+i))(l), 
                             (**(pppo_MaxModInLine+s*ai_Nl+i))(l+1),
                             prpo_ImagOut, pro_Mr2dRecData);
      }
    }

    // get a row
    for (int j=0; j<ai_Nr; j++) {

      // interpolate between max of line i at scale s
      for (int l=0; l<pro_NumberMaxInRow(s,j)-1; l++){

	Ortho_Proj_Operator (2*s+1, -1, j, 
                             (**(pppo_MaxModInRow+s*ai_Nr+j))(l), 
			     (**(pppo_MaxModInRow+s*ai_Nr+j))(l+1),
			     prpo_ImagOut, pro_Mr2dRecData);
      }
    }
  }
}



void init_last_scale (int pi_NbrPlan,                // number of plan
                      int pi_Nl,                     // number of line
                      int pi_Nr,                     // number of row
                      Ifloat *&prpo_ImagOut,         // Rec Data IN
                      MultiResol& po_Mr2dRecData) {  // Rec Data OUT
   for (int i=0; i<pi_Nl; i++)
     for (int j=0; j<pi_Nr; j++) 
       po_Mr2dRecData (2*pi_NbrPlan-1,i,j) = 
	 prpo_ImagOut[2*(pi_NbrPlan-1)](i,j);
}

void init_max (int pi_NbrPlan,                // number of plan
               int pi_Nl,                     // number of line
               int pi_Nr,                     // number of row
	       intarray& pro_NumberMaxInLine, // Number of max in line
               intarray** pppo_MaxModInLine,  // loc of max mod on line
               intarray& pro_NumberMaxInRow,  // Number of max in line
               intarray** pppo_MaxModInRow,   // loc of max mod on lin
               Ifloat *&prpo_ImagOut,         // Data IN
               MultiResol& pro_Mr2dRecData) { // max modulus of Data IN
  for (int s=0; s<pi_NbrPlan-1; s++) {

    for (int i=0; i<pi_Nl; i++) 
      for (int l=0; l<pro_NumberMaxInLine(s,i); l++)
	if (pro_NumberMaxInLine(s,i) > 0)
	  pro_Mr2dRecData (2*s, i, (**(pppo_MaxModInLine+s*pi_Nl+i))(l))
	    = prpo_ImagOut[2*s](i, (**(pppo_MaxModInLine+s*pi_Nl+i))(l));
  
	
    for (int j=0; j<pi_Nr; j++)
      for (int l=0; l<pro_NumberMaxInRow(s,j); l++)
	if (pro_NumberMaxInRow(s,j) > 0)
	  pro_Mr2dRecData (2*s+1, (**(pppo_MaxModInRow+s*pi_Nr+j))(l),j)
	    = prpo_ImagOut[2*s+1]((**(pppo_MaxModInRow+s*pi_Nr+j))(l),j);
  }
}

void search_number_max (int pi_NbrPlan,
			Ifloat*& prpo_ModData, 
                        intarray& pro_NumberMaxInLine,
                        intarray& pro_NumberMaxInRow) {
  int s;
  for (s=0; s<pi_NbrPlan-1; s++) 
    for (int i=0; i<prpo_ModData[s].nl(); i++)
      for (int j=0; j<prpo_ModData[s].nc(); j++) 
	if (prpo_ModData[s](i,j) < -EPS || prpo_ModData[s](i,j) > EPS) {
	  pro_NumberMaxInLine(s,i)++;
	  pro_NumberMaxInRow(s,j)++;
	}

  int ai_MaxLine=0, ai_MaxRow=0;
  for (s=0; s<pi_NbrPlan-1; s++) {
    for (int i=0; i<prpo_ModData[s].nl(); i++)
      ai_MaxLine += pro_NumberMaxInLine(s,i);
    for (int j=0; j<prpo_ModData[s].nc(); j++)
      ai_MaxRow += pro_NumberMaxInRow(s,j);
    cout << "Number max at scale " << s << " on line : " << ai_MaxLine <<
      endl;
    cout << "Number max at scale " << s << " on row  : " << ai_MaxRow << 
      endl;
  }
}

/*

void wave_max_2d_mallat_atrou ( Ifloat &pro_ImagIn, 
	                        Ifloat *&prpo_ImagOut,   
                                int pi_NbrPlan, 
	                        type_border ps_Border) {

  //---------------------------------------------------------------------
  //--------- IMAGE IN --------------------------------------------------
  //---------------------------------------------------------------------
  // length
  int ai_Nl = pro_ImagIn.nl(); 
  int ai_Nr = pro_ImagIn.nc();
  int s,b;

  // construct multiresolution data (DIADIC_MALLAT)
  wave_2d_mallat_atrou (pro_ImagIn, prpo_ImagOut, pi_NbrPlan, ps_Border);

  // construct Ifloat[] modulus and phase struct 
  // for calc_modulus_and_phase
  Ifloat* apo_ModData = new Ifloat [pi_NbrPlan-1];
  Ifloat* apo_PhaData = new Ifloat [pi_NbrPlan-1];
  for (s=0; s<pi_NbrPlan-1; s++) {
    apo_ModData[s] = pro_ImagIn; apo_PhaData[s] = pro_ImagIn;
  }

  // compute moduls and phase of image IN at all scale
  calc_mod_pha_wave_2d_mallat_atrou (prpo_ImagOut, apo_ModData, 
				     apo_PhaData, pi_NbrPlan);

  // compute local maximun wavelet of image IN 
  // we used class EDGE, method edge_from_1deriv
  for (s=0; s<pi_NbrPlan-1;s++) {
    EDGE ao_EdgeDetect (ED_CANNY);
    ao_EdgeDetect.edge_from_1deriv (apo_ModData[s], apo_PhaData[s]);
    //sprintf (atc_FileName, "maxmod_%d", s);
    //io_write_ima_float (atc_FileName, apo_ModData[s]);
  }


  // search coord (s,i,j) of all max max modulus 
  // vector ao_NumberMaxInLine    : number of max on line i at scale s
  // vector ao_NumberMaxInRow     : number of max on row j at scale s
  // intarray attpo_MaxModInLine  : coord of all max on line i at scale s
  // intarray attpo_MaxModInRow   : coord of all max on row j at scale s
  intarray ao_NumberMaxInLine (pi_NbrPlan, ai_Nl);
  intarray ao_NumberMaxInRow (pi_NbrPlan, ai_Nr);
  search_number_max (pi_NbrPlan, apo_ModData, ao_NumberMaxInLine, 
                     ao_NumberMaxInRow);
  intarray* attpo_MaxModInLine [pi_NbrPlan][ai_Nl];
  intarray* attpo_MaxModInRow [pi_NbrPlan][ ai_Nr];
  for (s=0; s<pi_NbrPlan-1; s++) {
    for (int i=0; i<ai_Nl; i++) {
      if (ao_NumberMaxInLine(s,i) > 0) {
	attpo_MaxModInLine [s][i] = new intarray (ao_NumberMaxInLine(s,i));
	int l=0;
	for (int j=0; j<ai_Nr; j++) 
	  if (apo_ModData[s](i,j) < -EPS || apo_ModData[s](i,j) > EPS)
	    (*attpo_MaxModInLine[s][i])(l++) = j;
      }
    }
    for (int j=0; j<ai_Nr; j++) {
      if (ao_NumberMaxInRow(s,j) > 0) {
	attpo_MaxModInRow [s][j] = new intarray (ao_NumberMaxInRow(s,j));
	int l=0;
	for (int i=0; i<ai_Nl; i++) 
	  if (apo_ModData[s](i,j) < -EPS || apo_ModData[s](i,j) > EPS)
	    (*attpo_MaxModInRow[s][j])(l++) = i;
      }
    }
  } 


  //---------------------------------------------------------------------
  //--------- RECONSTRUCT IMAGE -----------------------------------------
  //---------------------------------------------------------------------
  // init recons signal struct (format Ifloat)
  // used in rec of ao_Mr2dRecData
  Ifloat ao_RecData(ai_Nl, ai_Nr, "Recons");

  // create multiresol struct of reconstruct image (init)
  MultiResol ao_Mr2dRecData (ai_Nl, ai_Nr, pi_NbrPlan, TO_DIADIC_MALLAT, 
                             "Mr1d Reconstuct Signal");

  // construct multiresol (Ifloat (1 image) -> Multiresol (2s+1 images))
  ao_Mr2dRecData.transform (pro_ImagIn);

  // raz of all coef of recons signal (except last scale)
  for (b=0; b<2*(pi_NbrPlan-1); b++)
    for (int i=0; i<ai_Nl; i++) 
      for (int j=0; j<ai_Nr; j++)
  	ao_Mr2dRecData (b,i,j) = 0.0;

  // int last scale
  init_last_scale (pi_NbrPlan, ai_Nl, ai_Nr, prpo_ImagOut, 
		   ao_Mr2dRecData);


  // --------------------------------------------------------------------
  // --------- ITERATION FOR RECONSTRUCTION -----------------------------
  // --------------------------------------------------------------------
  // iteration for recons
  for (int iter=0; iter<5; iter++){

    //cout << "number iteration :" << iter << endl;;
	//if (iter==0) cout << " -- work in progress ...";
	//else cout << ".";

    // --------------------------------------------------- 
    // ---------------------------------------------------
    interpolate (pi_NbrPlan,
                 ao_NumberMaxInLine, (intarray**)attpo_MaxModInLine, 
                 ao_NumberMaxInRow,  (intarray**)attpo_MaxModInRow, 
                 prpo_ImagOut, ao_Mr2dRecData);

    // modify the modulus maximum of Mr1dRecDat : 
    // max (Mr1dRecData(s,i)) =  max (Mr1dData(s,i))
    init_max (pi_NbrPlan, ai_Nl, ai_Nr, 
	      ao_NumberMaxInLine, (intarray**)attpo_MaxModInLine,
	      ao_NumberMaxInRow,  (intarray**)attpo_MaxModInRow, 
	      prpo_ImagOut, ao_Mr2dRecData);

    // int last scale
    init_last_scale (pi_NbrPlan, ai_Nl, ai_Nr, prpo_ImagOut, 
		     ao_Mr2dRecData);

    // --------------------------------------------------- 
    // ---------------------------------------------------
    // invers and foward transform of Mr1dRecDat
    ao_Mr2dRecData.recons (ao_RecData);
    //sprintf (atc_FileName, "fic_%d", iter);
    //io_write_ima_float (atc_FileName, ao_RecData);
    //cout << "Relative quadartic error (%) :" 
    //   << calc_som_quad_error (ao_Dat, ao_RecData) <<endl;
    //cout << "Relative flux error (%) :" 
    //   << calc_som_error (ao_Dat, ao_RecData) <<endl;
    ao_Mr2dRecData.transform (ao_RecData);
  }
  //cout << ", work done" << endl;

  // --------------------------------------------------------------------
  // --------- COPY IMAGE REC IN IMAGE IN   -----------------------------
  // --------------------------------------------------------------------
  
  for (b=0; b<2*pi_NbrPlan-1; b++)
    for (int i=0; i<ai_Nl; i++) 
      for (int j=0; j<ai_Nr; j++)
  	prpo_ImagOut[b](i,j) = ao_Mr2dRecData (b,i,j);


}




void rec_wave_max_2d_mallat_atrou ( Ifloat *&prpo_ImagIn, 
	                            Ifloat &pro_ImagOut,   
                                    int pi_NbrPlan, 
	                            type_border ps_Border) {

  rec_wave_2d_mallat_atrou (prpo_ImagIn, pro_ImagOut, pi_NbrPlan,
                            ps_Border);
}
*/
