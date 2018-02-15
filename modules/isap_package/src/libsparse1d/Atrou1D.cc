 

#include "SB_Filter1D.h"
 
#define FIL_DEBUG 0


/*********************************************************************/
// ATROU,      
/*********************************************************************/
/****************************************************************************/

void ATROUS_1D_WT::b3spline_filtering (fltarray& Sig_in, 
                                       fltarray& Sig_out,  
                                       int Step_trou) 
{
    int Nx = Sig_in.nx();
    int i,Step;
    double Coeff_h0 = 3. / 8.;
    double Coeff_h1 = 1. / 4.;
    double Coeff_h2 = 1. / 16.;     
    Step = (int)(pow((double)2., (double) Step_trou) + 0.5);
    for (i=0; i<Nx; i ++) 
        Sig_out(i) = Coeff_h0 * Sig_in(i)
                  + Coeff_h1 * (Sig_in (i-Step, Bord)  + Sig_in (i+Step, Bord)) 
                  + Coeff_h2 * (Sig_in (i-2*Step, Bord)  + Sig_in (i+2*Step, Bord));
}


/****************************************************************************/

void ATROUS_1D_WT::transform (fltarray& Signal, fltarray* TabBand) 
{
   int s;
   TabBand[0] = Signal;
   fltarray Data_out;
   if (ModifiedAWT == True) Data_out.alloc(NbPts);
       
   for (s=0; s< NbScale-1; s++) 
   {
       b3spline_filtering(TabBand[s], TabBand[s+1], s);
       if (ModifiedAWT == True) 
       {
          b3spline_filtering(TabBand[s+1], Data_out, s);
          TabBand[s] -= Data_out;
       }
       else TabBand[s] -= TabBand[s+1];
   }
}

/****************************************************************************/

void ATROUS_1D_WT::recons(fltarray* TabBand, fltarray& Signal) 
{
   int s,i;
   
   if (Signal.nx() != NbPts) Signal.alloc(NbPts);
   else Signal.init();
   fltarray Data_out;
    
   if ((ModifiedAWT == False) && (AdjointRec == False))
           for (s=0; s<NbScale; s++) Signal += TabBand[s];
   else 
   {
      Data_out.alloc(NbPts);
      Signal = TabBand[NbScale-1];
      for (s=NbScale-2; s>= 0 ; s--) 
      {
 	  b3spline_filtering (Signal, Data_out, s);
	  for (i=0; i < NbPts; i++) Signal(i) = Data_out(i) + (TabBand[s])(i);
      }
   }
}

/****************************************************************************/

void ATROUS_1D_WT::free(fltarray* TabBand, int Nbr_Plan) 
{
    if (Nbr_Plan != 0) delete [] TabBand;
}

/****************************************************************************/

void ATROUS_1D_WT::alloc (fltarray* & TabBand, int Nx, int NbrBand) 
{
    NbScale = NbrBand;    
    FirstDetectScale = 0;
    NbPts = Nx;
    char ch[80];
    TabBand = new fltarray [NbrBand];
    for (int s=0; s<NbrBand; s++) 
    {
 	sprintf (ch, "band_%d", s+1);
        TabBand[s].alloc (Nx, ch);   
    }
}

/****************************************************************************/

float ATROUS_1D_WT::norm_band(int s) 
{
   static double TN[10] = 
    {0.721618,  0.286209, 0.178281, 0.121010, 0.083104, 0.0588543, 0.0446171,
     0.0289142, 0.0177307, 0.};
    static double TN_Modified[10] = 
    {0.721618,  0.286209, 0.178281, 0.121010, 0.083104, 0.0588543, 0.0446171,
     0.0289142, 0.0177307, 0.};
   if ((s < 0) || (s >= 10)) return 0;
   else if (ModifiedAWT == False) return TN[s];
   else return TN_Modified[s];
}

/****************************************************************************/
      
void ATROUS_1D_WT::KillScaleNotUsed (fltarray* AT_Trans, int FirstDetScale) 
{
   FirstDetectScale = FirstDetScale;
   if (FirstDetectScale != 0) 
   {
      for (int j=0; j < FirstDetectScale; j++) AT_Trans[j].init();
   }
}

/****************************************************************************/

void ATROUS_1D_WT::KillLastScale (fltarray* AT_Trans) 
{
   AT_Trans[NbScale-1].init();
}
/****************************************************************************/

void ATROUS_1D_WT::Threshold (fltarray* AT_Trans, float NSigma, int IterNumber) 
{
  int NumberDetectedCoef=0;
  for (int s=FirstDetectScale; s<NbScale-1; s++) 
  {
      float NSig = (s == 0) ? NSigma + 1: NSigma;
      if (NSigma == 0) NSig = 0;
      float Norm = (s == NbScale-1) ? norm_band(s-1): norm_band(s);
      float Noise =  SigmaNoise*Norm;
      float Level = NSig*Noise;
  
      // if (Verbose) 
      //cout << "Level atrou:" << Level/Norm << ", Noise:" 
      //     << Noise/Norm << endl;
      
      for (int i=0; i< NbPts; i++) 
      {
         float Coef = (AT_Trans[s])(i);  
	 (AT_Trans[s])(i) = update(Coef, Level, Noise);
         if ((SuppressIsolatedPixel) && (i > 0) && (i < NbPts-1))
	 {
	     if (   !is_on_support ((AT_Trans[s])(i-1), Level, Noise)
	         && !is_on_support ((AT_Trans[s])(i+1), Level, Noise))
	        (AT_Trans[s])(i) = 0.;
	 }
	 if (OnlyPositivDetect) 
	 {
	      if ((AT_Trans[s])(i) < 0) (AT_Trans[s])(i) = 0;
	 }
         if (fabs((AT_Trans[s])(i)) > 0) NumberDetectedCoef++;
      }
      if (Write)
      {
         char Name[256];
         sprintf (Name,"atrouWT_Thres_sc%d_iter%d",s+1,IterNumber);
         // fits_write_fltarr (Name, AT_Trans[s]);
      }
   }
   if (Verbose) 
      cout << "Number detected coef : " << NumberDetectedCoef << endl;
           
}


/****************************************************************************/

inline float ATROUS_1D_WT::update (float CoefSol, float Threshold,  float SoftLevel) 
{
  float NewCoef = hard_threshold(CoefSol, Threshold);
  if (UseNormL1 == True) 
  {
       float SoftL =  SoftLevel*NSigmaSoft;
       float T2 = Threshold/2.;
       if (SoftL < T2)
            NewCoef = soft_threshold(NewCoef, SoftL);
       else NewCoef = soft_threshold(NewCoef, T2);
   }
   return NewCoef;  
}

/****************************************************************************/

inline bool ATROUS_1D_WT::is_on_support (float CoefSol, float Threshold, 
                                         float SoftLevel) 
{					 
  float NewCoef = hard_threshold(CoefSol, Threshold);
  if (UseNormL1 == True) {
       float SoftL =  SoftLevel*NSigmaSoft;
       float T2 = Threshold/2.;
       if (SoftL < T2)
            NewCoef = soft_threshold(NewCoef, SoftL);
       else NewCoef = soft_threshold(NewCoef, T2);
   }
   return (!(NewCoef == 0));  
} 

/****************************************************************************/

float ATROUS_1D_WT::getAbsMaxTransf (fltarray* AT_Trans) 
{
   float absMaxLevel=0;
   for (int s=0; s<NbScale-1; s++) 
   {
      if (OnlyPositivDetect) (AT_Trans[s]).inf_threshold(0.0);
      
      float MaxLevel = fabs((AT_Trans[s]).maxfabs());
      float NormMaxLevel = MaxLevel/norm_band(s)/SigmaNoise;
                   
      if (Write)
         cout << "Lambda atrou (" << s+1 << ") : " << NormMaxLevel 
         << ", Max amplitude : " << MaxLevel << endl;  
                                                              
      if (absMaxLevel <= NormMaxLevel) absMaxLevel = NormMaxLevel;
   }
   if (Verbose) cout << "Atrou proj - Max abs level : " << absMaxLevel << endl;

   return absMaxLevel;
}

