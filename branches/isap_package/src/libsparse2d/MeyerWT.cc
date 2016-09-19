/*******************************************************************************
**
**    UNIT
**
**    Version: 1.0
**
**    Author: Jean-Luc Starck
**
**    Date:  21/01/2005
**    
**    File:  MeyerWT.cc
**
**    Modification history:
**       A. Woiselle 2010 : - Changed the norms to theoretical ones
**                          - Corrected the get_hfilter function (wrong floor(x))
**
*******************************************************************************
**
**    DESCRIPTION  Meyer wavelet decomposition
**    -----------  
**
******************************************************************************/
 

#include "MeyerWT.h"

/*********************************************************************/

// theoretical norms for non extended Meyer wavelets :
//  [sqrt(1-(3/8)^2), sqrt((3/4)^2-(3/8)^2) , 3/4 ]
static float Tab_Meyer2D[3] = 
	{ 0.9270, 0.6495, 0.75 };
 
/*********************************************************************/
 
void MEYER_WT::get_hfilter(Ifloat &H, double DNl, double DNc)
{
	int i,j;
	int Nl = H.nl();
	int Nc = H.nc();
	double Nl2 = (DNl > 0) ? DNl /2.: Nl / 2.;
	double Nc2 = (DNc > 0) ? DNc /2.: Nc / 2.;
	double Nl4 = Nl2 / 2.;
	double Nc4 = Nc2 / 2.;

	if (Isotrop == False)
	{
		dblarray H1(Nl);
		// fltarray HH1(Nl);
		dblarray H2(Nc);

		H1.init(1.);
		H2.init(1.);
		if (Extend == True)
		{
			for(i=0; i< Nl4; i++)
			{
				double x = (i -  int(Nl/2) + Nl2) / Nl4;
				H1(i) =  H1(Nl-i-1) =  lowpass_window(x);
			}
			for(j=0; j< Nc4; j++)
			{
				double x = (j - int(Nc/2) + Nc2) / Nc4;
				H2(j) =  H2(Nc-j-1) =  lowpass_window(x);
			}
		}
		else
		{
			for(i=0; i < Nl; i++)
			{
				double x = (i - int(Nl/2));
				double r = (ABS(x) - Nl4) / Nl4;
				if (r <= 0) H1(i) = 1.;
				else if (r <= 1) H1(i) = lowpass_window(1. - r);
				else H1(i) = 0.;
				// HH1(i) = H1(i);
			}
			// fits_write_fltarr("xx",  HH1);
			for(j=0; j< Nc; j++)
			{
				double x = (j - int(Nc/2));
				double r = (ABS(x) - Nc4) / Nc4;
				if (r <= 0) H2(j) = 1.;
				else if (r <= 1) H2(j) = lowpass_window(1. - r);
				else H2(j) = 0.;
				// HH1(j) = H2(j);
			}
			// fits_write_fltarr("xx1",  HH1);
		}
		for(i=0; i< Nl; i++)
		for(j=0; j< Nc; j++) H(i,j) = H1(i)*H2(j);
	}
	else
	{
		// cout << "ISOTROP" << endl;
		int N = MIN(Nl/4,Nc/4);
		H.init(1.);
		for(i=0; i< Nl; i++)
		for(j=0; j< Nc; j++)
		{
			double x = (i - Nl/2);
			double y = (j - Nc/2);
			double r = (sqrt(x*x+y*y) - N) / N;
			if (r <= 0) H(i,j) = 1.;
			else if (r <= 1) H(i,j) = lowpass_window(1. - r);
			else H(i,j) = 0.;
		}
	}
	// io_write_ima_float("low_pass.fits", H);
	// exit(0);
}

/*********************************************************************/

float MEYER_WT::get_norm(int s)
{
	// coef that takes into account the Fourier zero padding (rounded sizes and sometimes imposed odd size)
	// coef = sqrt( theoretical_redundancy / real_redundancy )
	double coef = sqrt( (float(D_ExtNl)/pow((double) 2,(double) s)*float(D_ExtNc)/pow((double) 2,(double) s))/(TabNl(s)*TabNc(s)) );
	if(s==0 && Extend==False) return Tab_Meyer2D[0]*coef;
	else if(s==NbrScale-1) return Tab_Meyer2D[2]*coef;
	else return Tab_Meyer2D[1]*coef;
}

void MEYER_WT::init(int Nbr_Scale, int Nl, int Nc, Bool ExtendWT, Bool IsotropWT, Bool WTNeedOddSize)
{
   NbrScale = Nbr_Scale;
   Nl_Imag = Nl;
   Nc_Imag = Nc;
   TabNl.alloc(NbrScale);
   TabNc.alloc(NbrScale);
   Extend = ExtendWT;
   Isotrop = IsotropWT;
   NeedOddSize = WTNeedOddSize;
   
   if (Isotrop  == True) Extend = False;
   if (Verbose == True)
   {
      cout << "INIT WT: " << "ImageSize = " << Nl << " " << Nc << " NbrScale = " << NbrScale << endl;
      if (Isotrop == True) cout << "   Use isotropic wavelets " << endl;
      else if (Extend == False) cout << "   Use Meyer's wavelets without image extension  " << endl;
      else cout << "   Use Meyer's wavelets " << endl;
   }
   
   if (Extend==True)
   {
       D_ExtNl =(4. / 3. * (float) Nl);
       D_ExtNc =(4. / 3. * (float) Nc);
   }
   else
   {
       D_ExtNl = (float) Nl;
       D_ExtNc = (float) Nc;
   }
   double DNL=D_ExtNl;
   double DNC=D_ExtNc;
   for (int s=0; s < Nbr_Scale; s++)
   {
      if ((Extend==True) || (NeedOddSize == True))
      {
         TabNl(s) = 2 * int(floor(DNL/2)) + 1; 
         TabNc(s) = 2 * int(floor(DNC/2)) + 1;  
      }
      else 
      {
          TabNl(s) = int( DNL + 0.5);
	  TabNc(s) = int( DNC + 0.5);
      }
      DNL /= 2.;
      DNC /= 2.;
   }
   ExtNl = TabNl(0);
   ExtNc = TabNc(0);
   
   // Memory allocation
   Tabcf_WT_Band = new Icomplex_f[NbrScale];
   for  (int s=0; s < NbrScale; s++) Tabcf_WT_Band[s].alloc(TabNl(s), TabNc(s));
   
   TF_ExtData.alloc(ExtNl,ExtNc,"Ext TFIMA");
   H.alloc(ExtNl,ExtNc,"H filter");
}

/*********************************************************************/

void MEYER_WT::get_extFourier(Ifloat &Data, Icomplex_f & TF_ExtData)
{
    int i,j;
    int Nl = Data.nl();
    int Nc = Data.nc(); 
    int Nl2 = Nl/2; 
    int Nc2 = Nc/2;  
    
    // cout << "get_extFourier" << endl;
     TF_ExtData.resize(ExtNl,ExtNc);
     if (Extend == True)
     {
        Icomplex_f TF_Ima(Nl,Nc,"TFIMA");
        int ExtNl2 = ExtNl/2;
        int ExtNc2 = ExtNc/2;
        FFT2D.fftn2d(Data, TF_Ima, False);
        float  Norm = sqrt(float(Nl*Nc));
        for(i=0; i<Nl; i++)
        for(j=0; j<Nc; j++)  TF_Ima(i,j) /= Norm;  
          
        for(i=0; i<ExtNl; i++)
        for(j=0; j<ExtNc; j++) 
        {
          int u = i - ExtNl2;
          int v = j - ExtNc2;
          if (u < -Nl2) u += Nl;
          else if (u >= Nl2) u -= Nl;
          if (v < -Nc2) v += Nc;
          else if (v >= Nc2) v -= Nc;
          TF_ExtData(i,j) = TF_Ima(u+Nl2,v+Nc2);
       }
    }
    else 
    {
       FFT2D.fftn2d(Data, TF_ExtData, False);
       float  Norm = sqrt(float(Nl*Nc));
       for(i=0; i<Nl; i++)
       for(j=0; j<Nc; j++)  TF_ExtData(i,j) /= Norm;
    }    
     // cout << "END get_extFourier" << endl;

}

/*********************************************************************/

void MEYER_WT::get_IFFT_Ima(Icomplex_f & TF_ExtData, Ifloat &Data)
{
    int i,j;
    int Nl = Data.nl();
    int Nc = Data.nc(); 

      // cout << "get_extFourier" << endl;
      
    if ( (ExtNl != TF_ExtData.nl()) || (ExtNc != TF_ExtData.nc()))
    {
       cout << "Error: bad image size in get_IFFT_Ima ... " << endl;
       cout << "ExtNl = " << ExtNl << "ExtNc = " << ExtNc << endl;
        cout << "InExtNl = " << TF_ExtData.nl() << "InExtNc = " << TF_ExtData.nc() << endl;
       exit(-1);
    }
    
    // cout << "alloc" << ExtNl << " " << ExtNc << endl;
    if (Extend == False)
    {
       float  Norm = sqrt(float(Nl*Nc));
       for(i=0; i<Nl; i++)
       for(j=0; j<Nc; j++)  TF_ExtData(i,j) *= Norm;
       FFT2D.ifftn2d(TF_ExtData);
       for(i=0; i<Nl; i++)
       for(j=0; j<Nc; j++)  Data(i,j) = TF_ExtData(i,j).real();    
    }
    else
    {
      Icomplex_f TF_Ima(Nl,Nc,"TFIMA"); 
      int ExtNl2 = ExtNl/2;
      int ExtNc2 = ExtNc/2;
      int Nl2 = Nl/2; 
      int Nc2 = Nc/2;  
      for(i=0; i<ExtNl; i++)
      for(j=0; j<ExtNc; j++) 
      {
       int u = i - ExtNl2;
       int v = j - ExtNc2;
       if (u < -Nl2) u += Nl;
       else if (u >= Nl2) u -= Nl;
       if (v < -Nc2) v += Nc;
       else if (v >= Nc2) v -= Nc;
       TF_Ima(u+Nl2,v+Nc2) += TF_ExtData(i,j);
      }
      // cout << "END get_extFourier" << endl;
      float  Norm = sqrt(float(Nl*Nc));
      for(i=0; i<Nl; i++)
      for(j=0; j<Nc; j++)  TF_Ima(i,j) *= Norm;
    
      FFT2D.ifftn2d(TF_Ima);
      for(i=0; i<Nl; i++)
      for(j=0; j<Nc; j++)  Data(i,j) = TF_Ima(i,j).real();
    }
 }

/*********************************************************************/

void MEYER_WT::transform_cf(Icomplex_f & TF_ExtData,  Icomplex_f * &TabWT)
{
   int i,j,s,Nls,Ncs;
   int Nl = TF_ExtData.nl();
   int Nc = TF_ExtData.nc();
   double DNL = D_ExtNl;
   double DNC = D_ExtNc;
   
   // TabWT = new Icomplex_f [NbrScale];
   if (TabWT == NULL)
   {
       cout << "Error in transform_cf: TabWT has not been initialized ... " << endl;
       exit(-1);
   }
   // io_write_ima_float("low_pass", H);
   // cout << "DATA NEW SIZE = " << TF_ExtData.nl() << "  " << TF_ExtData.nc() << endl;
   //cout << "H = " << H.nl() << "  " << H.nc() << endl;
   // cout << "MAXH = " << H.min() << " " << H.max() << endl;
   // TabWT[0].alloc(Nl,Nc,"scale");
   TabWT[0] = TF_ExtData;
   if (Extend == True)
   {
      get_hfilter(H, DNL, DNC);
      for (i=0; i < Nl; i++)
      for (j=0; j < Nc; j++) (TabWT[0])(i,j) *= H(i,j);
   }
   //cout << "ITER" << endl;
   // io_write_ima_complex_f("xx_in", TF_ExtData);
  
   for (s=0; s < NbrScale-1; s++)
   {
      Nl = TabNl(s);
      Nc = TabNc(s);   
      Nls = TabNl(s+1);
      Ncs = TabNc(s+1);
      if (Verbose == True)
             cout << "Scale " << s+1 << " Nl = " << Nl << " Nc = " << Nc <<  endl;
      H.resize(Nls,Ncs);
      TF_ExtData.resize(Nls,Ncs);
      if (Extend == True) get_hfilter(H,DNL/2.,DNC/2.);
      else get_hfilter(H);
//        {
//          char Name[256];
//          sprintf(Name, "low_%d" , s+1);
//          io_write_ima_float(Name, H);
//        }    
      // H.info("hfil");
      for (i=0; i < Nls; i++)
      for (j=0; j < Ncs; j++)
      {
         int Indi = i - Nls/2 + Nl/2;
	 int Indj = j - Ncs/2 + Nc/2;
	 TF_ExtData(i,j) = H(i,j) * (TabWT[s])(Indi, Indj);
	 (TabWT[s])(Indi, Indj) *= sqrt(1-H(i,j)*H(i,j));
      }
      // char Name[256];
      // sprintf(Name, "mband_%d" , s+1);
      // io_write_ima_complex_f(Name, TabWT[s]);
      TabWT[s+1].resize(Nl,Nc);
      TabWT[s+1] = TF_ExtData; 
      DNL =  DNL/2.;
      DNC =  DNC/2.;
   }
}
/*********************************************************************/

void MEYER_WT::recons_cf(Icomplex_f * &TabWT,  Icomplex_f & TF_ExtData)
{
   int i,j,s;
   int Nl = TabWT[0].nl();
   int Nc = TabWT[0].nc();
   int Nld2 = Nl/2;
   int Ncd2 = Nc/2;
   Ifloat H(Nl,Nc,"LowPass");
   Ifloat HB(Nl,Nc,"LowPass");

    // for  (s=0; s < NbrScale; s++) FFT2D.fftn2d(TabWT[s]);

   // cout << "REC ITER" << Nl << " " << Nc <<  endl;
   TF_ExtData.alloc(TabWT[0].nl(), TabWT[0].nc(), "rec");
   // cout << "DATA NEW SIZE = " << TF_ExtData.nl() << "  " << TF_ExtData.nc() << endl;
   double DNL=D_ExtNl;
   double DNC=D_ExtNc;
   
   if (Extend == True)
   {
     get_hfilter(HB,D_ExtNl,D_ExtNc);
     for (s=0; s < NbrScale-1; s++)
     {
      if (Verbose == True) cout << "Rec WT Scale " << s+1 << " " << TabNl(s) << " " << TabNc(s) << endl;
      Nl = TabNl(s);
      Nc = TabNc(s);
      int Nls = TabNl(s+1);
      int Ncs = TabNc(s+1);
      H.resize(Nls, Ncs);
      get_hfilter(H,DNL/2.,DNC/2.);
      for (i=0; i < Nl; i++)
      for (j=0; j < Nc; j++)   (TabWT[s])(i,j) *= HB(i,j);
      for (i=0; i < Nls; i++)
      for (j=0; j < Ncs; j++)  (TabWT[s])(i - Nls/2 + Nl/2, j - Ncs/2 + Nc/2) *= sqrt(1-H(i,j)*H(i,j));
      for (i=0; i < Nl; i++)
      for (j=0; j < Nc; j++)   TF_ExtData(i - Nl/2 + Nld2, j - Nc/2 + Ncd2) += (TabWT[s])(i, j);
 
      // char Name[256];
      // sprintf(Name, "mrband_%d" , s+1);
      // io_write_ima_complex_f(Name, TabWT[s]);
      if (s != NbrScale-2) HB = H;
      DNL =  DNL/2.;
      DNC =  DNC/2.;
     }
      s = NbrScale-1;
      Nl = TabWT[s].nl();
      Nc = TabWT[s].nc();
     for (i=0; i < Nl; i++)
     for (j=0; j < Nc; j++)  TF_ExtData(i - Nl/2 + Nld2, j - Nc/2 + Ncd2) += (TabWT[s])(i, j)*H(i,j);
  }
  else
  {
     s = NbrScale-1;
     Nl = TabNl(s);
     Nc = TabNc(s);
     TF_ExtData.resize(Nl, Nc);
     TF_ExtData = TabWT[s];
     for (s=NbrScale-2; s >= 0; s--)
     {
        if (Verbose == True) cout << "Rec WT Scale " << s+1 << " " << TabNl(s) << " " << TabNc(s) << endl;
        Nl = TabNl(s);
        Nc = TabNc(s);
        int Nls = TabNl(s+1);
        int Ncs = TabNc(s+1);
        H.resize(Nls, Ncs);
 	get_hfilter(H);
	// io_write_ima_float("low_pass1.fits", H);
	for (i=0; i < Nls; i++)
        for (j=0; j < Ncs; j++)
        {
           int Indi = i - Nls/2 + Nl/2;
	   int Indj = j - Ncs/2 + Nc/2;
	   (TabWT[s])(Indi, Indj) *= sqrt(1-H(i,j)*H(i,j));
 	    TF_ExtData(i,j) *= H(i,j);
 	   (TabWT[s])(Indi, Indj) += TF_ExtData(i,j);
	}   
	TF_ExtData.resize(Nl,Nc);
	TF_ExtData = TabWT[s];
	// DNL =  DNL/2.;
        // DNC =  DNC/2.;
     }
  }
  // io_write_ima_complex_f("reccf", TF_ExtData);
}
 
/*********************************************************************/

void MEYER_WT::ifft_tabima(Icomplex_f * & TabCF_Ima, Ifloat * & Tab_Ima, Bool Alloc)
{
    if (Alloc == True) Tab_Ima = new Ifloat[NbrScale];
    for  (int s=0; s < NbrScale; s++)
    {
       float Norm =  sqrt((float)(TabNl(s) * TabNc(s)));
       if (Alloc == True) Tab_Ima[s].alloc(TabNl(s), TabNc(s));
       for (int i=0; i < TabNl(s); i++)
       for (int j=0; j < TabNc(s); j++) (TabCF_Ima[s])(i,j) *= Norm;
       FFT2D.ifftn2d(TabCF_Ima[s]);
       for (int i=0; i < TabNl(s); i++)
       for (int j=0; j < TabNc(s); j++) (Tab_Ima[s])(i,j) = (TabCF_Ima[s])(i,j).real();
       
//        {
//          char Name[256];
//          sprintf(Name, "iband_%d" , s+1);
//          cout << "Scale " << s+1 << " " << Name << endl;
//          io_write_ima_complex_f(Name, TabCF_Ima[s]);
//        }    
    }
}

/*********************************************************************/

void MEYER_WT::fft_tabima(Ifloat * & Tab_Ima, Icomplex_f * & TabCF_Ima)
{
    // if (Alloc == True)  TabCF_Ima = new Icomplex_f[NbrScale];
    for  (int s=0; s < NbrScale; s++)
    {
       float Norm = sqrt( (float) TabNl(s) * TabNc(s));
       // if (Alloc == True) TabCF_Ima[s].alloc(TabNl(s), TabNc(s));
       FFT2D.fftn2d(Tab_Ima[s], TabCF_Ima[s]);
       for (int i=0; i < TabNl(s); i++)
       for (int j=0; j < TabNc(s); j++) (TabCF_Ima[s])(i,j) /= Norm;
    }
}

/*********************************************************************/
