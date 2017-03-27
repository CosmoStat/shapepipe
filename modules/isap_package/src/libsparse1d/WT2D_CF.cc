/*******************************************************************************
**
**    UNIT
**
**    Version: 1.0
**
**    Author: Jean-Luc Starck
**
**    Date:  25/09/2003 
**    
**    File:  WT2D_CF.cc
**
**    Modification history:
**
*******************************************************************************
**
**    DESCRIPTION  1D Sub-Band decomposition
**    -----------  
**                 
******************************************************************************/

#include "SB_Filter1D.h"
#include "SB_Filter.h"
#include "WT2D_CF.h"
#include "IM_IO.h"

#define Nbr_CF_H0_ODD 13
static float  F_cf_h0_odd[] = {-0.0017581,0,0.0222656,-0.0468750,
         -0.0482422,0.2968750,0.5554688,0.2968750,-0.0482422,
	 -0.0468750,0.0222656,0,-0.0017581};

#define Nbr_CF_G0_ODD 19
static float  F_cf_g0_odd[] = {-0.0000706,0,0.0013419,-0.0018834,-0.0071568,
0.0238560,0.0556431,-0.0516881,-0.2997576,0.5594308,-0.2997576,
-0.0516881,0.0556431,0.0238560,-0.0071568,-0.0018834,0.0013419,0,-0.0000706};

#define Nbr_CF_H0_EVEN 13
static float  F_cf_h0_even[] = {0, -0.0058109,0.0166977,-0.0000641,-0.0834914,
    0.0919537,0.4807151,0.4807151,0.0919537,
    -0.0834914,-0.0000641,0.0166977,-0.0058109};

#define Nbr_CF_G0_EVEN 17
static float  F_cf_g0_even[] = {-0.0004645 ,0.0013349,0.002206,-0.0130127,
0.0015360,0.0869008,0.0833552,-0.4885957,
0.4885957,-0.0833552,-0.0869008, -0.0015360,
0.0130127,-0.002206,-0.0013349,0.0004645, 0};

/***********************************************************************/

void WTCF2D::init()
{ 
cout << "INIT WTCF2D" << endl;
    float *Heven = F_cf_h0_even;
    float *Geven = F_cf_g0_even;
    int Nheven = Nbr_CF_H0_EVEN;
    int Ngeven = Nbr_CF_G0_EVEN;
    float *Hodd = F_cf_h0_odd;
    float *Godd = F_cf_g0_odd;
    int Nhodd = Nbr_CF_H0_ODD;
    int Ngodd = Nbr_CF_G0_ODD;
    sb_type_norm FilNorm = NORM_L2;
    SBF_EVEN = new SubBandFilter(Heven,Nheven,Geven,Ngeven,FilNorm);
    SBF_ODD = new SubBandFilter(Hodd,Nhodd,Godd,Ngodd,FilNorm);      
    TreeA = new Ortho_2D_WT(*SBF_EVEN,*SBF_EVEN);
    TreeC = new Ortho_2D_WT(*SBF_ODD,*SBF_EVEN);
    TreeB = new Ortho_2D_WT(*SBF_EVEN,*SBF_ODD);
    TreeD = new Ortho_2D_WT(*SBF_ODD,*SBF_ODD);   
    WT = new HALF_DECIMATED_2D_WT(*SBF_ODD);
    Verbose = False;
//     TreeA->Line_SubSample_H_Even = True;
//     TreeA->Line_SubSample_G_Odd = False;
//     TreeA->Col_SubSample_H_Even = True;
//     TreeA->Col_SubSample_G_Odd = False;
//     
//     TreeB->Line_SubSample_H_Even = True;
//     TreeB->Line_SubSample_G_Odd = False;
//     TreeB->Col_SubSample_H_Even = False;
//     TreeB->Col_SubSample_G_Odd = True;
//     
//     TreeC->Line_SubSample_H_Even = False;
//     TreeC->Line_SubSample_G_Odd = True;
//     TreeC->Col_SubSample_H_Even = True;
//     TreeC->Col_SubSample_G_Odd = False;
//     
//     TreeD->Line_SubSample_H_Even = False;
//     TreeD->Line_SubSample_G_Odd = True;
//     TreeD->Col_SubSample_H_Even = False;
//     TreeD->Col_SubSample_G_Odd = True;
}
 
/***********************************************************************/

void WTCF2D::transform(Ifloat & Data,  Ifloat * & TabTrans, int Nbr_Plan)
{ 
   int i,j,Nl = Data.nl();
   int Nc = Data.nc();
   int NlA = (Nl+1)/2;
   int NcA = (Nc+1)/2;
   Ifloat ImaA(NlA,NcA,"imaa");
   int NlB = Nl/2;
   int NcB = (Nc+1)/2;
   Ifloat ImaB(NlB,NcB,"imab");
   int NlC = (Nl+1)/2;
   int NcC = Nc/2;
   Ifloat ImaC(NlC,NcC,"imac");
   int NlD = Nl/2;
   int NcD = Nc/2;
   Ifloat ImaD(NlD,NcD,"imad");
   int NbrBand = WT->alloc(TabTrans, Nl, Nc, 2, 1);

   cout << "Transformation: NbrScale = " << Nbr_Plan << endl;
   WT->transform(Data,TabTrans, 2, 1);
      
   for (i=0; i < Nl; i+=2)
   for (j=0; j < Nc; j+=2)
   {
       ImaA(i/2,j/2) = (TabTrans[3])(i,j);
       if (i+1 < Nl) ImaB(i/2,j/2) = (TabTrans[3])(i+1,j);
       if (j+1 < Nc) ImaC(i/2,j/2) = (TabTrans[3])(i,j+1);
       if ((i+1 < Nl) && (j+1 < Nc)) ImaD(i/2,j/2) = (TabTrans[3])(i+1,j+1);
   }
   io_write_ima_float((char*)"xx_a.fits", ImaA); 
   io_write_ima_float((char*)"xx_b.fits", ImaB); 
   io_write_ima_float((char*)"xx_c.fits", ImaC);
   io_write_ima_float((char*)"xx_d.fits", ImaD);
   
   TreeA->transform(ImaA, Nbr_Plan-1);
   TreeB->transform(ImaB, Nbr_Plan-1);
   TreeC->transform(ImaC, Nbr_Plan-1);
   TreeD->transform(ImaD, Nbr_Plan-1);   
   
   io_write_ima_float((char*)"xx_ta.fits", ImaA); 
   io_write_ima_float((char*)"xx_tb.fits", ImaB); 
   io_write_ima_float((char*)"xx_tc.fits", ImaC);
   io_write_ima_float((char*)"xx_td.fits", ImaD);
   
   for (i=0; i < Nl; i+=2)
   for (j=0; j < Nc; j+=2)
   {
       (TabTrans[3])(i,j) =  ImaA(i/2,j/2);
       if (i+1 < Nl) (TabTrans[3])(i+1,j) = ImaB(i/2,j/2);
       if (j+1 < Nc) (TabTrans[3])(i,j+1) = ImaC(i/2,j/2);
       if ((i+1 < Nl) && (j+1 < Nc)) (TabTrans[3])(i+1,j+1) = ImaD(i/2,j/2);
   }
   // io_write_ima_float((char*)"xx_b4.fits", TabTrans[2]);
}
   
/***********************************************************************/

void WTCF2D::alloc_tabz(int Nl, int Nc, Ifloat * & TabZre, Ifloat * & TabZim, int Nbr_Plan, int & NzBand)
{ 
   int i,j,s,b;
   int NlA = (Nl+1)/2;
   int NcA = (Nc+1)/2;
   int NlB = Nl/2;
   int NcB = (Nc+1)/2;
   int NlC = (Nl+1)/2;
   int NcC = Nc/2;
   int NlD = Nl/2;
   int NcD = Nc/2;    
   int NbrBand_per_Resol = 6;
   char ch[80];
   NzBand = NbrBand_per_Resol*(Nbr_Plan-1)+2;
   if (Verbose == True) cout << " Nbr Z band = " << NzBand << endl;
   TabZre = new Ifloat [NzBand];
   TabZim = new Ifloat [NzBand];
   int Nls = NlA;
   int Ncs = NcA;
   int Cpt = 0;
   int IndBand = 0;
   for (s = 0; s < NzBand-2; s+=NbrBand_per_Resol)
   {
      for (b = 0; b < NbrBand_per_Resol; b++)
      {
        sprintf (ch, "band_%d", s+b);
        if (Verbose == True) cout << "Alloc band " << s << " " << s+b+1 << " Nls = " << Nls << " Ncs = " << Ncs  << endl;
        Cpt += 2*Nls*Ncs;
        TabZre[s+b].alloc(Nls, Ncs, ch);
        TabZim[s+b].alloc(Nls, Ncs, ch);
	
	IndBand++;
      }
      if (IndBand+NbrBand_per_Resol < NzBand)
      {
         Nls /= 2;
         Ncs /= 2;
      }
   }
   if (Verbose == True) cout << "Alloc band " << NzBand-1 << " Nls = " << Nls << " Ncs = " << Ncs  << endl;
   TabZre[NzBand-2].alloc(Nls, Ncs, ch);
   TabZim[NzBand-2].alloc(Nls, Ncs, ch); 
   if (Verbose == True) cout << "Alloc band " << NzBand << " Nls = " << Nls << " Ncs = " << Ncs  << endl;
   TabZre[NzBand-1].alloc(Nls, Ncs, ch);
   TabZim[NzBand-1].alloc(Nls, Ncs, ch);
   Cpt += 4*Nls*Ncs;
   // cout << "R = " << (float) Cpt / (float)(Nl*Nc) << " Cpt = " << Cpt << " Nls*Ncs = " <<  Nls*Ncs  << endl;
}

/***********************************************************************/

void WTCF2D::coeff_to_z(Ifloat *TabTrans, Ifloat *TabZre, Ifloat *TabZim, int Nbr_Plan, int NzBand)
{
   int i,j,s;
   int Nl = TabTrans[0].nl();
   int Nc = TabTrans[0].nc();
   int NlA = (Nl+1)/2;
   int NcA = (Nc+1)/2;
   int NlB = Nl/2;
   int NcB = (Nc+1)/2;
   int NlC = (Nl+1)/2;
   int NcC = Nc/2;
   int NlD = Nl/2;
   int NcD = Nc/2;    
   int NbrBand_per_Resol = 6;
   float z1r,z1i,z2r,z2i,a,b,c,d;
   char ch[80];
   
   for (s=0; s < 3; s++)
   {
       // cout << " Resol 1, band " << s+1 << endl;
       for (i=0; i < Nl; i+=2)
       for (j=0; j < Nc; j+=2)
       {
          a = (TabTrans[s])(i,j);
	  b = (TabTrans[s])(i+1,j, I_ZERO);
 	  c = (TabTrans[s])(i,j+1, I_ZERO);
	  d = (TabTrans[s])(i+1,j+1, I_ZERO);
	  abcd_to_z(a,b,c,d,z1r,z1i,z2r,z2i);
	  (TabZre[2*s])(i/2,j/2) = z1r;
	  (TabZim[2*s])(i/2,j/2) = z1i;
	  (TabZre[2*s+1])(i/2,j/2) = z2r;
	  (TabZim[2*s+1])(i/2,j/2) = z2i;
      }        
   }
   // cout << "First scale OK " << endl;
   int Nls = Nl;
   int Ncs = Nc;
   int Nl_2,Nc_2;
   s = 3;
   int IndZ1 = 6;
   int IndZ2 = 7;
   for (int k = 0; k < Nbr_Plan-2; k++)
   {
      // cout << "Nz =  " << NzBand << " Z1 =  " << IndZ1 <<  " Z2 =  " << IndZ2 <<endl;
      // cout << "Nls = " << TabTrans[s].nl() << " Ncs = " << TabTrans[s].nc() <<endl;
      Nl_2 = (Nls+1)/2;
      Nc_2 = (Ncs+1)/2;
      for (i=0; i < Nl_2; i+=2)
      for (j=0; j < Nc_2; j+=2)
      {
          a = (TabTrans[s])(i+Nl_2,j);
	  b = (TabTrans[s])(i+1+Nl_2,j, I_ZERO);
 	  c = (TabTrans[s])(i+Nl_2,j+1, I_ZERO);
	  d = (TabTrans[s])(i+1+Nl_2,j+1, I_ZERO);
	  abcd_to_z(a,b,c,d,z1r,z1i,z2r,z2i);
	  (TabZre[IndZ1])(i/2,j/2) = z1r;
	  (TabZim[IndZ1])(i/2,j/2) = z1i;
	  (TabZre[IndZ2])(i/2,j/2) = z2r;
	  (TabZim[IndZ2])(i/2,j/2) = z2i;
      }
      IndZ1+=2;
      IndZ2+=2;
      // cout << "Nz =  " << NzBand << " Z1 =  " << IndZ1 <<  " Z2 =  " << IndZ2 <<endl;
      for (i=0; i < Nl_2; i+=2)
      for (j=0; j < Nc_2; j+=2)
      {
          a = (TabTrans[s])(i,j+Nc_2);
	  b = (TabTrans[s])(i+1,j+Nc_2, I_ZERO);
 	  c = (TabTrans[s])(i,j+1+Nc_2, I_ZERO);
	  d = (TabTrans[s])(i+1,j+1+Nc_2, I_ZERO);
	  abcd_to_z(a,b,c,d,z1r,z1i,z2r,z2i);
	  (TabZre[IndZ1])(i/2,j/2) = z1r;
	  (TabZim[IndZ1])(i/2,j/2) = z1i;
	  (TabZre[IndZ2])(i/2,j/2) = z2r;
	  (TabZim[IndZ2])(i/2,j/2) = z2i;
      }
      IndZ1+=2;
      IndZ2+=2;
      // cout << "Nz =  " << NzBand << " Z1 =  " << IndZ1 <<  " Z2 =  " << IndZ2 <<endl;
      for (i=0; i < Nl_2; i+=2)
      for (j=0; j < Nc_2; j+=2)
      {
          a = (TabTrans[s])(i+Nl_2,j+Nc_2);
	  b = (TabTrans[s])(i+1+Nl_2,j+Nc_2, I_ZERO);
 	  c = (TabTrans[s])(i+Nl_2,j+1+Nc_2, I_ZERO);
	  d = (TabTrans[s])(i+1+Nl_2,j+1+Nc_2, I_ZERO);
	  abcd_to_z(a,b,c,d,z1r,z1i,z2r,z2i);
	  (TabZre[IndZ1])(i/2,j/2) = z1r;
	  (TabZim[IndZ1])(i/2,j/2) = z1i;
	  (TabZre[IndZ2])(i/2,j/2) = z2r;
	  (TabZim[IndZ2])(i/2,j/2) = z2i;
      }
      IndZ1+=2;
      IndZ2+=2;
      if (k != Nbr_Plan-3) 
      {
         Nls = (Nls+1)/2;
         Ncs = (Ncs+1)/2;
      }
   }
   // cout << "Last scale OK " << endl;
   // cout << "Nz =  " << NzBand << " Z1 =  " << IndZ1 <<  " Z2 =  " << IndZ2 <<endl;
  
   for (i=0; i < Nl_2; i+=2)
   for (j=0; j < Nc_2; j+=2)
   {
       a = (TabTrans[s])(i,j);
       b = (TabTrans[s])(i+1,j, I_ZERO);
       c = (TabTrans[s])(i,j+1, I_ZERO);
       d = (TabTrans[s])(i+1,j+1, I_ZERO);
       abcd_to_z(a,b,c,d,z1r,z1i,z2r,z2i);
       (TabZre[IndZ1])(i/2,j/2) = z1r;
       (TabZim[IndZ1])(i/2,j/2) = z1i;
       (TabZre[IndZ2])(i/2,j/2) = z2r;
       (TabZim[IndZ2])(i/2,j/2) = z2i;
    } 
}

/***********************************************************************/

void WTCF2D::threshold(Ifloat * & TabZre, Ifloat * & TabZim, int Nz, float SigmaNoise)
{
   int i,j,z;
   float NSigma=4.85;
   float Noise=SigmaNoise;
   for (z=0; z < Nz-2; z++)
   {
      float Level = (z < 6) ? (NSigma+0.5)*Noise: NSigma*Noise; 
      for (i=0; i < TabZre[z].nl(); i++)
      for (j=0; j < TabZre[z].nc(); j++)
      {
          float Re = (TabZre[z])(i,j);
	  float Im = (TabZim[z])(i,j);
          float Mod = sqrt(Re*Re+Im*Im);
	  if (Mod < Level) 
	  {
	     (TabZre[z])(i,j) = 0;
	     (TabZim[z])(i,j) = 0;
	  }
      }
   }
}

/***********************************************************************/

void WTCF2D::z_to_coeff(Ifloat * & TabZre, Ifloat * & TabZim, Ifloat *TabTrans, int Nbr_Plan)
{
   int i,j,s;
   int Nl = TabTrans[0].nl();
   int Nc = TabTrans[0].nc();
   int NlA = (Nl+1)/2;
   int NcA = (Nc+1)/2;
   int NlB = Nl/2;
   int NcB = (Nc+1)/2;
   int NlC = (Nl+1)/2;
   int NcC = Nc/2;
   int NlD = Nl/2;
   int NcD = Nc/2;    
   int NbrBand_per_Resol = 6;
   float z1r,z1i,z2r,z2i,a,b,c,d;
   char ch[80];
   
   for (s=0; s < 3; s++)
   {
       // cout << " Resol 1, band " << s+1 << endl;
       for (i=0; i < Nl; i+=2)
       for (j=0; j < Nc; j+=2)
       {
	  z1r = (TabZre[2*s])(i/2,j/2);
	  z1i = (TabZim[2*s])(i/2,j/2);
	  z2r = (TabZre[2*s+1])(i/2,j/2);
	  z2i = (TabZim[2*s+1])(i/2,j/2);
 	  z_to_abcd(z1r,z1i,z2r,z2i,a,b,c,d);
          (TabTrans[s])(i,j) = a;
	  (TabTrans[s])(i+1,j) = b;
	  (TabTrans[s])(i,j+1) = c;
	  (TabTrans[s])(i+1,j+1) = d;
      }        
   }
   // cout << "First scale OK " << endl;
   int Nls = Nl;
   int Ncs = Nc;
   int Nl_2,Nc_2;
   s = 3;
   int IndZ1 = 6;
   int IndZ2 = 7;
   for (int k = 0; k < Nbr_Plan-2; k++)
   {
      // cout <<  " Z1 =  " << IndZ1 <<  " Z2 =  " << IndZ2 <<endl;
      // cout << "Nls = " << TabTrans[s].nl() << " Ncs = " << TabTrans[s].nc() <<endl;
      Nl_2 = (Nls+1)/2;
      Nc_2 = (Ncs+1)/2;
      for (i=0; i < Nl_2; i+=2)
      for (j=0; j < Nc_2; j+=2)
      {
          z1r =  (TabZre[IndZ1])(i/2,j/2);
	  z1i =  (TabZim[IndZ1])(i/2,j/2);
	  z2r = (TabZre[IndZ2])(i/2,j/2);
	  z2i = (TabZim[IndZ2])(i/2,j/2);
	  z_to_abcd(z1r,z1i,z2r,z2i,a,b,c,d);
          (TabTrans[s])(i+Nl_2,j) = a;
	  (TabTrans[s])(i+1+Nl_2,j) = b;
 	  (TabTrans[s])(i+Nl_2,j+1) = c;
 	  (TabTrans[s])(i+1+Nl_2,j+1) = d;
      }
      IndZ1+=2;
      IndZ2+=2;
      // cout  << " Z1 =  " << IndZ1 <<  " Z2 =  " << IndZ2 <<endl;
      for (i=0; i < Nl_2; i+=2)
      for (j=0; j < Nc_2; j+=2)
      {
          z1r = (TabZre[IndZ1])(i/2,j/2);
	  z1i = (TabZim[IndZ1])(i/2,j/2);
	  z2r = (TabZre[IndZ2])(i/2,j/2);
	  z2i = (TabZim[IndZ2])(i/2,j/2);
	  z_to_abcd(z1r,z1i,z2r,z2i,a,b,c,d);
          (TabTrans[s])(i,j+Nc_2) = a ;
	  (TabTrans[s])(i+1,j+Nc_2) = b;
 	  (TabTrans[s])(i,j+1+Nc_2) = c;
	  (TabTrans[s])(i+1,j+1+Nc_2) = d;
      }
      IndZ1+=2;
      IndZ2+=2;
      // cout <<  " Z1 =  " << IndZ1 <<  " Z2 =  " << IndZ2 <<endl;
      for (i=0; i < Nl_2; i+=2)
      for (j=0; j < Nc_2; j+=2)
      {
          z1r = (TabZre[IndZ1])(i/2,j/2);
	  z1i = (TabZim[IndZ1])(i/2,j/2);
	  z2r = (TabZre[IndZ2])(i/2,j/2);
	  z2i = (TabZim[IndZ2])(i/2,j/2);
	  z_to_abcd(z1r,z1i,z2r,z2i,a,b,c,d);
          (TabTrans[s])(i+Nl_2,j+Nc_2) = a;
	  (TabTrans[s])(i+1+Nl_2,j+Nc_2) = b;
 	  (TabTrans[s])(i+Nl_2,j+1+Nc_2) = c;
	  (TabTrans[s])(i+1+Nl_2,j+1+Nc_2) = d;
      }
      IndZ1+=2;
      IndZ2+=2;
      if (k != Nbr_Plan-3) 
      {
         Nls = (Nls+1)/2;
         Ncs = (Ncs+1)/2;
      }
   }
   // cout << "Last scale OK " << endl;
   //cout <<  " Z1 =  " << IndZ1 <<  " Z2 =  " << IndZ2 <<endl;
  
   for (i=0; i < Nl_2; i+=2)
   for (j=0; j < Nc_2; j+=2)
   {
       z1r = (TabZre[IndZ1])(i/2,j/2);
       z1i = (TabZim[IndZ1])(i/2,j/2);
       z2r = (TabZre[IndZ2])(i/2,j/2);
       z2i = (TabZim[IndZ2])(i/2,j/2);
       z_to_abcd(z1r,z1i,z2r,z2i,a,b,c,d);
       (TabTrans[s])(i,j) = a;
       (TabTrans[s])(i+1,j) = b;
       (TabTrans[s])(i,j+1) = c;
       (TabTrans[s])(i+1,j+1) = d;
    } 
}

/***********************************************************************/

void WTCF2D::recons(Ifloat * & TabTrans, Ifloat & Data, int Nbr_Plan)
{ 
   int i,j,Nl = Data.nl();
   int Nc = Data.nc();
   int NlA = (Nl+1)/2;
   int NcA = (Nc+1)/2;
   Ifloat ImaA(NlA,NcA,"imaa");
   int NlB = Nl/2;
   int NcB = (Nc+1)/2;
   Ifloat ImaB(NlB,NcB,"imab");
   int NlC = (Nl+1)/2;
   int NcC = Nc/2;
   Ifloat ImaC(NlC,NcC,"imac");
   int NlD = Nl/2;
   int NcD = Nc/2;
   Ifloat ImaD(NlD,NcD,"imad");
   cout << "Reconstruction: NbrScale = " << Nbr_Plan << endl;
   for (i=0; i < Nl; i+=2)
   for (j=0; j < Nc; j+=2)
   {
       ImaA(i/2,j/2) = (TabTrans[3])(i,j);
       if (i+1 < Nl) ImaB(i/2,j/2) = (TabTrans[3])(i+1,j);
       if (j+1 < Nc) ImaC(i/2,j/2) = (TabTrans[3])(i,j+1);
       if ((i+1 < Nl) && (j+1 < Nc)) ImaD(i/2,j/2) = (TabTrans[3])(i+1,j+1);
   } 

   io_write_ima_float((char*)"xx_rta.fits", ImaA); 
   io_write_ima_float((char*)"xx_rtb.fits", ImaB); 
   io_write_ima_float((char*)"xx_rtc.fits", ImaC);
   io_write_ima_float((char*)"xx_rtd.fits", ImaD);
      
   TreeA->recons(ImaA, Nbr_Plan-1);
   TreeB->recons(ImaB, Nbr_Plan-1);
   TreeC->recons(ImaC, Nbr_Plan-1);
   TreeD->recons(ImaD, Nbr_Plan-1);
   
   io_write_ima_float((char*)"xx_ra.fits", ImaA); 
   io_write_ima_float((char*)"xx_rb.fits", ImaB); 
   io_write_ima_float((char*)"xx_rc.fits", ImaC);
   io_write_ima_float((char*)"xx_rd.fits", ImaD);
   
   for (i=0; i < Nl; i+=2)
   for (j=0; j < Nc; j+=2)
   {
       (TabTrans[3])(i,j) = ImaA(i/2,j/2);
       if (i+1 < Nl) (TabTrans[3])(i+1,j) = ImaB(i/2,j/2);
       if (j+1 < Nc) (TabTrans[3])(i,j+1) = ImaC(i/2,j/2);
       if ((i+1 < Nl) && (j+1 < Nc)) (TabTrans[3])(i+1,j+1) = ImaD(i/2,j/2);
   }         
   WT->recons(TabTrans, Data,2, 1);
}
/***********************************************************************/
