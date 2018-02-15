/******************************************************************************
**                   Copyright (C) 1999 by CEA
*******************************************************************************
**
**    UNIT
**
**    Version: 1.0
**
**    Author: J.L. Starck
**
**    Date:  1.03.99
**    
**    File:  MR_CorrNoise.cc
**
*******************************************************************************
**
**    DESCRIPTION  Class for correlated noise management
**    ----------- 
**                 
******************************************************************************/

#include "MR_Obj.h"
#include "IM_Math.h"
#include "MR_CorrNoise.h"

/****************************************************************************/

StatNoiseMap::StatNoiseMap(Ifloat &ImaNoise, int Nscale, type_transform Trans)
{
   // cout << "StatNoiseMap = " << "Nl = " << ImaNoise.nl() << " Nc = " << ImaNoise.nc() << " Nelem = " << ImaNoise.n_elem() << endl;
   alloc(ImaNoise, Nscale, Trans);
} 
/****************************************************************************/

StatNoiseMap::StatNoiseMap(Ifloat &NoiseIma, int Nscale, type_transform Trans, FilterAnaSynt *FAS, 
                         sb_type_norm Norm, int NbrUndec, type_undec_filter U_Filter)
{
   alloc(NoiseIma, Nscale, Trans, FAS,Norm,NbrUndec,U_Filter);
}
/****************************************************************************/
 
void StatNoiseMap::alloc(Ifloat &NoiseIma, int Nscale, type_transform Trans, FilterAnaSynt *FAS, 
                         sb_type_norm Norm, int NbrUndec, type_undec_filter U_Filter)
{
  int i,j,k,s;
  
  if (NoiseIma.n_elem() == 0)
  {
     cout << "Error in StatNoiseMap::alloc: incorrect NoiseIma parameter ... " << endl;
     cout << "      NoiseIma::Nl = " << NoiseIma.nl() << " Nc = " << NoiseIma.nc() << endl;
     exit(-1);
  }
  fltarray ThresholdMin, ThresholdMax ;
  Verbose = False;
  NbrScale = Nscale;
  Transform = Trans;
  
  // compute the wavelet transform of  noise model
  MultiResol MR_Data;
  MR_Data.alloc(NoiseIma.nl(), NoiseIma.nc(), NbrScale, Transform, FAS, Norm, NbrUndec, U_Filter);
  NbrBand = MR_Data.nbr_band();
  
  StepHisto.alloc(NbrBand-1);  
  MinHisto.alloc(NbrBand-1);  
  MaxHisto.alloc(NbrBand-1);       
  TabSigma.alloc(NbrBand-1);
  
  Histo.alloc(NbrBand-1, SIZE_CORR_HISTO);
  Repart.alloc(NbrBand-1, SIZE_CORR_HISTO);
  Tab_Bin.alloc(NbrBand-1, SIZE_CORR_HISTO);
  MR_Data.transform(NoiseIma);
   
  for (s=0; s < NbrBand-1; s++)
  {
      int Nlb = MR_Data.size_band_nl(s);
      int Ncb = MR_Data.size_band_nc(s);
      // Compute min max
      MinHisto(s)=min(MR_Data.band(s));
      MaxHisto(s)=max(MR_Data.band(s));
      TabSigma(s)=sigma(MR_Data.band(s));
      if (Verbose == True) cout << "band " << s << ":  Histogram [min,max] = [" << MinHisto(s) << "," << MaxHisto(s) << "]" << endl;
      
      // Compute histogram bin size  
      StepHisto(s) = (MaxHisto(s) - MinHisto(s))/(float)(SIZE_CORR_HISTO-1);

      // initialization
      for (k=0;k<SIZE_CORR_HISTO;k++)
      {
	  Tab_Bin(s,k) = MinHisto(s) + (float) k * StepHisto(s);
	  Histo(s,k)=0.;
	  Repart(s,k)= 0.;
      }
      
      // compute histogram 
      for(i=0; i< Nlb;i++)
      for(j=0; j< Ncb;j++)
      { 
	  k = (int)  ((MR_Data(s,i,j)-MinHisto(s))/StepHisto(s)) ;
	  if ((k < 0) || (k >= SIZE_CORR_HISTO))
	  {
	      cout << "Error: k = " << k << " in band " << s << endl;
              cout << "min = " << MinHisto(s) << endl;
	      cout << "max = " << MaxHisto(s) << endl;
	      cout << "step = " << StepHisto(s) << endl;
	      cout << "val = " <<  MR_Data(s,i,j) << endl;
	      exit(-1);
	  }
	  Histo(s,k)+=1; 
      }
 
      // histogram normalisation and repartition function
      for (k=0; k < SIZE_CORR_HISTO; k++)
      { 
	  Histo(s,k)/=(Nlb*Ncb);
	  if (k==0) Repart(s,k)=Histo(s,k);
	  else Repart(s,k) = Repart(s,k-1)+Histo(s,k);
      }
   }
}

/****************************************************************************/

float StatNoiseMap::prob(int Band, float Val)
{
  int k = (int)  ( (Val-Tab_Bin(Band,0)) / (Tab_Bin(Band,1)-Tab_Bin(Band,0)) );
  if ((k>=SIZE_CORR_HISTO)|| (k< 0))  return 0.;
  else return Histo(Band,k);
}

/****************************************************************************/

float StatNoiseMap::repartition(int Band, float Val)
{
  int k = (int)  ( (Val-Tab_Bin(Band,0)) / (Tab_Bin(Band,1)-Tab_Bin(Band,0)) );
  if (k>=SIZE_CORR_HISTO)  return 1.;
  else if (k<0)  return 0.;
  else return Repart(Band,k);
}

/****************************************************************************/ 

void StatNoiseMap::find_threshold(int NbrMaxBand, fltarray &Tab_Proba, 
                       fltarray &ThresholdMin, fltarray &ThresholdMax)
{
  int k;

  for (int s = 0; s < NbrMaxBand-1; s++)
  {
    k=0;
    while ((k < SIZE_CORR_HISTO) && (Repart(s,k) < Tab_Proba(s))) k++;
    if (k >= SIZE_CORR_HISTO) k = SIZE_CORR_HISTO-1; 
    ThresholdMin(s) = Tab_Bin(s,k);
    
    while ((k < SIZE_CORR_HISTO) && (Repart(s,k) < (1.-Tab_Proba(s)))) k++;
    if (k >= SIZE_CORR_HISTO) k = SIZE_CORR_HISTO-1; 
    ThresholdMax(s) = Tab_Bin(s,k);
  }
}  

/****************************************************************************/

void StatNoiseMap::write()
{
  fits_write_fltarr("Tab_Histo", Histo);
  fits_write_fltarr("Tab_Repart", Repart);
  fits_write_fltarr("Tab_Bin", Tab_Bin);
}

/****************************************************************************/

