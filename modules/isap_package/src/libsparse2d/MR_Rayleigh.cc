/******************************************************************************
**                   Copyright (C) 1997 by CEA
*******************************************************************************
**
**    UNIT
**
**    Version: 1.0
**
**    Author: Yanling Fang && Y. Bobichon and J.L. Starck
**
**    Date:  23.01.98
**    
**    File:  MR_Rayleigh.cc
**
*******************************************************************************
**
**    DESCRIPTION  Class for  Rayleigh noise management
**    ----------- 
**                 
******************************************************************************/

#include "MR_Obj.h"
#include "IM_Math.h"
#include "MR_Rayleigh.h"

/****************************************************************************/

StatRayleigh::StatRayleigh(type_noise_speckle TypeNoise, int Nscale, int Nima)
{
   alloc(TypeNoise, Nscale, Nima, TO_PAVE_BSPLINE);
}
/****************************************************************************/

StatRayleigh::StatRayleigh(type_noise_speckle TypeNoise, int Nscale, int Nima,
                           type_transform Trans)
{
   alloc(TypeNoise, Nscale, Nima, Trans);
} 
/****************************************************************************/

void StatRayleigh::alloc(type_noise_speckle TypeNoise, int Nscale, int Nima,
                           type_transform Trans)
{
  int N = 512; // size of image noise model
  int i,j,k,s;
  Ifloat NoiseIma(N,N,"noiseima");
  fltarray ThresholdMin, ThresholdMax ;
  //float sigma2;
  //float alpha=1.;
    
  Verbose = False;
  StatNoise=TypeNoise;
  NbrScale = Nscale;
  NbrImage = Nima;
  Transform = Trans;
  
  // compute the wavelet transform of  noise model
  MultiResol MR_Data (N,N, NbrScale, Transform,"MR_Data");
  NbrBand = MR_Data.nbr_band();
  
  StepHisto.alloc(NbrBand-1);  
  MinHisto.alloc(NbrBand-1);  
  MaxHisto.alloc(NbrBand-1);       
  
  Histo.alloc(NbrBand-1, SIZE_RAY_HISTO);
  Repart.alloc(NbrBand-1, SIZE_RAY_HISTO);
  Tab_Bin.alloc(NbrBand-1, SIZE_RAY_HISTO);
  
 
  switch (StatNoise)
  {
     case NOISE_LAPLACE: 
	    im_noise_laplace(NoiseIma,  NbrImage);
            break;
     case NOISE_RAYLEIGH: 
     case NOISE_LOG_RAYLEIGH:
            im_noise_rayleigh(NoiseIma,  NbrImage);
 	    break;
     default: cerr << "unknown Noise type" << endl;
            exit(-1);
            break;
  }
     
  // INFO(NoiseIma, "NoiseIma before");
  if (StatNoise == NOISE_LOG_RAYLEIGH)  
  {
        for(i=0; i < N;i++) 
	for(j=0; j < N;j++) NoiseIma(i,j)  = (NoiseIma(i,j) < FLOAT_EPSILON) ? 0: log (NoiseIma(i,j));
  }
  // INFO(NoiseIma, "NoiseIma");
 
  MR_Data.transform(NoiseIma);
   
  for (s=0; s < NbrBand-1; s++)
    {
      // Compute min max
      MinHisto(s)=min(MR_Data.band(s));
      MaxHisto(s)=max(MR_Data.band(s));
      
      if (Verbose == True) cout << "scale " << s << ":  Histogram [min,max] = [" << MinHisto(s) << "," << MaxHisto(s) << "]" << endl;
      
      // Compute histogram bin size  
      StepHisto(s) = (MaxHisto(s) - MinHisto(s))/(float)(SIZE_RAY_HISTO-1);

      // initialization
      for (k=0;k<SIZE_RAY_HISTO;k++)
	{
	  Tab_Bin(s,k) = MinHisto(s) + (float) k * StepHisto(s);
	  Histo(s,k)=0.;
	  Repart(s,k)= 0.;
	}
      
      // compute histogram 
      for(i=0;i<N;i++)
	{
	  for(j=0;j<N;j++)
	    { 
	      k = (int)  ((MR_Data(s,i,j)-MinHisto(s))/StepHisto(s)) ;
	      if ((k < 0) || (k >= SIZE_RAY_HISTO))
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
	}

      // histogram normalisation and repartition function
      for (k=0;k<SIZE_RAY_HISTO;k++)
	{
	  Histo(s,k)/=(N*N);
	  if (k==0) Repart(s,k)=Histo(s,k);
	  else Repart(s,k) = Repart(s,k-1)+Histo(s,k);
	}
      
    }
}

/****************************************************************************/

float StatRayleigh::prob(int Band, float Val)
{
  int k = (int)  ( (Val-Tab_Bin(Band,0)) / (Tab_Bin(Band,1)-Tab_Bin(Band,0)) );
  if ((k>=SIZE_RAY_HISTO)|| (k< 0))  return 0.;
  else return Histo(Band,k);
}

/****************************************************************************/

float StatRayleigh::repartition(int Band, float Val)
{
  int k = (int)  ( (Val-Tab_Bin(Band,0)) / (Tab_Bin(Band,1)-Tab_Bin(Band,0)) );
  if (k>=SIZE_RAY_HISTO)  return 1.;
  else if (k<0)  return 0.;
  else return Repart(Band,k);
}

/****************************************************************************/ 

void StatRayleigh::find_threshold(int NbrMaxBand, fltarray &Tab_Proba, 
                       fltarray &ThresholdMin, fltarray &ThresholdMax)
{
int k;

for (int s = 0; s < NbrMaxBand-1; s++)
  {
    k=0;
    while (Repart(s,k) < Tab_Proba(s)) k++;
    ThresholdMin(s) = Tab_Bin(s,k);
    while (Repart(s,k) < (1.-Tab_Proba(s))) k++;
    ThresholdMax(s) = Tab_Bin(s,k);
  }
}  

/****************************************************************************/

void StatRayleigh::write()
{
  fits_write_fltarr("Tab_Histo", Histo);
  fits_write_fltarr("Tab_Repart", Repart);
  fits_write_fltarr("Tab_Bin", Tab_Bin);
}

/****************************************************************************/

