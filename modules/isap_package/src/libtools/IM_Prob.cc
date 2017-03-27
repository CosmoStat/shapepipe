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
**    File:  IM_Prob.cc
**
*******************************************************************************
**
**    DESCRIPTION  Class for correlated noise management
**    ----------- 
**                 
******************************************************************************/


#include "IM_Prob.h"
#include "NR.h"
using namespace std;


/****************************************************************************/

CImaProb::CImaProb()
{
    Verbose = False;
    Histo.alloc(SIZE_PROB_HISTO);
    Repart.alloc(SIZE_PROB_HISTO);
    Tab_Bin.alloc(SIZE_PROB_HISTO);
    StepHisto=MinHisto=MaxHisto=0;
} 

/****************************************************************************/
 
void CImaProb::set(float *NoiseIma, int N, bool sym)
{
	int i,k;

	symetry=sym;
//cerr<<"Symetry = "<<symetry<<endl;
	if (N == 0)
	{
		cout << "Error in CImaProb::alloc: incorrect NoiseIma parameter ... " << endl;
		exit(-1);
	}
	MinHisto = *NoiseIma;
	MaxHisto = *NoiseIma;
	for (i=1; i < N; i++)
	{
		if (MinHisto > (double) NoiseIma[i]) MinHisto = (double) NoiseIma[i];
		if (MaxHisto < (double) NoiseIma[i]) MaxHisto = (double) NoiseIma[i];
	}
	
	if(symetry)
	{
		MaxHisto = max(abs(MaxHisto),abs(MinHisto));
		MinHisto= 0.0;
	}

	if (Verbose == True) cout << "Ima :  Histogram [min,max] = [" << MinHisto << "," << MaxHisto << "]" << endl;

	// Compute histogram bin size  
	StepHisto = (MaxHisto - MinHisto)/(double)(SIZE_PROB_HISTO-1);

	// initialization
	for (k=0;k<SIZE_PROB_HISTO;k++)
	{
		Tab_Bin(k) = MinHisto + (double) k * StepHisto;
		Histo(k)=0.;
		Repart(k)= 0.;
	}

	// compute histogram 
	if(!symetry)
		for(i=0; i< N;i++)
		{ 
			k = (int)  ((NoiseIma[i]-MinHisto)/StepHisto) ;
			if ((k < 0) || (k >= SIZE_PROB_HISTO))
			{
				cout << "Error: k = " << k << " in band " <<  endl;
				cout << "min = " << MinHisto << endl;
				cout << "max = " << MaxHisto << endl;
				cout << "step = " << StepHisto << endl;
				cout << "val = " <<  NoiseIma[i] << endl;
				exit(-1);
			}
			Histo(k) += 1; 
		}
	else// symetry
		for(i=0; i< N;i++)
		{ 
			k = (int)  (abs(NoiseIma[i])/StepHisto) ;
			if ((k < 0) || (k >= SIZE_PROB_HISTO))
			{
				cout << "Error: k = " << k << " in band " <<  endl;
				cout << "min = " << MinHisto << endl;
				cout << "max = " << MaxHisto << endl;
				cout << "step = " << StepHisto << endl;
				cout << "val = " <<  NoiseIma[i] << endl;
				exit(-1);
			}
			Histo(k) += 1; 
		}

	// histogram normalisation and repartition function
	for (k=0; k < SIZE_PROB_HISTO; k++)
	{ 
		//cerr<<"histo("<<k<<")="<<Tab_Bin(k)/3.<<", nb="<<Histo(k);
		Histo(k) /= (double) N;
		if (k==0) Repart(k)=Histo(k);
		else Repart(k) = Repart(k-1)+Histo(k);
		//cerr<<", repart %< ="<<Repart(k)<<", #above threshold="<<N*(1-Repart(k))<<endl;
	}
}

/****************************************************************************/

void CImaProb::set(Ifloat &NoiseIma)
{
   set (NoiseIma.buffer(),  NoiseIma.n_elem());
}

/****************************************************************************/

double CImaProb::prob(double Val)
{
  int k = (int)  ( (Val-Tab_Bin(0)) / (Tab_Bin(1)-Tab_Bin(0)) );
  if ((k>=SIZE_PROB_HISTO)|| (k< 0))  return 0.;
  else return Histo(k);
}

/****************************************************************************/

double CImaProb::repartition(double Val)
{
  int k = (int)  ( (Val-Tab_Bin(0)) / (Tab_Bin(1)-Tab_Bin(0)) );
  if (k>=SIZE_PROB_HISTO)  return 1.;
  else if (k<0)  return 0.;
  else return Repart(k);
}

/****************************************************************************/ 

void CImaProb::find_threshold(double Proba, double &ThresholdMin, double &ThresholdMax)
{
	int k=0;

	if(!symetry)
	{
		while ((k < SIZE_PROB_HISTO) && (Repart(k) < Proba)) k++;
		if (k >= SIZE_PROB_HISTO) k = SIZE_PROB_HISTO-1;
		// take the next bin, in case the current would contain many pixels 
		ThresholdMin = Tab_Bin(k)-StepHisto;
	}
	else ThresholdMin = 0.0;

	while ((k < SIZE_PROB_HISTO) && (Repart(k) < (1.-Proba))) k++;
	if (k >= SIZE_PROB_HISTO) k = SIZE_PROB_HISTO-1; 
	// take the next bin, in case the current would contain many pixels 
	ThresholdMax = Tab_Bin(k)+StepHisto;
}

/****************************************************************************/ 

void CImaProb::find_gthreshold(float N_Sigma, double &ThresholdMin, double &ThresholdMax)
{
	double Proba = (1. - erff((double) N_Sigma / sqrt((double) 2.)));

	if(!symetry) Proba/=2.0; // keep half on each side

	find_threshold(Proba, ThresholdMin, ThresholdMax);
}

/****************************************************************************/

void CImaProb::find_threshold(double Proba, double &ThresholdMax)
{
	int k=0;

	while ((k < SIZE_PROB_HISTO) && (Repart(k) < (1.-Proba))) k++;
	if (k >= SIZE_PROB_HISTO) k = SIZE_PROB_HISTO-1; 
	// take the next bin, in case the current would contain many pixels 
	ThresholdMax = Tab_Bin(k)+StepHisto;
}

/****************************************************************************/ 

void CImaProb::find_gthreshold(float N_Sigma, double &ThresholdMax)
{
	double Proba = (1. - erff((double) N_Sigma / sqrt((double) 2.)));
	
	find_threshold(Proba, ThresholdMax);
}

/****************************************************************************/

void CImaProb::write(char *Name)
{
  char FileName[255];

  sprintf(FileName, "%s_histo", Name);
//  fits_write_dblarr(FileName, Histo);
  sprintf(FileName, "%s_repart", Name);
//  fits_write_dblarr(FileName, Repart);
  sprintf(FileName, "%s_bin", Name);
//  fits_write_dblarr(FileName, Tab_Bin);
}

/****************************************************************************/

