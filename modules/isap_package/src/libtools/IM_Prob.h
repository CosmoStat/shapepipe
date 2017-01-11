
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
**    Date:  1.3.99
**    
**    File:  IM_Prob.h
**
*******************************************************************************
**
**    DESCRIPTION  Class for  correlated noise manadgement
**    ----------- 
**                 
******************************************************************************/

#ifndef _CIMAPROB_H_
#define _CIMAPROB_H_

#include "GlobalInc.h"

#define SIZE_PROB_HISTO 1024 /* number of bins for histogram */

class CImaProb {
  double StepHisto,MinHisto,MaxHisto;
  dblarray Histo;
  dblarray Repart;
  dblarray Tab_Bin;
  bool symetry;

public:
	Bool Verbose;  

	CImaProb();

	void set(Ifloat &NoiseMap);
	void set(float *NoiseIma, int N, bool sym=false);

	void find_threshold(double Proba, double & ThresholdMin,  double &ThresholdMax);
		// return the thresholds (min, max) corresponding 
		// to a Epsilon value 

	void find_gthreshold(float Nsigma, double & ThresholdMin,  double &ThresholdMax);
		// return the thresholds (min, max) corresponding 
		// to a Nsigma detection
	void find_threshold(double Proba, double &ThresholdMax);
	void find_gthreshold(float Nsigma, double &ThresholdMax);

	double prob(double Val); 
		// return the probability to have a  coefficient of value Val

	double  repartition(double Val);
		// return the  integrated probability to have a  coefficient
		// lower than Val  (Repartition Function)

	void write(char *Name); 
		// write histogram, repartition function and bin value 

	~CImaProb(){;}
};

#endif
