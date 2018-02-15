
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
**    File:  MR_CorrNoise.h
**
*******************************************************************************
**
**    DESCRIPTION  Class for  correlated noise manadgement
**    ----------- 
**                 
******************************************************************************/

#ifndef _CCORRNOISE_H_
#define _CCORRNOISE_H_

#define SIZE_CORR_HISTO 1024 /* number of bins for histogram */

class StatNoiseMap {
  fltarray StepHisto,MinHisto,MaxHisto;
  fltarray TabSigma;
  int NbrScale;
  int NbrBand;
  fltarray Histo;
  fltarray Repart;
  fltarray Tab_Bin;
  type_transform Transform;
public:
  Bool Verbose;
  int nbr_scale() const {return NbrScale;}
  int nbr_band() const {return NbrBand;}
  type_transform type_trans() const {return Transform;}
  float sigma_band(int Band)  const {return TabSigma(Band);}

  void alloc(Ifloat &NoiseMap, int Nscale, type_transform Trans, FilterAnaSynt *FAS=NULL, 
             sb_type_norm Norm=NORM_L1, int NbrUndec=-1, type_undec_filter U_Filter=DEF_UNDER_FILTER);

  StatNoiseMap(){Verbose=False;NbrBand=NbrScale=0;}

  StatNoiseMap(Ifloat &NoiseMap, int Nscale, type_transform Trans);
  StatNoiseMap(Ifloat &NoiseMap, int Nscale, type_transform Trans, FilterAnaSynt *FAS, sb_type_norm Norm,
               int NbrUndec, type_undec_filter U_Filter);
  
  void find_threshold(int NbrMaxBand, fltarray &Tab_Proba, fltarray &ThresholdMin, fltarray &ThresholdMax);
       // return the thresholds (min, max) corresponding 
       // to a Epsilon value for each scale
  
  float prob(int Band,  float Val); 
       // return the probability to have a wavelet
       // coefficient of value Val
   
  float repartition(int Band, float Val);
       // return the  integrated probability to have a wavelet coefficient
       // lower than Val  (Repartition Function)
    
  void write(); 
       // write histogram, repartition function and bin value 
  
  ~StatNoiseMap(){}
};

#endif
