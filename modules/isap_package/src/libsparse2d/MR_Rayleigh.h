
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
**    Date:  12.12.97
**    
**    File:  MR_Rayleigh.h
**
*******************************************************************************
**
**    DESCRIPTION  Class for  Rayleigh noise manadgement
**    ----------- 
**                 
******************************************************************************/

#ifndef _CRAYLEIGH_H_
#define _CRAYLEIGH_H_

#define SIZE_RAY_HISTO 1024 /* number of bins for histogram */


enum type_noise_speckle {NOISE_RAYLEIGH,
			 NOISE_LOG_RAYLEIGH,
			 NOISE_LAPLACE};

inline const char * StringTypeSpeckle (type_noise_speckle type)
{
    switch (type)
    {
    case NOISE_RAYLEIGH: 
      return ("Rayleigh noise");
      break;
    case NOISE_LOG_RAYLEIGH: 
      return ("Log transformed Rayleigh noise");
      break;
    case NOISE_LAPLACE: 
      return ("Square transformed Rayleigh noise");
      break;
    default:
      return ("Unknown noise");
    }
}

class StatRayleigh {
  type_noise_speckle StatNoise;
  fltarray StepHisto,MinHisto,MaxHisto;
  int NbrScale;
  int NbrBand;
  int NbrImage;
  fltarray Histo;
  fltarray Repart;
  fltarray Tab_Bin;
  type_transform Transform;
  void alloc(type_noise_speckle TypeNoise, int Nscale, 
             int Nima, type_transform Trans);
public:
  Bool Verbose;
  type_noise_speckle which_noise() {return StatNoise;}
  int nbr_scale() const {return NbrScale;}
  int nbr_band() const {return NbrBand;}
  type_transform type_trans() const {return Transform;}
  
  StatRayleigh(type_noise_speckle TypeNoise, int Nscale, 
               int Nima, type_transform Trans);
  StatRayleigh(type_noise_speckle TypeNoise, int Nscale, int Nima);
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
  
  ~StatRayleigh(){;}
};

#endif
