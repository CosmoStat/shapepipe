/******************************************************************************
**                   Copyright (C) 1998 by CEA
*******************************************************************************
**
**    UNIT
**
**    Version: 1.0
**
**    Author: Jean-Luc Starck
**
**    Date:  98/01/12 
**    
**    File:  Cerf.h
**
*******************************************************************************
**
**    DESCRIPTION: class MEM for signal restoration 
**    ----------- 
**                
******************************************************************************/

#ifndef _CMEM_H_
#define _CMEM_H_

const float DEF_MEM_STEP = 0.01;
const float DEF_MEM_MAX  = 5.;
const int   DEF_MEM_MAX_ITER = 100;
const float DEF_MEM_CVG = 0.001;

class CMemWave {
   double C1;    // internal constant = sqrt(2/PI)
   double C2;    // internal constant = sqrt(2)    
    
    
   double *TabHsGauss; // table for the derivative of Hs
   double *TabHnGauss; // table for the derivative of Hn
   int Np;       // size of the tables
      
   float Step;   // Step used for the probability integration
   float MaxVal; // Maximum value pre-calculated

   double grad_hn_sig1(double Val, CERF & ErTab);
   // calculate the derivative of hn in case of Gaussian noise (sigma=1)
   // hn = x/Sigma^2*erfc(x/(sqrt(2)sigma)) +
   //                        sqrt(2/PI)/sigma[1-exp(-x^2/(2sigma^2))]      

   double grad_hs_sig1(double Val, CERF & ErTab);
  // calculate the derivative of hs in case of Gaussian noise (sigma=1)
  // hs = -x/Sigma^2*erf(x/(sqrt(2)sigma)) +
  //                     sqrt(2/PI)/sigma[1-exp(-x^2/(2sigma^2))]   

   public:
   int MaxIter; // maximum number iteration for the dichotomy 
   float ConvgParam; // convergence parameter for the dichotomy
   
   CMemWave();  // initialization: set the table
 
   float grad_hs(float Val);  // calculate the derivative of hs at Val

   float grad_hn(float Val);  // calculate the derivative of hn at Val
 
    float filter (float CoefDat, float Alpha, float Sigma=1., float Model=0.);
   // calculate by the dichotomy method the solution wich minimize:
   // J = hs(y-x) + Alpha hn(x)
   // It is obtained when:  grad_hs(y-x) = Alpha grad_hn(x)
 
   float filter (float CoefDat, float Alpha, float SigmaS, float SigmaN, float Model);
   // calculate by the dichotomy method the solution wich minimize:
   // J = hs(y-x) + Alpha hn(x)
   // It is obtained when:  grad_hs(y-x) = Alpha grad_hn(x)
   
   void test_tab(fltarray &TabTest, float Nsig=3., Bool UseDataSNR=False);    
   // set a table TabTest(8,1000)  
   // TabTest(*,i) represents the filtering of a signal S, with 
   // S(k) = 0.01 * k : S is in [0,10[, with a step = 0.01
   // i = 0  ==> alpha = 0
   // i = 1  ==> alpha = 0.1
   // i = 2  ==> alpha = 0.5
   // i = 3  ==> alpha = 1
   // i = 4  ==> alpha = 2
   // i = 5  ==> alpha = 5
   // i = 6  ==> alpha = 10.;
   
   ~CMemWave();
};   

#endif
