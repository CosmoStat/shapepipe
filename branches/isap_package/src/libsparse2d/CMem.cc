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
**    File:  CMem.cc
**
*******************************************************************************
**
**    DESCRIPTION  class MEM for signal restoration
**    ----------- 
**                 
******************************************************************************/
 
#include "GlobalInc.h"
#include "CErf.h"
#include "CMem.h"

/****************************************************************************/
 
// calculate the derivative of hn in case of Gaussian noise
// hn = x/Sigma^2*erfc(x/(sqrt(2)sigma)) +
//                        sqrt(2/PI)/sigma[1-exp(-x^2/(2sigma^2))]
double CMemWave::grad_hn_sig1(double Val, CERF & ErTab)
{      
    int Sign = (Val >= 0) ? 1: -1;
    double Coef = Val*Sign;
    Coef = Coef*ErTab.get_erfc(Coef/C2) + C1 * (1-exp(-(Coef*Coef)/2.));
    return  Coef*Sign;
}
 
/****************************************************************************/
   
// calculate the derivative of hs in case of Gaussian noise
// hs = -x/Sigma^2*erf(x/(sqrt(2)sigma)) +
//                        sqrt(2/PI)/sigma[1-exp(-x^2/(2sigma^2))]
double CMemWave::grad_hs_sig1(double Val, CERF & ErTab)
{
     int Sign = (Val >= 0) ? 1: -1;
     double Coef = Val*Sign;
     Coef =  -Coef*ErTab.get_erf(Coef/C2) + C1 * (1-exp(-(Coef*Coef)/2.));
     return  Coef*Sign;
}
   
/***************************************/

   // initialization
CMemWave::CMemWave()
{       
   double Val=0.;
   int i=0;
   CERF ErTab; 

   MaxIter= DEF_MEM_MAX_ITER;
   ConvgParam= DEF_MEM_CVG;
   C1 = sqrt(2./PI);
   C2 = sqrt(2.);
   Step = DEF_MEM_STEP;
   MaxVal = DEF_MEM_MAX;
   
   Np = (int) (MaxVal / Step + 1.5);
   TabHsGauss = new double [Np];
   TabHnGauss = new double [Np];
      
   while (Val < MaxVal)
   {
      TabHsGauss[i] =  grad_hs_sig1(Val, ErTab);
      TabHnGauss[i++] = grad_hn_sig1(Val, ErTab);
      Val += Step;
   }         
}    
  
/***************************************/

float CMemWave::grad_hs(float Val)
{
    int Ind;
    double ValRet = 0.;
        
    if (Val >= MaxVal) ValRet = TabHsGauss[Np-1] - Val;
    else if (Val > 0)
    {
       Ind = (int) (Val / Step + 0.5);
       ValRet = TabHsGauss[Ind];
    }
    return ValRet;
}

/***************************************/

float CMemWave::grad_hn(float Val)
{
   int Ind;
   double ValRet = 0.;
        
   if (Val >= MaxVal) ValRet = TabHnGauss[Np-1];
   else if (Val > 0)
   {
       Ind = (int) (Val / Step + 0.5);
       ValRet = TabHnGauss[Ind];
   }
   return ValRet;
}
    
/***************************************/

// calculate by dichotomy method the solution wich minimize:
// J = hs(y-x) + Alpha hn(x)
// It is obtained when:  grad_hs(y-x) = Alpha grad_hn(x)
float CMemWave::filter (float CoefDat, float Alpha, float Sig, float Model)
{
   float ValRet=0;
   float Sigma=Sig;
   if (Sigma < FLOAT_EPSILON) Sigma = FLOAT_EPSILON ;
   
   if (ABS(Alpha) < FLOAT_EPSILON) ValRet = CoefDat;
   else if (Alpha < 0) ValRet = Model;
   else 
   {
      double ValDat = ABS(CoefDat-Model) / Sigma;
      double CoefMax = ValDat+3;
      double Coef,Diff,CoefMin = 0;
      int Iter=0;
      do
      {
        Iter++;
        Coef = (CoefMin+CoefMax)/2.;
        Diff =  grad_hs(ValDat-Coef) +  Alpha * grad_hn(Coef);
        if (Diff > 0.) CoefMax = Coef;
        else CoefMin = Coef;
       } while (( CoefMax-CoefMin > ConvgParam) && (Iter < MaxIter));
       Coef = (CoefDat >= Model) ? Coef : -Coef; 
       ValRet = Coef*Sigma+Model;
   }   
   return ValRet;
 }    
    
/***************************************/


// calculate by dichotomy method the solution wich minimize:
// J = hs(y-x) + Alpha hn(x)
// It is obtained when:  grad_hs(y-x) = Alpha grad_hn(x)
float CMemWave::filter (float CoefDat, float Alpha, float SigS, float SigN, float Model)
{
   float ValRet=0;
   float Sigma=SigS;
   float SigmaN=SigN;
    
   if (Sigma < FLOAT_EPSILON) Sigma = FLOAT_EPSILON ;
   if (SigmaN < FLOAT_EPSILON) SigmaN = FLOAT_EPSILON ;

   if (ABS(Alpha) < FLOAT_EPSILON) ValRet = CoefDat;
   else if (Alpha < 0) ValRet = Model;
   else 
   {
      double ValDat = ABS(CoefDat-Model);
      double CoefMax = ValDat+3*MAX(Sigma,SigmaN);
      double Coef,Diff,CoefMin = 0;
      int Iter=0;
      do
      {
        Iter++;
        Coef = (CoefMin+CoefMax)/2.;
        Diff =  Sigma*grad_hs( (ValDat-Coef)/Sigma) +  Alpha * SigmaN * grad_hn(Coef/SigmaN);
        if (Diff > 0.) CoefMax = Coef;
        else CoefMin = Coef;
       } while (( (CoefMax-CoefMin)/Sigma > ConvgParam) && (Iter < MaxIter));
       Coef = (CoefDat >= Model) ? Coef : -Coef; 
       ValRet = Coef+Model;
   }   
   return ValRet;
 }    
   /***************************************/

void CMemWave::test_tab(fltarray &TabTest, float NSig, Bool UseDataSNR)
{
   int k,U=1000;
   fltarray Tab(U,8);
   float AlphaP, Rg=1., RegulParam; 
   float TabReg[7];
   TabReg[0] = 0.;
   TabReg[1] = 0.1;
   TabReg[2] = 0.5;
   TabReg[3] = 1.;
   TabReg[4] = 2.;
   TabReg[5] = 5.;
   TabReg[6] = 10.;
   
   for (int i=0; i< U; i++) 
   {
       Tab(i,0) = 0.01*i;
       if (UseDataSNR == True)
       {
          AlphaP = ABS( Tab(i,0) / NSig);
          if (AlphaP > 1) AlphaP = 1;
          if (AlphaP < FLOAT_EPSILON) Rg = -1;
          else Rg = (1-AlphaP)/AlphaP;
      }
      for (k=1; k < 7; k++)
      {
         RegulParam = TabReg[k]*Rg;
 	 if (ABS(RegulParam) < FLOAT_EPSILON)  Tab(i,k)=  Tab(i,0);
         else if (RegulParam < 0)   Tab(i,k) = 0;
         else  Tab(i, k) = filter( Tab(i,0), RegulParam);   
      }
           
       // cout << "Val = " << i*0.5 << " Dicho = " << filter(0.01*i, RegulParam) << endl;
    }
    TabTest = Tab;
 }

/***************************************/

CMemWave::~CMemWave()  
{
   if (TabHnGauss  != NULL) delete  [] TabHnGauss;
   if ( TabHsGauss != NULL) delete  [] TabHsGauss;
   Np = 0;
}

/***************************************/


