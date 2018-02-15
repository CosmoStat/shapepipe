/*******************************************************************************
**
**    UNIT
**
**    Version: 1.0
**
**    Author: Jean-Luc Starck
**
**    Date:  5/03/2000 
**    
**    File:  Lifting.cc
**
**    Modification history:
**
*******************************************************************************
**
**    DESCRIPTION  1D lifting Sub-Band decomposition
**    -----------  
**                 
******************************************************************************/
 
#include "SB_Filter1D.h"

const double F79_Alpha = -1.586134342;
const double F79_Beta = -0.05298011854;
const double F79_Gamma = 0.8829110762;
const double F79_Delta = 0.4435068522;
const double F79_Xsi = 1.149604398;

/*************************************************************/
// LIFTING SCHEME
/*************************************************************/

const char * StringLSTransform (type_lift  type)
{
    switch (type)
    {
        case TL_MEDIAN:
	      return ("Lifting scheme: median prediction"); 
        case  TL_INT_HAAR: 
              return ("Lifting scheme: integer Haar WT"); 
        case  TL_INT_CDF: 
              return ("Lifting scheme: integer CDF WT"); 
        case  TL_INT_4_2:
	      return ("Lifting scheme: integer (4,2) interpolating transform");
	case  TL_CDF: 
              return ("Lifting scheme: CDF WT"); 
        case  TL_F79: 
              return ("Lifting scheme: 7/9 WT"); 
	case  TL_INT_F79: 
              return ("Lifting scheme: integer 7/9 WT"); 	      
        default: 
	      return ("Undefined transform");
     }
}

/*************************************************************/

void lifting_usage(type_lift Transform)
{
    fprintf(OUTMAN, "        [-l type_of_lifting_transform]\n");
    for (int i = 1; i <= NBR_LIFT; i++)
    fprintf(OUTMAN, "              %d: %s \n",i,
                                           StringLSTransform  ((type_lift )i));
    fprintf(OUTMAN, "             default is %s\n", StringLSTransform((type_lift) Transform));
}

/*************************************************************/

float Lifting::lift_predict(int i, int N, float *High, int Step)
{
   float Predict=0.;
   int Ind1,Ind2,Ind3,Ind4;
   int IndHigh = test_index(i-Step, N);
   int Index = test_index(i+Step, N);
   float TabVal[5];
   
   if ((IndHigh < 0) || (IndHigh >= N))
   {
      cout << "PB IndHigh: " << IndHigh << " N = " << N << " i = " << i << " Step = " << Step <<endl;
      exit(-1);
   }
   if ((Index < 0) || (Index >= N))
   {
      cout << "PB IndHigh: " << IndHigh << " N = " << N << " i = " << i << " Step = " << Step <<endl;
      exit(-1);
   }
    switch(TypeTrans)
    {
        case TL_MEDIAN:
           Ind1 = test_index(IndHigh+2*Step, N); 
           Ind2 = test_index(IndHigh-2*Step, N);      
           Ind3 = test_index(IndHigh+4*Step, N); 
           Ind4 = test_index(IndHigh-4*Step, N);
           TabVal[0] = High[IndHigh];
           TabVal[1] = High[Ind1];
           TabVal[2] = High[Ind2];
	   TabVal[3] = High[Ind3];
	   TabVal[4] = High[Ind4];
	   Predict = opt_med5(TabVal);
 	   break;
        case TL_INT_HAAR:
	    Predict = (int) (High[IndHigh] + 0.5);
	   break;
	case TL_INT_CDF: 
	   Predict = (int)(0.5*(High[IndHigh]+High[Index]) + 0.5);
           break;
	case TL_CDF: 
	   Predict = 0.5*(High[IndHigh]+High[Index]);
	   break;
	case TL_INT_4_2:
	   Ind1 = test_index(IndHigh+2*Step, N); 
           Ind2 = test_index(IndHigh-2*Step, N);      
	   Ind3 = test_index(IndHigh+4*Step, N); 
	   Predict = (int)(9./16.*(High[IndHigh]+High[Ind1]) -
		               (1./16.*(High[Ind2]+High[Ind3])) + 0.5);
	break;
	default:
	        cerr << "Error: unknown lifting transform ... " << endl;
	        exit(-1);
	        break;
   }
   return Predict;
}

/*************************************************************/

float Lifting::lift_update(int i, int N, float *Det, int Step)
{
   float Update=0.;
   int Ind1, IndDet = test_index(i, N);
   int Index = test_index(IndDet - Step, N);
   switch(TypeTrans)
   {
       case TL_MEDIAN:
	     Update = Det[IndDet];
             break;
       case TL_INT_HAAR:
	     Update =  (int) (Det[IndDet]/2. + 0.5);
             break;
       case TL_INT_CDF: 
	     Update = (int)(0.25*( Det[IndDet]+Det[Index]) + 0.5);
	     break;
       case TL_CDF: 
	     Update = 0.25*( Det[IndDet]+Det[Index]);
             break;
       case TL_INT_4_2: 
	     Ind1 = test_index(IndDet-Step, N);
  	     Update =  (int)(1./4*(Det[IndDet]+Det[Ind1]) + 0.5);
             break;
       default: 
	     cerr << "Error: unknown lifting transform ... " << endl;
	     exit(-1);
	     break;
    }
    return  Update;
}
   
/*************************************************************/

void Lifting::transform(int N, float *High, float *Low, float *Det)
{
   int i;
      // linear prediction
   if  (TypeTrans == TL_F79) transform_f79(N, High, Low,Det);
   else if (TypeTrans == TL_INT_F79) transform_int_f79(N, High, Low,Det);
   else
   {
      for (i = 1; i < N; i += 2)
        Det[i/2] =  High[i] - lift_predict(i, N, High, DistPix);
      for (i = 0; i < N; i += 2)
        Low[i/2] =  High[i] + lift_update(i/2, N/2, Det, DistPix);
   }
}

/*************************************************************/

void Lifting::transform_f79(int N, float *High, float *Low, float *Det)
{
      int i;
      int Nd = N/2;
      int Ns = (N+1)/2;
       
      // Split
      for (i = 0; i < N; i += 2) Low[i/2] = High[i];
      for (i = 1; i < N; i += 2) Det[i/2] = High[i];
      
      for (i = 0; i < Nd; i ++)  
          Det[i] += F79_Alpha*(Low[i]+Low[test_index(i+1,Ns)]);
   
      for (i = 0; i < Ns; i ++)  
          Low[i] += F79_Beta*(Det[test_index(i,Nd)]+Det[test_index(i-1,Nd)]);

      for (i = 0; i < Nd; i ++)  
          Det[i] += F79_Gamma*(Low[i]+Low[test_index(i+1,Ns)]);
   
      for (i = 0; i < Ns; i ++)  
          Low[i] += F79_Delta*(Det[test_index(i,Nd)]+Det[test_index(i-1,Nd)]);
   
      for (i = 0; i < Ns; i ++)  Low[i] *= F79_Xsi;     
      for (i = 0; i < Nd; i ++)  Det[i] /= F79_Xsi;
}

/*************************************************************/

void Lifting::transform_int_f79(int N, float *High, float *Low, float *Det)
{
      int i;
      int Nd = N/2;
      int Ns = (N+1)/2;
       
      // Split
      for (i = 0; i < N; i += 2) Low[i/2] = High[i];
      for (i = 1; i < N; i += 2) Det[i/2] = High[i];
      
      for (i = 0; i < Nd; i ++)  
          Det[i] += (int)(F79_Alpha*(Low[i]+Low[test_index(i+1,Ns)])+0.5);
   
      for (i = 0; i < Ns; i ++)  
          Low[i] += (int)(F79_Beta*(Det[test_index(i,Nd)]+Det[test_index(i-1,Nd)])+0.5);

      for (i = 0; i < Nd; i ++)  
          Det[i] += (int)(F79_Gamma*(Low[i]+Low[test_index(i+1,Ns)])+0.5);
   
      for (i = 0; i < Ns; i ++)  
          Low[i] += (int)(F79_Delta*(Det[test_index(i,Nd)]+Det[test_index(i-1,Nd)])+0.5);   
}

/*************************************************************/

void Lifting::recons_int_f79(int N, float *Low, float *Det, float *High)
{
    int i;    
    int Nd = N/2;
    int Ns = (N+1)/2;

    for (i = 0; i < Ns; i ++)  
          Low[i] -= (int) (F79_Delta*(Det[test_index(i,Nd)]+Det[test_index(i-1,Nd)])+0.5);

    for (i = 0; i < Nd; i ++)  
          Det[i] -= (int) (F79_Gamma*(Low[i]+Low[test_index(i+1,Ns)])+0.5);
	  
    for (i = 0; i < Ns; i ++)  
          Low[i] -= (int) (F79_Beta*(Det[test_index(i,Nd)]+Det[test_index(i-1,Nd)])+0.5);	  
	  
    for (i = 0; i < Nd; i ++)  
          Det[i] -= (int) (F79_Alpha*(Low[i]+Low[test_index(i+1,Ns)])+0.5);
	  	  
    for (i = 0; i < N; i += 2) High[i] = Low[i/2];
    for (i = 1; i < N; i += 2) High[i] = Det[i/2];
    
}

/*************************************************************/

void Lifting::recons_f79(int N, float *Low, float *Det, float *High)
{
    int i;    
    int Nd = N/2;
    int Ns = (N+1)/2;

    for (i = 0; i < Ns; i ++)  Low[i] /= F79_Xsi;     
    for (i = 0; i < Nd; i ++)  Det[i] *= F79_Xsi;

    for (i = 0; i < Ns; i ++)  
          Low[i] -= F79_Delta*(Det[test_index(i,Nd)]+Det[test_index(i-1,Nd)]);

    for (i = 0; i < Nd; i ++)  
          Det[i] -= F79_Gamma*(Low[i]+Low[test_index(i+1,Ns)]);
	  
    for (i = 0; i < Ns; i ++)  
          Low[i] -= F79_Beta*(Det[test_index(i,Nd)]+Det[test_index(i-1,Nd)]);	  
	  
    for (i = 0; i < Nd; i ++)  
          Det[i] -= F79_Alpha*(Low[i]+Low[test_index(i+1,Ns)]);
	  	  
    for (i = 0; i < N; i += 2) High[i] = Low[i/2];
    for (i = 1; i < N; i += 2) High[i] = Det[i/2];
    
}

/*************************************************************/

void Lifting::recons(int N, float *Low, float *Det, float *High)
{
    int i;    
   if  (TypeTrans == TL_F79) recons_f79(N, Low,Det,High);
   else if (TypeTrans == TL_INT_F79) recons_int_f79(N,Low,Det,High);
   else
   {
    for (i = 0; i < N; i += 2)
           High[i] = Low[i/2] - lift_update(i/2,N/2,Det, DistPix);
    for (i = 1; i < N; i += 2)
           High[i] = Det[i/2] + lift_predict(i, N, High, DistPix);
   }
}

/*************************************************************/

void Lifting::transform(int N, float *High, float *Low, float *Det, int Step)
{
   int i;
   if  (TypeTrans == TL_F79) 
   {
     cout << "not implemented ... " << endl;
     exit (-1);
     // transform_f79(N, High, Low,Det);
   }
   else  if (TypeTrans == TL_INT_F79) 
   {
     cout << "not implemented ... " << endl;
     exit (-1);
      // transform_int_f79(N, High, Low,Det);
   }
   else
   {
      for (i = 0; i < N; i ++) 
            Det[i] =  High[i] - lift_predict(i, N, High, Step);
      for (i = 0; i < N; i ++) 
          Low[i] =  High[i] + lift_update(i+Step, N, Det, 2*Step);
   }
}

/*************************************************************/

void Lifting::recons(int N, float *Low, float *Det, float *High, int Step)
{
    int i;    
    float *Temp = new float [N];
    
    for (i = 0; i < N; i ++)  High[i] = Temp[i] = 0;
    
    for (i = 0; i < N; i ++)
           High[i] = Low[i] - lift_update(i+Step, N, Det, 2*Step);
    for (i = 1; i < N; i ++)
            High[i] = Det[i] + lift_predict(i, N, High, Step);

    for (i = 1; i < N; i ++)
           Temp[i] = Low[i] - lift_update(i+Step, N, Det, 2*Step);
    for (i = 0; i < N; i ++)
           Temp[i] = Det[i] + lift_predict(i, N, Temp, Step);

    for (i = 0; i < N; i++) High[i] =  0.5 *(Temp[i]+High[i]);
	   
    delete [] Temp;
}

/*************************************************************/
