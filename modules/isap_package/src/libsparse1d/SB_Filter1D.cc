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
**    File:  SB_Filter1D.cc
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
 
#define FIL_DEBUG 0


static void filter_complement(float *Filter, float *ComFilter, int N, int Signe)
{
    int i,Sgn=Signe;
    for (i = 0; i < N; i++)
    {
        ComFilter[N-i-1] = Sgn * Filter[i];
        Sgn *= -1;
    }       
}

// static void filter_transpose(float *Filter, int N)
// {
//     int i;
//     float *Aux = new float[N];
//     for (i = 0; i < N; i++) Aux[i] = Filter[i];
//     for (i = 0; i < N; i++) Filter[N-i-1] = Aux[i];
//     delete Aux;
// }

/***********************************************************************/

void SubBandFilter::reset_param()
{
   H0 = G0 = H1 = G1 = NULL;
   Size_H0 = Size_H1 = Size_G0 = Size_G1 = 0;
   Start_H0 = Start_H1 = Start_G0 = Start_G1 =0;   
   TypeFilter = SB_UNKNOWN;
   NormCoef = 0;
   Verbose = False;

}

/***********************************************************************/

SubBandFilter::SubBandFilter(FilterAnaSynt & FILTER)
{
   reset_param();
   init(FILTER, DEF_SB_NORM);
}

/***********************************************************************/

SubBandFilter::SubBandFilter(FilterAnaSynt & FILTER, sb_type_norm Norm)
{
   reset_param();
   init(FILTER, Norm);
}

/***********************************************************************/

SubBandFilter::SubBandFilter(type_sb_filter T_Filter)
{
   reset_param();
   FilterAnaSynt FILTER(T_Filter);
   init(FILTER, DEF_SB_NORM);
}

/***********************************************************************/

SubBandFilter::SubBandFilter(type_sb_filter T_Filter, sb_type_norm Norm)
{
   reset_param();
   FilterAnaSynt FILTER(T_Filter);
   init(FILTER, Norm);
}

/***********************************************************************/

SubBandFilter::SubBandFilter(char *FileName)
{
   reset_param();
   FilterAnaSynt FILTER(FileName);
   init(FILTER, DEF_SB_NORM);
}

/***********************************************************************/

SubBandFilter::SubBandFilter(char *FileName, sb_type_norm Norm)
{
   reset_param();
   FilterAnaSynt FILTER(FileName);
   init(FILTER, Norm);
}


/***********************************************************************/
 
SubBandFilter::SubBandFilter(float *H, int Nh, float *G, int Ng, sb_type_norm FilNorm)
{
   sb_type_norm FilterNorm = NORM_L1; // filters h and g are supposed to be L1 normalized
   int i,Sgn,FirstPix;
   TypeNorm = FilNorm;
   NormCoef = (TypeNorm == NORM_L1) ? 2 : 1;
      
   Size_H0 = Nh;
   Size_G0 = Ng;
   Size_H1 = Size_G0;
   Size_G1 = Size_H0;
   H0 = new float[Size_H0];
   H1 = new float[Size_H1];
   G1 = new float[Size_G1];
   G0 = new float[Size_G0];
   float NormCorrection=1;
   if (FilterNorm == TypeNorm) NormCorrection = 1.;
   else if (TypeNorm == NORM_L1)  NormCorrection = 1. / sqrt(2.);
   else if (TypeNorm == NORM_L2)  NormCorrection = sqrt(2.);
   for (i=0; i < Size_H0; i++) H0[i] = H[i]*NormCorrection;
   for (i=0; i < Size_G0; i++) G0[i] = G[i]*NormCorrection;
   // for (i=0; i < Size_G0; i++) H1[i] = G[i]*NormCorrection;
    
   Start_H0 = - Size_H0 / 2;
   Start_H1 = - Size_H1 / 2;
   Start_G0 = - Size_G0 / 2;
   Start_G1 = - Size_G1 / 2; 

    FirstPix = - Size_H0 / 2;
    if (FirstPix % 2 == 0) Sgn = 1;
    else Sgn = -1;    
    filter_complement(H0, G1, Size_H0, Sgn);
    FirstPix = - Size_H1 / 2;
    if (FirstPix % 2 == 0) Sgn = 1;
    else Sgn = -1;
    filter_complement(G0, H1, Size_G0, Sgn);  
    // filter_complement(H1, G0, Size_H1, Sgn);
    float T=0;
    for (i=0; i < Size_H0; i++) T += H0[i];
    printf("H0 = %f\n", T);
    for (i=0; i < Size_H0; i++) printf("%f\n", H0[i]);

    T=0;
    for (i=0; i < Size_G0; i++) T+= G0[i];
    printf("G0 = %f\n", T);
    for (i=0; i < Size_G0; i++) printf("%f\n", G0[i]);
    
    T=0;
    for (i=0; i < Size_H1; i++) T+= H1[i];
    printf("H1 = %f\n", T);
    for (i=0; i < Size_H1; i++) printf("%f\n", H1[i]);
    
    T=0;
    for (i=0; i < Size_G1; i++) T+= G1[i];
    printf("G1 = %f\n", T);
    for (i=0; i < Size_G1; i++) printf("%f\n", G1[i]);
}

/***********************************************************************/

void SubBandFilter::init(FilterAnaSynt & FILTER, sb_type_norm Norm)
{
   int i,Sgn,FirstPix;
   TypeNorm = Norm;
   NormCoef = (TypeNorm == NORM_L1) ? 2 : 1;
  
   TypeFilter = FILTER.type_filter();
   float *F_H0 = FILTER.analysis();
   float *F_H1 = FILTER.synthesis();
   Size_H0 = FILTER.size_analysis();
   Size_H1 = FILTER.size_synthesis();
   Size_G1 = Size_H0;
   Size_G0 = Size_H1;
   H0 = new float[Size_H0];
   H1 = new float[Size_H1];

   // If the wanted normalization is not this of the filter
   // we need to correct the coefficient
   float NormCorrection=1;
   if (FILTER.norm () == TypeNorm) NormCorrection = 1.;
   else if (TypeNorm == NORM_L1)  NormCorrection = 1. / sqrt(2.);
   else if (TypeNorm == NORM_L2)  NormCorrection = sqrt(2.);

   for (i=0; i < Size_H0; i++) H0[i] = F_H0[i]*NormCorrection;
   for (i=0; i < Size_H1; i++) H1[i] = F_H1[i]*NormCorrection;
 
   // float SumH=0.,SumG=0.;
   // for (i=0; i < Size_H0; i++) SumH += H0[i]*H0[i];
   // for (i=0; i < Size_H1; i++) SumG += H1[i]*H1[i];
   // cout << " BAND SUM = " << Size_H0 <<  " " << Size_H1 <<  " " << SumH << " " << SumG << endl;
   // cout << " H0 " << endl;
   // for (i=0; i < Size_H0; i++) printf("%f\n", H0[i]);
   // cout << " H1 " << endl;
   // for (i=0; i < Size_H1; i++) printf("%f\n", H1[i]);
   
   G1 = new float[Size_G1];
   G0 = new float[Size_G0];

    FirstPix = - Size_H0 / 2;
    if (FirstPix % 2 == 0) Sgn = 1;
    else Sgn = -1;    
    filter_complement(H0, G1, Size_H0, Sgn);
    FirstPix = - Size_H1 / 2;
    if (FirstPix % 2 == 0) Sgn = 1;
    else Sgn = -1;
    filter_complement(H1, G0, Size_H1, Sgn);    
    Start_H0 = - Size_H0 / 2;
    Start_H1 = - Size_H1 / 2;
    Start_G0 = - Size_G0 / 2;
    Start_G1 = - Size_G1 / 2; 

#if FIL_DEBUG
    printf("Filter \n");
    for (i = 0; i <  Size_H0; i++)
                printf("H0(%d) = %1.11f\n", i,  H0[i]);
    printf("\n");
    for (i = 0; i <  Size_G0; i++)
                printf("G0(%d) = %1.11f\n", i, G0[i]);
    printf("\n");
    for (i = 0; i <   Size_H1; i++)
              printf("H1(%d) = %1.11f\n", i,  H1[i]);
    printf("\n");
    for (i = 0; i <  Size_G1; i++)
              printf("G1(%d) = %1.11f\n", i,  G1[i]);
    printf("ENDFilter \n");
#endif
};

/*************************************************************/

SubBandFilter::~SubBandFilter() 
{
     if (H0 != NULL) delete [] H0;
     if (H1 != NULL) delete [] H1;
     if (G0 != NULL) delete [] G0;
     if (G1 != NULL) delete [] G1;
     TypeFilter = SB_UNKNOWN;
};
 
/*************************************************************/

void SubBandFilter::convol_h0 (int N, float *Input, float *Output )
/* convolves the data with the filter h0 */
{
    int i, Index, p;
    int SizeFilter = Size_H0;    
    float *F = H0;
    int Dep = 0;
//      static int Pas=0;
//     if (Pas == 0)
//     {
//        if (SubSample_H_Even  == False) cout << "SubSample_H_Even False" << endl;
//        else cout << "SubSample_H_Even True" << endl;
//        Pas=1; 
//     }   
    if (SubSample_H_Even == False) Dep = 1;
    if(Border== I_ZERO) { //FCS ADDED
       for (i = Dep; i < N; i += 2){
           Output[i/2] = 0;
           for (p = 0; p < SizeFilter; p++) {
              Index = i + Start_H0 + p;
              if((Index >= 0) && (Index < N)) Output[i/2] += Input[Index] * F[SizeFilter-p-1];
	   } 
       }
       if ((Dep == 1) && (N % 2 == 1)){
          i = N-1;
          Output[N/2] = 0;
          for (p  = 0; p < SizeFilter; p ++){
              Index = i + p + Start_H0 ;
              if((Index >= 0) && (Index < N)) Output[N/2] +=  Input[Index]*F[SizeFilter-p-1];
          }
      }
    } else {
        for (i = Dep; i < N; i += 2)
        {
           Output[i/2] = 0;
           for (p = 0; p < SizeFilter; p++)
	   {
              Index = i + Start_H0 + p;
              Index = test_index(Index, N);
              Output[i/2] += Input[Index] * F[SizeFilter-p-1];
	   } 
       }
       if ((Dep == 1) && (N % 2 == 1))
       {
          i = N-1;
          Output[N/2] = 0;
          for (p  = 0; p < SizeFilter; p ++)
          {
              Index = i + p + Start_H0 ;
              Index = test_index(Index, N);
              Output[N/2] +=  Input[Index]*F[SizeFilter-p-1];
          }
      }
   }
}

/*************************************************************/

void SubBandFilter::convol_g0 (int N, float *Input, float *Output)
{
    int i, Index, p;
    int SizeFilter = Size_G0;
    float *F = G0;
    int Dep = 1;
//     static int Pas=0;
//     if (Pas == 0)
//     {
//         if (SubSample_G_Odd == False) cout << "SubSample_G_Odd False" << endl;
//        else cout << "SubSample_G_Odd True" << endl;
//        Pas=1;
//     }
    
    if (SubSample_G_Odd == False) Dep = 0;
    if(Border== I_ZERO) { //FCS ADDED

       for (i = Dep; i < N; i += 2){        
           Output[i/2] = 0;
           for (p  = 0; p < SizeFilter; p ++){
              Index = i + Start_G0 + p;
              if((Index >= 0) && (Index < N)) Output[i/2] +=  Input[Index]*F[SizeFilter-p-1];
	   }
       }
       if (N % 2 == 1){
          i = N-1;
          Output[N/2] = 0;
          for (p  = 0; p < SizeFilter; p ++){
              Index = i + p + Start_G0;
              if((Index >= 0) && (Index < N)) Output[N/2] +=  Input[Index]*F[SizeFilter-p-1];
          }
       }

    } else {
       for (i = Dep; i < N; i += 2)
       {        
           Output[i/2] = 0;
           for (p  = 0; p < SizeFilter; p ++)
	   {
              Index = i + Start_G0 + p;
              Index = test_index(Index, N);
              Output[i/2] +=  Input[Index]*F[SizeFilter-p-1];
	   }
       }
       if (N % 2 == 1)
       {
          i = N-1;
          Output[N/2] = 0;
          for (p  = 0; p < SizeFilter; p ++)
          {
              Index = i + p + Start_G0;
              Index = test_index(Index, N);
              Output[N/2] +=  Input[Index]*F[SizeFilter-p-1];
          }
      }
   }
}

/*************************************************************/

void SubBandFilter::noise_convol_h0 (int N, float *Input, float *Output )
/* convolves the data with the filter h0 */
{
    int i, Index, p;
    int SizeFilter = Size_H0;    
    float *F = H0;
    if(Border== I_ZERO) { //FCS ADDED
       for (i = 0; i < N; i += 2)
       {
           Output[i/2] = 0;
           for (p = 0; p < SizeFilter; p++)
           {
              Index = i + (p + Start_H0)*DistPix;
              if((Index >= 0) && (Index < N)) Output[i/2] += Input[Index] * pow((double) F[SizeFilter-p-1], 2.0);
	   }  
       }

    } else {
       for (i = 0; i < N; i += 2)
       {
           Output[i/2] = 0;
           for (p = 0; p < SizeFilter; p++)
	   {
              Index = i + (p + Start_H0)*DistPix;
              Index = test_index(Index, N);
              Output[i/2] += Input[Index] * pow((double) F[SizeFilter-p-1], 2.0);
	   } 
       }
   }
}

/*************************************************************/

void SubBandFilter::noise_convol_g0 (int N, float *Input, float *Output)
{
    int i, Index, p;
    int SizeFilter = Size_G0;
    float *F = G0;
    
    if(Border== I_ZERO) { //FCS ADDED
       for (i = 1; i < N; i += 2) {        
           Output[i/2] = 0;
           for (p  = 0; p < SizeFilter; p ++)
	   {
              Index = i+ (p + Start_G0)*DistPix;
              if((Index >= 0) && (Index < N)) Output[i/2] +=  Input[Index]*pow((double) F[SizeFilter-p-1],2.0);
	   }
       }
       if (N % 2 == 1)
       {
          i = N-1;
          Output[N/2] = 0;
          for (p  = 0; p < SizeFilter; p ++)
          {
              Index = i+ (p + Start_G0)*DistPix;
              if((Index >= 0) && (Index < N)) Output[N/2] +=  Input[Index]*pow((double) F[SizeFilter-p-1],2.0);
          }
      }
    } else {
       for (i = 1; i < N; i += 2)
       {        
          Output[i/2] = 0;
           for (p  = 0; p < SizeFilter; p ++)
	   {
              Index = i+ (p + Start_G0)*DistPix;
              Index = test_index(Index, N);
              Output[i/2] +=  Input[Index]*pow((double) F[SizeFilter-p-1],2.0);
	   }
       }
       if (N % 2 == 1)
       {
          i = N-1;
          Output[N/2] = 0;
          for (p  = 0; p < SizeFilter; p ++)
          {
              Index = i+ (p + Start_G0)*DistPix;
              Index = test_index(Index, N);
              Output[N/2] +=  Input[Index]*pow((double) F[SizeFilter-p-1],2.0);
          }
      }
   }
}

/*************************************************************/
 
void SubBandFilter::convol_h1 (int N, float *Input, float *Output)
/* convolves the data with the filter h1 */
{
    int i, Index, p;
    float   *Temp = new float [N];
    int SizeFilter = Size_H1;
    float *F = H1;
    int Dep = 0;
    for (i=0; i < N; i++) Temp[i]=0;

//      static int Pas=0;
//      if (Pas == 0)
//      {
//        if (SubSample_H_Even  == False) cout << "SubSample_H_Even False" << endl;
//        else cout << "SubSample_H_Even True" << endl;
//        Pas=1;    
//     }
    if (SubSample_H_Even == False) Dep = 1;
    for (i = Dep; i < N; i += 2)  Temp[i] = Input[i/2];
    if(Border== I_ZERO) { //FCS ADDED
       for (i = 0; i < N; i++) {
           Output[i] = 0;
           for (p = 0; p < SizeFilter; p++){
               Index = i + Start_H1 + p;
               if((Index >= 0) && (Index < N)) Output[i] += Temp[Index]*F[p];
	   }
       }
   
    } else {
       for (i = 0; i < N; i++)
       {
           Output[i] = 0;
           for (p = 0; p < SizeFilter; p++)
	   {
               Index = i + Start_H1 + p;
               Index = test_index(Index, N);
               Output[i] += Temp[Index]*F[p];
	   }
       }
    }
    delete [] Temp;
}

/****************************************************************************/

void SubBandFilter::convol_g1 (int N, float *Input, float *Output)
/* convolves the data with the filter g1 */
{
    int i,Index, p;
    float *Temp = new float [N];
    int SizeFilter = Size_G1;
    float *F = G1;    
    int Dep = 1;
//     static int Pas=0;
//     if (Pas == 0)
//     {
//       if (SubSample_G_Odd == False) cout << "SubSample_G_Odd False" << endl;
//       else cout << "SubSample_G_Odd True" << endl;
//       Pas=1;
//     }    
    for (i=0; i < N; i++) Temp[i]=0;
 
    if (SubSample_G_Odd == False) Dep = 0;
    for (i = Dep; i<N; i += 2)  Temp[i] = Input[i/2];
    if(Border== I_ZERO) { //FCS ADDED
        for (i = 0; i< N; i++) { 
            Output[i] = 0;
            for (p = 0; p < SizeFilter; p++) {
               Index = i + Start_G1 + p;
               if((Index >= 0) && (Index < N)) Output[i] += Temp[Index]*F[p];
	    }
        }

    } else {
       for (i = 0; i< N; i++)
       {  
           Output[i] = 0;
           for (p = 0; p < SizeFilter; p++)
	   {
              Index = i + Start_G1 + p;
              Index = test_index(Index, N);
              Output[i] += Temp[Index]*F[p];
	   }
       }
    }
    delete [] Temp;
}

/****************************************************************************/

void SubBandFilter::transform (int N, float *SignalIn, float *SignalOut, float *DetailOut)
{
   if (DistPix == 1)
   {
       convol_h0 (N, SignalIn, SignalOut);
       convol_g0 (N, SignalIn, DetailOut);
   }
   else
   {
      int i,s,Dim;
      fltarray TabSignal(N);
      // TabSignal.NameObj = "TabSignal";
      fltarray TabSignalOut(N);
      // TabSignalOut.NameObj = "TabSignalOut";
      fltarray TabDetailOut(N);
      // TabDetailOut.NameObj = "TabDetailOut";


      for (i = 0; i < N/2; i++) DetailOut[i] = 0;
      for (s = 0; s < DistPix; s++)
      {
         Dim = 0;
         for (i = s; i < N; i += DistPix) TabSignal(Dim++) = SignalIn[i];
         TabSignalOut.reform((N+1)/2);
	 TabDetailOut.reform(N/2);
         convol_h0 (Dim, TabSignal.buffer(), TabSignalOut.buffer());
         convol_g0 (Dim, TabSignal.buffer(), TabDetailOut.buffer());
         for (i = 0; i < (Dim+1)/2; i ++)
         {
            int Ind = i * DistPix + s;
	    if (Ind < (N+1)/2)
	    {
 	       SignalOut[Ind] = TabSignalOut(i);
	    }
  	 }
         for (i = 0; i < Dim/2; i ++)
         {
            int Ind = i * DistPix + s;
            DetailOut[Ind] = TabDetailOut(i);
	 }
      } 
   }
}
       
/*************************************************************/

void SubBandFilter::recons (int N, float *SignalIn,  float *DetailIn, 
                           float *SignalOut)
{
    int i;
    float *Temp = new float [N];
    for (i=0; i < N; i++) Temp[i]=0;
    
    if (DistPix == 1)
    {
      convol_h1 (N, SignalIn, SignalOut);
      convol_g1 (N, DetailIn, Temp);
      for (i = 0; i< N; i++) 
             SignalOut[i] = NormCoef*(SignalOut[i]+Temp[i]);
    }
    else
    {
      int i,s,Dim,Ind;
      int v = N / DistPix;
      int r = N - v * DistPix;
      fltarray TabSignal(N);
      // TabSignal.NameObj = "REC TabSignal";
      fltarray TabSignalIn(N);
      // TabSignalIn.NameObj = "REC TabSignalIn";
      fltarray TabDetailIn(N);
      // TabDetailIn.NameObj = "REC  TabDetailIn";

      for (s = 0; s < DistPix; s++)
      {
         Dim = (s < r) ? v + 1: v;
 	 TabSignalIn.reform((Dim+1)/2);
	 TabDetailIn.reform(Dim/2);
	 TabSignal.reform(Dim);
         for (i = 0; i < (Dim+1)/2; i ++)
         {
            int Ind = i * DistPix + s;
	    if (Ind >= (N+1)/2) 
	    {
	       // cout << "Error: ind = " << Ind << endl;
	       // exit(-1);
	       Ind = (N+1)/2 - 1;
	    }
            TabSignalIn(i) = SignalIn[Ind];
	 } 
	 convol_h1 (Dim, TabSignalIn.buffer(), TabSignal.buffer());
         for (i = 0; i < Dim/2; i ++)
         {
            int Ind = i * DistPix + s;
            if (Ind >= N) Ind = N - 1;
	    TabDetailIn(i) = DetailIn[Ind];
 	 }	
	 convol_g1 (Dim, TabDetailIn.buffer(), Temp);
	 Ind = 0;
	 for (i = s; i < N; i += DistPix) 
	 {
	     if (Ind >= Dim)
	     {
	        cout << "Error: dim = " << Dim << " ind = " << Ind << endl;
	        exit(-1);
	     }
             SignalOut[i] = NormCoef*(TabSignal(Ind)+Temp[Ind]);
 	     Ind ++;
	 }	 
      }
    }
    delete [] Temp;
}

/*************************************************************/

void SubBandFilter::noise_transform (int N, float *SignalIn, float *SignalOut, 
                   float *DetailOut)
{
       noise_convol_h0 (N, SignalIn, SignalOut);
       noise_convol_g0 (N, SignalIn, DetailOut);
}

/*************************************************************/
// UNDECIMATED PART
/*************************************************************/

void SubBandFilter::convol_filter(int N, float *Input, 
          float *Output, float *F, int SizeFilter, int Start_Filter, int Step)
{
    int i, Index, p;    
    for (i = 0; i < N; i ++)
    {        
        Output[i] = 0;
        for (p  = 0; p < SizeFilter; p ++)
	{
           Index = i + (p  + Start_Filter) * Step;
           Index = test_index(Index, N);
           Output[i] +=  Input[Index]*F[SizeFilter-p-1];
	}
    }
}

/*************************************************************/

void SubBandFilter::rec_convol_filter(int N, float *Input, float *Output, 
                       float *F, int SizeFilter, int Start_Filter, int Step)
{
    int i, Index, p;
    double Val;
    double *Temp1 = new double [N];
    double *Temp2 = new double [N];
    
    for (i=0; i < N; i++) 
    {
       Temp1[i]=0;
       Temp2[i]=0;
    }
    for (i = 0; i < N; i += 2)  Temp1[i] = Input[i];
    for (i = 1; i < N; i += 2)  Temp2[i] = Input[i];
    
    for (i = 0; i < N; i++)
    {
        Val = 0;
        for (p = 0; p < SizeFilter; p++)
	{
            Index = i+ (p+Start_Filter)*Step;
            Index = test_index(Index, N);
            Val +=  Temp1[Index] * F[p] + Temp2[Index] * F[p];
 	}
	Output[i] = (float) Val;
    }
    delete [] Temp1;
    delete [] Temp2;    
}

/*************************************************************/

void SubBandFilter::transform (int N, float *SignalIn, float *SignalOut, 
                   float *DetailOut, int Step)
{
       convol_filter(N, SignalIn, SignalOut, H0, Size_H0, Start_H0, Step);
       convol_filter(N, SignalIn, DetailOut, G0, Size_G0, Start_G0, Step);
}


/*************************************************************/

void SubBandFilter::recons(int N, float *SignalIn,  float *DetailIn, 
                           float *SignalOut, int Step)
{
    int i;
    float *Temp = new float [N];
    
    for (i=0; i < N; i++) Temp[i]=0;
    rec_convol_filter(N, SignalIn, SignalOut, H1, Size_H1, Start_H1, Step);
    rec_convol_filter(N, DetailIn, Temp, G1, Size_G1, Start_G1, Step);
		           
    for (i = 0; i< N; i++) 
             SignalOut[i] = 0.5*NormCoef*(SignalOut[i]+Temp[i]);
    delete [] Temp;
}

/*************************************************************/

void SubBandFilter::transform (fltarray& Sig_in, fltarray& Sig_out, 
                           fltarray* Smooth) {
    int Nx = Sig_in.nx();
float *SigHigh = Sig_in.buffer();
    int Nx_2;
    int i;
    float* Sig_low, * Sig_det;

    Nx_2 = (Nx+1)/2;

    Sig_low = new float [Nx_2];
    Sig_det = new float [Nx_2];
   transform(Nx, SigHigh, Sig_low, Sig_det);


   (*Smooth).resize(Nx_2);
    
    for (i = 0; i< Nx_2; i++) {
       (*Smooth)(i) = Sig_low[i];
       if (i+Nx_2 < Nx) Sig_out(i) = Sig_det[i];
    }
    
    
    delete [] Sig_low;
    delete [] Sig_det;

}

/*************************************************************/

void SubBandFilter::recons (fltarray& Sig_rec, fltarray& Sig_in, 
                            fltarray& Smooth) {

    int Nx = Sig_rec.nx();
       
    int  Nx_2 = (Nx+1)/2;
          
    /* Allocation  */
    float* signal_rec = new float [Nx];
    float* signal_w= new float [Nx_2];
    float* signal_c= new float [Nx_2];

    for (int i = 0; i< Nx_2; i++) {     
       signal_c[i] = (Smooth)(i);
       signal_w[i] = (Sig_in)(i);
    }
    recons (Nx, signal_c, signal_w, signal_rec);    
    for (int i = 0; i< Nx; i++) {     
       Sig_rec(i) = signal_rec[i];
    }
    delete [] signal_rec;
    delete [] signal_w;
    delete [] signal_c;
    
}




/****************************************************************************/
// 1D full decimated decomposition
/****************************************************************************/

void Ortho_1D_WT::transform(fltarray &Data, fltarray &Mallat, int Nbr_Plan)
{			  
    int i,s;
    int Np_2,Nps;
    float *ImagLow, *ImagHigh, *DataResol;
    int Np = Data.nx();
      
   // cout << "mallat : " << Nbr_Plan << " np = " << Np << endl;
    // Data.NameObj = "transform data";
    // Mallat.NameObj = "transform mallat";
    
    if (Np !=  Mallat.n_elem()) Mallat.alloc(Np);
    
    /* Allocation  */
    Nps = size_resol(1, Np);
    DataResol    = new float[Np];
    ImagHigh     = new float[Nps];
    ImagLow      = new float[Nps];
 
    for (i=0; i< Np; i++) DataResol[i] = Data(i);
  
    /* Compute the wavelet coefficients */
    for (s = 0; s < Nbr_Plan-1; s++)
    {  
  //  cout << " step " << s << endl;
        Nps = size_resol(s, Np);
        Np_2 = (Nps+1)/2;
	
        Ptr_SB1D->transform(Nps,DataResol,ImagLow,ImagHigh);
        for (i=0; i < Nps/2; i++)  Mallat(Np_2+i) = ImagHigh[i];
 
  	if (s != Nbr_Plan-1)
             for (i = 0; i < Np_2; i++)  
	               DataResol[i] = ImagLow[i];
     }
     Np_2 = size_resol(Nbr_Plan-1, Np);
     for (i=0; i< Np_2; i++) Mallat(i) = ImagLow[i];
 // cout << "exit" << endl;
    delete [] ImagHigh;
    delete [] ImagLow;
    delete [] DataResol;
}

/****************************************************************************/

void Ortho_1D_WT::recons (fltarray &Mallat, fltarray &Data, int Nbr_Plan)
{
    int  Nps;
    register int i,s;
    float *image_h1;
    float *image;
    float *image_g1;
    int Np = Mallat.nx();
    int Np_2=Np;
    
    // Data.NameObj = "rec data";
    // Mallat.NameObj = "rec mallat";
    // cout << "rec mallat : " << Nbr_Plan << " np = " << Np << endl;
    
    if (Data.n_elem() != Np) Data.alloc(Np);

    /* Allocation */
    image    = new float [Np];
    image_h1 = new float [Np];
    image_g1 = new float [Np];
  
    /* last scale size */
    Np_2 = size_resol( Nbr_Plan-1, Np);
         
     /* initial image construction : image_h1_h1. */
    for (i=0; i< Np_2; i++) image_h1[i] = Mallat(i);

    for (s = Nbr_Plan-2; s >= 0; s--)
    {
    // cout << " step " << s << endl;
        Nps = size_resol(s, Np);
        Np_2 = (Nps+1)/2;
        for (i=0; i< Nps/2; i++) image_g1[i] = Mallat(i+Np_2);
        Ptr_SB1D->recons (Nps, image_h1, image_g1, image);
        /* Next iteration */
        for (i = 0; i < Nps; i++) image_h1[i] = image[i];
    }
   //  cout << "out" << endl;
    
    for (i=0; i< Np; i++) Data(i) = image[i];
    delete [] image;
    delete [] image_h1;
    delete [] image_g1;
}

/****************************************************************************/
// 1D full undecimated  decomposition
/****************************************************************************/

void PAVE_1D_WT::transform(fltarray &Data, fltarray &Mallat, int Nbr_Plan)
{
    int i,s;
    float *ImagLow, *ImagHigh, *DataResol;
    int Np = Data.nx();
       
   // cout << "mallat : " << Nbr_Plan << " np = " << Np << endl;
    
    if ((Np !=  Mallat.nx()) || (Nbr_Plan != Mallat.ny())) 
                                        Mallat.alloc(Np, Nbr_Plan);
     
    /* Allocation  */
    DataResol    = new float [Np];
    ImagHigh     = new float [Np];
    ImagLow      = new float [Np];
 
    for (i=0; i< Np; i++) DataResol[i] = Data(i);
  
    /* Compute the wavelet coefficients */
    for (s = 0; s < Nbr_Plan-1; s++)
    {  
       int Step = POW2(s);
       // cout << " step " << Step  << endl;
 	
       Ptr_SB1D->transform(Np,DataResol,ImagLow,ImagHigh, Step);
       for (i=0; i < Np; i++)  Mallat(i, s) = ImagHigh[i];
 
  	if (s != Nbr_Plan-1)
             for (i = 0; i < Np; i++) DataResol[i] = ImagLow[i];
     }
     for (i=0; i< Np; i++) Mallat(i, Nbr_Plan-1) = ImagLow[i];
    delete [] ImagHigh;
    delete [] ImagLow;
    delete [] DataResol;
}

/****************************************************************************/

void PAVE_1D_WT::recons (fltarray &Mallat, fltarray &Data, int Nbr_Plan)
{
    register int i,s;
    float *image_h1;
    float *image;
    float *image_g1;
    int Np = Mallat.nx();
 
    // cout << "rec mallat : " << Nbr_Plan << " np = " << Np << endl;
    
    if (Data.n_elem() != Np) Data.alloc(Np);

    /* Allocation */
    image    = new float [Np];
    image_h1 = new float [Np];
    image_g1 = new float [Np];
  
    /* initial image construction : image_h1_h1. */
    for (i=0; i< Np; i++) image_h1[i] = Mallat(i, Nbr_Plan-1);

    for (s = Nbr_Plan-2; s >= 0; s--)
    {
       int Step = POW2(s);
       // cout << " step " << Step << endl;
       for (i=0; i< Np; i++) image_g1[i] = Mallat(i, s);
       Ptr_SB1D->recons (Np, image_h1, image_g1, image, Step);
       /* Next iteration */
       for (i = 0; i < Np; i++) image_h1[i] = image[i];
    }
   //  cout << "out" << endl;
    
    for (i=0; i< Np; i++) Data(i) = image[i];
     
    delete [] image;
    delete [] image_h1;
    delete [] image_g1;
}
/****************************************************************************/
void PAVE_1D_WT::transform (fltarray& SignalIn, fltarray& Transf_out, 
                        fltarray& Smooth, int step) {
                        
    int Nx = SignalIn.nx();
    float *SignalHigh = SignalIn.buffer();
    int i;
    float *Signal_low, *Signal_det;
    Signal_low = new float[Nx];
    Signal_det = new float[Nx];
    
    Ptr_SB1D->transform(Nx, SignalHigh, Signal_low, Signal_det, step);
    
    for (i = 0; i<Nx; i++) {    
    Smooth(i) = Signal_low[i];
    Transf_out(i) = Signal_det[i];
    }
    delete [] Signal_low;
    delete [] Signal_det;
}
    
/****************************************************************************/        
void PAVE_1D_WT::recons (fltarray &Det, fltarray &Smooth, 
                     fltarray &Sig_Rec, int step) {
                     
    int i;
    float* signal_h1, * signal_g1;
    int Nx = Sig_Rec.nx();
       
    /* Allocation  */
    signal_h1 = new float [Nx];
    signal_g1 = new float [Nx];
    
    for (i = 0; i< Nx; i++){
        signal_h1[i] = Smooth(i);
        signal_g1[i] = Det(i);
        
    }
    float* Sig_recons = Sig_Rec.buffer();    
    Ptr_SB1D->recons (Nx, signal_h1, signal_g1, Sig_recons, step);
 
    delete [] signal_h1;
    delete [] signal_g1;

                 
}

float PAVE_1D_WT::getAbsMaxTransf (fltarray & WT_Trans, float SigmaNoise, Bool OnlyPositivDetect, Bool Verbose) 
{
   fltarray WTab1DSignifLevel(100);
   WTab1DSignifLevel(0) = 0.723490;
   WTab1DSignifLevel(1) = 0.285450;
   WTab1DSignifLevel(2) = 0.177948 ;
   WTab1DSignifLevel(3) = 0.122223;
   WTab1DSignifLevel(4) = 0.0858113;
   WTab1DSignifLevel(5) = 0.0605703;
   WTab1DSignifLevel(6) = 0.0428107;
   WTab1DSignifLevel(7) = 0.023;
   WTab1DSignifLevel(8) = 0.017;
   for (int i = 9; i < 100; i++) WTab1DSignifLevel(i) = 0.017;
	        
   float absMaxLevel=0;
   fltarray Scale(WT_Trans.nx());
   for (int s=0; s< WT_Trans.ny()-1; s++) 
   {
      float NoiseLevel = SigmaNoise * WTab1DSignifLevel(s);
      for (int i = 0; i< WT_Trans.nx(); i++) Scale(i) = WT_Trans(i,s);
      if (OnlyPositivDetect) Scale.inf_threshold(0.0);
      float MaxLevel = fabs(Scale.maxfabs());
      float NormMaxLevel = MaxLevel/NoiseLevel;
      if (absMaxLevel <= NormMaxLevel) absMaxLevel = NormMaxLevel;
   }
   if (Verbose) cout << "Lambda  Half Decimated: " << absMaxLevel << endl;

   return absMaxLevel;
}

/****************************************************************************/
// 1D full half decimated decomposition
/****************************************************************************/

void HALF_1D_WT::transform(fltarray &Data, fltarray &Mallat, int Nbr_Plan)
{			  
    int i,s;
    int Np_2,Nps;
    float *ImagLow, *ImagHigh, *DataResol;
    int Np = Data.nx();

         
   // cout << "mallat : " << Nbr_Plan << " np = " << Np << endl;
    // Data.NameObj = "transform data";
    // Mallat.NameObj = "transform mallat";
    
    if (Np !=  Mallat.n_elem()) Mallat.alloc(Np);
    
    /* Allocation  */
    Nps = size_resol(1, Np);
    DataResol    = new float[Np];
    ImagHigh     = new float[Nps];
    ImagLow      = new float[Nps];
 
    for (i=0; i< Np; i++) DataResol[i] = Data(i);
  
    /* Compute the wavelet coefficients */
    for (s = 0; s < Nbr_Plan-1; s++)
    {  
  //  cout << " step " << s << endl;
        Nps = size_resol(s, Np);
        Np_2 = (Nps+1)/2;
	
        Ptr_SB1D->transform(Nps,DataResol,ImagLow,ImagHigh);
        for (i=0; i < Nps/2; i++)  Mallat(Np_2+i) = ImagHigh[i];
 
  	if (s != Nbr_Plan-1)
             for (i = 0; i < Np_2; i++)  
	               DataResol[i] = ImagLow[i];
     }
     Np_2 = size_resol(Nbr_Plan-1, Np);
     for (i=0; i< Np_2; i++) Mallat(i) = ImagLow[i];
 // cout << "exit" << endl;
    delete [] ImagHigh;
    delete [] ImagLow;
    delete [] DataResol;
}


/****************************************************************************/
/*    Undecimated non-orthgonal filter bank (one step transform)
*
****************************************************************************/
            
const char * StringUndecFilter (type_undec_filter TypeUnderFilter)
{
    switch (TypeUnderFilter)
    {
        case U_B3SPLINE:
	      return ("SplineB3-Id+H:  H=[1,4,6,4,1]/16, Ht=H, G=Id-H, Gt=Id+H"); 
        case U_B3SPLINE_2:
	      return ("SplineB3-Id:  H=[1,4,6,4,1]/16, Ht=H, G=Id-H*H, Gt=Id");
        case U_B2SPLINE:
	      return ("SplineB2-Id: H=4[1,2,1]/4, Ht=H, G=Id-H*H, Gt=Id"); 
        case U_HAAR_B3S: 
             return ("Harr/Spline: [H,G]=Haar,Ht=[1,4,6,4,1]/16,Gt=[1,6,16,-6,-1]/16"); 
        case U_HAAR_B3S_POS:
	     return ("Harr/Spline POS: H=Haar,G=[-1/4,1/2,-1/4],Ht=[1,3,3,1]/8,Gt=[1,6,1]/4");
	default: 
	      return ("Undefined filters");
    }
} 

/***********************************************************************/

void usb_usage(type_undec_filter Filter)
{
    fprintf(OUTMAN, "         [-U type_of_non_orthog_filters]\n");
    for (int i = 1; i <= NBR_UNDEC_FILTER; i++)
    fprintf(OUTMAN, "              %d: %s \n",i,
                                           StringUndecFilter((type_undec_filter  ) (i-1)));
    fprintf(OUTMAN, "             default is %s\n\n", StringUndecFilter ((type_undec_filter) Filter));
     
}

/***********************************************************************/

void UndecSubBandFilter::init()
{
   Verbose = False;
   switch (TypeUnderFilter)
   {
      case U_B2SPLINE:
           {
	       // PHI = bi-linear interpolation
	      // H = B, B=b3 spline
	      // Ht = H
	      // G = Id - H*H
	      // Ht = Id
	      // it works fine
              float Coeff_h0 = 3. / 8.;
              float Coeff_h1 = 1. / 4.;
              float Coeff_h2 = 1. / 16.;
              SizeFilterH=SizeFilterTH=3;
	      SizeFilterG=5;
	      SizeFilterTG=1;
	      FilterH = new float[SizeFilterH];
	      FilterG = new float[SizeFilterG];
	      FilterTG = new float[SizeFilterTG];
 	      FilterTH = FilterH;
	      FilterH[0]=FilterH[2]=1./4.;
	      FilterH[1]=1./2.;
 	      FilterG[0]=FilterG[4]=-Coeff_h2;
	      FilterG[1]=FilterG[3]=-Coeff_h1;
	      FilterG[2]= 1.-Coeff_h0;
	      FilterTG[0] = 1;	      
// 	      for (int i=0; i < SizeFilterG; i++) printf("%5.2f, ", FilterG[i]); 
// 	      printf("\n");
// 	      for (int i=0; i < SizeFilterH; i++) printf("%5.2f, ", FilterH[i]); 
// 	      printf("\n");
// 	      for (int i=0; i < SizeFilterTH; i++) printf("%5.2f, ", FilterTH[i]);
// 	      printf("\n"); 
// 	      for (int i=0; i < SizeFilterTG; i++) printf("%5.2f, ", FilterTG[i]); 
// 	      printf("\n");
// 	      cout << " SPLINE 2" << endl;
	    }
	    break;
      case U_B3SPLINE_2:
           {
	      // H = B, B=b3 spline
	      // Ht = H
	      // G = Id - H*H
	      // Gt = Id
	      // it works fine
	      double Norm = 16.*16.;
              float Coeff_h0 = 70. / Norm;
              float Coeff_h1 = 56. / Norm;
              float Coeff_h2 = 28. / Norm;
              float Coeff_h3 =  8. / Norm;
              float Coeff_h4 =  1. / Norm;
              SizeFilterH=SizeFilterTH=5;
	      SizeFilterG=9;
	      SizeFilterTG=3;
	      FilterH = new float[SizeFilterH];
	      FilterG = new float[SizeFilterG];
	      FilterTG = new float[SizeFilterTG];
 	      FilterTH = FilterH;
	      FilterH[0]=FilterH[4]= 1. / 16.;
	      FilterH[1]=FilterH[3]= 1. / 4.;
	      FilterH[2]= 3. / 8;
 	      FilterG[0]=FilterG[8]=-Coeff_h4;
	      FilterG[1]=FilterG[7]=-Coeff_h3;
	      FilterG[2]=FilterG[6]=-Coeff_h2;
	      FilterG[3]=FilterG[5]=-Coeff_h1;
	      FilterG[4]= 1.-Coeff_h0;
	      FilterTG[1] = 1;
	      FilterTG[0] = FilterTG[2] =0;
	      
// 	      for (int i=0; i < SizeFilterG; i++) printf("%5.2f, ", FilterG[i]*Norm); 
// 	      printf("\n");
// 	      for (int i=0; i < SizeFilterH; i++) printf("%5.2f, ", FilterH[i]); 
// 	      printf("\n");
// 	      for (int i=0; i < SizeFilterTH; i++) printf("%5.2f, ", FilterTH[i]);
// 	      printf("\n"); 
// 	      for (int i=0; i < SizeFilterTG; i++) printf("%5.2f, ", FilterTG[i]); 
// 	      printf("\n");
// 	      cout << " SPLINE 3" << endl;
	    }
	    break;
       case U_B3SPLINE:
           {
               // H = B, B=b3 spline
	       // Ht = H
	       // G = Id - H
	       // Gt = Id + H
	       // it works fine:
              float Coeff_h0 = 3. / 8.;
              float Coeff_h1 = 1. / 4.;
              float Coeff_h2 = 1. / 16.;
              SizeFilterH=SizeFilterG=SizeFilterTH=SizeFilterTG=5;
	      FilterH = new float[SizeFilterH];
	      FilterG = new float[SizeFilterG];
	      FilterTH = FilterH;
	      FilterTG = new float[SizeFilterTG];
	      FilterH[0]=FilterH[4]=Coeff_h2;
	      FilterH[1]=FilterH[3]=Coeff_h1;
	      FilterH[2]=Coeff_h0;
	      FilterG[0]=FilterG[4]=-Coeff_h2;
	      FilterG[1]=FilterG[3]=-Coeff_h1;
	      FilterG[2]=1.-Coeff_h0;
	      float Coeff_Tg0 = 22./16.;
              float Coeff_Tg1 = 4./16.;
              float Coeff_Tg2 = 1./16.;
	      FilterTG[0]=FilterTG[4]=Coeff_Tg2;
	      FilterTG[1]=FilterTG[3]=Coeff_Tg1;
	      FilterTG[2]=Coeff_Tg0;
	    }
	    break;
      case U_HAAR_B3S:
            {
              SizeFilterH=SizeFilterG=3;
	      SizeFilterTH=SizeFilterTG=5;
	      FilterH = new float[SizeFilterH];
	      FilterG = new float[SizeFilterG];
	      FilterTH = new float[SizeFilterH];;
	      FilterTG = new float[SizeFilterTG];
	      FilterH[0]=FilterH[1]=1./2.;
	      FilterH[2]=0;
 	      FilterG[0]=-1./2;
	      FilterG[1]=1./2.;
	      FilterG[2]=0;
   	      FilterTG[0]= 1./ 16.;
	      FilterTG[1]= 6./ 16.;
	      FilterTG[2]= 1.;
	      FilterTG[3]= -6./ 16.;
	      FilterTG[4]= -1./ 16.;
              float Coeff_h0 = 3. / 8.;
              float Coeff_h1 = 1. / 4.;
              float Coeff_h2 = 1. / 16.;
 	      FilterTH[0]=FilterTH[4]=Coeff_h2;
	      FilterTH[1]=FilterTH[3]=Coeff_h1;
	      FilterTH[2]=Coeff_h0;
 	    }
	    break;
      case U_HAAR_B3S_POS:
            {
	      //  H=Haar,G=[-1/4,1/2,-1/4],Ht=[1,3,3,1]/8,Gt=[1,6,1]/4");
              SizeFilterH=3;
	      SizeFilterG=3;
	      SizeFilterTH=5;
	      SizeFilterTG=3;
	      FilterH = new float[SizeFilterH];
	      FilterG = new float[SizeFilterG];
	      FilterTH = new float[SizeFilterH];;
	      FilterTG = new float[SizeFilterTG];
	      FilterH[0]=FilterH[1]=1./2.;
	      FilterH[2]=0;
 	      FilterG[0]=-1./4.;
	      FilterG[1]=1./2.;
	      FilterG[2]=-1./4.;
   	      FilterTG[0]= 1./ 4.;
	      FilterTG[1]= 6./ 4.;
	      FilterTG[2]= 1./ 4.;
	      
	      FilterTH[0]= 0.;
	      FilterTH[1]= 1./ 8.;
	      FilterTH[2]= 3./ 8.;
	      FilterTH[3]= 3./ 8.;
	      FilterTH[4]= 1./ 8.;
  	    }
	    break;
     }
     if (Verbose == True)
     {
     	int i;
	 cout << "Undec. Filter Bank = " << StringUndecFilter(TypeUnderFilter) << endl;
    	 printf("\nH = [");
    	 for (i=0; i < SizeFilterH; i++) printf("%5.6f, ", FilterH[i]);
    	 printf("]\nG = [");
    	 for (i=0; i < SizeFilterG; i++) printf("%5.6f, ", FilterG[i]);
    	 printf("]\nTH = [");
    	 for (i=0; i < SizeFilterTH; i++) printf("%5.6f, ", FilterTH[i]);
    	 printf("]\nTG = [");
    	 for (i=0; i < SizeFilterTG; i++) printf("%5.6f, ", FilterTG[i]);
    	 printf("]\n");
     }


//     if (TypeNorm == NORM_L2) // L2 NORMALIZATION DOES NOT WORK WITH NON ORTHOGONAL FILTERS
//     {
//        float A=0.;
//        float B=0.;
//        for (int i=0; i < SizeFilterH; i++) A+= FilterH[i]*FilterH[i];
//        A = sqrt(A);
//        for (int i=0; i < SizeFilterG; i++) B+= FilterG[i]*FilterG[i];
//        B = sqrt(B);
//        for (i=0; i < SizeFilterH; i++)  FilterH[i] /= A;
//        if (FilterH != FilterTH)
//               for (i=0; i < SizeFilterTH; i++)  FilterTH[i] *= A;
//        for (i=0; i < SizeFilterG; i++)  FilterG[i] /= B;
//        if (FilterG != FilterTG)
//            for (i=0; i < SizeFilterTG; i++)  FilterTG[i] *= B;
//  
//          A=0.;
//          B=0.;
//        for (int i=0; i < SizeFilterH; i++) A+= FilterH[i]*FilterH[i];
//        for (int i=0; i < SizeFilterG; i++) B+= FilterG[i]*FilterG[i];
//        cout << "SUMH2 = " << A << " SIMG2 " << B << endl;
//     }
}

/****************************************************************************/

void UndecSubBandFilter::convol_filter(int N, float *Input, float *Output, float *F, int SizeFilter,  int Step)
{
    int i, Index, p;
    int Start_Filter = - SizeFilter / 2; 
    for (i = 0; i < N; i ++)
    {        
        Output[i] = 0;
        for (p  = 0; p < SizeFilter; p ++)
	{
           Index = i + (p  + Start_Filter) * Step;
           Index = test_index(Index, N);
           Output[i] +=  Input[Index]*F[p];
	}
    }
}

/****************************************************************************/

void UndecSubBandFilter::transform(int Nx, float *High, float *Low, float *Det, int Step)
{
   convol_filter(Nx, High, Low, FilterH, SizeFilterH, Step);
   convol_filter(Nx, High, Det, FilterG, SizeFilterG, Step);
   // for (int i=0;i<Nx; i++) Det[i] = High[i] - Low[i];
}

/****************************************************************************/

void UndecSubBandFilter::recons(int Nx, float *Low, float *Det, float *High, int Step)
{
    int i;
    fltarray NewDet(Nx);
    convol_filter(Nx, Low, High, FilterTH, SizeFilterTH, Step);
    convol_filter(Nx, Det, NewDet.buffer(), FilterTG, SizeFilterTG, Step);
    for (i=0;i<Nx; i++) High[i] += NewDet(i);
}


/****************************************************************************/

/****************************************************************************/
//                   Half Decimated Orthogonal ALGORITHM
/****************************************************************************/

void HALF_DECIMATED_1D_WT::free (fltarray *TabBand, int Nbr_Plan) {

    if (Nbr_Plan != 0) delete [] TabBand;
}

void HALF_DECIMATED_1D_WT::set_tabdec (int NumUndec, Bool* & TabDec, 
                                       int Nbr_Plan) {
                                       
    int i;
    int NU = (NumUndec < 0) ? Nbr_Plan: NumUndec;
    TabDec = new Bool [Nbr_Plan];
    for (i=0; i < MIN(NU,Nbr_Plan); i++) TabDec[i] = False;
    for (i=MIN(NU,Nbr_Plan); i < Nbr_Plan; i++) TabDec[i] = True;
}


int HALF_DECIMATED_1D_WT::alloc (fltarray * & TabBand, int NSig,
	                          int Nbr_Plan,int NumUndec)
{
    Bool *TabDec;
    set_tabdec(NumUndec, TabDec, Nbr_Plan);
    int NbrBand = alloc (TabBand,NSig,Nbr_Plan,TabDec);
    delete [] TabDec;
    return NbrBand;
}


int HALF_DECIMATED_1D_WT::alloc (fltarray * & TabBand, int NSig, 
                                 int Nbr_Plan, Bool *TabDec) {
                                 
    char ch[80];
    int Nx = NSig;
    int NbrBand=Nbr_Plan;
    FirstDetectScale = 0;
    TabBand = new fltarray [Nbr_Plan];
    
    for (int s=0; s<NbrBand-1; s++) {
    
        Bool Dec = TabDec[s];
	int Nlx = Nx;
	
        if (Dec == True)
           Nlx = (Nx+1)/2;
           
	sprintf (ch, "band_%d", s+1);
        TabBand[s].alloc (Nlx, ch);   
            	
	if (Dec == True)
            Nx = (Nx+1)/2;
     }
     sprintf (ch, "band_%d", NbrBand);
     TabBand[NbrBand-1].alloc(Nx, ch);
     return NbrBand;
}

void HALF_DECIMATED_1D_WT::transform (fltarray& Signal, fltarray* TabTrans) {
                                        
    Bool *TabDec;
    set_tabdec (NbrUndecimatedScale, TabDec, NbScale);
    transform (Signal, TabTrans, TabDec);
    delete [] TabDec;
}

void HALF_DECIMATED_1D_WT::transform (fltarray& Signal, fltarray* TabTrans, 
	                              Bool *TabDec) {
                                      
   int Step = 1;
   PAVE_1D_WT PWT(*Ptr_SB1D);
   SubBandFilter SBT(*Ptr_SB1D);
   int Nx = Signal.nx();
   fltarray Sig_Aux(Nx,"aux");
   for (int s=0; s<NbScale-1; s++) {
#if SBDEBUG
       cout << "Scale " << s+1 << " Step = " << Step << endl;
       if (s == 0) cout << "  signal: Nx = " << Signal.nx() 
                        << " sigma = " <<  Signal.sigma()<<  endl;
       else cout << "  signal: Nx = " <<   Sig_Aux.nx() 
                 << " sigma = " << Sig_Aux.sigma() << endl;
#endif   
      if (s == 0) {
      
         if (TabDec[s] == False){ 
	    PWT.transform (Signal, TabTrans[s], Sig_Aux, Step);
	 }
	 else
  	    SBT.transform (Signal, TabTrans[s], &Sig_Aux);
      }    
      else  {
      
         Ptr_SB1D->DistPix = Step;
         if (TabDec[s] == False) 
	    PWT.transform (Sig_Aux, TabTrans[s], Sig_Aux, Step);
         else {
	    SBT.transform (Sig_Aux, TabTrans[s], &Sig_Aux);
      }
      }
      
      if (TabDec[s] == False) Step *= 2;
#if SBDEBUG      
      cout << "  signal: Nx = " << TabTrans[s].nx() 
           << " sigma = " << TabTrans[s].sigma() << endl;
      cout << "  Smooth: Nx = " << Sig_Aux.nx() 
           << " sigma = " << Ima_Aux.sigma() << endl << endl;
#endif
   }
   
   TabTrans[NbScale-1] = Sig_Aux;
}

void HALF_DECIMATED_1D_WT::recons (fltarray* TabTrans, fltarray& Signal) {
                                   
    Bool *TabDec;
    set_tabdec (NbrUndecimatedScale, TabDec, NbScale);
    recons (TabTrans, Signal, TabDec);
    delete []TabDec;
}

void HALF_DECIMATED_1D_WT::recons (fltarray* TabTrans, fltarray& Signal, 
                                   Bool *TabDec) {
                                   
   int Step = 1;
   PAVE_1D_WT PWTR(*Ptr_SB1D);
   //SubBandFilter SBTR(*Ptr_SB1D);
   FilterAnaSynt *pFilter = new FilterAnaSynt(F_HAAR);
   SubBandFilter *psb1d = new SubBandFilter (*pFilter, NORM_L2);
   psb1d->Border = I_MIRROR;
   
   int Nx = Signal.nx();
   fltarray Signal_Aux(Nx,"aux");
   int *TabNx = new int [NbScale];
   
   for (int s=0; s<NbScale-1; s++)
      if (TabDec[s] == False) Step *= 2;
      
   TabNx[0] = Nx;
   
   for (int s=1; s<NbScale; s++) {
     
      if (TabDec[s-1] == True)
         TabNx[s] = (TabNx[s-1]+1)/2;
      
      else 
         TabNx[s] = TabNx[s-1];
         
   }
   
   for (int s=NbScale-2; s>=0; s--) {  
       
      fltarray Buff(TabNx[s], "Buff");
      // cout << "resize to " <<  TabNx[s]  << endl;
      if (TabDec[s] == False) Step /= 2;

#if SBDEBUG      
      cout << "Scale " << s+1 << " Step = " << Step << endl;
      cout << "  Signal: Nx = " << TabTrans[s].nx() 
           << " sigma = " << TabTrans[s].sigma() << endl;
      if (s == Nbr_Plan-2) 
         cout << "  Smooth: Nx = " << TabTrans[s+1].nx() 
              << " sigma = " << TabTrans[s+1].sigma() << endl;
      else 
         cout << "  Smooth: Nx = " <<  Signal_Aux.nx() 
              << " sigma = " <<  Signal_Aux.sigma() << endl;
#endif
      if (s == NbScale-2) {
      
         Ptr_SB1D->DistPix = Step;

         if (TabDec[s] == False)
            PWTR.recons (TabTrans[s], TabTrans[s+1], Buff, Step);
         else {
            //SBTR.recons (Buff, TabTrans[s], TabTrans[s+1]);
	    psb1d->recons (Buff, TabTrans[s], TabTrans[s+1]);
	    
         }
      }
      else  {

         Ptr_SB1D->DistPix = Step;
	 if (TabDec[s] == False) 
	    PWTR.recons (TabTrans[s], Signal_Aux, Buff, Step);
         else{
	    //SBTR.recons (Buff, TabTrans[s], Signal_Aux);
	    psb1d->recons(Buff, TabTrans[s], Signal_Aux);
         }
      }
      if (s==0) Signal = Buff;
      else {
         Signal_Aux.resize (Buff.nx());
         Signal_Aux = Buff;
      }
   }
   
}


void HALF_DECIMATED_1D_WT::threshold (fltarray* WT_Trans, float NSigma, 
                                      int IterNumber) {
                                                                         
   for (int s=0; s < NbScale-1; s++) {
     
      float NSig = (s < 3) ? NSigma + 1: NSigma;
      float Noise = SigmaNoise;
      float Level = NSig*Noise;
          
      for (int i=0; i <  (WT_Trans[s]).nx(); i++){
        
         float Coef = (WT_Trans[s])(i);  
	 (WT_Trans[s])(i) = update(Coef, Level, Noise);
         if (SuppressIsolatedPixel){
	     if (   !is_on_support ((WT_Trans[s])(i-1), Level, Noise)
	         && !is_on_support ((WT_Trans[s])(i+1), Level, Noise))
	        (WT_Trans[s])(i) = 0.;
	 }
	 if (OnlyPositivDetect) {
	      if ((WT_Trans[s])(i) < 0) (WT_Trans[s])(i) = 0;
	 }    
      }
     
      if (Write) {
         char Name[256];
         sprintf (Name,"HDWT_sc%d_iter%d", s, IterNumber);
         // fits_write_fltarr (Name, WT_Trans[s]);
      }
   }    
      
   if (Write) {
      char Name[256];
      sprintf (Name,"HDWT_sc%d_iter%d", NbScale-1, IterNumber);
      // fits_write_fltarr (Name, WT_Trans[NbScale-1]);
   }  
}



/****************************************************************************/
inline float HALF_DECIMATED_1D_WT::update (float CoefSol, float Threshold, 
                                  float SoftLevel) {
					 
  float NewCoef = hard_threshold(CoefSol, Threshold);
  if (UseNormL1 == True) {
       float SoftL =  SoftLevel*NSigmaSoft;
       float T2 = Threshold/2.;
       if (SoftL < T2)
            NewCoef = soft_threshold(NewCoef, SoftL);
       else NewCoef = soft_threshold(NewCoef, T2);
   }
   return NewCoef;  
} 

inline bool HALF_DECIMATED_1D_WT::is_on_support (float CoefSol, float Threshold, 
                                  float SoftLevel) {
					 
  float NewCoef = hard_threshold(CoefSol, Threshold);
  if (UseNormL1 == True) {
       float SoftL =  SoftLevel*NSigmaSoft;
       float T2 = Threshold/2.;
       if (SoftL < T2)
            NewCoef = soft_threshold(NewCoef, SoftL);
       else NewCoef = soft_threshold(NewCoef, T2);
   }
   return (!(NewCoef == 0));  
} 

      
void HALF_DECIMATED_1D_WT::KillScaleNotUsed (fltarray* WT_Trans, int FirstDetScale) {

   FirstDetectScale = FirstDetScale;
   if (FirstDetectScale != 0) {
      for (int j=0; j < FirstDetectScale; j++) {
         WT_Trans[j].init();
      }
   }
}

void HALF_DECIMATED_1D_WT::KillLastScale (fltarray* WT_Trans) {
   
   WT_Trans[NbScale-1].init();
}


float HALF_DECIMATED_1D_WT::getAbsMaxTransf (fltarray* WT_Trans) {

   float absMaxLevel=0;
   
   for (int s=0; s<NbScale-1; s++) {

      if (OnlyPositivDetect) (WT_Trans[s]).inf_threshold(0.0);
      float MaxLevel = fabs((WT_Trans[s]).maxfabs());
      float NormMaxLevel = MaxLevel/SigmaNoise;
                   
      if (Write)
         cout << "Lambda HDWT (" << s+1 << ") : " << NormMaxLevel 
              << ", Max amplitude : " << MaxLevel << endl;  
                                                              
      if (absMaxLevel <= NormMaxLevel) absMaxLevel = NormMaxLevel;
   }
   if (Verbose) cout << "Lambda  Half Decimated: " << absMaxLevel << endl;

   return absMaxLevel;
}

