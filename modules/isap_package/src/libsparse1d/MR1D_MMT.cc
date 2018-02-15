/*******************************************************************************
**
**    UNIT
**
**    Version: 3.2
**
**    Author: Jean-Luc Starck
**
**    Date:  96/06/13 
**    
**    File:  MR1D_Transf.cc
**
*****************************************************************************
**
**    DESCRIPTION 1D wavelet transform and reconstruction for the 
**    ----------- multiresolution median transform   
**   
*****************************************************************************
**
** void mr1d_median (fltarray &Signal, fltarray &W_1D, int N, int Nbr_Plan,
**                   int MedianWindowSize, type_border Border)
** Signal = IN: data
** W_1D = OUT: wavelet transform
** N = IN: number of data in the Signal
** Nbr_Plan = OUT: number of scales
** MedianWindowSize = IN: window size for median filtering
** Border = IN: type of border to use
**
****************************************************************************/

 
#include "MR1D_Obj.h"
// #include "NR.h"
// #include "IM_IO.h"

 
/***************************************************************************/

void smooth1d_b3spline(fltarray &Signal, fltarray &Result, int Step, type_border Border)
{
   int N = Signal.nx();
   int i,indi1,indi2,indi3,indi4;

   for (i = 0; i < N; i ++)
   {
      indi1 = border_ind_test (i - Step, N, Border);
      indi2 = border_ind_test (i + Step, N, Border);
      indi3 = border_ind_test (i - 2 * Step, N, Border);
      indi4 = border_ind_test (i + 2 * Step, N, Border);

      Result(i) = 0.0625 * ( Signal(indi3)+Signal(indi4))
                  + 0.25 * ( Signal(indi1)+Signal(indi2)) + 0.375 * Signal(i);
   }
}

/***************************************************************************/
/*
void smooth1d_median(fltarray &Signal, fltarray &Result,
                     int Window_Size, int Step, type_border Border)
{
   int N = Signal.nx();
   int i,k,ind_fen,ind;
   int Dim2,Dim4;
   float *fenetre;
   extern float hmedian(float *, int);

   fenetre = new float [Window_Size+1];
   Dim2 = Window_Size / 2;
   Dim4 = Dim2 / 2;
   for (i = 0; i < N; i ++)
   {
      ind_fen = 1;
      for (k = i - Dim2*Step; k <= i + Dim2*Step; k+=Step)
      {
         ind = border_ind_test(k, N, Border);
         fenetre[ind_fen++] = Signal(ind);
      }
      sort((unsigned long) (Window_Size), fenetre);
      Result(i) = fenetre[Dim2+1];
      
      //Result(i) = fenetre[Dim2+1];
      // Result(i) = hmedian(fenetre, Window_Size);
   }

   delete [] fenetre;
}
*/
/***************************************************************************/
/*
void mr1d_medianspline3_bis (fltarray &Signal, fltarray & W_1D, 
                         int Nbr_Plan, int WindowSize, type_border Border)
  multiresolution transform with a median and a b3-spline 
{
    int i,Num_Plan, D2n;
    int Step, N = Signal.nx();
    fltarray Data(N);
    fltarray SData(N);
    fltarray Temp(N);
    fltarray SignalRec(N);

    int D2 = WindowSize / 2;
    Data = Signal;
    SignalRec.init(0);
    for (Num_Plan = 0; Num_Plan < Nbr_Plan - 1; Num_Plan++)
    {
          float Sigm;

          Step = iround (pow((double)2., (double) Num_Plan));
          D2n = D2*Step;
          smooth1d_median(Data, SData, 2*D2n+1, 1, Border);
          Temp = Data - SData;
          Sigm = sigma(Temp);
cout << "Scale " << Num_Plan+1 << " Window Size = " << 2*D2n+1 << "  Sigma = " << Sigm << endl;

          for (i = 0; i < N; i++)
          {
              if (ABS(Temp(i)) < 3*Sigm) Temp(i) = 0.;
              else W_1D(i,Num_Plan) = Temp(i);
          }
          SignalRec += Temp;
          Data = SData;
    }
    for (i = 0; i < N; i++) W_1D(i,Nbr_Plan-1) = Data (i);
    SignalRec += Data;
    Data = Signal - SignalRec;

    fits_write_fltarr("xx_rec.fits",SignalRec);

    for (Num_Plan = 0; Num_Plan < Nbr_Plan - 1; Num_Plan++)
    {
          Step = iround (pow((double)2., (double) Num_Plan));
          D2n = D2*Step;
          smooth1d_b3spline(Data, SData, Step, Border);

          for (i = 0; i < N; i++) W_1D(i,Num_Plan) += Data(i) - SData(i);
          Data = SData;
    }
     for (i = 0; i < N; i++) W_1D(i,Nbr_Plan-1) += Data (i);

}
*/
/***************************************************************************/
/*
void mr1d_medianspline3 (fltarray &Signal, fltarray & W_1D, 
                         int Nbr_Plan, int WindowSize, type_border Border)
 multiresolution transform with a median and a b3-spline 
{
    int i,Num_Plan, D2n;
    int Step, N = Signal.nx();
    fltarray Data(N);
    fltarray SData(N);
    fltarray Temp(N);
    int D2 = WindowSize / 2;
    Data = Signal;
    for (Num_Plan = 0; Num_Plan < Nbr_Plan - 1; Num_Plan++)
    {
          float Sigm;

          Step = iround (pow((double)2., (double) Num_Plan));
          D2n = D2*Step;
cout << "Scale " << Num_Plan+1 << " Window Size = " << 2*D2n+1 << endl;
          smooth1d_median(Data, SData, 2*D2n+1, 1,Border);
          Temp = Data - SData;
          Sigm = sigma(Temp);

          for (i = 0; i < N; i++)
          {
              W_1D(i,Num_Plan) = Data(i) -  SData(i);
              if (ABS(W_1D(i,Num_Plan)) < FLOAT_EPSILON)
              {
                 int k;
                 double Dist, Wei, Val, TWei;
                 TWei = 0.;
                 Val = 0.;
                 for (k = i-D2n; k <= i+D2n; k++)
                    if ((k >= 0) && (k < N))
                    {
                        if (ABS(Data(k) - Data(i)) < 3*Sigm)
                        {
                            Dist = ABS(i-k) / (float) D2n * 2.;
                            Wei = im_b3_spline(Dist);
                            TWei += Wei;
                            Val += Wei*Data(k);
                        }
                    }
                  SData(i) = (float) (Val / TWei);
                 W_1D(i,Num_Plan) =  Data(i) - SData(i);
              }
          }


          Data = SData;
    }
     for (i = 0; i < N; i++) W_1D(i,Nbr_Plan-1) = Data (i);
}
*/

 /*********************  MMT ****************************/

static void tri_1d (float *fenetre, int Window_Size)
{
    int i,j;
    float Aux;

    for (i = 1; i < Window_Size; i++)
    for (j = 0; j < Window_Size - i; j++)
        if (fenetre[j] > fenetre[j+1])
        {
           Aux = fenetre[j+1];
           fenetre[j+1] = fenetre[j];
           fenetre[j] = Aux;
        }
}
void filt1d_mediane (const fltarray &S1, fltarray &S2, int N, int Window_Size,
                     type_border Border)
{
    int i,k,ind_fen,ind;
    int Dim,Dim2;
    float *fenetre;

    fenetre = new float [Window_Size];

    Dim = Window_Size;
    Dim2 = Dim / 2;
    for (i = 0; i < N; i++) 
    {
        ind_fen = 0;
        
        for (k = i -Dim2; k <= i+Dim2; k++)
        {
            ind = border_ind_test (k, N, Border);
            fenetre[ind_fen++] = S1(ind);
        }
        tri_1d (fenetre, Window_Size);
        S2(i) = fenetre[Dim2];
    }
    delete[] fenetre;
}

/***************************************************************************/

void mr1d_median (fltarray &Signal, fltarray &W_1D, int N, int Nbr_Plan, int MedianWindowSize, type_border Border)
/* multiresolution transform by the median */
{
    int i,Num_Plan;
    int Window_Size=MedianWindowSize;
    fltarray Data (N);
    fltarray Data_Sub (N);
    int Step = Window_Size/2; 
    for (i = 0; i < N; i++) Data(i) = Signal(i);

    for (Num_Plan = 0; Num_Plan < Nbr_Plan-1; Num_Plan++)
    {
          for (i = 0; i < N; i++) W_1D(i, Num_Plan) = Data (i);

          Window_Size = 2 * Step + 1;
          Step *= 2;
          filt1d_mediane (Data, Data_Sub, N, Window_Size, Border);

          /* Calcule the multiresolutions coefficients */
          for (i = 0; i < N; i++) W_1D(i,Num_Plan) -= Data_Sub (i);
          Data = Data_Sub;
    }
    
    /* copy the low resolution */
    for (i = 0; i < N; i++) W_1D(i,Nbr_Plan-1) = Data (i);
}

/***************************************************************************/

 
