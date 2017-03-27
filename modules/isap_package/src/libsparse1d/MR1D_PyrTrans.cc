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
**    DESCRIPTION 1D wavelet transform and reconstruction for pyramidal 
**    ----------- transform   
**   
****************************************************************************/

 
#include "MR1D_Obj.h"
// #include "NR.h"
 
#define C1 .00245
#define C2 -.00915
#define C3 .03414
#define C4 -0.12720
#define C5 0.60048

#define CS ((C1+C2+C3+C4+C5)*2.)
#define C1n (C1/CS)
#define C2n (C2/CS)
#define C3n (C3/CS)
#define C4n (C4/CS)
#define C5n (C5/CS)

static float Tab_Coef_Inter[10] = {C1n,C2n,C3n,C4n,C5n,C5n,C4n,C3n,C2n,C1n};


 
 
/***************************************************************************/


static void extand_sample (fltarray &Data, fltarray &Data_Ext, int Ns)
{
   int i,k,ind;
   int N = Ns/2+Ns%2;

    for (i = 0; i < N; i++)
    {
        Data_Ext(2*i) = Data(i);
        Data_Ext(2*i+1) = 0.;

        for (k = 0; k < 10; k++)
        {
            ind = i - 4 + k;
            if (ind < 0) ind = 0;
            else if (ind >= N) ind = N - 1;
            Data_Ext(2*i+1) += Data(ind)*Tab_Coef_Inter [k];
        } 
       /* if (i < N-1) Data_Ext(2*i+1) = 0.5 * Data(i) + 0.5 * Data(i+1);
        else Data_Ext(2*i+1) = Data(i);*/
    }

    Data_Ext(2*(N-1)) = Data(N-1);
    Data_Ext(2*N-1) = Data(N-1);
    Data_Ext(Ns-1) = Data(N-1);
}

/***************************************************************************/

static void reduc_sample (fltarray &Data, fltarray &Data_Sub, int N)
{
   int i;

    for (i = 0; i < N/2+N%2; i++) Data_Sub (i)  = Data (2*i);
    for (i = N/2+N%2; i < N; i++) Data_Sub (i)  = 0.;
}

/***************************************************************************/

void mr1d_pyr_median (fltarray &Signal, fltarray &W_1D, int Np, int Nbr_Plan,
                      int WinSize, type_border Border)
/* pyramidal multiresolution transform by the median */
{
    int i,Num_Plan;
    fltarray Data (Np);
    fltarray Data_Sub (Np);
    int N = Np;

    for (i = 0; i < N; i++) Data_Sub(i) = Signal(i);

    for (Num_Plan = 0; Num_Plan < Nbr_Plan; Num_Plan++)
            for (i = 0; i < N; i++) W_1D(i,Num_Plan) = 0.;
    
    for (Num_Plan = 0; Num_Plan < Nbr_Plan - 1; Num_Plan++)
    {
          /* copy the image in a cube */
          for (i = 0; i < N; i++) W_1D(i,Num_Plan) = Data_Sub (i);
          filt1d_mediane (Data_Sub, Data, N, WinSize, Border);
          reduc_sample (Data, Data_Sub, N);
          extand_sample (Data_Sub, Data, N);

          /* computes the multiresolutions coefficients */
          for (i = 0; i < N; i++) W_1D(i,Num_Plan) -= Data (i);
          N =  N/2 + N %2;
    }
    /* copy the low resolution */
    for (i = 0; i < N; i++) W_1D(i,Nbr_Plan-1) = Data_Sub (i);
}

/***************************************************************************/

void pyr_1d_spline3 (fltarray &Signal, fltarray & W_1D, int Np, int Nbr_Plan,
                    type_border Border)
/* pyramidal wavelet transform with a b3-spline */
{
    int i,indi1,indi2,indi3,indi4,Num_Plan;
    int Step;
    fltarray Data (Np);
    fltarray Data_Sub (Np);
    int N = Np;

    for (i = 0; i < N; i++) Data_Sub(i) = Signal(i);

    for (Num_Plan = 0; Num_Plan < Nbr_Plan - 1; Num_Plan++)
    {
          /* copy the image in a cube */
          for (i = 0; i < N; i++) W_1D (i,Num_Plan) = Data_Sub (i);

          Step = 1;
          for (i = 0; i < N; i ++)
          {
              indi1 = border_ind_test (i - Step, N, Border);
              indi2 = border_ind_test (i + Step, N, Border);
              indi3 = border_ind_test (i - 2 * Step, N, Border);
              indi4 = border_ind_test (i + 2 * Step, N, Border);

              Data(i) =
                  0.0625 * ( W_1D(indi3,Num_Plan)+W_1D(indi4,Num_Plan))
                  + 0.25 * ( W_1D(indi1,Num_Plan)+W_1D(indi2,Num_Plan))
                          + 0.375 * W_1D(i,Num_Plan);
          }

          /* Calcule the wavelet coefficients */
          for (i = 0; i < N; i++) W_1D(i,Num_Plan) -= Data (i);

          reduc_sample (Data, Data_Sub, N);
          N /= 2;
    }
    
    /* copy the low resolution */
    for (i = 0; i < N; i++) W_1D(i,Nbr_Plan-1) = Data_Sub (i);
}

/***************************************************************************/

void pyr_1d_linear (fltarray &Signal, fltarray & W_1D, int Np, int Nbr_Plan, type_border Border)
/* pyramidal wavelet transform with a linear wavelet */
{
    int i,indi1,indi2,Num_Plan;
    int Step;
    fltarray Data (Np);
    fltarray Data_Sub (Np);
    int N = Np;

    for (i = 0; i < N; i++) Data_Sub(i) = Signal(i);

    for (Num_Plan = 0; Num_Plan < Nbr_Plan - 1; Num_Plan++)
    {
          /* copy the image in a cube */
          for (i = 0; i < N; i++) W_1D (i,Num_Plan) = Data_Sub (i);

          Step = 1;
          for (i = 0; i < N; i ++)
          {
              indi1 = border_ind_test (i - Step, N, Border);
              indi2 = border_ind_test (i + Step, N, Border);
              Data(i) = 0.25 * (W_1D(indi1,Num_Plan) + W_1D(indi2,Num_Plan)) 
                       + 0.5 * W_1D(i,Num_Plan);
          }

          /* Calcule the wavelet coefficients */
          for (i = 0; i < N; i++) W_1D(i,Num_Plan) -= Data (i);

          reduc_sample (Data, Data_Sub, N);
          N /= 2;
    }
    
    /* copy the low resolution */
    for (i = 0; i < N; i++) W_1D(i,Nbr_Plan-1) = Data_Sub (i);
}

/***************************************************************************/

void mr1d_pyr_rec (fltarray &Signal, fltarray &W_1D, int Np, int Nbr_Plan)
{
    int i,Num_Plan;
    int N=Np;
    fltarray Ext (Np);

    for (i = 0; i < Np; i++) Signal(i) = W_1D(i,Nbr_Plan-1); 
    for (Num_Plan = Nbr_Plan - 2; Num_Plan >= 0; Num_Plan--)
    {
          N = Np;
          for (i = 0; i < Num_Plan; i++) N = N/2+N%2;

          for (i = 0; i < N/2+N%2; i++)  Ext(i) = Signal(i);
          extand_sample (Ext, Signal, N);
          for (i = 0; i < N; i++) Signal(i) += W_1D(i,Num_Plan);
    }
}

/***************************************************************************/

 
