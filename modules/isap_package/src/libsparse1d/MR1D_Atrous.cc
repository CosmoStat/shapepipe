/*******************************************************************************
**
**    UNIT
**
**    Version: 1.0
**
**    Author: Jean-Luc Starck
**
**    Date:  98/04/20 
**    
**    File:  MR1D_Atrous.cc
**
*****************************************************************************
**
**    DESCRIPTION 1D wavelet transform and reconstruction
**    ----------- using wavelet derived from a Spline   
**   
*****************************************************************************
** 
** wave_1d_spline1 (fltarray &Signal, fltarray & W_1D, int N, int Nbr_Plan)
**
** wavelet transform, a trous algorithm,  with a b1-spline 
** Signal = IN: data
** W_1D = OUT: wavelet transform
** N = IN: number of data in the Signal
** Nbr_Plan = OUT: number of scales
**
****************************************************************************
** 
** wave_1d_spline3 (fltarray &Signal, fltarray & W_1D, int N, int Nbr_Plan)
**
** wavelet transform with a b3-spline, a trous algorithm
**
** Signal = IN: data
** W_1D = OUT: wavelet transform
** N = IN: number of data in the Signal
** Nbr_Plan = OUT: number of scales
** 
***************************************************************************
** 
** wave_1d_algo_trou_rec (fltarray &W_1D, fltarray & Signal, int N, int Nbr_Plan)
**
** reconstruction from a a-trou algorithm 
**
** Signal = IN: data
** W_1D = OUT: wavelet transform
** N = IN: number of data in the Signal
** Nbr_Plan = OUT: number of scales
**
****************************************************************************
** 
** wave_1d_linear (fltarray &Signal, fltarray & W_1D, int N, int Nbr_Plan)
**
** linear wavelet transform, a trous algorithm
**
** Signal = IN: data
** W_1D = OUT: wavelet transform
** N = IN: number of data in the Signal
** Nbr_Plan = OUT: number of scales
**
****************************************************************************
**
** void wave_1d_B3deriv_atrou (fltarray &Signal, fltarray & W_1D,   
**                           int Nbr_Plan, type_border Border)
**  a trous wavelet transform with the mallat algorithm
**  the wavelet is the derivative of a spline 
** 
** Signal = IN: data
** W_1D = OUT: wavelet transform
** Nbr_Plan = OUT: number of scales
** Border = IN: type of border to use
**
****************************************************************************
**
** void wave_1d_rec_B3deriv_atrou (fltarray &W_1D, fltarray & Signal, 
**                                int Nbr_Plan, type_border Border)
**
** Signal = IN: data
** W_1D = OUT: wavelet transform
** Nbr_Plan = OUT: number of scales
** Border = IN: type of border to use
** 				
****************************************************************************/

#include "MR1D_Obj.h"
// #include "NR.h"

#define SIZE_FILTER 7

static double Tab_H1[SIZE_FILTER] = {0,0,0.125,0.375, 0.375, 0.125,0};
static double Tab_H2[SIZE_FILTER] = {0,0.125,0.375, 0.375, 0.125,0, 0};

static double Tab_G1[SIZE_FILTER] = {0,0,0,-0.5, 0.5, 0,0};
static double Tab_G2[SIZE_FILTER] = {0.03125, 0.21875, 0.6875, 
                                  -0.6875,  -0.21875,  -0.03125, 0};
				     

/***************************************************************************/

static void convol1d (fltarray &Signal, double *Filter, int SizeFilter,
                      fltarray &Result, int Step, type_border Border)
{
   int N = Signal.nx();
   int i,j,Ind;
   int DemiSize = SizeFilter/2;
   double Val;
   for (i = 0; i < N; i ++)
   {
      Val = 0;
      for (j= -DemiSize; j <= DemiSize; j++)
      {
         Ind = border_ind_test (i+j*Step, N, Border);
         Val += (double) Signal(Ind) * Filter[j+DemiSize];
      }
      Result(i) = (float) Val;
   }
}
  

/***************************************************************************/

void wave_1d_B3deriv_atrou (fltarray &Signal, fltarray & W_1D,   
                           int Nbr_Plan, type_border Border)
/* a trous wavelet transform with the mallat algorithm
   the wavelet is the derivative of a spline */
{
    int i,Num_Plan;
    int Step;
    int N = Signal.nx();
    fltarray Data(N);
    fltarray DataH(N);
    fltarray DataG(N);
    
    Data = Signal;
    for (Num_Plan = 0; Num_Plan < Nbr_Plan - 1; Num_Plan++)
    {
       Step = iround (pow((double)2., (double) Num_Plan));
       convol1d(Data,Tab_H1,SIZE_FILTER,DataH, Step, Border);
       convol1d(Data,Tab_G1,SIZE_FILTER,DataG, Step, Border);
       Data = DataH;
       for (i = 0; i < N; i ++) W_1D(i,Num_Plan) = DataG(i);  
    }
    
    /* copy the low resolution */
    for (i = 0; i < N; i++) W_1D(i,Nbr_Plan-1) = DataH(i);
}

/***************************************************************************/

void wave_1d_rec_B3deriv_atrou (fltarray &W_1D, fltarray & Signal, 
                                int Nbr_Plan, type_border Border)
/* reconstruction from a a-trou algorithm */
{
    int i,Num_Plan;
    int Step;
    int N = Signal.nx();
    fltarray DataG1(N);
    fltarray DataG2(N);
    fltarray DataH(N);
    for (i = 0; i < N; i++)  Signal(i) = W_1D(i,Nbr_Plan-1);
    for (Num_Plan = Nbr_Plan - 2; Num_Plan >= 0; Num_Plan--)
    {
       Step = iround (pow((double)2., (double) Num_Plan));
       for (i = 0; i < N; i++)  DataG1(i) =  W_1D(i, Num_Plan);
       convol1d(Signal,Tab_H2,SIZE_FILTER,DataH, Step, Border);
       convol1d(DataG1,Tab_G2,SIZE_FILTER,DataG2, Step, Border);
       for (i = 0; i < N; i++) Signal(i) = DataH(i) + DataG2(i);
    }
}
 

/*************************************************************************/

void wave_1d_spline3 (fltarray &Signal, fltarray & W_1D, int N, int Nbr_Plan,
                          type_border Border)
/* wavelet transform with a b3-spline */
{
    int i,indi1,indi2,indi3,indi4,Num_Plan;
    int Step;
    float *Data;

    Data = new float [N];
    for (i = 0; i < N; i++) Data[i] = Signal(i);

    for (Num_Plan = 0; Num_Plan < Nbr_Plan - 1; Num_Plan++)
    {
          /* copy the image in a cube */
          for (i = 0; i < N; i++) W_1D (i,Num_Plan) = Data [i];

          Step = iround (pow((double)2., (double) Num_Plan));

          for (i = 0; i < N; i ++)
          {
              indi1 = border_ind_test (i - Step, N, Border);
              indi2 = border_ind_test (i + Step, N, Border);
              indi3 = border_ind_test (i - 2 * Step, N, Border);
              indi4 = border_ind_test (i + 2 * Step, N, Border);

              Data[i] =
                  0.0625 * ( W_1D(indi3,Num_Plan)+W_1D(indi4,Num_Plan))
                  + 0.25 * ( W_1D(indi1,Num_Plan)+W_1D(indi2,Num_Plan))
                          + 0.375 * W_1D(i,Num_Plan);
          }

          /* Calcule the wavelet coefficients */
          for (i = 0; i < N; i++) W_1D(i,Num_Plan) -= Data [i];
    }
    
    /* copy the low resolution */
    for (i = 0; i < N; i++) W_1D(i,Nbr_Plan-1) = Data [i];

    delete[] Data;
}

/*************************************************************************/

void wave_1d_spline3_gen2 (fltarray &Signal, fltarray & W_1D, int N, int Nbr_Plan,
                      type_border Border)
/* wavelet transform with a b3-spline */
{
    int i,indi1,indi2,indi3,indi4,Num_Plan;
    int Step;
    float *Data;
    float *Data_aux;
    
    Data = new float [N];
    Data_aux = new float [N];
    for (i = 0; i < N; i++) Data[i] = Signal(i);
    
    for (Num_Plan = 0; Num_Plan < Nbr_Plan - 1; Num_Plan++)
    {
        /* copy the image in a cube */
        for (i = 0; i < N; i++) W_1D (i,Num_Plan) = Data [i];
        
        Step = iround (pow((double)2., (double) Num_Plan));
        
        
        for (i = 0; i < N; i ++)
        {
            indi1 = border_ind_test (i - Step, N, Border);
            indi2 = border_ind_test (i + Step, N, Border);
            indi3 = border_ind_test (i - 2 * Step, N, Border);
            indi4 = border_ind_test (i + 2 * Step, N, Border);
            
            Data_aux[i] =
            0.0625 * ( Data[indi3]+Data[indi4])
            + 0.25 * ( Data[indi1]+Data[indi2])
            + 0.375 * Data[i];
        }
        
        
        
        
        for (i = 0; i < N; i ++)
        {
            indi1 = border_ind_test (i - Step, N, Border);
            indi2 = border_ind_test (i + Step, N, Border);
            indi3 = border_ind_test (i - 2 * Step, N, Border);
            indi4 = border_ind_test (i + 2 * Step, N, Border);
            
            Data[i] =
            0.0625 * ( Data_aux[indi3] + Data_aux[indi4])
            + 0.25 * ( Data_aux[indi1] + Data_aux[indi2])
            + 0.375 * Data_aux[i];
        }
        
        
        /* Calcule the wavelet coefficients */
        for (i = 0; i < N; i++) W_1D(i,Num_Plan) -= Data [i];    
        
        for (i = 0; i < N; i++) Data[i] = Data_aux[i];
    }
    
    /* copy the low resolution */
    for (i = 0; i < N; i++) W_1D(i,Nbr_Plan-1) = Data [i];
    
    delete[] Data;
}


/***************************************************************************/

void wave_1d_spline3_gen2_rec (fltarray &W_1D, fltarray & Signal, int N, int Nbr_Plan)
/* reconstruction from a a-trou algorithm */
{
    int i,j,indi1,indi2,indi3,indi4,Step;
    
    float *Data;
    Data = new float [N];
    
    for (i = 0; i < N; i++) Signal(i) = W_1D(i,Nbr_Plan-1); 
    
    for (j = Nbr_Plan-2; j >= 0; j--)
    {
        
         Step = iround (pow((double)2., (double) j));
        
        for (i = 0; i < N; i ++)
        {
            indi1 = border_ind_test (i - Step, N, DEFAULT_BORDER);
            indi2 = border_ind_test (i + Step, N, DEFAULT_BORDER);
            indi3 = border_ind_test (i - 2 * Step, N, DEFAULT_BORDER);
            indi4 = border_ind_test (i + 2 * Step, N, DEFAULT_BORDER);
            
            Data[i] =
            0.0625 * ( Signal(indi3) + Signal(indi4))
            + 0.25 * ( Signal(indi1) + Signal(indi2))
            + 0.375 * Signal(i);
        }
        
        for (i = 0; i < N; i++) Signal(i) = Data[i] + W_1D(i,j); 
    }
}

/***************************************************************************/

void wave_1d_algo_trou_rec (fltarray &W_1D, fltarray & Signal, int N, int Nbr_Plan)
/* reconstruction from a a-trou algorithm */
{
    int i,j;

    for (i = 0; i < N; i++)
    {
        Signal(i) = 0.;
        for (j = 0; j < Nbr_Plan; j++) Signal(i) += W_1D(i,j);
    }
}

/***************************************************************************/

void wave_1d_spline1 (fltarray &Signal, fltarray & W_1D, int N, int Nbr_Plan, type_border Border)
/* wavelet transform with a b1-spline */
{
    int i,indi1,indi2,Num_Plan;
    int Step;
    float *Data;

//    printf ("Nbr_Plan = %d\n", Nbr_Plan);

    Data = new float [N];
    for (i = 0; i < N; i++) Data[i] = Signal(i);

    for (Num_Plan = 0; Num_Plan < Nbr_Plan - 1; Num_Plan++)
    {
          /* copy the image in the cube */
          for (i = 0; i < N; i++) W_1D(i,Num_Plan) = Data [i];

          Step = iround (pow((double)2., (double) Num_Plan));

          for (i = 0; i < N; i ++)
          {
              indi1 = border_ind_test (i - Step, N, Border);
              indi2 = border_ind_test (i + Step, N, Border);
              Data[i] =(0.5*(W_1D(indi1,Num_Plan)+W_1D(indi2,Num_Plan))
                          + 2. * W_1D (i,Num_Plan)) / 3.;
          }

          /* Calcul des coefficients d'ondelettes */
          for (i = 0; i < N; i++) W_1D (i,Num_Plan) -= Data [i];
    }
    
    /* copy the law resolution signal in the cube */
    for (i = 0; i < N; i++) W_1D (i,Nbr_Plan-1) = Data [i];

    delete [] Data;
}

/***************************************************************************/

void wave_1d_linear (fltarray &Signal, fltarray & W_1D, int N, int Nbr_Plan, type_border Border)
/* linear wavelet transform */
{
    int i,indi1,indi2,Num_Plan;
    int Step;
    float *Data;

//    printf ("Nbr_Plan = %d\n", Nbr_Plan);

    Data = new float [N];
    for (i = 0; i < N; i++) Data[i] = Signal(i);
    for (Num_Plan = 0; Num_Plan < Nbr_Plan - 1; Num_Plan++)
    {
          /* copy the data in the cube */
          for (i = 0; i < N; i++) W_1D(i, Num_Plan) = Data [i];

          Step = iround (pow((double)2., (double) Num_Plan));
          for (i = 0; i < N; i ++)
          {
              indi1 = border_ind_test (i - Step, N, Border);
              indi2 = border_ind_test (i + Step, N, Border);
              Data[i] = 0.25 * (W_1D(indi1,Num_Plan) + W_1D(indi2,Num_Plan)) 
                       + 0.5 * W_1D(i,Num_Plan);
          }

          /* calcul the wavelet coefficients */
          for (i = 0; i < N; i++) W_1D(i,Num_Plan) -= Data [i];
    }
    
    /* copy the low resolution signal in the cube */
    for (i = 0; i < N; i++) W_1D(i,Nbr_Plan-1) = Data [i];

    delete [] Data;
}
/***************************************************************************/

void wave_1d_haar (fltarray &Signal, fltarray & W_1D, int N, int Nbr_Plan, type_border Border)
/* linear wavelet transform */
{
    int i,indi1,indi2,Num_Plan;
    int Step;
    float *Data;

//    printf ("Nbr_Plan = %d\n", Nbr_Plan);

    Data = new float [N];
    for (i = 0; i < N; i++) Data[i] = Signal(i);
    for (Num_Plan = 0; Num_Plan < Nbr_Plan - 1; Num_Plan++)
    {
          /* copy the data in the cube */
          for (i = 0; i < N; i++) W_1D(i, Num_Plan) = Data [i];

          Step = iround (pow((double)2., (double) Num_Plan));
          for (i = 0; i < N; i ++)
          {
              indi1 = border_ind_test (i , N, Border);
              indi2 = border_ind_test (i - Step, N, Border);
              Data[i] = 0.5 * (W_1D(indi1,Num_Plan) + W_1D(indi2,Num_Plan)); 
	  }

          /* calcul the wavelet coefficients */
          for (i = 0; i < N; i++) W_1D(i,Num_Plan) -= Data [i];
    }
    
    /* copy the low resolution signal in the cube */
    for (i = 0; i < N; i++) W_1D(i,Nbr_Plan-1) = Data [i];

    delete [] Data;
}

/*********************************************************************/

void wp1d_atrous_transform (fltarray &Data, fltarray &WP,
                            int Nstep, type_border Border,
			    int StartPos, int NumStep)
{
   double Tab[5] = {0.0625, 0.25, 0.375, 0.25, 0.0625};
   int i,Np = Data.nx();
   fltarray  LowResol(Np);
   fltarray  HighResol(Np);
   int Nb = (int) pow(2., (double) (Nstep-NumStep));
   int  Step = (int) (pow((double)2., (double) NumStep));

// cout << "Step = " << Step << endl;

   convol1d (Data, Tab, 5, LowResol, Step, Border);
   HighResol = Data - LowResol;
   
   if (NumStep < Nstep-1)	
   {
      wp1d_atrous_transform (HighResol, WP, Nstep, Border, StartPos, 
                             NumStep+1);
      wp1d_atrous_transform (LowResol, WP, Nstep, Border,
                             StartPos+Nb/2, NumStep+1);
   }
   else
   {
// cout << "pos " << StartPos << " " << StartPos+1 << endl;
      for (i=0; i< Np; i++) 
      {
         WP(i, StartPos) = HighResol(i);
	 WP(i, StartPos+1) = LowResol(i);
      }
   }
}
