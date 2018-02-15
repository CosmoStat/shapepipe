/*****************************************************************************
**                   Copyright (C) 2003 by CEA
******************************************************************************
**
**    UNIT
**
**    Version: 1.0
**
**    Author: Jean-Luc Starck
**
**    Date:  30/03/03 
**    
**    File:  DefFunc.cc
**
*****************************************************************************/

#include "DefFunc.h"

static double MSVST_T1[] = {1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1};
static double MSVST_T2[] = {1, 0.2734375, 0.12347412109375, 0.0603594779968262, 0.0300147198140621, \
                   0.014986945054261, 0.00749092722048772, 0.00374514565088013, \
                   0.00187253308688793, 0.00093626157632392, 0.000468130167278172 \
                  };
static double MSVST_T3[] = {1, 0.08447265625, 0.017338752746582, 0.00414866860955954, 0.00102616434651281, \
                   0.00025586256661736, 6.39233749606246e-5, 1.59782042631771e-5, \
                   3.99438613268932e-6, 9.9858622538761e-7, 2.49645912118706e-7 \
                  };
static double MSVST_sumhh[] = {
         0, // dummy number
         0.375,
             0.15771484375,
        0.0760841369628906,
        0.0377178490161896,
        0.0188190417829901,
       0.00940455525596917,
       0.00470165753587537,
       0.00235075127553241,
       0.00117536595181247,
      0.000587681765180673,
      0.000293840731250224
};
            		   
/****************************************************************************/

void make_gaussian3d(fltarray & Result, float Sigma, int Indi, int Indj, int Indk)
// Create a 3D Gaussian at position  Indi,Indj,Indk), with a standard deviation
// equals to Sigma
{
   int Nx = Result.nx();
   int Ny = Result.ny();
   int Nz = Result.nz();
   int Dist,Delta_i,Delta_j;
   float sigma2,Energ = 0.;
   register int i,j,k;
   int Ci=Indi, Cj=Indj, Ck=Indk;
    if (Ci < 0) Ci = Nx / 2;
    if (Cj < 0) Cj = Ny / 2;
    if (Ck < 0) Ck = Nz / 2;
    sigma2 = -2. * Sigma * Sigma;
    for (i=0; i < Nx; i++)
    {
        Delta_i = (i - Ci) * (i - Ci);
        for (j = 0; j < Ny; j++)
	{
	    Delta_j = Delta_i + (j - Cj) * (j - Cj);
	    for (k = 0; k < Nz; k++)
	    {
	       Dist = Delta_j + (k - Ck) * (k - Ck);
	       Result(i,j,k) = exp((double)((float) Dist / sigma2));
	       Energ +=  Result(i,j,k);
	    }
	 }
    }
    for (i=0; i < Nx; i++)
    for (j=0; j < Ny; j++)
    for (k=0; k < Nz; k++)  Result(i,j,k) /= Energ ; 
}

/****************************************************************************/

void make_gaussian2d(fltarray & Result, float Sigma, int Indi, int Indj)
// Create a 2D Gaussian at position  Indi,Indj, with a standard deviation
// equals to Sigma
{
   int Nx = Result.nx();
   int Ny = Result.ny();
   int Dist,Delta_i;
   float sigma2,Energ = 0.;
   register int i,j;
   int Ci=Indi, Cj=Indj;
    if (Ci < 0) Ci = Nx / 2;
    if (Cj < 0) Cj = Ny / 2;
    sigma2 = -2. * Sigma * Sigma;
    for (i=0; i < Nx; i++)
    {
        Delta_i = (i - Ci) * (i - Ci);
        for (j = 0; j < Ny; j++)
	{
	    Dist = Delta_i + (j - Cj) * (j - Cj);
            Result(i,j) = exp((double)((float) Dist / sigma2));
	    Energ +=  Result(i,j);
	 }
    }
    for (i=0; i < Nx; i++)
    for (j=0; j < Ny; j++) Result(i,j) /= Energ ; 
}

/****************************************************************************/

void make_gaussian1d(fltarray & Result, float Sigma, int Indi)
// Create a 1D Gaussian at position  Indi, with a standard deviation
// equals to Sigma
{
   int Nx = Result.nx();
   int Dist;
   float sigma2,Energ = 0.;
   register int i;
   int Ci=Indi;
   if (Ci < 0) Ci = Nx / 2;
   sigma2 = -2. * Sigma * Sigma;
   for (i=0; i < Nx; i++)
   {
        Dist = (i - Ci) * (i - Ci);
        Result(i) = exp((double)((float) Dist / sigma2));
	Energ += Result(i);
   }
   for (i=0; i < Nx; i++) Result(i) /= Energ ; 
}

/****************************************************************************/

void simu_gaussian_noise(fltarray & NoiseData, float Sigma, unsigned int InitRnd)
// Generate a  Gaussian white noise
{
    int i;
    double x;
    if (InitRnd > 0) init_random (InitRnd);
    for (i=0; i < NoiseData.n_elem(); i++)
    {
       x = get_random();
       NoiseData(i) = xerf(x) * Sigma;
    } 
}

/***************************************************************************/

void simu_gaussian_noise3d (fltarray & NoiseData, float Sigma, unsigned int InitRnd)
// Generate a  Gaussian white noise
{
    int i,j,k;
    double x;
    int Nx = NoiseData.nx();
    int Ny = NoiseData.ny();
    int Nz = NoiseData.nz();
        
    if (InitRnd > 0) init_random (InitRnd);
    for (i=0; i < Nx; i++)
    for (j=0; j < Ny; j++)
    for (k=0; k < Nz; k++)
    {
       x = get_random();
       NoiseData(i,j,k) = xerf(x) * Sigma;
    }
}


/***************************************************************************/

void make_gaussian(fltarray & Result, float Sigma)
{
    switch (Result.naxis())
    {
       case 1: make_gaussian1d(Result, Sigma); break;
       case 2: make_gaussian2d(Result, Sigma); break;
       case 3: make_gaussian3d(Result, Sigma); break;
       default:
            cout << "Error in make_gaussian: naxis = " << Result.naxis() << endl;
	    exit(-1);
    }
}

/***************************************************************************/


void ms_vst(int dim, int scale, double &b, double &c)
{    
    if (scale > 10) 
    {
        cout << "Error:VST coefs not precalculated for scale > 10" << endl;
	exit(-1);
    }
    double tt1 = pow(MSVST_T1[scale], dim);
    double tt2 = pow(MSVST_T2[scale], dim);
    double tt3 = pow(MSVST_T3[scale], dim);
    c = 7*tt2/(8*tt1) - tt3 / (2*tt2);
    b = 2 * sqrt(tt1/tt2);
}

/***************************************************************************/

void msvst_transform (fltarray &data, fltarray &vstdata, int scale, Bool Coupled)
{
    double b, c;
    ms_vst(data.naxis(), scale, b, c);
    vstdata.resize(data.nx(), data.ny(), data.nz());
    int n = data.n_elem();
    for (int i=0; i<n; i++)
    {
        if (Coupled == True) vstdata(i) = sqrt(data(i) + c); // tau_1 = 1
        else  vstdata(i) = b * sqrt(data(i) + c);
    }
}

/***************************************************************************/

void msvst_inv_transform(fltarray &vstdata, fltarray &data, int scale, Bool CBBiais)
{
    double b, c;
    ms_vst(vstdata.naxis(), scale, b, c);
    double cb = (CBBiais == True) ? 1. / (b * b) : 0.;
    data.resize(vstdata.nx(), vstdata.ny(), vstdata.nz());
    int n = vstdata.n_elem();
    for (int i=0; i<n; i++)
        data(i) = vstdata(i)*vstdata(i) + cb - c; // tau_1 = 1    
}
		  
/***************************************************************************/

void msvst_transform (Ifloat &data, Ifloat &vstdata, int scale, Bool Coupled)
{
    double b, c;
    ms_vst(2, scale, b, c);
    vstdata.resize(data.nl(), data.nc());
    int n = data.n_elem();
    for (int i=0; i<n; i++)
    {
        if (Coupled == True) vstdata(i) = sqrt(data(i) + c); // tau_1 = 1
        else  vstdata(i) = b * sqrt(data(i) + c);
    }
}

/***************************************************************************/

void msvst_inv_transform(Ifloat &vstdata, Ifloat &data, int scale, Bool CBBiais)
{
    double b, c;
    ms_vst(2, scale, b, c);
    double cb = (CBBiais == True) ? 1. / (b * b) : 0.;

    data.resize(vstdata.nl(), vstdata.nc());
    int n = vstdata.n_elem();
    for (int i=0; i<n; i++)
        data(i) = vstdata(i)*vstdata(i) + cb - c; // tau_1 = 1    
}

/***************************************************************************/
	  
double msvst_var (int dim, int scale)
{
    if (scale > 10)  
    {
        cout << "Error:VST coefs Var not precalculated for scale > 10" << endl;
	exit(-1);
    }
    double tt11 = pow(MSVST_T1[scale-1], dim);
    double tt1  = pow(MSVST_T1[scale], dim);
    double tt21 = pow(MSVST_T2[scale-1], dim);
    double tt2  = pow(MSVST_T2[scale], dim);
    double shh  = pow(MSVST_sumhh[scale], dim);
    
    return (tt21 / (4. * tt11 * tt11) + tt2 / (4. * tt1 * tt1) - shh / (2. * tt11 * tt1));
}

/***************************************************************************/
 

void survival (fltarray & Data, fltarray &Survival, int NpixSurv, float NSigma, Bool Norm, float SigmaNorm, float MeanNorm)
{
   void sort(unsigned long n, double *arr);
   int i,t,N = Data.n_elem();
   dblarray X(N);
   if (Norm == True)
   {
      float Mean = (MeanNorm==-1.) ? Data.mean(): MeanNorm;
      float Sigma = (SigmaNorm==0) ? Data.sigma() : SigmaNorm;
      for (i=0; i < N; i++) X(i) = (Data(i) - Mean) / Sigma;
   } 
   else for (i=0; i < N; i++) X(i) = Data(i);
   
   int Np = NpixSurv;
   float Nsig = NSigma;
   Survival.alloc(Np,2);
   for (i=0; i < Np; i++) Survival(i,0) =  2. * double(i) / double(Np-1) * Nsig - Nsig;
   
   double *Ptr = X.buffer() - 1;
   sort((unsigned long) N, Ptr);
   
   int Cpt=0;
   i=0;
   // for (i=0; i < 10; i++) cout << " " << X(i) << endl;
   
   for (t=0; t < Np/2; t++) 
   {
      double Threshold=Survival(t,0);
      while ( (i < N) && (X(i) < Threshold) && (X(i) < 0))
      {
         Cpt ++;
	 i++;
      }
      Survival(t,1) = Cpt;
   }
   Cpt=0;
   i=N-1;
   for (t=Np-1; t >= Np/2; t--) 
   {
      double Threshold=Survival(t,0);
      while ((i >=0) && (X(i) > Threshold) && (X(i) > 0))
      {
         Cpt ++;
	 i--;
      }
      Survival(t,1) = Cpt;
   }
   // cout << endl <<  " CptP =  " << Cpt << endl;
}

/***************************************************************************/


void get_stat(double *x, int N, dblarray &TabStat)
{
    int NStat = 9;
    TabStat.alloc(NStat);
    double Mean, Sigma, Skew, Curt, Min, Max;
    moment4(x, N,  Mean,  Sigma, Skew,  Curt,  Min,  Max);
    TabStat(0) = Sigma;
    TabStat(1) = Skew;
    TabStat(2) = Curt;
    TabStat(3) = Min;
    TabStat(4) = Max;
    double HC1, HC2;
    hc_test(x, N, HC1, HC2, Sigma, Mean);
    TabStat(5) = HC1;
    TabStat(6) = HC2;
    
    int NbrCumulant=6;
    double *TabCumulant = new double [NbrCumulant];
    
    double *y = new double [N];
    for (int i=0; i < N; i++) y[i] = (x[i] - Mean) / Sigma;
    //   if not keyword_set(norm) then x = x / sigma(x)
    for (int c=0; c <= 5; c++) TabCumulant[c] = cumulants(y, N, c+1);
    TabStat(7) = TabCumulant[4];
    TabStat(8) = TabCumulant[5];
    delete TabCumulant;    
    delete y;
}


/****************************************************************************/

