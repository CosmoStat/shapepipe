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
**    Date:  3/12/98 
**    
**    File:  DefMath.h
**
*******************************************************************************
**
**    DESCRIPTION  
**    ----------- 
**                 
**    PARAMETRES    
**    ----------    
** 
**    RESULTS      
**    -------  
**
**
******************************************************************************/

#ifndef _DEF_MATH_H_
#define _DEF_MATH_H_

// #include "DefComplex_f.h"
// #include "DefComplex_d.h"

#include "GlobalInc.h"

#define MIN(a,b) (((a) < (b) ? (a):(b)))
#define MAX(a,b) (((a) > (b) ? (a):(b)))

/*
template<class T> inline T MAX(T a, T b)
    { if (a > b) return a; else return b; } 
template<class T> inline T MIN(T a, T b)   
    { if (a > b) return b; else return a; }
*/

// #ifdef WINDOWS
// #define MAXFLOAT 1e20
// #endif

#ifndef MAXFLOAT
#define MAXFLOAT ((float)3.40282346638528860e+38)
#endif

#define FLOAT_EPSILON 5.96047e-08
#define DOUBLE_EPSILON 1.11077e-16
#define Maxfloat MAXFLOAT

#ifndef INFINITY
#define	INFINITY 1.0e+20
#endif

#ifndef PI
// #define PI 3.1415926536 
#define		PI	((double)3.14159265358979323846264338327950288419716939937510)
#endif

#define	ZERO	1.0e-20

inline int iround(float point)
{
  int result;
  if (point >= 0.0) result = (int) (point+0.5);
  else result = (int) (point-0.5);
  return result;
}
inline int iround(double point)
{
   int result;
   if (point >= 0.0) result = (int) (point+0.5);
   else result = (int) (point - 0.5);
   return result;
}

inline int ifloor(float point)
{
   int result;
   if (point >= 0.0) result = (int) (point);
   else result = (int) (point - 1.0);
   return result;
}

inline int ifloor(double point)
{
   int result;
   if (point >= 0.0) result = (int) (point);
   else result = (int) (point - 1.0);
   return result;
}

inline int ABS(int f) {return ( (f < 0) ? -f : f);}
inline short ABS(short arg)  {return (arg < 0)? -arg : arg;}
inline long ABS(long arg) {  return (arg < 0)? -arg : arg;}
inline double ABS(double arg)  {return (arg < 0.0)? -arg : arg;}
inline float ABS(float arg)  {return (arg < 0.0)? -arg : arg;}
/*
#ifdef IABS
inline int abs(int f) {return ( (f < 0) ? -f : f);}
#endif
#ifdef SABS
inline short abs(short arg)  {return (arg < 0)? -arg : arg;}
#endif
#ifdef LABS
inline long abs(long arg) {  return (arg < 0)? -arg : arg;}
#endif
#ifdef DABS
inline double abs(double arg)  {return (arg < 0.0)? -arg : arg;}
#endif
#ifdef FABS
inline float abs(float arg)  {return (arg < 0.0)? -arg : arg;}
#endif
*/
inline int sign(long arg){  return (arg == 0) ? 0 : ( (arg > 0) ? 1 : -1 );}
inline int sign(double arg){return (arg == 0.0) ? 0 : ( (arg > 0.0) ? 1 : -1);}
inline long sqr(long arg){  return arg * arg;}

// psqrt :  sqrt with positive projection
inline double psqrt(double arg) { return (arg <= 0 ? 0 : sqrt(arg)); }
// iilog2 : integer part of the log2. of an integer
inline int iilog2(int arg) { int s = 0; while (arg > 1) {arg >>= 1; s++;} return s;}

inline double sqr(double arg){return arg * arg;}
inline int even(long arg){  return !(arg & 1);}
inline int odd(long arg){return (arg & 1);}
inline void (setbit)(long& x, long b){  x |= (1 << b);}
inline void clearbit(long& x, long b){  x &= ~(1 << b);}
inline int testbit(long x, long b){  return ((x & (1 << b)) != 0);}

inline float POW(float f1, float f2) 
                              {return ( (float) pow (double(f1), double(f2)));}
inline int POW(int f1,int f2) {return (iround(pow (double(f1), double(f2))));}

inline double POW2(double f2) {return pow (double(2.), double(f2));}
inline float POW2(float f2) {return ((float) pow (double(2.), double(f2)));}
inline int POW2(int f2) {return (iround(pow (double(2.), double(f2))));}
inline int IPOW(int x, int y) 
{ int z,l;  
  for (l=0,z=1; l<y; ++l, z*=x);
  return z;
} 

inline double gauss2poisson(float N_Sigma)
{
   double EpsilonPoisson = (1. - erf((double) N_Sigma / sqrt((double) 2.)));
   return EpsilonPoisson;
}
inline double TTgauss2poisson(float N_Sigma)
{
   double EpsilonPoisson = (1. - erf((double) N_Sigma / sqrt((double) 2.)));
   return EpsilonPoisson;
}
/*#define ARG(a,b,Arg) \
   { \
      float Val,Va,Vb;\
      Va = (float) a; Vb = (float) b;\
      if (fabs(Va) < FLOAT_EPSILON) \
      {\
          if (fabs(Vb) < FLOAT_EPSILON)  Arg = 0.; \
          else if (Vb < 0.) Arg = PI / 2.; \
               else Arg =  - PI / 2.; \
      }\
      else \
      {\
          Val = Vb / Va; \
          Arg = atan(Val);\
      }\
   }
*/

#define ARG(a,b,Arg) \
   { \
      double Va,Vb;\
      Va = (double) a; Vb = (double) b;\
      Arg = atan2(Vb,Va);\
   }     
      

/* =============== is power of two ===============================*/

inline Bool is_power_of_2(int  length)
{
   Bool Val;
   int len_exp = (int)(0.3+log((double)(length))/(log(2.0)));
   Val = (length == IPOW(2,len_exp)) ? True: False;
   return Val;
}

#define INT_POW(x,y,z) { int l,xx,yy; xx = (x) ; yy = (y);  for (l=0,(z)=1;l<yy;++ l,z *= xx); }
inline int next_power_of_2(int N) 
{
    int len_exp,temp;

    len_exp = (int)(0.3+log((double)(N))/(log(2.0)));
    INT_POW(2,len_exp,temp);
    if (temp < N) temp *= 2;
    return temp;
}

double xerf (double X);
double xerfc (double X);

/***********************************************************************/

inline float soft_threshold(float Val, float T)
{
   float Coef = Val;
   if (ABS(Coef) < T) Coef = 0.;
   else if (Coef > 0) Coef -= T;
        else Coef += T;
   return Coef;
}

/***********************************************************************/

inline float hard_threshold(float Val, float T)
{
   float Coef = Val;
   if (ABS(Coef) < T) Coef = 0.;
   return Coef;
}

/***********************************************************************/

/* =============== Randonm value ===============================*/
void  init_random (unsigned int Init=100);
float get_random (float Min, float Max);
float get_random();
double b3_spline (double x);
double entropy (float *Data, int Npix, float StepHisto=1.);
float get_sigma_mad(float *Data, int N);
float get_sigma_clip(float *Data, int N, int Nit=3, Bool Average_Non_Null=True, 
                    Bool UseBadPixel=False, float BadPVal=0.);
double skewness(float *Dat, int N);
double curtosis(float *Dat, int N);
void moment4(float *Dat, int N, double &Mean, double &Pow, 
             double &M3, double & M4, float & Min, float & Max, bool not_centered=false);
void moment4(double *Dat, int N, double &Mean, double &Sigma, 
             double &Skew, double & Curt, double & Min, double & Max);
void moment4_center(int N, double x1, double x2, double x3, double x4, double &Sigma, double &Skew, double &Curt);
 
void hc_test(float *Dat, int N, float & HC1, float & HC2, float Sigma, float Mean);

void hc_test(double *Dat, int N, double & HC1, double & HC2, double Sigma, double Mean);

// Higher Criticism Test
void hc_ima_test(float *Dat, float *HCout, int N, float & HC1, float & HC2);

void hc_test(float *Dat, int N, float & HC1, float & HC2, float Mean);
// Higher Criticism Test
// Sigma = MAD(Dat)
void hc_test(float *Dat, int N, float & HC1, float & HC2);
// Higher Criticism Test
// Sigma = MAD(Dat)
// Mean = mean(Dat)
void gausstest(float *Band, int N, float &T1, float &T2);
double cumulants(double *x, int N, int p);


// inline float sqrt(float x) {return (float)sqrt((double)x);}
// inline float log (float x) {return (float)log((double) x);}
// inline float exp (float x) {return (float)exp((double) x);}
// inline float pow (float x, float y) {return (float) pow((double) x, (double) y);}
// inline float pow (double x, float y) {return (float) pow((double) x, (double) y);}
// inline float pow (float x, double y) {return (float) pow((double) x, (double) y);}

// float  fdr_pvalue(float  *TabPValue, int N, float Alpha);
// double fdr_pvalue(double *TabPVal,   int N, float Alpha);
float  fdr_pvalue(float  *TabPValue, int N, double Alpha, bool indep=true);
double fdr_pvalue(double *TabPVal,   int N, double Alpha, bool indep=true);
// return the pvalue detection limit using the False Discovery Rate method

// return the pvalue detection limit using the False Discovery Rate method
// 
// TabPValue = p-values array,
// N = number of pavalues
// Alpha = FDR method parameter

float fdr_gauss_threshold(float *TabVal, int N, float Alpha, float SigmaNoise=1.);
// return the detection threshold value using the False Discovery Rate method
//   TabVal = data
//   N = number of data values
//   Alpha = FDR method parameter
//   SigmaNoise = noise standard deviation

enum filter_type {FT_UNKNOWN, FT_HARD, FT_SOFT, FT_WIENER, FT_FDR, FT_SBT, FT_CONTRAST};// FilterType
char* string_filter_type(filter_type t);

// to solve polynomial equations
int quintic(double [], double [], double [], int*, double);
int quartic(double[], double[], double[], int* );
int cubic(double[], double[], int*);
int signR(double);
double CBRT(double);

#endif
