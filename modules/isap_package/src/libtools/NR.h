/******************************************************************************
**                   Copyright (C) 1996 by CEA
*******************************************************************************
**
**    UNIT
**
**    Version: 2.1
**
**    Author: C. Delattre & R. Gastaud
**
**    Date:  98/05/12
**    
**    File:  NR.h
**
*******************************************************************************
**
**    DESCRIPTION  Definition of the numerical recipies routines
**    -----------  
**              
******************************************************************************/




#ifndef NR_H
#define NR_H

#include <math.h>

#ifndef _UTIL_H
#include "NR_util.h"
#endif

#define EPSILON 1.0e-7

float bessi(int n, float x);
float bessi1(float x);
float bessi0(float x);

void fit(float x[],float y[],int ndata,float sig[],int mwt,float *a,float *b,float *siga,float *sigb,float *chi2,float *q);
float gammq(float a,float x);
float rofunc(float b);
void medfit(float *x,float *y,int ndata,float *a,float *b,float *abdev);
void covsrt(float **covar,int ma,int lista[],int mfit);
void fgauss(float x,float a[],float *y,float dyda[],int na);
void old_mrqcof(float x[],float y[],float sig[],int ndata,float a[],int ma,int lista[],int mfit,float **alpha,float beta[],float *chisq,void (*funcs)(float,float *,float *,float *,int));
void mrqcof(float x[], float y[], float sig[], int ndata, float a[], int ia[],
	int ma, float **alpha, float beta[], float *chisq,
	void (*funcs)(float, float [], float *, float [], int));

void old_mrqmin(float x[],float y[],float sig[],int ndata,float a[],int ma,int lista[],int mfit,float **covar,float **alpha,float *chisq,void (*funcs)(float,float *,float *,float *,int),float *alamda);

void mrqmin(float x[], float y[], float sig[], int ndata, float a[], int ia[],
	int ma, float **covar, float **alpha, float *chisq,
	void (*funcs)(float, float [], float *, float [], int), float *alamda);

float select(unsigned long k, unsigned long n, float arr[]);

float amotry(float **p,float *y,float *psum,int ndim,float ( *funk)(float [],float [],float [],float [],int ,char []),float time[], float data[], float mask[],int size, char *method,int ihi,float fac);
float amotry_overshoot(float **p,float *y,float *psum,int ndim,float ( *funk)(float [],float [],float [],float [],int ,char []),float time[], float data[], float mask[],int size, char *method,int ihi,float fac);

void amoeba(float **p,float y[],int ndim,float ftol,float (*funk)(float [],float [],float [],float [],int ,char []),float time[], float data[],float mask[], int size, char *method,int *nfunk, int NMAX);
void amoeba_overshoot(float **p,float y[],int ndim,float ftol,float (*funk)(float [],float [],float [],float [],int ,char []),float time[], float data[],float mask[], int size, char *method,int *nfunk, int NMAX);

void gcf(float *gammcf,float a,float x,float *gln);
void gser(float *gamser,float a,float x,float *gln);
void sortb(int n,float *ra);
void sort(unsigned long n, float *arr);
void sort(unsigned long n, double *arr);
void sort(unsigned long n, int *arr);
void hpsort(unsigned long n, float *ra);

void gaussj(float **a,int n,float **b,int m);
float gammln(float xx);
float gammp(float a, float x);
float erffc(float x);
float erfcc(float x);
float erff(float x);
float gammln(float xx);
void indexx(unsigned long n, float arr[], unsigned long indx[]);

void eigsrt(float *d, float **v, int n);
/* void eigsrt(float *d, float **v, int n)
**
** Given the eigenvalues d[1..n] and the eigen vectors v[1..][1..n] as output
** from the jacobi function, this routine sorts the eigenvalues into descending
** order, and rearranges the columns of v corresponding.
**
** From numerical recipes in C p.468
*/

void jacobi(float **a, int n, float *d, float **v, int *nrot);
/* void jacobi(float **a, int n, float *d, float **v, int *nrot)
** 
** return the eigenvalues and the eignenvectors of a real symmetric matrix
** 
** a: input real symmetric matrix
** n: input matrix size (a is square)
** d: output vector of eigenvalues
** v: output matrix whose columns contains the normalized  eigenvectors
** nrot: number of Jacobi rotations 
**
** From numerical recipes in C p.467. The procedure is modified to keep 
** the input matrix a unchange. The procedure operate on a copy a1 of matrix a
*/
#endif /* NR_H */
