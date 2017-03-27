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
**    File:  NR_util.h
**
*******************************************************************************
**
**    DESCRIPTION  Definition of the numerical recipies memory  allocation
**    -----------  and deallocation routine.
**              
******************************************************************************/


#ifndef _UTIL_H
#define _UTIL_H

void nrerror(const char error_text[]);
float *vector(int nl,int nh);
int *ivector(int nl,int nh);
double *dvector(int nl,int nh);
float **matrix(int nrl,int nrh,int ncl,int nch);
double **dmatrix(int nrl,int nrh,int ncl,int nch);
int **imatrix(int nrl,int nrh,int ncl,int nch);
float **submatrix(float **a,int oldrl,int oldrh,int oldcl,int oldch,int newrl,int newcl);
void free_vector(float *v,int nl,int nh);
void free_ivector(int *v, int nl,int nh);
void free_dvector(double *v,int nl,int nh);
void free_matrix(float **m,int nrl,int nrh,int ncl,int nch);
void free_dmatrix(double **m,int nrl,int nrh,int ncl,int nch);
void free_imatrix(int **m,int nrl,int nrh,int ncl,int nch);
void free_submatrix(float **b,int nrl,int nrh,int ncl,int nch);
float **convert_matrix(float *a,int nrl,int nrh,int ncl,int nch);
void free_convert_matrix(float **b,int nrl,int nrh,int ncl,int nch);
float ***f3tensor(long nrl, long nrh, long ncl, long nch, long ndl, long ndh);
unsigned long *lvector(long nh,long nl);


#endif /* _UTIL_H */
