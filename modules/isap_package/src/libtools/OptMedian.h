/******************************************************************************
**                   Copyright (C) 1998 CEA
*******************************************************************************
**
**    UNIT
**
**    Version: 1.0
**
**    Author: J.L. Starck
**
**    Date:  8/12/98 
**    
**    File:  OptMedian.h
**
*******************************************************************************
**
**    DECRIPTION  Optimized median declaration
**    ---------- 
**
*****************************************************************************/


#ifndef _OPT_MEDIAN_
#define _OPT_MEDIAN_

int opt_med3(int  *p);
int opt_med5(int  *p);
int opt_med7(int  *p);
int opt_med9(int  *p);
int kth_smallest(int a[], int n, int k);
int get_median(int a[], int n);
int abs_kth_smallest(int a[], int n, int k);
int get_abs_median(int a[], int n);

float opt_med3(float  *p);
float opt_med5(float  *p);
float opt_med7(float  *p);
float opt_med9(float  *p);
float kth_smallest(float a[], int n, int k);
float get_median(float a[], int n);
float abs_kth_smallest(float a[], int n, int k);
float get_abs_median(float a[], int n);

int hmedian(int  *ra, int n);
float hmedian(float  *ra, int n);

#endif
