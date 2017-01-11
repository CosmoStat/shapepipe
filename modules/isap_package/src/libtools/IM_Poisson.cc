/******************************************************************************
**                   Copyright (C) 1995 by CEA
*******************************************************************************
**
**    UNIT
**
**    Version: 3.1
**
**    Author: Jean-Luc Starck
**
**    Date:  96/05/02 
**    
**    File:  IM_Poisson.cc
**
*******************************************************************************
**
**    DESCRIPTION : multiresolution transform procedures
**    -----------  
**
*******************************************************************************
**
**  void float poidev(float xm, int *idum)
**
**  add a poisson noise to xm. idum is ramdom value.
** 
*******************************************************************************
**
**  void im_poisson_noise(Ifloat &Data)
**  
**  add a poisson noise to an image
** 
******************************************************************************/

// static char sccsid[] = "@(#)IM_Poisson.cc 3.1 96/05/02 CEA 1995 @(#)";

#include "IM_Obj.h"


#define M1 259200
#define IA1 7141
#define IC1 54773
#define RM1 (1.0/M1)
#define M2 134456
#define IA2 8121
#define IC2 28411
#define RM2 (1.0/M2)
#define M3 243000
#define IA3 4561
#define IC3 51349

static float ran1(int *idum);
static float nr_gammln(float xx);

/*********************************************************/

float ran1(int *idum)
{
	static long ix1,ix2,ix3;
	static float r[98];
	float temp;
	static int iff=0;
	int j;

	if (*idum < 0 || iff == 0) {
		iff=1;
		ix1=(IC1-(*idum)) % M1;
		ix1=(IA1*ix1+IC1) % M1;
		ix2=ix1 % M2;
		ix1=(IA1*ix1+IC1) % M1;
		ix3=ix1 % M3;
		for (j=1;j<=97;j++) {
			ix1=(IA1*ix1+IC1) % M1;
			ix2=(IA2*ix2+IC2) % M2;
			r[j]=(ix1+ix2*RM2)*RM1;
		}
		*idum=1;
	}
	ix1=(IA1*ix1+IC1) % M1;
	ix2=(IA2*ix2+IC2) % M2;
	ix3=(IA3*ix3+IC3) % M3;
	j=1 + ((97*ix3)/M3);
	if (j > 97 || j < 1) cout << "ERROR in RAN1: This cannot happen." << endl;
	temp=r[j];
	r[j]=(ix1+ix2*RM2)*RM1;
	return temp;
}

/***********************************************************************/

float nr_gammln(float xx)
{
        double x,y,tmp,ser;
        static double cof[6]={76.18009172947146,-86.50532032941677,
                24.01409824083091,-1.231739572450155,
                0.1208650973866179e-2,-0.5395239384953e-5};
        int j;

        y=x=xx;
        tmp=x+5.5;
        tmp -= (x+0.5)*log(tmp);
        ser=1.000000000190015;
        for (j=0;j<=5;j++) ser += cof[j]/++y;
        return -tmp+log(2.5066282746310005*ser/x);
}


/***********************************************************************/

float poidev(float xm,int *idum)
{
	static float sq,alxm,g,oldm=(-1.0);
	float em,t,y;


	if (xm < 12.0) {
		if (xm != oldm) {
			oldm=xm;
			g=exp(-xm);
		}
		em = -1;
		t=1.0;
		do {
			em += 1.0;
			t *= ran1(idum);
		} while (t > g);
	} else {
		if (xm != oldm) {
			oldm=xm;
			sq=sqrt(2.0*xm);
			alxm=log(xm);
			g=xm*alxm-nr_gammln(xm+1.0);
		}
		do {
			do {
				y=tan(PI*ran1(idum));
				em=sq*y+xm;
			} while (em < 0.0);
			em=floor(em);
			t=0.9*(1.0+y*y)*exp(em*alxm-nr_gammln(em+1.0)-g);
		} while (ran1(idum) > t);
	}
	return em;
}

/***********************************************************************/

void im_poisson_noise(Ifloat &Data)
{
  int Nc,Nl;
  int i,j;
  int number=-1;

  Nc = Data.nc();
  Nl = Data.nl();

  for (i=0;i<Nl;i++)
  for (j=0;j<Nc;j++)
   	   Data(i,j) = poidev(Data(i,j),&number);
     
}

/***********************************************************************/

#undef M1
#undef IA1
#undef IC1
#undef RM1
#undef M2
#undef IA2
#undef IC2
#undef RM2
#undef M3
#undef IA3
#undef IC3



