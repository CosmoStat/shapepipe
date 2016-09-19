/******************************************************************************
**
**    UNIT
**
**    Version: 1.0
**
**    Author: Jalal FADILI 
**
**    Date:  20/10/03
**    
**    File:  BesselBayes.cc
**           
******************************************************************************/
 
#include "GlobalInc.h"

#define MAXP 1E1
#define MAXS 2E1

#define EPS 1E-10
#define eps 2.2204e-16
#define pi 3.141592653589793E0

static int SizeDV=0;
static int SizeDP=0;

void dvla(double va,double x,double* pd);

int FIX(double x){return((int)x);}

double SIGN(double x)
{
if(x>=0) return(1);
else return(-1);
}

void vvla(double va,double x,double* pv)
{
/*
%       ===================================================
%       Purpose: Compute parabolic cylinder function Vv(x)
%                for large argument
%       Input:   x  --- Argument
%                va --- Order
%       Output:  PV --- Vv(x)
%       Routines called:
%            (1) DVLA for computing Dv(x) for large |x|
%            (2) GAMMA for computing ג(x)
%       ===================================================
*/


double pdl,gl,Neps,qe,a0,r,x1,dsl;
int k;

pdl=gl=0;
Neps=1.0E-12;
qe=exp(0.25*x*x);
a0=pow(fabs(x),(-va-1.0E0))*sqrt(2.0E0/pi)*qe;
r=1.0E0;
*pv=1.0E0;
for(k=1;k<=18;k++){
	r=0.5E0*r*(2.0*k+va-1.0)*(2.0*k+va)/(k*x*x);
	*pv=(*pv)+r;
	if(fabs(r/(*pv)) < Neps) break;
}
*pv=a0*(*pv);
if(x < 0.0E0){
	x1=-x;
	dvla(va,x1,&pdl);
	gl=exp(lgamma(-va));
	dsl=sin(pi*va)*sin(pi*va);
	*pv=dsl*gl/pi*pdl-cos(pi*va)*(*pv);
}
}


void dvsa(double va,double x,double* pd)
/*
%       ===================================================
%       Purpose: Compute parabolic cylinder function Dv(x)
%                for small argument
%       Input:   x  --- Argument
%                va --- Order
%       Output:  PD --- Dv(x)
%       Routine called: GAMMA for computing ג(x)
%       ===================================================
*/
{

double ga0,g1,g0,gm,Neps,sq2,ep,va0,a0,vt,r,m,vm,r1;

ga0=g1=g0=gm=0;
Neps=1.0E-15;
sq2=sqrt(2.0E0);
ep=exp(-.25E0*x*x);
va0=0.5E0*(1.0E0-va);
if(va == 0.0) *pd=ep;
else{
	if(x == 0.0){
		if(va0 <= 0.0&va0 == FIX(va0)) *pd=0.0E0;
		else{
			ga0=exp(lgamma(va0));
			*pd=sqrt(pi)/(pow(2.0E0,(-.5E0*va))*ga0);
		}
  	}
	else{
		g1=exp(lgamma(-va));
		a0=pow(2.0E0,(-0.5E0*va-1.0E0))*ep/g1;
		vt=-.5E0*va;
		g0=exp(lgamma(vt));
		*pd=g0;
		r=1.0E0;
		for(m=1;m<=250;m++){
		 vm=.5E0*(m-va);
		 gm=exp(lgamma(vm));
		 r=-r*sq2*x/m;
		 r1=gm*r;
		 *pd=*pd+r1;
		 if(fabs(r1) < fabs(*pd)*Neps) {break;}
		}
            	*pd=a0*(*pd);
	}
}
}

void dvla(double va,double x,double* pd)
{
/*
%       ====================================================
%       Purpose: Compute parabolic cylinder functions Dv(x)
%                for large argument
%       Input:   x  --- Argument
%                va --- Order
%       Output:  PD --- Dv(x)
%       Routines called:
%            (1) VVLA for computing Vv(x) for large |x|
%            (2) GAMMA for computing ג(x)
%       ====================================================
*/


double vl,gl,Neps,a0,r,k,ep,x1;

vl=gl=0;
Neps=1.0E-12;
ep=exp(-.25*x*x);
a0=pow(fabs(x),va)*ep;
r=1.0E0;
*pd=1.0E0;
for(k=1;k<=16;k++){
	r=-0.5E0*r*(2.0*k-va-1.0)*(2.0*k-va-2.0)/(k*x*x);
	*pd=*pd+r;
	if(fabs(r/(*pd)) < Neps) break;
}
*pd=a0*(*pd);
if(x < 0.0E0){
x1=-x;
vvla(va,x1,&vl);
gl=exp(lgamma(-va));
*pd=pi*vl/gl+cos(pi*va)*(*pd);
}
}

void err_dv(int N, char *Mes)
{
   if (N >= SizeDV)
   {
      cout << Mes << ", Error in pbdv: ind = " << N << " SizeDV = " << SizeDV << endl;
      exit(-1);
   }
}
void err_dp(int N, char *Mes)
{
   if (N >= SizeDP)
   {
      cout << Mes << "Error in pbdv: ind = " << N << " SizeDP = " << SizeDP << endl;
      exit(-1);
   }
}

void pbdv(double v,double x,double* dv,double* dp,double* pdf,double* pdd)
{
/*
%       ====================================================
%       Purpose: Compute parabolic cylinder functions Dv(x)
%                and their derivatives
%       Input:   x --- Argument of Dv(x)
%                v --- Order of Dv(x)
%       Output:  DV(na) --- Dn+v0(x)
%                DP(na) --- Dn+v0'(x)
%               ( na = |n|, v0 = v-n, |v0| < 1,
%                  n = 0,ס1,ס2,תתת )
%                PDF --- Dv(x)
%                PDD --- Dv'(x)
%       Routines called:
%            (1) DVSA for computing Dv(x) for small |x|
%            (2) DVLA for computing Dv(x) for large |x|
%       ====================================================
*/



double f1,f0,pd0,pd1,xa,v0,ep,v1,vh,pd,v2,f,s0;
int l,k,m,na,ja,nv,nk;

f1=f0=pd0=pd1=0;

xa=fabs(x);
vh=v;
v=v+(fabs(1.0E0)*SIGN(v));
nv=FIX(v);
v0=v-nv;
na=ABS(nv);
ep=exp(-.25E0*x*x);
if(na >= 1) {ja=1;}
if(v >= 0.0){
        if(v0 == 0.0){
		pd0=ep;
		pd1=x*ep;
  	}
	else {
	 for(l=0;l<=ja;l++){
		v1=v0+l;
	 	if(xa <= 5.8) {dvsa(v1,x,&pd1);}
	 	if(xa > 5.8) {dvla(v1,x,&pd1);}
	 	if(l == 0) {pd0=pd1;}
	}
	}
	dv[1]=pd0;
	dv[2]=pd1;
	for(k=2;k<=na;k++) {
		*pdf=x*pd1-(k+v0-1.0E0)*pd0;
		dv[k+1]=*pdf;
		pd0=pd1;
		pd1=*pdf;
	}
 }
else{
	if(x <= 0.0){
		if(xa <= 5.8E0){
		 dvsa(v0,x,&pd0);
		 v1=v0-1.0E0;
		 dvsa(v1,x,&pd1);
   		}
		else{
		 dvla(v0,x,&pd0);
		 v1=v0-1.0E0;
		 dvla(v1,x,&pd1);
		}
		err_dv(2,(char *) "dv2");
		dv[1]=pd0;
		dv[2]=pd1;
		for(k=2;k<=na;k++)
		{
		  pd=(-x*pd1+pd0)/(k-1.0E0-v0);		
		  err_dv(k+1, (char *)"dvk+1");
		  dv[k+1]=pd;
		  pd0=pd1;
		  pd1=pd;
		}
	}
	else{
        	if(x <= 2.0){
			v2=nv+v0;
			if(nv == 0) {v2=v2-1.0E0;}
			nk=FIX(-v2);
			dvsa(v2,x,&f1);
			v1=v2+1.0E0;
			dvsa(v1,x,&f0);		  
			err_dv(nk+1, (char *)"dvnk+1");
                        dv[nk+1]=f1;
                        dv[nk]=f0;
                        for(k=nk-2;k>=0;k--) 
			{
                         f=x*f0+(k-v0+1.0E0)*f1;
			 err_dv(k+1, (char *)"dvk+1-f");
			 dv[k+1]=f;
			 f1=f0;
			 f0=f;
			}
  		}
		else{
			if(xa <= 5.8) dvsa(v0,x,&pd0);
			if(xa > 5.8) dvla(v0,x,&pd0);
			dv[1]=pd0;
			m=100+na;
			f1=0.0E0;
			f0=1.0E-30;
			for(k=m;k>=0;k--)
			{
			 f=x*f0+(k-v0+1.0E0)*f1;
			 err_dv(k+1, (char *)"dvk+1-f2");
			 if(k <= na) dv[k+1]=f;
			 f1=f0;
			 f0=f;
			}
			s0=pd0/f;
			for(k=0;k<=na;k++)
			{
			 err_dv(k+1, (char *)"dvk+1-dv");
			 dv[k+1]=s0*dv[k+1];
			}
		}
	}
	for(k=0;k<=na-1;k++) {
		v1=fabs(v0)+k;
		if(v >= 0.0E0)
			dp[k+1]=0.5E0*x*dv[k+1]-dv[k+1+1];
		else
		dp[k+1]=-0.5E0*x*dv[k+1]-v1*dv[k+1+1];
	}
}
err_dv(na, (char *)"dvna end");
*pdf=dv[na];
*pdd=dp[na];
v=vh;
}

void mpbdv(double v,double x,double* pdf,double* pdd,double *dv,double* dp)
{
/*This program is a direct conversion of the corresponding Fortran program in
%S. Zhang & J. Jin "Computation of Special Functions"(Wiley, 1996).
%       =========================================================
%       Purpose: This program computes the parabolic cylinder
%                functions Dv(x) and their derivatives using
%                subroutine PBDV
%       Input:   x --- Argument of Dv(x)
%                v --- Order of Dv(x)
%       Output:  DV(na) --- Dn+v0(x)
%                DP(na) --- Dn+v0'(x)
%               ( na = |n|, n = int(v), v0 = v-n, |v0| < 1
%                  n = 0,ס1,ס2,תתת, |n| ף 100 )
%                PDF --- Dv(x)
%                PDD --- Dv'(x)
%       Example: v = 5.5,  x =10.0,  v0 = 0.5,  n = 0,1,...,5

%                  n+v0      Dv(x)           Dv'(x)
%                ---------------------------------------
%                  0.5   .43971930D-10  -.21767183D-09
%                  1.5   .43753148D-09  -.21216995D-08
%                  2.5   .43093569D-08  -.20452956D-07
%                  3.5   .41999741D-07  -.19491595D-06
%                  4.5   .40491466D-06  -.18355745D-05
%                  5.5   .38601477D-05  -.17073708D-04
%
%                Dv(x)= .38601477D-05,  Dv'(x)=-.17073708D-04
%       =========================================================*/

double nv,v0,na,vk;
int k;

*dv=0;*dp=0;*pdf=0;*pdd=0;

nv=FIX(v);
v0=v-nv;
na=fabs(nv);
pbdv(v,x,dv,dp,pdf,pdd);

for(k=0;k<=na;k++)
	vk=k*(SIGN(nv))+v0;
}

/*int main(int argc,char* argv[])
{double pdf,pdd,dv[101],dp[101];
mpbdv(atof(argv[1]),atof(argv[2]),&pdf,&pdd,dv,dp);
printf("%g\n",pdf);
}*/





double moment_centered(double data[], int n, int order)
{
	int j;
	double s,ave,mom;

	if (n <= 1) 
	{
	   cout << "Error: n must be at least 2 in moment" << endl;
	   exit(-1);
	}
	s=0.0;
	for (j=1;j<=n;j++) s += data[j];
	ave=s/n;
	if(order == 1){
		mom=ave;
		return(mom);
	}
	mom=0.0;
	for (j=1;j<=n;j++)
		mom += pow(data[j]-ave,(double)order);
	
	mom /= n;
	return(mom);
}

double moment_non_centered(double data[], int n, int order)
{
	int j;
	double s,ave,mom;

	if (n <= 1) 
	{
	   cout << "Error: n must be at least 2 in moment" << endl;
	   exit(-1);
	}
	s=0.0;
	for (j=1;j<=n;j++) s += data[j];
	ave=s/n;
	if(order == 1){
		mom=ave;
		return(mom);
	}
	mom=0.0;
	for (j=1;j<=n;j++)
		mom += pow(data[j],(double)order);
	
	mom /= n;
	return(mom);
}

double cumulant(double* x, int n, int p)
/* Computes the unbiased estimator of the p-th cumulant up to 4th order of vector
 x using k-statistics, except for orders 5 and 6 where biased estimators are used.

*/
{
double mc[7],m[7],kp,nn;
int k;

if(p>6 || p<1)
{ 
	cout << "Error: Order of kn must be <=6 and >=1" << endl;
	exit(-1);
}

nn=(double)n;
mc[1]=moment_centered(x,n,1);
for(k=2;k<=6;k++) 
	mc[k]=moment_centered(x,n,k);

switch(p){
 case 1:
	kp=mc[1];
	break;
 case 2:
	kp=nn*mc[2]/(nn-1);
	break;
 case 3:
	kp=nn*nn*mc[3]/((nn-1)*(nn-2));
	break;
 case 4:
	kp=(nn*nn*((nn+1)*mc[4]-3*(nn-1)*mc[2]*mc[2]))/((nn-1)*(nn-2)*(nn-3));
	break;
 case 5: /* Biased */
	kp=mc[5]-10*mc[2]*mc[3];
	break;
 case 6: /* Biased */
 	m[1]=mc[1];
 	for(k=2;k<=6;k++) 
		m[k]=moment_non_centered(x,n,k);
		
 	kp=-120*pow(m[1],6.)+360*pow(m[1],4)*m[2]-270*m[1]*m[1]*m[2]*m[2]+30*pow(m[2],3)
	   -120*pow(m[1],3.)*m[3]+120*m[1]*m[2]*m[3]-10*m[3]*m[3]+30*m[1]*m[1]*m[4]
	   -15*m[2]*m[4]-6*m[1]*m[5]+m[6];
	break;
}
return(kp);
}


void estim_params_kbesselcont(double* x,int n,double Sigma,double* P,double* C)
{
   double k2,k4;

   k2=cumulant(x,n,2);
   k4=cumulant(x,n,4);

   *P=3*pow(k2-Sigma*Sigma,2.)/k4;
   if(*P<=0) *P=eps;
   *C=(k2-Sigma*Sigma)/(*P);
   if(*C<=0) *C=eps;
}


// ===========================================================================
// ===========================================================================

/********************************************************/
/* 	Parabolic Cylinder Function prototype		*/
/********************************************************/
// void mpbdv(double v, double x, double* pdf, double* pdd, double *dv, double* dp);


/********************************************************/
/* 	Numerator of the PCM Bayesian estimator		*/
/********************************************************/

double numerator(double x, double P, double C)
{
   double pdf1,pdf2,pdd,*dv,*dp;
   double Ret;
   int D=1000;
   int N= (int) P + D;
   if (N < 1)
   {
      cout << " Error: x = " << x << " P = " << P << " C = " << C << " N = " << N << endl;
      exit(-1);
   } 
   dv= new double [N];
   dp= new double [N];
   SizeDV=N;
   SizeDP=N;
   mpbdv(-P-1,x+sqrt(2/C),&pdf1,&pdd,dv,dp);
   mpbdv(-P-1,-x+sqrt(2/C),&pdf2,&pdd,dv,dp);
   Ret = P*(-exp(pow((x+sqrt(2./C)),2)/4.)*pdf1 + exp(pow((-x+sqrt(2./C)),2.)/4.)*pdf2);
   delete [] dv;
   delete [] dp;
   return Ret;
}

/********************************************************/
/* 	Denominator of the PCM Bayesian estimator	*/
/********************************************************/

double denominator(double x, double P, double C)
{
   double pdf1,pdf2,pdd,*dv,*dp;
   double Ret;
   int D=1000;
   int N= (int) P + D;

   if (N < 1)
   {
      cout << " Error: x = " << x << " P = " << P << " C = " << C << " N = " << N << endl;
      exit(-1);
   }    
   SizeDV=N;
   SizeDP=N;
   dv= new double [N];
   dp= new double [N];
   mpbdv(-P,x+sqrt(2/C),&pdf1,&pdd,dv,dp);
   mpbdv(-P,-x+sqrt(2/C),&pdf2,&pdd,dv,dp);
   Ret = exp(pow((x+sqrt(2./C)),2.)/4.)*pdf1 + exp(pow((-x+sqrt(2./C)),2.)/4.)*pdf2;
   delete [] dv;
   delete [] dp;

   return Ret;
}

/********************************************************/
/* 			BKF PCM estimator		*/
/********************************************************/

void estim_bkf_pcm(double* input, double *output, int data_size,double Sigma,double P,double C)
{
   int i;

   // output=(double*)calloc(data_size+1,sizeof(double));

   if(P>EPS && C>EPS && P<MAXP) 
   {
       for(i=1;i<=data_size;i++) 
       {
  	if(fabs(input[i])<MAXS*Sigma)
		output[i] = Sigma*numerator(input[i]/Sigma,P,C/(Sigma*Sigma))/denominator(input[i]/Sigma,P,C/(Sigma*Sigma));
	else
		output[i] = input[i];
       }
    }
    else
    {
      for(i=1;i<=data_size;i++)
  	output[i] = input[i]*(P*C)/(Sigma*Sigma+P*C); /* Gaussian case.*/
    }
}

/********************************************************/
/* 			Main				*/
/********************************************************/

void dbl_bessel_estimate_bayes(double *input, int data_size, double *result, double Sigma, Bool Verbose)
{
   double  P,C; 
   if (Verbose == True) cout << "BESSEL: N = " << data_size << " Sigma = " << Sigma << endl;
   estim_params_kbesselcont(input,data_size,Sigma,&P,&C);
   if (Verbose == True) cout << " P = " << P << " C = " << C << endl;
   estim_bkf_pcm(input,result, data_size, Sigma,P,C);
   if (Verbose == True) cout << "OK" << endl;
}

/*************************************************************************/

void bessel_estimate_bayes(float *Data, int N, float Sigma, Bool Verbose)
{
   double *input = new double[N+1];
   double *result = new double[N+1];
   input[0] = result[0] = 0;
   for (int i=0; i < N; i++) input[i+1] = Data[i];
   dbl_bessel_estimate_bayes(input, N, result, (double) Sigma, Verbose);
//     for (int i=0; i < N; i++) 
//     {
//       if ( ABS(input[i+1]) < 3*Sigma) result[i+1] = 0.;
//       else result[i+1] = input[i+1];
//    }
   for (int i=0; i < N; i++) Data[i] = result[i+1];
   delete [] input;
   delete [] result;
}

/*************************************************************************/

