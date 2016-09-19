
#ifndef _FISTA_H
#define _FISTA_H

#include <fstream>
#include "GlobalInc.h"
#include "FloatTrans.h"
#include "writefits3d.h"

extern Bool Verbose;

class Fista_params
{
public:
	void reset();			// Initialise the class with default parameters
	char NameOut[256];
	int MaxNiter;			// Maximum number of iterations
	float Threshold;		// Thresholding level
	float TolVar;			// Stopping criterium : ||xt+1 - xt||_infty < TolVar
	float Mu;				// Gradient regularization parameter
	bool Fast;				// (F)ISTA
	bool Decreasing;		// if true : linearily decreasing threshold
	bool Positivity;		// The reconstruction must be >0
	bool No_coarse;			// if true : kills the coarse scale before apaplying fista, then thresholds it
};
void Fista_params::reset()
{
	MaxNiter=10;
	Threshold = 1;
	TolVar = 1e-4;
	Mu=1;
	Fast = false;
	Decreasing = false;
	Positivity = false;
	No_coarse = false;
}


class Fista : public Fista_params
{
public:
	Fista(FloatTrans *domain);
	~Fista(){};
	
	FloatTrans* Domain;
	
	void run(fltarray &b, fltarray &z, void (*_degrade)(fltarray &,bool,bool));
	// b : observed signal, in degraded domain
	// z : recovered signal, in degraded domain
	// _degrade : degradation operator. 1bool to apply the degradation, 2bool to transpose
	void run(cfarray &b, cfarray &z, void (*_degrade)(cfarray &,bool,bool));
	// b : observed signal, in degraded domain
	// z : recovered signal, in degraded domain
	// _degrade : degradation operator. 1bool to apply the degradation, 2bool to transpose
	// _degrade(z,*,true) is assumed real (imaginary part ignored)
};


Fista::Fista(FloatTrans *domain)
{
	reset();
	Domain = domain;
}

void Fista::run(fltarray &b, fltarray &z, void (*_degrade)(fltarray &,bool,bool))
{
	char filename[64];
	float *x, *xold, *xtmp; // current estimate, previous estimate, and temp variable

// Other variables
	float tk=1,told;
	float lvl = Threshold;
	int i=0; bool done=false;
	float speed=0; float old_speed=0;
	fltarray y;// previous reconstruction
	bool allocD=true; // Domain allocation
	
//cerr<<"z nx ny nl nc "<<z.nx()<<","<<z.ny()<<","<<z.nl()<<","<<z.nc()<<endl;
//cerr<<"b nx ny nl nc "<<b.nx()<<","<<b.ny()<<","<<b.nl()<<","<<b.nc()<<endl;
// Initial point z=b
	float* CS;
	if(No_coarse) // Substract the coarse scale
	{
		_degrade(b, false, true);// degraded to direct space, no degradation made
		Domain->transform(b, xtmp, allocD); allocD=false;// allocate xtmp and localTB
		Domain->substract_coarse_scale(xtmp,CS,true);
		Domain->recons(xtmp, b);// allocate xtmp and localTB
		_degrade(b, false, false);// direct to degraded space, no degradation made
	}
	z = b; // degraded space
	_degrade(z, false, true); // degraded to direct space, no degradation made
	Domain->transform(z, xtmp, allocD);// allocate xtmp and localTB
	int n=Domain->size_transform();
	
	
	x = new float[n]; for(int k=0;k<n;k++) x[k] = xtmp[k];
	xold = new float[n];

	cerr<<"##########\nBegin ";
	if(Decreasing) cerr<<"MCA";
	else if(Fast) cerr<<"FISTA";
	else cerr<<"ISTA";
	cerr<<" with mu="<<Mu<<", and "<<MaxNiter<<" iterations.\n##########"<<endl;
	
// Fista algorithm. 
	while( i<MaxNiter && (!done || Decreasing) )
	{
//cerr<<"z i nx ny nl nc "<<i<<":"<<z.nx()<<","<<z.ny()<<","<<z.nl()<<","<<z.nc()<<endl;
	// Threshold update
		if(Decreasing) 
		{
			//lvl = Threshold * (Niter-i-1.)/(Niter-1.);
			lvl = 1-i*1.0/float(MaxNiter-1);			// 1..0
			lvl = (pow((double)lvl,5)+lvl/100.)/1.01;
			lvl = Threshold*(1+ lvl*(1000-1));
		}
		
	// Save current solution
		if(i!=0) Domain->recons(xtmp,z);
		if(Positivity) for(int k=0;k<z.n_elem();k++) z(k) = z(k) *(z(k)>0);
		sprintf(filename,"%s_%05d.fits",NameOut,i);writefltarr(filename, z);
		
	// Evaluate the evolution
		y = y-z; // error between the last two estimates
		speed = abs(y.maxfabs()/z.maxfabs());
		done = (speed < TolVar) && (i>0) && (speed<old_speed);
		old_speed=speed;
		if(Verbose) cerr<<" Step "<<i<<", lvl="<<lvl<<", || (z - zt)/z ||_infty = "<<speed<<", TolVar="<<TolVar<<endl;
		y=z; // Save the new solution
		
	// Gradient step
//cerr<<"zx i nx ny nl nc "<<i<<":"<<z.nx()<<","<<z.ny()<<","<<z.nl()<<","<<z.nc()<<endl;
		_degrade(z, true, false);
//cerr<<"zy i nx ny nl nc "<<i<<":"<<z.nx()<<","<<z.ny()<<","<<z.nl()<<","<<z.nc()<<endl;
		z = b - z ;
//cerr<<"zz i nx ny nl nc "<<i<<":"<<z.nx()<<","<<z.ny()<<","<<z.nl()<<","<<z.nc()<<endl;
		_degrade(z, true, true);
		Domain->transform(z, xold);
		for(int k=0;k<n;k++) xtmp[k] = xtmp[k] + Mu*xold[k];
		
	// Save estimate
		if(Fast) for(int k=0;k<n;k++) xold[k] = x[k];
		
	// Proximal operator : soft thresholding
		for(int k=0;k<n;k++) x[k] = xtmp[k];
		Domain->soft_threshold(x,lvl,No_coarse);
		
	// New point
		if(Fast) 
		{
			told=tk;
			tk = (1.+sqrt(1.+4.*told*told))/2.;
			for(int k=0;k<n;k++) xtmp[k] = x[k] + (told-1)/tk * (x[k]-xold[k]);
		}
		else
			for(int k=0;k<n;k++) xtmp[k] = x[k];
		
//cerr<<"z i nx ny nl nc "<<i<<":"<<z.nx()<<","<<z.ny()<<","<<z.nl()<<","<<z.nc()<<endl;
		i++;
	}
	if(No_coarse)
		Domain->add_coarse_scale(x,CS);
	Domain->recons(x,z);
	sprintf(filename,"%s_recons.fits",NameOut);writefltarr(filename, z);
	cerr<<"##########\nEnd.\n##########"<<endl;
}













void Fista::run(cfarray &b, cfarray &z, void (*_degrade)(cfarray&,bool,bool))
{
	char filename[64];
	float *x, *xold, *xtmp; // current estimate, previous estimate, and temp variable
	fltarray zreal;
	
// Other variables
	float tk=1,told;
	float lvl = Threshold;
	int iter=0; bool done=false;
	float speed=0, old_speed=0;
	fltarray yreal;// previous reconstruction
	bool allocD=true; // Domain allocation
	
// Initial point
	float* CS;
	if(No_coarse)
	{
		_degrade(b, false, true);// degraded to direct space, no degradation made
		zreal.alloc(b.nx(),b.ny(),b.nz());
		if(b.nz()) for(int k=0;k<b.nz();k++) for(int j=0;j<b.ny();j++) for(int i=0;i<b.nx();i++) zreal(i,j,k) = b(i,j,k).real();
		else for(int j=0;j<b.ny();j++) for(int i=0;i<b.nx();i++) zreal(i,j) = b(i,j).real();
		Domain->transform(zreal, xtmp, true); allocD=false; // allocate xtmp and localTB
		Domain->substract_coarse_scale(xtmp,CS,true);
		//Domain->add_coarse_scale(xtmp,CS);
		Domain->recons(xtmp, zreal);// allocate xtmp and localTB
		if(b.nz()) for(int k=0;k<b.nz();k++) for(int j=0;j<b.ny();j++) for(int i=0;i<b.nx();i++) b(i,j,k) = complex_f(zreal(i,j,k),0.);
		else for(int j=0;j<b.ny();j++) for(int i=0;i<b.nx();i++) b(i,j) = complex_f(zreal(i,j),0.);
		_degrade(b, false, false);// direct to degraded space, no degradation made
	}
	z = b; // degraded space
	_degrade(z, false, true); // degraded to direct space, no degradation made
	zreal.alloc(z.nx(),z.ny(),z.nz());
	if(z.nz()) for(int k=0;k<z.nz();k++) for(int j=0;j<z.ny();j++) for(int i=0;i<z.nx();i++) zreal(i,j,k) = z(i,j,k).real();
	else for(int j=0;j<z.ny();j++) for(int i=0;i<z.nx();i++) zreal(i,j) = z(i,j).real();
//sprintf(filename,"%sx_%05d.fits",NameOut,iter);writefltarr(filename, zreal);
	Domain->transform(zreal, xtmp, allocD);
	int n=Domain->size_transform();
	x = new float[n]; for(int k=0;k<n;k++) x[k] = xtmp[k];
	xold = new float[n];

	cerr<<"##########\nBegin ";
	if(Decreasing) cerr<<"MCA";
	else if(Fast) cerr<<"FISTA";
	else cerr<<"ISTA";
	cerr<<" with mu="<<Mu<<", and "<<MaxNiter<<" iterations.\n##########"<<endl;

// Fista algorithm. 
	while( iter<MaxNiter && (!done || Decreasing) )
	{
	// Threshold update
		if(Decreasing) 
		{
			//lvl = Threshold * (Niter-i-1.)/(Niter-1.);
			lvl = 1-iter*1.0/float(MaxNiter-1);			// 1..0
			lvl = (pow((double)lvl,3)+lvl/25.)/1.04;
			lvl = Threshold*(1+ lvl*(10-1));
		}
		
	// Save current solution
//cerr<<"zreal nx ny nl nc "<<zreal.nx()<<","<<zreal.ny()<<","<<zreal.nl()<<","<<zreal.nc()<<endl;
		if(iter!=0) Domain->recons(xtmp,zreal);
//cerr<<"zreal nx ny nl nc "<<zreal.nx()<<","<<zreal.ny()<<","<<zreal.nl()<<","<<zreal.nc()<<endl;
		if(Positivity) for(int k=0;k<zreal.n_elem();k++) zreal(k) = zreal(k) *(zreal(k)>0);
		if(z.nz()) for(int k=0;k<z.nz();k++) for(int j=0;j<z.ny();j++) for(int i=0;i<z.nx();i++) z(i,j,k) = complex_f(zreal(i,j,k),0.);
		else for(int j=0;j<z.ny();j++) for(int i=0;i<z.nx();i++) z(i,j) = complex_f(zreal(i,j),0.);
		sprintf(filename,"%s_%05d.fits",NameOut,iter);writefltarr(filename, zreal);
		
	// Evaluate the evolution
		yreal = yreal - zreal; // error between the last two estimates
		speed = abs(yreal.maxfabs());
		done = (speed < TolVar) && (iter>0) && (speed<old_speed);
		old_speed=speed;
		if(Verbose) cerr<<" Step "<<iter<<", lvl="<<lvl<<", || z - zt ||_infty = "<<abs(speed)<<", TolVar="<<TolVar<<endl;
		yreal = zreal; // Save the new solution
		
	// Gradient step
		_degrade(z, true, false);// direct to degraded space
		z = b - z ;// degraded space
		_degrade(z,true,true); // degraded to direct space
		if(z.nz()) for(int k=0;k<z.nz();k++) for(int j=0;j<z.ny();j++) for(int i=0;i<z.nx();i++) zreal(i,j,k) = z(i,j,k).real();
		else for(int j=0;j<z.ny();j++) for(int i=0;i<z.nx();i++) zreal(i,j) = z(i,j).real();
		Domain->transform(zreal, xold);
		for(int k=0;k<n;k++) xtmp[k] = xtmp[k] + Mu*xold[k];
		
	// Save estimate
		if(Fast) for(int k=0;k<n;k++) xold[k] = x[k];
		
	// Proximal operator : soft thresholding
		for(int k=0;k<n;k++) x[k] = xtmp[k];
		Domain->soft_threshold(x,lvl,No_coarse);
		
	// New point
		if(Fast) 
		{
			told=tk;
			tk = (1.+sqrt(1.+4.*told*told))/2.;
			for(int k=0;k<n;k++) xtmp[k] = x[k] + (told-1)/tk * (x[k]-xold[k]);
		}
		else
			for(int k=0;k<n;k++) xtmp[k] = x[k];
		
		iter++;
	}
	if(No_coarse)
		Domain->add_coarse_scale(x,CS);
	Domain->recons(x,zreal);
	if(Positivity) for(int k=0;k<zreal.n_elem();k++) zreal(k) = zreal(k) *(zreal(k)>0);
	sprintf(filename,"%s_recons.fits",NameOut);writefltarr(filename, zreal);
	for(int i=0;i<z.nx();i++) for(int j=0;j<z.ny();j++) for(int k=0;k<z.nz();k++) z(i,j,k) = complex_f(zreal(i,j,k),0.);
	cerr<<"##########\nEnd.\n##########"<<endl;
}


#endif


