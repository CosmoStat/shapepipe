/******************************************************************************
**                   Copyright (C) 2005 by CEA
*******************************************************************************
**
**    UNIT
**
**    Version: 1.0
**
**    Author: Laurent Demanet && Jean-Luc Starck
**
**    Date:  21/01/2005 
**    
**    File:  FCur.h
**
**    Modification history :
**
*******************************************************************************
**
**    DESCRIPTION:  Fast Curvelet Transform 
**    ----------- 
******************************************************************************/

#include "FCur.h"
using namespace std;
/***************************************************************************/

FCUR::~FCUR()
{
	for  (int s=0; s < nbr_scale(); s++) 
		delete [] TabCF_Band[s];
	delete [] TabCF_Band;
	delete [] TabSizeNl;
	delete [] TabSizeNc;
}
/***************************************************************************/

static int rangecompute(double XL1, double XL2, int& XS1, int& XS2, int& XF1, int& XF2, double& XR1, double& XR2)
{
	XS1 = 2*int(floor(XL1/2))+1;  XS2 = 2*int(floor(XL2/2))+1; //number of samples
	XF1 = int(floor(XL1/2));  XF2 = int(floor(XL2/2)); //offset on either side
	XR1 = XL1/2;  XR2 = XL2/2;
	// printf("XL1 = %f, XL2 = %f, XS1 = %d, XS2 = %d, XF1 = %d, XF2 = %d, XR1  = %f, XR2 = %f\n",XL1,XL2,XS1,XS2,XF1,XF2,XR1,XR2);
	return 0;
}

/*********************************************************************/

inline int window(double x, double& l, double& r)
{
	if(x<=0) {
		l = 0;	 r = 1;} 
	else if(x>=1) {
		l = 1;	 r = 0;} 
	else
	{
		l = exp(1-1/(1-exp(1-1/(1-x))));
		r = exp(1-1/(1-exp(1-1/x)));
		double norm = sqrt(l*l+r*r);
		l /= norm;
		r /= norm;
	}
	return 0;
}

/*********************************************************************/
// recons exacte sans les modif de fenetre sur wpdata
void static scale_into_wedge(Icomplex_f & Scale, Icomplex_f * & TabBand, int NbrAngles, double XL2, double XL1, bool Undecimated)
{//cerr<<"SIW"<<Undecimated<<endl;
	int nbangles_perquad = NbrAngles / 4;
	int nbquadrants = 4;
	int wcnt = 0;
	int XS1, XS2;  
	int XF1, XF2;  
	double XR1, XR2;

	rangecompute(XL1, XL2, XS1, XS2, XF1, XF2, XR1, XR2);

	for(int q=0; q<nbquadrants; q++) 
	{
		double XW1 = XL1/nbangles_perquad;
		double XW2 = XL2/nbangles_perquad;

		if(1)
		{
			double xs = -XR1;
			double xe = -XR1/4 + (XW1/2)/4; //a quarter of half wedge length
			double ys = -XR2 - (XW2/2);		  
			double ye = -XR2 + 1.5*XW2;
			//how large is the local grid
			int xn,yn;
			if(Undecimated)
			{
				xn = Scale.nc();
				yn = Scale.nl();
			}
			else
			{
				xn = int(ceil(xe-xs));
				yn = int(ceil(ye-ys)); 
			}
			int xf = int(ceil(xs));    int yf = int(ceil(ys)); //where the local grid starts
			double s0 = (-XR2)/(-XR1+XW1/2);
			//double s1 = (-XR2)/(-XR1);
			double s2 = (-XR2+XW2/2)/(-XR1);
			//double s3 = (-XR2+XW2)/(-XR1);
			double s4 = (-XR2+3*XW2/2)/(-XR1);
			double tu = (-XR2 + XW2/2 - (-XR1)) / (-XR2 + XW2/2 + (-XR1));
			double td = (-XR2 - (-XR1+XW1/2)) / (-XR2 + (-XR1+XW1/2));

			cfarray wpdata(xn, yn);
			TabBand[wcnt].resize(yn,xn);
			//Ifloat XXr(yn,xn);
			//Ifloat XXi(yn,xn);

			for(int xcur=xf; xcur<xe; xcur++) {
				int xid = xcur - xf;
				int fm = (int)ceil( max(-XR2, xcur*s0) );
				int to = (int)floor(xcur*s4);
				for(int ycur=fm; ycur<=to; ycur++) {
					int tmp = (ycur-yf) % yn; if(tmp<0) tmp+=yn; assert(tmp>=0);
					if(Undecimated)
					{
						xid = xcur+XF1;
						tmp = ycur+XF2;
					}
					wpdata(xid, tmp) = Scale(ycur+XF2, xcur+XF1);
					//partition of unity
					if(ycur<s2*xcur) { //below s2
						double l,r; window( (double(ycur-xcur)/double(ycur+xcur)-td)/(tu-td), l, r);
						wpdata(xid, tmp) *= l;
					}
					if(ycur>s2*xcur) { //above s2
						double l,r; window( (double(ycur)/double(xcur)-s2)/(s4-s2), l, r);
						wpdata(xid, tmp) *= r;
					}
				}
			}

			// c[sc][wcnt].resize(xn, yn);
			// fftwnd_one(p, (fftw_complex*)wpdata.data(), (fftw_complex*)(c[sc][wcnt].data()));
			// double sqrtprod = sqrt(double(xn*yn));
			for(int i=0; i<xn; i++)
			for(int j=0; j<yn; j++) {
				//c[sc][wcnt](i,j).re /= sqrtprod;				
				//c[sc][wcnt](i,j).im /= sqrtprod;
				(TabBand[wcnt])(j,i) =   wpdata(i,j);
				//XXr(j,i) =  wpdata(i,j).real();
				//XXi(j,i) =  wpdata(i,j).imag() ;
			}
			{
//				cerr<<"WRITE "<<q<<endl;
//				char filename1[64],filename2[64];
//				sprintf(filename1,"xx%d1_re",q);
//				sprintf(filename2,"xx%d1_im",q);
//				io_write_ima_float(filename1, XXr);
//				io_write_ima_float(filename2, XXi);
//				//XXr.info("MXXr");		
			}


			wcnt ++;
		}
		// cout << "REG " << endl;
		//regular
		for(int w=1; w<nbangles_perquad-1; w++)
		{
			double xs = -XR1;		  double xe = -XR1/4;
			double ys = -XR2 + w*XW2-XW2/2;		  
			double ye = -XR2 + (w+1)*XW2+XW2/2;
			int xn,yn;
			if(Undecimated)
			{
				xn = Scale.nc();
				yn = Scale.nl();
			}
			else
			{
				xn = int(ceil(xe-xs));
				yn = int(ceil(ye-ys)); 
			}
			int xf = int(ceil(xs));      int yf = int(ceil(ys));
			double s0 = ys/(-XR1);
			//double s1 = (ys+XW2/2)/(-XR1);
			double s2 = (ys+XW2)/(-XR1);
			//double s3 = (ys+XW2*3/2)/(-XR1);
			double s4 = (ys+XW2*2)/(-XR1);
			cfarray wpdata(xn, yn);
			TabBand[wcnt].resize(yn,xn);
			for(int xcur=xf; xcur<xe; xcur++) {
				int xid = xcur - xf;
				int fm = (int)ceil(xcur*s0);
				int to = (int)floor(xcur*s4);
				for(int ycur=fm; ycur<=to; ycur++) {
					int tmp = (ycur-yf) % yn; if(tmp<0) tmp+=yn; //assert(tmp>=0);
					if(Undecimated)
					{
						xid = xcur+XF1;
						tmp = ycur+XF2;
					}
					wpdata(xid, tmp) = Scale(ycur+XF2, xcur+XF1);
					//partition of unity
					if(ycur<s2*xcur) { //below s2
						double l,r; window( (double(ycur)/double(xcur)-s0)/(s2-s0), l, r);
						wpdata(xid, tmp) *= l;
					}
					if(ycur>s2*xcur) { //above s2
						double l,r; window( (double(ycur)/double(xcur)-s2)/(s4-s2), l, r);
						wpdata(xid, tmp) *= r;
					}
				}
			}

			// double sqrtprod = sqrt(double(xn*yn));
			for(int i=0; i<xn; i++)
			for(int j=0; j<yn; j++) {
				// c[sc][wcnt](i,j).re /= sqrtprod;				
				// c[sc][wcnt](i,j).im /= sqrtprod;
				(TabBand[wcnt])(j,i) =  wpdata(i,j);
			}
			wcnt ++;
		}
		//right
		// cout << "Right " << endl;
		if(1)
		{
			double xs = -XR1;		  double xe = -XR1/4 + (XW1/2)/4; //a quarter of half wedge length
			double ys = XR2 - 1.5*XW2;		  double ye = XR2+(XW2/2);
			int xn,yn;
			if(Undecimated)
			{
				xn = Scale.nc();
				yn = Scale.nl();
			}
			else
			{
				xn = int(ceil(xe-xs));
				yn = int(ceil(ye-ys)); 
			}
			int xf = int(ceil(xs));      int yf = int(ceil(ys));
			double s0 = (XR2-3*XW2/2)/(-XR1);
			//double s1 = (XR2-XW2)/(-XR1);
			double s2 = (XR2-XW2/2)/(-XR1);
			//double s3 = (XR2)/(-XR1);
			double s4 = (XR2)/(-XR1+XW1/2);
			double tu = (XR2+(-XR1+XW1/2))/(XR2-(-XR1+XW1/2));
			double td = (XR2-XW2/2 + (-XR1))/(XR2-XW2/2-(-XR1));

			cfarray wpdata(xn, yn);
			TabBand[wcnt].resize(yn,xn);
			for(int xcur=xf; xcur<xe; xcur++) {
				int xid = xcur - xf;
				int fm = (int)ceil(xcur*s0);
				int to = (int)floor(min(XR2, xcur*s4));
				for(int ycur=fm; ycur<=to; ycur++) {
					int tmp = (ycur-yf) % yn; if(tmp<0) tmp+=yn; //assert(tmp>=0);
					if(Undecimated)
					{
						xid = xcur+XF1;
						tmp = ycur+XF2;
					}
					wpdata(xid, tmp) = Scale(ycur+XF2, xcur+XF1);
					//partition of unity
					if(ycur<s2*xcur) { //below s2
						double l,r; window( (double(ycur)/double(xcur)-s0)/(s2-s0), l, r);
						wpdata(xid, tmp) *= l;
					}
					if(ycur>s2*xcur) {
						double l,r; window( (double(ycur+xcur)/double(ycur-xcur)-td)/(tu-td), l, r);
						wpdata(xid, tmp) *= r;
					}
				}
			}

			// double sqrtprod = sqrt(double(xn*yn));
			for(int i=0; i<xn; i++)
			for(int j=0; j<yn; j++) {
				// c[sc][wcnt](i,j).re /= sqrtprod;				
				// c[sc][wcnt](i,j).im /= sqrtprod;
				(TabBand[wcnt])(j,i) =  wpdata(i,j);
			}
			wcnt ++;
		}
		// cout << "Rot" << XS1 << " " << XS2 << endl;
		//rotate 90 degrees
		// CpxOffMat tmp(XS2, XS1, -XF2, -XF1);
		Icomplex_f TempIma(XS1, XS2);
		for(int i=-XF1; i<-XF1+XS1; i++)
		for(int j=-XF2; j<-XF2+XS2; j++)
			TempIma(i+XF1, -j+XF2) = Scale(j+XF2,i+XF1);
		//io_write_ima_complex_f("xx", Scale);
		//io_write_ima_complex_f("xxr", TempIma);
		Scale.resize(XS1, XS2);
		Scale = TempIma;
		// cout << "EERot" << endl;
		{ double tmp = XL1;    XL1 = XL2;    XL2 = tmp; }
		{ int    tmp = XS1;    XS1 = XS2;    XS2 = tmp; }
		{ int    tmp = XF1;    XF1 = XF2;    XF2 = tmp; }
		{ double tmp = XR1;    XR1 = XR2;    XR2 = tmp; }

	} //for loop quadrant
	 // cout << "END wedge " << endl;
}



void static wedge_into_scale(Icomplex_f * & TabBand, Icomplex_f & Scale, int NbrAngles, double XL2, double XL1, bool Undecimated)
{//cerr<<"WIS"<<Undecimated<<endl;
	int nbangles_perquad = NbrAngles / 4;
	int nbquadrants = 4;
	int wcnt = 0;
	int XS1, XS2;  
	int XF1, XF2;  
	double XR1, XR2;
//Ifloat TB(TabBand->nl(), TabBand->nc());
//for(int i=0;i<TabBand->nl();i++)
//for(int j=0;j<TabBand->nc();j++)
//	TB(i,j) = TabBand[2](i,j).real();

	rangecompute(XL1, XL2, XS1, XS2, XF1, XF2, XR1, XR2);
	Scale.init();
	for(int q=0; q<nbquadrants; q++) {
		double XW1 = XL1/nbangles_perquad;
		double XW2 = XL2/nbangles_perquad;
		//left
		if(1)
		{
			double xs = -XR1;		  
			double xe = -XR1/4 + (XW1/2)/4; //a quarter of half wedge length
			double ys = -XR2 - (XW2/2);		  
			double ye = -XR2 + 1.5*XW2;
			int xn,yn;
			if(Undecimated)
			{
				xn = Scale.nc();
				yn = Scale.nl();
			}
			else
			{
				xn = int(ceil(xe-xs));
				yn = int(ceil(ye-ys)); 
			}
			int xf = int(ceil(xs)); //where the local grid starts   
			int yf = int(ceil(ys)); 
			double s0 = (-XR2)/(-XR1+XW1/2);
			double s2 = (-XR2+XW2/2)/(-XR1);
			double s4 = (-XR2+3*XW2/2)/(-XR1);

			cfarray wpdata(xn, yn);
			for(int i=0; i<xn; i++)
			for(int j=0; j<yn; j++) 
				wpdata(i,j)  = (TabBand[wcnt])(j,i);
			wcnt ++;

			for(int xcur=xf; xcur<xe; xcur++)
			{
				int xid = xcur - xf;
				int fm = (int)ceil( max(-XR2, xcur*s0) );
				int to = (int)floor(xcur*s4);
				for(int ycur=fm; ycur<=to; ycur++) {
					int tmp = (ycur-yf) % yn; if(tmp<0) tmp+=yn; assert(tmp>=0);
					if(Undecimated)
					{
						xid = xcur+XF1;
						tmp = ycur+XF2;
					}
					//partition of unity
					if(ycur<s2*xcur) { //below s2
						double tu = (-XR2 + XW2/2 - (-XR1)) / (-XR2 + XW2/2 + (-XR1));
						double td = (-XR2 - (-XR1+XW1/2)) / (-XR2 + (-XR1+XW1/2));
						double l,r; window( (double(ycur-xcur)/double(ycur+xcur)-td)/(tu-td), l, r);
						wpdata(xid, tmp) *= l;
					}
					if(ycur>s2*xcur) { //above s2
						double l,r; window( (double(ycur)/double(xcur)-s2)/(s4-s2), l, r);
						wpdata(xid, tmp) *= r;
					}
					Scale(ycur+XF2,xcur+XF1) +=  wpdata(xid, tmp);
/*
//					cerr<<"xn yn="<<Undecimated<<xn<<":"<<yn<<endl;
					if(ycur+XF2<0)
						cerr<<"scaleY<0:"<<ycur+XF2<<endl;
					if(xcur+XF1<0)
						cerr<<"scaleX<0:"<<ycur+XF2<<endl;
					if(ycur+XF2>xn)
						cerr<<"scaleY>0:"<<ycur+XF2<<endl;
					if(xcur+XF1>yn)
						cerr<<"scaleX>0:"<<ycur+XF2<<endl;
					if(ycur+XF2>yn)
						cerr<<"scaleY>0:"<<ycur+XF2<<endl;
					if(xcur+XF1>xn)
						cerr<<"scaleX>0:"<<ycur+XF2<<endl;
					Scale(ycur+XF2,xcur+XF1) +=  wpdata(xcur+XF1,ycur+XF2);
*/
				}
			}
		}
		//regular
		for(int w=1; w<nbangles_perquad-1; w++)
		{
			double xs = -XR1;		  double xe = -XR1/4;
			double ys = -XR2 + w*XW2-XW2/2;		  double ye = -XR2 + (w+1)*XW2+XW2/2;
			int xn,yn;
			if(Undecimated)
			{
				xe=XR1;
				ye=XR2;
				xn = Scale.nc();
				yn = Scale.nl();
			}
			else
			{
				xn = int(ceil(xe-xs));
				yn = int(ceil(ye-ys)); 
			}
			int xf = int(ceil(xs));      int yf = int(ceil(ys));
			double s0 = ys/(-XR1);
			//double s1 = (ys+XW2/2)/(-XR1);
			double s2 = (ys+XW2)/(-XR1);
			//double s3 = (ys+XW2*3/2)/(-XR1);
			double s4 = (ys+XW2*2)/(-XR1);

			cfarray wpdata(xn, yn);
			for(int i=0; i<xn; i++)
			for(int j=0; j<yn; j++) 
				wpdata(i,j) = (TabBand[wcnt])(j,i);
			wcnt ++;

			for(int xcur=xf; xcur<xe; xcur++)
			{
				int xid = xcur - xf;
				int fm = (int)ceil(xcur*s0);
				int to = (int)floor(xcur*s4);
				for(int ycur=fm; ycur<=to; ycur++)
				{
					int tmp = (ycur-yf) % yn; if(tmp<0) tmp+=yn; //assert(tmp>=0);
					if(Undecimated)
					{
						xid = xcur+XF1;
						tmp = ycur+XF2;
					}
					//partition of unity
					if(ycur<s2*xcur) { //below s2
						double l,r; window( (double(ycur)/double(xcur)-s0)/(s2-s0), l, r);
						wpdata(xid, tmp) *= l;
					}
					if(ycur>s2*xcur) { //above s2
						double l,r; window( (double(ycur)/double(xcur)-s2)/(s4-s2), l, r);
						wpdata(xid, tmp) *= r;
					}
					Scale(ycur+XF2,xcur+XF1) +=  wpdata(xid, tmp);
				}
			}
		}
		//right
		if(1) {
			double xs = -XR1;		  double xe = -XR1/4 + (XW1/2)/4; //a quarter of half wedge length
			double ys = XR2 - 1.5*XW2;		  double ye = XR2+(XW2/2);
			int xn,yn;
			if(Undecimated)
			{
				xe=XR1;
				ye=XR2;
				xn = Scale.nc();
				yn = Scale.nl();
			}
			else
			{
				xn = int(ceil(xe-xs));
				yn = int(ceil(ye-ys)); 
			}
			int xf = int(ceil(xs));      int yf = int(ceil(ys));
			double s0 = (XR2-3*XW2/2)/(-XR1);
			//double s1 = (XR2-XW2)/(-XR1);
			double s2 = (XR2-XW2/2)/(-XR1);
			//double s3 = (XR2)/(-XR1);
			double s4 = (XR2)/(-XR1+XW1/2);

			cfarray wpdata(xn, yn);
			for(int i=0; i<xn; i++)
			for(int j=0; j<yn; j++) 
				wpdata(i,j)  = (TabBand[wcnt])(j,i);
			wcnt ++;

			for(int xcur=xf; xcur<xe; xcur++) {
				int xid = xcur - xf;
				int fm = (int)ceil(xcur*s0);
				int to = (int)floor(min(XR2, xcur*s4));
				for(int ycur=fm; ycur<=to; ycur++) {
					int tmp = (ycur-yf) % yn; if(tmp<0) tmp+=yn; //assert(tmp>=0);
					if(Undecimated)
					{
						xid = xcur+XF1;
						tmp = ycur+XF2;
					}
					//partition of unity
					if(ycur<s2*xcur) { //below s2
						double l,r; window( (double(ycur)/double(xcur)-s0)/(s2-s0), l, r);
						wpdata(xid, tmp) *= l;
					}
					if(ycur>s2*xcur) {
						double tu = (XR2+(-XR1+XW1/2))/(XR2-(-XR1+XW1/2));
						double td = (XR2-XW2/2 + (-XR1))/(XR2-XW2/2-(-XR1));
						double l,r; window( (double(ycur+xcur)/double(ycur-xcur)-td)/(tu-td), l, r);
						wpdata(xid, tmp) *= r;
					}
					Scale(ycur+XF2,xcur+XF1) +=  wpdata(xid, tmp);
				}
			}
		} //right
		//rotate 90 degrees
		Icomplex_f TempIma(XS1, XS2);
		for(int i=-XF1; i<-XF1+XS1; i++)
		for(int j=-XF2; j<-XF2+XS2; j++)
			TempIma(i+XF1, -j+XF2) = Scale(j+XF2,i+XF1);
		Scale.resize(XS1, XS2);
		Scale = TempIma;

		{ double tmp = XL1;    XL1 = XL2;    XL2 = tmp; }
		{ int    tmp = XS1;    XS1 = XS2;    XS2 = tmp; }
		{ int    tmp = XF1;    XF1 = XF2;    XF2 = tmp; }
		{ double tmp = XR1;    XR1 = XR2;    XR2 = tmp; }
	} //quad
}

/***************************************************************************/

void FCUR::threshold(float SigmaNoise, float N_Sigma)
{
    int i,j;
    Ifloat Band;
    if (Verbose == True) cout << "Noise thresholding: N_Sigma = " << N_Sigma << " SigmaNoise = " << SigmaNoise << endl;
    for (int s=0; s < nbr_scale()-1; s++)
    for (int b=0; b < nbr_band(s); b++) 
    {
		int cnt=0;
        float Nsig = (s==0) ? N_Sigma+1: N_Sigma;
        float Level = Nsig * SigmaNoise * TabCurSigma(s,b);  
        // Level = (s == 0) ? 0.35: 0.25;
	// Level *=  Nsig * SigmaNoise;
        // cout << " Scale " << s+1 <<  " Band " << b+1 << "  N_Sigma = " << Nsig << " Level = " << Level << endl;
        for (i=0; i < size_band_nl(s,b); i++)
        for (j=0; j < size_band_nc(s,b); j++)
        {
           if (ABS((*this)(s,b,i,j))  < Level)
		   {
			   (*this)(s,b,i,j) = 0.;
			   cnt++;
		   }
        }
		float tot=float(size_band_nl(s,b)*size_band_nc(s,b));
	    if (Verbose == True) cout << "Scale ("<<s<<","<<b<<") = thresholded " << cnt << "#:" << (tot-float(cnt))/tot << endl;
    }
} 
 
/****************************************************************************/

void FCUR::wiener_filter(float SigmaNoise, int WienerBlockSize)
{
    int i,j;
    int B2 = WienerBlockSize / 2;
    Ifloat Band;
    if (Verbose == True) cout << "Wiener Noise Filtering ..." << WienerBlockSize <<  " Noise = " <<  SigmaNoise << endl;
    for (int s=0; s < nbr_scale()-1; s++)
    for (int b=0; b < nbr_band(s); b++) 
    {
        float NoiseLevel = SigmaNoise * TabCurSigma(s,b);  
        get_band(s, b, Band);
        float  VarianceSignal=0.,VarianceData=0., VarianceNoise;
        for (i=0; i < size_band_nl(s,b); i++)
        for (j=0; j < size_band_nc(s,b); j++)
        {
		   VarianceData = 0.;
           for (int k=i-B2 ; k <= i+B2 ; k++) 
           for (int l=j-B2 ; l <= j+B2 ; l++) 
		       VarianceData += Band(k,l,I_MIRROR)*Band(k,l,I_MIRROR);
           VarianceData /= (float)(WienerBlockSize*WienerBlockSize); 
	 	   VarianceNoise = NoiseLevel*NoiseLevel;
           VarianceSignal = MAX(0,VarianceData-VarianceNoise);
           Band(i,j) *= VarianceSignal  / VarianceData;
        }
		put_band( s, b,Band);
    }
} 
/****************************************************************************/

void FCUR::get_stat(fltarray &TabStat)
{
    int NbrStatPerBand = 5; // moment of order 2,3,4  
    TabStat.alloc(nbr_tot_band(), NbrStatPerBand);
    Ifloat Band;
    int NumBand = 0;
    for (int s=0; s < nbr_scale(); s++)
    for (int b=0; b < nbr_band(s); b++) 
    {
        int N;
	double Mean, Sigma, Skew, Curt;
	float  Min, Max;
 	get_band(s, b, Band);
	N = Band.nl()*Band.nc();
        moment4(Band.buffer(), N, Mean, Sigma, Skew, Curt, Min, Max);
        TabStat(NumBand, 0) = (float) Sigma;
	TabStat(NumBand, 1) = (float) Skew;
	TabStat(NumBand, 2) = (float) Curt;
        TabStat(NumBand, 3) = (float) Min;
        TabStat(NumBand, 4) = (float) Max; 
        if (Verbose == True)
            printf("  Band %d (%d,%d): Nl = %d, Nc = %d, Sigma = %5.3f, Skew = %5.3f, Curt = %5.3f, Min = %5.3f, Max = %5.3f\n",
	           NumBand+1, s+1, b+1,  Band.nl(), Band.nc(), TabStat(NumBand,0), 
		   TabStat(NumBand, 1),TabStat(NumBand,2),TabStat(NumBand, 3), TabStat(NumBand, 4));
    }  
}
     
/*********************************************************************/
 
void FCUR::get_norm_coeff(Ifloat &ImaNoise, float N_Sigma)
{
   if (Verbose == True)
       cout << "get_norm_coeff " << N_Sigma << endl;
   cur_trans(ImaNoise);
   set_noise_level(N_Sigma);
}

/****************************************************************************/

void FCUR::get_norm_coeff(float N_Sigma)
{
   if (Verbose == True)
       cout << "get_norm_coeff " << N_Sigma << endl;
   float Sig=1.;
   unsigned int InitRnd = 10;

   Ifloat ImSimu(nlima(),ncima(),"ImSimu");
   im_noise_gaussian (ImSimu, Sig, InitRnd);
   cur_trans(ImSimu);   
   set_noise_level(N_Sigma);
}

/****************************************************************************/

void FCUR::set_noise_level(float N_Sigma)
{
    CImaProb CP;
    double LMin, LMax;
    Ifloat Band;

    if (Verbose == True) cout << "Curvelet transform of the noise " << endl;
    for (int s=0; s < nbr_scale(); s++)
    for (int b=0; b <  nbr_band(s); b++) 
    {
        get_band(s, b, Band);
        CP.set(Band);
        CP.find_gthreshold(N_Sigma, LMin, LMax);
        TabCurSigma(s,b) = MAX(ABS(LMin), LMax) /  N_Sigma;
	// TabCurSigma(s,b) = Band.sigma();
        if (Verbose == True) 
             printf("Band(%2d,%2d): Noise Level = %f,  SigBand = %f\n", 
	                 s+1, b+1,  TabCurSigma(s,b),  Band.sigma());
     }
	 isset_tabsigma=true;
}

/*********************************************************************/

void FCUR::import_norm_coef(fltarray *TabSigma)
{
	for (int s=0; s < nbr_scale(); s++)
	for (int b=0; b <  nbr_band(s); b++) 
		TabCurSigma(s,b) = TabSigma[0](s,b);
	isset_tabsigma=true;
}
void FCUR::export_norm_coef(fltarray * &TabSigma)
{
	if(TabSigma==NULL) TabSigma = new fltarray(nbr_scale(),TabNbrBandPerScale.max());
	for (int s=0; s < nbr_scale(); s++)
	for (int b=0; b <  nbr_band(s); b++) 
		TabSigma[0](s,b) = TabCurSigma(s,b);
}

/*********************************************************************/

void FCUR::get_band(int s, int br, Ifloat &Band)
{
   int b = (real() == False) ? br / 2: br;
   Bool GetRe = (br %2 ==0) ? True: False;
   
   int Nls = TabCF_Band[s][b].nl();
   int Ncs = TabCF_Band[s][b].nc();
   Band.resize(Nls,Ncs);
   if ((GetRe == True) || (real() == True))
   {
      for (int i=0;  i < Nls; i++)
      for (int j=0;  j < Ncs; j++) Band(i,j) = (TabCF_Band[s][b])(i,j).real();
   }
   else
   {  
      for (int i=0;  i < Nls; i++)
      for (int j=0;  j < Ncs; j++) Band(i,j) = (TabCF_Band[s][b])(i,j).imag();
   }
}

/*********************************************************************/

void FCUR::put_band(int s, int br, Ifloat &Band)
{
   int b = (real() == False) ? br / 2:  br;
   Bool GetRe = (br %2 ==0) ? True: False;
   
   int Nls = TabCF_Band[s][b].nl();
   int Ncs = TabCF_Band[s][b].nc();
   
   if ((Nls != Band.nl()) || (Ncs != Band.nc()))
   {
      cout << "Error: incorrect image size in FCUR::put_band ... " << endl;
      cout << "       Band size = " << Band.nl() << " " << Band.nc() << endl;
      exit(-1);
   }
   if ((GetRe == True) || (real() == True))
   {
      for (int i=0;  i < Nls; i++)
      for (int j=0;  j < Ncs; j++) 
      {
          complex_f Val = complex_f(Band(i,j), (TabCF_Band[s][b])(i,j).imag());
          (TabCF_Band[s][b])(i,j) = Val;
      }
   }
   else
   {  
      for (int i=0;  i < Nls; i++)
      for (int j=0;  j < Ncs; j++)
      {
 	  complex_f Val = complex_f((TabCF_Band[s][b])(i,j).real(), Band(i,j));
          (TabCF_Band[s][b])(i,j) = Val;
      }
   }
}

/*********************************************************************/

int fct_real_get_band_size(int Nbr_Scale, int Nl, int Nc, int NDir, intarray &TabSizeNl, intarray &TabSizeNc)
{ 
	int NewNl, NewNc, NbrBand;
	if (Nl %2 == 0) NewNl = Nl + 1;
	else NewNl = Nl;
	if (Nc %2 == 0) NewNc = Nc + 1;
	else NewNc = Nc;

	intarray TabDir;
	TabDir.alloc(Nbr_Scale);
	TabDir(0) = NDir;
	int Pas=0;
	for (int s=1; s < Nbr_Scale-1; s++) 
	{
		if (Pas == 0) 
		{
			TabDir(s) = TabDir(s-1);
			Pas = 1;
		}
		else
		{
			TabDir(s) = TabDir(s-1)/2;
			if (TabDir(s) < 8) TabDir(s) = 8;
			Pas = 0;
		}
	}

	double DNL=NewNl;
	double DNC=NewNc;
	intarray TabNl(Nbr_Scale);
	intarray TabNc(Nbr_Scale);
	for (int s=0; s < Nbr_Scale; s++)
	{
		TabNl(s) = 2 * int(floor(DNL/2)) + 1; 
		TabNc(s) = 2 * int(floor(DNC/2)) + 1;  
		DNL /= 2.;
		DNC /= 2.;
	}
	DNL=NewNl;
	DNC=NewNc;

	NbrBand =  0;
	intarray TabNbrBandPerScale(Nbr_Scale);
	for (int s=0; s < Nbr_Scale-1; s++)
	{
		TabNbrBandPerScale(s) = TabDir(s);
		if (TabNbrBandPerScale(s) < 8)  TabNbrBandPerScale(s) = 8;
	}
	TabNbrBandPerScale(Nbr_Scale-1) = 1;
	NbrBand = TabNbrBandPerScale.total();

	TabSizeNl.alloc(NbrBand);
	TabSizeNc.alloc(NbrBand);

	int IndBand=0;
	for  (int s=0; s < Nbr_Scale; s++) 
	{
		int NDir = TabNbrBandPerScale(s);
		int nbangles_perquad = NDir / 4;
		int nbquadrants = 4;
		int wcnt = 0;
		int XS1, XS2;  
		int XF1, XF2;  
		double XR1, XR2;
		double XL1 = DNC;
		double XL2 = DNL;

		rangecompute(XL1, XL2, XS1, XS2, XF1, XF2, XR1, XR2);
		if (NDir <= 1)
		{
			TabSizeNc(IndBand) = TabNc(s);
			TabSizeNl(IndBand++) = TabNl(s);
		}
		else
		{
			for(int q=0; q<nbquadrants; q++) 
			{
				double XW1 = XL1/nbangles_perquad;
				double XW2 = XL2/nbangles_perquad;
				double xs = -XR1;		  
				double xe = -XR1/4 + (XW1/2)/4;  
				double ys = -XR2 - (XW2/2);		  
				double ye = -XR2 + 1.5*XW2;
				int xn = int(ceil(xe-xs));    
				int yn = int(ceil(ye-ys));  
				TabSizeNc(IndBand) = xn;
				TabSizeNl(IndBand++) = yn;
				wcnt ++;
				for(int w=1; w<nbangles_perquad-1; w++) 
				{ 
					xs = -XR1;		  
					xe = -XR1/4;
					ys = -XR2 + w*XW2-XW2/2;		  
					ye = -XR2 + (w+1)*XW2+XW2/2;
					xn = int(ceil(xe-xs));      
					yn = int(ceil(ye-ys));
					TabSizeNc(IndBand) = xn;
					TabSizeNl(IndBand++) = yn;
					wcnt ++;
				}
				xs = -XR1;		  
				xe = -XR1/4 + (XW1/2)/4;  
				ys = XR2 - 1.5*XW2;		  
				ye = XR2+(XW2/2);
				xn = int(ceil(xe-xs));      
				yn = int(ceil(ye-ys));
				TabSizeNc(IndBand) = xn;
				TabSizeNl(IndBand++) = yn;	  
				wcnt ++;
				{ double tmp = XL1;    XL1 = XL2;    XL2 = tmp; }
				{ int    tmp = XS1;    XS1 = XS2;    XS2 = tmp; }
				{ int    tmp = XF1;    XF1 = XF2;    XF2 = tmp; }
				{ double tmp = XR1;    XR1 = XR2;    XR2 = tmp; }
			}
		} 	
		DNL =  DNL/2.;
		DNC =  DNC/2.;
	}
	return NbrBand;
}

/*********************************************************************/

void FCUR::alloc_with_tab(int Nbr_Scale, int Nl, int Nc, intarray & TabDir, Bool ExtendWT, Bool IsotropWT, Bool RealCur)
{
   DataNl = Nl;
   DataNc = Nc;
   if ((Nl %2 == 0) || (Nc %2 == 0)) ModifSize = True;
   else ModifSize = False;
   if (Nl %2 == 0) NewNl = Nl + 1;
   else NewNl = Nl;
   if (Nc %2 == 0) NewNc = Nc + 1;
   else NewNc = Nc;
   
   init(Nbr_Scale,NewNl,NewNc,ExtendWT,IsotropWT,True);
   RealBand = RealCur;
   
   if (TabDir.nx() < Nbr_Scale-1)
   {
       cout << "Error: TabDir size incorrect in FCUR::alloc_with_tab ... " << endl;
       exit(-1);
   }
   TabNbrBandPerScale.alloc(Nbr_Scale);
   TabNbrAnglePerScale.alloc(Nbr_Scale);
   NbrBand = NbrBandCF = 0;
   if (Verbose == True)
    cout << " Number of Scales = " << Nbr_Scale << " " << nbr_scale() << endl;
    
   for (int s=0; s < Nbr_Scale-1; s++)
   {
      TabNbrAnglePerScale(s) = TabDir(s);
      if (TabNbrAnglePerScale(s) < 8) TabNbrAnglePerScale(s) = 8;
      if (real() == False) TabNbrBandPerScale(s) = TabNbrAnglePerScale(s)*2;
      else TabNbrBandPerScale(s) = TabNbrAnglePerScale(s);
      if (Verbose == True) cout << "  Scale " << s+1 << " NbrDir = " << TabNbrAnglePerScale(s) << " NbrBandPerScale = " << TabNbrBandPerScale(s) << endl;
   }
   TabNbrAnglePerScale(Nbr_Scale-1) = 1;
   TabNbrBandPerScale(Nbr_Scale-1) = 1;
   
   TabCF_Band = new Icomplex_f * [nbr_scale()];
   TabSizeNc = new intarray [nbr_scale()];
   TabSizeNl = new intarray [nbr_scale()];
   for  (int s=0; s < nbr_scale(); s++) 
   {
       TabCF_Band[s] = new Icomplex_f [TabNbrBandPerScale(s)];   
       TabSizeNl[s].alloc(TabNbrBandPerScale(s));
       TabSizeNc[s].alloc(TabNbrBandPerScale(s));
   }
   get_size();
   TabCurSigma.alloc(Nbr_Scale,TabNbrBandPerScale.max());
   TabCurSigma.init(1.);
   isset_tabsigma = false;
}

/*********************************************************************/

void FCUR::alloc_from_coarse(int Nbr_Scale, int Nl, int Nc, int NbrDir, Bool ExtendWT, Bool IsotropWT, Bool Real)
{
   intarray TabDir;
   TabDir.alloc(Nbr_Scale);
   TabDir(Nbr_Scale-2) = NbrDir;
   int Pas=0;
   for (int s=Nbr_Scale-3; s >=0; s--) 
   {
       if (Pas == 0) 
       {
           TabDir(s) = TabDir(s+1);
	   Pas = 1;
       }
       else
       {
           TabDir(s) = 2*TabDir(s+1);
	   Pas = 0;
       }
   }
   alloc_with_tab(Nbr_Scale, Nl, Nc, TabDir, ExtendWT, IsotropWT, Real);
}


/*********************************************************************/

void FCUR::alloc_from_fine(int Nbr_Scale, int Nl, int Nc, int NbrDir, Bool ExtendWT, Bool IsotropWT, Bool Real)
{
   intarray TabDir;
   TabDir.alloc(Nbr_Scale);
   TabDir(0) = NbrDir;
   int Pas=0;
   for (int s=1; s < Nbr_Scale-1; s++) 
   {
       if (Pas == 0) 
       {
           TabDir(s) = TabDir(s-1);
	   Pas = 1;
       }
       else
       {
           TabDir(s) = TabDir(s-1)/2;
	   if (TabDir(s) < 8) TabDir(s) = 8;
	   Pas = 0;
       }
   }
   alloc_with_tab(Nbr_Scale, Nl, Nc, TabDir, ExtendWT, IsotropWT, Real);
}

/*********************************************************************/

void FCUR::cur_write(char *Name)
{
    char FileName[256];
    for  (int s=0; s < nbr_scale(); s++) 
    {
       cout << "Scale " << s+1 << " NbrDir = " <<  TabNbrAnglePerScale(s) << endl;
       for (int b=0; b < TabNbrAnglePerScale(s); b++) 
       {
          sprintf(FileName, "%s_%d_%d", Name, s+1, b+1);
          io_write_ima_complex_f(FileName, TabCF_Band[s][b]);
	  Ifloat Temp(TabCF_Band[s][b].nl(), TabCF_Band[s][b].nc());
	  Ifloat TempI(TabCF_Band[s][b].nl(), TabCF_Band[s][b].nc());
	  for (int i=0;  i < Temp.nl(); i++)
	  for (int j=0;  j < Temp.nc(); j++) Temp(i,j) = (TabCF_Band[s][b])(i,j).real();
	  for (int i=0;  i < Temp.nl(); i++)
	  for (int j=0;  j < Temp.nc(); j++) TempI(i,j) = (TabCF_Band[s][b])(i,j).imag();
          cout << "Band  " << b+1 << ": " << FileName  << " Nl = " << TabCF_Band[s][b].nl() << " Nc = " << TabCF_Band[s][b].nc();
	  cout << " SigRe = " << Temp.sigma() << " SigIm = " << TempI.sigma();
	  cout << " EnergRe = " << Temp.energy() << " ErnergIm = " << TempI.energy() << endl;
       }
    }
}

/*********************************************************************/

void FCUR::put_wedges(Icomplex_f * &TabWT)
{
	int b,i,j;
	double DNL =  nld();
	double DNC =  ncd();
	float NormVal;

	for  (int s=0; s < nbr_scale(); s++) 
	{
		int NDir = TabNbrAnglePerScale(s);
		int UseNDir = (real() == True) ? NDir / 2: NDir;
		if (s == nbr_scale()-1) 
			UseNDir=1;
		
		if ((s != (nbr_scale()-1)) && (real() == True))
		{
			for  (b=0; b < UseNDir; b++)
			{
				// cout << s+1 << " "  << b+1 << ": " <<  (TabCF_Band[s][b]).nl() << " " << (TabCF_Band[s][b]).nc() << " ==> " << (TabCF_Band[s][b+NDir/2]).nl() << " " << (TabCF_Band[s][b+NDir/2]).nc() <<  endl;
			// Theoretical norm used : 
				NormVal = sqrt( 0.5/nbr_dir(s)*nls(s)*ncs(s)/(size_band_nl(s,b)*size_band_nc(s,b)) ) * get_norm(s);
				NormVal = 1/NormVal;
				for (i=0; i < (TabCF_Band[s][b]).nl(); i ++)
				for (j=0; j < (TabCF_Band[s][b]).nc(); j ++)
				{
					float Re = (TabCF_Band[s][b])(i,j).real()  / NormVal;
					float Im = (TabCF_Band[s][b+NDir/2])(i,j).real() / NormVal;
					(TabCF_Band[s][b])(i,j) =  complex_f(Re,Im);
					(TabCF_Band[s][b+NDir/2])(i,j) = complex_f(Re,-Im);
				}
			}
		}
		if ((s == (nbr_scale()-1)) && (real() == True))
		{
			b=0;
			NormVal = 1/get_norm(s);
			for (i=0; i < (TabCF_Band[s][b]).nl(); i ++)
			for (j=0; j < (TabCF_Band[s][b]).nc(); j ++)
			{
				float Re = (TabCF_Band[s][b])(i,j).real()  / NormVal;
				// float Im = (TabCF_Band[s][b])(i,j).imag()  / NormVal;
				(TabCF_Band[s][b])(i,j) =  complex_f(Re,0.);
			}
		}

		for  (b=0; b < UseNDir; b++)
		{
			float Norm = sqrt((float) (TabCF_Band[s][b]).nl() * (TabCF_Band[s][b]).nc());
			FFT2D.fftn2d(TabCF_Band[s][b]);
			for (i=0; i < (TabCF_Band[s][b]).nl(); i ++)
			for (j=0; j < (TabCF_Band[s][b]).nc(); j ++) (TabCF_Band[s][b])(i,j) /= Norm;
		}

		if ((s != (nbr_scale()-1)) && (real() == True))
		{
			for  (b=0; b < NDir/2; b++)
			{
				for (i=0; i < (TabCF_Band[s][b]).nl(); i ++)
				for (j=0; j < (TabCF_Band[s][b]).nc(); j ++)
					(TabCF_Band[s][b+NDir/2])(i,j) = conj(  (TabCF_Band[s][b])(i,j));
			}
		}


		TabWT[s].resize(nls(s), ncs(s));
		if (TabNbrAnglePerScale(s) > 1)
		{
			if (Verbose == True)
				cout << "Put Wedges Scale " << s+1 << ", NbrDir = " << TabNbrAnglePerScale(s) << endl;
			wedge_into_scale(TabCF_Band[s], TabWT[s], TabNbrAnglePerScale(s), DNL, DNC, Undecimated);
		}
		else TabWT[s] = TabCF_Band[s][0];
		//cout << "FFT " << endl;
		DNL =  DNL/2.;
		DNC =  DNC/2.;
	}  
    
}

/*********************************************************************/

void FCUR::get_size()
{
     double DNL = nld();
     double DNC = ncd();
//    cerr<<"GS"<<Undecimated<<endl; 
    for  (int s=0; s < nbr_scale(); s++) 
    {
		int NDir = TabNbrAnglePerScale(s);
		int nbangles_perquad = NDir / 4;
		int nbquadrants = 4;
		int wcnt = 0;
		int XS1, XS2;  
		int XF1, XF2;  
		double XR1, XR2;
		double XL1 = DNC;
		double XL2 = DNL;

		rangecompute(XL1, XL2, XS1, XS2, XF1, XF2, XR1, XR2);
		if (NDir <= 1)
		{
			TabSizeNc[s](0) = ncs(s);
			TabSizeNl[s](0) = nls(s);
		}
		else
		{
		if(Undecimated)
		{
			for(int q=0; q<nbquadrants; q++) 
			for(int w=0; w<nbangles_perquad; w++) 
			{
				if(q%2)
				{
					TabSizeNc[s](wcnt) = nls(s);
					TabSizeNl[s](wcnt++) = ncs(s);
				}
				else
				{
					TabSizeNc[s](wcnt) = ncs(s);
					TabSizeNl[s](wcnt++) = nls(s);
				}
//cerr<<q<<"x"<<w<<":"<<ncs(s)<<"x"<<nls(s)<<endl;
			}
		}
		else
		{
			for(int q=0; q<nbquadrants; q++) 
			{
				double XW1 = XL1/nbangles_perquad;
				double XW2 = XL2/nbangles_perquad;
				double xs = -XR1;		  
				double xe = -XR1/4 + (XW1/2)/4;  
				double ys = -XR2 - (XW2/2);		  
				double ye = -XR2 + 1.5*XW2;
				int xn = int(ceil(xe-xs));    
				int yn = int(ceil(ye-ys));  
				TabSizeNc[s](wcnt) = xn;
				TabSizeNl[s](wcnt) = yn;
				wcnt ++;
				for(int w=1; w<nbangles_perquad-1; w++) 
				{ 
					xs = -XR1;		  
					xe = -XR1/4;
					ys = -XR2 + w*XW2-XW2/2;		  
					ye = -XR2 + (w+1)*XW2+XW2/2;
					xn = int(ceil(xe-xs));      
					yn = int(ceil(ye-ys));
					TabSizeNc[s](wcnt) = xn;
					TabSizeNl[s](wcnt) = yn;
					wcnt ++;
				}
				xs = -XR1;		  
				xe = -XR1/4 + (XW1/2)/4;  
				ys = XR2 - 1.5*XW2;		  
				ye = XR2+(XW2/2);
				xn = int(ceil(xe-xs));      
				yn = int(ceil(ye-ys));
				TabSizeNc[s](wcnt) = xn;
				TabSizeNl[s](wcnt) = yn;		  
				wcnt ++;
				{ double tmp = XL1;    XL1 = XL2;    XL2 = tmp; }
				{ int    tmp = XS1;    XS1 = XS2;    XS2 = tmp; }
				{ int    tmp = XF1;    XF1 = XF2;    XF2 = tmp; }
				{ double tmp = XR1;    XR1 = XR2;    XR2 = tmp; }
			}
		}
		} 	
		DNL =  DNL/2.;
		DNC =  DNC/2.;
    }
     
//     for (int s=0; s < nbr_scale(); s++)
//     {
//      cout << " Scale " << s+1 << " nDir " << TabNbrAnglePerScale(s) <<  " NBand = " <<  TabNbrBandPerScale(s) << endl;
//       for (int d=0; d < TabNbrBandPerScale(s); d++)  
//          cout << "   ==> Scale " << s+1 << " Dir " << d+1 << " Nl = " << TabSizeNl[s](d) <<  " Nc = " << TabSizeNc[s](d) << endl;
//     }
}

/*********************************************************************/

void FCUR::get_wedges(Icomplex_f * &TabWT)
{
	int b,i,j;
	double DNL = nld();
	double DNC = ncd();
	float NormVal;
	// FFT2D.CenterZeroFreq=False;
	for  (int s=0; s < nbr_scale(); s++) 
	{
		int NDir = TabNbrAnglePerScale(s);
		int UseNDir = (real() == True) ? NDir / 2: NDir;
		if (s == nbr_scale()-1) 
			UseNDir=1;
		if (NDir > 1)
		{
			if (Verbose == True) cout << "Get Wedges Scale " << s+1 << ", NbrDir = " << TabNbrAnglePerScale(s) <<  endl;
			scale_into_wedge(TabWT[s], TabCF_Band[s], TabNbrAnglePerScale(s), DNL, DNC, Undecimated);
			// cout << "ret scale_into_wedge " << endl;
		}
		else TabCF_Band[s][0] = TabWT[s];
		//cout << "FFT " << endl;

		// FFT2D.CenterZeroFreq = False;
		for  (b=0; b < UseNDir; b++)  
		{
			float Norm = sqrt((float) (TabCF_Band[s][b]).nl() * (TabCF_Band[s][b]).nc());
			for (i=0; i < (TabCF_Band[s][b]).nl(); i ++)
			for (j=0; j < (TabCF_Band[s][b]).nc(); j ++) (TabCF_Band[s][b])(i,j) *= Norm;
			FFT2D.ifftn2d(TabCF_Band[s][b]);

		// Theoretical norm used : 
			NormVal = sqrt( 0.5/nbr_dir(s)*nls(s)*ncs(s)/(size_band_nl(s,b)*size_band_nc(s,b)) ) * get_norm(s);
			NormVal = 1/NormVal;
			if ((s != nbr_scale()-1)  && (real() == True))
			{
				for (i=0; i < (TabCF_Band[s][b]).nl(); i ++)
				for (j=0; j < (TabCF_Band[s][b]).nc(); j ++)
				{
					float Re = (TabCF_Band[s][b])(i,j).real()*NormVal;
					float Im = (TabCF_Band[s][b])(i,j).imag()*NormVal;
					(TabCF_Band[s][b])(i,j) = complex_f(Re, 0.);
					(TabCF_Band[s][b+NDir/2])(i,j) = complex_f(Im, 0.);
				}
			}
			else if ((s == nbr_scale()-1)  && (real() == True))
			{
				NormVal = 1/get_norm(s);
				for (i=0; i < (TabCF_Band[s][0]).nl(); i ++)
				for (j=0; j < (TabCF_Band[s][0]).nc(); j ++)  
				{
					float Re = (TabCF_Band[s][0])(i,j).real()*NormVal;
					// float Im = (TabCF_Band[s][0])(i,j).imag()*NormVal;
					(TabCF_Band[s][0])(i,j) = complex_f(Re, 0.);
				}
			}
		}
		// FFT2D.CenterZeroFreq = True;
		DNL =  DNL/2.;
		DNC =  DNC/2.;
	}
	// FFT2D.CenterZeroFreq=True;
	// cur_write("xxc");
}

/*********************************************************************/

void FCUR::cur_trans(Ifloat &Data)
{
	if (Verbose == True) cout << "Transform WT ... " << endl;
	if (ModifSize == False) transform(Data);
	else
	{
		if (Verbose == True) cout << " NewNl... " << NewNl << endl;
		Ifloat ModIma(NewNl, NewNc, "New");
		for (int i=0; i < Data.nl(); i++)
		for (int j=0; j < Data.nc(); j++) ModIma(i,j) = Data(i,j);

		if (NewNl != DataNl)
			for (int j=0; j < Data.nc(); j++)
				ModIma(DataNl,j) = Data(DataNl-1,j);
		if (NewNc != DataNc)
			for (int i=0; i < Data.nl(); i++)
				ModIma(i,DataNc) = Data(i,DataNc-1);
		if ((NewNl != DataNl) && (NewNc != DataNc))
			ModIma(DataNl,DataNc) = Data(DataNl-1,DataNc-1);
		transform(ModIma);
	}
	if (Verbose == True) cout << "Get wedges ..." << endl;
	get_wedges(Tabcf_WT_Band);
}

/*********************************************************************/

void FCUR::cur_recons(Ifloat &Data)
{
   Data.resize(DataNl, DataNc);
   if (Verbose == True) cout << "Put wedges ..." << endl;
   put_wedges(Tabcf_WT_Band);
   
   if (Verbose == True) cout << "Recons WT ... " << endl;
   if (ModifSize == False) recons(Data);
   else
   {
       Ifloat ModIma(NewNl, NewNc, "New");
       recons(ModIma);
       for (int i=0; i < Data.nl(); i++)
       for (int j=0; j < Data.nc(); j++) Data(i,j) =  ModIma(i,j);
   }
}

/*********************************************************************/

/*************************************************************************/

static void mr_io_name (char *File_Name_In, char *File_Name_Out)
{
	int L;

	strcpy (File_Name_Out, File_Name_In);

	L = strlen (File_Name_In);
	if ((L < 4) || (File_Name_In[L-1] != 't')
				|| (File_Name_In[L-2] != 'c')
				|| (File_Name_In[L-3] != 'f')
				|| (File_Name_In[L-4] != '.'))
		strcat (File_Name_Out, ".fct");
}

/****************************************************************************/
/****************************************************************************/

/*--------------------------------------------------------------------------*/
static void PrintError( int status)
{
    /*****************************************************/
    /* Print out cfitsio error messages and exit program */
    /*****************************************************/

    char status_str[FLEN_STATUS], errmsg[FLEN_ERRMSG];
  
    if (status)
      fprintf(stderr, "\n*** Error occurred during program execution ***\n");

    ffgerr(status, status_str);        /* get the error status description */
    fprintf(stderr, "\nstatus = %d: %s\n", status, status_str);

    if ( ffgmsg(errmsg) )  /* get first message; null if stack is empty */
    {
         fprintf(stderr, "\nError message stack:\n");
         fprintf(stderr, " %s\n", errmsg);

         while ( ffgmsg(errmsg) )  /* get remaining messages */
             fprintf(stderr, " %s\n", errmsg);
    }

    exit( status );       /* terminate the program, returning error status */
}

/*--------------------------------------------------------------------------*/

/*********************************************************************/

void FCUR::mr_io_fill_header(fitsfile *fptr)
{
  int status = 0; // this means OK for cfitsio !!!
    /*****************************************************/
     /* write optional keyword to the header */
    /*****************************************************/
  if ( ffpkyj(fptr, (char*)"Nl", (long) DataNl,(char*)"NlIma",&status))
     PrintError( status );  
  if ( ffpkyj(fptr,(char*)"Nc",(long) DataNc, (char*)"NcIma",&status))
     PrintError( status );  
  if ( ffpkyj(fptr,(char*)"NbrDir",(long) TabNbrAnglePerScale(0),(char*)"NbrDir",&status))
     PrintError( status ); 
  Bool Ext = extend();
  if ( ffpkyj(fptr,(char*)"Extend",(long) ((Ext == True) ? 1: 0),(char*)"Extend WT",&status))
     PrintError( status );
  Bool Iso=isotrop();
  if ( ffpkyj(fptr,(char*)"Isotrop",(long) ((Iso == True) ? 1: 0),(char*)"Isotropic WT",&status))
     PrintError( status );
  if ( ffpkyj(fptr,(char*)"NbrScale",(long) nbr_scale(),(char*)"Number of scales",&status))
     PrintError( status );    
   if ( ffpkyj(fptr,(char*)"Real",(long) ((real() == True) ? 1: 0),(char*)"Real Curv",&status))
     PrintError( status );  
} 

/*********************************************************************/

void FCUR::write(char *Name)
{
 char filename[256];
 fitsfile *fptr;    
 int status;
 //int i,j,s;
 //float *Ptr;
 int simple;
 int bitpix;
 long naxis=0;
 long naxes[3];
 long nelements;
 long group = 1;  /* group to write in the fits file, 1= first group */
 long firstpixel = 1;    /* first pixel to write (begin with 1) */
 Ifloat Aux;
 //long fpixels[3];
 //long lpixels[3];

/* we keep mr as extension even if its fits ! */
 mr_io_name (Name, filename);
 // inc[0]=1;  inc[1]=1; inc[2]=1;

#if DEBUG_IO  
    cout << "Write on " << filename << endl;
#endif

 FILE *FEXIST = fopen(filename, "rb");
 if (FEXIST)
 {
    fclose(FEXIST);
    remove(filename);               /* Delete old file if it already exists */
 }
 
 
  status = 0;         /* initialize status before calling fitsio routines */
   /* open the file */
 if ( ffinit(&fptr, filename, &status) )     /* create the new FITS file */
     PrintError( status );           /* call PrintError if error occurs */
                                                                              
 /* write  the header */
 simple   = True;
 bitpix   =  -32;   /* 32-bit real pixel values      */
 long pcount   =   0;  /* no group parameters */
 long gcount   =   1;  /* only a single image/group */
 int  extend   =   False;
 naxis = 1;


  int Nelem=1+nbr_scale();
  for (int s=0; s < nbr_scale(); s++) 
  for (int b=0; b < TabNbrAnglePerScale(s); b++) 
  {
     if ((s == nbr_scale()-1) || (real() == False)) Nelem += 2 +  TabCF_Band[s][b].nl()*TabCF_Band[s][b].nc()*2;
     else Nelem += 2 +  TabCF_Band[s][b].nl()*TabCF_Band[s][b].nc();
  } 
  
  fltarray Data(Nelem);
  int ind=0;
  Data(ind++) = nbr_scale();
  for  (int s=0; s < nbr_scale(); s++) Data(ind++) = TabNbrAnglePerScale(s);

  // Ifloat Band;
  for (int s=0; s < nbr_scale(); s++) 
  for (int b=0; b < TabNbrAnglePerScale(s); b++)
  {
     Data(ind++) = TabCF_Band[s][b].nl();
     Data(ind++) = TabCF_Band[s][b].nc();
    //  cout << " s = " << s+1 << " b = " << b+1 << " Nl = " << TabCF_Band[s][b].nl() << " Nc = " << TabCF_Band[s][b].nc() << endl;
    //get_band( s,  b,  Band);
    //cout << "   sig = " << Band.sigma() << endl;
     
     for (int i=0;  i < TabCF_Band[s][b].nl(); i++)
     for (int j=0;  j < TabCF_Band[s][b].nc(); j++) Data(ind++) = (TabCF_Band[s][b])(i,j).real();
     // if ((s == nbr_scale()-1) || (real() == False))
     if (real() == False)
     {
        for (int i=0;  i < TabCF_Band[s][b].nl(); i++)
        for (int j=0;  j < TabCF_Band[s][b].nc(); j++) Data(ind++) = (TabCF_Band[s][b])(i,j).imag();
     }
  }
  
  naxes[0] = ind;
  
 // write first header part (parameters)
 if ( ffphpr(fptr,simple,bitpix,naxis,naxes,pcount,gcount,extend,&status) )
     PrintError( status );          /* call PrintError if error occurs */
  
   // write the header of the multiresolution file
  mr_io_fill_header(fptr);
 
 nelements = ind;
 if ( ffppre(fptr, group, firstpixel, nelements, Data.buffer(), &status) )
              PrintError( status );  
     
  /* close the FITS file */
  if ( ffclos(fptr, &status) )  PrintError( status ); 
 }

/*********************************************************************/

void FCUR::read (char *Name)
{
    char filename[256];
    fitsfile *fptr;           /* pointer to the FITS file */
    int status=0, hdutype ;
    long hdunum;
    char comment[FLEN_COMMENT]; 
    int naxis;
    long naxes[3];
    long mon_long;
    int anynul = 0;
    long nulval = 0;
    long inc[3];
    void PrintError( int status);
    long nelements = 0 ; // naxes[0] * naxes[1] in the image
    //long fpixels[3];
    //long int lpixels[3];
    
     // for multiresol
    float *Ptr;
    //int my_logical; // sais pas...

     mr_io_name (Name, filename);

    inc[0]=1;  inc[1]=1; inc[2]=1;
 
#if DEBUG_IO
    cout << "Read in " << filename << endl;
#endif
   
    /* open the file */
    status = 0;         /* initialize status before calling fitsio routines */
    if ( ffopen(&fptr, filename, (int) READONLY, &status) ) 
         PrintError( status );
                                    
    hdunum = 1;  /*read  table */
    if ( ffmahd(fptr, hdunum, &hdutype, &status) ) /* move to the HDU */
           PrintError( status );

    int simple, bitpix, extend;
    long pcount, gcount;
    if ( ffghpr(fptr, 3, &simple, &bitpix, &naxis, naxes, &pcount,
            &gcount, &extend, &status)) /* move to the HDU */
           PrintError( status );

    nelements = naxes[0];
     
    if (ffgkyj(fptr,(char*)"Nl", &mon_long, comment, &status)) PrintError( status );
    int Nli = (int) mon_long;  
    if (ffgkyj(fptr,(char*)"Nc", &mon_long, comment, &status)) PrintError( status );
    int Nci = (int) mon_long; 
    if (ffgkyj(fptr,(char*)"NbrDir", &mon_long, comment, &status)) PrintError( status );
    int Ndir  = (int) mon_long;
    if (ffgkyj(fptr,(char*)"NbrScale", &mon_long, comment, &status)) PrintError( status );
    int NbrScale2D = (int) mon_long;
    if ( ffgkyj(fptr,(char*)"Extend",&mon_long, comment, &status))
    PrintError( status );    
    Bool Ext  = (mon_long == 0) ? False : True;
    if ( ffgkyj(fptr,(char*)"Isotrop",&mon_long, comment, &status))
    PrintError( status );
    Bool Iso  = (mon_long == 0) ? False : True;
    if ( ffgkyj(fptr,(char*)"Real",&mon_long, comment, &status))
    PrintError( status );    
    Bool Real  = (mon_long == 0) ? False : True;
     //cout << "        NbrScale2D = " << NbrScale2D << " Nli = " << Nli <<  " Nci = " << Nci <<   " Ndir = " << Ndir << endl;
    alloc_from_fine(NbrScale2D, Nli,  Nci, Ndir, Ext, Iso, Real);
   //cout << "        Nl = " << (TabCF_Band[0][1]).nl() << " Nc = " << TabCF_Band[0][1].nc() << endl;
    fltarray Tab(nelements);
   //cout << "        Nl = " << (TabCF_Band[0][1]).nl() << " Nc = " << TabCF_Band[0][1].nc() << endl;

    Ptr = Tab.buffer();
    if ( ffgpve(fptr, 1, 1, nelements, nulval, Ptr, &anynul, &status))
             PrintError( status );

    int ind = nbr_scale()+1;
    // Ifloat Band;
    for (int s=0; s < nbr_scale(); s++) 
    for (int b=0; b < TabNbrAnglePerScale(s); b++)
    {
       int Ny = (int) Tab(ind++);
       int Nx = (int) Tab(ind++);
       TabCF_Band[s][b].alloc(Ny ,Nx);
       //cout << " s = " << s+1 << " b = " << b+1 << " Nl = " << Ny << " Nc = " << Nx << endl;
       //cout << "        Nl = " << (TabCF_Band[s][b]).nl() << " Nc = " << TabCF_Band[s][b].nc() << endl;
       for (int i=0;  i < Ny; i++)
       for (int j=0;  j < Nx; j++)
            (TabCF_Band[s][b])(i,j) = complex_f(Tab(ind++), 0.);
       // if ((s == nbr_scale()-1) || (real() == False))
       if (real() == False)
       {
          for (int i=0;  i < Ny; i++)
          for (int j=0;  j < Nx; j++) 
                (TabCF_Band[s][b])(i,j) = complex_f( (TabCF_Band[s][b])(i,j).real(), Tab(ind++));
       }
       // get_band( s,  b,  Band);
       //cout << "   sig = " << Band.sigma() << endl;
   }
      
   if ( ffclos(fptr, &status) ) PrintError( status );
// cout << " end of read fits file " << endl;
}


/*********************************************************************/
