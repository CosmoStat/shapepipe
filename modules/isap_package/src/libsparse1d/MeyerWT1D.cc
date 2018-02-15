/******************************************************************************
**                   Copyright (C) 2008 by CEA
*******************************************************************************
**
**    UNIT
**
**    Version: 1.0
**
**    Author: Arnaud Woiselle
**
**    Date:  30/09/2008
**    
**    File:  MeyerWT.cc
**
**    Modification history :
**
*******************************************************************************
**
**    DESCRIPTION:  1D MEYER WAVELET TRANSFORM 
**    ----------- 
******************************************************************************/
 

#include "MeyerWT1D.h"
//extern bool Verbose;

using namespace std;

static float Tab_1D_Meyer[10] = 
	{	0.790166, 0.611871, 0.613524, 0.613336, 0.611089, 0.6, 0.6, 0.6, 0.6, 0.6 };

/*********************************************************************/
static inline void PrintError( int status)
{
    // ***************************************************** 
    // * Print out cfitsio error messages and exit program * 
    // ***************************************************** 

    char status_str[FLEN_STATUS], errmsg[FLEN_ERRMSG];

    if (status)
      fprintf(stderr, "\n*** Error occurred during program execution ***\n");

    ffgerr(status, status_str);        // get the error status description 
    fprintf(stderr, "\nstatus = %d: %s\n", status, status_str);

    if ( ffgmsg(errmsg) )  // get first message; null if stack is empty 
    {
         fprintf(stderr, "\nError message stack:\n");
         fprintf(stderr, " %s\n", errmsg);

         while ( ffgmsg(errmsg) )  // get remaining messages 
             fprintf(stderr, " %s\n", errmsg);
    }

    exit( status );
}

/****************************************************************************/

static inline void mr_io_name (char *File_Name_In, char *File_Name_Out)
{
    int L;

    strcpy (File_Name_Out, File_Name_In);

    L = strlen (File_Name_In);
    if ((L < 3) || (File_Name_In[L-1] != 'r')
                || (File_Name_In[L-2] != 'm')
                || (File_Name_In[L-3] != '.'))
    {
        strcat (File_Name_Out, ".mr");
    }
}

/****************************************************************************/

MEYER_WT1D::MEYER_WT1D()
{
	FFT1D.CenterZeroFreq = True;
	NbrScale=0;
	Tabcf_WT_Band=NULL;
    Verbose = false;
}

MEYER_WT1D::~MEYER_WT1D()
{
	if (Tabcf_WT_Band !=NULL) delete [] Tabcf_WT_Band;
}

void MEYER_WT1D::get_hfilter(fltarray &H, double DNx, double DNy, double DNz)
{
//	cerr<<"get_hfilter"<<endl;
	int i;
	int Nx = H.nx();

//	cerr<<" Nx Nx2 "<<Nx<<" "<<Ny<<" "<<Nz<<"   "<<Nx2<<" "<<Ny2<<" "<<Nz2<<endl;
	H.init(1.);
	for(i=0; i < Nx; i++)
	{
		int N = Nx/4;
		double x = (i - Nx/2);
		double r = (ABS(x) - N) / N;
		if (r <= 0) H(i) = 1.;
		else if (r < 1) H(i) = lowpass_window(1. - r);
		else H(i) = 0.;
	}
	// io_write_cube_float("low_pass.fits", H);
	// exit(0);
//	cerr<<"end get_hfilter"<<endl;
}

/*********************************************************************/

void MEYER_WT1D::init(int Nbr_Scale, int Nx, Bool WTNeedOddSize)
{
	bool LocVerbose = false && Verbose;
	if(Verbose) cerr<<"MEYER_WT1D::init("<<Nbr_Scale<<","<<Nx<<","<<WTNeedOddSize<<")"<<endl;

	NbrScale = Nbr_Scale;
	Nx_Cube = Nx;
	TabNx.alloc(NbrScale);
	NeedOddSize = WTNeedOddSize;

	if (LocVerbose)
	{
		cout << " INIT WT: " << "CubeSize = " << Nx << " NbrScale = " << NbrScale << endl;
		cout << "   Use Meyer's wavelets without cube extension  " << endl;
	}

	D_ExtNx = (float) Nx;
	if (LocVerbose) cerr<<" Nx "<<Nx<<endl;
	if (LocVerbose) cerr<<" D_ExtNx "<<D_ExtNx<<endl;
	
	double DNX=D_ExtNx;
	
	for (int s=0; s < Nbr_Scale; s++)
	{
		if ((NeedOddSize == True))
			TabNx(s) = 2 * int(floor(DNX/2)) + 1; 
		else 
			TabNx(s) = int( DNX + 0.5);
		DNX /= 2.;
		if (LocVerbose) cerr<<" TabNx("<<s<<") "<<TabNx(s)<<endl;
	}
	ExtNx = TabNx(0);
	
// Memory allocation
	Tabcf_WT_Band = new cfarray[NbrScale];
	for  (int s=0; s < NbrScale; s++)
		Tabcf_WT_Band[s].alloc(TabNx(s));
	TF_ExtData.alloc(ExtNx);
	H.alloc(ExtNx);
	
	if(Verbose) cerr<<"...end MEYER_WT1D::init"<<endl;
}

/*********************************************************************/

void MEYER_WT1D::get_extFourier(fltarray &Data, cfarray &TF_ExtData)
{
	if(Verbose) cerr<<"WT::get_extFourier..."<<endl;
	
	int i;
	int Nx = Data.nx();

	TF_ExtData.resize(ExtNx);
	FFT1D.fftn1d(Data, TF_ExtData, False);

	float  Norm = sqrt(float(Nx));
	for(i=0; i<Nx; i++)
		TF_ExtData(i) /= Norm;

	if(Verbose) cerr<<"...End WT::get_extFourier"<<endl;
}

/*********************************************************************/

void MEYER_WT1D::get_IFFT_Cube(cfarray & TF_ExtData, fltarray &Data)
{
	if(Verbose) cerr<<"Get_IFFT_Cube..."<<endl;
	int i;
	int Nx = Data.nx();

	if (ExtNx != TF_ExtData.nx())
	{
		cout << "Error: bad cube size in get_IFFT_Cube ... " << endl;
		cout << "ExtNx = " << ExtNx << endl;
		cout << "InExtNx = " << TF_ExtData.nx() << endl;
		exit(-1);
	}

	float  Norm = sqrt(float(Nx));
	for(i=0; i<Nx; i++)
		TF_ExtData(i) *= Norm;
	FFT1D.fftn1d(TF_ExtData,True);
	for(i=0; i<Nx; i++)
		Data(i) = TF_ExtData(i).real();    

	if(Verbose) cerr<<"...End Get_IFFT_Cube"<<endl;
}

/*********************************************************************/

void MEYER_WT1D::transform_cf(cfarray & TF_ExtData,  cfarray * &TabWT)
{
	if(Verbose) cerr<<"MEYER_WT1D::transform_cf..."<<endl;
	int i,s,Nxs;
	int Nx = TF_ExtData.nx();
	double DNX = D_ExtNx;

	TabWT[0] = TF_ExtData;
	for (s=0; s < NbrScale-1; s++) 
	{
		Nx = TabNx(s);
		Nxs = TabNx(s+1);
		H.resize(Nxs);
		TF_ExtData.resize(Nxs);
		get_hfilter(H);
		for (i=0; i < Nxs; i++)
		{
			int Indi = i - Nxs/2 + Nx/2;
			TF_ExtData(i) = H(i) * (TabWT[s])(Indi);
			(TabWT[s])(Indi) *= sqrt(1-H(i)*H(i));
		}
		TabWT[s+1].resize(Nx);
		TabWT[s+1] = TF_ExtData; 
		DNX =  DNX/2.;
	}
	if(Verbose) cerr<<"...End MEYER_WT1D::transform_cf"<<endl;
}
/*********************************************************************/

void MEYER_WT1D::recons_cf(cfarray * &TabWT,  cfarray & TF_ExtData)
{
	bool LocVerbose = False && Verbose;
	if(Verbose) cerr<<"MEYER_WT1D::recons_cf..."<<endl;

	int i,s;
	int Nx = TabWT[0].nx();
	fltarray H(Nx);
	fltarray HB(Nx);

	TF_ExtData.alloc(TabWT[0].nx());

	s = NbrScale-1;
	Nx = TabNx(s);
	TF_ExtData.resize(Nx);
	TF_ExtData = TabWT[s];
	for (s=NbrScale-2; s >= 0; s--)
	{
		if (LocVerbose == True) cout << "Rec WT Scale " << s+1 << " " << TabNx(s) << endl;
		Nx = TabNx(s);
		int Nxs = TabNx(s+1);
		H.resize(Nxs);
		get_hfilter(H);
		for (i=0; i < Nxs; i++)
		{
			int Indi = i - Nxs/2 + Nx/2;
			(TabWT[s])(Indi) *= sqrt(1-H(i)*H(i));
			TF_ExtData(i) *= H(i);
			(TabWT[s])(Indi) += TF_ExtData(i);
		}   
		TF_ExtData.resize(Nx);
		TF_ExtData = TabWT[s];
	}
	if(Verbose) cerr<<"...End MEYER_WT1D::recons_cf"<<endl;
}
 
/*********************************************************************/

void MEYER_WT1D::ifft_tabcube(cfarray * & TabCF_Cube, fltarray * & Tab_Cube, Bool Alloc)
{
	if(Verbose) cerr<<"MEYER_WT1D::ifft_tabcube, alloc="<<(bool)Alloc<<endl;
    if (Alloc == True) Tab_Cube = new fltarray[NbrScale];
    for  (int s=0; s < NbrScale; s++)
    {
       float Norm =  sqrt((float)(TabNx(s)));
       if (Alloc == True) Tab_Cube[s].alloc(TabNx(s));
       for (int i=0; i < TabNx(s); i++) (TabCF_Cube[s])(i) *= Norm;
       FFT1D.fftn1d(TabCF_Cube[s], True);
       for (int i=0; i < TabNx(s); i++) (Tab_Cube[s])(i) = (TabCF_Cube[s])(i).real();
    }
	if(Verbose) cerr<<"...End MEYER_WT1D::ifft_tabcube"<<endl;
}

/*********************************************************************/

void MEYER_WT1D::fft_tabcube(fltarray * & Tab_Cube, cfarray * & TabCF_Cube)
{
    // if (Alloc == True)  TabCF_Cube = new cfarray[NbrScale];
    for  (int s=0; s < NbrScale; s++)
    {
       float Norm = sqrt( (float) TabNx(s));
       FFT1D.fftn1d(Tab_Cube[s], TabCF_Cube[s]);
       for (int i=0; i < TabNx(s); i++) (TabCF_Cube[s])(i) /= Norm;
    }
}

/*********************************************************************/

// The same as extract_stat except that there is no normalisation 
void MEYER_WT1D::noise_calibration(fltarray *TabBand, char* Outname)
{
	bool LocVerbose = true & Verbose;
	if(Verbose) cerr<<"Noise_calibration..."<<endl;
	
	TabStat = new double[NbrScale];
	
// Output stat file
	char Statname[250];
	strcpy(Statname, Outname);
	strcat(Statname, "_noise_calib.dat");
	// fstream cstat;
	// cstat.open (Statname, fstream::out);
	
	double mean;
	
	for(int s=0;s<NbrScale-1;s++)
	{
		mean=(TabBand[s]).mean();
		TabStat[s]=(TabBand[s]).sigma();
		// if(LocVerbose)
		// 	cerr<< "\tStat echelle "<<s<<" : mean,sigma = "<<mean<<" "<<TabStat[s]<<endl;
		// cstat <<mean<<"\t"<<TabStat[s]<<"\t"<<endl;
	}
	
// coarse scale : N-1
	// No calibration for the coarsest scale : no "/TabSigma(N-1,0)" defined, only stats

	// cstat.close();
	if(Verbose) cerr<<"...End Noise_calibration"<<endl;
}

/*********************************************************************/

void MEYER_WT1D::extract_stat(fltarray *TabBand, char* Outname)
{
	bool LocVerbose = false & Verbose;
	if(Verbose) cerr<<"Extract_stat..."<<endl;
	
	TabStat = new double[NbrScale];
	MaxCoef = new double[NbrScale-1];//we ommit the coarse scale
	
// Output stat file
	char Statname[250];
	strcpy(Statname, Outname);
	strcat(Statname, "_stat.dat");
	fstream cstat;
	cstat.open (Statname, fstream::out);
	char MaxCoefname[250];
	strcpy(MaxCoefname, Outname);
	strcat(MaxCoefname, "_maxcoef.dat");
	fstream ccoef;
	ccoef.open (MaxCoefname, fstream::out);
	
	
// fine scales : 0..N-2
	for(int s=0;s<NbrScale-1;s++)
	{
		double mean;
		int maxx;
		
		mean=(TabBand[s]).mean();
		TabStat[s]=(TabBand[s]).sigma()/Tab_1D_Meyer[s];
		MaxCoef[s]=abs((TabBand[s]).maxfabs(maxx))/Tab_1D_Meyer[s];
		
		if(LocVerbose)
		{
			cerr<< "\tStat echelle "<<s<<" : mean,sigma = "<<mean<<" "<<TabStat[s]<<endl;
			cerr << "\t MaxCoef Value = "<< MaxCoef[s] << " at "<< maxx << endl;
		}
		cstat << mean <<"\t"<<TabStat[s]<<"\t"<<endl;
		ccoef << s  << "\t" << MaxCoef[s] << "\t" << maxx << endl;
	}
	
// coarse scale : N-1
// No calibration for the coarsest scale : no "/TabSigma(N-1,0)" defined, only stats
	{
		int s=NbrScale-1;
		double mean;
		
		mean=(TabBand[s]).mean();
		TabStat[s]=(TabBand[s]).sigma();
		
		if(LocVerbose)
			cerr<< "\tStat echelle "<<s<<" : mean,sigma = "<<mean<<" "<<TabStat[s]<<endl;
		cstat << mean <<"\t"<<TabStat[s]<<"\t"<<endl;
	}

	cstat.close();
	ccoef.close();
	
	if(Verbose) cerr<<"...End Extract_stat"<<endl;
}

/*********************************************************************/

void MEYER_WT1D::normalize(fltarray *TabBand, fltarray *TabBandNorm)
{
	if(Verbose) cerr<<"RCurvelet3D::normalize..."<<endl;
	
	// Fine scales
	for(int s=0;s<NbrScale-1;s++)
	{
		for (int i=0; i < TabNx(s); i++)
			(TabBandNorm[s])(i)=(TabBand[s])(i)/Tab_1D_Meyer[s];
	}
	
	// Coarse Scale : no normalisation
	{
		int s=NbrScale-1;
		for (int i=0; i < TabNx(s); i++)
			(TabBandNorm[s])(i)=(TabBand[s])(i);
	}

	if(Verbose) cerr<<"...End RCurvelet3D::normalize"<<endl;
}

/*********************************************************************/

void MEYER_WT1D::normalize_self(fltarray *TabBand, bool inverse)
{
	if(Verbose) cerr<<"RCurvelet3D::normalize_self(.,"<<inverse<<")"<<endl;
	
	// Fine scales
	for(int s=0;s<NbrScale-1;s++)
	{
		for (int i=0; i < TabNx(s); i++)
			if(inverse)
				(TabBand[s])(i) *= Tab_1D_Meyer[s];
			else
				(TabBand[s])(i) /= Tab_1D_Meyer[s];
	}
	
	// Coarse Scale : no normalisation
	
	if(Verbose) cerr<<"...End RCurvelet3D::normalize_self"<<endl;
}

/*********************************************************************/

void MEYER_WT1D::threshold(fltarray *TabBand, float SigmaNoise, float NSigma, filter_type FilterType, bool force3)
{
	if(Verbose) cerr<<"MEYER_WT1D::threshold(.,"<<","<<SigmaNoise<<","<<NSigma<<","<<FilterType<<","<<force3<<")"<<endl;

	float lvl;
	// Fine scales
	for(int s=0;s<NbrScale-1;s++)
	{
		float cnt=0;
		int nx = TabBand[s].nx();
		
		for(int i = 0; i < nx; i++)
		{
			float Nsig=NSigma;
			if(s==0) Nsig=(force3 ? NSigma : (NSigma+1));
			lvl = SigmaNoise * Nsig * Tab_1D_Meyer[s];

			if( abs((TabBand[s])(i)) < lvl )
			{
				cnt++;
				(TabBand[s])(i)=0; // hard
			}
			else if(FilterType==FT_SOFT) (TabBand[s])(i) -= (2*int( (TabBand[s])(i)>0 )-1)*lvl;
		}
		if(Verbose) cerr<<" Scale "<<s<<", n#proportion non seuillee ("<<lvl<<")="<<nx-cnt<<"#"<<(nx-cnt)/nx<<endl;
	}

	// Coarse scale
	// Do not denoise the coarsest scale

	if(Verbose) cerr<<"...End MEYER_WT1D::threshold"<<endl;
}

/*********************************************************************/

void MEYER_WT1D::wiener(fltarray *TabBand, float noise_lvl, int LocalBS)
{
	if(Verbose) cerr<<"MEYER_WT1D::wiener("<<noise_lvl<<","<<LocalBS<<")..."<<endl;
//	bool LocVerbose = false & Verbose;
	
	// Fine scales
	for(int s=0;s<NbrScale-1;s++)
	{
		int nx = nxs(s);
		int Nx = nx/LocalBS + 1 ;
		fltarray coef_wiener(Nx);

		coef_wiener.init(-2);

		// Wiener blocks
		for(int kx=0 ; kx < Nx ; kx++)
		{
			double sigma2 = 0.0;
			float cnt=0;

		// Sigma calculation
			// Pixels in a wiener block
			for(int bx = 0; bx < LocalBS; bx++)
			{
				cnt+=1;
				// spatial position = position in the block + block_position
				int x = (bx + kx*LocalBS) % nx;
				sigma2+=pow((TabBand[s])(x) / Tab_1D_Meyer[s],2);
			}
			float sig = sqrt(max( 0.0, sigma2/cnt - pow(noise_lvl,2) ));
			float norm = sig / (sig+noise_lvl);
			coef_wiener(kx) = norm;
		}


	// Apply the wiener coefficient
		// Wiener blocks (int angle2 and space) for mainly vertical Ridgelets
		for(int kx=0 ; kx < Nx ; kx++)
			for(int bx = 0; bx < LocalBS; bx++)
			{
				int x = (bx + kx*LocalBS);
				if( x<nx )
					(TabBand[s])(x) *= coef_wiener(kx);
			}
	}// end scale

	// Coarse scale
	// Do not denoise the coarsest scale

	if(Verbose) cerr<<"...End MEYER_WT1D::wiener"<<endl;
}

/*********************************************************************/

void MEYER_WT1D::values_at(fltarray *TabBand, char * filename, char* Outname)
{
	int s,x;
	float m;
	
	// Output stream
	char MaxCoefname[250];
	strcpy(MaxCoefname, Outname);
	strcat(MaxCoefname, "_values_at.dat");
	fstream cval;
	cval.open (MaxCoefname, fstream::out);
	
	FILE* fic;
	fic=fopen(filename, "r");

	if(fic)
	{
		for(int s3=0;s3<NbrScale-1;s3++)
		{
			if(!fscanf(fic,"%d %f %d",&s,&m,&x)) cerr<<"Error while reading coordinates in "<<filename<<endl;
			cval<<s<<"\t"<<abs((TabBand[s])(x))<<"\t"<<x<<endl;
		}
		fclose(fic);
	}
	else cerr<<"Warnig: File "<<filename<<" not found, for use in BCurvelet3D::values_at"<<endl;
	cval.close();
}

/*********************************************************************/

void MEYER_WT1D::transform(fltarray &Data)
{
	if(Verbose) cerr<<"MEYER_WT1D::transform(.)..."<<endl;
	bool LocVerbose = false & Verbose;
	get_extFourier(Data, TF_ExtData);

	if(LocVerbose) 
	{
		char filename[64];
		fltarray TabWavelet;
		cfarray TabWaveletF;
		TabWaveletF=TF_ExtData;

		cerr<<"size TF = "<<TF_ExtData.nx()<<endl;
		TabWavelet.alloc(nxs(0));

		// save cube in fourier space
		for (int i=0; i < nxs(0); i++) TabWavelet(i) = (TabWaveletF)(i).real();
		//sprintf(filename,"%s_F.fits","out");
		//writefltarr(filename, TabWavelet);
	}

	transform_cf(TF_ExtData, Tabcf_WT_Band);
	
	if(LocVerbose) 
	{
		char filename[64];
		for  (int s=0; s < NbrScale; s++)
		{
			fltarray TabWavelet;
			cfarray TabWaveletF;
			TabWaveletF=Tabcf_WT_Band[s];

//			float Norm =  sqrt((float)(TabNx(s)));
			TabWavelet.alloc(TabNx(s));

			// save module in fourier space
			for (int i=0; i < TabNx(s); i++) TabWavelet(i) = (TabWaveletF)(i).real();
			//sprintf(filename,"%s_FW%d.fits","out",s);
			//writefltarr(filename, TabWavelet);
		}
	}
	if(Verbose) cerr<<"...End MEYER_WT1D::transform"<<endl;
}

/*********************************************************************/

void MEYER_WT1D::transform(fltarray &Data, fltarray * & Tab_Cube, Bool Alloc)
{
	if(Verbose) cerr<<"MEYER_WT1D::transform(.,.,"<<Alloc<<")..."<<endl;
	transform(Data);
	ifft_tabcube(Tabcf_WT_Band, Tab_Cube, Alloc);
	if(Verbose) cerr<<"End MEYER_WT1D::transform"<<endl;
}

/*********************************************************************/

void MEYER_WT1D::recons(fltarray * & Tab_Cube, fltarray &Data, Bool Alloc)
{
   if (Alloc == True) Data.resize(Nx_Cube);
   fft_tabcube(Tab_Cube, Tabcf_WT_Band);
   recons(Data);
}

/*********************************************************************/

void MEYER_WT1D::recons(fltarray &Data)
{
	recons_cf(Tabcf_WT_Band, TF_ExtData);

/*	if(Verbose) 
	{
		char filename[64];
		fltarray TabWavelet;
		cfarray TabWaveletF;
		TabWaveletF=TF_ExtData;

		cerr<<"size rTF = "<<TF_ExtData.nx()<<endl;
		TabWavelet.alloc(nxs(0), nys(0), nzs(0));

		// save module in fourier space
		for (int i=0; i < nxs(0); i++)
		for (int j=0; j < nys(0); j++)
		for (int k=0; k < nzs(0); k++) TabWavelet(i) = (TabWaveletF)(i).real();
		sprintf(filename,"%s_rF.fits","out");
		writefltarr(filename, TabWavelet);
	}*/

	get_IFFT_Cube(TF_ExtData, Data);
}

/*********************************************************************/

void MEYER_WT1D::test(fltarray *TabBand)
{
/*
	for(int s3=0;s3<NbrScale-1;s3++)
	{
		cerr<<"s nx"<<s3<<" "<<TabBand[s3].nx()<<endl;
		for(int i=0;i<TabBand[s3].nx();i++)
			TabBand[s3](i)=0;
	}
*/
}

/*********************************************************************/

void MEYER_WT1D::write (char *Name, fltarray * & data, bool Normalize)
{
	if(Verbose) cerr<<"MEYER_WT1D::write_multi..."<<endl;
	
	char filename[256];
	fitsfile *fptr;    
	int status;
	int simple;
	int bitpix;
	long naxis=0;
	long naxes[3];
	long group = 1; 

	// .mr extention
	mr_io_name (Name, filename);

	FILE *FEXIST = fopen(filename, "rb");
	if (FEXIST)
	{
		fclose(FEXIST);
		remove(filename);               // Delete old file if it already exists 
	}
	status = 0;         // initialize status before calling fitsio routines 

// open the file
	if ( ffinit(&fptr, filename, &status) )	// create the new FITS file 
		PrintError( status );					// call PrintError if error occurs 
	
// write  the header 
	simple   = True;
	bitpix   =  -32;   // 32-bit real pixel values      
	long pcount   =   0;  // no group parameters 
	long gcount   =   1;  // only a single image/group 

// write first header part (parameters)
	naxis=0;
	if (ffphps(fptr, bitpix, naxis, naxes, &status))
		PrintError( status );  
	// write optional keyword to the header 
		if ( ffpkyj(fptr, (char*)"Type_Tra", (long) 0, (char*)"Meyer 3D Wavele Transform", &status))
			PrintError( status );  
		if ( ffpkyj(fptr, (char*)"Normaliz", (long) Normalize, (char*)"1 if the transform is normalized, else 0", &status))
			PrintError( status );  
		if ( ffpkyj(fptr, (char*)"NbrScale", (long) NbrScale, (char*)"Number of bands 3D", &status))
			PrintError( status );  
		if ( ffpkyj(fptr, (char*)"Nx_Cube", (long) Nx_Cube, (char*)"x size of the original sinal", &status))
			PrintError( status );  
	
// write other headers and associated data	
	for (int s=0; s < nbr_scale(); s++)
	{
		naxis=1;
		naxes[0] = nxs(s);

		if(ffcrhd(fptr,&status))
			PrintError( status );
		if ( ffphpr(fptr,simple,bitpix,naxis,naxes,pcount,gcount,False,&status) )
			PrintError( status );
		
	// save the data
		if ( ffppre(fptr, group, 1, nxs(s), (data[s]).buffer(), &status) )
			PrintError( status ); 
	}

// close the FITS file 
	if ( ffclos(fptr, &status) )  PrintError( status );  
	
	if(Verbose) cerr<<"...end MEYER_WT1D::write_multi"<<endl;
}

/*********************************************************************/

void MEYER_WT1D::read(char *Name, fltarray * & data, bool *Normalize)
{
	if(Verbose) cerr<<"MEYER_WT1D::read_multi..."<<endl;
	char filename[256];
	fitsfile *fptr;           // pointer to the FITS file 
	int status=0, hdutype ;
	char comment[FLEN_COMMENT];
	long mon_long;
	int anynul = 0;
	long nulval = 0;
	void PrintError( int status);

	mr_io_name (Name, filename);

// open the file 
	status = 0;         // initialize status before calling fitsio routines 
	if ( ffopen(&fptr, filename, (int) READONLY, &status) ) 
		PrintError( status );

// get number of Headers
	int nhead;
	fits_get_num_hdus(fptr, &nhead, &status);
	
// read primary header
	if ( ffmahd(fptr, 1, &hdutype, &status) ) PrintError( status );
	
	// Read params
	int _NbrScale,_Nx;
	if (ffgkyj(fptr,(char*)"Type_Tra", &mon_long, comment, &status)) PrintError( status );
	if (ffgkyj(fptr,(char*)"Normaliz", &mon_long, comment, &status)) PrintError( status );
	*Normalize = (bool)mon_long;
	if (ffgkyj(fptr,(char*)"NbrScale", &mon_long, comment, &status)) PrintError( status );
	_NbrScale = (int)mon_long;
	if (ffgkyj(fptr,(char*)"Nx_Cube", &mon_long, comment, &status)) PrintError( status );
	_Nx = (int)mon_long;
	
	init(_NbrScale, _Nx);
	if(nhead!=_NbrScale+1) { cerr<<"Wrong header number in fits file : Hdr="<<nhead<<", NbrScale="<<NbrScale<<endl; exit(0); }
	
//	init the structure and the output vector
	data = new fltarray[NbrScale];
	for(int s=0; s < NbrScale; s++)
		data[s].alloc(TabNx(s));
	
// read data
	for(int s=0;s<_NbrScale;s++)
	{
		if (fits_movabs_hdu(fptr, s+2, NULL, &status)) PrintError( status );
		if (ffgpve(fptr, 1, 1, nxs(s), nulval, (data[s]).buffer(), &anynul, &status)) PrintError( status );
	}
	
// close the FITS file 
	if ( ffclos(fptr, &status) ) PrintError( status );
	
	if(Verbose) cerr<<"...end MEYER_WT1D::read_multi"<<endl;
}



/*********************************************************************/

