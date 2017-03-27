/******************************************************************************
**                   Copyright (C) 1994 by CEA
*******************************************************************************
**
**    UNIT
**
**    Version: 3.2
**
**    Author: Jean-Luc Starck
**
**    Date:  96/06/13
**    
**    File:  mr_info.cc
**
*******************************************************************************
**
**    DESCRIPTION  image multiresolution information  
**    ----------- 
**                 
**    PARAMETRES    
**    ----------    
** 
**    USAGE: mr_info option in_image 
**        where options = 
**
**          [-t type_of_multiresolution_transform]
**                  1: linear wavelet transform: a trous algorithm 
**                  2: bspline wavelet transform: a trous algorithm 
**                  3: wavelet transform in Fourier space 
**                  4: morphological median transform 
**                  5: morphological minmax transform 
**                  6: pyramidal linear wavelet transform 
**                  7: pyramidal bspline wavelet transform 
**                  8: pyramidal wavelet transform in Fourier space: 
**                     wavelet =  between two resolutions 
**                  9: pyramidal wavelet transform in Fourier space: 
**                     wavelet = difference  between the square 
**                                                of two resolutions
**                 10: pyramidal median transform 
**                 11: morphological pyramidal minmax transform 
**                 12: pyramidal laplacian 
**                 13: decomposition on scaling function 
**                 18: Feauveau's wavelet transform without undersampling 
**
**                 default is 2
**                 17 is not yet implemented
**
**           [-p]
**                Poisson noise
**                default is gaussian noise
**
**           [-g sigma]
**                Gaussian noise
**                  sigma = noise standard deviation 
**                by default, the noise is gaussian, and the standard 
**                devaition is automatically estimated. 
**
**           [-c gain,sigma,mean]
**                case of a CCD: noise = Poisson noise + read-out noise
**                  gain = CCD gain 
**                  sigma = standard deviation of the read-out noise
**                  mean = mean of the read-out noise
**                if this option is set, 
**                           Noise = Poisson + Gaussian read-out Noise
**                it is generally the case with the CCD.
**                Attention, these parameters must be separated by a comma 
**                without space. example: -c 0.133,1.733,0.
**                If mean, or sigma and mean are omitted, default values are 0.
**                gain can not be omitted. 
**
**           [-n number_of_scales]
**                number of scales used in the multiresolution transform
**                default is 4
**
**           [-s NSigma]
**                Thresolding at NSigma * SigmaNoise at each scale
**                default is 3
**
**           [-k]
**                kill isolated pixels in the multiresolution support
**                if the PSF is large, then isolated pixels are certainly 
**                residual noise, or cosmic rays, or artefacts. If this
**                option is set, we suppress these pixels in the support
**                default is no.
**
**           [-l]
**                dilate the multiresolution support
**                if this option is set, each scale is dilated by using 
**                the mathematical morphology operator. Usefull if artefacts
**                remain arround objects. 
**                default is no
**
**           [-w support_file_name]
**                if this option is set, an image is created from the 
**                multiresolution support. Default is no.
**
**
******************************************************************************/
 
#include "IM_Obj.h"
#include "IM_IO.h"
#include "MR_Obj.h"
#include "MR_Support.h"
#include "MR_NoiseModel.h"
#include "IM_Math.h"

char Name_Imag_In[256];  /* input file image */
char Name_Imag_Out[256]; //  output file name  
int Nbr_Plan=DEFAULT_NBR_SCALE;  /* number of scales */
float N_Sigma=DEFAULT_N_SIGMA;   /* number of sigma (for the noise) */
float Noise_Ima=0.;              /* noise standard deviation */
float Fwhm=0.5;                   /* Full width at half maximum */
type_noise Stat_Noise = DEFAULT_STAT_NOISE;   /* type of noise */
type_transform Transform = TO_PAVE_BSPLINE; /* type of transform */
Bool Dil=False;                 /* dilate the support */
Bool SupIsol=False;             /* suppress isolated pixel in the support */
char Name_Write_Sup[256];          /* output support file name */
Bool WriteSup = False;          /* write the support on the disk */
Bool ObjAna = False;
extern float PasCodeur;  /* CCD gain */
extern float SigmaGauss; /* CCD read-out noise standard deviation */
extern float MeanGauss;  /* CCD read-out noise mean */

extern int  OptInd;
extern char *OptArg;

extern int  GetOpt(int argc, char *const*argv, char *opts);
Bool UseNSigma=False;
Bool Verbose = False;

sb_type_norm Norm = NORM_L1;
type_sb_filter SB_Filter = F_MALLAT_7_9;
int NbrUndec = -1;                     /*number of undecimated scale */
Bool WriteStat = False;
int BorderSize=0;
Bool AllStat = False;
Bool NormAllBand = False;

/*********************************************************************/

static void usage(char *argv[])
{
 
    fprintf(OUTMAN, "Usage: %s options in_image [out_file]\n\n", argv[0]);
    fprintf(OUTMAN, "   where options =  \n");
    
    manline();
    transform_usage(Transform);
    manline();
    sb_usage(SB_Filter);
    manline();
    nbr_nbr_undec_usage();
    manline();
    nbr_scale_usage(Nbr_Plan);
    manline();
    poisson_noise_usage();
    manline();
    gauss_usage();
    manline();
    ccd_usage();
    manline();
    nsigma_usage(N_Sigma);
    manline();
    kill_isol_pix_usage();
    manline();
    dilate_sup_usage();
    manline();
    analyse_struct_usage();
    fprintf(OUTMAN, "         [-B BorderSize]\n");
    fprintf(OUTMAN, "             Do not take into account coeffs on the borders. \n");    
    manline();
    vm_usage();
    manline();
    verbose_usage();    
    manline();
    manline();
    exit(-1);
}
 
/*********************************************************************/

/* GET COMMAND LINE ARGUMENTS */
static void infinit(int argc, char *argv[])
{
    Bool OptL = False, Optf = False;
    int c;
#ifdef LARGE_BUFF
    int VMSSize=-1;
    Bool OptZ = False;
    char VMSName[1024] = "";
#endif    
    /* get options */
    while ((c = GetOpt(argc,argv,"NSB:u:pklt:T:Lg:c:n:s:avzZ:")) != -1) 
    {
	switch (c) 
        {
	   case 'N': NormAllBand = True; break;
	   case 'S': AllStat = True; break;
	   case 'B': if (sscanf(OptArg,"%d",&BorderSize) != 1) 
                {
		    fprintf(OUTMAN, "Error: bad number of scales: %s\n", OptArg);
	            exit(-1);
                    
		}
   		break;
	   case 'u':
 		if (sscanf(OptArg,"%d",&NbrUndec) != 1) 
                {
		    fprintf(OUTMAN, "Error: bad number of scales: %s\n", OptArg);
	            exit(-1);
                    
		}
  		break;
	   case 'v': Verbose = True; break;
	   case 't':
		/* -t <type> type of transform */
		if (sscanf(OptArg,"%d",&c ) != 1) 
                {
		    fprintf(OUTMAN, "Error: bad type of multiresolution transform: %s\n", OptArg);
	            exit(-1);
                    
		}
                if ((c > 0) && (c <= NBR_TRANSFORM)) 
                                        Transform = (type_transform) (c-1);
                else  
                {
		    fprintf(OUTMAN, "Error: bad type of transform: %s\n", OptArg);
	            exit(-1);
 		}
		break;
	   case 'T': 
		Optf = True;
		SB_Filter = get_filter_bank(OptArg);
		break;
	   case 'L': Norm = NORM_L2; OptL = True; break;
           case 'p':
                /* Poisson noise */
                Stat_Noise = NOISE_POISSON;
                Noise_Ima = 1.;
                ObjAna = True;
               break;
	    case 'g':
		/* -g <sigma_noise> */
		if (sscanf(OptArg,"%f",&Noise_Ima) != 1) 
                {
		    fprintf(OUTMAN, "Error: bad sigma noise: %s\n", OptArg);
		    exit(-1);
		}
                Stat_Noise = NOISE_GAUSSIAN;
                ObjAna = True;
		break;
             case 'c':
		/* -c <gain sigma mean> */
                printf("OptArg = %s\n", OptArg);
		if (sscanf(OptArg,"%f,%f,%f", &PasCodeur,
                                              &SigmaGauss, &MeanGauss) <= 0) 
                {
		    fprintf(OUTMAN, "Error: bad noise parameter: %s\n", OptArg);
		    exit(-1);
		}
                Stat_Noise = NOISE_POISSON;
                Noise_Ima = 1.;
                ObjAna = True;
		break;
	   case 'n':
		/* -n <Nbr_Plan> */
		if (sscanf(OptArg,"%d",&Nbr_Plan) != 1) 
                {
		    fprintf(OUTMAN, "Error: bad number of scales: %s\n", OptArg);
		    exit(-1);
		}
                if ((Nbr_Plan <= 1) || (Nbr_Plan > MAX_SCALE)) 
                {
		    fprintf(OUTMAN, "Error: bad number of scales: %s\n", OptArg);
		    fprintf(OUTMAN, "       1 < Nbr Scales <= %d\n", MAX_SCALE);
		    exit(-1);
		}
		break;
	   case 's':
		/* -s <nsigma> */
		if (sscanf(OptArg,"%f",&N_Sigma) != 1) 
                {
		    fprintf(OUTMAN, "Error: bad N_Sigma: %s\n", OptArg);
		    exit(-1);
		}
                if (N_Sigma <= 0.)  N_Sigma = DEFAULT_N_SIGMA;
                ObjAna = True;
                UseNSigma= True;
		break;
            case 'k':
               /* kill i */
               SupIsol = True;
               ObjAna = True;
               break;
            case 'l':
               /* dilate the support */
                Dil = True;
                ObjAna = True;
               break;
 	   case 'a':
                ObjAna = True;
                break;
#ifdef LARGE_BUFF
	    case 'z':
	        if (OptZ == True)
		{
                   fprintf(OUTMAN, "Error: Z option already set...\n");
                   exit(-1);
                }
	        OptZ = True;
	        break;
            case 'Z':
	        if (sscanf(OptArg,"%d:%s",&VMSSize, VMSName) < 1)
		{
		   fprintf(OUTMAN, "Error: syntaxe is Size:Directory ... \n");
		   exit(-1);
		}
	        if (OptZ == True)
		{
                   fprintf(OUTMAN, "Error: z option already set...\n");
                   exit(-1);
                }
		OptZ = True;
                break;
#endif
            case '?': usage(argv); break;
	    default: usage(argv); break;
 		}
	}

       /* get optional input file names from trailing 
          parameters and open files */
	if (OptInd < argc) strcpy(Name_Imag_In, argv[OptInd++]);
         else usage(argv);

	if (OptInd < argc)
	{
	   strcpy(Name_Imag_Out, argv[OptInd++]);
	   WriteStat = True;
	}
	 
	/* make sure there are not too many parameters */
	if (OptInd < argc)
        {
		fprintf(OUTMAN, "Too many parameters: %s ...\n", argv[OptInd]);
		exit(-1);
	}

	if ( (isotrop(Transform) == False) && (ObjAna == True))
	{
	   cerr << "Error: significant coefficient analysis cannot be done " << endl;
	   cerr << "       with this transform." << endl;
	   exit(0);
	}

	if ((Transform != TO_UNDECIMATED_MALLAT) && (Transform != TO_MALLAT) && ((OptL == True) || (Optf == True)))
	{
	   fprintf(OUTMAN, "Error: option -T and -L are only valid with Mallat transform ... \n");
           exit(0);
	}

#ifdef LARGE_BUFF
    if (OptZ == True) vms_init(VMSSize, VMSName, Verbose);
#endif

}

/*********************************************************************/

float get_new_stat(float *Band, int N, float KSig)
{
   int i;
   float Sigma = get_sigma_mad(Band, N);
   float T = KSig*Sigma;
   double Sum = 0.;
   for (i=0; i < N; i++) if (ABS(Band[i]) > T) Sum += ABS(Band[i]); 
   Sum /= (Sigma*sqrt( (float) N));
   return (float) Sum;
}

/*********************************************************************/
// function hs, val
// x = abs(val)
// step = 0.01
// n = fix(x / step + 0.5)
// if x lt step then h = 0 $
// else begin
//   a = findgen(n) / float(n) * x + step
//   s2 = sqrt(2.)
//   u=a
//   E = u*errorf( (x-u) /S2)*step
//   h = total(E)
// end
// return, h
// end

float hs(float val)
{
   float u,h=0;
   float x = ABS(val);
   float step = 0.01;
   int n = (int)(x / step + 0.5);
   if (x >= step)
   {
      float s2 = sqrt(2.);
      for (int i=0; i < n; i++)
      {
         u = (float) i / (float) n * x + step;
 	 h += u*erf( (x-u) /s2)*step;
      }
   }
   return h;
}


/*********************************************************************/

double neg_entropy(float *TabVal, int N, double Mean, double Sigma, double  Param=2);

double neg_entropy(float *TabVal, int N, double Mean, double Sigma, double  Par)
{
    double Val, Neg=0.;

    for (int i=0; i < N; i++) 
    {
       Val = (TabVal[i] - Mean)/ Sigma; 
       Neg += log(cosh(Par*Val));
    }
    Neg /= (Par*N);
    return Neg;
}

/*********************************************************************/

int main(int argc, char *argv[]) 
{
    int i,j,s,Nls,Ncs;
    Ifloat Dat;
    int Nl, Nc;
    fitsstruct Header;
    float *x_ef, *y_ef, *f_ef, Coef;
    int nmax,Cpt;
    Ifloat Imag, Imag_Detect;
    float Level = 0.5;
    int *Tab_Histo;
    int NbrVal;
 
     /* Get command line arguments, open input file(s) if necessary */
    lm_check(LIC_MR1);
    infinit(argc, argv);
 
    if (Verbose == True)
    {
       cout << endl << endl << "PARAMETERS: " << endl << endl;
       cout << "File Name in = " << Name_Imag_In << endl;
       cout << "Transform = " << StringTransform(Transform) << endl;
       if ((Transform == TO_MALLAT) || (Transform == TO_UNDECIMATED_MALLAT))
       {
          cout << StringSBFilter(SB_Filter) << endl;
          if (Norm == NORM_L2) cout << "L2 normalization" << endl;
       }
       cout << "Number of scales = " << Nbr_Plan << endl;
       if (Transform == TO_UNDECIMATED_MALLAT) cout << "Number of undecimated scales = " <<  NbrUndec << endl;
        if (BorderSize != 0) cout << "BorderSize = " << BorderSize << endl;
    }
    io_read_ima_float(Name_Imag_In, Dat, &Header);
    Nl = Dat.nl();
    Nc = Dat.nc();
    check_scale(Dat.nl(), Dat.nc(), Nbr_Plan);

    MultiResol MR_Data;
    FilterAnaSynt FAS;
    FilterAnaSynt *PtrFAS = NULL;
    if ((Transform == TO_MALLAT) || (Transform == TO_UNDECIMATED_MALLAT))
    {
        FAS.Verbose = Verbose;
        FAS.alloc(SB_Filter);
        PtrFAS = &FAS;
    }
    MR_Data.alloc (Nl, Nc, Nbr_Plan, Transform, PtrFAS, Norm, NbrUndec);

    int NbrBand = MR_Data.nbr_band();
    MRNoiseModel ModelData;

    if (ObjAna == True)
    {
       MR_Data.Border = I_ZERO;
       ModelData.alloc(Stat_Noise, Nl,Nc,Nbr_Plan, Transform, PtrFAS, Norm, NbrUndec);
       if (Noise_Ima > FLOAT_EPSILON) ModelData.SigmaNoise = Noise_Ima;
       if (SupIsol == True) ModelData.SupIsol = True;
       if (Dil == True) ModelData.DilateSupport = True;    
       if (UseNSigma  == True)
                   for (i=0; i < NbrBand; i++) ModelData.NSigma[i]=N_Sigma;
       ModelData.CCD_Gain = PasCodeur;
       ModelData.CCD_ReadOutSigma = SigmaGauss;
       ModelData.CCD_ReadOutMean = MeanGauss;       ModelData.model(Dat, MR_Data);
       
       // if TransImag == true MR_Data contains the wavelet coefficient
       // of the transform image, and not from the image
       // then we compute the WT of the image
       if (ModelData.TransImag == True) MR_Data.transform (Dat, I_ZERO);
    }
    else MR_Data.transform (Dat);
    
    int NbrStatPerBand = 5;
    if (AllStat == True) NbrStatPerBand = 9;
    fltarray TabStat(NbrBand-1, NbrStatPerBand);
    for (s = 0; s < NbrBand-1; s++)
    {
       double Mean, Sigma, Skew, Curt;
       float  Min, Max;	
       Ifloat Band; 
       Band = MR_Data.band(s);
       Sigma = Band.sigma();
       if (NormAllBand == True)
       {
          for (i=0; i < Band.nl(); i++)
	  for (j=0; j < Band.nc(); j++)
	    Band(i,j) /= Sigma;
       }
       
       // moment4(MR_Data.band(s).buffer(), N, Mean, Sigma, Skew, Curt, Min, Max);
       im_moment4(MR_Data.band(s), Mean, Sigma, Skew, Curt, Min, Max, BorderSize);
       
       if (Verbose == True)
       {
          cout << "Band " << s+1 << ": Min = " << Min << ", Max = " << Max << ", Mean = " << Mean  <<  ", Sigma = " << Sigma << endl; 
          cout << "   Skewness = " << Skew  << ", Kurtosis = " << Curt << endl;
       }
       TabStat(s, 0) = (float) Sigma;
       TabStat(s, 1) = (float) Skew;
       TabStat(s, 2) = (float) Curt;
       TabStat(s, 3) = (float) Min;
       TabStat(s, 4) = (float) Max; 
       
       if (AllStat == True) 
       {
          float HC1, HC2, TC1, TC2, Q1, Q2;
          fltarray Buff;
          float *PtrBuff = NULL;
       	  int N = MR_Data.size_band_nl(s) * MR_Data.size_band_nc(s);
          if (BorderSize > 0)
          {
             int p=0;
             Buff.resize(N);
	     PtrBuff = Buff.buffer();
	     for (i=BorderSize; i < MR_Data.size_band_nl(s)-BorderSize; i++)
             for (j=BorderSize; j < MR_Data.size_band_nc(s)-BorderSize; j++)
                 Buff(p++) = MR_Data(s,i,j);
 	     N = p;
          }
          else PtrBuff = MR_Data.band(s).buffer();
	  
	  // Multiscale EntropyCMemWave
	  //float SigMad = get_sigma_mad(PtrBuff, N);
	  //float ME=0.;
	  //float C1=2.*SigMad*SigMad;
	  //float C2=sqrt(2.)*SigMad;
	  //for (i=0;i < N; i++) ME += hs(Buff(i)/SigMad);
	  //             // C1*Buff(i)*Buff(i)*erf( ABS(Buff(i)) / C2);
	  //TabStat(s, 3) = ME;
	         	
          gausstest(PtrBuff, N, TC1, TC2);
          if (s != NbrBand-1) im_hc_test(MR_Data.band(s), HC1, HC2, BorderSize, True);
          else im_hc_test(MR_Data.band(s), HC1, HC2, BorderSize);
	  // Q1 = get_new_stat(PtrBuff, N, 2.);
	  // Q2 = get_new_stat(PtrBuff, N, 2.5);
      
          TabStat(s, 5) = (float) HC1;
          TabStat(s, 6) = (float) HC2;
	  TabStat(s, 7) = (float) TC1;
          TabStat(s, 8) = (float) TC2;
	  TabStat(s, 8) = (float) neg_entropy(PtrBuff, N, MR_Data.band(s).mean(), MR_Data.band(s).sigma(), 1.);
	  // TabStat(s, 9) = (float) Q1;
          // TabStat(s, 10) = (float) Q2;
          if (Verbose == True)
          {
              cout << "    HC1 = " << HC1 << ", HC2 = " << HC2 << endl;
              cout << "    T1 = " << TC1  << ", T2 = " << TC2 << endl;
	  }
       }
       
       if ((BorderSize > 0) &&
           ((MR_Data.Set_Transform == TRANSF_PAVE) || 
	    (MR_Data.Set_Transform == TRANSF_UNDECIMATED_MALLAT) ||
	    (MR_Data.Set_Transform == TRANSF_DIADIC_MALLAT))) BorderSize *= 2;

 
       if (ObjAna == True)
       {
          Nls = MR_Data.size_band_nl(s);
          Ncs = MR_Data.size_band_nc(s);
          Imag.resize(Nls, Ncs);
          Imag_Detect.resize(Nls, Ncs);
          Imag.init();

          // Thresholding
          Cpt = 0;
          for (i = 0; i < Nls; i++)
          for (j = 0; j < Ncs; j++)
          if (ModelData(s,i,j) == True)
          {
             Imag(i,j)= MR_Data(s,i,j);
             Cpt++;
          }
          Coef = (float) Cpt / (float)(Nls*Ncs) * 100.;
          printf("    Significant Coefficients = %5.2f %%\n\n",Coef);
                   
          im_detect(Imag, Imag_Detect, &x_ef, &y_ef,  &f_ef, &nmax);
          cout << "    Number of maxima = " << nmax << endl;

          // Binarization
          for (i = 0; i < Nls; i++)
          for (j = 0; j < Ncs; j++) if (ModelData(s,i,j) == True) Imag(i,j)=1.;
          im_segment (Imag, Imag_Detect, nmax, Level);
          cout << "    Number of structures = " << nmax << endl;  
          im_histo (Imag_Detect, &Tab_Histo, NbrVal);
          Max = 0;
          for (i = 1; i < NbrVal; i++)
                       if (Tab_Histo [i] > Max) Max = Tab_Histo [i];
          delete [] Tab_Histo;
          cout << "    Bigger structure (in pixel) = " << Max << endl; 
          cout << endl; 
       }
    }
    
    if (WriteStat == True) fits_write_fltarr(Name_Imag_Out, TabStat); 
    exit(0);
}

/*********************************************************************/
