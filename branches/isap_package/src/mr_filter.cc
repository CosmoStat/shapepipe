/******************************************************************************
**                   Copyright (C) 1995 by CEA
*******************************************************************************
**
**    UNIT
**
**    Version: 3.3
**
**    Author: Jean-Luc Starck
**
**    Date:  96/07/03
**    
**    File:  mr_filter.cc
**
*******************************************************************************
**
**    DESCRIPTION  filter an image
**    ----------- 
**                 
**    PARAMETRES    
**    ----------    
** 
**    USAGE: mr_filter option image output
**        where options = 
**          [-f type_of_filtering]
**                 1: Multiresolution thresholding 
**                 2: Iterative multiresolution thresholding 
**                 3: Hierarchical thresholding 
**                 4: Multiresolution Filtering
**                 5: Hierarchical Wiener filtering 
**                 6: Multiresolution Wiener filtering 
**                 7: Median filtering 
**                 8: Average filtering 
**                 9: Bspline filtering 
**                 10: Anisotropic diffusion filtering 
**
**                 default is 2
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
**                 14: Mallat's wavelet transform 
**                 15: G transform (morphological min-max algorithm 
**                 16: Feauveau's wavelet transform 
**                 17: Haar's wavelet transform 
**                 18: Feauveau's wavelet transform without undersampling 
**
**                 default is 2
**                 17 is not yet implemented
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
**           [-e Epsilon]
**                Convergence parameter
**                default is 0.001
**
**           [-i number_of_iterations] 
**                Maximum number of iterations 
**
**           [-w support_file_name]
**                if this option is set, an image is created from the 
**                multiresolution support. Default is no.
**
**           t,p,g,c,n, options are only used if
**           type_of_filtering in [1..5]
**
**           s option is only used with type_of_filtering in [1,3]
**
**           e,i,w are used with type_of_filtering = 2
**    
**  
**
**
**
******************************************************************************/
 
#include "IM_Obj.h"
#include "IM_IO.h"
#include "MR_Obj.h"
#include "MR_Filter.h"
#include "MR_NoiseModel.h"
#include "MR_Abaque.h"
#include "MR_Psupport.h"
#include "NR.h"
#include "MR_Threshold.h"
#include <iomanip>

char Name_Imag_In[1024]; /* input file image */
char Name_Imag_Out[1024]; /* output file name */
char Name_Write_Sup[1024];          /* output support file name */
int Nbr_Plan=DEFAULT_NBR_SCALE;  /* number of scales */
float N_Sigma=DEFAULT_N_SIGMA;   /* number of sigma (for the noise) */
float Noise_Ima=0.;              /* noise standard deviation */
type_noise Stat_Noise = DEFAULT_STAT_NOISE;   /* type of noise */
type_transform Transform = DEFAULT_MR_TRANS; /* type of transform */
int Max_Iter = DEFAULT_MAX_ITER_FILTER; /* Maximum number of iteration */
type_filter Filter = FILTER_THRESHOLD;  /* default filtering method */
float Epsilon = DEFAULT_EPSILON_FILTERING; /* convergence parameter */
Bool UseSupport = True; /* Use the multiresolution support */
Bool WriteSup = False;          /* write the support on the disk */
Bool SupIsol=False;             /* suppress isolated pixel in the support */
Bool SupDil=False;              /* dilate the support */
Bool UseNSigma =False;
Bool PositivIma = DEF_POSITIV_CONSTRAINT;
Bool MaxIma = DEF_MAX_CONSTRAINT;
Bool KillLastScale =False;
Bool WriteInfoOnProbMap = False;

char Name_RMSMap[1024]; 
Bool UseRMSMap=False;
 
double EpsilonPoisson = DEFAULT_EPSILON;
Bool KeepPositivSup = False;
int SizeBlock = DEFAULT_SIZE_BLOCK_SIG;
int NiterClip = 1;
int FirstScale = DEF_FIRST_DETECT_SCALE;

extern float PasCodeur;  /* CCD gain */
extern float SigmaGauss; /* CCD read-out noise standard deviation */
extern float MeanGauss;  /* CCD read-out noise mean */

extern int  OptInd;
extern char *OptArg;

extern int  GetOpt(int argc, char **argv, char *opts);

Bool EpsOpt=False;
Bool MaxIterOpt=False;
Bool PosOpt=False;
Bool KillLastOpt = False;
Bool WindowOpt=False;    

Bool Verbose=False;
Bool GetEgde=False;

sb_type_norm Norm = NORM_L1;
type_sb_filter SB_Filter = F_MALLAT_7_9;
type_border Bord = I_CONT;
int NbrUndec = -1;                     /*number of undecimated scale */

Bool BackgroundImage = False;
char Name_Background[1024];
Bool MissingData=False;
float RegulParam = 0.1;
Bool OptRegul = False;
Bool PoissonFisz = False;
type_threshold TypeThreshold = DEF_THRESHOLD;
type_undec_filter U_Filter = DEF_UNDER_FILTER;
Bool UseFlatField = False;
char Name_FlatField[1024]; 
Bool PositiveReconsFilter = False;
int MinEventNumber=0;

Bool OldPoisson = False;
Bool WriteThreshold = False;

Bool ProbMR =False;
char ProbMRFile[1024];
Bool MadCenter=False;
float TabNSigma[15];
int NbrTabNSigma=0;

extern "C" {
   double sap_erf( double );
   double sap_erfc( double );
};
double dierfc(double y);


extern float BadPixalVal; 
extern Bool BadPixel;

/****************************************************************************/

static void usage(char *argv[])
{
    fprintf(OUTMAN, "Usage: %s options in_image out_image\n\n", argv[0]);
    fprintf(OUTMAN, "   where options =  \n");
    
    filter_usage(Filter);
     
    get_threshold_usage(TypeThreshold);
         
    transform_usage(Transform);
 
    fprintf(OUTMAN, "         [-T type_of_filters]\n");
    for (int i = 1; i <= NBR_SB_FILTER; i++)
    fprintf(OUTMAN, "              %d: %s \n",i,
                                           StringSBFilter((type_sb_filter  )i));
    fprintf(OUTMAN, "             default is %s\n", StringSBFilter ((type_sb_filter)  SB_Filter));
 
    usb_usage(U_Filter);
     
    nbr_nbr_undec_usage(NbrUndec);
          
    gauss_usage();
     
    ccd_usage();
     
    noise_usage();
 
    nbr_scale_usage(Nbr_Plan);
    nbr_scalep_usage(DEF_N_SCALE);
     
    nsigma_usage(N_Sigma);
    fprintf(OUTMAN, "             Default is 2 for FDR detection method.\n"); 
     
    max_iter_usage(Max_Iter);
     
    converg_param_usage(Epsilon);
    convergp_param_usage(DEF_EPS_EVENT_FILTERING);
     
    support_file_usage();
 
    kill_isol_pix_usage();
  
    kill_last_scale_usage();
     
    detect_pos_usage();
  
    prec_eps_poisson_usage(DEFAULT_EPSILON);
 
    size_block_usage(SizeBlock);
     
    sigma_clip_block_usage(NiterClip);
  
    first_detect_scale_usage();
  
    rms_noise_usage();
     
    fprintf(OUTMAN, "         [-P]\n");
    fprintf(OUTMAN, "             Suppress the positivity constraint.\n");
     
    fprintf(OUTMAN, "         [-b]\n");
    fprintf(OUTMAN, "             Add the maximum level constraint.\n");
    fprintf(OUTMAN, "             Max value is 255. Default is no.\n");
     
    fprintf(OUTMAN, "         [-B BackgroundModelImage]\n");
    fprintf(OUTMAN, "             Background Model Image: the background image is  \n");
    fprintf(OUTMAN, "             subtracted during the filtering.\n");
    fprintf(OUTMAN, "             Default is no.\n");
 	
	fprintf(OUTMAN, "         [-M Flat_Image]\n");
    fprintf(OUTMAN, "             Flat Image: The solution is corrected from the flat (i.e. Sol = Input / Flat)  \n");
    fprintf(OUTMAN, "             Default is no.\n");
    fprintf(OUTMAN, "         [-A]\n");
    fprintf(OUTMAN, "             If set, Use the second generation filter (positive reconstruction filter) in the starlet transform (i.e. a trous algorithm)..\n");
    fprintf(OUTMAN, "             Only valid with -t2 option. \n");
    fprintf(OUTMAN, "             Default is no.\n");
 	fprintf(OUTMAN, "         [-H]\n");
    fprintf(OUTMAN, "             If set, all pxiels zero values in the input image are considered as missing data.\n");
    fprintf(OUTMAN, "             Default is no.\n");
    fprintf(OUTMAN, "         [-h]\n");
    fprintf(OUTMAN, "             write info used for computing the probability map.\n");
    fprintf(OUTMAN, "             Default is no.\n");
 
    fprintf(OUTMAN, "         [-G RegulParam]\n");
    fprintf(OUTMAN, "              Regularization parameter for the TV method.\n");
    fprintf(OUTMAN, "              default is %f\n", RegulParam);
      
    vm_usage();
    verbose_usage();
     
    manline();
    manline();
    exit(-1);
}
  
/*********************************************************************/

/* GET COMMAND LINE ARGUMENTS */
static void filtinit(int argc, char *argv[])
{
    int c;
    Bool TransfOpt=False;
    Bool FilterOpt=False;
    Bool NscaleOpt=False;
    Bool EpsOptPoisson=False;
    Bool Optf=False;
#ifdef LARGE_BUFF
    int VMSSize=-1;
    Bool OptZ = False;
    char VMSName[1024] = "";
#endif     
    
    /* get options */
    while ((c = GetOpt(argc,argv, (char *) "HQ:AD:M:U:C:B:u:pbm:f:t:T:g:c:n:s:i:e:w:kE:S:N:F:KPR:vzZ:G:odh")) != -1) 
    { 
	switch (c) 
        {
	  case 'H': MissingData = True;  break; // MadCenter = True;
	       break;
	  case 'Q': if (sscanf(OptArg,"%s", ProbMRFile) != 1) 
                {
		   fprintf(OUTMAN, "Error: bad file name: %s\n", OptArg);
		   exit(-1);
		}
                ProbMR = True;
	       break;
	  case 'D': if (sscanf(OptArg,"%d",&MinEventNumber) != 1) 
                    {
		       fprintf(OUTMAN, "Error: bad type of detection: %s\n", OptArg);
	               exit(-1);
 		    }
		    break;
	   case 'A': PositiveReconsFilter = True; break;
	   case 'X': break; // PoissonFisz = True; break;
           case 'C':
		/* -f <type> type of filtering */
		if (sscanf(OptArg,"%d",&c ) != 1) 
                {
		    fprintf(OUTMAN, "Error: bad type of detection: %s\n", OptArg);
	            exit(-1);
 		}
                if ((c > 0) && (c <= NBR_THRESHOLD)) TypeThreshold = (type_threshold) (c-1);
                else  
                {
		    fprintf(OUTMAN, "Error: bad type of detection: %s\n", OptArg);
	            exit(-1);
 		}
 		break;
	  case 'B':
	       if (sscanf(OptArg,"%s", Name_Background) != 1) 
                {
		   fprintf(OUTMAN, "Error: bad file name: %s\n", OptArg);
		   exit(-1);
		}
                BackgroundImage = True;
	       break;
	   case 'u':
 		if (sscanf(OptArg,"%d",&NbrUndec) != 1) 
                {
		    fprintf(OUTMAN, "Error: bad number of scales: %s\n", OptArg);
	            exit(-1);
                    
		}
  		break;
 	   case 'G':
		/* -s <nsigma> */
		if (sscanf(OptArg,"%f",&RegulParam) != 1) 
                {
		    fprintf(stderr, "bad Regularization Parameter: %s\n", OptArg);
		    exit(-1);
		}
                if (RegulParam  < 0.) RegulParam = 0.1;
		OptRegul = True;
 		break;	   
	   case 'f':
		/* -f <type> type of filtering */
		if (sscanf(OptArg,"%d",&c ) != 1) 
                {
		    fprintf(OUTMAN, "Error: bad type of filtering: %s\n", OptArg);
	            exit(-1);
                    
		}
                if ((c > 0) && (c <= NBR_FILTERING)) Filter = (type_filter) (c-1);
                else  
                {
		    fprintf(OUTMAN, "Error: bad type of filtering: %s\n", OptArg);
	            exit(-1);
 		}
                FilterOpt = True;
		break;
	   case 't':
		/* -d <type> type of transform */
		if (sscanf(OptArg,"%d",&c ) != 1) 
                {
		    fprintf(OUTMAN, "ErrorL bad type of multiresolution transform: %s\n", OptArg);
	            exit(-1);
                    
		}
                if ((c > 0) && (c <= NBR_TRANSFORM)) 
                                        Transform = (type_transform) (c-1);
                else  
                {
		    fprintf(OUTMAN, "Error: bad type of transform: %s\n", OptArg);
	            exit(-1);
 		}
                TransfOpt = True;
 		break;
	   case 'U': 
		/* -d <type> type of transform */
		if (sscanf(OptArg,"%d",&c ) != 1) 
                {
		    fprintf(OUTMAN, "ErrorL bad type of filters: %s\n", OptArg);
	            exit(-1);
                    
		}
                if ((c > 0) && (c <= NBR_UNDEC_FILTER)) 
                                        U_Filter = (type_undec_filter) (c-1);
                else  
                {
		    fprintf(OUTMAN, "Error: bad type of filters: %s\n", OptArg);
	            exit(-1);
 		}
  		break;
	   case 'T': 
		Optf = True;
		SB_Filter = get_filter_bank(OptArg);
		break;
            case 'm':
		if (sscanf(OptArg,"%d",&c ) != 1) 
                {
		    fprintf(OUTMAN, "Error: bad type of noise: %s\n", OptArg);
	            exit(-1);
                    
		}
                if ((c > 0) && (c <= NBR_NOISE+1)) 
                                        Stat_Noise = (type_noise) (c-1);
                else  
                {
		    fprintf(OUTMAN, "Error: bad type of noise: %s\n", OptArg);
	            exit(-1);
 		}
		break;
	    case 'g':
		/* -g <sigma_noise> */
		if (sscanf(OptArg,"%f",&Noise_Ima) != 1) 
                {
		    fprintf(OUTMAN, "Error: bad sigma noise: %s\n", OptArg);
		    exit(-1);
		}
                Stat_Noise = NOISE_GAUSSIAN;
		break;
             case 'c':
		/* -c <gain sigma mean> */
		if (sscanf(OptArg,"%f,%f,%f", &PasCodeur,
                                              &SigmaGauss, &MeanGauss) <= 0) 
                {
		    fprintf(OUTMAN, "Error: bad noise parameter: %s\n", OptArg);
		    exit(-1);
		}
                Stat_Noise = NOISE_POISSON;
                Noise_Ima = 1.;
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
		NscaleOpt = True;
		break;
	   case 's':
		/* -s <nsigma> */
		if (sscanf(OptArg,"%f,%f,%f,%f,%f,%f,%f,%f,%f,%f", &TabNSigma[0], &TabNSigma[1], &TabNSigma[2],&TabNSigma[3] ,&TabNSigma[4],&TabNSigma[5],&TabNSigma[6],&TabNSigma[7],&TabNSigma[8],&TabNSigma[9]) == 10) NbrTabNSigma=10;
		else  if (sscanf(OptArg,"%f,%f,%f,%f,%f,%f,%f,%f,%f", &TabNSigma[0], &TabNSigma[1], &TabNSigma[2],&TabNSigma[3] ,&TabNSigma[4],&TabNSigma[5],&TabNSigma[6],&TabNSigma[7],&TabNSigma[8]) == 9) NbrTabNSigma=9;
        	else if ( sscanf(OptArg,"%f,%f,%f,%f,%f,%f,%f,%f", &TabNSigma[0], &TabNSigma[1], &TabNSigma[2],&TabNSigma[3] ,&TabNSigma[4],&TabNSigma[5],&TabNSigma[6],&TabNSigma[7]) == 8) NbrTabNSigma=8;
        	else if (sscanf(OptArg,"%f,%f,%f,%f,%f,%f,%f", &TabNSigma[0], &TabNSigma[1], &TabNSigma[2],&TabNSigma[3] ,&TabNSigma[4],&TabNSigma[5],&TabNSigma[6]) == 7) NbrTabNSigma=7;
        	else if (sscanf(OptArg,"%f,%f,%f,%f,%f,%f", &TabNSigma[0], &TabNSigma[1], &TabNSigma[2],&TabNSigma[3] ,&TabNSigma[4],&TabNSigma[5]) == 6) NbrTabNSigma=6;
        	else if (sscanf(OptArg,"%f,%f,%f,%f,%f", &TabNSigma[0], &TabNSigma[1], &TabNSigma[2],&TabNSigma[3] ,&TabNSigma[4]) == 5) NbrTabNSigma=5;
    	else if (sscanf(OptArg,"%f,%f,%f,%f", &TabNSigma[0], &TabNSigma[1], &TabNSigma[2],&TabNSigma[3] ) == 4) NbrTabNSigma=4;
       else if (sscanf(OptArg,"%f,%f,%f", &TabNSigma[0], &TabNSigma[1], &TabNSigma[2] ) == 3) NbrTabNSigma=3;
       else if (sscanf(OptArg,"%f,%f", &TabNSigma[0], &TabNSigma[1] ) == 2) NbrTabNSigma=2;
       else if (sscanf(OptArg,"%f",&N_Sigma) != 1) 
        	{
   		    fprintf(OUTMAN, "Error: bad N_Sigma: %s\n", OptArg);
		    exit(-1);
         	}
         UseNSigma =True;
          if (N_Sigma <= 0.)  N_Sigma = DEFAULT_N_SIGMA;
		break;
 	   case 'i':
		/* -i < Number of iterations> */
		if (sscanf(OptArg,"%d",&Max_Iter) != 1) 
                {
		   fprintf(OUTMAN, "Error: bad Max_Iter: %s\n", OptArg);
				exit(-1);
		}
                if (Max_Iter <= 0)   Max_Iter = DEFAULT_MAX_ITER_FILTER;
                MaxIterOpt=True;
		break;
	   case 'e': 
		/* -e < Convergence parameter> */
		if (sscanf(OptArg,"%f",&Epsilon) != 1) 
                {
		   fprintf(OUTMAN, "Error: bad convergence parameter: %s\n", OptArg);
		   exit(-1);
		}
                if ((Epsilon < 0) || (Epsilon > 1.)) 
                           Epsilon = DEFAULT_EPSILON_FILTERING;
		EpsOpt=True;
		break;
 	   case 'w':
		/* -w < support file name> */
		if (sscanf(OptArg,"%s", Name_Write_Sup) != 1) 
                {
		   fprintf(OUTMAN, "Error: bad file name: %s\n", OptArg);
		   exit(-1);
		}
                WriteSup = True;
 		break;
	   case 'M':	 
	   	/* -w < support file name> */
		if (sscanf(OptArg,"%s", Name_FlatField) != 1) 
                {
		   fprintf(OUTMAN, "Error: bad file name: %s\n", OptArg);
		   exit(-1);
		}
                UseFlatField = True;
 		break;
 	   case 'R':
		/* -w < support file name> */
		if (sscanf(OptArg,"%s", Name_RMSMap) != 1) 
                {
		   fprintf(OUTMAN, "Error: bad file name: %s\n", OptArg);
		   exit(-1);
		}
                UseRMSMap = True;
 		break;
            case 'k':
               /* kill i */
               SupIsol = True;
               break;
            case 'v': Verbose = True;break;
         case 'K':
               /* kill last scale */
               KillLastScale = True;
               KillLastOpt=True;
               break; 
            case 'P':
                PositivIma=False;
                PosOpt=True;
               break;            
	    case 'b':
                MaxIma=True;
               break;
	 case 'p':
		 KeepPositivSup=True;
		break; 
	 case 'h':
		 WriteInfoOnProbMap=True;
		break; 
	 case 'E':
	      { float f;
		if (sscanf(OptArg,"%f",&f) != 1) 
                {
		    fprintf(OUTMAN, "Error: bad precision number: %s\n", OptArg);
		    exit(-1);
		}
		EpsilonPoisson = (double) f;
		if ((EpsilonPoisson <= 0.) || (EpsilonPoisson > MAX_EPSILON)) 
                {
		    fprintf(OUTMAN, "Error: bad precision number: %s\n", OptArg);
		    fprintf(OUTMAN, "       %f <= Precision <= %f\n",MIN_EPSILON, MAX_EPSILON);
 		    exit(-1);
		}
		EpsOptPoisson=True;
	        }
		break; 
	 case 'S':
		if (sscanf(OptArg,"%d",&SizeBlock) != 1) 
                {
		    fprintf(OUTMAN, "Error: bad block size: %s\n", OptArg);
		    exit(-1);
		}
		if (SizeBlock  < 2)
		{
		   fprintf(OUTMAN, "Error: bad  SizeBlock parameter: %s\n", OptArg);
		   fprintf(OUTMAN, "              SizeBlock > 1\n");
		   exit(-1);
		}
		break; 
	 case 'N':
		if (sscanf(OptArg,"%d",&NiterClip) != 1) 
                {
		    fprintf(OUTMAN, "Error: bad it. number for the 3sigma clipping: %s\n", OptArg);
		    exit(-1);
		}
		if (NiterClip < 1)
		{
		   fprintf(OUTMAN, "Error: bad NiterClip parameter: %s\n", OptArg);
		   fprintf(OUTMAN, "             NiterClip > 0\n");
		   exit(-1);
		}
		break; 
	 case 'F':
		if (sscanf(OptArg,"%d", &FirstScale) != 1) 
                {
		   fprintf(OUTMAN, "Error: bad first detection scale number: %s\n", OptArg);
		   exit(-1);
		}
		FirstScale --;
		if (FirstScale < 0)
		{
		   fprintf(OUTMAN, "Error: bad FirstScale parameter: %s \n", OptArg);
		   fprintf(OUTMAN, "           FirstScale > 0\n");
		   exit(-1);
		}
		break;
            case 'o' : OldPoisson = True; break;
            case 'd' : WriteThreshold = True; break;
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

	if (OptInd < argc) strcpy(Name_Imag_Out, argv[OptInd++]);
         else usage(argv);

	/* make sure there are not too many parameters */
	if (OptInd < argc)
        {
		fprintf(OUTMAN, "Too many parameters: %s ...\n", argv[OptInd]);
		exit(-1);
	}

       if ((UseNSigma == False) && (TypeThreshold == T_FDR)) 
       {
          N_Sigma = 2.;
	  UseNSigma = True;
       }
              
 	// Test inconsistencies in the option call 
 	if (Stat_Noise == NOISE_EVENT_POISSON)
 	{
 	   if ((TransfOpt == True) && (Transform != TO_PAVE_BSPLINE))
 	   {
 	      cerr << "WARNING: with this noise model, only the BSPLINE A TROUS can be used ... " << endl;
 
 	      cerr << "        Type transform is set to: BSPLINE A TROUS ALGORITHM " << endl;
 	   }
 	   Transform = TO_PAVE_BSPLINE;
 	   if (NscaleOpt != True) Nbr_Plan = DEF_N_SCALE;
 	   if (EpsOpt == False) Epsilon = DEF_EPS_EVENT_FILTERING;
 	}
	
	if ((Stat_Noise == NOISE_CORREL) && (UseRMSMap != True))
	{
 	   cerr << endl << endl;
 	   cerr << "  Error: this noise model need a noise map (-R option) " << endl;
           exit(-1);
  	}

 	if (UseRMSMap == True)
 	{
 	   if ((Stat_Noise != NOISE_NON_UNI_ADD) && (Stat_Noise !=  NOISE_CORREL))
 	   {
 	      cerr << "Error: this noise model is not correct when RMS map option is set." << endl;
 	      cerr << "       Valid models are: " << endl;
 	      cerr << "        " << StringNoise(NOISE_NON_UNI_ADD) << endl;
	      cerr << "        " << StringNoise(NOISE_CORREL) << endl;
	      exit(-1);
   	   }
       }
 
	if ( ((SetTransform(Transform) == TRANSF_MALLAT) ||
	       (SetTransform(Transform) == TRANSF_FEAUVEAU))
 	      && ((Stat_Noise == NOISE_NON_UNI_ADD) ||
 	         (Stat_Noise  == NOISE_NON_UNI_MULT)))
  	{
 	   cerr << endl << endl;
 	   cerr << "  Error: with this transform, non stationary noise models are not valid : " << StringFilter(FILTER_THRESHOLD) << endl;
           exit(-1);
  	}
  	      
       // isolated pixel removal is not valid with
       // non isotropic transform.
       if ((isotrop(Transform) != True) && (SupIsol == True) && (Transform != TO_UNDECIMATED_MALLAT) && (Transform != TO_UNDECIMATED_NON_ORTHO))
       {
          fprintf(OUTMAN, "Error: option -k and -l are not valid with non isotropic transform. ...\n");
		exit(-1);
       }
       

       // Eps option and MaxIter option are valid only with iterative methods
       if ((EpsOpt == True) || (MaxIterOpt == True))
       {
          if (   (Filter != FILTER_ITER_THRESHOLD) 
              && (Filter != FILTER_ITER_ADJOINT_COEFF)
	      && (Filter != FILTER_WAVELET_CONSTRAINT)
	      && (Filter != FILTER_TV_CONSTRAINT))
          {
             fprintf(OUTMAN, "Error: option -e and -i are not valid with non iterative filtering methods. ...\n");
		exit(-1);
	  }
       }
       
       if ((EpsOptPoisson == False) && (UseNSigma == True)) {
          /*EpsilonPoisson = (1. - erf((double) N_Sigma / sqrt((double) 2.))) / 2;
          std::cout << EpsilonPoisson << std::endl;
	  EpsilonPoisson = sap_erfc( (double) N_Sigma / sqrt((double) 2.) ) / 2.;
          std::cout << EpsilonPoisson << std::endl;
	  EpsilonPoisson = dierfc( (double) N_Sigma / sqrt((double) 2.) ) / 2.;
          std::cout << EpsilonPoisson << std::endl;
	  EpsilonPoisson = erffc( (double) N_Sigma / sqrt((double) 2.) ) / 2.;
          std::cout << EpsilonPoisson << std::endl;*/
          if( OldPoisson ) 
	     EpsilonPoisson = (1. - erf((double) N_Sigma / sqrt((double) 2.))) / 2;
	  else
	     EpsilonPoisson = erffc( (double) N_Sigma / sqrt((double) 2.) ) / 2.;
          
	  if ((Stat_Noise == NOISE_EVENT_POISSON) && ( N_Sigma > 12 )) {
	     fprintf( OUTMAN, "Error: NSigma (option -s...) must be set to a lower value (<12) ..\n");
             exit(-1);
	  }
       } 


 	if ((Transform != TO_UNDECIMATED_MALLAT) && (Transform != TO_MALLAT) && (Optf == True))
	{
	   fprintf(OUTMAN, "Error: option -T is only valid with Mallat transform ... \n");
           exit(0);
	}
	     
#ifdef LARGE_BUFF
    if (OptZ == True) vms_init(VMSSize, VMSName, Verbose);
#endif  
}

 
/****************************************************************************/
  
double dierfc(double y)
{
    double s, t, u, w, x, z;

    z = y;
    if (y > 1) {
        z = 2 - y;
    }
    w = 0.916461398268964 - log(z);
    u = sqrt(w);
    s = (log(u) + 0.488826640273108) / w;
    t = 1 / (u + 0.231729200323405);
    x = u * (1 - s * (s * 0.124610454613712 + 0.5)) - 
        ((((-0.0728846765585675 * t + 0.269999308670029) * t + 
        0.150689047360223) * t + 0.116065025341614) * t + 
        0.499999303439796) * t;
    t = 3.97886080735226 / (x + 3.97886080735226);
    u = t - 0.5;
    s = (((((((((0.00112648096188977922 * u + 
        1.05739299623423047e-4) * u - 0.00351287146129100025) * u - 
        7.71708358954120939e-4) * u + 0.00685649426074558612) * u + 
        0.00339721910367775861) * u - 0.011274916933250487) * u - 
        0.0118598117047771104) * u + 0.0142961988697898018) * u + 
        0.0346494207789099922) * u + 0.00220995927012179067;
    s = ((((((((((((s * u - 0.0743424357241784861) * u - 
        0.105872177941595488) * u + 0.0147297938331485121) * u + 
        0.316847638520135944) * u + 0.713657635868730364) * u + 
        1.05375024970847138) * u + 1.21448730779995237) * u + 
        1.16374581931560831) * u + 0.956464974744799006) * u + 
        0.686265948274097816) * u + 0.434397492331430115) * u + 
        0.244044510593190935) * t - 
        z * exp(x * x - 0.120782237635245222);
    x += s * (x * s + 1);
    if (y > 1) {
        x = -x;
    }
    return x;
}

 
/****************************************************************************/

int main(int argc, char *argv[])
{
    int s,i,j,k;
    Ifloat DataB;
    //double e=0.00027;
     
     /* support image creation */
    fitsstruct Header;
    char Cmd[4096];
 
 
    //cout << "T 0.001 ==> ", << sqrt(2.) *ABS(dierfc( 0.001 )) << endl;;
    //cout << "T 0.0001 ==> ", << sqrt(2.) *ABS(dierfc( 0.0001 )) << endl;
    //cout << "T 0.00001 ==> ", << sqrt(2.) *ABS(dierfc( 0.00001 )) << endl;
 
 
    Cmd[0] = '\0';
    for (k =0; k < argc; k++) sprintf(Cmd, "%s %s", Cmd, argv[k]);

     /* Get command line arguments, open input file(s) if necessary */
#ifndef MRPOL  
    lm_check(LIC_MR1);
#else
    lm_check(LIC_POL);
#endif
    filtinit(argc, argv);
if (Verbose == True)
{
    cout << endl << endl << "PARAMETERS: " << endl << endl;
    cout << "File Name in = " << Name_Imag_In << endl;
    cout << "File Name Out = " << Name_Imag_Out << endl;
    cout << "Transform = " << StringTransform(Transform) << endl;
    cout << "Number of scales = " << Nbr_Plan << endl;
    cout << "Noise Type = " << StringNoise(Stat_Noise) << endl;
    if ((Transform == TO_MALLAT) || (Transform == TO_UNDECIMATED_MALLAT))
    {
       cout << StringSBFilter(SB_Filter) << endl;
       if (Norm == NORM_L2) cout << "L2 normalization" << endl;
       if (Transform == TO_UNDECIMATED_MALLAT) cout << "Number of undecimated scales = " <<  NbrUndec << endl;
    }
    if (Transform == TO_UNDECIMATED_NON_ORTHO)  cout << "Undec. Filter Bank = " << StringUndecFilter(U_Filter) << endl;
    if (Stat_Noise == NOISE_GAUSS_POISSON)
    {
           cout << "Type of Noise = POISSON" << endl;
           cout << "  Gain = " << PasCodeur << endl;
           cout << "  Read-out Noise Sigma  = " << SigmaGauss << endl;
           cout << "  Read-out Mean = " << MeanGauss << endl;
    }
    if ((Stat_Noise == NOISE_GAUSSIAN) && (Noise_Ima > FLOAT_EPSILON))
                cout << "Sigma Noise = " << Noise_Ima << endl;
    if (Stat_Noise ==  NOISE_EVENT_POISSON) 
        cout << "Epsilon Poisson = " <<  EpsilonPoisson << endl;
    if  (Stat_Noise == NOISE_SPECKLE)
        cout << "Epsilon speckle = " <<  EpsilonPoisson << endl;
    if ((Stat_Noise !=  NOISE_EVENT_POISSON)
         && (Stat_Noise != NOISE_SPECKLE))
    cout << "N_Sigma = " << N_Sigma << endl;
    
    cout << "Filter = " << StringFilter(Filter) << endl;
    cout << "Epsilon = " << Epsilon << endl;
    cout << "Max_Iter = " << Max_Iter << endl;
    if (PositivIma == True) cout << "Positivity constraint" << endl;
    else cout << "No positivity constraint" << endl;    
    if (MaxIma == True) cout << "Maximum level constraint" << endl;
    else cout << "No maximum level constraint" << endl;
    if (KeepPositivSup == True) 
       cout << "Only positive wavelet coefficients are detected" << endl;
    if (FirstScale > 0)
       cout << "Start the detect at scale " << FirstScale+1 << endl;
    if (KillLastScale == True)
       cout << "Suppress the last scale" << endl;
       
    if (WriteSup == True)
      cout << "Support Image file name : " << Name_Write_Sup << endl;
}
 
    io_read_ima_float(Name_Imag_In, DataB, &Header);
    Header.origin = Cmd; 
    Header.bitpix = BP_FLOAT;
    Ifloat Result (DataB.nl(), DataB.nc(), (char *) "Result Filtering");
    check_scale(DataB.nl(),  DataB.nc(), Nbr_Plan);

    // noise model class initialization
    FilterAnaSynt FAS;
    FilterAnaSynt *PtrFAS = NULL;
    if ((Transform == TO_MALLAT) || (Transform == TO_UNDECIMATED_MALLAT))
    {
        FAS.Verbose = Verbose;
        FAS.alloc(SB_Filter);
        PtrFAS = &FAS;
    }
    MRNoiseModel ModelData;
    ModelData.set_old_poisson( OldPoisson );
    ModelData.write_threshold( WriteThreshold ); 
    ModelData.alloc(Stat_Noise, DataB.nl(), DataB.nc(), 
                    Nbr_Plan, Transform, PtrFAS, Norm, NbrUndec);

    int NbrBand = ModelData.nbr_band();
    if (Noise_Ima > FLOAT_EPSILON) ModelData.SigmaNoise = Noise_Ima;
    if (UseNSigma  == True)
    {
    	if(Verbose) cout << " NbrTabNSigma = " << NbrTabNSigma << endl;
    	 if (NbrTabNSigma == 0) for (i=0; i < NbrBand; i++) ModelData.NSigma[i]=N_Sigma;
        else 
        {
        	  for (i=0; i < NbrTabNSigma ; i++)   ModelData.NSigma[i]=TabNSigma[i];
        	  for (i=NbrTabNSigma ; i < NbrBand; i++)   ModelData.NSigma[i]=TabNSigma[NbrTabNSigma-1];
        }
    }
         
    if (N_Sigma  == 111)
    {
        for (i=0; i < 4; i++)  ModelData.NSigma[i] = 5;
	ModelData.NSigma[4] = 4.;
	ModelData.NSigma[5] = 3.5;
 	for (i=6; i < NbrBand; i++) ModelData.NSigma[i] = 3.;
    }
    // for (i=0; i < NbrBand; i++) cout << " Scale " << i+1 << " NSig = " <<  ModelData.NSigma[i] << endl;
//     ModelData.NSigma[0] += 1.;
//     ModelData.NSigma[1] += 1.;
//     ModelData.NSigma[2] += 1.;
    if (SupIsol == True) ModelData.SupIsol = True;
    if (SupDil == True) ModelData.DilateSupport = True;
     for (s=0; s <  NbrBand; s++) ModelData.TabEps[s] = EpsilonPoisson;
    ModelData.OnlyPositivDetect = KeepPositivSup;
    ModelData.NiterSigmaClip = NiterClip;
    ModelData.SizeBlockSigmaNoise = SizeBlock;
    ModelData.FirstDectectScale = FirstScale;
    ModelData.CCD_Gain = PasCodeur;
    ModelData.CCD_ReadOutSigma = SigmaGauss;
    ModelData.CCD_ReadOutMean = MeanGauss;
    ModelData.GetEgde = GetEgde;
    ModelData.PoissonFisz = PoissonFisz;
    ModelData.TypeThreshold=TypeThreshold;
    ModelData.U_Filter=U_Filter;
    ModelData.MadEstimFromCenter = MadCenter;
    if (MinEventNumber > 0) ModelData.MinEventNumber = MinEventNumber;
    
    if (UseRMSMap == True) 
    {
       ModelData.UseRmsMap = True;
       io_read_ima_float(Name_RMSMap, ModelData.RmsMap);
    }


 
//    MultiResol MR_Data(DataB.nl(), DataB.nc(), Nbr_Plan, Transform, "MRNoiseModel");
//    ModelData.model(DataB, MR_Data);
//    INFO(MR_Data.band(0), "band 0");
//    ModelData.threshold(MR_Data);
//    INFO(MR_Data.band(0), "band 0");
//    ModelData.write_support_mr("sup.mr");
//    if (ModelData.signif(1000., 0, 100, 100) == True) cout << "OK" << endl;
//    else cout << "KO" << endl;
   
    // Filtering class initialization
    MRFiltering CFilter(ModelData, Filter);
    if (BackgroundImage == True)
    {
       io_read_ima_float(Name_Background, CFilter.BackgroundData);
       if (( CFilter.BackgroundData.nl() != DataB.nl()) || 
           ( CFilter.BackgroundData.nc() != DataB.nc()))
       {
          cout << "Error: the background image must have the same size as the input data ... " << endl;
	      exit(-1);
       }
       CFilter.BackgroundImage = True;
    }
        
    if (MissingData == True)
    {
        CFilter.MissingData = MissingData;
        (CFilter.MaskIma).alloc(DataB.nl(), DataB.nc(), "mask");
        for (i=0; i < DataB.nl(); i++) 
        for (j=0; j < DataB.nc(); j++) (CFilter.MaskIma)(i,j) = (DataB(i,j) == 0) ? 0: 1;
        BadPixel = True;
        BadPixalVal = 0.;
    }
    
    
    if (MaxIterOpt == True)  CFilter.Max_Iter = Max_Iter;
    if (KillLastOpt == True) CFilter.KillLastScale = KillLastScale;
    if (EpsOpt == True) CFilter.Epsilon = Epsilon;
    if (PosOpt == True) CFilter.PositivIma = PositivIma;
    if (MaxIma == True) CFilter.MaxIma = True;
    CFilter.Verbose = Verbose;
    if ((OptRegul == True) || (Filter == FILTER_TV_CONSTRAINT)) CFilter.RegulParam = RegulParam;
    
    CFilter.UseFlatField = UseFlatField;
    if (UseFlatField == True)  io_read_ima_float(Name_FlatField, CFilter.FlatField);
    if (PositiveReconsFilter == True) CFilter.Sup_Set = False;
    CFilter.filter(DataB, Result,PositiveReconsFilter);
    if (Stat_Noise == NOISE_EVENT_POISSON) Header.bitpix = BP_FLOAT;
    io_write_ima_float(Name_Imag_Out, Result, &Header);

    // support image creation 
    if ((WriteSup == True) && (Nbr_Plan > 1))
    {
        ModelData.write_support_mr(Name_Write_Sup);
        //ModelData.write_support_ima(Name_Write_Sup);
    }    
    if ((Verbose == True) && (Stat_Noise == NOISE_GAUSSIAN)  && (ABS(Noise_Ima) < FLOAT_EPSILON))
         cout << "Noise standard deviation = " <<  ModelData.SigmaNoise << endl;





    if (ProbMR == True)
     {
        int b,i,j;
	
	if( WriteInfoOnProbMap == True ) ModelData.write_in_few_event( True ) ;
	   
	
       //  cout << "ProbMR " << endl;
       MultiResol MR_Data;
       MR_Data.alloc(DataB.nl(), DataB.nc(),Nbr_Plan,ModelData.type_trans(), 
                 ModelData.filter_bank(), ModelData.TypeNorm, ModelData.nbr_undec_scale(), ModelData.U_Filter);
   
       if (ModelData.TransImag==True)  ModelData.im_transform(DataB);
       int LastScale=MR_Data.nbr_band()-1;
       if (BackgroundImage == True)
       {
           DataB -= CFilter.BackgroundData;
           MR_Data.transform(DataB);
       } 
       else MR_Data.transform(DataB);  
       if (KeepPositivSup == True)
       {
          for (b=0; b< LastScale; b++) 
          for (i=0; i< MR_Data.size_band_nl(b); i++)
          for (j=0; j< MR_Data.size_band_nc(b); j++)
             if (MR_Data(b,i,j) < 0) MR_Data(b,i,j) = 0.;
       }
       
       //ModelData.prob_noise(MR_Data,  True);
       //ModelData.threshold(MR_Data,  False);
       //MR_Data.band(LastScale).init();
  
 //       INFO(MR_Data.band(1), "scale 1");     
//       MR_Data.write("xx_prog.mr");
//        cout << "xerf(0.0027) = " << ABS(xerfc(0.5+(1-0.0027)/2.)) << endl;
//        cout << "xerf(0.0027) = " << ABS(dierfc(0.5+(1-0.0027)/2.)) << endl;
//        cout << "xerf(0.0027) = " << sqrt(2.) * ABS(dierfc((0.0027))) << endl;
//        cout << "xerf(0.00) = " << ABS(dierfc((0.0))) << endl;
//        cout << "xerf(1.) = " << ABS(dierfc((1.0))) << endl;
//         cout << "xerf(0.0027) = " << ABS(xerfc(0.5+(1-0.05)/2.)) << " " << sqrt(2.) * ABS(dierfc((0.05))) <<  endl;
//        cout << "xerf(0.0001) = " <<  sqrt(2.) * ABS(dierfc((0.0001))) << endl;
//        cout << "xerf(0.00001) = " <<  sqrt(2.) * ABS(dierfc((0.00001))) << endl;
//        cout << "xerf(0.000001) = " <<  sqrt(2.) * ABS(dierfc((0.000001))) << endl;
//        cout << "xerf(0.0000001) = " <<  sqrt(2.) * ABS(dierfc((0.0000001))) << endl;
//         cout << "xerf(1e-9) = " <<  sqrt(2.) * ABS(dierfc((1e-9))) << endl;
// 	cout << "xerf(1e-10) = " <<  sqrt(2.) * ABS(dierfc((1e-10))) << endl;
// 	cout << "xerf(1e-15) = " <<  sqrt(2.) * ABS(dierfc((1e-15))) << endl;
// 	cout << "xerf(1e-20) = " <<  sqrt(2.) * ABS(dierfc((1e-20))) << endl;
// 	cout << "xerf(1e-25) = " <<  sqrt(2.) * ABS(dierfc((1e-25))) << endl;
// 	cout << "xerf(1e-30) = " <<  sqrt(2.) * ABS(dierfc((1e-30))) << endl;
       for (b=0; b< LastScale; b++) 
       for (i=0; i< MR_Data.size_band_nl(b); i++)
       for (j=0; j< MR_Data.size_band_nc(b); j++)
       {
           //double Coef = 1. - MR_Data(b,i,j);
           //std::cout << setprecision( 16 ) << "prob(" << b << "," << i 
	   //          << "," << j << ")=" << Coef << std::endl;

           double Coef = ModelData.prob_signal_few_event( MR_Data(b,i,j), 
	                                                  b, i, j ) ;
           //std::cout << setprecision( 16 ) << "prob(" << b << "," << i 
	   //          << "," << j << ")=" << prob << std::endl;

           if (b < LastScale) 
	   {
 	      if (Coef < 1e-35) 
	         Coef = (float) (sqrt((double) 2.) * fabs(dierfc((1e-35))));
 	      else Coef = sqrt(2.) *ABS(dierfc( Coef ));
	      MR_Data(b,i,j) = Coef;
 	   }
	   if (MR_Data(b,i,j) > MR_Data(LastScale,i,j)) MR_Data(LastScale,i,j) = MR_Data(b,i,j);
       }
       // INFO(MR_Data.band(1), "scale 1");     
       MR_Data.write(ProbMRFile);
        
     }
    exit(0);
}



