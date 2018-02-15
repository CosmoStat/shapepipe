/******************************************************************************
**                   Copyright (C) 1995 by CEA
*******************************************************************************
**
**    UNIT
**
**    Version: 3.4
**
**    Author: Jean-Luc Starck
**
**    Date:  96/07/02
**    
**    File:  mr_deconv.cc
**
*******************************************************************************
**
**    DESCRIPTION  deconvolution of an image  
**    ----------- 
**                 
**    PARAMETRES    
**    ----------    
**
**    USAGE: mr_deconv option image psf output
**        where options = 
**          [-d type_of_deconvolution]
**               type_of_deconvolution = 
**                  1: Deconvolution by Van Cittert's algorithm
**                  2: Deconvolution by Gradient's algorithm
**                  3: Deconvolution by division in Fourier space
**                  4: Deconvolution by Lucy's algorithm
**                  5: Deconvolution by CLEAN algorithm
**                  6: Deconvolution by multiresolution Van Cittert's algorithm
**                  7: Deconvolution by multiresolution Gradient's algorithm
**                  8: Deconvolution by multiresolution Lucy's algorithm 
**                  9: Deconvolution by multiresolution CLEAN algorithm 
** 
**                  default is 8
**                  5 and 9 are not yet implemented
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
**           [-e Epsilon]
**                Convergence parameter
**                default is 0.0001
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
**           [-a Boolean_image]
**                add a Boolean image to the multiresolution support
**                Boolean_image = image file name
**                if we have information about positions of stars, ...,
**                we can add this position to the support
**                default is no addition
**
**           [-f Fwhm]
**                Fwhm = full width at half maximum
**                only used if type_of_deconvolution in [3,5,9]
**
**           [-w support_file_name]
**                if this option is set, an image is created from the 
**                multiresolution support. Default is no.
**
**           [-r residual_file_name]
**                if this option is set, the residual is written to 
**                the file of name residual_file_name. By default, the
**                residual is not written.
**
**           t,p,g,c,n,s,e,k,l,a,w options are only used if
**           type_of_deconvolution in [6..9]
**  
**
**
******************************************************************************/

#include "IM_Obj.h"
#include "IM_IO.h"
#include "MR_Obj.h"
#include "IM_Deconv.h"
#include "MR_Deconv.h"
#include "MR_Sigma.h"
#include "MR_Sigma.h"

char Name_Imag_In[1024];  /* input file image */
char Name_Psf_In[1024];   /* PSF */
char Name_Imag_Out[1024]; /* output file name */
char Name_Imag_Start[1024]; // First guess input solution 
char Name_Imag_ICF[1024];   // ICF file name

Bool WriteResi = False;         /* write the residual */
Bool PsfMaxShift = True;        /* shift the max of the PSF to the center of the image */
char Name_Resi[1024];  

Bool UseNSigma =False;
Bool GaussConv=False;
Bool UseICF=False;
Bool UseGuess=False;

float Fwhm = 0.;
float Converg = 1.;                      // convergence parameter 
int Max_Iter=DEFAULT_MAX_ITER_DECONV;    // Maximum number of iteration 

int Nbr_Plan=DEFAULT_NBR_SCALE;   /* number of scales */
float N_Sigma=DEFAULT_N_SIGMA;   /* number of sigma (for the noise) */
float Noise_Ima=0.;              /* noise standard deviation */
type_noise Stat_Noise = DEFAULT_STAT_NOISE;   /* type of noise */
type_transform Transform = DEFAULT_TRANSFORM; /* type of transform */
type_deconv Deconv = DEFAULT_DECONV; /* type of deconvolution */
float Epsilon=DEFAULT_EPSILON_DECONV;/* convergence parameter */
int SizeBlock = DEFAULT_SIZE_BLOCK_SIG;
int NiterClip = 1;
Bool KillLastScale =False;
char Name_RMSMap[1024]; 
Bool UseRMSMap=False;
Bool KeepPositivSup = False;
float RegulParam = 0.;
 
extern float PasCodeur;  /* CCD gain */
extern float SigmaGauss; /* CCD read-out noise standard deviation */
extern float MeanGauss;  /* CCD read-out noise mean */

extern int  OptInd;
extern char *OptArg;

extern int  GetOpt(int argc, char *const*argv, char *opts);

Bool Verbose;
Bool Optim = False;
Bool PositivIma = True;
Bool SupIsol = False;

const int NBR_OK_METHOD = 5;
static type_deconv TabMethod[NBR_OK_METHOD] = {DEC_MR_CITTERT,DEC_MR_GRADIENT,
                            DEC_MR_LUCY, DEC_MR_MAP,DEC_MR_VAGUELET};

sb_type_norm Norm = NORM_L1;
type_sb_filter SB_Filter = F_MALLAT_7_9;
int NbrUndec = -1;                     /*number of undecimated scale */

/*********************************************************************/

static void usage(char *argv[])
{
    int i;
    fprintf(OUTMAN, "Usage: %s options in_image in_psf out_image\n\n", argv[0]);
    fprintf(OUTMAN, "   where options =  \n");
    manline();
    
    fprintf(OUTMAN, "         [-d type_of_deconvolution]\n");
    for (i = 0; i < NBR_OK_METHOD; i++)
    fprintf(OUTMAN, "              %d: %s \n",i+1,StringDeconv((type_deconv) TabMethod[i]));
    fprintf(OUTMAN, "              default is %s\n", StringDeconv((type_deconv)  Deconv));
    manline();
   
    transform_usage(Transform);
    manline();

    fprintf(OUTMAN, "         [-T type_of_filters]\n");
    for (int i = 1; i <= NBR_SB_FILTER; i++)
    fprintf(OUTMAN, "              %d: %s \n",i,
                                           StringSBFilter((type_sb_filter  )i));
    fprintf(OUTMAN, "             default is %s\n\n", StringSBFilter ((type_sb_filter) SB_Filter));
    manline();

    nbr_nbr_undec_usage(NbrUndec);
    manline();
    
    gauss_usage();
    manline();
    
    ccd_usage();
    manline();
    
    fprintf(OUTMAN, "         [-m type_of_noise]\n");
    for (i = 0; i < NBR_NOISE-1; i++)
    fprintf(OUTMAN, "              %d: %s \n",i+1,
                                            StringNoise((type_noise )i));
    fprintf(OUTMAN, "             default is Gaussian noise\n");
    manline();
      
    nbr_scale_usage(Nbr_Plan);
    manline();
 
    nsigma_usage(N_Sigma);
    manline();
     
    max_iter_usage(Max_Iter);
    manline();
  
    converg_param_usage(Epsilon);
    manline();
    
    // sigma_clip_block_usage(NiterClip);
    // manline(); 
    
    rms_noise_usage();
    manline();
        
    fprintf(OUTMAN, "         [-f ICF_Fwhm]\n");
    fprintf(OUTMAN, "              Intrinsic correlation function.\n");
    fprintf(OUTMAN, "              Fwhm = Full Width at Half Maximum.\n");
    manline();
   
    fprintf(OUTMAN, "         [-P]\n");
    fprintf(OUTMAN, "             Suppress the positivity constraint.\n");
    manline();
    
    fprintf(OUTMAN, "         [-I ICF_FileName]\n");
    fprintf(OUTMAN, "              Intrinsic correlation function file.\n");
    manline();
    
    fprintf(OUTMAN, "         [-F First_Guess]\n");
    fprintf(OUTMAN, "              Input solution file.\n");
    manline();
    write_residual_usage();
    manline();    
    psf_not_center_usage();
    manline();
    detect_pos_usage();
    manline();
    kill_isol_pix_usage();
    manline();
    kill_last_scale_usage();
    manline();
    fprintf(OUTMAN, "         [-O]\n");
    fprintf(OUTMAN, "              Optimization.\n");
    manline();
    fprintf(OUTMAN, "         [-G RegulParam]\n");
    fprintf(OUTMAN, "              Regularization parameter \n");
    fprintf(OUTMAN, "              default is %f\n", RegulParam);
    manline();
    vm_usage();
    manline();    
    verbose_usage();
    manline();
    
    manline();
    manline();
    exit(-1);
}
 
/*********************************************************************/

/* GET COMMAND LINE ARGUMENTS */
static void decinit(int argc, char *argv[])
{
    Bool Optf=False;
    Bool OptC=False;
    int c;
#ifdef LARGE_BUFF
    int VMSSize=-1;
    Bool OptZ = False;
    char VMSName[1024] = "";
#endif    
    
    /* get options */
    while ((c = GetOpt(argc,argv,"C:G:ku:T:F:I:R:KPm:d:t:g:c:n:s:i:e:r:f:SvzZ:Op")) != -1) 
    {
	switch (c) 
        {
	   case 'C':
		/* -e < Convergence parameter> */
		if (sscanf(OptArg,"%f",&Converg) != 1) 
                {
		   fprintf(OUTMAN, "Error: bad convergence parameter: %s\n", OptArg);
		   exit(-1);
		}
		OptC = True;
 		break;
           case 'G':
		/* -s <nsigma> */
		if (sscanf(OptArg,"%f",&RegulParam) != 1) 
                {
		    fprintf(stderr, "bad Regularization Parameter: %s\n", OptArg);
		    exit(-1);
		}
                if (RegulParam  < 0.) RegulParam = 0.1;
 		break;	   
         case 'u':
 		if (sscanf(OptArg,"%d",&NbrUndec) != 1) 
                {
		    fprintf(OUTMAN, "Error: bad number of scales: %s\n", OptArg);
	            exit(-1);
                    
		}
  		break;
	  case 'F':
 		if (sscanf(OptArg,"%s",Name_Imag_Start) != 1) 
                {
		   fprintf(OUTMAN, "Error: bad file name: %s\n", OptArg);
		   exit(-1);
		}
                UseGuess = True;
 		break; 
	  case 'I':
 		if (sscanf(OptArg,"%s",Name_Imag_ICF) != 1) 
                {
		   fprintf(OUTMAN, "Error: bad file name: %s\n", OptArg);
		   exit(-1);
		}
                UseICF = True;
	        if (GaussConv == True)
		{
		   fprintf(OUTMAN, "Error: -I and -f options are not compatible .. \n ");
		   exit(-1);
		}
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
            case 'K': KillLastScale = True; break;
            case 'P': PositivIma=False; break; 
	    case 'p': KeepPositivSup = True;break; 
	    case 'k': SupIsol = True; break;   	   
	    case 'd':
		/* -d <type> type of deconvolution */
		if (sscanf(OptArg,"%d",&c ) != 1) 
                {
		    fprintf(OUTMAN, "Error: bad type of deconvolution: %s\n", OptArg);
	            exit(-1);
                    
		}        
		if ((c <= 0) || (c > NBR_OK_METHOD))
		{
		   fprintf(OUTMAN, "Error: bad type of deconvolution: %s\n", OptArg);
	           exit(-1);
                    
		}                             
		Deconv  = (type_deconv) (TabMethod[c-1]);
                break;
            case 'm':
		if (sscanf(OptArg,"%d",&c ) != 1) 
                {
		    fprintf(OUTMAN, "Error: bad type of noise: %s\n", OptArg);
	            exit(-1);
                    
		}
                if ((c > 0) && (c < NBR_NOISE)) 
                                        Stat_Noise = (type_noise) (c-1);
                else  
                {
		    fprintf(OUTMAN, "Error: bad type of noise: %s\n", OptArg);
	            exit(-1);
 		}
		break;
	   case 't':
		/* -d <type> type of transform */
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
                printf("OptArg = %s\n", OptArg);
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
 		break;
	   case 's':
		/* -s <nsigma> */
		if (sscanf(OptArg,"%f",&N_Sigma) != 1) 
                {
		    fprintf(OUTMAN, "Error: bad N_Sigma: %s\n", OptArg);
		    exit(-1);
		}
                if (N_Sigma <= 0.) N_Sigma = DEFAULT_N_SIGMA;
		UseNSigma =True;
                break;
	   case 'i':
		/* -i < Number of iterations> */
		if (sscanf(OptArg,"%d",&Max_Iter) != 1) 
                {
		   fprintf(OUTMAN, "Error: bad Max_Iter: %s\n", OptArg);
				exit(-1);
		}
                if (Max_Iter <= 0) Max_Iter = DEFAULT_MAX_ITER_DECONV;
		break;
	   case 'e':
		/* -e < Convergence parameter> */
		if (sscanf(OptArg,"%f",&Epsilon) != 1) 
                {
		   fprintf(OUTMAN, "Error: bad convergence parameter: %s\n", OptArg);
		   exit(-1);
		}
                if ((Epsilon < 0) || (Epsilon > 1.)) 
                           Epsilon = DEFAULT_EPSILON_DECONV;
		break;
            case 'S':
               PsfMaxShift = False;
               break;
           case 'v': Verbose = True;break;
 	   case 'O' : Optim = True; break;
	   case 'r':
		/* -r < residual file name> */
		if (sscanf(OptArg,"%s",Name_Resi) != 1) 
                {
		   fprintf(OUTMAN, "Error: bad file name: %s\n", OptArg);
		   exit(-1);
		}
                WriteResi = True;
 		break;
	   case 'f':
		/* -f < Fwhm parameter> */
		if (sscanf(OptArg,"%f",&Fwhm) != 1) 
                {
		   fprintf(OUTMAN, "Error: bad Fwhm: %s\n", OptArg);
		   exit(-1);
		} 
		if (Fwhm > 0) GaussConv = True;
		if (UseICF == True)
		{
		   fprintf(OUTMAN, "Error: -I and -f options are not compatible .. \n ");
		   exit(-1);
		}
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
  	   // Stat_Noise = NOISE_NON_UNI_ADD;
 	}
	if ((Stat_Noise ==  NOISE_CORREL) && (UseRMSMap != True))
	{
 	   cerr << endl << endl;
 	   cerr << "  Error: this noise model need a noise map (-R option) " << endl;
           exit(-1);
  	}
	
	if ((isotrop(Transform) == False)
 	      && ((Stat_Noise == NOISE_NON_UNI_ADD) ||
 	         (Stat_Noise  == NOISE_NON_UNI_MULT)))
  	{
 	   cerr << endl << endl;
 	   cerr << "  Error: with this transform, non stationary noise models are not valid : " << StringFilter(FILTER_THRESHOLD) << endl;
           exit(-1);
  	}
	if (KillLastScale == True)
	{
	   if (Optim == True)
	   {
	      cerr << endl << endl;
 	      cerr << "  Error: -K and -O options are not compatible ... " << endl;
	      exit(-1);
	   }
	   if ((Deconv == DEC_MR_LUCY) || (Deconv ==DEC_MR_MAP))
	   {
	      cerr << endl << endl;
 	      cerr << "  Error: -K option cannot be used with this deconvolution method  ... " << endl;
	      exit(-1);
	   }
  	}	
	
	if ((Max_Iter == DEFAULT_MAX_ITER_DECONV) && (Deconv == DEC_MR_VAGUELET))
	   Max_Iter = 10;
 
        if ((Transform != TO_UNDECIMATED_MALLAT) && (Transform != TO_MALLAT) && (Optf == True))
	{
	   fprintf(OUTMAN, "Error: option -T is only valid with Mallat transform ... \n");
           exit(0);
	}
        if ((OptC == False) && (RegulParam > 1.))
                             Converg = 1. / (2. * RegulParam);
			     
       /* get optional input file names from trailing 
          parameters and open files */
	if (OptInd < argc) strcpy(Name_Imag_In, argv[OptInd++]);
         else usage(argv);

	if (OptInd < argc) strcpy(Name_Psf_In, argv[OptInd++]);
         else usage(argv);

	if (OptInd < argc) strcpy(Name_Imag_Out, argv[OptInd++]);
         else usage(argv);

	/* make sure there are not too many parameters */
	if (OptInd < argc)
        {
		fprintf(OUTMAN, "Error: Too many parameters: %s ...\n", argv[OptInd]);
		exit(-1);
	}

#ifdef LARGE_BUFF
    if (OptZ == True) vms_init(VMSSize, VMSName, Verbose);
#endif
}

/*********************************************************************/

int main(int argc, char *argv[])
{
    int b, k;
    Ifloat Result;
    Ifloat Resi;
    fitsstruct Header;
    char Cmd[8192];
    Ifloat Guess, Ima_ICF;
   
    Cmd[0] = '\0';
    for (k =0; k < argc; k++) sprintf(Cmd, "%s %s", Cmd, argv[k]);
    
    /* Get command line arguments, open input file(s) if necessary */
    lm_check(LIC_MR1);
    decinit(argc, argv);

if (Verbose == True)
{
    cout << endl << endl << "PARAMETERS: " << endl << endl;
    cout << "File Name in = " << Name_Imag_In << endl;
    cout << "File Name Out = " << Name_Imag_Out << endl;
    cout << "Transform = " << StringTransform(Transform) << endl;
    cout << "Number of scales = " << Nbr_Plan << endl;
    if (Stat_Noise == NOISE_GAUSSIAN)
    {
        cout << "Type of Noise = GAUSSIAN" << endl;
	if (Noise_Ima > 0) cout << "Sigma Noise = " << Noise_Ima << endl;  
    }
    else
    {
           cout << "Type of Noise = POISSON" << endl;
           cout << "  Gain = " << PasCodeur << endl;
           cout << "  Read-out Noise Sigma  = " << SigmaGauss << endl;
           cout << "  Read-out Mean = " << MeanGauss << endl;
    }
    cout << "Deconv = " << StringDeconv(Deconv) << endl;
    cout << "N_Sigma = " << N_Sigma << endl;
    cout << "Epsilon = " << Epsilon << endl;
    cout << "Max_Iter = " << Max_Iter << endl;
    cout << "Convergence paramter = " << Converg << endl;
    if (KillLastScale == True) cout << "Kill last scale " << endl;
    cout << "Fwhm = " << Fwhm << endl;
    if (WriteResi == True)
      cout << "Image file name : " << Name_Resi << endl;
}

    /* read input image */
    MRDeconv CDec;
    io_read_ima_float(Name_Imag_In, CDec.Imag, &Header);
    io_read_ima_float(Name_Psf_In, CDec.Psf);
    Header.origin = Cmd;
    if (UseGuess == True) io_read_ima_float(Name_Imag_Start, Guess);
    if (UseICF == True) io_read_ima_float(Name_Imag_ICF, Ima_ICF);
    CDec.KillLastScale = KillLastScale;
    CDec.PositivConstraint = PositivIma;
    CDec.DecMethod = Deconv;
    CDec.PsfMaxShift = PsfMaxShift;
    CDec.Noise_Ima = Noise_Ima;
    CDec.MaxIter = Max_Iter;
    CDec.EpsCvg = Epsilon;
    CDec.IterCvg = Converg;
    CDec.GaussConv = GaussConv;
    CDec.Fwhm = Fwhm;
    CDec.OptimParam = Optim;
    CDec.Verbose = Verbose;
    CDec.RegulParam = RegulParam;
    Ifloat *Pt_G = NULL;
    if (UseGuess == True) Pt_G = &Guess;
    Ifloat *Pt_ICF = NULL;
    if (UseICF == True) Pt_ICF = &Ima_ICF;
    
    //DECONVOLUTION
    if (Verbose == TRUE) cout << " Start the deconvolution ... " << endl;
    CDec.StatNoise = Stat_Noise;

    // noise model class initialization
    // MRNoiseModel ModelData(Stat_Noise, CDec.Imag.nl(), CDec.Imag.nc(), 
    //                       Nbr_Plan, Transform);
    MRNoiseModel ModelData;
   // noise model class initialization
    FilterAnaSynt FAS;
    FilterAnaSynt *PtrFAS = NULL;
    if ((Transform == TO_MALLAT) || (Transform == TO_UNDECIMATED_MALLAT))
    {
        FAS.Verbose = Verbose;
        FAS.alloc(SB_Filter);
        PtrFAS = &FAS;
    }
    ModelData.alloc(Stat_Noise, CDec.Imag.nl(), CDec.Imag.nc(), 
                    Nbr_Plan, Transform, PtrFAS, Norm, NbrUndec);

    int NbrBand = ModelData.nbr_band();
    ModelData.OnlyPositivDetect = KeepPositivSup;
    if (Noise_Ima > FLOAT_EPSILON) ModelData.SigmaNoise = Noise_Ima;
    if (UseNSigma  == True)
          for (b=0; b < NbrBand; b++) ModelData.NSigma[b]=N_Sigma;
    ModelData.NiterSigmaClip = NiterClip;
    ModelData.SizeBlockSigmaNoise = SizeBlock;
    ModelData.CCD_Gain = PasCodeur;
    ModelData.CCD_ReadOutSigma = SigmaGauss;
    ModelData.CCD_ReadOutMean = MeanGauss;
    if (SupIsol == True) ModelData.SupIsol = True;
    if (UseRMSMap == True)
    {
       ModelData.UseRmsMap = True;
       io_read_ima_float(Name_RMSMap, ModelData.RmsMap);
    }
    CDec.ModelData = &ModelData;    
    
    CDec.im_deconv(Pt_G, Pt_ICF);

    io_write_ima_float(Name_Imag_Out, CDec.Obj, &Header);
    if (WriteResi == True) io_write_ima_float(Name_Resi, CDec.Resi, &Header);   
    exit(0);
} 

