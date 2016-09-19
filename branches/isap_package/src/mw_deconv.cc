/******************************************************************************
**                   Copyright (C) 1998 by CEA
*******************************************************************************
**
**    UNIT
**
**    Version: 1.0
**
**    Author: Jean-Luc Starck
**
**    Date:  98/01/01
**    
**    File:  mw_deconv.cc
**
*******************************************************************************
**
**    DESCRIPTION  deconvolution of an image by the multiscale entropy method 
**    ----------- 
**                 
**
******************************************************************************/
 
#include "IM_Obj.h"
#include "IM_Math.h"
#include "IM_IO.h"
#include "MR_Obj.h"
#include "IM_Deconv.h"
#include "MR_Sigma.h"
#include "MR_Filter.h"
#include "MR_NoiseModel.h"
#include "MR_Abaque.h"
#include "MR_Psupport.h"
#include "CErf.h"
#include "CMem.h"
#include "MW_Filter.h"
#include "MR_Deconv.h"
#include "MW_Deconv.h"

char Name_Imag_In[256];    // input file image 
char Name_Psf_In[256];     // PSF 
char Name_Imag_Out[256];   // output file name 
char Name_Imag_Start[256]; // First guess input solution 
char Name_Imag_ICF[256];   // ICF file name

int Nbr_Plan=DEFAULT_NBR_SCALE;  // number of scales 
float N_Sigma=DEFAULT_N_SIGMA;   // number of sigma (for the noise) 
float Noise_Ima=0.;              // noise standard deviation 
type_noise Stat_Noise = DEFAULT_STAT_NOISE;   // type of noise 
type_transform Transform = DEFAULT_TRANSFORM; // type of transform 
Bool WriteResi = False;         // write the residual 
Bool PsfMaxShift = True;        // shift the max of the PSF to the center of the image
char Name_Resi[256];            // residual file name 
extern float PasCodeur;  // CCD gain 
extern float SigmaGauss; // CCD read-out noise standard deviation 
extern float MeanGauss;  // CCD read-out noise mean 

Bool UseNSigma =False;
Bool GaussConv=False;
Bool UseICF=False;
Bool UseGuess=False;
Bool KeepPositivSup=False;

float RegulParam = 1.;
float Fwhm = 0.;
float Converg = 1.;       // convergence parameter 
float Epsilon = 1e-4;    // convergence parameter 
// float Epsilon=DEFAULT_EPSILON_DECONV;/* convergence parameter */
int Max_Iter=500;        // Maximum number of iteration 

extern int  OptInd;
extern char *OptArg;

extern int  GetOpt(int argc, char *const*argv, char *opts);
Bool Verbose;
Bool Optim = False;

type_deconv     DecMethod  = DEF_DEC_ENTROP;
type_mem_weight TypeWeight = DEF_DEC_WEIGHT;
type_mem_regul  TypeRegul  = DEF_DEC_REGUL;


Bool PositivIma = True;
float EpsilonPoisson = DEFAULT_EPSILON;
int SizeBlock = DEFAULT_SIZE_BLOCK_SIG;
int NiterClip = 1;
Bool KillLastScale =False;
char Name_RMSMap[256]; 
Bool UseRMSMap=False;
Bool GetMR_FirstGuess=False;
float NSigmaObj = 3.;

const int NBR_OK_METHOD = 2;
static type_deconv TabMethod[NBR_OK_METHOD] = {DEC_MR_MEM,DEC_MR_MEM_NOISE};
			    
Bool DataSNR = False;
Bool SupIsol = False;

/*********************************************************************/

static void usage(char *argv[])
{
    int i;
    fprintf(OUTMAN, "Usage: %s options in_image in_psf out_image\n\n", argv[0]);
    fprintf(OUTMAN, "   where options =  \n");
    manline();
    
    transform_usage(Transform);
    manline();
    
    fprintf(OUTMAN, "         [-H EntropyFunction]\n");
    fprintf(OUTMAN, "             1: H = Wavelet coef. Energy (N1-MSE).\n");
    fprintf(OUTMAN, "             2: H = Noise Information (N2-MSE, for Gaussian noise only).\n");
    fprintf(OUTMAN, "                Default is 1.\n");
    fprintf(OUTMAN, "                For Gaussian noise, default 2.\n");
    manline();

    fprintf(OUTMAN, "         [-m type_of_noise]\n");
    for (i = 0; i < NBR_NOISE-1; i++)
    fprintf(OUTMAN, "              %d: %s \n",i+1,
                                            StringNoise((type_noise )i));
    fprintf(OUTMAN, "             default is Gaussian noise\n");
    manline();
    
    gauss_usage();
    manline();  
     
    ccd_usage();
    manline();    
    
    nbr_scale_usage(Nbr_Plan);
    manline();
 
    nsigma_usage(N_Sigma);
    manline();
     
    max_iter_usage(Max_Iter);
    manline();

    converg_param_usage(Epsilon);
    manline();
    
    
    kill_last_scale_usage();
    manline();
              
    //size_block_usage(SizeBlock);
    //manline();
    
    //sigma_clip_block_usage(NiterClip);
    //manline();

    rms_noise_usage();
    manline();
    
    fprintf(OUTMAN, "         [-P]\n");
    fprintf(OUTMAN, "             Suppress the positivity constraint.\n");
    manline();
            
    fprintf(OUTMAN, "         [-f ICF_Fwhm]\n");
    fprintf(OUTMAN, "             Intrinsic correlation function.\n");
    fprintf(OUTMAN, "             Fwhm = Full Width at Half Maximum.\n");
    manline();
    
    fprintf(OUTMAN, "         [-I ICF_FileName]\n");
    fprintf(OUTMAN, "             Intrinsic correlation function file.\n");
    manline();
    
    fprintf(OUTMAN, "         [-F First_Guess]\n");
    fprintf(OUTMAN, "             Input solution file.\n");
    manline();
    
    fprintf(OUTMAN, "         [-W DataWeightingType]\n"); 
    fprintf(OUTMAN, "             1: no weighting  \n"); 
    fprintf(OUTMAN, "             2: soft weighting \n"); 
    fprintf(OUTMAN, "             3: hard weighting \n");
    fprintf(OUTMAN, "             default is %1d.\n", (int) DEF_DEC_WEIGHT+1);
    manline();
    
    fprintf(OUTMAN, "         [-A RegulProtectType]\n"); 
    fprintf(OUTMAN, "             0: no regularization (all protected) \n");
    fprintf(OUTMAN, "             1: no protection from regularization \n"); 
    fprintf(OUTMAN, "             2: soft protection \n"); 
    fprintf(OUTMAN, "             3: hard protection \n");
    fprintf(OUTMAN, "             4: soft + hard protection \n");
    fprintf(OUTMAN, "             default is %1d.\n", (int) DEF_DEC_REGUL);
    manline();       
    fprintf(OUTMAN, "         [-G RegulParam]\n");
    fprintf(OUTMAN, "             Regularization parameter \n");
    fprintf(OUTMAN, "             default is 1\n");
    fprintf(OUTMAN, "             with hard weighting (-W 3), default is 0.1\n");
    
    manline();       
    fprintf(OUTMAN, "         [-M NSigmaObj]\n");
    fprintf(OUTMAN, "             NSigma level for object multiresolution determination\n");
    fprintf(OUTMAN, "             default is 3\n");
    fprintf(OUTMAN, "             with hard weighting (-W 3), default is 1\n");
     
    manline();
    fprintf(OUTMAN, "         [-C ConvergParam]\n");
    fprintf(OUTMAN, "             Convergence parameter. \n");
    fprintf(OUTMAN, "             default is %f\n", Converg);
    manline();    
    
//     fprintf(OUTMAN, "         [-d]\n");
//     fprintf(OUTMAN, "             Wavelet-Waguelet deconvolution method. \n");
//     fprintf(OUTMAN, "             default is no.\n");
//     manline(); 
//     
//     fprintf(OUTMAN, "         [-D]\n");
//     fprintf(OUTMAN, "             Wavelet-Waguelet deconvolution method \n");
//     fprintf(OUTMAN, "             with protection of high sigma to noise ration coefficients\n");
//     fprintf(OUTMAN, "             default is no.\n");
//     manline();
    
    write_residual_usage();
    manline();
    
    psf_not_center_usage();
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
static void decinit(int argc, char *argv[])
{
    int c;    
    Bool UseOptEntrop = False;
    Bool UseOptG = False;
    Bool UseSigmaObj = False;
    
 #ifdef LARGE_BUFF
    int VMSSize=-1;
    Bool OptZ = False;
    char VMSName[1024] = "";
#endif  
   
    /* get options */
    while ((c = GetOpt(argc,argv,"pkVf:t:g:c:n:s:i:e:r:vC:G:OH:W:A:Pm:SKR:I:F:zZ:M:")) != -1) 
    {
	switch (c) 
        {
	   case 'k': SupIsol = True; break;
	   case 'p': KeepPositivSup = True; break;
	   case 'd': DecMethod = DEC_MR_VAGUELET; break;
	   case 'D': DecMethod = DEC_MR_VAGUELET; DataSNR = True; break;         
	   case 'V': GetMR_FirstGuess = True; break;
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
 	 case 'g':
		/* -g <sigma_noise> */
		if (sscanf(OptArg,"%f",&Noise_Ima) != 1) 
                {
		    fprintf(OUTMAN, "Error: bad sigma noise: %s\n", OptArg);
		    exit(-1);
		}
                Stat_Noise = NOISE_GAUSSIAN;
 		break;
          case 'f':
 		if (sscanf(OptArg,"%f",&Fwhm) != 1) 
                {
		    fprintf(OUTMAN, "Error: bad ICF Fwhm: %s\n", OptArg);
		    exit(-1);
		}
		if (Fwhm > 0) GaussConv = True;
		if (UseICF == True)
		{
		   fprintf(OUTMAN, "Error: -I and -f options are not compatible .. \n ");
		   exit(-1);
		}
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
		UseNSigma =True;
                if (N_Sigma <= 0.) N_Sigma = DEFAULT_N_SIGMA;
                 break;
	   case 'M':
		/* -s <nsigma> */
		if (sscanf(OptArg,"%f",&NSigmaObj) != 1) 
                {
		    fprintf(OUTMAN, "Error: bad N_Sigma: %s\n", OptArg);
		    exit(-1);
		}
		UseSigmaObj = True;
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
	   case 'C':
		/* -e < Convergence parameter> */
		if (sscanf(OptArg,"%f",&Converg) != 1) 
                {
		   fprintf(OUTMAN, "Error: bad convergence parameter: %s\n", OptArg);
		   exit(-1);
		}
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
           case 'S': PsfMaxShift = False; break;
           case 'v': Verbose = True;break;
 	   case 'O' : Optim = True; break;
	   case 'H' : 
	        if (sscanf(OptArg,"%d",&c) != 1) 
                {
		   fprintf(OUTMAN, "Error: bad entropy function: %s\n", OptArg);
		   exit(-1);
		}
 	        if ((c <= 0) || (c > NBR_OK_METHOD))
		{
                   fprintf(OUTMAN, "Error: bad entropy function ... \n");
	           exit(-1);
                    
		}      
 		UseOptEntrop = True;
		DecMethod = (type_deconv) (TabMethod[c-1]);
                break;
           case 'W':
 		if (sscanf(OptArg,"%d",&c) != 1) 
                {
		   fprintf(OUTMAN, "Error: bad weighting method: %s\n", OptArg);
		   exit(-1);
		}
                if ((c < 1) || (c > DEC_NBR_WEIGHT))
		{
		  fprintf(OUTMAN, "Error: bad weighting method: %s\n", OptArg);
		  exit(-1);
		}	
		TypeWeight = (type_mem_weight) (c-1);	
		break;
	   case 'A':
 		if (sscanf(OptArg,"%d",&c) != 1) 
                {
		   fprintf(OUTMAN, "Error: bad regularization method: %s\n", OptArg);
		   exit(-1);
		}
                if ((c < 0) || (c  >= DEC_NBR_REGUL))
		{
		  fprintf(OUTMAN, "Error: bad regularization method: %s\n", OptArg);
		  exit(-1);
		}	
		TypeRegul = (type_mem_regul) c;
		break;
	   case 'G':
		/* -s <nsigma> */
		if (sscanf(OptArg,"%f",&RegulParam) != 1) 
                {
		    fprintf(stderr, "bad Regularization Parameter: %s\n", OptArg);
		    exit(-1);
		}
		UseOptG = True;
                // if (RegulParam  < 0.) RegulParam = 1.;
		break;	   
	  case 'r':
		/* -r < residual file name> */
		if (sscanf(OptArg,"%s",Name_Resi) != 1) 
                {
		   fprintf(OUTMAN, "Error: bad file name: %s\n", OptArg);
		   exit(-1);
		}
                WriteResi = True;
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

	if (OptInd < argc) strcpy(Name_Psf_In, argv[OptInd++]);
         else usage(argv);

	if (OptInd < argc) strcpy(Name_Imag_Out, argv[OptInd++]);
         else usage(argv);

	/* make sure there are not too many parameters */
	if (OptInd < argc)
        {
		fprintf(OUTMAN, "Error: too many parameters: %s ...\n", argv[OptInd]);
		exit(-1);
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
           // Stat_Noise = NOISE_NON_UNI_ADD;
        }
        if ((isotrop(Transform) == False)
              && ((Stat_Noise == NOISE_NON_UNI_ADD) ||
                 (Stat_Noise  == NOISE_NON_UNI_MULT)))
        {
           cerr << endl << endl;
           cerr << "  Error: with this transform, non stationary noise models are not valid : " << StringFilter(FILTER_THRESHOLD) << endl;
           exit(-1);
        }
	
	
        if ((UseOptEntrop == False) && (DecMethod != DEC_MR_VAGUELET))
	{
	    if (Stat_Noise == NOISE_GAUSSIAN) DecMethod = DEF_DEC_ENTROP_GAUSS;
	    else DecMethod = DEF_DEC_ENTROP;
	}
	else
	{
	   if ((Stat_Noise != NOISE_GAUSSIAN) && (DecMethod ==  DEF_DEC_ENTROP_GAUSS))
	   {
	      cerr << "Error: this entropy function can only be used with Gaussian noise ... " << endl;
	   }
	}
	
	if ((TypeWeight == DEC_WEIGHT_SUPPORT) && (DecMethod != DEC_MR_VAGUELET))
	{
	   if (UseOptG == False) RegulParam = 0.1;
	   if (UseSigmaObj == False) NSigmaObj = 1;
	}
	
	
#ifdef LARGE_BUFF
    if (OptZ == True) vms_init(VMSSize, VMSName, Verbose);
#endif	    
}


/****************************************************************************/

main(int argc, char *argv[])
{
    Ifloat Data;
    int b,k;
    Ifloat Result;
    Ifloat Resi;
    Ifloat Psf, Psf1, Guess, Ima_ICF;
    Icomplex_f Psf_cf;
    fitsstruct Header;
    char Cmd[256];
    extern softinfo Soft;

    Soft.mr2();
    Cmd[0] = '\0';
    for (k =0; k < argc; k++) sprintf(Cmd, "%s %s", Cmd, argv[k]);
    lm_check(LIC_MR2);
    
     /* Get command line arguments, open input file(s) if necessary */
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
       cout << "Deconv = " << StringDeconv(DecMethod) << endl;
       cout << "N_Sigma = " << N_Sigma << endl;
       cout << "Epsilon = " << Epsilon << endl;
       cout << "Max_Iter = " << Max_Iter << endl;
       cout << "RegulParam = " << RegulParam << endl;
       cout << "Convergence paramter = " << Converg << endl;
       if (KillLastScale == True) cout << "Kill last scale " << endl;
       if (Fwhm > 0) cout << "Fwhm = " << Fwhm << endl;
       if (WriteResi == True)
          cout << "Image residual file name : " << Name_Resi << endl;
       if (TypeWeight == DEC_WEIGHT_PROBA)
              cout << "Soft data weighting " << endl;
       else if (TypeWeight == DEC_WEIGHT_SUPPORT)
              cout << "Hard data weighting " << endl;
       else cout << "No data weighting " << endl;
       
        switch(TypeRegul)
	{
	   case DEC_NO_REGUL:  cout << "No regul " << endl; break;
	   case DEC_REGUL_NO_PROTECT:  cout << "Regul " << endl; break;
	   case DEC_REGUL_PROBA:  cout << "Regul SNR" << endl; break;
	   case DEC_REGUL_SUPPORT:  cout << "Regul + SUPPORT" << endl; break;
 	   case DEC_REGUL_PROBA_SUPPORT: cout << "Regul SNR + SUPPORT" << endl; break;
           default:
	          cerr << "Error: unknown regularization method ... " << endl;
		  exit(-1);
		  break;
	 }     
    }

    /* read input image */
    MWDeconv CDec;
    io_read_ima_float(Name_Imag_In, CDec.Imag, &Header);
    
    io_read_ima_float(Name_Psf_In, CDec.Psf);
    if (UseGuess == True) io_read_ima_float(Name_Imag_Start, Guess);
    if (UseICF == True) io_read_ima_float(Name_Imag_ICF, Ima_ICF);
    
    Header.origin = Cmd;    
    CDec.TypeRegul = TypeRegul;
    CDec.TypeWeight = TypeWeight;
    CDec.RegulParam = RegulParam;
           
    CDec.KillLastScale = KillLastScale;
    CDec.PositivConstraint = PositivIma;
    CDec.DecMethod = DecMethod;
    CDec.PsfMaxShift = PsfMaxShift;
    CDec.Noise_Ima = Noise_Ima;
    CDec.MaxIter = Max_Iter;
    CDec.EpsCvg = Epsilon;
    CDec.IterCvg = Converg;
    CDec.GaussConv = GaussConv;
    CDec.Fwhm = Fwhm;
    CDec.OptimParam = Optim;
    CDec.GetMR_FirstGuess = GetMR_FirstGuess;
    CDec.NSigmaObj = NSigmaObj;
    
    CDec.Verbose = Verbose;
    Ifloat *Pt_G = NULL;
    if (UseGuess == True) Pt_G = &Guess;
    Ifloat *Pt_ICF = NULL;
    if (UseICF == True) Pt_ICF = &Ima_ICF;
   
    // noise model class initialization
    MRNoiseModel ModelData(Stat_Noise, CDec.Imag.nl(), CDec.Imag.nc(), 
                           Nbr_Plan, Transform);
    int NbrBand = ModelData.nbr_band();
    if (Noise_Ima > FLOAT_EPSILON) ModelData.SigmaNoise = Noise_Ima;

    if (UseNSigma  == True) 
    {
        EpsilonPoisson = (1. - erff((double) N_Sigma / sqrt((double) 2.)));
        for (b=0; b < NbrBand; b++) 
	{
	   ModelData.NSigma[b]=N_Sigma;
	   ModelData.TabEps[b] = EpsilonPoisson;
	}
    }
    else for (b=0; b < NbrBand; b++) ModelData.TabEps[b] = EpsilonPoisson;

    ModelData.NiterSigmaClip = NiterClip;
    ModelData.SizeBlockSigmaNoise = SizeBlock;
    ModelData.CCD_Gain = PasCodeur;
    ModelData.CCD_ReadOutSigma = SigmaGauss;
    ModelData.CCD_ReadOutMean = MeanGauss;
    ModelData.OnlyPositivDetect = KeepPositivSup;
    if (SupIsol == True) ModelData.SupIsol = True;
    if (UseRMSMap == True)
    {
       ModelData.UseRmsMap = True;
       io_read_ima_float(Name_RMSMap, ModelData.RmsMap);
    }
    CDec.ModelData = &ModelData;
     
    //DECONVOLUTION
    if (Verbose == TRUE) cout << " Start the deconvolution ... " << endl;
    if (CDec.DecMethod != DEC_MR_VAGUELET) CDec.mw_deconv(Pt_G, Pt_ICF);
    else CDec.mw_vaguelet(DataSNR, Pt_G, Pt_ICF);
    
    io_write_ima_float(Name_Imag_Out, CDec.Obj, &Header);
    if (WriteResi == True) io_write_ima_float(Name_Resi, CDec.Resi, &Header);   

    exit(0);
}

