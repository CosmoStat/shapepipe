/******************************************************************************
**                   Copyright (C) 1994 by CEA
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
**    File:  mr_transform.cc
**
*******************************************************************************
**
**    DESCRIPTION  multiresolution transform  an image
**    ----------- 
**                 
**    PARAMETRES    
**    ---------- 
**
**    USAGE: mr_transform option image output
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
**                 14: Mallat's wavelet transform 
**                 15: G transform (morphological min-max algorithm 
**                 16: Feauveau's wavelet transform 
**                 17: Haar's wavelet transform 
**                 18: Feauveau's wavelet transform without undersampling 
**
**                 default is 2
**                 17 is not yet implemented
**   
**           [-n number_of_scales]
**                number of scales used in the multiresolution transform
**                default is 4
**
**           [-x] 
**                write all scales separately as images with prefix "scale_j"
**                (j being the scale number)
**
**           [-b] 
**                same as x option, but interpolate by block the scales.
**                This option is valid only if the choosen multiresolution 
**                transform is pyramidal (6,7,8,9,10,11,12)
**
**           [-i] 
**                same as x option, but interpolate by B3-spline the scales.
**                This option is valid only if the choosen multiresolution 
**                transform is pyramidal (6,7,8,9,10,11,12)
**
**           [-c iter] 
**                iterative transformation. Iter = number of iterations.
**                This option is valid only if the choosen multiresolution 
**                transform is pyramidal (6,7,8,9,10,11). The reconstruction
**                is not exact and we need few iterations. Generally 3
**                iterations are enough.
**
**
******************************************************************************/
 
#include "IM_Obj.h"
#include "IM_IO.h"
#include "MR_Obj.h"
// #include "CPU_Time.h"
#include "IM_Prob.h"

char NamePre[256];
char Name_Imag_In[256]; /* input file image */
char Name_Imag_Out[256]; /* output file name */
int Nbr_Plan=DEFAULT_NBR_SCALE;  /* number of scales */
type_transform Transform = DEFAULT_TRANSFORM; /* type of transform */
type_lift LiftingTrans = DEF_LIFT;

Bool WriteScalex = False;          /* write each scale separately on the disk 
                                      with x option */
Bool WriteScaleb = False;          /* write each scale separately on the disk 
                                      with b option */
Bool WriteScalei = False;          /* write each scale separately on the disk 
                                      with i option */
int Iterc = 0;                     /* Iterative option */
int NbrUndec = -1;                     /*number of undecimated scale */

extern int  OptInd;
extern char *OptArg;

extern int  GetOpt(int argc, char *const*argv, char *opts);
Bool Verbose = False;
 
sb_type_norm Norm = NORM_L1;
type_sb_filter SB_Filter = F_MALLAT_7_9;
type_border Bord = I_CONT;
Bool OptS = False;
type_undec_filter U_Filter = DEF_UNDER_FILTER;

#define TR_DEBUG 0 

/*********************************************************************/

static void usage(char *argv[])
{
    fprintf(OUTMAN, "Usage: %s options in_image out_mr_file\n\n", argv[0]);
    fprintf(OUTMAN, "   where options =  \n");

    manline();
    all_transform_usage(Transform);
    manline();
    
    lifting_usage(LiftingTrans);
    manline();
    
    nbr_scale_usage(Nbr_Plan);
    manline();
    
    write_scales_x_band_usage();
    manline();
    
    write_scales_b_usage();
    manline();
 
    iter_transform_usage();
    manline();
    sb_usage(SB_Filter);
    manline();
    usb_usage(U_Filter);
    manline();
    nbr_nbr_undec_usage(NbrUndec);
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
static void transinit(int argc, char *argv[])
{
    int c;
    Bool OptL = False, Optf = False;
#ifdef LARGE_BUFF
    int VMSSize=-1;
    Bool OptZ = False;
    char VMSName[1024] = "";
#endif  
    
    /* get options */
    while ((c = GetOpt(argc,argv,"U:u:t:n:c:xBl:T:LvzZ:S")) != -1) 
    {
	switch (c) 
        {
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
           case 'u':
 		if (sscanf(OptArg,"%d",&NbrUndec) != 1) 
                {
		    fprintf(OUTMAN, "Error: bad number of scales: %s\n", OptArg);
	            exit(-1);
                    
		}
  		break;
           case 'S': OptS = True; break;
	   case 'L': Norm = NORM_L2; OptL = True; break;
	   case 'v': Verbose = True; break;
	   case 't':
		/* -d <type> type of transform */
		if (sscanf(OptArg,"%d",&c ) != 1) 
                {
		    fprintf(OUTMAN, "bad type of multiresolution transform: %s\n", OptArg);
	            exit(-1);
                    
		}
                if ((c > 0) && (c <= NBR_TOT_TRANSFORM+1)) 
                                        Transform = (type_transform) (c-1);
                else  
                {
		    fprintf(OUTMAN, "bad type of transform: %s\n", OptArg);
	            exit(-1);
 		}                
 		break;
	   case 'l':
 		if (sscanf(OptArg,"%d",&c ) != 1) 
                {
		    fprintf(OUTMAN, "Error: bad lifting method: %s\n", OptArg);
	            exit(-1);
                    
		}
                if ((c > 0) && (c <= NBR_LIFT)) LiftingTrans = (type_lift) (c);
                else  
                {
		    fprintf(OUTMAN, "Error: bad lifting method: %s\n", OptArg);
	            exit(-1);
 		}
  		break;
	   case 'T': 
		Optf = True;
		SB_Filter = get_filter_bank(OptArg);
		break;
	   case 'n':
		/* -n <Nbr_Plan> */
		if (sscanf(OptArg,"%d",&Nbr_Plan) != 1) 
                {
		    fprintf(OUTMAN, "bad number of scales: %s\n", OptArg);
		    exit(-1);
		}
                if ((Nbr_Plan <= 1) || (Nbr_Plan > MAX_SCALE)) 
                 {
		    fprintf(OUTMAN, "bad number of scales: %s\n", OptArg);
		    fprintf(OUTMAN, "1 < Nbr Scales <= %d\n", MAX_SCALE);
 		    exit(-1);
		}
		break;
 	   case 'c':
		/* -c < nbr of iterations> */
		if (sscanf(OptArg,"%d",&Iterc) != 1) 
                {
		    fprintf(OUTMAN, "bad number of iterations: %s\n", OptArg);
		    exit(-1);
		}
                if ((Iterc <= 1) || (Iterc > 20)) 
                 {
		    fprintf(OUTMAN, "bad number of iterations: %s\n", OptArg);
		    fprintf(OUTMAN, "1 < Nbr Iter <= 20\n");
 		    exit(-1);
		}
		break;
            case 'x':
                /* write scales */
                WriteScalex = True;
                if (WriteScaleb == True)
                {
		    fprintf(OUTMAN, "Error: options -x -B are exclusive ...\n");
 		    exit(-1);
		}

                break;
            case 'B':
                /* write scales */
                WriteScaleb = True;
                if (WriteScalex == True)
                 {
		    fprintf(OUTMAN, "Error: options -x -B are exclusive ...\n");
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
		usage(argv);
	}

        if ((Iterc > 1) && (SetTransform(Transform) != TRANSF_PYR))
        {
           fprintf(OUTMAN, "Error: option -c is only available with pyramidal transform ... \n");
            exit(0);
	}
	
	if ((Transform != TO_LIFTING) && (LiftingTrans != DEF_LIFT))
	{
           fprintf(OUTMAN, "Error: option -l is only available with lifting transform ... \n");
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

int main(int argc, char *argv[])
{
    int s;
    Ifloat Dat;
    char Prefix[256];
    char Name_Imag[256];

    /* Get command line arguments, open input file(s) if necessary */
    lm_check(LIC_MR1);
    transinit(argc, argv);

if (Verbose == True)
{
    cout << endl << endl << "PARAMETERS: " << endl << endl;
    cout << "File Name in = " << Name_Imag_In << endl;
    cout << "File Name Out = " << Name_Imag_Out << endl;
    cout << "Transform = " << StringTransform(Transform) << endl;
    if ((Transform == TO_MALLAT) || (Transform == TO_UNDECIMATED_MALLAT))
    {
       cout << StringSBFilter(SB_Filter) << endl;
       if (Norm == NORM_L2) cout << "L2 normalization" << endl;
    }
    if (Transform == TO_LIFTING)
       cout << StringLSTransform(LiftingTrans) << endl;
    cout << "Number of scales = " << Nbr_Plan << endl;
     if (Transform == TO_UNDECIMATED_MALLAT) cout << "Number of undecimated scales = " <<  NbrUndec << endl;
     if (Transform == TO_UNDECIMATED_NON_ORTHO)  cout << "Undec. Filter Bank = " << StringUndecFilter(U_Filter) << endl;
}

    io_read_ima_float(Name_Imag_In, Dat);
    // check_scale(Dat.nl(), Dat.nc(), Nbr_Plan);

    MultiResol MR_Data;
    FilterAnaSynt FAS;
    FilterAnaSynt *PtrFAS = NULL;
    if ((Transform == TO_MALLAT) || (Transform == TO_UNDECIMATED_MALLAT))
    {
        FAS.Verbose = Verbose;
        FAS.alloc(SB_Filter);
        PtrFAS = &FAS;
    }
    MR_Data.alloc (Dat.nl(), Dat.nc(), Nbr_Plan, Transform, PtrFAS, Norm, NbrUndec,U_Filter);
          
    if (Verbose == True)
        cout << "Number of bands = " << MR_Data.nbr_band()  << endl;
           
    if ((WriteScalei==True) && (MR_Data.Set_Transform != TRANSF_PYR))
    {
         cerr << "Error: illegal interpolation option " << endl;
         cerr << "       interpolation options are not allowed with a non-pyramidal transform. " << endl;
          exit (0);
    }
    if (Transform == TO_LIFTING) MR_Data.LiftingTrans = LiftingTrans;
//     if ((Transform == TO_MALLAT) || (Transform == TO_UNDECIMATED_MALLAT))
//     {
//         MR_Data.SBFilter = SB_Filter;
// 	MR_Data.TypeNorm = Norm;
//     }
    MR_Data.Border = Bord;
    MR_Data.Verbose = Verbose;

    // Perform the transformation
    {
      // CPUTime CPU;
      MR_Data.transform (Dat);
    }




    // write the results
    MR_Data.write (Name_Imag_Out);

    /* correct the reconstruction error */
    if ((Iterc > 1) && (MR_Data.Set_Transform == TRANSF_PYR))
                              mr_correct_pyr (Dat, MR_Data, Iterc);

    /* write separately each image */
    io_strcpy_prefix(Prefix,  Name_Imag_Out);
         
    if (WriteScaleb==True)
    {
        Ifloat Dats (MR_Data.size_ima_nl(), MR_Data.size_ima_nc(), "Interp");
        for (s = 0; s < MR_Data.nbr_band(); s++)
        {
            sprintf (Name_Imag, "%s_band_%d", Prefix, s+1);
            cout << "Write " << Name_Imag << endl;
            im_block_extend(MR_Data.band(s), Dats);
            io_write_ima_float(Name_Imag, Dats);
        }
    }
    else if (WriteScalex==True)
    {
        for (s = 0; s < MR_Data.nbr_band(); s++)
        {
            sprintf (Name_Imag, "%s_band_%d", Prefix, s+1);
            cout << "Write " << Name_Imag << endl;
            MR_Data.write_band(Name_Imag, s);
        }
    }


#if TR_DEBUG
//     for (int b = 0; b < MR_Data.nbr_band(); b+=2)
//     {
//        int  Nls = MR_Data.size_band_nl(b);
//        int Ncs = MR_Data.size_band_nc(b);
//             for (int i = 0; i < Nls;i++)
//             for (int j = 0; j < Ncs; j++)
//             {
//                 MR_Data(b,i,j) = sqrt(POW(MR_Data(b,i,j),2.)+POW(MR_Data(b+1,i,j), 2.));
//             }
//        INFO(MR_Data.band(b), "band ");
//      }
//      exit(-1);


    cout << "Reconstruction ... " << endl;
    Ifloat Rec(Dat.nl(), Dat.nc(), "rec");
 //    MR_Data.recons(Rec);
//     Dat -= Rec;
//     io_write_ima_float("x1.fits", Rec);
//     io_write_ima_float("x2.fits", Dat);
//     INFO(Dat,"Resi");
    MultiResol  MM;
    MM.read(Name_Imag_Out);
    MM.recons(Rec);
    io_write_ima_float("x1.fits", Rec);
    io_write_ima_float("x2.fits", Dat);
    Dat -= Rec;
    INFO(Dat,"Resi");
    
    if (OptS == True)
    {
       for (s = 0; s < MR_Data.nbr_band(); s++)
       {   
          double TMin, TMax;
          CImaProb CIP;
          CIP.Verbose = True;
          CIP.set(MR_Data.band(s));
          CIP.find_gthreshold(3., TMin, TMax);
          float Simu = MAX(ABS((float)TMin), ABS((float) TMax)) / 3.;
          cout << " Band " << s+1 << " Min = ";
          cout << TMin << " Max = " <<  TMax << " Sigma = " << Simu << " Norm = " << MR_Data.band_norm(s) << endl;
       }
    }
#endif

    exit(0);
} 

