/******************************************************************************
**                   Copyright (C) 2007 by CEA
*******************************************************************************
**
**    UNIT
**
**    Version: 1.0
**
**    Author: Jean-Luc Starck - Jerome Bobin
**
**    Date:  11/30/07
**    
**    File:  mr_gmca.cc
**
*******************************************************************************
**
**    DESCRIPTION  Generalized Morphological Component Analysis
**    ----------- 
**                 
**    Usage: mr_gmca options cube output
**
******************************************************************************/

#include "MR2D1D.h"
#include "DefMath.h"
#include "MatrixOper.h"  // defined in $TOOLS
#include <cmath>
#include "Array.h"
#include "NR.h"
#include "IM_Obj.h"
#include "GMCA.h"

 
/****************************************************************************/


class GMCA_2D: public MR2D1D, public GMCA 
{
    public:
    
    GMCA_2D():GMCA() {};
    void run(fltarray &TabCannels, fltarray & TabSource, fltarray & InvMixingMat);
    void inpainting_run(fltarray &TabCannels, fltarray & InpData);
    void transrecons_sources(fltarray &TabVect,fltarray &Recdata);
    void transform_sources(fltarray &Data,fltarray &TabVect);
   ~GMCA_2D() {} ;
};

/****************************************************************************/

void GMCA_2D::run(fltarray &TabCannels, fltarray & TabSource, fltarray & InvMixingMat)
{
cout << "TO BE DONE" << endl;
}

/****************************************************************************/

void GMCA_2D::transform_sources(fltarray &Data,fltarray &TabVect)
{
   // TabVect.alloc(nbr_coef_2d(), nbr_coef_1d());   // Pas necessaire
   int i,j,b,z;
   int Nxx = Data.nx();
   int Nyy = Data.ny();
   int Nzz = Data.nz();
   Ifloat Frame(Nyy, Nxx);
   fltarray Vect(Nzz);
   intarray TabPosBand(nbr_band_2d());
    
   // 2D wt transform per frame
   for (z=0; z < Nzz; z++)
   {
      int Pix=0;
      for (i=0; i < Nyy; i++)
      for (j=0; j < Nxx; j++) Frame(i,j) = Data(j,i,z);
      WT2D.transform(Frame);
      for (b=0; b < nbr_band_2d(); b++)
      {
           for (i=0; i < WT2D.size_band_nl(b); i++)
	       for (j=0; j < WT2D.size_band_nc(b); j++) TabVect(Pix++,z) = WT2D(b,i,j);
	       if (z == 0) TabPosBand(b) = Pix;
       }
   }


}

/****************************************************************************/

void GMCA_2D::transrecons_sources(fltarray &TabVect,fltarray &Recdata)
{
    
   int Nxx = nbr_pix_nx();
   int Nyy = nbr_pix_ny();
   int Nzz = TabVect.ny();
   //int Nsx = nbr_coef_2d();
      
   Recdata.alloc(Nxx,Nyy,Nzz);
   
   // RECONSTRUCTION
   
   int i,j,b,z;
   Ifloat Frame(Nyy, Nxx);
   intarray TabPosBand(nbr_band_2d());
    
   // 2D wt transform per frame
   for (z=0; z < Nzz; z++)
   {
      int Pix=0;

      for (b=0; b < nbr_band_2d(); b++)
      {
           for (i=0; i < WT2D.size_band_nl(b); i++)
	       for (j=0; j < WT2D.size_band_nc(b); j++) WT2D(b,i,j) = TabVect(Pix++,z);
	       if (z == 0) TabPosBand(b) = Pix;
       }
       WT2D.recons(Frame);
       for (i=0; i < Nyy; i++)
       for (j=0; j < Nxx; j++) Recdata(j,i,z) = Frame(i,j);
   }

    

}

/****************************************************************************/

void GMCA_2D::inpainting_run(fltarray &TabCannels, fltarray & InpData)
{
cout << "TO BE DONE" << endl;
}

/****************************************************************************/

char Name_Cube_In[256];
char Name_Out[256];
char Name_Mask[256];
char Name_KnowColumn[256];
char Name_Out_2[256];

int NbrScale2d = 5;
int Nbr_Plan=1;
type_transform  Transform=TO_MALLAT;

extern int  OptInd;
extern char *OptArg;

extern int  GetOpt(int argc, char *const*argv, char *opts);

// sb_type_norm Norm = NORM_L2;
// type_sb_filter SB_Filter = F_HAAR; // F_MALLAT_7_9;
// type_lift LiftingTrans = DEF_LIFT;
Bool Verbose=False;
Bool GetMax = False;
float N_Sigma=3;                // number of sigma (for the noise)  
Bool Normalize=False;       // normalize data in
 
type_border Bord = I_MIRROR;
Bool Reverse = False;

int NbrCannels=0;  // Number of cannels 
int NbrSources=0;  // Number of Sources
float KMin = 3; // Last K-Mad Thresholding
int Max_GMCA_Iter = DEF_GMCA_MAX_ITER; // Maximum number of iterations in the GMCA algorithm
Bool UseKnownColomn = False;
Bool UseMask = False;
int NbrKnownColumn = 0;   // Number of column which are known in the mixing matrix (a priori)
fltarray MatKnownColumn; // Matrix containing the known column
type_gmca_threshold_decrease TypeThresholdDecrease = GMCA_THRES_MAD; // Select the decreasing threshold strategy 
Bool Inpainting = False;            // If true, GMCA consider missing data, and apply an inpainting technique
fltarray Mask;              // if Inpainting==True, Mask contains the available data Mask(i) = 1 or 0 (0 for no data).
Bool PositiveMatrix = False;        // if True, the mixing matrice is assume to be positive
Bool PositiveSource = False;        // if True, the sources are assumed to be positive
Bool L1_Matrix_Constraint = False;  // if True, a sparsity constraint is applied on the matrix
Bool EstimNbSources = False;  // if True, jointly estimate the number of sources
Bool WriteChannels = False;
Bool DisjSpec = False;
Bool OrthoSpectra = False;
Bool GThrd = False;
float RelErrorMax = 40;
Bool UsePCA = False;
float MadStop = 0;
     
/*********************************************************************/

/*static int max_scale_number (int Nc)
{
    int Nmin = Nc;
    int ScaleMax;

    ScaleMax=iround(log((float)Nmin/(float)MAX_SIZE_LAST_SCALE) / log(2.)+ 1.);
    return (ScaleMax);
}*/

/*********************************************************************/
 
static void usage(char *argv[])
{
    fprintf(OUTMAN, "Usage: %s options input_data output_source [output_channel] \n\n", argv[0]);
    fprintf(OUTMAN, "   where options =  \n");

    
    fprintf(OUTMAN, "         [-i NbrIter]\n");
    fprintf(OUTMAN, "             Number of iterations. Default is no %d. \n", Max_GMCA_Iter);
    fprintf(OUTMAN, "         [-N number_of_scales_1D]\n");
    fprintf(OUTMAN, "              number of scales used in the 1D multiresolution transform.\n");
    all_transform_usage(Transform);
    manline();
    fprintf(OUTMAN, "         [-n number_of_scales_2D]\n");
    fprintf(OUTMAN, "              number of scales used in the 2D multiresolution transform\n");
    fprintf(OUTMAN, "              by default, no 1D transform is applied. \n");
    fprintf(OUTMAN, "         [-M]\n");
    fprintf(OUTMAN, "             Normalize the data. Default is no. \n");    
    fprintf(OUTMAN, "         [-S Number_of_Sources]\n");
    fprintf(OUTMAN, "             Number of sources. \n");    
    //  fprintf(OUTMAN, "             It is automatically estimated by default.\n");
    fprintf(OUTMAN, "         [-p]\n");
    fprintf(OUTMAN, "             Positivity contraint on the sources. Default is no. \n"); 
    fprintf(OUTMAN, "         [-U]\n");
    fprintf(OUTMAN, "             Use PCA to constrain the subspace spaned by the columns of MixingMat. Default is no. \n");    
    fprintf(OUTMAN, "         [-P]\n");
    fprintf(OUTMAN, "             Positivity contraint on the mixing matrix. Default is no. \n");        
    fprintf(OUTMAN, "         [-I MissingDataMaskFileName]\n");
    fprintf(OUTMAN, "             If set, the file name contains the position of the available data.  \n");  
    fprintf(OUTMAN, "             An inpainting teachnique is applied with GMCA. Default is no. \n");   
    fprintf(OUTMAN, "         [-A KnowColumnFileName]\n");
    fprintf(OUTMAN, "             File which contains the known column (a priori information).  \n");    
    fprintf(OUTMAN, "             By default, no a priori information \n");
    fprintf(OUTMAN, "         [-l]\n");
    fprintf(OUTMAN, "             Apply a l_1 constraint also on the mixing matrix.  Default is no. \n");  
    fprintf(OUTMAN, "         [-d]\n");
    fprintf(OUTMAN, "             Estimate the number of sources.  Default is no. \n"); 
    fprintf(OUTMAN, "         [-d]\n");
    fprintf(OUTMAN, "             Estimate the number of sources.  Default is no. \n"); 
    fprintf(OUTMAN, "         [-m]\n");
    fprintf(OUTMAN, "             Mad-based stopping criterion when the number of sources is estimated.  Default is 5 - default criterion is l2-based. \n"); 
    fprintf(OUTMAN, "         [-L]\n");
    fprintf(OUTMAN, "             L2-based stopping criterion when the number of sources is estimated.  Default is 40 (in dB). \n"); 
    fprintf(OUTMAN, "         [-D] \n");
    fprintf(OUTMAN, "             Spectra with disjoint supports for thresholds higher than 7 Mad.  \n");   // Should be an option
    fprintf(OUTMAN, "         [-K Last K-Mad]\n");
    fprintf(OUTMAN, "             Last value of K for K-Mad Thresholding.  \n");
    fprintf(OUTMAN, "         [-G Global Thresholding]\n");
    fprintf(OUTMAN, "         [-O]\n");        
    fprintf(OUTMAN, "                 Orthogonalization of the spectra\n");       
    verbose_usage();    
    vm_usage();
     manline();
    exit(-1);
}

/*********************************************************************/

/* GET COMMAND LINE ARGUMENTS */

static void transinit(int argc, char *argv[])
{
    int c;
#ifdef LARGE_BUFF
    int VMSSize=-1;
    Bool OptZ = False;
    char VMSName[1024] = "";
#endif   

    /* get options */
    while ((c = GetOpt(argc,argv,"i:N:t:n:MS:PpI:A:ldm:L:DK:GUOvzZ:")) != -1) 
    {
	switch (c) 
        {
	  case 'I':
		/* -n <Nbr_Plan> */
		if (sscanf(OptArg,"%s",Name_Mask) != 1) 
                {
		    fprintf(OUTMAN, "bad file name: %s\n", OptArg);
		    exit(-1);
		}
		UseMask = True;
 		break;
	  case 'A':
		/* -n <Nbr_Plan> */
		if (sscanf(OptArg,"%s", Name_KnowColumn) != 1) 
                {
		    fprintf(OUTMAN, "bad file name: %s\n", OptArg);
		    exit(-1);
		}
		UseKnownColomn = True;
 		break;
	  case 'l': L1_Matrix_Constraint = True; break;
	  case 'd': EstimNbSources = True; break;
	  case 'D': DisjSpec= True; break;
	  case 'p': PositiveSource  = True; break;
	  case 'P': PositiveMatrix = True; break;
	  case 'O': OrthoSpectra = True; break;
	  case 'G': GThrd = True; break;
	  case 'U': UsePCA = True; break;
	  case 'S':
		/* -n <Nbr_Plan> */
		if (sscanf(OptArg,"%d",&NbrSources) != 1) 
                {
		    fprintf(OUTMAN, "bad number of sources: %s\n", OptArg);
		    exit(-1);
		}
 		break;
          case 'i':
		/* -n <Nbr_Plan> */
		if (sscanf(OptArg,"%d",&Max_GMCA_Iter) != 1) 
                {
		    fprintf(OUTMAN, "bad number of iterations: %s\n", OptArg);
		    exit(-1);
		}
 		break;
 	  case 'm':
		/* -n <Nbr_Plan> */
		if (sscanf(OptArg,"%g",&MadStop) != 1) 
                {
		    fprintf(OUTMAN, "bad k-mad value: %s\n", OptArg);
		    exit(-1);
		}
 		break;
       case 'L':
		/* -n <Nbr_Plan> */
		if (sscanf(OptArg,"%g",&RelErrorMax) != 1) 
                {
		    fprintf(OUTMAN, "bad l2 maximal error in dB: %s\n", OptArg);
		    exit(-1);
		}
 		break;
	  case 't':
		/* -d <type> type of transform */
		if (sscanf(OptArg,"%d",&c ) != 1) 
                {
		    fprintf(OUTMAN, "bad type of multiresolution transform: %s\n", OptArg);
	            exit(-1);
                    
		}
                if ((c > 0) && (c <= NBR_TRANSFORM+1)) 
                                        Transform = (type_transform) (c-1);
                else  
                {
		    fprintf(OUTMAN, "bad type of transform: %s\n", OptArg);
	            exit(-1);
 		}                
 		break;
       case 'M': Normalize = (Normalize == True) ? False: True; break;
 	   case 'v': Verbose = True; break;
	   case 'n':
		/* -n <Nbr_Plan> */
		if (sscanf(OptArg,"%d",&NbrScale2d) != 1) 
                {
		    fprintf(OUTMAN, "bad number of scales: %s\n", OptArg);
		    exit(-1);
		}
                if ((NbrScale2d <= 1) || (NbrScale2d > MAX_SCALE_1D)) 
                 {
		    fprintf(OUTMAN, "bad number of scales: %s\n", OptArg);
  		    exit(-1);
		}
		break;
	    case 'N':
		/* -n <Nbr_Plan> */
		if (sscanf(OptArg,"%d",&Nbr_Plan) != 1) 
                {
		    fprintf(OUTMAN, "bad number of scales: %s\n", OptArg);
		    exit(-1);
		}
                if ((Nbr_Plan <= 0) || (Nbr_Plan > MAX_SCALE_1D)) 
                 {
		    fprintf(OUTMAN, "bad number of scales: %s\n", OptArg);
  		    exit(-1);
		}
		break;
		case 'K':
		/* -n <Nbr_Plan> */
		if (sscanf(OptArg,"%g",&KMin) < 0) 
                {
		    fprintf(OUTMAN, "bad value of last K \n");
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
 	   case '?':
			usage(argv);
		}
	}
      
	
       /* get optional input file names from trailing 
          parameters and open files */
	if (OptInd < argc) strcpy(Name_Cube_In, argv[OptInd++]);
         else usage(argv);

	if (OptInd < argc) strcpy(Name_Out, argv[OptInd++]);
        else usage(argv);

	if (OptInd < argc) 
	{
	   strcpy(Name_Out_2, argv[OptInd++]);
           WriteChannels = True;
	}
	
	/* make sure there are not too many parameters */
	if (OptInd < argc)
        {
		fprintf(OUTMAN, "Error: too many parameters: %s ...\n", argv[OptInd]);
		exit(-1);
	}
#ifdef LARGE_BUFF
    if (OptZ == True) vms_init(VMSSize, VMSName, Verbose);
#endif
}

/*********************************************************************/
  
int main(int argc, char *argv[]) 
{
    fltarray Dat;
    /* Get command line arguments, open input file(s) if necessary */
   fitsstruct Header;
   char Cmd[512];
   Cmd[0] = '\0';
   for (int k =0; k < argc; k++) sprintf(Cmd, "%s %s", Cmd, argv[k]);
   
    lm_check(LIC_MR3);
    transinit(argc, argv);
 
    if (Verbose == True)
    {
           cout << "Filename in = " << Name_Cube_In << endl;
           cout << "Filename out = " << Name_Out  << endl;
           cout << "Transform = " << StringTransform((type_transform) Transform) << endl;
           cout << "NbrScale2d = " << NbrScale2d<< endl;
           cout << "NbrScale1d = " << Nbr_Plan<< endl;
	       cout << "NbrGMCA_Iter = " <<  Max_GMCA_Iter << endl;
	       if (PositiveSource  == True) cout << "Positive source constraint " << endl;
	       if (PositiveMatrix  == True) cout << "Positive matrix constraint " << endl;
           if (L1_Matrix_Constraint  == True) cout << "l_1 norm constraint on the  matrix" << endl;
           if (UseKnownColomn == True) cout << "Use a priori information on the matrix. File name = " << Name_KnowColumn << endl;
	       if (UseMask == True) cout << "Use a mask for missing data. File name = " <<  Name_Mask << endl;
	       if (OrthoSpectra == True) cout << "Orthogonalization of the spectra" << endl;
	       if (DisjSpec == True) cout << " Spectra with disjoint supports "  << endl;
	       if (GThrd == True) cout << " Use Global Thresholding "  << endl;
	       if (UsePCA == True) cout << " Use PCA for subspace constraint "  << endl;
     }
	
       if (Verbose == True) cout << "\n Reading the data"<< endl;
   
       io_3d_read_data(Name_Cube_In, Dat, &Header);
 
       int Nx = Dat.nx();
       int Ny = Dat.ny();
       int Nz = Dat.nz();
       if (Verbose == True) cout << "Nx = " << Dat.nx() << " Ny = " << Dat.ny() << " Nz = " << Dat.nz() << endl;
     
       if (Normalize == True)
       {
         double Mean = Dat.mean();
         double Sigma = Dat.sigma();
         for (int i=0;i<Nx;i++)
         for (int j=0;j<Ny;j++)
         for (int k=0;k<Nz;k++) Dat(i,j,k) = (Dat(i,j,k)-Mean)/Sigma;
       }    
    
   
       // MR2D1D WT;
       GMCA_2D WT;
       
       if (Verbose == True) cout << "Alloc ...  " << endl;
       WT.alloc(Nx, Ny, Nz, Transform, NbrScale2d, 1); // On ne regularise que la MixingMat
 
       if (Verbose == True)	cout << "2d1d_trans ... "<< endl;
 
       //WT.transform (Dat);  // Pas utile
       
       // WT.write(Name_Out);
       // Compute the 2D1D transform
       fltarray TabVect;
       WT.transform_to_vectarray(Dat, TabVect);
       // fits_write_fltarr ("xx_tabvect.fits", TabVect);

       // Initalize the class for GMCA
       
       WT.NbrCannels = Nz;
       // cout << Nz << endl;
       if (NbrSources == 0) NbrSources = Nz;
       WT.NbrSources = NbrSources; 
       WT.Max_GMCA_Iter = Max_GMCA_Iter;
       WT.TypeThresholdDecrease = TypeThresholdDecrease;
       WT.Inpainting = UseMask;
       WT.PositiveMatrix = PositiveMatrix;
       WT.PositiveSource  = PositiveSource;
       WT.L1_Matrix_Constraint = L1_Matrix_Constraint;
       WT.EstimNbSources = EstimNbSources;
       WT.KMin = KMin;
       WT.DisjSpec = DisjSpec;
       WT.OrthoSpectra =OrthoSpectra;
       WT.GlobThrd = GThrd;
       WT.SVConst = UsePCA;
       WT.MatNbrScale1D = Nbr_Plan;
       fltarray QSVec;
       
       if (UsePCA == True)
       {
       	WT.pca(TabVect,QSVec); // Compute the Singular Vectors
       	WT.SingVec = QSVec;
       }
       
       if (UseMask == True)
       {
          fits_read_fltarr (Name_Mask, WT.Mask);
       }
       
       if (UseKnownColomn  == True)
       {
          fits_read_fltarr (Name_KnowColumn, WT.MatKnownColumn);
	      if (WT.MatKnownColumn.naxis() == 1) WT.NbrKnownColumn = 1;
	      else WT.NbrKnownColumn = WT.MatKnownColumn.axis(2);
       }
       
       //########## HERE THE NUMBER OF SOURCES MAY CHANGE run the GMCA method
       if (Verbose == True) cout << "Running GMCA ... "<< endl;
	   
	  if (EstimNbSources == False)  // THE NUMBER OF SOURCES IS FIXED
       {
       	fltarray TabSource;
       	int NbrCoef = TabVect.nx();
      	TabSource.alloc(NbrCoef,NbrSources);
       	WT.GMCA::Verbose = Verbose;
       	WT.run_gmca(TabVect, TabSource);

       }
       
       if (EstimNbSources == TRUE)  // THE NUMBER OF SOURCES IS ESTIMATED
       {
       	int NbrSourcesMax = NbrSources;
       	float RelError = 0;
       	//float SigmaData=TabVect.sigma();
       	NbrSources = 1;
       	bool ExitCriterion = False; /// CHANGEDDDDD
       	float OldRelError=0;
       	float DiffRelError;
       	
       	while (NbrSources <= NbrSourcesMax && ExitCriterion == False)
       	{
       		NbrSources++;
       		if (Verbose == True) cout << "Running GMCA ... Number of Estimated Sources : " << NbrSources << endl;
       		WT.NbrSources = NbrSources;
       		fltarray TabSource;
       		int NbrCoef = TabVect.nx();
	      	TabSource.alloc(NbrCoef,NbrSources);
       		WT.GMCA::Verbose = Verbose;
	       	WT.run_gmca(TabVect, TabSource);
	       	
	        WT.RetrieveSourcesFromMixmat(TabVect,TabSource);
	       	
	       	if (MadStop <= 0)
	       	{
	       		RelError = WT.L2err;
	       		if (Verbose == True) cout << "Running GMCA ... Global Relative Error : " << RelError << " / " << RelErrorMax << " dB" << endl;
	       		if (RelError > RelErrorMax) ExitCriterion = True;
	       	}
	       	
	        if (MadStop > 0)
	       	{
	       		RelError = WT.L2err;
	       		if (Verbose == True) cout << "Running GMCA ... Global Relative Error : " << RelError << " / " << RelErrorMax << " dB" << endl;
	       		DiffRelError = abs(RelError - OldRelError)/RelError;
	       		
	       		if (RelError > RelErrorMax) ExitCriterion = True;
	       		if (MadStop > DiffRelError) ExitCriterion = True;
	       	}
	       	
	       	
	       	WT.NbrColumnInit = NbrSources; // Keep the last mixmat
	       	WT.MatColumnInit = WT.MixingMat;
	       	
	       		       	
       	}
       	
       	if (Verbose == True) cout << "Stopping GMCA" << endl;
       }
       
       
       
       // Reconstruction :
       if (Verbose == True) cout << "Reconstruction ... "<< endl;
       fltarray EstSources;
       WT.recons_sources(Dat,EstSources);
       // WT.Sort_Sources(EstSources);

       if (Verbose == True) cout << "Write results ... "<< endl;
       // fits_write_fltarr ("xx_EstSources.fits", EstSources);
       fits_write_fltarr ("xx_EstMixmat.fits", WT.MixingMat);
       fits_write_fltarr ("xx_InvMixingMat.fits", WT.InvMixingMat);

    // Header.origin = Cmd;	 
    fits_write_fltarr(Name_Out, EstSources);
 
    exit(0);
}
