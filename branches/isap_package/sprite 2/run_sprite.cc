#include <stdio.h>
#include <stdlib.h>
#include "sr_util.h"
#include "sprite.h"
#include "IM_Obj.h"
#include "IM_IO.h"
#include <iomanip>


extern int  OptInd;
extern char *OptArg;

#define OUTMAN stdout

char Name_Cube_In[256]; /* input file image */
char Name_Imag_Out[256]; /* output file name */
char Name_Resi_Out[256]; /* Residual file name */
char outputs_Dir[256]; /* Outputs directory file name */

int Nbr_Plan=DEFAULT_NBR_SCALE;  /* number of scales */
float N_Sigma=DEFAULT_N_SIGMA;   /* number of sigma (for the noise) */
Bool WriteResi;
Bool WriteParam;
bool Verbose=False;
float ResolUp = 2.;
int DEFAULT_SRITE_FILTER=50;
int Max_Iter = DEFAULT_SRITE_FILTER; /* Maximum number of iteration */
int Max_Iter_Noise_Sim = 200;
int Noise_est_meth = 0; /* Noise estimation method for thresholding */
int Mod_denoise_meth = MR_SUPP_DENOISING; /* Model denoising method */
type_transform transform = TO_UNDECIMATED_MALLAT; /* Dictionary for inverse problem constraint */
type_transform transform_mod = TO_UNDECIMATED_MALLAT; /* Dictionary for first guess denoising */
Bool Thresh_type=True;
int nb_rw=2;
Bool eq_noise = True;
Bool eq_flux = True;
int cv_win = 50;
float cv_tol = 0.1;
char* shift_filename = NULL;
/****************************************************************************/

static void usage(char *argv[])
{
  fprintf(OUTMAN, "Usage: %s options in_Cube out_Image outputs_Dir\n\n", argv[0]);
  fprintf(OUTMAN, "   where options =  \n");
  all_transform_usage(transform);
  manline();
  fprintf(OUTMAN, "         [-T type_of_multiresolution_transform for model denoising.]\n");
  fprintf(OUTMAN, "             Same optional transforms as -t optional.\n");
  fprintf(OUTMAN, "             Default is %f.\n", ResolUp);
  thresh_type_usage();
  nbr_scale_usage(Nbr_Plan);
     
  nsigma_usage(N_Sigma);
  
  max_iter_usage(Max_Iter);

  fprintf(OUTMAN, "         [-r ResolUp]\n");
  fprintf(OUTMAN, "             Default is %f.\n", ResolUp);
  
  fprintf(OUTMAN, "         [-R nb_rw]\n");
  fprintf(OUTMAN, "             Number of reweighted passes.\n");
  fprintf(OUTMAN, "             Default is %d.\n", nb_rw);
    
  fprintf(OUTMAN, "         [-w ResiduFileName]\n");
  fprintf(OUTMAN, "             By default, the residu is not written.\n");
  
  fprintf(OUTMAN, "         [-S ShiftFileName]\n");
  fprintf(OUTMAN, "             If not provided, the shifts are estimated.\n");
  
  noise_est_usage();
  
  fprintf(OUTMAN, "         [-I number_of_realisations]\n");
  fprintf(OUTMAN, "             Number of realisations in the case noise simulation is chosen.\n");
  fprintf(OUTMAN, "             Default is %d.\n", Max_Iter_Noise_Sim);
  
  mod_denoise_usage();


  fprintf(OUTMAN, "         [-C cv_win]\n");
  fprintf(OUTMAN, "             Number of iterates for convergence assessment.\n");
  fprintf(OUTMAN, "             Default is %d.\n", cv_win);

  fprintf(OUTMAN, "         [-a cv_tol]\n");
  fprintf(OUTMAN, "             Tolerance parameter for convergence assessment.\n");
  fprintf(OUTMAN, "             Default is %f.\n", cv_tol);
  
  fprintf(OUTMAN, "         [-W]\n");
  fprintf(OUTMAN, "             If set, different diagnostic parameters are saved on the disk (centroids, flux estimates...).\n");

  fprintf(OUTMAN, "         [-F]\n");
  fprintf(OUTMAN, "             If set, the photometric flux are estimated, otherwise they are assumed to be equal.\n");

  fprintf(OUTMAN, "         [-N]\n");
  fprintf(OUTMAN, "             If not set, the noise level is assumed to be the same in the low resolution images.\n");
    
  vm_usage();
  verbose_usage();
    
  manline();
  manline();
  exit(-1);
} 

/*********************************************************************/

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
  while ((c = GetOpt(argc,argv, (char *) "t:T:c:n:s:i:r:R:w:e:I:m:C:a:S:WFNvzZ:")) != -1) 
    { 
      switch (c) 
        {
	case 't':
	  /* -t <type> type of transform */
	  if (sscanf(OptArg,"%d",&c ) != 1) 
	    {
	      fprintf(OUTMAN, "bad type of multiresolution transform: %s\n", OptArg);
	      exit(-1);
                    
	    }
	  if ((c > 0) && (c <= NBR_TOT_TRANSFORM+1)) 
	    transform = (type_transform) (c-1);
	  else  
	    {
	      fprintf(OUTMAN, "bad type of transform: %s\n", OptArg);
	      exit(-1);
	    }                
	  break;
	case 'T':
	  /* -T <type> type of transform */
	  if (sscanf(OptArg,"%d",&c ) != 1) 
	    {
	      fprintf(OUTMAN, "bad type of multiresolution transform for model denoising: %s\n", OptArg);
	      exit(-1);
                    
	    }
	  if ((c > 0) && (c <= NBR_TOT_TRANSFORM+1)) 
	    transform_mod = (type_transform) (c-1);
	  else  
	    {
	      fprintf(OUTMAN, "bad type of transform: %s\n", OptArg);
	      exit(-1);
	    }                
	  break;
	case 'c':
	  /* -c <type> type of thresholding */
	  if (sscanf(OptArg,"%d",&c ) != 1) 
	    {
	      fprintf(OUTMAN, "bad type of thresholding: %s\n", OptArg);
	      exit(-1);
                    
	    }
	  if ((c >= 0) && (c <2)) 
	    Thresh_type = (Bool)c;
	  else  
	    {
	      fprintf(OUTMAN, "bad type of transform: %s\n", OptArg);
	      exit(-1);
	    }                
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
	  if (sscanf(OptArg,"%f",&N_Sigma) != 1) 
	    {
	      fprintf(OUTMAN, "Error: bad N_Sigma: %s\n", OptArg);
	      exit(-1);
	    }
	  if (N_Sigma <= 0.)  N_Sigma = DEFAULT_N_SIGMA;
	  break;
	case 'i':
	  /* -i < Number of iterations> */
	  if (sscanf(OptArg,"%d",&Max_Iter) != 1) 
	    {
	      fprintf(OUTMAN, "Error: bad Max_Iter: %s\n", OptArg);
	      exit(-1);
	    }
	  if (Max_Iter <= 0)   Max_Iter = DEFAULT_SRITE_FILTER;
	  break;
	case 'r':
	  if (sscanf(OptArg,"%f",&ResolUp) != 1) 
	    {
	      fprintf(OUTMAN, "Error: bad ResolUp: %s\n", OptArg);
	      exit(-1);
	    }
	  if (ResolUp <= 0.)  ResolUp = 2;
	  break;    
	case 'R':
	  if (sscanf(OptArg,"%d",&nb_rw) != 1) 
	    {
	      fprintf(OUTMAN, "Error: nb_rw: %s\n", OptArg);
	      exit(-1);
	    }
	  break;
	case 'C':
	  if (sscanf(OptArg,"%d",&cv_win) != 1) 
	    {
	      fprintf(OUTMAN, "Error: cv_win: %s\n", OptArg);
	      exit(-1);
	    }
	  break;
	case 'a':
	  
	  if (sscanf(OptArg,"%f",&cv_tol) != 1) 
	    {
	      fprintf(OUTMAN, "Error: cv_tol: %s\n", OptArg);
	      exit(-1);
	    }
	  break;
	     
	  
	case 'w':
	  /* -w < residual file name> */
	  if (sscanf(OptArg,"%s", Name_Resi_Out) != 1) 
	    {
	      fprintf(OUTMAN, "Error: bad file name: %s\n", OptArg);
	      exit(-1);
	    }
	  WriteResi = True;
	  break;
	case 'S':
	  /* -S < shifts file name> */
	  shift_filename = new char[1024];
	  if (sscanf(OptArg,"%s", shift_filename) != 1) 
	    {
	      fprintf(OUTMAN, "Error: bad file name: %s\n", OptArg);
	      exit(-1);
	    }
	  WriteResi = True;
	  break;
	case 'e':
	    
	  if (sscanf(OptArg,"%d",&Noise_est_meth) != 1) 
	    {
	      fprintf(OUTMAN, "Error: bad Noise_est_meth: %s\n", OptArg);
	      exit(-1);
	      
	    }
	  break;
	case 'I':
	  /* -I < Number of iterations in noise simulation> */
	  if (sscanf(OptArg,"%d",&Max_Iter_Noise_Sim) != 1) 
	    {
	      fprintf(OUTMAN, "Error: bad Max_Iter_Noise_Sim: %s\n", OptArg);
	      exit(-1);
	    }
	  break;
	case 'm':
	  /* -m < Model denoising method> */
	  if (sscanf(OptArg,"%d",&Mod_denoise_meth) != 1) 
	    {
	      fprintf(OUTMAN, "Error: bad Mod_denoise_meth: %s\n", OptArg);
	      exit(-1);
	    }
	  break;
	case 'W':WriteParam = True;break;
	case 'F':eq_flux = False;break;
	case 'N':eq_noise = False;break;
	case 'v': Verbose = True;break;
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
  if (OptInd < argc) strcpy(Name_Cube_In, argv[OptInd++]);
  else usage(argv);
    
  if (OptInd < argc) strcpy(Name_Imag_Out, argv[OptInd++]);
  else usage(argv);

  if (OptInd < argc) strcpy(outputs_Dir, argv[OptInd++]);
  else usage(argv);
    
  	     
       
#ifdef LARGE_BUFF
  if (OptZ == True) vms_init(VMSSize, VMSName, Verbose);
#endif  
}


/****************************************************************************/


int main(int argc, char *argv[]) 
{
  int s,i,j,k;
  Ifloat DataB;
  fitsstruct Header;
  char Cmd[256];
    
  Cmd[0] = '\0';
  for (k =0; k < argc; k++) sprintf(Cmd, "%s %s", Cmd, argv[k]);
    
  /* Get command line arguments, open input file(s) if necessary */
  filtinit(argc, argv);
  if (Verbose == True)
    {
      cout << endl << endl << "PARAMETERS: " << endl << endl;
      cout << "File Cube in = " << Name_Cube_In << endl;
      cout << "File Ima Out = " << Name_Imag_Out << endl;
      cout << "Outputs directory = "<< outputs_Dir << endl;
      cout << "Type of transform for sparse constraint = "<< transform <<endl;
      cout << "Type of transform for first guess denoising = "<< transform_mod <<endl;
      cout << "Thresholding method = "<< Thresh_type <<endl;
      cout << "Number of scales = " << Nbr_Plan << endl;
      cout << "Max_Iter = " << Max_Iter << endl;
      cout << "Number of reweighted passes = " << nb_rw << endl;
      cout << "Zooming factor = " << ResolUp << endl;
      cout << "N_Sigma = "<< N_Sigma << endl; 
      cout << "Noise_est_meth: " << Noise_est_meth << endl;
      cout << "Number of iterates used for convergence assessment = " <<cv_win << endl;
      cout << "Associated tolerance parameter: " <<cv_tol << endl;
      if(Noise_est_meth ==1)
	cout << "Max_Iter_Noise_Sim = " << Max_Iter_Noise_Sim << endl;
      if (WriteResi == True)
	cout << "Residual Cube file name : " << Name_Resi_Out << endl;
      if (shift_filename!=NULL)
	cout << "Shift filename : "<< shift_filename << endl;
      if (eq_flux == False)
	cout << "Photometric flux estimation enabled." <<endl;
      if (eq_noise == False)
	cout << "The noise level is not assumed to be the same in the low resolution images." <<endl;
      if(WriteParam==True)
	cout << "Diagnostic files will be saved in the outputs directory"<< endl;
    }
  
  string output_file = Name_Imag_Out;
  string rep_path = outputs_Dir;
  string resi_file = Name_Resi_Out;
  int D0  =12/ResolUp; // Scale parameter for weighted centroid estimation
  double sig = 89.108911; // Weighting function standard deviation   
  double r = 3; // Aperture radius for flux estimation
  fltarray Dat;
  io_3d_read_data(Name_Cube_In, Dat, &Header);
  
  int Nz = Dat.nz(),Nx = Dat.nx(),Ny=Dat.ny();
  mr_opt sr_opt = mr_opt_init(Nx*ResolUp,transform,-1,-1,NORM_L1,DEF_UNDER_FILTER,F_MALLAT_7_9,I_CONT,False);
  int siz = Nx*ResolUp;
  if (Mod_denoise_meth==LR_DENOISING)
    siz=Nx;
  mr_opt mod_opt = mr_opt_init(siz,transform_mod,Nbr_Plan,Nbr_Plan,NORM_L1,DEF_UNDER_FILTER,F_MALLAT_7_9,I_CONT,False);
  Sprite sprite_test(ResolUp,&Dat,D0,mod_opt,sr_opt,Max_Iter,r,sig,Mod_denoise_meth,N_Sigma,rep_path,output_file,resi_file,Noise_est_meth,0,eq_flux,eq_noise,cv_win,cv_tol,Verbose,shift_filename);
  Bool Wr_en=True;
  sprite_test.do_it(Wr_en,WriteResi,WriteParam,Max_Iter_Noise_Sim,Thresh_type,nb_rw);
  if (shift_filename!=NULL)
    delete [] shift_filename;
  free(sr_opt.FAS);
  free(mod_opt.FAS);
  Dat.free();
  return 0;
}


