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
**    Date:  05/12/1999
**    
**    File:  im3d_deconv.cc
**
*******************************************************************************
**
**    DESCRIPTION: Deconvolution of a set of images using a set of PSF   
**    -----------  We search O such that:
**                    P1 * O = I1
**                     ...
**                    Pn * O = In
**
**                 Shift between images I is taken into account.
**                 The multiresolution support is used for the regularization
**
******************************************************************************/

#include "IM_Obj.h"
#include "IM_IO.h"
#include "IM3D_IO.h"
#include "CoaddCorrel.h"
#include "IM_Deconv.h"
#include "MR_Obj.h"
#include "MR_Filter.h"

char Name_Imag_In[256];
char Name_Imag_Out[256];
char Name_Offset[256];
char Name_PSF[256];

extern int  OptInd;
extern char *OptArg;

extern int  GetOpt(int argc, char **argv, char *opts);

#define OUTMAN stdout
 
float Zoom = 1.;
int MaxDist=4;
Bool OptOffset = False;
Bool OptNoOffset = False;

int MaxIter = 50;
Bool Pos = True;
Bool Verbose = False;
int PosX = -1;
int PosY = -1;
int Surface = 10;

type_interp_corr TypeInterp = ICF_TANH;
 
// int   Nbr_Plan=DEFAULT_NBR_SCALE;  /* number of scales */
int   Nbr_Plan=DEFAULT_NBR_SCALE;
float N_Sigma=DEFAULT_N_SIGMA;   /* number of sigma (for the noise) */
float Noise_Ima=0.;              /* noise standard deviation */
// type_noise Stat_Noise = DEFAULT_STAT_NOISE;   /* type of noise */
type_noise Stat_Noise = NOISE_UNI_UNDEFINED;
type_transform Transform = DEFAULT_MR_TRANS; /* type of transform */
Bool SupIsol=False;             /* suppress isolated pixel in the support */
Bool SupDil=False;              /* dilate the support */
Bool UseNSigma =False;  
Bool KeepPositivSup=False; 
Bool NoiseRegul=False;

#define NBR_DEC3D 4
enum T_DEC3D {DEC3D_CITTERT, DEC3D_LAND, DEC3D_LUCY, DEC3D_MAP};
T_DEC3D TypeDeconv = DEC3D_MAP;

inline char * StringDec3D (T_DEC3D type)
{
    switch (type)
    {
        case DEC3D_CITTERT:
             return ("Van-Citter iteration");break;
        case DEC3D_LAND:
	     return ("Landweber iteration");break;
        case DEC3D_LUCY:
	     return ("Lucy iteration");break;
        case DEC3D_MAP:
	     return ("MAP iteration");break;
   }
   return ("Error: bad type of deconvolution");
}
inline void dec3d_usage(T_DEC3D typeDec)
{
    fprintf(OUTMAN, "         [-d type_of_deconvolution]\n");
    fprintf(OUTMAN, "              %d: %s \n",1,StringDec3D(DEC3D_CITTERT));
    fprintf(OUTMAN, "              %d: %s \n",2,StringDec3D(DEC3D_LAND));
    fprintf(OUTMAN, "              %d: %s \n",3,StringDec3D(DEC3D_LUCY));
    fprintf(OUTMAN, "              %d: %s \n",4,StringDec3D(DEC3D_MAP));
    fprintf(OUTMAN, "              default is %s.\n", StringDec3D(typeDec));
}

/****************************************************************************/

static void usage(char *argv[])
{
   int i;
   
    fprintf(stdout, "Usage: %s options in_data_cube in_psf_image out_image\n\n", argv[0]);
     
    fprintf(stdout, "   where options =  \n");
        
    dec3d_usage(TypeDeconv);
    manline();
	 
    fprintf(OUTMAN, "         [-f type_of_interpolation]\n");
    for (i = 0; i < NBR_CORR_INTERP; i++)
    fprintf(OUTMAN, "              %d: %s \n",i+1,
                                       StringCorrInterp ((type_interp_corr)i));
    fprintf(OUTMAN, "              Default is %s\n", StringCorrInterp (TypeInterp));
    manline();

    fprintf(stdout, "         [-r ZoomFactor]\n");
    fprintf(stdout, "              Rebin the reconstructed image.\n\n");

    fprintf(stdout, "         [-N]\n");
    fprintf(stdout, "              No offset.\n\n");
    
    fprintf(stdout, "         [-o InputOffsetFileName]\n");
    fprintf(stdout, "              Offset Fits array (2,Nz).\n");
    fprintf(stdout, "              Array(0,z) = x offset from frame z to first frame.\n");
    fprintf(stdout, "              Array(1,z) = y offset from frame z to first frame.\n\n");
    
    fprintf(stdout, "         [-i MaxIter]\n");
    fprintf(stdout, "              Number of iterations. Default is %d\n", MaxIter);
    manline();
      
    fprintf(stdout, "         [-p]\n");
    fprintf(stdout, "             Suppress the positivity constraint.\n");
    manline();

    fprintf(stdout, "         [-W]\n");
    fprintf(stdout, "             Regularization by the Wavelet Transform.\n");
    manline();

    transform_usage(Transform);
    manline();
    nbr_scale_usage(Nbr_Plan);
    manline();
    
//    gauss_usage();
//    manline();
    
    nsigma_usage(N_Sigma);
    manline();
    
    vm_usage();
    manline();    
    verbose_usage();
    manline();
    exit(-1);
}

/*********************************************************************/

  
/* GET COMMAND LINE ARGUMENTS */

static void transinit(int argc, char *argv[])
{
    int c;
    Bool Optwave = False;
 #ifdef LARGE_BUFF
    int VMSSize=-1;
    Bool OptZ = False;
    char VMSName[1024] = "";
#endif    
    /* get options */
    while ((c = GetOpt(argc,argv,"d:NWn:t:s:f:r:o:i:pvzZ:")) != -1) 
    {
	switch (c) 
        {
           case 'N': OptNoOffset = (OptNoOffset == True) ? False: True;break;
           case 'W': NoiseRegul = (NoiseRegul == True) ? False: True;break;
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
 		Optwave = True;
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
		Optwave = True;
 		break;
           case 's':
                /* -s <nsigma> */
                if (sscanf(OptArg,"%f",&N_Sigma) != 1) 
                {
                    fprintf(OUTMAN, "Error: bad N_Sigma: %s\n", OptArg);
                    exit(-1);
                }
                Optwave =True;
                if (N_Sigma <= 0.)  N_Sigma = DEFAULT_N_SIGMA;
                break;
	   case 'd':
		/* -d <type> type of transform */
		if (sscanf(OptArg,"%d",&c ) != 1) 
                {
		    fprintf(OUTMAN, "Error: bad type of deconvolution: %s\n", OptArg);
	            exit(-1);
                    
		}
                if ((c > 0) && (c <= NBR_DEC3D)) TypeDeconv = (T_DEC3D) (c-1);
                else  
                {
		    fprintf(OUTMAN, "Error: bad type of deconvolution: %s\n", OptArg);
	            exit(-1);
 		}
 		break;		
	 case 'f':
	        if (sscanf(OptArg,"%d",&c ) != 1) 
                {
                    fprintf(OUTMAN, "Error: bad interpolation method: %s\n", OptArg);
                    exit(-1);
                }
		if ((c > 0) && (c <=  NBR_CORR_INTERP))  
		       TypeInterp = (type_interp_corr) (c-1);
                else  
                {
                    fprintf(OUTMAN, "Error: bad interpolation method: %s\n", OptArg);
                    exit(-1);
                }
		break;
 	 case 'v': Verbose = True; break;
        case 'i':
	       if (sscanf(OptArg,"%d",&MaxIter) != 1) 
                {
		    fprintf(stdout, "bad number of iterations: %s\n", OptArg);
		    usage(argv);
		}
		break; 	
         case 'r':
	       if (sscanf(OptArg,"%f",&Zoom) != 1) 
                {
		    fprintf(stdout, "bad zoom factor: %s\n", OptArg);
		    usage(argv);
		}
		break; 
         case 'o':
	       if (sscanf(OptArg,"%s", Name_Offset) != 1) 
                {
		    fprintf(stdout, "bad offset file name: %s\n", OptArg);
		    usage(argv);
		}
		OptOffset = True;
		break;  
 	case 'p': Pos = False; break;
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
	
	if ((Optwave == True) && (NoiseRegul == False))
	{
	   cout << "Error: t,s and n options are not valid when -W is not set ... " << endl;
	   exit(-1);
	}
       /* get optional input file names from trailing 
          parameters and open files */
	if (OptInd < argc) strcpy(Name_Imag_In, argv[OptInd++]);
         else usage(argv);

	if (OptInd < argc) strcpy(Name_PSF, argv[OptInd++]);
         else usage(argv);
	 
	if (OptInd < argc) strcpy(Name_Imag_Out, argv[OptInd++]);
        else usage(argv);

	/* make sure there are not too many parameters */
	if (OptInd < argc)
        {
		fprintf(stdout, "Too many parameters: %s ...\n", argv[OptInd]);
		usage(argv);
	}
	
#ifdef LARGE_BUFF
    if (OptZ == True) vms_init(VMSSize, VMSName, Verbose);
#endif
}

/*********************************************************************/

void cube_to_ima(fltarray &Data, Ifloat &Ima)
{
   int Nx, Ny, Nz;
   Nx = Data.nx();
   Ny = Data.ny();
   Nz = Data.nz();
   int x,y,z;

   for (x=0; x < Nx; x++)
   for (y=0; y < Ny; y++)
   {	
       Ima(y,x) = 0.;
       for (z=0; z < Nz; z++) Ima(y,x) += Data(x,y,z);
       Ima(y,x) /= (float) Nz;
   }  
}

/*********************************************************************/

void cube_to_ima(fltarray &Dat, Ifloat &Ima, xcorr_def *Xcor, float Zoom,
                 char *interp_method)
{
   int ValRet, i, nz = Dat.nz();
   int Nx1 = (int)(Zoom * Dat.nx());
   int Ny1 = (int)(Zoom * Dat.ny());
   fltarray DatOut(Nx1, Ny1, Dat.nz());   

   float **TabFrame;
   float **TabFrameOut;
   TabFrame = new float * [nz];
   TabFrameOut = new float * [nz];   
   for (i = 0; i < nz; i++) 
        TabFrame[i] = Dat.buffer() + i*Dat.nx()*Dat.ny();   
   for (i = 0; i < nz; i++) 
        TabFrameOut[i] = DatOut.buffer() + i*DatOut.nx()*DatOut.ny();

   ValRet = cube_resample(TabFrame, TabFrameOut, Xcor, Dat.nx(), Dat.ny(), Dat.nz(),
                           DatOut.nx(), DatOut.ny(), interp_method); 
   cube_to_ima(DatOut, Ima);
}

/****************************************************************************/

void ima_to_cube(Ifloat &Ima, fltarray &Dat, xcorr_def *Xcor, float Zoom,
                 char *interp_method)
// Project an image to a cube taking into account the zoom factor and
// shift.		 
{
   int i, nz = Dat.nz();
   int Nx1 = Dat.nx();
   int Ny1 = Dat.ny();
   double dx, dy;
   
   float *Frame;
   for (i = 0; i < nz; i++) 
   {
      Frame = Dat.buffer() + i*Nx1*Ny1;     
      dx = (Xcor->xc_x)[i] * Zoom;
      dy = (Xcor->xc_y)[i] * Zoom;
      warp_image(Ima.buffer(), Ima.nc(), Ima.nl(),
                 Frame, Nx1, Ny1, -dx, -dy, (double) (1./Zoom), interp_method);
   }
}


/****************************************************************************/

void calcul_residu(Ifloat & Ima, fltarray & Data, fltarray & Resi,  
                   xcorr_def *Xcor, float Zoom, char *interp_method)
{
   int i;
   int nz = Data.nz();
      
   ima_to_cube(Ima, Resi, Xcor, Zoom, interp_method);
   for (i = 0; i < nz; i++) 
   {
      for (int x=0; x < Data.nx(); x++)
      for (int y=0; y < Data.ny(); y++) Resi(x,y,i) = Data(x,y,i) - Resi(x,y,i);
   }
}

/****************************************************************************/

void calcul_gradient(Ifloat &Obj, Ifloat & Gradient, fltarray & Resi, 
                  Icomplex_f *TabPsf_cf, xcorr_def *Xcor, 
		  float Zoom, char *interp_method)
{
   int i;
   int nz = Resi.nz();
   Ifloat Image;
   
   for (i = 0; i < nz; i++) 
   {
      Image.alloc(Resi.buffer()+Resi.nx()*Resi.ny()*i, Resi.ny(), Resi.nx());
      psf_convol_conj(Image, TabPsf_cf[i]);
   }
   cube_to_ima(Resi, Gradient, Xcor, Zoom, interp_method);
}

/****************************************************************************/

void dec_cube(Ifloat &Obj, fltarray & Data, fltarray & Resi, 
 //             Icomplex_f *TabPsf_cf,
              Icomplex_f & Psf_cf,
              int MaxIter, xcorr_def *Xcor, float Zoom, 
	      Bool Positiv, char *interp_method)
{
   int i,k,l,NbrBand;
   Ifloat Gradient(Obj.nl(), Obj.nc(), "gradient");
   Ifloat ImaConv(Data.ny(), Data.nx(), "ima");

   MRNoiseModel ModelData;
   MultiResol MR_Data;
   
   // noise model class initialization
   if (NoiseRegul == True)
   {
      ModelData.alloc (Stat_Noise, Obj.nl(),Obj.nc(),Nbr_Plan, Transform);
      NbrBand = ModelData.nbr_band();
      if (Noise_Ima > FLOAT_EPSILON) ModelData.SigmaNoise = Noise_Ima;
      if (UseNSigma  == True)
          for (i=0; i < NbrBand; i++) ModelData.NSigma[i]=N_Sigma;
      if (SupIsol == True) ModelData.SupIsol = True;
      if (SupDil == True) ModelData.DilateSupport = True;
      ModelData.OnlyPositivDetect = KeepPositivSup;
  
      MR_Data.alloc(Obj.nl(), Obj.nc(), Nbr_Plan, Transform, "MRNoiseModel");
      ModelData.model(Obj, MR_Data);
      ModelData.threshold(MR_Data);
      MR_Data.recons(Obj);
      if (Positiv == True) threshold(Obj);
  }
//    INFO(MR_Data.band(0), "band 0");
//    ModelData.threshold(MR_Data);
//    INFO(MR_Data.band(0), "band 0");
//    ModelData.write_support_mr("sup.mr");

   for (i = 0; i < MaxIter; i++)
   {
       ImaConv = Obj;
       psf_convol(ImaConv, Psf_cf);
       calcul_residu(ImaConv, Data, Resi, Xcor, 
                     Zoom, interp_method);
       // Resi.fits_write("resi.fits");
       
       // calcul_gradient(Obj, Gradient, Resi, TabPsf_cf,
       //               Xcor, Zoom, interp_method);
       cube_to_ima(Resi, Gradient, Xcor, Zoom, interp_method);

       if (NoiseRegul == True)
       {
          MR_Data.transform( Gradient);
          ModelData.threshold(MR_Data);
          MR_Data.recons( Gradient);
       }
       
       switch (TypeDeconv)
       {
          case DEC3D_LAND:
	          psf_convol_conj(Gradient, Psf_cf);
          case DEC3D_CITTERT:
               for (k=0; k < ImaConv.nl(); k++)
               for (l=0; l < ImaConv.nc(); l++) Obj(k,l) += Gradient(k,l); 
               break;
	  case DEC3D_LUCY:
               for (k=0; k < ImaConv.nl(); k++)
               for (l=0; l < ImaConv.nc(); l++)
               {
                  if (ImaConv(k,l) > FLOAT_EPSILON)
                  Gradient(k,l) = ( ImaConv(k,l) + Gradient(k,l)) /  ImaConv(k,l);
	          else Gradient(k,l) = FLOAT_EPSILON;
               }
               psf_convol_conj(Gradient, Psf_cf);
               for (k=0; k < ImaConv.nl(); k++)
               for (l=0; l < ImaConv.nc(); l++) Obj(k,l) *= Gradient(k,l);	       break;
          case DEC3D_MAP: 
               for (k=0; k < ImaConv.nl(); k++)
               for (l=0; l < ImaConv.nc(); l++)
               {
                  if (ImaConv(k,l) > FLOAT_EPSILON)
                  Gradient(k,l) = Gradient(k,l) /  ImaConv(k,l);
	          else Gradient(k,l) = 0.;
               }
               psf_convol_conj(Gradient, Psf_cf);
               for (k=0; k < ImaConv.nl(); k++)
               for (l=0; l < ImaConv.nc(); l++)  Obj(k,l) *= exp(Gradient(k,l));
               break;
       }
//        // Resi.info("RESIDUAL"); 
//        // Gradient.info("GRADIENT");
       if (Verbose == True)
       {
          cout << i+1 << ": Sigma(Resi) = " << Resi.sigma();
          cout <<  "  Sigma grad = " << sigma(Gradient) << endl;
       }
       // Positivity constraint
       if (Positiv == True) threshold(Obj);
   }
}

/****************************************************************************/

int main(int argc, char *argv[])
{
    fltarray Dat;
    fltarray Mallat;
    int i,ValRet;
    fltarray TabPsf;
    Ifloat Psf;
    
    /* Get command line arguments, open input file(s) if necessary */
    lm_check(LIC_MR3);
    transinit(argc, argv);
    
    if (Verbose == True)
    {
      cout << endl << endl << "PARAMETERS: " << endl << endl;
      cout << "Name File in = " << Name_Imag_In << endl ;
      cout << "Name File out = " << Name_Imag_Out << endl ;
      // cout << "nelem = " << Dat.n_elem() << endl;
      // cout << "Sigma = " << Dat.sigma();
      // cout << " Min   = " << Dat.min() << " Max   = "  << Dat.max() << endl;
      cout << "Interpolation method = " << StringCorrInterp(TypeInterp) << endl;
      if (OptNoOffset == True) cout << "No offset calculation " << endl;
      cout << "Deconvolution Method = " << StringDec3D(TypeDeconv) << endl;
      if (NoiseRegul  == True) 
      {
         cout << "Regularization using the wavelet transform " << endl;
         cout << "Transform = " << StringTransform(Transform) << endl;
         cout << "Number of scales = " << Nbr_Plan << endl;
         cout << "N_Sigma = " << N_Sigma << endl;
      }
      if (Zoom != 1) cout << "Zoom = " << Zoom << endl;
   }
   
      io_3d_read_data(Name_Imag_In, Dat);
   // io_3d_read_data(Name_PSF, TabPsf);
   io_read_ima_float(Name_PSF, Psf);
   if (Verbose == True)
       cout << "nx = " << Dat.nx() << " ny = " << Dat.ny() << " nz = " << Dat.nz() << endl;

//    if (TabPsf.nz() != Dat.nz())
//    {
//        cerr << "Error: bad number of frames in the PSF file ..." << endl;
//        exit(-1);
//    }   
   int nz = Dat.nz();
   xcorr_def Xcor;
   Xcor.xc_np = nz;
   Xcor.xc_x = new double [nz];
   Xcor.xc_y = new double [nz];
   Xcor.xc_level = new double [nz];
   for (i=0; i < nz; i++) Xcor.xc_x[i] = Xcor.xc_y[i] = Xcor.xc_level[i] = 0.;
   Xcor.xc_dx = MaxDist;
   Xcor.xc_dy = MaxDist;
   if ((Surface < 1) || (Surface > Dat.nx()/2))
   {
      Xcor.xc_hx = MAX(0, Dat.nx()/2 -1 - MaxDist);
      Xcor.xc_hy = MAX(0, Dat.ny()/2 -1 - MaxDist);
   }
   else
   {
      Xcor.xc_hx = Surface;
      Xcor.xc_hy = Surface;
   }
   if (PosX < 0) PosX =  Dat.nx() / 2;
   if (PosY < 0) PosY =  Dat.ny() / 2;
       
   //Xcor.xc_hx = 10;
   //Xcor.xc_hy = 10;
   Xcor.xc_method = 0;
   const char *interp_method= StringCorrInterp(TypeInterp);
   
   float **TabFrame = new float * [nz];
   float * TabMean = new float [nz];
   for (i = 0; i < nz; i++) TabFrame[i] = Dat.buffer() + i*Dat.nx()*Dat.ny();
   
   // If an input offset file is available, we used it
   if (OptOffset == True)
   {
      fltarray TabOffset;
      fits_read_fltarr(Name_Offset, TabOffset);
      for (i = 0; i < nz; i++) 
      {
         Xcor.xc_x[i] = (double) TabOffset(0,i);
	 Xcor.xc_y[i] = (double) TabOffset(1,i);
      }
      Xcor.xc_init = 1;
   }
   else if (OptNoOffset == False) // otherwise we calculte the offset frame by correlation
   {  
      Xcor.xc_init = 0;
      for (i = 0; i < nz; i++)  Xcor.xc_x[i] = Xcor.xc_y[i] = 0.;
      float *Pattern = TabFrame[0];
      
      ValRet = cube_get_offset(TabFrame, Pattern, 
                               Dat.nx(), Dat.ny(), Dat.nz(), 
			       PosX, PosY, &Xcor);
      if (Verbose == True)
        for (i = 0; i < nz; i++)
           cout << "Offset Frame: " << i+1 << ", dx = " << Xcor.xc_x[i] << " dy = " << Xcor.xc_y[i] << endl;
   }
   
   int Nx1 = (int)(Zoom * Dat.nx());
   int Ny1 = (int)(Zoom * Dat.ny());
   
//    if ((Ny1 != Psf.nl()) || (Nx1 != Psf.nc()))
//    {
//       cout << "Error: bad PSF size: " << Psf.nl() << "x" << Psf.nc() << endl;
//       cout << "       Expected size: " <<  Ny1 << "x" << Nx1 << endl;
//       cout << "       Zoom = " << Zoom << endl;
//       exit(-1);
//    }
   Ifloat Ima(Ny1, Nx1, "ResIma");
   cube_to_ima(Dat, Ima, &Xcor, Zoom, (char *) interp_method);

   // The multiresolution support must be computed from Ima

   // Deconvolution Part
   fltarray Resi(Dat.nx(), Dat.ny(), Dat.nz());
   // Icomplex_f *TabPsf_cf;
   // Ifloat Psf;
   // TabPsf_cf = new Icomplex_f [nz];
   // for (i = 0; i < nz; i++) 
   // {
   //     Psf.alloc(TabPsf.buffer()+TabPsf.nx()*TabPsf.ny()*i, 
   //              TabPsf.ny(), TabPsf.nx());
   //    psf_get(Psf, TabPsf_cf[i], Ny1, Nx1, True);
   //}
   Icomplex_f Psf_cf;
   psf_get(Psf, Psf_cf, Ny1, Nx1, True);
    
   if (Verbose == True) cout << "run deconv ... " <<  Ny1 << "x" << Nx1 << endl;
   // dec_cube(Ima, Dat, Resi, TabPsf_cf, MaxIter, &Xcor,  
   //             Zoom, Pos, interp_method);
   dec_cube(Ima, Dat, Resi, Psf_cf, MaxIter, &Xcor,
            Zoom, Pos, (char *) interp_method);
	    
   io_write_ima_float(Name_Imag_Out, Ima);
   delete [] TabMean;
   delete [] TabFrame;
   exit(0);
}


