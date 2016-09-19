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
**    Date:  9/2/99
**    
**    File:  im_coadd.cc
**
*******************************************************************************
**
**    DESCRIPTION   
**    ----------- 
**                 
**
**
******************************************************************************/

#include "IM_Obj.h"
#include "IM_IO.h"
#include "IM3D_IO.h"
#include "CoaddCorrel.h"
#include "IM_Deconv.h"
#include "OptMedian.h"

char Name_Imag_In[512];
char Name_Imag_Out[512];
char Name_Offset[512];
char Name_Cube_Out[512];
char Name_RMS_Out[512];

extern int  OptInd;
extern char *OptArg;

extern int  GetOpt(int argc, char **argv, char *opts);

#define OUTMAN stdout

float Zoom = 1.;
int MaxDist=4;
Bool OptOffset = False;
Bool OptWOffset = False;
Bool OptNoOffset = False;
Bool Noise_STD = False;
Bool GetMedian= False;
Bool WriteCube = False;
Bool MeanSub = False;

Bool Verbose = False;
int PosX = -1;
int PosY = -1;
int Surface = 10;

type_interp_corr TypeInterp = ICF_TANH;
Ifloat SigmaIma;

/****************************************************************************/

static void usage(char *argv[])
{
   int i;
   
    fprintf(stdout, "Usage: %s options in_cube out_image\n\n", argv[0]);
     
    fprintf(stdout, "   where options =  \n");
    fprintf(OUTMAN, "         [-f type_of_interpolation]\n");
    for (i = 0; i < NBR_CORR_INTERP; i++)
    fprintf(OUTMAN, "              %d: %s \n",i+1,
                                       StringCorrInterp ((type_interp_corr)i));
    fprintf(OUTMAN, "              Default is %s\n", StringCorrInterp (TypeInterp));
    manline();

    fprintf(stdout, "         [-r ZoomFactor]\n");
    fprintf(stdout, "              Rebin the reconstructed image.\n\n");

    fprintf(stdout, "         [-d MaxDist]\n");
    fprintf(stdout, "              Maximum offset.\n");
    fprintf(stdout, "              Default is %d\n\n", MaxDist);

    fprintf(stdout, "         [-a Surface]\n");
    fprintf(stdout, "              Surface measurement parameter.\n");
    fprintf(stdout, "              Default is %d\n\n", Surface);

    fprintf(stdout, "         [-x XPos]\n");
    fprintf(stdout, "              X position. Default is image center.\n\n");

    fprintf(stdout, "         [-y YPos]\n");
    fprintf(stdout, "              Y position. Default is image center.\n\n");

    fprintf(stdout, "         [-m]\n");
    fprintf(stdout, "              Subtract to each frame its mean value.\n\n");
 
    fprintf(stdout, "         [-o InputOffsetFileName]\n");
    fprintf(stdout, "              Offset Fits array (2,Nz).\n");
    fprintf(stdout, "              Array(0,z) = x offset from frame z to first frame.\n");
    fprintf(stdout, "              Array(1,z) = y offset from frame z to first frame.\n\n");
    
    fprintf(stdout, "         [-w OutputOffsetFileName]\n");
    fprintf(stdout, "              Offset Fits array (2,Nz).\n\n");

    fprintf(stdout, "         [-W RegisterCubeFileName]\n");
    fprintf(stdout, "             Write on the disk the registered cube.\n\n");
    
    fprintf(stdout, "         [-M]\n");
    fprintf(stdout, "             Take the median instead of the mean when coadding the frame.\n\n");
    
    fprintf(stdout, "         [-R OutputRMSFileName]\n");
    fprintf(stdout, "              Write on the disk the Root Mean Square Map.\n\n");
    
    fprintf(stdout, "         [-N]\n");
    fprintf(stdout, "              No offset.\n\n");
    
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
#ifdef LARGE_BUFF
    int VMSSize=-1;
    Bool OptZ = False;
    char VMSName[1024] = "";
#endif    
    /* get options */
    while ((c = GetOpt(argc,argv,"R:MW:Na:x:y:mw:f:r:d:o:vzZ:")) != -1) 
    {
	switch (c) 
        {
	case 'R':
               if (sscanf(OptArg,"%s", Name_RMS_Out) != 1) 
                {
                    fprintf(stdout, "bad cube file name: %s\n", OptArg);
                    usage(argv);
                }
                Noise_STD = True;
                break; 
	case 'M': GetMedian = (GetMedian == True) ? False: True;break;
        case 'W':
               if (sscanf(OptArg,"%s", Name_Cube_Out) != 1) 
                {
                    fprintf(stdout, "bad cube file name: %s\n", OptArg);
                    usage(argv);
                }
                WriteCube = True;
                break; 
         case 'N': OptNoOffset = (OptNoOffset == True) ? False: True;break;
	 case 'm': MeanSub = True; break;
	 case 'a':
	        if (sscanf(OptArg,"%d",&Surface) != 1) 
                {
		    fprintf(stdout, "bad surface parameter: %s\n", OptArg);
		    usage(argv);
		}
		break; 
	 case 'x':
	        if (sscanf(OptArg,"%d",&PosX) != 1) 
                {
		    fprintf(stdout, "bad x position: %s\n", OptArg);
		    usage(argv);
		}
		break; 
	case 'y':
	        if (sscanf(OptArg,"%d",&PosY) != 1) 
                {
		    fprintf(stdout, "bad y position: %s\n", OptArg);
		    usage(argv);
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
         case 'r':
	       if (sscanf(OptArg,"%f",&Zoom) != 1) 
                {
		    fprintf(stdout, "bad zoom factor: %s\n", OptArg);
		    usage(argv);
		}
		break; 
        case 'd':
	       if (sscanf(OptArg,"%d",&MaxDist) != 1) 
                {
		    fprintf(stdout, "bad maximum distance: %s\n", OptArg);
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
	case 'w':
	       if (sscanf(OptArg,"%s", Name_Offset) != 1) 
                {
		    fprintf(stdout, "bad offset file name: %s\n", OptArg);
		    usage(argv);
		}
		OptWOffset = True;
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
        if ((OptOffset == True) && (OptWOffset == True))
	{
           fprintf(OUTMAN, "\nError: w and o options cannot are exclusive ...\n\n");
           exit(-1);
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
		fprintf(stdout, "Too many parameters: %s ...\n", argv[OptInd]);
		usage(argv);
	}
#ifdef LARGE_BUFF
    if (OptZ == True) vms_init(VMSSize, VMSName, Verbose);
#endif
}

/*********************************************************************/

void coadd_all_frames(fltarray &Data, Ifloat &Ima)
{
   int Nx, Ny, Nz;
   Nx = Data.nx();
   Ny = Data.ny();
   Nz = Data.nz();
   int x,y,z;
   double Val;
   
   if (GetMedian == False)
   {
     for (x=0; x < Nx; x++)
     for (y=0; y < Ny; y++)
     {	
         Val = 0.;
         for (z=0; z < Nz; z++) Val += Data(x,y,z);
         if (Nz > 0) Ima(y,x) = (float) (Val / (double) Nz);
	 else Ima(y,x) = 0;
     }
   }
   else
   {
      float *V = new float [Nz];
      for (x=0; x < Nx; x++)
      for (y=0; y < Ny; y++)
      {
         if (Nz > 0) 
	 {
            for (z=0; z < Nz; z++)  V[z] = Data(x,y,z);
	    Ima(y,x) = get_median(V, Nz);
	 }
	 else Ima(y,x) = 0.;
      }
      delete [] V;
   }
   
   if (Noise_STD == True)
   {
      double RMS;
      SigmaIma.resize(Ny,Nx);
      for (x=0; x < Nx; x++)
      for (y=0; y < Ny; y++)
      {
         RMS = 0.;
         for (z=0; z < Nz; z++)  
	   RMS +=  (Data(x,y,z) - Ima(y,x))*(Data(x,y,z) - Ima(y,x));
	 if (Nz > 0)  SigmaIma(y,x) = (float) (sqrt(RMS / (float) Nz));
         else SigmaIma(y,x) = 0;
      }
      io_write_ima_float(Name_RMS_Out, SigmaIma);
   }
}

/*********************************************************************/

void coadd_all_frames(fltarray &Dat, Ifloat &Ima, xcorr_def *Xcor, float Zoom,
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

   if (WriteCube == True) io_3d_write_data(Name_Cube_Out, DatOut);

   coadd_all_frames(DatOut, Ima);
}
 
/****************************************************************************/

int main(int argc, char *argv[])
{
    fitsstruct Header;
    fltarray Dat;
    fltarray Mallat;
    int i,k,l, ValRet;
    //Bool OnePSF = True;
    fltarray TabPSF;
    char Cmd[1024];
    Cmd[0] = '\0';
    for (k =0; k < argc; k++) sprintf(Cmd, "%s %s", Cmd, argv[k]);
    Header.origin = Cmd;
   
    /* Get command line arguments, open input file(s) if necessary */
    lm_check(LIC_MR3);
    transinit(argc, argv);
 
    if (Verbose == True)
    {
      cout << endl << endl << "PARAMETERS: " << endl << endl;
      cout << "Name File in = " << Name_Imag_In << endl ;
      cout << "Name File out = " << Name_Imag_Out << endl ;
      cout << "nx = " << Dat.nx() << " ny = " << Dat.ny() << " nz = " << Dat.nz() << endl;
      // cout << "nelem = " << Dat.n_elem() << endl;
      // cout << "Sigma = " << Dat.sigma();
      // cout << " Min   = " << Dat.min() << " Max   = "  << Dat.max() << endl;
      cout << "Interpolation method = " << StringCorrInterp(TypeInterp) << endl;
      if (MaxDist != 4) cout << "Max dist = " << MaxDist  << endl;
      if ((Surface >= 1) || (Surface <= Dat.nx()/2))
      {
         cout << "X Measurement surface parameter = " <<  Surface << endl;
         cout << "Y Measurement surface parameter = " << Surface << endl;
      }
      if ((PosX > 0) && (PosX > 0))
           cout << "XPos = " << PosX << " YPos = " << PosY  << endl;
      if (MeanSub == True) cout << "Subtract to each frame its mean value " << endl;
      if (GetMedian == True) cout << "Take the median instead of the mean " << endl;
      if (OptNoOffset == True) cout << "No offset calculation " << endl;
      if (Noise_STD == True) cout << "Compute the RMS Map " << endl;
      if (WriteCube  == True) cout << "Write the registered cube " << endl;
   }
  
   io_3d_read_data(Name_Imag_In, Dat);
   
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
   int D1 = MIN(PosX-MaxDist,PosY-MaxDist);
   int D2 = MIN(Dat.nx()-PosX-1-MaxDist,Dat.ny()-PosY-1-MaxDist);
   int DM = MIN(D1,D2);
   if (DM < 2)
   {
      cout << "Error: search position cannot be at the border of the image ... " << endl;
      exit(-1);
   } 
   if (Surface >= DM) 
   {
      Surface = DM-1;
      Xcor.xc_hx = Surface;
      Xcor.xc_hy = Surface;
      cout << "Warning: surface must be decreased: new value = " << Surface << endl;
   } 
       
   //Xcor.xc_hx = 10;
   //Xcor.xc_hy = 10;
   Xcor.xc_method = 0;
   char *interp_method= (char*)StringCorrInterp(TypeInterp);
   
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
	 if (Verbose == True) printf("OFFSETS[%d] = %f, %f\n",i,Xcor.xc_x[i],Xcor.xc_y[i]);
      }
      Xcor.xc_init = 1;
   }
   else if (OptNoOffset == False) // otherwise we calculte the offset frame by correlation
   {  
      Ifloat ImaFrame;
      Xcor.xc_init = 0;
      
      // Test position and surface
      
      for (i = 0; i < nz; i++) 
      {
          Xcor.xc_x[i] = Xcor.xc_y[i] = 0.;
	  if (MeanSub == True)
	  {
	     ImaFrame.alloc(TabFrame[i], Dat.ny(), Dat.nx());
	     TabMean[i] = (float) average(ImaFrame);
 	     for (k=0; k < ImaFrame.nl(); k++)
	     for (l=0; l < ImaFrame.nc(); l++)  ImaFrame(k,l) -= TabMean[i];
	  }
      }
      float *Pattern = TabFrame[0];
      
      ValRet = cube_get_offset(TabFrame, Pattern, 
                               Dat.nx(), Dat.ny(), Dat.nz(), 
			       PosX, PosY, &Xcor);
      for (i = 0; i < nz; i++) 
      {
 	  ImaFrame.alloc(TabFrame[i], Dat.ny(), Dat.nx());
	  if (MeanSub == True)
	  {
 	     for (k=0; k < ImaFrame.nl(); k++)
	     for (l=0; l < ImaFrame.nc(); l++)  ImaFrame(k,l) += TabMean[i];
	  }
      }
      // cout <<"ValRet = " << ValRet << endl; 
      if (Verbose == True)
        for (i = 0; i < nz; i++)
           cout << "Offset Frame: " << i+1 << ", dx = " << Xcor.xc_x[i] << " dy = " << Xcor.xc_y[i] << endl;     

      if (OptWOffset == True)
      {
         fltarray TabOffset(2,nz);
	 for (i = 0; i < nz; i++) 
         {
             TabOffset(0,i) = Xcor.xc_x[i];
	     TabOffset(1,i) = Xcor.xc_y[i];
         }
         fits_write_fltarr(Name_Offset, TabOffset);
      }
   }
   
   int Nx1 = (int)(Zoom * Dat.nx());
   int Ny1 = (int)(Zoom * Dat.ny());
   Ifloat Ima(Ny1, Nx1, "ResIma");
   if ((OptNoOffset == False) || (Zoom != 1))
   {
       if (Verbose == True) cout << "Resample cube ... " << endl;
       coadd_all_frames(Dat, Ima, &Xcor, Zoom, interp_method);
   }
   else 
   {
      coadd_all_frames(Dat, Ima);
      if (WriteCube == True)
      {
         cout << "Warning: the registered cube is identical to the input cube ... " << endl;
	 cout << "         the registred cube is not written." << endl;
      }
   }
   io_write_ima_float(Name_Imag_Out, Ima);
   delete [] TabMean;
   delete [] TabFrame;
   exit(0);
}


