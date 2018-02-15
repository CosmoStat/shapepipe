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
**    File:  mr_recons.cc
**
*******************************************************************************
**
**    DESCRIPTION  reconstruction of  an image
**    ----------- 
**                 
**    PARAMETRES    
**    ----------    
** 
**    USAGE: mr_recons  multiresolution_file  output_image
**        multiresolution_file = file (.mr) which contains the multiresolution
**                               transform.
**        output_image = output scale file name
**
**
******************************************************************************/
 
#include "IM_Obj.h"
#include "IM_IO.h"
#include "MR_Obj.h"

char Name_MR_In[256];
char Name_Imag_Out[256];

extern int  OptInd;
extern char *OptArg;

extern int  GetOpt(int argc, char *const*argv, char *opts);
Bool Verbose = False;
 
/*********************************************************************/

static void usage(char *argv[])
{
    fprintf(OUTMAN, "Usage: %s options in_mr_file out_image\n\n", argv[0]);
    fprintf(OUTMAN, "   where options =  \n");
    vm_usage();
    manline();
    verbose_usage();    
    manline();
    exit(-1);
}

/*********************************************************************/

/* GET COMMAND LINE ARGUMENTS */
static void recinit(int argc, char *argv[])
{
    int c;
#ifdef LARGE_BUFF
    int VMSSize=-1;
    Bool OptZ = False;
    char VMSName[1024] = "";
#endif    
    
    /* get options */
    while ((c = GetOpt(argc,argv,"vzZ:")) != -1) 
    {
	switch (c) 
        {
	   case 'v': Verbose = True; break;
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
	if (OptInd < argc) strcpy(Name_MR_In, argv[OptInd++]);
         else usage(argv);

	if (OptInd < argc) strcpy(Name_Imag_Out, argv[OptInd++]);
         else usage(argv);

	/* make sure there are not too many parameters */
	if (OptInd < argc)
        {
		fprintf(stderr, "Too many parameters: %s ...\n", argv[OptInd]);
		usage(argv);
	}
#ifdef LARGE_BUFF
    if (OptZ == True) vms_init(VMSSize, VMSName, Verbose);
#endif
}

/*********************************************************************/

int  main(int argc, char *argv[])
{
    MultiResol MR_Data;

    /* Get command line arguments, open input file(s) if necessary */
    lm_check(LIC_MR1);
    recinit(argc, argv);
    
    if (Verbose == True)
    {
       cout << endl << endl << "PARAMETERS: " << endl << endl;
       cout << "File Name in = " << Name_MR_In << endl;
       cout << "File Name Out = " << Name_Imag_Out << endl;
    }
    MR_Data.Verbose = Verbose;
    MR_Data.read (Name_MR_In);
    if (Verbose == True) MR_Data.print_info();

    Ifloat Dat (MR_Data.size_ima_nl(), MR_Data.size_ima_nc(), "Reconstruct.");
    MR_Data.recons(Dat); 
    
    type_format FormatData = io_which_format(Name_Imag_Out);
    if ((FormatData == F_UNKNOWN) && (MR_Data.FormatInputImag != F_UNKNOWN))  
                    io_set_format(MR_Data.FormatInputImag);
                    
    io_write_ima_float (Name_Imag_Out, Dat);
    exit(0);
} 

