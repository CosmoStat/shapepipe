


#include "writefits3d.h"

static void PrintError( int status)
{
	char status_str[FLEN_STATUS], errmsg[FLEN_ERRMSG];
	if (status)
		fprintf(stderr, "\n*** Error occurred during program execution ***\n");
	ffgerr(status, status_str);    
	fprintf(stderr, "\nstatus = %d: %s\n", status, status_str);
	if ( ffgmsg(errmsg) )  
	{
		fprintf(stderr, "\nError message stack:\n");
		fprintf(stderr, " %s\n", errmsg);
		while ( ffgmsg(errmsg) )  
			fprintf(stderr, " %s\n", errmsg);
	}
	exit( status );
}

void writefltarr(char * filename, fltarray &Mat)
{
// Header
	long naxis=3;
	long naxes[3]; naxes[0]=Mat.nx(); naxes[1]=Mat.ny(); naxes[2]=Mat.nz();
	if(naxes[2]==0) naxes[2]=1; // Ifloat
	if(naxes[1]==0) naxes[1]=1; // just in case : may not happen
	int status;
	int simple;
	int bitpix=-32;
	simple   = True;
	long pcount   =   0;  // no group parameters 
	long gcount   =   1;  // only a single image/group 
	int  extend   =   False;
	long group = 1;  // group to write in the fits file, 1= first group 
	fitsfile *fptr;    
	
// File exists
	FILE *FEXIST = fopen(filename, "rb");
	if (FEXIST)
	{
		fclose(FEXIST);
		remove(filename);
	}
	status = 0; 

// open file
	if (ffinit(&fptr, filename, &status))
		PrintError( status );
	
// write the header
	if (ffphpr(fptr,simple,bitpix,naxis,naxes,pcount,gcount,extend,&status) )
		PrintError( status );
	
// write the data
	if (ffppre(fptr, group, 1, naxes[0]*naxes[1]*naxes[2], Mat.buffer(), &status) )
		PrintError( status );
	
// close file
	if (ffclos(fptr, &status) )
		PrintError( status );  
}

void writefltarr(char * filename, dblarray &Mat)
{
// Header
	long naxis=3;
	long naxes[3]; naxes[0]=Mat.nx(); naxes[1]=Mat.ny(); naxes[2]=Mat.nz();
	if(naxes[2]==0) naxes[2]=1; // Ifloat
	if(naxes[1]==0) naxes[1]=1; // just in case : may not happen
	int status;
	int simple;
	int bitpix=-64;
	simple   = True;
	long pcount   =   0;  // no group parameters 
	long gcount   =   1;  // only a single image/group 
	int  extend   =   False;
	long group = 1;  // group to write in the fits file, 1= first group 
	fitsfile *fptr;    
	
// File exists
	FILE *FEXIST = fopen(filename, "rb");
	if (FEXIST)
	{
		fclose(FEXIST);
		remove(filename);
	}
	status = 0; 

// open file
	if (ffinit(&fptr, filename, &status))
		PrintError( status );
	
// write the header
	if (ffphpr(fptr,simple,bitpix,naxis,naxes,pcount,gcount,extend,&status) )
		PrintError( status );
	
// write the data
	if (ffpprd(fptr, group, 1, naxes[0]*naxes[1]*naxes[2], Mat.buffer(), &status) )
		PrintError( status );
	
// close file
	if (ffclos(fptr, &status) )
		PrintError( status );  
}

// write a complex array as a 2times larger (in x) float array
void writecfarr(char * filename, cfarray &Mat)
{
// Header
	long naxis=3;
	long naxes[3]; naxes[0]=2*Mat.nx(); naxes[1]=Mat.ny(); naxes[2]=Mat.nz();
	if(naxes[2]==0) naxes[2]=1; // Ifloat
	if(naxes[1]==0) naxes[1]=1; // just in case : may not happen
	int status;
	int simple;
	int bitpix=-32;
	simple   = True;
	long pcount   =   0;  // no group parameters 
	long gcount   =   1;  // only a single image/group 
	int  extend   =   False;
	long group = 1;  // group to write in the fits file, 1= first group 
	fitsfile *fptr;    
	
// File exists
	FILE *FEXIST = fopen(filename, "rb");
	if (FEXIST)
	{
		fclose(FEXIST);
		remove(filename);
	}
	status = 0; 

// open file
	if (ffinit(&fptr, filename, &status))
		PrintError( status );
	
// write the header
	if (ffphpr(fptr,simple,bitpix,naxis,naxes,pcount,gcount,extend,&status) )
		PrintError( status );
	
// write the data
	if (ffppre(fptr, group, 1, naxes[0]*naxes[1]*naxes[2], (float*)Mat.buffer(), &status) )
		PrintError( status );
	
// close file
	if (ffclos(fptr, &status) )
		PrintError( status );  
}
